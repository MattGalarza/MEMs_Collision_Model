# ------------------------------------------ Libraries --------------------------------------

using Sundials, ForwardDiff, DifferentialEquations, OrdinaryDiffEq, Plots

# --------------------------------------- Analytical Model ----------------------------------

# Define the AnalyticalModel Module
module AnalyticalModel
using DifferentialEquations
using ForwardDiff
using Parameters
using SpecialFunctions
export Params, p, electrostatic, CoupledSystem!

@with_kw mutable struct Params{T<:Real}
    # Fundamental parameters
    m1::T = 2.0933e-6        # Shuttle mass (kg)
    E::T = 170e9             # Young's modulus (Pa) - Updated to match paper
    eta::T = 1.849e-5        # Dynamic viscosity of air (Pa·s)
    c::T = 0.015             # Damping scaling factor
    g0::T = 14e-6            # Electrode gap (m)
    Tp::T = 120e-9           # Thickness of parylene layer (m)
    Tf::T = 25e-6            # Device thickness (m)
    gss::T = 14e-6           # Soft-stopper initial gap (m)
    rho::T = 2330.0          # Density of silicon (kg/m³) - Updated to match paper
    cp::T = 5e-12            # Capacitance of parylene (F)
    wt::T = 9e-6             # Electrode width, top (m)
    wb::T = 30e-6            # Electrode width, bottom (m)
    ws::T = 14.7e-6          # Suspension spring width (m)
    Lss::T = 1000e-6         # Soft-stopper length (m) - Updated to match paper
    Lf::T = 450e-6           # Full electrode length (m) - Using Lf as in paper
    Leff::T = 400e-6         # Electrode effective overlap length (m)
    Lsp::T = 1400e-6         # Suspension spring length (m) - Added from paper
    e::T = 8.85e-12          # Permittivity of free space (F/m)
    ep::T = 3.2              # Permittivity of parylene
    Vbias::T = 3.0           # Bias voltage (V)
    Rload::T = 0.42e6        # Load resistance (Ω)
    N::Int = 160             # Number of electrodes
    kss::T = 6.0             # Soft-stopper spring force (N/m)
    lambda::T = 70e-9        # Mean free path of air molecules (m) - Added for rarefaction effects
    sigmap::T = 1.016        # Slip coefficient for rarefaction - Added from paper

    # Dependent parameters
    gp::T = :($(g0 - 2 * Tp))                      # Initial gap including parylene layer (m)
    a::T = :($((wb - wt) / Leff))                  # Aspect ratio (tilt factor) - Renamed to match paper
    k1::T = :($((2 / 3) * ((E * Tf * (ws^3)) / (Lsp^3))))  # Linear spring constant (N/m) - Updated per paper
    k3::T = :($((18 / 25) * ((E * Tf * ws) / (Lsp^3))))    # Cubic spring constant (N/m³) - Updated per paper
    I::T = :($((1 / 12) * Tf * (((wt + wb) / 2)^3)))        # Electrode moment of inertia (m⁴) - Using average width
    m2::T = :($((33 / 140) * rho * Tf * Lf * ((wt + wb) / 2)))  # Mass of electrode (kg) - Using average width
    ke::T = :($((E * I) / (Lf^3 / 3)))          # Electrode spring constant (N/m) - Modified calculation
end

function Params{T}(p::Params{S}) where {T<:Real, S<:Real}
    # Extract field names
    fnames = fieldnames(typeof(p))
    # Get field values and convert to type T
    fvalues = T.(getfield.(Ref(p), fnames))
    # Construct a NamedTuple of field values
    params_nt_T = NamedTuple{fnames}(fvalues)
    # Return a new Params{T} instance
    return Params{T}(; params_nt_T...)
end

# Initialize a default Params instance
p = Params{Float64}()

# Suspension spring force, Fsp
function spring(x1, k1, k3, gss, kss)
    # From paper: suspension spring force + soft-stopper force
    # Fsp = -k1 * x1  # Suspension beam force with only linear term
    Fsp = -k1 * x1 - 0.25 * k3 * (x1^3)  # Suspension beam force with cubic term
    
    # Soft-stopper force (piecewise linear as in paper)
    if abs(x1) <= gss
        Fss = 0.0
    else
        Fss = -kss * (abs(x1) - gss) * sign(x1)
    end
    
    # Total spring force
    Fs = Fsp + Fss
    return Fs
end

# Electrode collision force, Fc
function collision(x1, x2, m2, ke, gp)
    # Updated to match paper's piecewise electrode force model
    if abs(x2) < gp  # Non-collision regime
        Fc = -ke * (x1 - x2)  # Simple spring force between shuttle and electrode
        collision_state = "translational"
    else  # Collision regime
        Fc = -ke * (gp * sign(x2) - x1)  # Force acts between shuttle and boundary
        collision_state = "rotational"
    end
    
    return Fc, collision_state
end

# Knudsen number for rarefaction effects
function calculate_knudsen(y, lambda, a)
    # Knudsen number calculation for rarefaction effects
    # Only used in rotational regime
    Kn = lambda / (a * y)
    return Kn
end

# Pressure flow factor for rarefaction
function pressure_flow_factor(Kn, sigmap)
    # From paper: accounts for rarefaction effects
    Qp = 1 / (1 + 6 * sigmap * Kn)
    return Qp
end

# Viscous damping, Fd
function damping(x1, x1dot, x2, x2dot, a, c, gp, Leff, Tf, eta, lambda, sigmap)
    # Determine if in translational or rotational regime
    if abs(x2) < gp  # Translational regime
        # Implementation based on Moy et al. paper equation (9)
        # For MEMS energy harvester with tapered gap geometry
            
        # For our system:
        # - Static bottom plate: ωb = 0, θb = 0
        # - Top plate has no rotation: ωa = 0
        # - x2dot represents dh0/dt (velocity of displacement)
            
        # Calculate parameters for both sides
        # Left side
        h0_l = gp + x2
        hL_l = h0_l + a * Leff
            
        # Right side 
        h0_r = gp - x2
        hL_r = h0_r + a * Leff
            
        # Width (thickness in our case)
        W = Tf
            
        # Calculate integration constants from equations (7) and (8)
        # Left side constants (simplified for ωa = ωb = 0)
        C1_l = (h0_l/(12*eta*a^3*Leff*(h0_l + 2*a))) * (a*Leff*12*eta*x2dot)
        C2_l = (1/(a^4*Leff*(2*h0_l + a*Leff))) * (a*Leff*12*eta*x2dot)
            
            # Right side constants (simplified for ωa = ωb = 0, with opposite velocity)
        C1_r = (h0_r/(12*eta*a^3*Leff*(h0_r + 2*a))) * (a*Leff*12*eta*x2dot)
        C2_r = (1/(a^4*Leff*(2*h0_r + a*Leff))) * (a*Leff*12*eta*x2dot)
            
        # Implement the full force equation (9) from Moy paper for left side
        F_l = W * (
            Leff * C2_l - 
            (6 * eta * Leff) / (a * h0_l * hL_l) * C1_l + 
            # Main α^-4 terms:
            (6 * eta * Leff) / (a^4 * hL_l) * x2dot + 
            (6 * eta / a^4) * (2 * a * x2dot) * log(h0_l) -
            (a * (2 * x2dot)) * log(hL_l)
        )
            
        # Implement the full force equation (9) from Moy paper for right side
        # Note: same x2dot is used but the effect is opposite due to gap geometry
        F_r = W * (
            Leff * C2_r - 
            (6 * eta * Leff) / (a * h0_r * hL_r) * C1_r + 
            # Main α^-4 terms:
            (6 * eta * Leff) / (a^4 * hL_r) * x2dot + 
            (6 * eta / a^4) * (2 * a * x2dot) * log(h0_r) -
            (a * (2 * x2dot)) * log(hL_r)
        )
            
        # Total damping force with damping coefficient
        Fd = -c * (F_l + F_r)
    else  # Rotational regime (collision)
        # For the colliding electrode (using approximation with rarefaction effects)
        y_effective = Leff / 2  # Representative point for average Knudsen calculation
        Kn = calculate_knudsen(y_effective, lambda, a)
        Qp = pressure_flow_factor(Kn, sigmap)
        
        # From paper: rotational damping with rarefaction correction
        Fd_c = (3 * pi * eta * Tf^2 * x1dot * Leff / (8 * x1)) * (
            1 + log((x1 / Tf)^2)
        ) * Qp
        
        # For the non-colliding electrode (use translational formula but with x1 instead of x2)
        non_colliding_side = (x2 > 0) ? "left" : "right"
        if non_colliding_side == "left"
            ul = gp + abs(x2)
            ll = ul + a * Leff
            Fd_nc = (12 * eta * Tf * x1dot / a^3) * (
                (2 * a * Leff) / (2 * (gp + abs(x2)) + a * Leff) + 
                log(abs((gp + abs(x2)) / (gp + abs(x2) + a * Leff)))
            )
        else
            ur = gp + abs(x2)
            lr = ur + a * Leff
            Fd_nc = (12 * eta * Tf * x1dot / a^3) * (
                (2 * a * Leff) / (2 * (gp + abs(x2)) + a * Leff) + 
                log(abs((gp + abs(x2)) / (gp + abs(x2) + a * Leff)))
            )
        end
        
        # Total damping
        Fd = -c * (Fd_c + Fd_nc)
    end
    
    return Fd
end

# Electrostatic coupling, Fe
function electrostatic(x1, x2, Qvar, g0, gp, a, e, ep, cp, wt, wb, ke, E, I, Leff, Tf, Tp, N)
    # Function for variable capacitance, Cvar
    function Cvar_func(x2)
        # Constant dielectric layer capacitance
        crl = (e * ep * Leff * Tf) / Tp
        
        if abs(x2) < gp  # Non-collision regime (translational)
            # Right side air capacitance
            Cair_r = (e * Tf / a) * log((gp - x2 + a * Leff) / (gp - x2))
            Cvar_r = 1 / (2/crl + 1/Cair_r) 
            
            # Left side air capacitance
            Cair_l = (e * Tf / a) * log((gp + x2 + a * Leff) / (gp + x2))
            Cvar_l = 1 / (2/crl + 1/Cair_l)  
            
            # Total variable capacitance
            Cvar_value = (N / 2) * (Cvar_r + Cvar_l)
        else  # Collision regime (rotational)
            # Colliding electrode (use logarithmic approximation from paper)
            # Calculate scaling factor k
            k = 2 * Tp / (a * Leff)
            
            # Determine collision side
            collision_side = (x2 > 0) ? "right" : "left"
            
            # Approximate min and max capacitance values for rotational regime
            Cmin = 7.105299639935359*10^-14
            Cmax = 9.95200974248769*10^-11
            
            if collision_side == "right"
                # Colliding side (right)
                Cair_c = Cmin + (Cmax - Cmin) * (
                    log(1 + k * (abs(x2) - gp)) / log(1 + k * (a * Leff))
                )
                Cvar_c = 1 / (2/crl + 1/Cair_c)
                
                # Non-colliding side (left)
                Cair_nc = ((e * Tf) / (2 * a)) * log((gp + abs(x2) + a * Leff) / (gp + abs(x2)))
                Cvar_nc = 1 / (2/crl + 1/Cair_nc)
            else
                # Colliding side (left)
                Cair_c = Cmin + (Cmax - Cmin) * (
                    log(1 + k * (abs(x2) - gp)) / log(1 + k * (a * Leff))
                )
                Cvar_c = 1 / (2/crl + 1/Cair_c)
                
                # Non-colliding side (right)
                Cair_nc = (e * Tf / a) * log((gp + abs(x2) + a * Leff) / (gp + abs(x2)))
                Cvar_nc = 1 / (2/crl + 1/Cair_nc)
            end
            
            # Total variable capacitance
            Cvar_value = (N / 2) * (Cvar_c + Cvar_nc)
        end
        
        return Cvar_value
    end
    
    # Compute Cvar and its derivative
    Cvar = Cvar_func(x2)
    dC = ForwardDiff.derivative(Cvar_func, x2)
    
    # Total capacitance, Ctotal
    Ctotal = Cvar + cp
    
    # Electrostatic force, Fe 
    Fe = -((Qvar^2 / (2 * Ctotal^2)) * dC) / (N / 2)
    
    return Ctotal, Fe
end

function CoupledSystem!(dz, z, p, t, current_acceleration)
    # Unpack state variables
    z1, z2, z3, z4, z5, Vout = z  # x1, x1dot, x2, x2dot, Qvar, Vout
    
    # Compute forces
    Fs = spring(z1, p.k1, p.k3, p.gss, p.kss)
    Fc, collision_state = collision(z1, z3, p.m2, p.ke, p.gp)
    
    # Determine which velocity to use for damping based on collision state
    if collision_state == "translational"
        x_dot_for_damping = z4  # Use electrode velocity for non-collision
        Fd = damping(z1, z2, z3, z4, p.a, p.c, p.gp, p.Leff, p.Tf, p.eta, p.lambda, p.sigmap)
    else
        x_dot_for_damping = z2  # Use shuttle velocity for collision (rotation)
        Fd = damping(z1, z2, z3, z2, p.a, p.c, p.gp, p.Leff, p.Tf, p.eta, p.lambda, p.sigmap)
    end
    
    Ctotal, Fe = electrostatic(z1, z3, z5, p.g0, p.gp, p.a, p.e, p.ep, p.cp, p.wt, p.wb, p.ke, p.E, p.I, p.Leff, p.Tf, p.Tp, p.N)
    
    # Use current_acceleration as the external force
    Fext = current_acceleration
    
    # Compute derivatives
    dz[1] = z2  
    dz[2] = (Fs + (p.N / 2) * Fc) / p.m1 - Fext
    dz[3] = z4  
    dz[4] = (-Fc + Fd + Fe) / p.m2 - Fext
    dz[5] = (p.Vbias - (z5 / Ctotal)) / p.Rload  
    dz[6] = (p.Vbias - z5 / Ctotal - Vout) / (p.Rload * Ctotal)
end
end

# Import the AnalyticalModel module without bringing all exports into Main
import .AnalyticalModel  

# Function to convert Params to NamedTuple
function params_to_namedtuple(p)
    param_names = fieldnames(typeof(p))
    param_values = getfield.(Ref(p), param_names)
    return NamedTuple{param_names}(param_values)
end

# --------------------------------------- External Force ------------------------------------

# Sine Wave External Force
f = 20.0 # Frequency (Hz)
alpha = 3.1 # Applied acceleration constant (dimensionless)
g = 9.81 # Gravitational constant (m/s²)
A = alpha * g
t_ramp = 0.2 # Ramp-up duration (s)
# Define the ramp function (linear ramp)
ramp(t) = t < t_ramp ? t / t_ramp : 1.0
Fext_sine = t -> A * ramp(t) * sin(2 * π * f * t)

# ------------------------------------- Set Input Force ------------------------------------

# Set to `true` to use sine wave, `false` for displaced IC
use_sine = true

# Define Fext_input based on your choice
if use_sine
    Fext_input = Fext_sine
else
    Fext_input = t -> 0.0
end

# ------------------------------------ Initialize Parameters --------------------------------

# Create a new Params instance by copying the default
p_new = deepcopy(AnalyticalModel.p)

# Initial conditions
x10 = 0.0 # Initial displacement
x10dot = 0.0 # Initial velocity
x20 = 0.0 # Initial displacement
x20dot = 0.0 # Initial velocity

# Compute initial electrostatic parameters
Ctotal0, Fe0 = AnalyticalModel.electrostatic(
    x10, x20, 0.0, 
    p_new.g0, p_new.gp, p_new.a, 
    p_new.e, p_new.ep, p_new.cp, 
    p_new.wt, p_new.wb, p_new.ke, 
    p_new.E, p_new.I, p_new.Leff, 
    p_new.Tf, p_new.Tp, p_new.N
)
Q0 = p_new.Vbias * Ctotal0 # Initial charge
Vout0 = p_new.Vbias - (Q0 / Ctotal0) # Initial voltage
z0 = [x10, x10dot, x20, x20dot, Q0, Vout0]
tspan = (0.0, 0.5) # simulation length
# teval = () # evaluation steps
abstol = 1e-9 # absolute solver tolerance
reltol = 1e-6 # relative solver tolerance

# ---------------------------------- Solve Analytical Model ---------------------------------

# Define a wrapper function
function CoupledSystem_wrapper!(dz, z, p, t)
    current_acceleration = Fext_input(t)
    AnalyticalModel.CoupledSystem!(dz, z, p, t, current_acceleration)
end

# Define and solve the ODE problem
eqn = ODEProblem(CoupledSystem_wrapper!, z0, tspan, p_new)

# Solve the system using Rosenbrock23 solver
sol = solve(eqn, Rosenbrock23(); abstol=abstol, reltol=reltol, maxiters=1e7)
# If the system is too stiff, use CVODE_BDF from Sundials
# sol = solve(eqn, CVODE_BDF(), abstol=abstol, reltol=reltol, maxiters=Int(1e9))

# Verify the solution structure
println("Type of sol.u: ", typeof(sol.u))
println("Size of sol.u: ", size(sol.u))
println("Solver status: ", sol.retcode)

# ----------------------------------------- Plotting -----------------------------------------

# Plot the related terms
x1 = [u[1] for u in sol.u]
x1dot = [u[2] for u in sol.u]
x2 = [u[3] for u in sol.u]
x2dot = [u[4] for u in sol.u]
Qvar = [u[5] for u in sol.u]
V = [u[6] for u in sol.u]

p3 = plot(sol.t, x1, xlabel = "Time (s)", ylabel = "x1 (m)", title = "Shuttle Mass Displacement (x1)")
display(p3)
p4 = plot(sol.t, x1dot, xlabel = "Time (s)", ylabel = "x1dot (m/s)", title = "Shuttle Mass Velocity (x1dot)")
display(p4)
p5 = plot(sol.t, x2, xlabel = "Time (s)", ylabel = "x2 (m)", title = "Mobile Electrode Displacement (x2)")
display(p5)
p6 = plot(sol.t, x2dot, xlabel = "Time (s)", ylabel = "x2dot (m/s)", title = "Mobile Electrode Velocity (x2)")
display(p6)
p7 = plot(sol.t, Qvar, xlabel = "Time (s)", ylabel = "Qvar (C)", title = "Charge (Qvar)")
display(p7)
p8 = plot(sol.t, V, xlabel = "Time (s)", ylabel = "Vout (V)", title = "Output Voltage")
display(p8)

# Generate forces during the simulation
# Initialize arrays to store forces
Fs_array = Float64[] # Suspension spring force
Fc_array = Float64[] # Collision force
Fd_array = Float64[] # Damping force
Fe_array = Float64[] # Electrostatic force
CollisionState_array = String[] # Track collision state

# Iterate over each solution point to compute forces
for (i, t) in enumerate(sol.t)
    # Extract state variables at time t
    z = sol.u[i]
    z1, z2, z3, z4, z5, Vout = z
    
    # Compute Fs (Suspension spring force)
    Fs = AnalyticalModel.spring(z1, p_new.k1, p_new.k3, p_new.gss, p_new.kss)
    push!(Fs_array, Fs)
    
    # Compute Fc (Collision force)
    Fc, collision_state = AnalyticalModel.collision(z1, z3, p_new.m2, p_new.ke, p_new.gp)
    push!(Fc_array, Fc)
    push!(CollisionState_array, collision_state)
    
    # Compute Fd (Damping force)
    # Use appropriate velocity based on collision state
    if collision_state == "translational"
        Fd = AnalyticalModel.damping(z1, z2, z3, z4, p_new.a, p_new.c, p_new.gp, p_new.Leff, p_new.Tf, p_new.eta, p_new.lambda, p_new.sigmap)
    else
        Fd = AnalyticalModel.damping(z1, z2, z3, z2, p_new.a, p_new.c, p_new.gp, p_new.Leff, p_new.Tf, p_new.eta, p_new.lambda, p_new.sigmap)
    end
    push!(Fd_array, Fd)
    
    # Compute Fe (Electrostatic force)
    Ctotalx, Fe = AnalyticalModel.electrostatic(z1, z3, z5, p_new.g0, p_new.gp, p_new.a, 
                                            p_new.e, p_new.ep, p_new.cp, p_new.wt, 
                                            p_new.wb, p_new.ke, p_new.E, p_new.I, 
                                            p_new.Leff, p_new.Tf, p_new.Tp, p_new.N)
    push!(Fe_array, Fe)
end

# Plotting respective forces
p9 = plot(sol.t, Fs_array, xlabel = "Time (s)", ylabel = "Fs (N)", title = "Suspension + Soft-stopper Spring Force")
display(p9)
p10 = plot(sol.t, Fc_array, xlabel = "Time (s)", ylabel = "Fc (N)", title = "Mobile Electrode Collision Force")
display(p10)
p11 = plot(sol.t, Fd_array, xlabel = "Time (s)", ylabel = "Fd (N)", title = "Viscous Damping Force")
display(p11)
p12 = plot(sol.t, Fe_array, xlabel = "Time (s)", ylabel = "Fe (N)", title = "Electrostatic Force")
display(p12)
p13 = plot(sol.t, Fext_input, xlabel = "Time (s)", ylabel = "Fext (N)", title = "Applied External Force")
display(p13)

# Plot collision state transitions
p14 = plot(sol.t, [s == "translational" ? 0 : 1 for s in CollisionState_array], 
          xlabel = "Time (s)", ylabel = "State", yticks = ([0, 1], ["Translational", "Rotational"]), 
          title = "Collision State Transitions")
display(p14)

function plot_combined_results(sol, x1, x1dot, x2, x2dot, Qvar, V,
    Fs_array, Fc_array, Fd_array, Fe_array, CollisionState_array)
# Create a plot for all state variables
p_states = plot(layout=(3,2), size=(1000, 800), dpi=150)

# Shuttle displacement
plot!(p_states[1], sol.t, x1 .* 1e6, 
xlabel = "", ylabel = "x₁ (μm)", 
title = "Shuttle Displacement", 
legend = false, linewidth = 2, color = :blue)

# Shuttle velocity
plot!(p_states[2], sol.t, x1dot .* 1e3, 
xlabel = "", ylabel = "x₁dot (mm/s)", 
title = "Shuttle Velocity", 
legend = false, linewidth = 2, color = :red)

# Electrode displacement
plot!(p_states[3], sol.t, x2 .* 1e6, 
xlabel = "", ylabel = "x₂ (μm)", 
title = "Electrode Displacement", 
legend = false, linewidth = 2, color = :green)

# Electrode velocity
plot!(p_states[4], sol.t, x2dot .* 1e3, 
xlabel = "", ylabel = "x₂dot (mm/s)", 
title = "Electrode Velocity", 
legend = false, linewidth = 2, color = :purple)

# Charge
plot!(p_states[5], sol.t, Qvar .* 1e12, 
xlabel = "Time (s)", ylabel = "Q (pC)", 
title = "Charge", 
legend = false, linewidth = 2, color = :brown)

# Output voltage
plot!(p_states[6], sol.t, V .* 1e3, 
xlabel = "Time (s)", ylabel = "Vout (mV)", 
title = "Output Voltage", 
legend = false, linewidth = 2, color = :orange)

# Create a plot for all forces
p_forces = plot(layout=(3,2), size=(1000, 800), dpi=150)

# Spring force
plot!(p_forces[1], sol.t, Fs_array .* 1e6, 
xlabel = "", ylabel = "Fs (μN)", 
title = "Spring Force", 
legend = false, linewidth = 2, color = :blue)

# Collision force
plot!(p_forces[2], sol.t, Fc_array .* 1e6, 
xlabel = "", ylabel = "Fc (μN)", 
title = "Collision Force", 
legend = false, linewidth = 2, color = :red)

# Damping force
plot!(p_forces[3], sol.t, Fd_array .* 1e6, 
xlabel = "", ylabel = "Fd (μN)", 
title = "Damping Force", 
legend = false, linewidth = 2, color = :green)

# Electrostatic force
plot!(p_forces[4], sol.t, Fe_array .* 1e6, 
xlabel = "", ylabel = "Fe (μN)", 
title = "Electrostatic Force", 
legend = false, linewidth = 2, color = :purple)

# Collision state
collision_numeric = [s == "translational" ? 0 : 1 for s in CollisionState_array]
plot!(p_forces[6], sol.t, collision_numeric, 
xlabel = "Time (s)", ylabel = "State", 
yticks = ([0, 1], ["Trans.", "Rot."]), 
title = "Collision State", 
legend = false, linewidth = 2, color = :black)

return p_states, p_forces
end

# Create the combined plots
p_states, p_forces = plot_combined_results(
sol, x1, x1dot, x2, x2dot, Qvar, V,
Fs_array, Fc_array, Fd_array, Fe_array, CollisionState_array
)

# Display the combined plots
display(p_states)
display(p_forces)