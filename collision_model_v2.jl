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
    # Fundamental geometric parameters
    g0::T = 14e-6  # Initial gap
    Tp::T = 120e-9  # Parylene-C thickness
    Tf::T = 25e-6  # Electrode thickness
    wt::T = 9e-6  # Electrode width, top
    wb::T = 30e-6  # Electrode width, bottom
    ws::T = 14.7e-6  # Suspension spring width
    wss::T = 14e-6  # Soft-stopper width
    Leff::T = 400e-6  # Effective electrode length
    Lf::T = 450e-6  # Full electrode length
    Lsp::T = 1400e-6  # Suspension spring length
    Lss::T = 1000e-6  # Soft-stopper length
    gss::T = 14e-6  # Soft-stopper position
    
    # Mass and material properties
    m1::T = 2.0933e-6  # Shuttle mass
    rho::T = 2330.0  # Density of silicon
    E::T = 170e9  # Young's modulus
    e::T = 8.85e-12  # Permittivity of free space
    ep::T = 3.2  # Relative permittivity of Parylene-C
    eta::T = 1.849e-5  # Viscosity of air
    c::T = 0.015  # Damping coefficient
    lambda::T = 70e-9        # Mean free path of air molecules (m)
    sigmap::T = 1.016        # Slip coefficient for rarefaction
    
    # Electrical parameters
    N::Int = 160  # Number of electrodes
    cp::T = 5e-12  # Capacitance of Parylene-C
    Vbias::T = 3.0  # Bias voltage
    Rload::T = 0.42e6  # Load resistance
    
    # Derived parameters - these will be calculated by create_params()
    gp::T = :($(g0 - 2 * Tp))  # Initial electrode gap
    a::T = :($((wb - wt) / Leff))  # Taper ratio
    m2::T = 0.0  # Modal mass of electrode
    I::T = 0.0  # Electrode moment of inertia
    ke::T = 0.0  # Electrode spring constant
    k1::T = 0.0  # Linear spring constant
    k3::T = 0.0  # Cubic spring constant
    kss::T = 6.0  # Soft-stopper spring constant
end

# Function to create and update dependent parameters for a Params instance
function create_params(p::Params{T}; verbose=true) where T<:Real
    # Calculate electrode moment of inertia
    p.I = (1/48) * p.Lf * p.Tf * (p.wb + p.wt) * (p.wb^2 + p.wt^2)
    
    # Calculate and set modal mass, m2
    modalcoeff = 0.236 + 0.045 * (1.0 - p.wt / p.wb)
    m2Physical = 0.5 * p.Lf^2 * p.rho * p.Tf * (p.wb + p.wt)
    p.m2 = modalcoeff * m2Physical
    
    # Calculate electrode spring constant, ke
    F = 1 # Normalized force used in ke derivation
    num = (p.E * p.Tf * p.wt^2 * ((p.wb - p.wt)^3))
    dem = 6 * F * p.Lf^3 * ((p.wb - 3 * p.wt) * (p.wb - p.wt) + 2 * p.wt^2 * (log(p.Lf * p.wb) - log(p.Lf * p.wt)))
    p.ke = num / dem
    
    # Calculate spring constants, k1/k3/kss
    p.k1 = (4.0/6.0) * ((p.E * p.Tf * (p.ws^3)) / (p.Lsp^3))
    p.k3 = (18.0/25.0) * ((p.E * p.Tf * p.ws) / (p.Lsp^3))
    p.kss = (p.E * p.Tf * (p.wss^3)) / (p.Lss^3)

    # Print calculation details if verbose is true
    if verbose
        println("\n--- Modal Mass Calculation ---")
        println("Modal coefficient terms:  ", modalcoeff)
        println("Physical mass terms:  ", m2Physical)
        println("Modal mass, m2:  ", p.m2)
        
        println("\n--- Electrode Spring Constant Calculation ---")
        println("Tilt angle, a:  ", p.a)
        println("Moment of inertia, I:  ", p.I)
        println("Electrode spring constant, ke:  ", p.ke)
        
        println("\n--- Spring Constants ---")
        println("Linear suspension spring constant, k1:  ", p.k1)
        println("Cubic suspension spring constant, k3:  ", p.k3)
        println("Soft-stopper spring constant, kss:  ", p.kss)
    end
    
    return p
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

# Initialize a default Params instance and calculate dependent parameters
p = Params{Float64}()
p = create_params(p; verbose=false)

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

# --------------------------------------- Fancy Plots ---------------------------------------
# Journal-Quality Individual Plots for State Variables and Forces
using Plots, LaTeXStrings, Measures

# Define high-quality theme for journal publication
function set_journal_theme()
    default(
        fontfamily="Computer Modern",  # LaTeX-like font
        linewidth=2.5,                 # Thicker lines
        foreground_color_legend=nothing, # Transparent legend background
        background_color_legend=nothing, # Transparent legend background
        legendfontsize=10,             # Legend font size
        guidefontsize=12,              # Axis label font size
        tickfontsize=10,               # Tick label font size
        titlefontsize=14,              # Title font size
        size=(800, 600),               # Figure size
        dpi=600,                       # High DPI for print quality
        margin=10mm,                   # Margins around the plot
        grid=false,                    # No grid by default
        framestyle=:box,               # Box-style frame
        foreground_color_axis=:black,  # Black axes
        tick_direction=:out,           # Ticks pointing outward
        palette=:default               # Color palette
    )
end

# Apply the theme
set_journal_theme()

# 1. Shuttle Displacement (x1)
p_x1 = plot(
    sol.t, 
    x1 .* 1e6,  # Convert to μm for better visualization
    linewidth=1.0,
    color=:royalblue3,
    xlabel=L"Time $t$ (s)",
    ylabel=L"Shuttle Displacement $x_1$ ($\mu$m)",
    title="",  # Journal plots often don't have titles
    framestyle=:box,
    legend=false,
    dpi=600,
    size=(800, 500),
    left_margin=10mm,
    bottom_margin=10mm,
    top_margin=5mm,
    right_margin=5mm
)
savefig(p_x1, "shuttle_displacement.pdf")  # PDF for publication
savefig(p_x1, "shuttle_displacement.png")  # PNG for quick viewing

# 2. Shuttle Velocity (x1dot)
p_x1dot = plot(
    sol.t, 
    x1dot .* 1e3,  # Convert to mm/s
    linewidth=1.5,
    color=:firebrick3,
    xlabel=L"Time $t$ (s)",
    ylabel=L"Shuttle Velocity $\dot{x}_1$ (mm/s)",
    title="",
    framestyle=:box,
    legend=false,
    dpi=600,
    size=(800, 500),
    left_margin=10mm,
    bottom_margin=10mm,
    top_margin=5mm,
    right_margin=5mm
)
savefig(p_x1dot, "shuttle_velocity.pdf")
savefig(p_x1dot, "shuttle_velocity.png")

# 3. Electrode Displacement (x2)
p_x2 = plot(
    sol.t, 
    x2 .* 1e6,  # Convert to μm
    linewidth=1.5,
    color=:green4,
    xlabel=L"Time $t$ (s)",
    ylabel=L"Electrode Displacement $x_2$ ($\mu$m)",
    title="",
    framestyle=:box,
    legend=false,
    dpi=600,
    size=(800, 500),
    left_margin=10mm,
    bottom_margin=10mm,
    top_margin=5mm,
    right_margin=5mm
)
savefig(p_x2, "electrode_displacement.pdf")
savefig(p_x2, "electrode_displacement.png")

# 4. Electrode Velocity (x2dot)
p_x2dot = plot(
    sol.t, 
    x2dot .* 1e3,  # Convert to mm/s
    linewidth=1.5,
    color=:purple3,
    xlabel=L"Time $t$ (s)",
    ylabel=L"Electrode Velocity $\dot{x}_2$ (mm/s)",
    title="",
    framestyle=:box,
    legend=false,
    dpi=600,
    size=(800, 500),
    left_margin=10mm,
    bottom_margin=10mm,
    top_margin=5mm,
    right_margin=5mm
)
savefig(p_x2dot, "electrode_velocity.pdf")
savefig(p_x2dot, "electrode_velocity.png")

# 5. Charge (Qvar)
p_qvar = plot(
    sol.t, 
    Qvar .* 1e12,  # Convert to pC
    linewidth=1.5,
    color=:darkorange,
    xlabel=L"Time $t$ (s)",
    ylabel=L"Charge $Q_{var}$ (pC)",
    title="",
    framestyle=:box,
    legend=false,
    dpi=600,
    size=(800, 500),
    left_margin=10mm,
    bottom_margin=10mm,
    top_margin=5mm,
    right_margin=5mm
)
savefig(p_qvar, "charge.pdf")
savefig(p_qvar, "charge.png")

# 6. Output Voltage (V)
p_vout = plot(
    sol.t, 
    V .* 1e3,  # Convert to mV
    linewidth=1.5,
    color=:darkorange2,
    xlabel=L"Time $t$ (s)",
    ylabel=L"Output Voltage $V_{out}$ (mV)",
    title="",
    framestyle=:box,
    legend=false,
    dpi=600,
    size=(800, 500),
    left_margin=10mm,
    bottom_margin=10mm,
    top_margin=5mm,
    right_margin=5mm
)
savefig(p_vout, "output_voltage.pdf")
savefig(p_vout, "output_voltage.png")

# 7. Relative Displacement (x1-x2)
rel_disp = x1 - x2  # Calculate relative displacement
p_rel = plot(
    sol.t, 
    rel_disp .* 1e6,  # Convert to μm
    linewidth=1.5,
    color=:teal,
    xlabel=L"Time $t$ (s)",
    ylabel=L"Relative Displacement $x_1 - x_2$ ($\mu$m)",
    title="",
    framestyle=:box,
    legend=false,
    dpi=600,
    size=(800, 500),
    left_margin=10mm,
    bottom_margin=10mm,
    top_margin=5mm,
    right_margin=5mm
)
savefig(p_rel, "relative_displacement.pdf")
savefig(p_rel, "relative_displacement.png")

# 1. Suspension Spring Force (Fs)
p_fs = plot(
    sol.t, 
    Fs_array .* 1e6,  # Convert to μN
    linewidth=1.5,
    color=:royalblue3,
    xlabel=L"Time $t$ (s)",
    ylabel=L"Spring Force $F_s$ ($\mu$N)",
    title="",
    framestyle=:box,
    legend=false,
    dpi=600,
    size=(800, 500),
    left_margin=10mm,
    bottom_margin=10mm,
    top_margin=5mm,
    right_margin=5mm
)
savefig(p_fs, "spring_force.pdf")
savefig(p_fs, "spring_force.png")

# 2. Collision Force (Fc)
p_fc = plot(
    sol.t, 
    Fc_array .* 1e6,  # Convert to μN
    linewidth=1.5,
    color=:firebrick3,
    xlabel=L"Time $t$ (s)",
    ylabel=L"Collision Force $F_c$ ($\mu$N)",
    title="",
    framestyle=:box,
    legend=false,
    dpi=600,
    size=(800, 500),
    left_margin=10mm,
    bottom_margin=10mm,
    top_margin=5mm,
    right_margin=5mm
)
savefig(p_fc, "collision_force.pdf")
savefig(p_fc, "collision_force.png")

# 3. Damping Force (Fd)
p_fd = plot(
    sol.t, 
    Fd_array .* 1e6,  # Convert to μN
    linewidth=1.5,
    color=:green4,
    xlabel=L"Time $t$ (s)",
    ylabel=L"Damping Force $F_d$ ($\mu$N)",
    title="",
    framestyle=:box,
    legend=false,
    dpi=600,
    size=(800, 500),
    left_margin=10mm,
    bottom_margin=10mm,
    top_margin=5mm,
    right_margin=5mm
)
savefig(p_fd, "damping_force.pdf")
savefig(p_fd, "damping_force.png")

# 4. Electrostatic Force (Fe)
p_fe = plot(
    sol.t, 
    Fe_array .* 1e6,  # Convert to μN
    linewidth=1.5,
    color=:purple3,
    xlabel=L"Time $t$ (s)",
    ylabel=L"Electrostatic Force $F_e$ ($\mu$N)",
    title="",
    framestyle=:box,
    legend=false,
    dpi=600,
    size=(800, 500),
    left_margin=10mm,
    bottom_margin=10mm,
    top_margin=5mm,
    right_margin=5mm
)
savefig(p_fe, "electrostatic_force.pdf")
savefig(p_fe, "electrostatic_force.png")

# 5. External Force (Fext)
p_fext = plot(
    sol.t, 
    [Fext_input(t) for t in sol.t] .* 1e6,  # Calculate and convert to μN
    linewidth=1.5,
    color=:darkorange2,
    xlabel=L"Time $t$ (s)",
    ylabel=L"External Force $F_{ext}$ ($\mu$N)",
    title="",
    framestyle=:box,
    legend=false,
    dpi=600,
    size=(800, 500),
    left_margin=10mm,
    bottom_margin=10mm,
    top_margin=5mm,
    right_margin=5mm
)
savefig(p_fext, "external_force.pdf")
savefig(p_fext, "external_force.png")

# 6. Collision State
collision_numeric = [s == "translational" ? 0 : 1 for s in CollisionState_array]
p_cs = plot(
    sol.t, 
    collision_numeric,
    linewidth=1.5,
    color=:black,
    xlabel=L"Time $t$ (s)",
    ylabel="Collision State",
    yticks=([0, 1], ["Translational", "Rotational"]),
    title="",
    framestyle=:box,
    legend=false,
    dpi=600,
    size=(800, 500),
    left_margin=10mm,
    bottom_margin=10mm,
    top_margin=5mm,
    right_margin=5mm
)
savefig(p_cs, "collision_state.pdf")
savefig(p_cs, "collision_state.png")

println("All individual journal-quality plots created successfully!")


# ----------------------------------------- Short Fancy Plots -----------------------------------------
# Journal-Quality Individual Plots for State Variables and Forces (Optimized Oscillations)
using Plots, LaTeXStrings, Measures, Statistics

# Define high-quality theme for journal publication
function set_journal_theme()
    default(
        fontfamily="Computer Modern",  # LaTeX-like font
        linewidth=1.25,                 # Thicker lines
        foreground_color_legend=nothing, # Transparent legend background
        background_color_legend=nothing, # Transparent legend background
        legendfontsize=10,             # Legend font size
        guidefontsize=12,              # Axis label font size
        tickfontsize=10,               # Tick label font size
        titlefontsize=14,              # Title font size
        size=(800, 600),               # Figure size
        dpi=600,                       # High DPI for print quality
        grid=false,                    # No grid by default
        framestyle=:box,               # Box-style frame
        foreground_color_axis=:black,  # Black axes
        tick_direction=:out,           # Ticks pointing outward
        palette=:default               # Color palette
    )
end

# Apply the theme
set_journal_theme()

# Fixed time window for exactly 4 oscillations based on the image you shared
start_time = 0.275    # Same start time as in the image
end_time = 0.327     # End time calibrated from the image you provided

# Find indices corresponding to the desired time range
start_idx = findfirst(t -> t >= start_time, sol.t)
end_idx = findfirst(t -> t >= end_time, sol.t)

if start_idx === nothing
    start_idx = 1
end
if end_idx === nothing
    end_idx = length(sol.t)
end

# Extract data for the specified time window
plot_time = sol.t[start_idx:end_idx]
plot_x1 = x1[start_idx:end_idx]
plot_x1dot = x1dot[start_idx:end_idx]
plot_x2 = x2[start_idx:end_idx]
plot_x2dot = x2dot[start_idx:end_idx]
plot_Qvar = Qvar[start_idx:end_idx]
plot_V = V[start_idx:end_idx]
plot_rel_disp = plot_x1 - plot_x2
plot_Fs = Fs_array[start_idx:end_idx]
plot_Fc = Fc_array[start_idx:end_idx]
plot_Fd = Fd_array[start_idx:end_idx]
plot_Fe = Fe_array[start_idx:end_idx]
plot_Fext = [Fext_input(t) for t in plot_time]
plot_collision_state = CollisionState_array[start_idx:end_idx]
plot_collision_numeric = [s == "translational" ? 0 : 1 for s in plot_collision_state]

# Print the selected time range
println("Selected time range: $(plot_time[1]) to $(plot_time[end]) seconds")
println("Number of data points: $(length(plot_time))")

# 1. Shuttle Displacement (x1)
p_x1 = plot(
    plot_time, 
    plot_x1 .* 1e6,  # Convert to μm for better visualization
    linewidth=1.25,
    color=:royalblue3,
    xlabel=L"Time $t$ (s)",
    ylabel=L"Shuttle Displacement $x_1$ ($\mu$m)",
    title="",  # Journal plots often don't have titles
    framestyle=:box,
    legend=false,
    dpi=600,
    size=(800, 500),
    left_margin=10mm,
    bottom_margin=10mm,
    top_margin=5mm,
    right_margin=5mm
)
savefig(p_x1, "shuttle_displacement_opt.pdf")
savefig(p_x1, "shuttle_displacement_opt.png")

# 2. Shuttle Velocity (x1dot)
p_x1dot = plot(
    plot_time, 
    plot_x1dot .* 1e3,  # Convert to mm/s
    linewidth=1.25,
    color=:firebrick3,
    xlabel=L"Time $t$ (s)",
    ylabel=L"Shuttle Velocity $\dot{x}_1$ (mm/s)",
    title="",
    framestyle=:box,
    legend=false,
    dpi=600,
    size=(800, 500),
    left_margin=10mm,
    bottom_margin=10mm,
    top_margin=5mm,
    right_margin=5mm
)
savefig(p_x1dot, "shuttle_velocity_opt.pdf")
savefig(p_x1dot, "shuttle_velocity_opt.png")

# 3. Electrode Displacement (x2)
p_x2 = plot(
    plot_time, 
    plot_x2 .* 1e6,  # Convert to μm
    linewidth=1.25,
    color=:green4,
    xlabel=L"Time $t$ (s)",
    ylabel=L"Electrode Displacement $x_2$ ($\mu$m)",
    title="",
    framestyle=:box,
    legend=false,
    dpi=600,
    size=(800, 500),
    left_margin=10mm,
    bottom_margin=10mm,
    top_margin=5mm,
    right_margin=5mm
)
savefig(p_x2, "electrode_displacement_opt.pdf")
savefig(p_x2, "electrode_displacement_opt.png")

# 4. Electrode Velocity (x2dot)
p_x2dot = plot(
    plot_time, 
    plot_x2dot .* 1e3,  # Convert to mm/s
    linewidth=1.25,
    color=:purple3,
    xlabel=L"Time $t$ (s)",
    ylabel=L"Electrode Velocity $\dot{x}_2$ (mm/s)",
    title="",
    framestyle=:box,
    legend=false,
    dpi=600,
    size=(800, 500),
    left_margin=10mm,
    bottom_margin=10mm,
    top_margin=5mm,
    right_margin=5mm
)
savefig(p_x2dot, "electrode_velocity_opt.pdf")
savefig(p_x2dot, "electrode_velocity_opt.png")

# 5. Charge (Qvar)
p_qvar = plot(
    plot_time, 
    plot_Qvar .* 1e12,  # Convert to pC
    linewidth=1.25,
    color=:darkorange,
    xlabel=L"Time $t$ (s)",
    ylabel=L"Charge $Q_{var}$ (pC)",
    title="",
    framestyle=:box,
    legend=false,
    dpi=600,
    size=(800, 500),
    left_margin=10mm,
    bottom_margin=10mm,
    top_margin=5mm,
    right_margin=5mm
)
savefig(p_qvar, "charge_opt.pdf")
savefig(p_qvar, "charge_opt.png")

# 6. Output Voltage (V)
p_vout = plot(
    plot_time, 
    plot_V .* 1e3,  # Convert to mV
    linewidth=1.25,
    color=:darkorange2,
    xlabel=L"Time $t$ (s)",
    ylabel=L"Output Voltage $V_{out}$ (mV)",
    title="",
    framestyle=:box,
    legend=false,
    dpi=600,
    size=(800, 500),
    left_margin=10mm,
    bottom_margin=10mm,
    top_margin=5mm,
    right_margin=5mm
)
savefig(p_vout, "output_voltage_opt.pdf")
savefig(p_vout, "output_voltage_opt.png")

# 7. Relative Displacement (x1-x2)
p_rel = plot(
    plot_time, 
    plot_rel_disp .* 1e6,  # Convert to μm
    linewidth=1.25,
    color=:teal,
    xlabel=L"Time $t$ (s)",
    ylabel=L"Relative Displacement $x_1 - x_2$ ($\mu$m)",
    title="",
    framestyle=:box,
    legend=false,
    dpi=600,
    size=(800, 500),
    left_margin=10mm,
    bottom_margin=10mm,
    top_margin=5mm,
    right_margin=5mm
)
savefig(p_rel, "relative_displacement_opt.pdf")
savefig(p_rel, "relative_displacement_opt.png")

# 1. Suspension Spring Force (Fs)
p_fs = plot(
    plot_time, 
    plot_Fs .* 1e6,  # Convert to μN
    linewidth=1.25,
    color=:royalblue3,
    xlabel=L"Time $t$ (s)",
    ylabel=L"Spring Force $F_s$ ($\mu$N)",
    title="",
    framestyle=:box,
    legend=false,
    dpi=600,
    size=(800, 500),
    left_margin=10mm,
    bottom_margin=10mm,
    top_margin=5mm,
    right_margin=5mm
)
savefig(p_fs, "spring_force_opt.pdf")
savefig(p_fs, "spring_force_opt.png")

# 2. Collision Force (Fc)
p_fc = plot(
    plot_time, 
    plot_Fc .* 1e6,  # Convert to μN
    linewidth=1.25,
    color=:firebrick3,
    xlabel=L"Time $t$ (s)",
    ylabel=L"Collision Force $F_c$ ($\mu$N)",
    title="",
    framestyle=:box,
    legend=false,
    dpi=600,
    size=(800, 500),
    left_margin=10mm,
    bottom_margin=10mm,
    top_margin=5mm,
    right_margin=5mm
)
savefig(p_fc, "collision_force_opt.pdf")
savefig(p_fc, "collision_force_opt.png")

# 3. Damping Force (Fd)
p_fd = plot(
    plot_time, 
    plot_Fd .* 1e6,  # Convert to μN
    linewidth=1.25,
    color=:green4,
    xlabel=L"Time $t$ (s)",
    ylabel=L"Damping Force $F_d$ ($\mu$N)",
    title="",
    framestyle=:box,
    legend=false,
    dpi=600,
    size=(800, 500),
    left_margin=10mm,
    bottom_margin=10mm,
    top_margin=5mm,
    right_margin=5mm
)
savefig(p_fd, "damping_force_opt.pdf")
savefig(p_fd, "damping_force_opt.png")

# 4. Electrostatic Force (Fe)
p_fe = plot(
    plot_time, 
    plot_Fe .* 1e6,  # Convert to μN
    linewidth=1.25,
    color=:purple3,
    xlabel=L"Time $t$ (s)",
    ylabel=L"Electrostatic Force $F_e$ ($\mu$N)",
    title="",
    framestyle=:box,
    legend=false,
    dpi=600,
    size=(800, 500),
    left_margin=10mm,
    bottom_margin=10mm,
    top_margin=5mm,
    right_margin=5mm
)
savefig(p_fe, "electrostatic_force_opt.pdf")
savefig(p_fe, "electrostatic_force_opt.png")

# 5. External Force (Fext)
p_fext = plot(
    plot_time, 
    plot_Fext .* 1e6,  # Convert to μN
    linewidth=1.25,
    color=:darkorange2,
    xlabel=L"Time $t$ (s)",
    ylabel=L"External Force $F_{ext}$ ($\mu$N)",
    title="",
    framestyle=:box,
    legend=false,
    dpi=600,
    size=(800, 500),
    left_margin=10mm,
    bottom_margin=10mm,
    top_margin=5mm,
    right_margin=5mm
)
savefig(p_fext, "external_force_opt.pdf")
savefig(p_fext, "external_force_opt.png")

# 6. Collision State
p_cs = plot(
    plot_time, 
    plot_collision_numeric,
    linewidth=1.25,
    color=:black,
    xlabel=L"Time $t$ (s)",
    ylabel="Collision State",
    yticks=([0, 1], ["Translational", "Rotational"]),
    title="",
    framestyle=:box,
    legend=false,
    dpi=600,
    size=(800, 500),
    left_margin=10mm,
    bottom_margin=10mm,
    top_margin=5mm,
    right_margin=5mm
)
savefig(p_cs, "collision_state_opt.pdf")
savefig(p_cs, "collision_state_opt.png")

# Create a combined state variables plot (2x3)
p_states = plot(
    layout=(2,3),
    size=(1200, 800),
    dpi=600,
    left_margin=10mm,
    bottom_margin=10mm,
    top_margin=5mm,
    right_margin=5mm,
    plot_title="MEMS Energy Harvester: State Variables",
    plot_titlelocation=:center,
    plot_titlefontsize=14
)

# Add each state variable to the combined plot
plot!(p_states[1], plot_time, plot_x1 .* 1e6, linewidth=2, color=:royalblue3, 
      xlabel="", ylabel=L"$x_1$ ($\mu$m)", title="Shuttle Displacement", legend=false)
      
plot!(p_states[2], plot_time, plot_x1dot .* 1e3, linewidth=2, color=:firebrick3, 
      xlabel="", ylabel=L"$\dot{x}_1$ (mm/s)", title="Shuttle Velocity", legend=false)
      
plot!(p_states[3], plot_time, plot_x2 .* 1e6, linewidth=2, color=:green4, 
      xlabel="", ylabel=L"$x_2$ ($\mu$m)", title="Electrode Displacement", legend=false)
      
plot!(p_states[4], plot_time, plot_x2dot .* 1e3, linewidth=2, color=:purple3, 
      xlabel=L"Time $t$ (s)", ylabel=L"$\dot{x}_2$ (mm/s)", title="Electrode Velocity", legend=false)
      
plot!(p_states[5], plot_time, plot_V .* 1e3, linewidth=2, color=:darkorange2, 
      xlabel=L"Time $t$ (s)", ylabel=L"$V_{out}$ (mV)", title="Output Voltage", legend=false)
      
plot!(p_states[6], plot_time, plot_rel_disp .* 1e6, linewidth=2, color=:teal, 
      xlabel=L"Time $t$ (s)", ylabel=L"$x_1 - x_2$ ($\mu$m)", title="Relative Displacement", legend=false)

savefig(p_states, "combined_states_opt.pdf")
savefig(p_states, "combined_states_opt.png")

# Create a combined forces plot (2x3)
p_forces = plot(
    layout=(2,3),
    size=(1200, 800),
    dpi=600,
    left_margin=10mm,
    bottom_margin=10mm,
    top_margin=5mm,
    right_margin=5mm,
    plot_title="MEMS Energy Harvester: Forces",
    plot_titlelocation=:center,
    plot_titlefontsize=14
)

# Add each force to the combined plot
plot!(p_forces[1], plot_time, plot_Fs .* 1e6, linewidth=2, color=:royalblue3, 
      xlabel="", ylabel=L"$F_s$ ($\mu$N)", title="Spring Force", legend=false)
      
plot!(p_forces[2], plot_time, plot_Fc .* 1e6, linewidth=2, color=:firebrick3, 
      xlabel="", ylabel=L"$F_c$ ($\mu$N)", title="Collision Force", legend=false)
      
plot!(p_forces[3], plot_time, plot_Fd .* 1e6, linewidth=2, color=:green4, 
      xlabel="", ylabel=L"$F_d$ ($\mu$N)", title="Damping Force", legend=false)
      
plot!(p_forces[4], plot_time, plot_Fe .* 1e6, linewidth=2, color=:purple3, 
      xlabel=L"Time $t$ (s)", ylabel=L"$F_e$ ($\mu$N)", title="Electrostatic Force", legend=false)
      
plot!(p_forces[5], plot_time, plot_Fext .* 1e6, linewidth=2, color=:darkorange2, 
      xlabel=L"Time $t$ (s)", ylabel=L"$F_{ext}$ ($\mu$N)", title="External Force", legend=false)
      
plot!(p_forces[6], plot_time, plot_collision_numeric, linewidth=2, color=:black, 
      xlabel=L"Time $t$ (s)", ylabel="State", yticks=([0, 1], ["Trans.", "Rot."]), 
      title="Collision State", legend=false)

savefig(p_forces, "combined_forces_opt.pdf")
savefig(p_forces, "combined_forces_opt.png")

# Create a comprehensive figure for the paper with just 4 key plots
p_paper = plot(
    layout=(2,2),
    size=(1200, 1000),
    dpi=600,
    left_margin=12mm,
    bottom_margin=12mm,
    top_margin=8mm,
    right_margin=8mm,
    plot_title="MEMS Electrostatic Energy Harvester Dynamics",
    plot_titlelocation=:center,
    plot_titlefontsize=16
)

# Add key plots: Displacements, Voltage, Collision Force, and Damping Force
plot!(p_paper[1], plot_time, plot_x1 .* 1e6, linewidth=2.5, color=:royalblue3, 
      xlabel="", ylabel=L"Displacement ($\mu$m)", title="Shuttle and Electrode Displacement", 
      label=L"$x_1$ (Shuttle)")
plot!(p_paper[1], plot_time, plot_x2 .* 1e6, linewidth=2.5, color=:green4, 
      linestyle=:dash, label=L"$x_2$ (Electrode)")
      
plot!(p_paper[2], plot_time, plot_V .* 1e3, linewidth=2.5, color=:darkorange2, 
      xlabel="", ylabel=L"Voltage $V_{out}$ (mV)", title="Output Voltage", legend=false)
      
plot!(p_paper[3], plot_time, plot_Fc .* 1e6, linewidth=2.5, color=:firebrick3, 
      xlabel=L"Time $t$ (s)", ylabel=L"Force $F_c$ ($\mu$N)", title="Collision Force", legend=false)
      
plot!(p_paper[4], plot_time, plot_Fd .* 1e6, linewidth=2.5, color=:green4, 
      xlabel=L"Time $t$ (s)", ylabel=L"Force $F_d$ ($\mu$N)", title="Damping Force", legend=false)

savefig(p_paper, "paper_figure_opt.pdf")
savefig(p_paper, "paper_figure_opt.png")

println("All optimized plots created successfully!")



# --------------------------------------- Alpha Comparison --------------------------------------
# Define high-quality theme for journal publication
function set_journal_theme()
    default(
        fontfamily="Computer Modern",  # LaTeX-like font
        linewidth=2.5,                 # Thicker lines
        foreground_color_legend=nothing, # Transparent legend background
        background_color_legend=nothing, # Transparent legend background
        legendfontsize=10,             # Legend font size
        guidefontsize=12,              # Axis label font size
        tickfontsize=10,               # Tick label font size
        titlefontsize=14,              # Title font size
        size=(800, 600),               # Figure size
        dpi=600,                       # High DPI for print quality
        grid=false,                    # No grid by default
        framestyle=:box,               # Box-style frame
        foreground_color_axis=:black,  # Black axes
        tick_direction=:out,           # Ticks pointing outward
        palette=:default               # Color palette
    )
end

# Apply the theme
set_journal_theme()

# Main function to generate plots for different alpha values
function generate_alpha_plots(alpha_values; num_oscillations=4, start_time=0.2, prefix="")
    # Store results for each alpha value
    all_results = Dict()
    
    # Create a directory for alpha comparison plots if it doesn't exist
    alpha_dir = "alpha_comparison"
    if !isdir(alpha_dir)
        mkdir(alpha_dir)
    end
    
    # Generate and save plots for each alpha value
    for (i, alpha) in enumerate(alpha_values)
        println("Processing α = $alpha ($(i)/$(length(alpha_values)))")
        
        # Run the simulation with this alpha value
        result = run_simulation(alpha)
        
        # Store the result
        all_results[alpha] = result
        
        # Generate plots for this alpha value
        generate_plots_for_alpha(result, alpha, num_oscillations, start_time, prefix)
    end
    
    # Generate comparison plots across all alpha values
    generate_comparison_plots(all_results, alpha_values, num_oscillations, start_time)
    
    return all_results
end

# Function to run a simulation with a given alpha value
function run_simulation(alpha)
    # Set to `true` to use sine wave, `false` for displaced IC
    use_sine = true

    # Set the external force parameters
    f = 20.0  # Frequency (Hz)
    g = 9.81  # Gravitational constant (m/s²)
    A = alpha * g
    t_ramp = 0.2  # Ramp-up duration (s)
    
    # Define Fext_input based on your choice
    if use_sine
        Fext_input = t -> A * (t < t_ramp ? t / t_ramp : 1.0) * sin(2 * π * f * t)
    else
        Fext_input = t -> 0.0
    end
    
    # Initialize parameters
    p_new = deepcopy(AnalyticalModel.p)
    
    # Initial conditions
    x10 = 0.0      # Initial displacement
    x10dot = 0.0   # Initial velocity
    x20 = 0.0      # Initial displacement
    x20dot = 0.0   # Initial velocity
    
    # Compute initial electrostatic parameters
    Ctotal0, Fe0 = AnalyticalModel.electrostatic(
        x10, x20, 0.0, 
        p_new.g0, p_new.gp, p_new.a, 
        p_new.e, p_new.ep, p_new.cp, 
        p_new.wt, p_new.wb, p_new.ke, 
        p_new.E, p_new.I, p_new.Leff, 
        p_new.Tf, p_new.Tp, p_new.N
    )
    
    Q0 = p_new.Vbias * Ctotal0  # Initial charge
    Vout0 = p_new.Vbias - (Q0 / Ctotal0)  # Initial voltage
    z0 = [x10, x10dot, x20, x20dot, Q0, Vout0]
    tspan = (0.0, 0.5)  # simulation length
    abstol = 1e-9  # absolute solver tolerance
    reltol = 1e-6  # relative solver tolerance
    
    # Define the wrapper function for ODE solver
    function CoupledSystem_wrapper!(dz, z, p, t)
        current_acceleration = Fext_input(t)
        AnalyticalModel.CoupledSystem!(dz, z, p, t, current_acceleration)
    end
    
    # Solve the ODE system
    eqn = ODEProblem(CoupledSystem_wrapper!, z0, tspan, p_new)
    sol = solve(eqn, Rosenbrock23(); abstol=abstol, reltol=reltol, maxiters=1e7)
    
    # Extract results
    x1 = [u[1] for u in sol.u]      # Shuttle displacement
    x1dot = [u[2] for u in sol.u]   # Shuttle velocity
    x2 = [u[3] for u in sol.u]      # Electrode displacement
    x2dot = [u[4] for u in sol.u]   # Electrode velocity
    Qvar = [u[5] for u in sol.u]    # Charge
    V = [u[6] for u in sol.u]       # Output voltage
    
    # Initialize arrays to store forces
    Fs_array = Float64[]  # Suspension spring force
    Fc_array = Float64[]  # Collision force
    Fd_array = Float64[]  # Damping force
    Fe_array = Float64[]  # Electrostatic force
    CollisionState_array = String[]  # Collision state
    
    # Calculate forces at each time point
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
        if collision_state == "translational"
            Fd = AnalyticalModel.damping(z1, z2, z3, z4, p_new.a, p_new.c, p_new.gp, p_new.Leff, p_new.Tf, p_new.eta, p_new.lambda, p_new.sigmap)
        else
            Fd = AnalyticalModel.damping(z1, z2, z3, z2, p_new.a, p_new.c, p_new.gp, p_new.Leff, p_new.Tf, p_new.eta, p_new.lambda, p_new.sigmap)
        end
        push!(Fd_array, Fd)
        
        # Compute Fe (Electrostatic force)
        Ctotalx, Fe = AnalyticalModel.electrostatic(
            z1, z3, z5, p_new.g0, p_new.gp, p_new.a, 
            p_new.e, p_new.ep, p_new.cp, p_new.wt, 
            p_new.wb, p_new.ke, p_new.E, p_new.I, 
            p_new.Leff, p_new.Tf, p_new.Tp, p_new.N
        )
        push!(Fe_array, Fe)
    end
    
    # Return the simulation results as a dictionary
    return Dict(
        "time" => sol.t,
        "x1" => x1,
        "x1dot" => x1dot,
        "x2" => x2,
        "x2dot" => x2dot,
        "Qvar" => Qvar,
        "V" => V,
        "Fs" => Fs_array,
        "Fc" => Fc_array,
        "Fd" => Fd_array,
        "Fe" => Fe_array,
        "CollisionState" => CollisionState_array,
        "Fext" => [Fext_input(t) for t in sol.t],
        "alpha" => alpha,
        "frequency" => f
    )
end

# Function to generate plots for a single alpha value
function generate_plots_for_alpha(result, alpha, num_oscillations, start_time, prefix)
    # Calculate the approximate end time based on the frequency and number of oscillations
    f = result["frequency"]
    period = 1.0 / f
    end_time = start_time + num_oscillations * period
    
    # Find indices for the specified time window
    time = result["time"]
    start_idx = findfirst(t -> t >= start_time, time)
    end_idx = findfirst(t -> t >= end_time, time)
    
    if start_idx === nothing
        start_idx = 1
    end
    if end_idx === nothing
        end_idx = length(time)
    end
    
    # Extract data for the specified time window
    plot_time = time[start_idx:end_idx]
    plot_x1 = result["x1"][start_idx:end_idx]
    plot_x1dot = result["x1dot"][start_idx:end_idx]
    plot_x2 = result["x2"][start_idx:end_idx]
    plot_x2dot = result["x2dot"][start_idx:end_idx]
    plot_Qvar = result["Qvar"][start_idx:end_idx]
    plot_V = result["V"][start_idx:end_idx]
    plot_Fs = result["Fs"][start_idx:end_idx]
    plot_Fc = result["Fc"][start_idx:end_idx]
    plot_Fd = result["Fd"][start_idx:end_idx]
    plot_Fe = result["Fe"][start_idx:end_idx]
    plot_Fext = result["Fext"][start_idx:end_idx]
    plot_collision_state = result["CollisionState"][start_idx:end_idx]
    plot_collision_numeric = [s == "translational" ? 0 : 1 for s in plot_collision_state]
    
    # Directory for this alpha value
    alpha_dir = "alpha_$(alpha)"
    if !isdir(alpha_dir)
        mkdir(alpha_dir)
    end
    
    # Create a combined state variables plot (2x3) with updated layout
    p_states = plot(
        layout=(2,3),
        size=(1200, 800),
        dpi=600,
        left_margin=10mm,
        bottom_margin=10mm,
        top_margin=5mm,
        right_margin=5mm,
        plot_title="MEMS Energy Harvester: State Variables (α = $alpha)",
        plot_titlelocation=:center,
        plot_titlefontsize=14
    )

    # Add each state variable to the combined plot with the new layout
    plot!(p_states[1], plot_time, plot_x1 .* 1e6, linewidth=2, color=:royalblue3, 
          xlabel="", ylabel=L"$x_1$ ($\mu$m)", title="Shuttle Displacement", legend=false)
          
    plot!(p_states[2], plot_time, plot_x1dot .* 1e3, linewidth=2, color=:firebrick3, 
          xlabel="", ylabel=L"$\dot{x}_1$ (mm/s)", title="Shuttle Velocity", legend=false)
          
    plot!(p_states[3], plot_time, plot_x2 .* 1e6, linewidth=2, color=:green4, 
          xlabel="", ylabel=L"$x_2$ ($\mu$m)", title="Electrode Displacement", legend=false)
          
    plot!(p_states[4], plot_time, plot_x2dot .* 1e3, linewidth=2, color=:purple3, 
          xlabel=L"Time $t$ (s)", ylabel=L"$\dot{x}_2$ (mm/s)", title="Electrode Velocity", legend=false)
          
    plot!(p_states[5], plot_time, plot_Qvar .* 1e12, linewidth=2, color=:darkorange, 
          xlabel=L"Time $t$ (s)", ylabel=L"$Q_{var}$ (pC)", title="Charge", legend=false)
          
    plot!(p_states[6], plot_time, plot_V .* 1e3, linewidth=2, color=:darkorange2, 
          xlabel=L"Time $t$ (s)", ylabel=L"$V_{out}$ (mV)", title="Output Voltage", legend=false)

    savefig(p_states, joinpath(alpha_dir, "$(prefix)combined_states_alpha_$(alpha).pdf"))
    savefig(p_states, joinpath(alpha_dir, "$(prefix)combined_states_alpha_$(alpha).png"))

    # Create a combined forces plot (2x3)
    p_forces = plot(
        layout=(2,3),
        size=(1200, 800),
        dpi=600,
        left_margin=10mm,
        bottom_margin=10mm,
        top_margin=5mm,
        right_margin=5mm,
        plot_title="MEMS Energy Harvester: Forces (α = $alpha)",
        plot_titlelocation=:center,
        plot_titlefontsize=14
    )

    # Add each force to the combined plot
    plot!(p_forces[1], plot_time, plot_Fs .* 1e6, linewidth=2, color=:royalblue3, 
          xlabel="", ylabel=L"$F_s$ ($\mu$N)", title="Spring Force", legend=false)
          
    plot!(p_forces[2], plot_time, plot_Fc .* 1e6, linewidth=2, color=:firebrick3, 
          xlabel="", ylabel=L"$F_c$ ($\mu$N)", title="Collision Force", legend=false)
          
    plot!(p_forces[3], plot_time, plot_Fd .* 1e6, linewidth=2, color=:green4, 
          xlabel="", ylabel=L"$F_d$ ($\mu$N)", title="Damping Force", legend=false)
          
    plot!(p_forces[4], plot_time, plot_Fe .* 1e6, linewidth=2, color=:purple3, 
          xlabel=L"Time $t$ (s)", ylabel=L"$F_e$ ($\mu$N)", title="Electrostatic Force", legend=false)
          
    plot!(p_forces[5], plot_time, plot_Fext .* 1e6, linewidth=2, color=:darkorange2, 
          xlabel=L"Time $t$ (s)", ylabel=L"$F_{ext}$ ($\mu$N)", title="External Force", legend=false)
          
    plot!(p_forces[6], plot_time, plot_collision_numeric, linewidth=2, color=:black, 
          xlabel=L"Time $t$ (s)", ylabel="State", yticks=([0, 1], ["Trans.", "Rot."]), 
          title="Collision State", legend=false)

    savefig(p_forces, joinpath(alpha_dir, "$(prefix)combined_forces_alpha_$(alpha).pdf"))
    savefig(p_forces, joinpath(alpha_dir, "$(prefix)combined_forces_alpha_$(alpha).png"))

    # Create a comprehensive figure for the paper with just 4 key plots
    p_paper = plot(
        layout=(2,2),
        size=(1200, 1000),
        dpi=600,
        left_margin=12mm,
        bottom_margin=12mm,
        top_margin=8mm,
        right_margin=8mm,
        plot_title="MEMS Electrostatic Energy Harvester Dynamics (α = $alpha)",
        plot_titlelocation=:center,
        plot_titlefontsize=16
    )

    # Add key plots: Displacements, Voltage, Collision Force, and Damping Force
    plot!(p_paper[1], plot_time, plot_x1 .* 1e6, linewidth=2.5, color=:royalblue3, 
          xlabel="", ylabel=L"Displacement ($\mu$m)", title="Shuttle and Electrode Displacement", 
          label=L"$x_1$ (Shuttle)")
    plot!(p_paper[1], plot_time, plot_x2 .* 1e6, linewidth=2.5, color=:green4, 
          linestyle=:dash, label=L"$x_2$ (Electrode)")
          
    plot!(p_paper[2], plot_time, plot_V .* 1e3, linewidth=2.5, color=:darkorange2, 
          xlabel="", ylabel=L"Voltage $V_{out}$ (mV)", title="Output Voltage", legend=false)
          
    plot!(p_paper[3], plot_time, plot_Fc .* 1e6, linewidth=2.5, color=:firebrick3, 
          xlabel=L"Time $t$ (s)", ylabel=L"Force $F_c$ ($\mu$N)", title="Collision Force", legend=false)
          
    plot!(p_paper[4], plot_time, plot_Fd .* 1e6, linewidth=2.5, color=:green4, 
          xlabel=L"Time $t$ (s)", ylabel=L"Force $F_d$ ($\mu$N)", title="Damping Force", legend=false)

    savefig(p_paper, joinpath(alpha_dir, "$(prefix)paper_figure_alpha_$(alpha).pdf"))
    savefig(p_paper, joinpath(alpha_dir, "$(prefix)paper_figure_alpha_$(alpha).png"))
    
    return Dict(
        "states_plot" => p_states,
        "forces_plot" => p_forces,
        "paper_plot" => p_paper
    )
end

# Function to generate comparison plots across all alpha values
function generate_comparison_plots(all_results, alpha_values, num_oscillations, start_time)
    # Create a directory for alpha comparison if it doesn't exist
    alpha_dir = "alpha_comparison"
    if !isdir(alpha_dir)
        mkdir(alpha_dir)
    end
    
    # Calculate the time range based on the first alpha value
    f = all_results[alpha_values[1]]["frequency"]
    period = 1.0 / f
    end_time = start_time + num_oscillations * period
    
    # For each alpha value, extract data for the specified time window
    trimmed_data = Dict()
    
    for alpha in alpha_values
        result = all_results[alpha]
        time = result["time"]
        
        # Find indices for the specified time window
        start_idx = findfirst(t -> t >= start_time, time)
        end_idx = findfirst(t -> t >= end_time, time)
        
        if start_idx === nothing
            start_idx = 1
        end
        if end_idx === nothing
            end_idx = length(time)
        end
        
        # Extract and store the trimmed data
        trimmed_data[alpha] = Dict(
            "time" => time[start_idx:end_idx],
            "x1" => result["x1"][start_idx:end_idx],
            "x1dot" => result["x1dot"][start_idx:end_idx],
            "x2" => result["x2"][start_idx:end_idx],
            "x2dot" => result["x2dot"][start_idx:end_idx],
            "Qvar" => result["Qvar"][start_idx:end_idx],
            "V" => result["V"][start_idx:end_idx],
            "Fs" => result["Fs"][start_idx:end_idx],
            "Fc" => result["Fc"][start_idx:end_idx],
            "Fd" => result["Fd"][start_idx:end_idx],
            "Fe" => result["Fe"][start_idx:end_idx],
            "Fext" => result["Fext"][start_idx:end_idx],
            "CollisionState" => result["CollisionState"][start_idx:end_idx]
        )
    end
    
    # Create comparison plots for key variables
    variables = ["x1", "x1dot", "x2", "x2dot", "V", "Qvar", "Fc", "Fd"]
    titles = ["Shuttle Displacement", "Shuttle Velocity", "Electrode Displacement", 
              "Electrode Velocity", "Output Voltage", "Charge", "Collision Force", "Damping Force"]
    units = [L"$\mu$m", L"mm/s", L"$\mu$m", L"mm/s", L"mV", L"pC", L"$\mu$N", L"$\mu$N"]
    scale_factors = [1e6, 1e3, 1e6, 1e3, 1e3, 1e12, 1e6, 1e6]  # Convert to appropriate units
    
    for (var, title, unit, scale) in zip(variables, titles, units, scale_factors)
        p = plot(
            size=(800, 600),
            dpi=600,
            left_margin=10mm,
            bottom_margin=10mm,
            top_margin=5mm,
            right_margin=5mm,
            title="Effect of Acceleration on $title",
            titlelocation=:center,
            titlefontsize=14,
            xlabel=L"Time $t$ (s)",
            ylabel="$title $unit",
            legend=:topright
        )
        
        # Add a line for each alpha value
        for alpha in alpha_values
            data = trimmed_data[alpha]
            plot!(p, data["time"], data[var] .* scale, 
                  linewidth=2, 
                  label=L"\alpha = %$alpha")
        end
        
        savefig(p, joinpath(alpha_dir, "comparison_$(var).pdf"))
        savefig(p, joinpath(alpha_dir, "comparison_$(var).png"))
    end
    
    # Create RMS voltage vs. alpha plot
    rms_voltages = Float64[]
    for alpha in alpha_values
        data = trimmed_data[alpha]
        rms_voltage = sqrt(mean(data["V"].^2)) * 1e3  # Convert to mV
        push!(rms_voltages, rms_voltage)
    end
    
    p_rms = plot(
        alpha_values, 
        rms_voltages,
        linewidth=1.15,
        marker=:circle,
        markersize=6,
        color=:darkorange2,
        xlabel=L"Acceleration Parameter $\alpha$",
        ylabel=L"RMS Voltage $V_{rms}$ (mV)",
        title="RMS Voltage vs. Acceleration",
        legend=false,
        dpi=600,
        size=(800, 600),
        left_margin=10mm,
        bottom_margin=10mm,
        top_margin=5mm,
        right_margin=5mm
    )
    
    savefig(p_rms, joinpath(alpha_dir, "rms_voltage_vs_alpha.pdf"))
    savefig(p_rms, joinpath(alpha_dir, "rms_voltage_vs_alpha.png"))
    
    # Create max displacement vs. alpha plot
    max_displacements = Float64[]
    for alpha in alpha_values
        data = trimmed_data[alpha]
        max_displacement = maximum(abs.(data["x1"])) * 1e6  # Convert to μm
        push!(max_displacements, max_displacement)
    end
    
    p_disp = plot(
        alpha_values, 
        max_displacements,
        linewidth=1.15,
        marker=:circle,
        markersize=6,
        color=:royalblue3,
        xlabel=L"Acceleration Parameter $\alpha$",
        ylabel=L"Maximum Displacement $x_{1,max}$ ($\mu$m)",
        title="Maximum Displacement vs. Acceleration",
        legend=false,
        dpi=600,
        size=(800, 600),
        left_margin=10mm,
        bottom_margin=10mm,
        top_margin=5mm,
        right_margin=5mm
    )
    
    savefig(p_disp, joinpath(alpha_dir, "max_displacement_vs_alpha.pdf"))
    savefig(p_disp, joinpath(alpha_dir, "max_displacement_vs_alpha.png"))
    
    # Create summary plot with both RMS voltage and max displacement
    p_summary = plot(
        layout=(2,1),
        size=(800, 800),
        dpi=600,
        left_margin=10mm,
        bottom_margin=10mm,
        top_margin=5mm,
        right_margin=5mm,
        plot_title="Effect of Acceleration on MEMS Energy Harvester Performance",
        plot_titlelocation=:center,
        plot_titlefontsize=14
    )
    
    plot!(p_summary[1], alpha_values, rms_voltages, linewidth=2.5, marker=:circle, markersize=6, 
          color=:darkorange2, xlabel="", ylabel=L"RMS Voltage $V_{rms}$ (mV)", 
          title="RMS Voltage vs. Acceleration", legend=false)
    
    plot!(p_summary[2], alpha_values, max_displacements, linewidth=2.5, marker=:circle, markersize=6, 
          color=:royalblue3, xlabel=L"Acceleration Parameter $\alpha$", 
          ylabel=L"Maximum Displacement $x_{1,max}$ ($\mu$m)", 
          title="Maximum Displacement vs. Acceleration", legend=false)
    
    savefig(p_summary, joinpath(alpha_dir, "performance_summary.pdf"))
    savefig(p_summary, joinpath(alpha_dir, "performance_summary.png"))
end

# Example usage:
# Define alpha values for the parameter sweep
alpha_values = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5]

# Generate plots for all alpha values
all_results = generate_alpha_plots(alpha_values, num_oscillations=4, start_time=0.2, prefix="")

println("All alpha parameter plots created successfully!")



# --------------------------------------- Optimized Plots for Journal Publication --------------------------------------
# Automated Journal-Quality Plotting System with Flexible Time Window Selection
using Plots, LaTeXStrings, Measures, Statistics

# Define high-quality theme for journal publication
function set_journal_theme()
    default(
        fontfamily="Computer Modern",  # LaTeX-like font
        linewidth=1.15,                 # Thicker lines
        foreground_color_legend=nothing, # Transparent legend background
        background_color_legend=nothing, # Transparent legend background
        legendfontsize=10,             # Legend font size
        guidefontsize=12,              # Axis label font size
        tickfontsize=10,               # Tick label font size
        titlefontsize=14,              # Title font size
        size=(800, 600),               # Figure size
        dpi=600,                       # High DPI for print quality
        grid=false,                    # No grid by default
        framestyle=:box,               # Box-style frame
        foreground_color_axis=:black,  # Black axes
        tick_direction=:out,           # Ticks pointing outward
        palette=:default               # Color palette
    )
end

# Apply the theme
set_journal_theme()

# Main function to generate plots for different alpha values
function generate_alpha_plots(alpha_values; 
                              time_window=:auto,   # Can be :auto or a tuple (start_time, end_time)
                              num_oscillations=4,  # Used if time_window is :auto
                              start_time=0.2,      # Default start time if time_window is :auto
                              prefix="")
    # Store results for each alpha value
    all_results = Dict()
    
    # Create a directory for alpha comparison plots if it doesn't exist
    alpha_dir = "alpha_comparison"
    if !isdir(alpha_dir)
        mkdir(alpha_dir)
    end
    
    # Generate and save plots for each alpha value
    for (i, alpha) in enumerate(alpha_values)
        println("Processing α = $alpha ($(i)/$(length(alpha_values)))")
        
        # Run the simulation with this alpha value
        result = run_simulation(alpha)
        
        # Store the result
        all_results[alpha] = result
        
        # Generate plots for this alpha value
        generate_plots_for_alpha(result, alpha, time_window, num_oscillations, start_time, prefix)
    end
    
    # Only generate comparison plots if we have more than one alpha value
    if length(alpha_values) > 1
        # Generate comparison plots across all alpha values
        generate_comparison_plots(all_results, alpha_values, time_window, num_oscillations, start_time)
        
        # Generate combined grid plot for all alphas
        generate_combined_grid_plot(all_results, alpha_values, time_window, num_oscillations, start_time)
    else
        println("Only one alpha value provided. Skipping comparison plots.")
    end
    
    return all_results
end

# Function to run a simulation with a given alpha value
function run_simulation(alpha)
    # Set to `true` to use sine wave, `false` for displaced IC
    use_sine = true

    # Set the external force parameters
    f = 20.0  # Frequency (Hz)
    g = 9.81  # Gravitational constant (m/s²)
    A = alpha * g
    t_ramp = 0.2  # Ramp-up duration (s)
    
    # Define Fext_input based on your choice
    if use_sine
        Fext_input = t -> A * (t < t_ramp ? t / t_ramp : 1.0) * sin(2 * π * f * t)
    else
        Fext_input = t -> 0.0
    end
    
    # Initialize parameters
    p_new = deepcopy(AnalyticalModel.p)
    
    # Initial conditions
    x10 = 0.0      # Initial displacement
    x10dot = 0.0   # Initial velocity
    x20 = 0.0      # Initial displacement
    x20dot = 0.0   # Initial velocity
    
    # Compute initial electrostatic parameters
    Ctotal0, Fe0 = AnalyticalModel.electrostatic(
        x10, x20, 0.0, 
        p_new.g0, p_new.gp, p_new.a, 
        p_new.e, p_new.ep, p_new.cp, 
        p_new.wt, p_new.wb, p_new.ke, 
        p_new.E, p_new.I, p_new.Leff, 
        p_new.Tf, p_new.Tp, p_new.N
    )
    
    Q0 = p_new.Vbias * Ctotal0  # Initial charge
    Vout0 = p_new.Vbias - (Q0 / Ctotal0)  # Initial voltage
    z0 = [x10, x10dot, x20, x20dot, Q0, Vout0]
    tspan = (0.0, 0.5)  # simulation length
    abstol = 1e-9  # absolute solver tolerance
    reltol = 1e-6  # relative solver tolerance
    
    # Define the wrapper function for ODE solver
    function CoupledSystem_wrapper!(dz, z, p, t)
        current_acceleration = Fext_input(t)
        AnalyticalModel.CoupledSystem!(dz, z, p, t, current_acceleration)
    end
    
    # Solve the ODE system
    eqn = ODEProblem(CoupledSystem_wrapper!, z0, tspan, p_new)
    sol = solve(eqn, Rosenbrock23(); abstol=abstol, reltol=reltol, maxiters=1e7)
    
    # Extract results
    x1 = [u[1] for u in sol.u]      # Shuttle displacement
    x1dot = [u[2] for u in sol.u]   # Shuttle velocity
    x2 = [u[3] for u in sol.u]      # Electrode displacement
    x2dot = [u[4] for u in sol.u]   # Electrode velocity
    Qvar = [u[5] for u in sol.u]    # Charge
    V = [u[6] for u in sol.u]       # Output voltage
    
    # Initialize arrays to store forces
    Fs_array = Float64[]  # Suspension spring force
    Fc_array = Float64[]  # Collision force
    Fd_array = Float64[]  # Damping force
    Fe_array = Float64[]  # Electrostatic force
    CollisionState_array = String[]  # Collision state
    
    # Calculate forces at each time point
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
        if collision_state == "translational"
            Fd = AnalyticalModel.damping(z1, z2, z3, z4, p_new.a, p_new.c, p_new.gp, p_new.Leff, p_new.Tf, p_new.eta, p_new.lambda, p_new.sigmap)
        else
            Fd = AnalyticalModel.damping(z1, z2, z3, z2, p_new.a, p_new.c, p_new.gp, p_new.Leff, p_new.Tf, p_new.eta, p_new.lambda, p_new.sigmap)
        end
        push!(Fd_array, Fd)
        
        # Compute Fe (Electrostatic force)
        Ctotalx, Fe = AnalyticalModel.electrostatic(
            z1, z3, z5, p_new.g0, p_new.gp, p_new.a, 
            p_new.e, p_new.ep, p_new.cp, p_new.wt, 
            p_new.wb, p_new.ke, p_new.E, p_new.I, 
            p_new.Leff, p_new.Tf, p_new.Tp, p_new.N
        )
        push!(Fe_array, Fe)
    end
    
    # Return the simulation results as a dictionary
    return Dict(
        "time" => sol.t,
        "x1" => x1,
        "x1dot" => x1dot,
        "x2" => x2,
        "x2dot" => x2dot,
        "Qvar" => Qvar,
        "V" => V,
        "Fs" => Fs_array,
        "Fc" => Fc_array,
        "Fd" => Fd_array,
        "Fe" => Fe_array,
        "CollisionState" => CollisionState_array,
        "Fext" => [Fext_input(t) for t in sol.t],
        "alpha" => alpha,
        "frequency" => f
    )
end

# Function to determine time window based on user preferences
function determine_time_window(result, time_window, num_oscillations, start_time)
    if time_window == :auto
        # Calculate time window based on frequency and desired number of oscillations
        f = result["frequency"]
        period = 1.0 / f
        end_time = start_time + num_oscillations * period
        println("Using frequency-based time window: $(start_time) to $(round(end_time, digits=4)) seconds ($(num_oscillations) oscillations)")
    else
        # Use explicit time window
        start_time, end_time = time_window
        f = result["frequency"]
        period = 1.0 / f
        num_oscillations = (end_time - start_time) / period
        println("Using explicit time window: $(start_time) to $(end_time) seconds (approximately $(round(num_oscillations, digits=2)) oscillations)")
    end
    
    return start_time, end_time
end

# Function to generate plots for a single alpha value
function generate_plots_for_alpha(result, alpha, time_window, num_oscillations, start_time, prefix)
    # Determine the time window
    start_time, end_time = determine_time_window(result, time_window, num_oscillations, start_time)
    
    # Find indices for the specified time window
    time = result["time"]
    start_idx = findfirst(t -> t >= start_time, time)
    end_idx = findfirst(t -> t >= end_time, time)
    
    if start_idx === nothing
        start_idx = 1
    end
    if end_idx === nothing
        end_idx = length(time)
    end
    
    # Extract data for the specified time window
    plot_time = time[start_idx:end_idx]
    plot_x1 = result["x1"][start_idx:end_idx]
    plot_x1dot = result["x1dot"][start_idx:end_idx]
    plot_x2 = result["x2"][start_idx:end_idx]
    plot_x2dot = result["x2dot"][start_idx:end_idx]
    plot_Qvar = result["Qvar"][start_idx:end_idx]
    plot_V = result["V"][start_idx:end_idx]
    plot_Fs = result["Fs"][start_idx:end_idx]
    plot_Fc = result["Fc"][start_idx:end_idx]
    plot_Fd = result["Fd"][start_idx:end_idx]
    plot_Fe = result["Fe"][start_idx:end_idx]
    plot_Fext = result["Fext"][start_idx:end_idx]
    plot_collision_state = result["CollisionState"][start_idx:end_idx]
    plot_collision_numeric = [s == "translational" ? 0 : 1 for s in plot_collision_state]
    
    # Directory for this alpha value
    alpha_dir = "alpha_$(alpha)"
    if !isdir(alpha_dir)
        mkdir(alpha_dir)
    end
    
    # Create a combined state variables plot (2x3) with updated layout
    p_states = plot(
        layout=(2,3),
        size=(1200, 800),
        dpi=600,
        left_margin=10mm,
        bottom_margin=10mm,
        top_margin=5mm,
        right_margin=5mm,
        plot_title="MEMS Energy Harvester: State Variables (α = $alpha)",
        plot_titlelocation=:center,
        plot_titlefontsize=14
    )

    # Add each state variable to the combined plot with the new layout
    plot!(p_states[1], plot_time, plot_x1 .* 1e6, linewidth=2, color=:royalblue3, 
          xlabel="", ylabel=L"$x_1$ ($\mu$m)", title="Shuttle Displacement", legend=false)
          
    plot!(p_states[2], plot_time, plot_x1dot .* 1e3, linewidth=2, color=:firebrick3, 
          xlabel="", ylabel=L"$\dot{x}_1$ (mm/s)", title="Shuttle Velocity", legend=false)
          
    plot!(p_states[3], plot_time, plot_x2 .* 1e6, linewidth=2, color=:green4, 
          xlabel="", ylabel=L"$x_2$ ($\mu$m)", title="Electrode Displacement", legend=false)
          
    plot!(p_states[4], plot_time, plot_x2dot .* 1e3, linewidth=2, color=:purple3, 
          xlabel=L"Time $t$ (s)", ylabel=L"$\dot{x}_2$ (mm/s)", title="Electrode Velocity", legend=false)
          
    plot!(p_states[5], plot_time, plot_Qvar .* 1e12, linewidth=2, color=:darkorange, 
          xlabel=L"Time $t$ (s)", ylabel=L"$Q_{var}$ (pC)", title="Charge", legend=false)
          
    plot!(p_states[6], plot_time, plot_V .* 1e3, linewidth=2, color=:darkorange2, 
          xlabel=L"Time $t$ (s)", ylabel=L"$V_{out}$ (mV)", title="Output Voltage", legend=false)

    savefig(p_states, joinpath(alpha_dir, "$(prefix)combined_states_alpha_$(alpha).pdf"))
    savefig(p_states, joinpath(alpha_dir, "$(prefix)combined_states_alpha_$(alpha).png"))

    # Create a combined forces plot (2x3)
    p_forces = plot(
        layout=(2,3),
        size=(1200, 800),
        dpi=600,
        left_margin=10mm,
        bottom_margin=10mm,
        top_margin=5mm,
        right_margin=5mm,
        plot_title="MEMS Energy Harvester: Forces (α = $alpha)",
        plot_titlelocation=:center,
        plot_titlefontsize=14
    )

    # Add each force to the combined plot
    plot!(p_forces[1], plot_time, plot_Fs .* 1e6, linewidth=2, color=:royalblue3, 
          xlabel="", ylabel=L"$F_s$ ($\mu$N)", title="Spring Force", legend=false)
          
    plot!(p_forces[2], plot_time, plot_Fc .* 1e6, linewidth=2, color=:firebrick3, 
          xlabel="", ylabel=L"$F_c$ ($\mu$N)", title="Collision Force", legend=false)
          
    plot!(p_forces[3], plot_time, plot_Fd .* 1e6, linewidth=2, color=:green4, 
          xlabel="", ylabel=L"$F_d$ ($\mu$N)", title="Damping Force", legend=false)
          
    plot!(p_forces[4], plot_time, plot_Fe .* 1e6, linewidth=2, color=:purple3, 
          xlabel=L"Time $t$ (s)", ylabel=L"$F_e$ ($\mu$N)", title="Electrostatic Force", legend=false)
          
    plot!(p_forces[5], plot_time, plot_Fext .* 1e6, linewidth=2, color=:darkorange2, 
          xlabel=L"Time $t$ (s)", ylabel=L"$F_{ext}$ ($\mu$N)", title="External Force", legend=false)
          
    plot!(p_forces[6], plot_time, plot_collision_numeric, linewidth=2, color=:black, 
          xlabel=L"Time $t$ (s)", ylabel="State", yticks=([0, 1], ["Trans.", "Rot."]), 
          title="Collision State", legend=false)

    savefig(p_forces, joinpath(alpha_dir, "$(prefix)combined_forces_alpha_$(alpha).pdf"))
    savefig(p_forces, joinpath(alpha_dir, "$(prefix)combined_forces_alpha_$(alpha).png"))

    # Create a comprehensive figure for the paper with just 4 key plots
    p_paper = plot(
        layout=(2,2),
        size=(1200, 1000),
        dpi=600,
        left_margin=12mm,
        bottom_margin=12mm,
        top_margin=8mm,
        right_margin=8mm,
        plot_title="MEMS Electrostatic Energy Harvester Dynamics (α = $alpha)",
        plot_titlelocation=:center,
        plot_titlefontsize=16
    )

    # Add key plots: Displacements, Voltage, Collision Force, and Damping Force
    plot!(p_paper[1], plot_time, plot_x1 .* 1e6, linewidth=2.5, color=:royalblue3, 
          xlabel="", ylabel=L"Displacement ($\mu$m)", title="Shuttle and Electrode Displacement", 
          label=L"$x_1$ (Shuttle)")
    plot!(p_paper[1], plot_time, plot_x2 .* 1e6, linewidth=2.5, color=:green4, 
          linestyle=:dash, label=L"$x_2$ (Electrode)")
          
    plot!(p_paper[2], plot_time, plot_V .* 1e3, linewidth=2.5, color=:darkorange2, 
          xlabel="", ylabel=L"Voltage $V_{out}$ (mV)", title="Output Voltage", legend=false)
          
    plot!(p_paper[3], plot_time, plot_Fc .* 1e6, linewidth=2.5, color=:firebrick3, 
          xlabel=L"Time $t$ (s)", ylabel=L"Force $F_c$ ($\mu$N)", title="Collision Force", legend=false)
          
    plot!(p_paper[4], plot_time, plot_Fd .* 1e6, linewidth=2.5, color=:green4, 
          xlabel=L"Time $t$ (s)", ylabel=L"Force $F_d$ ($\mu$N)", title="Damping Force", legend=false)

    savefig(p_paper, joinpath(alpha_dir, "$(prefix)paper_figure_alpha_$(alpha).pdf"))
    savefig(p_paper, joinpath(alpha_dir, "$(prefix)paper_figure_alpha_$(alpha).png"))
    
    return Dict(
        "states_plot" => p_states,
        "forces_plot" => p_forces,
        "paper_plot" => p_paper,
        "plot_data" => Dict(
            "time" => plot_time,
            "x1" => plot_x1,
            "x1dot" => plot_x1dot,
            "x2" => plot_x2,
            "x2dot" => plot_x2dot,
            "Qvar" => plot_Qvar,
            "V" => plot_V,
            "Fs" => plot_Fs,
            "Fc" => plot_Fc,
            "Fd" => plot_Fd,
            "Fe" => plot_Fe,
            "Fext" => plot_Fext,
            "CollisionState" => plot_collision_numeric
        )
    )
end

# Function to generate comparison plots across all alpha values
function generate_comparison_plots(all_results, alpha_values, time_window, num_oscillations, start_time)
    # Check if we have enough alpha values for comparisons
    if length(alpha_values) < 2
        println("Need at least 2 alpha values for comparison plots. Skipping.")
        return Dict()
    end
    
    # Create a directory for alpha comparison if it doesn't exist
    alpha_dir = "alpha_comparison"
    if !isdir(alpha_dir)
        mkdir(alpha_dir)
    end
    
    # Store trimmed data for all alpha values
    trimmed_data = Dict()
    
    # Process each alpha value to get the time window and extract data
    for alpha in alpha_values
        result = all_results[alpha]
        
        # Determine the time window
        start_t, end_t = determine_time_window(result, time_window, num_oscillations, start_time)
        
        # Find indices for the specified time window
        time = result["time"]
        start_idx = findfirst(t -> t >= start_t, time)
        end_idx = findfirst(t -> t >= end_t, time)
        
        if start_idx === nothing
            start_idx = 1
        end
        if end_idx === nothing
            end_idx = length(time)
        end
        
        # Extract and store the trimmed data
        trimmed_data[alpha] = Dict(
            "time" => time[start_idx:end_idx],
            "x1" => result["x1"][start_idx:end_idx],
            "x1dot" => result["x1dot"][start_idx:end_idx],
            "x2" => result["x2"][start_idx:end_idx],
            "x2dot" => result["x2dot"][start_idx:end_idx],
            "Qvar" => result["Qvar"][start_idx:end_idx],
            "V" => result["V"][start_idx:end_idx],
            "Fs" => result["Fs"][start_idx:end_idx],
            "Fc" => result["Fc"][start_idx:end_idx],
            "Fd" => result["Fd"][start_idx:end_idx],
            "Fe" => result["Fe"][start_idx:end_idx],
            "Fext" => result["Fext"][start_idx:end_idx],
            "CollisionState" => result["CollisionState"][start_idx:end_idx]
        )
    end
    
    # Create comparison plots for key variables
    variables = ["x1", "x1dot", "x2", "x2dot", "V", "Qvar", "Fc", "Fd"]
    titles = ["Shuttle Displacement", "Shuttle Velocity", "Electrode Displacement", 
              "Electrode Velocity", "Output Voltage", "Charge", "Collision Force", "Damping Force"]
    units = [L"$\mu$m", L"mm/s", L"$\mu$m", L"mm/s", L"mV", L"pC", L"$\mu$N", L"$\mu$N"]
    scale_factors = [1e6, 1e3, 1e6, 1e3, 1e3, 1e12, 1e6, 1e6]  # Convert to appropriate units
    
    for (var, title, unit, scale) in zip(variables, titles, units, scale_factors)
        p = plot(
            size=(800, 600),
            dpi=600,
            left_margin=10mm,
            bottom_margin=10mm,
            top_margin=5mm,
            right_margin=5mm,
            title="Effect of Acceleration on $title",
            titlelocation=:center,
            titlefontsize=14,
            xlabel=L"Time $t$ (s)",
            ylabel="$title $unit",
            legend=:topright
        )
        
        # Add a line for each alpha value
        for alpha in alpha_values
            data = trimmed_data[alpha]
            plot!(p, data["time"], data[var] .* scale, 
                  linewidth=2, 
                  label=L"\alpha = %$alpha")
        end
        
        savefig(p, joinpath(alpha_dir, "comparison_$(var).pdf"))
        savefig(p, joinpath(alpha_dir, "comparison_$(var).png"))
    end
    
    # Create RMS voltage vs. alpha plot
    rms_voltages = Float64[]
    for alpha in alpha_values
        data = trimmed_data[alpha]
        rms_voltage = sqrt(mean(data["V"].^2)) * 1e3  # Convert to mV
        push!(rms_voltages, rms_voltage)
    end
    
    p_rms = plot(
        alpha_values, 
        rms_voltages,
        linewidth=1.15,
        marker=:circle,
        markersize=6,
        color=:darkorange2,
        xlabel=L"Acceleration Parameter $\alpha$",
        ylabel=L"RMS Voltage $V_{rms}$ (mV)",
        title="RMS Voltage vs. Acceleration",
        legend=false,
        dpi=600,
        size=(800, 600),
        left_margin=10mm,
        bottom_margin=10mm,
        top_margin=5mm,
        right_margin=5mm
    )
    
    savefig(p_rms, joinpath(alpha_dir, "rms_voltage_vs_alpha.pdf"))
    savefig(p_rms, joinpath(alpha_dir, "rms_voltage_vs_alpha.png"))
    
    # Create max displacement vs. alpha plot
    max_displacements = Float64[]
    for alpha in alpha_values
        data = trimmed_data[alpha]
        max_displacement = maximum(abs.(data["x1"])) * 1e6  # Convert to μm
        push!(max_displacements, max_displacement)
    end
    
    p_disp = plot(
        alpha_values, 
        max_displacements,
        linewidth=1.15,
        marker=:circle,
        markersize=6,
        color=:royalblue3,
        xlabel=L"Acceleration Parameter $\alpha$",
        ylabel=L"Maximum Displacement $x_{1,max}$ ($\mu$m)",
        title="Maximum Displacement vs. Acceleration",
        legend=false,
        dpi=600,
        size=(800, 600),
        left_margin=10mm,
        bottom_margin=10mm,
        top_margin=5mm,
        right_margin=5mm
    )
    
    savefig(p_disp, joinpath(alpha_dir, "max_displacement_vs_alpha.pdf"))
    savefig(p_disp, joinpath(alpha_dir, "max_displacement_vs_alpha.png"))
    
    # Create summary plot with both RMS voltage and max displacement
    p_summary = plot(
        layout=(2,1),
        size=(800, 800),
        dpi=600,
        left_margin=10mm,
        bottom_margin=10mm,
        top_margin=5mm,
        right_margin=5mm,
        plot_title="Effect of Acceleration on MEMS Energy Harvester Performance",
        plot_titlelocation=:center,
        plot_titlefontsize=14
    )
    
    plot!(p_summary[1], alpha_values, rms_voltages, linewidth=2.5, marker=:circle, markersize=6, 
          color=:darkorange2, xlabel="", ylabel=L"RMS Voltage $V_{rms}$ (mV)", 
          title="RMS Voltage vs. Acceleration", legend=false)
    
    plot!(p_summary[2], alpha_values, max_displacements, linewidth=2.5, marker=:circle, markersize=6, 
          color=:royalblue3, xlabel=L"Acceleration Parameter $\alpha$", 
          ylabel=L"Maximum Displacement $x_{1,max}$ ($\mu$m)", 
          title="Maximum Displacement vs. Acceleration", legend=false)
    
    savefig(p_summary, joinpath(alpha_dir, "performance_summary.pdf"))
    savefig(p_summary, joinpath(alpha_dir, "performance_summary.png"))
    
    return trimmed_data
end

# Function to generate a combined grid plot showing all alphas together on each state variable subplot
function generate_combined_grid_plot(all_results, alpha_values, time_window, num_oscillations, start_time)
    # Check if we have enough alpha values for comparison
    if length(alpha_values) < 2
        println("Need at least 2 alpha values for comparison plots. Skipping.")
        return
    end
    
    # Create a directory for alpha comparison if it doesn't exist
    alpha_dir = "alpha_comparison"
    if !isdir(alpha_dir)
        mkdir(alpha_dir)
    end
    
    # Store trimmed data for all alpha values
    trimmed_data = Dict()
    
    # Process each alpha value to get the time window and extract data
    for alpha in alpha_values
        result = all_results[alpha]
        
        # Determine the time window
        start_t, end_t = determine_time_window(result, time_window, num_oscillations, start_time)
        
        # Find indices for the specified time window
        time = result["time"]
        start_idx = findfirst(t -> t >= start_t, time)
        end_idx = findfirst(t -> t >= end_t, time)
        
        if start_idx === nothing
            start_idx = 1
        end
        if end_idx === nothing
            end_idx = length(time)
        end
        
        # Extract and store the trimmed data
        trimmed_data[alpha] = Dict(
            "time" => time[start_idx:end_idx],
            "x1" => result["x1"][start_idx:end_idx],
            "x1dot" => result["x1dot"][start_idx:end_idx],
            "x2" => result["x2"][start_idx:end_idx],
            "x2dot" => result["x2dot"][start_idx:end_idx],
            "Qvar" => result["Qvar"][start_idx:end_idx],
            "V" => result["V"][start_idx:end_idx]
        )
    end
    
    # Define variables for the grid plot
    variables = ["x1", "x1dot", "x2", "x2dot", "Qvar", "V"]
    titles = ["Shuttle Displacement", "Shuttle Velocity", "Electrode Displacement", 
              "Electrode Velocity", "Charge", "Output Voltage"]
    units = [L"$\mu$m", L"mm/s", L"$\mu$m", L"mm/s", L"pC", L"mV"]
    scale_factors = [1e6, 1e3, 1e6, 1e3, 1e12, 1e3]  # Convert to appropriate units
    
    # Create a 2×3 grid where each subplot shows all alpha values
    combined_plot = plot(
        layout=(2, 3),
        size=(1200, 800),
        dpi=600,
        left_margin=10mm,
        bottom_margin=10mm,
        top_margin=5mm,
        right_margin=5mm,
        plot_titlelocation=:center,
        plot_titlefontsize=16
    )
    
    # For each state variable (subplot position), plot all alpha values
    for (i, (var, title, unit, scale)) in enumerate(zip(variables, titles, units, scale_factors))
        # Create a subplot for this state variable
        subplot_idx = i
        
        # Initialize the subplot
        plot!(
            combined_plot[subplot_idx],
            xlabel=(i > 3) ? L"Time $t$ (s)" : "",  # Only show x-label for bottom row
            ylabel="$title $unit",
            title=title,
            legend=(i == 1) ? :topright : false,  # Only show legend in first subplot
            framestyle=:box,
            grid=false
        )
        
        # Plot each alpha value as a separate line on this subplot
        for (j, alpha) in enumerate(alpha_values)
            data = trimmed_data[alpha]
            
            plot!(
                combined_plot[subplot_idx],
                data["time"],
                data[var] .* scale,
                linewidth=1,
                label=(i == 1) ? L"\gamma = %$alpha" : "",  # Only add labels in first subplot
                color=j  # Use different colors for different alpha values
            )
        end
    end
    
    # Save the combined grid plot
    savefig(combined_plot, joinpath(alpha_dir, "combined_state_comparison_grid.pdf"))
    savefig(combined_plot, joinpath(alpha_dir, "combined_state_comparison_grid.png"))
    
    # Create a separate set of individual state plots (one plot per state variable, all alphas)
    for (i, (var, title, unit, scale)) in enumerate(zip(variables, titles, units, scale_factors))
        p_var = plot(
            size=(800, 600),
            dpi=600,
            left_margin=10mm,
            bottom_margin=10mm,
            top_margin=5mm,
            right_margin=5mm,
            title="Effect of Acceleration on $title",
            xlabel=L"Time $t$ (s)",
            ylabel="$title $unit",
            legend=:topright,
            framestyle=:box,
            grid=false
        )
        
        # Plot each alpha value
        for (j, alpha) in enumerate(alpha_values)
            data = trimmed_data[alpha]
            
            plot!(
                p_var,
                data["time"],
                data[var] .* scale,
                linewidth=1,
                label=L"\gamma = %$alpha",
                color=j
            )
        end
        
        # Save this individual comparison plot
        savefig(p_var, joinpath(alpha_dir, "comparison_$(var)_overlay.pdf"))
        savefig(p_var, joinpath(alpha_dir, "comparison_$(var)_overlay.png"))
    end
end

# Function to create a special plot with voltage and external force overlay
function plot_voltage_with_force(result, alpha, time_range=(0.275, 0.327); prefix="")
    # Create directory for this alpha value if it doesn't exist
    alpha_dir = "alpha_$(alpha)"
    if !isdir(alpha_dir)
        mkdir(alpha_dir)
    end
    
    # Extract the time range
    time = result["time"]
    start_time, end_time = time_range
    
    # Find indices for the specified time window
    start_idx = findfirst(t -> t >= start_time, time)
    end_idx = findfirst(t -> t >= end_time, time)
    
    if start_idx === nothing
        start_idx = 1
    end
    if end_idx === nothing
        end_idx = length(time)
    end
    
    # Extract data for the specified time window
    plot_time = time[start_idx:end_idx]
    plot_V = result["V"][start_idx:end_idx]
    plot_Fext = result["Fext"][start_idx:end_idx]
    
    # Create plot
    p_voltage_force = plot(
        size=(900, 600),
        dpi=600,
        left_margin=15mm,
        bottom_margin=10mm,
        top_margin=10mm,
        right_margin=10mm,
        titlefontsize=14,
        titlelocation=:center
    )
    
    # Plot the voltage on the primary y-axis
    plot!(
        p_voltage_force,
        plot_time,
        plot_V .* 1e3,  # Convert to mV
        linewidth=2.5,
        color=:royalblue3,
        xlabel=L"Time $t$ (s)",
        ylabel=L"Output Voltage, $V_{\rm out}$ (mV)",
        label="Voltage",
        framestyle=:box,
        legend=:topright,
        grid=false
    )
    
    # Calculate a scaling factor for the external force to fit it within the voltage range
    # This ensures the force plot is visible but doesn't dominate
    v_range = maximum(plot_V .* 1e3) - minimum(plot_V .* 1e3)
    f_range = maximum(plot_Fext) - minimum(plot_Fext)
    f_scale = v_range / f_range * 0.6  # Scale to ~60% of voltage range
    
    # Add the external force as a dashed line - no right y-axis
    # Use broadcasting (.+) when combining scalar and vector
    plot!(
        p_voltage_force,
        plot_time,
        minimum(plot_V .* 1e3) .+ (plot_Fext .- minimum(plot_Fext)) .* f_scale,  # Fix: added dot for broadcasting
        linewidth=2.5,
        linestyle=:dash,
        color=:black,
        label="External Force",
        grid=false
    )
    
    # Save the plot
    savefig(p_voltage_force, joinpath(alpha_dir, "$(prefix)voltage_with_force_alpha_$(alpha).pdf"))
    savefig(p_voltage_force, joinpath(alpha_dir, "$(prefix)voltage_with_force_alpha_$(alpha).png"))
    
    return p_voltage_force
end

# Function to plot relative displacement (x1-x2) for multiple alpha values
function plot_relative_displacement(all_results, alpha_values, time_range=(0.2, 0.45); prefix="")
    # Create directory for alpha comparison if it doesn't exist
    alpha_dir = "alpha_comparison"
    if !isdir(alpha_dir)
        mkdir(alpha_dir)
    end
    
    # Initialize the plot
    p_rel_disp = plot(
        size=(900, 600),
        dpi=600,
        left_margin=15mm,
        bottom_margin=10mm,
        top_margin=10mm,
        right_margin=10mm,
        titlefontsize=14,
        titlelocation=:center,
        xlabel=L"Time $t$ (s)",
        ylabel=L"Relative Displacement, $x_1 - x_2$ ($\mu$m)",
        legend=:topright,
        framestyle=:box,
        grid=false
    )
    
    # Process each alpha value
    for (i, alpha) in enumerate(alpha_values)
        result = all_results[alpha]
        time = result["time"]
        start_time, end_time = time_range
        
        # Find indices for the specified time window
        start_idx = findfirst(t -> t >= start_time, time)
        end_idx = findfirst(t -> t >= end_time, time)
        
        if start_idx === nothing
            start_idx = 1
        end
        if end_idx === nothing
            end_idx = length(time)
        end
        
        # Extract data for the specified time window
        plot_time = time[start_idx:end_idx]
        plot_x1 = result["x1"][start_idx:end_idx]
        plot_x2 = result["x2"][start_idx:end_idx]
        
        # Calculate relative displacement and convert to μm
        rel_disp = (plot_x1 .- plot_x2) .* 1e6
        
        # Plot with different color/style for each alpha
        plot!(
            p_rel_disp,
            plot_time,
            rel_disp,
            linewidth=2.5,
            color=i,  # Different color for each alpha
            label=L"\gamma = %$alpha",
            grid=false
        )
    end
    
    # Add a horizontal line at zero to show when shuttle and electrode are at the same position
    hline!(
        p_rel_disp,
        [0],
        linestyle=:dash,
        color=:black,
        alpha=0.5,
        linewidth=1.5,
        label=""
    )
    
    # Save the plot
    savefig(p_rel_disp, joinpath(alpha_dir, "$(prefix)relative_displacement_comparison.pdf"))
    savefig(p_rel_disp, joinpath(alpha_dir, "$(prefix)relative_displacement_comparison.png"))
    
    # Create an alternative version with normalized time (showing phase differences more clearly)
    p_rel_disp_norm = plot(
        size=(900, 600),
        dpi=600,
        left_margin=15mm,
        bottom_margin=10mm,
        top_margin=10mm,
        right_margin=10mm,
        titlefontsize=14,
        titlelocation=:center,
        xlabel=L"Normalized Time $t/T$",
        ylabel=L"Relative Displacement $x_1 - x_2$ ($\mu$m)",
        legend=:topright,
        framestyle=:box,
        grid=false
    )
    
    # Assume the frequency is the same for all alpha values (get from first result)
    base_freq = all_results[alpha_values[1]]["frequency"]
    period = 1.0 / base_freq
    
    # Process each alpha value again for normalized time plot
    for (i, alpha) in enumerate(alpha_values)
        result = all_results[alpha]
        time = result["time"]
        start_time, end_time = time_range
        
        # Find indices for the specified time window
        start_idx = findfirst(t -> t >= start_time, time)
        end_idx = findfirst(t -> t >= end_time, time)
        
        if start_idx === nothing
            start_idx = 1
        end
        if end_idx === nothing
            end_idx = length(time)
        end
        
        # Extract data for the specified time window
        plot_time = time[start_idx:end_idx]
        plot_x1 = result["x1"][start_idx:end_idx]
        plot_x2 = result["x2"][start_idx:end_idx]
        
        # Calculate relative displacement and convert to μm
        rel_disp = (plot_x1 .- plot_x2) .* 1e6
        
        # Normalize time by dividing by the period
        norm_time = (plot_time .- plot_time[1]) ./ period
        
        # Only keep data for a few complete cycles for clarity
        max_cycles = 4
        cycle_indices = findall(t -> t <= max_cycles, norm_time)
        
        if !isempty(cycle_indices)
            norm_time = norm_time[cycle_indices]
            rel_disp = rel_disp[cycle_indices]
            
            # Plot with different color/style for each alpha
            plot!(
                p_rel_disp_norm,
                norm_time,
                rel_disp,
                linewidth=1.15,
                color=i,  # Different color for each alpha
                label=L"\alpha = %$alpha",
                grid=false
            )
        end
    end
    
    # Add a horizontal line at zero
    hline!(
        p_rel_disp_norm,
        [0],
        linestyle=:dash,
        color=:black,
        alpha=0.5,
        linewidth=1.15,
        label=""
    )
    
    # Save the normalized time plot
    savefig(p_rel_disp_norm, joinpath(alpha_dir, "$(prefix)relative_displacement_normalized_time.pdf"))
    savefig(p_rel_disp_norm, joinpath(alpha_dir, "$(prefix)relative_displacement_normalized_time.png"))
    
    return p_rel_disp, p_rel_disp_norm
end

p_rel, p_rel_norm = plot_relative_displacement(all_results, alpha_values, (0.275, 0.327))

alpha = 3.25
result = all_results[alpha]
plot_voltage_with_force(result, alpha, (0.275, 0.327), prefix="special_")

# Option 1: Frequency-based time window (auto-calculated based on number of oscillations)
alpha_values = [3.25]
# Option 1: Frequency-based time window (auto-calculated based on number of oscillations)
# all_results = generate_alpha_plots(alpha_values, time_window=:auto, num_oscillations=2, start_time=0.275)
# Option 2: Explicit time window (directly specify start and end times)
all_results = generate_alpha_plots(alpha_values, time_window=(0.275, 0.327), prefix="fixed_")      