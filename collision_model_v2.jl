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


 # --------------------------------------- Frequency Sweep ---------------------------------------
using Statistics, Plots, Printf, LaTeXStrings, StatsBase

# Define the frequency range and increments
freq_min = 5.0
freq_max = 150.0
num_points = 75  # Increase for finer resolution if needed
freq_values = range(freq_min, freq_max, length=num_points)

# Initialize arrays to store outputs
ss_times = Float64[]
max_vout_values = Float64[]
max_z1_values = Float64[]
max_z1dot_values = Float64[]
max_z2_values = Float64[]
max_z2dot_values = Float64[]
rms_vout_values = Float64[]

# For storing forces at resonance
resonance_forces = Dict()
resonance_states = Dict()
resonance_time = nothing

# Parameters for detection (similar to MATLAB code)
peak_threshold = 5e-9      # Threshold for peak detection
steady_state_threshold = 2.5e-9  # Threshold for steady state detection
min_steady_state_time = 0.2      # Minimum time to start considering steady state
num_peaks_to_analyze = 35        # Number of peaks to analyze in each iteration

println("Starting frequency sweep from $(freq_min) Hz to $(freq_max) Hz...")
println("-----------------------------------------------------")
println("| Frequency (Hz) | RMS Voltage (V) | Max Disp. (µm) |")
println("-----------------------------------------------------")

# Run the model for each frequency value
for (i, current_freq) in enumerate(freq_values)
    # Set the frequency for this iteration
    f = current_freq  # Update the frequency variable
    
    # Set the external force parameters
    alpha = 1.0  # 1g acceleration
    g = 9.81     # Gravitational constant 
    A = alpha * g
    t_ramp = 0.2  # Ramp-up time
    
    # Update the external force function (sine wave with frequency)
    Fext_input = t -> A * (t < t_ramp ? t / t_ramp : 1.0) * sin(2 * π * current_freq * t)
    
    # Initialize parameters 
    p_new = deepcopy(AnalyticalModel.p)
    
    # Set simulation time (matching MATLAB's simulation time)
    tspan = (0.0, 75.0)  # Long simulation for steady state
    
    # Initial conditions
    x10 = 0.0      # Initial displacement
    x10dot = 0.0   # Initial velocity
    x20 = 0.0      # Initial electrode displacement
    x20dot = 0.0   # Initial electrode velocity
    
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
    
    # ODE solver settings
    abstol = 1e-9  # Absolute tolerance
    reltol = 1e-6  # Relative tolerance
    
    # Define the wrapper function for the ODE solver
    function CoupledSystem_wrapper!(dz, z, p, t)
        current_acceleration = Fext_input(t)
        AnalyticalModel.CoupledSystem!(dz, z, p, t, current_acceleration)
    end
    
    # Solve the ODE system (similar to MATLAB's ode45 solver)
    eqn = ODEProblem(CoupledSystem_wrapper!, z0, tspan, p_new)
    sol = solve(eqn, Rosenbrock23(); abstol=abstol, reltol=reltol, maxiters=1e7)
    
    # Extract the results
    time = sol.t
    z1 = [u[1] for u in sol.u]      # Shuttle displacement
    z1dot = [u[2] for u in sol.u]   # Shuttle velocity
    z2 = [u[3] for u in sol.u]      # Electrode displacement
    z2dot = [u[4] for u in sol.u]   # Electrode velocity
    qvar = [u[5] for u in sol.u]    # Charge
    vout = [u[6] for u in sol.u]    # Output voltage
    
    # Find the index corresponding to the minimum time
    min_time_index = findfirst(t -> t >= min_steady_state_time, time)
    if min_time_index === nothing
        min_time_index = 1
    end
    
    # Trim the data to start from the specified minimum time
    z1_trimmed = z1[min_time_index:end]
    z1dot_trimmed = z1dot[min_time_index:end]
    z2_trimmed = z2[min_time_index:end]
    z2dot_trimmed = z2dot[min_time_index:end]
    qvar_trimmed = qvar[min_time_index:end]
    vout_trimmed = vout[min_time_index:end]
    time_trimmed = time[min_time_index:end]
    
    # Find the indices of all peaks in the displacement data (similar to findpeaks)
    peak_indices = Int[]
    for i in 2:(length(z1_trimmed)-1)
        if z1_trimmed[i] > z1_trimmed[i-1] && z1_trimmed[i] > z1_trimmed[i+1] && z1_trimmed[i] > peak_threshold
            push!(peak_indices, i)
        end
    end
    
    # Initialize steady state detected flag
    steady_state_detected = false
    ss_index = length(time_trimmed)  # Default to end of simulation
    
    # Continue analyzing peaks until steady state is detected or no more peaks
    if !isempty(peak_indices)
        # Analyze peaks in chunks of num_peaks_to_analyze
        for start_idx in 1:num_peaks_to_analyze:length(peak_indices)
            end_idx = min(start_idx + num_peaks_to_analyze - 1, length(peak_indices))
            
            if end_idx <= length(peak_indices)
                current_peaks = z1_trimmed[peak_indices[start_idx:end_idx]]
                
                if length(current_peaks) >= 2
                    # Separate odd and even peaks (like MATLAB code)
                    odd_idxs = collect(1:2:length(current_peaks))
                    even_idxs = collect(2:2:length(current_peaks))
                    
                    odd_values = length(odd_idxs) > 0 ? current_peaks[odd_idxs] : Float64[]
                    even_values = length(even_idxs) > 0 ? current_peaks[even_idxs] : Float64[]
                    
                    # Check if peak variation is within threshold
                    if !isempty(odd_values) && !isempty(even_values)
                        diff_odd = maximum(odd_values) - minimum(odd_values)
                        diff_even = maximum(even_values) - minimum(even_values)
                        
                        if diff_odd < steady_state_threshold || diff_even < steady_state_threshold
                            steady_state_detected = true
                            ss_index = peak_indices[start_idx]
                            push!(ss_times, time_trimmed[ss_index])
                            break  # Exit the loop if steady state is detected
                        end
                    end
                end
            end
        end
    end
    
    # Record values based on steady state detection
    if steady_state_detected
        # Calculate RMS voltage from steady state onward
        vout_ss = vout_trimmed[ss_index:end]
        rms_vout = sqrt(mean(vout_ss.^2))
        
        # Record maximum values
        push!(max_vout_values, maximum(abs.(vout_ss)))
        push!(max_z1_values, maximum(abs.(z1_trimmed[ss_index:end])))
        push!(max_z1dot_values, maximum(abs.(z1dot_trimmed[ss_index:end])))
        push!(max_z2_values, maximum(abs.(z2_trimmed[ss_index:end])))
        push!(max_z2dot_values, maximum(abs.(z2dot_trimmed[ss_index:end])))
        push!(rms_vout_values, rms_vout)
        
        # Print progress
        @printf("| %12.2f | %14.6f | %13.2f |\n", 
                current_freq, rms_vout, maximum(abs.(z1_trimmed[ss_index:end])) * 1e6)
                
        # If this is close to resonance (highest RMS voltage so far), store the forces
        if length(rms_vout_values) == 1 || rms_vout > maximum(rms_vout_values[1:end-1])
            # Store the full solution at resonance for later force plots
            resonance_forces = Dict(
                "time" => time,
                "z1" => z1,
                "z1dot" => z1dot,
                "z2" => z2,
                "z2dot" => z2dot,
                "vout" => vout,
                "qvar" => qvar,
                "frequency" => current_freq
            )
            resonance_time = time
        end
    else
        # If steady state not detected, use the last section of simulation
        last_quarter_idx = round(Int, length(time_trimmed) * 0.75)
        
        # Calculate values from the last quarter
        vout_last = vout_trimmed[last_quarter_idx:end]
        rms_vout = sqrt(mean(vout_last.^2))
        
        # Add values
        push!(ss_times, NaN)
        push!(max_vout_values, maximum(abs.(vout_last)))
        push!(max_z1_values, maximum(abs.(z1_trimmed[last_quarter_idx:end])))
        push!(max_z1dot_values, maximum(abs.(z1dot_trimmed[last_quarter_idx:end])))
        push!(max_z2_values, maximum(abs.(z2_trimmed[last_quarter_idx:end])))
        push!(max_z2dot_values, maximum(abs.(z2dot_trimmed[last_quarter_idx:end])))
        push!(rms_vout_values, rms_vout)
        
        # Print progress with warning
        @printf("| %12.2f | %14.6f | %13.2f | (SS not detected)\n", 
                current_freq, rms_vout, maximum(abs.(z1_trimmed[last_quarter_idx:end])) * 1e6)
    end
end

println("-----------------------------------------------------")
println("Frequency sweep completed.")

# Now extract the resonance data for force calculation
if !isempty(resonance_forces)
    # Calculate forces at resonance frequency
    res_time = resonance_forces["time"]
    res_z1 = resonance_forces["z1"]
    res_z1dot = resonance_forces["z1dot"]
    res_z2 = resonance_forces["z2"]
    res_z2dot = resonance_forces["z2dot"]
    res_qvar = resonance_forces["qvar"]
    res_freq = resonance_forces["frequency"]
    
    # Initialize arrays for forces
    Fs_array = Float64[]  # Suspension spring force
    Fc_array = Float64[]  # Collision force
    Fd_array = Float64[]  # Damping force
    Fe_array = Float64[]  # Electrostatic force
    Fext_array = Float64[] # External force
    CollisionState_array = String[]  # Collision state
    
    # Calculate forces at each time point
    p_new = deepcopy(AnalyticalModel.p)
    
    for i in 1:length(res_time)
        # Extract state variables
        z1 = res_z1[i]
        z1dot = res_z1dot[i]
        z2 = res_z2[i]
        z2dot = res_z2dot[i]
        qvar = res_qvar[i]
        
        # Calculate external force at this time
        t = res_time[i]
        alpha = 1.0
        g = 9.81
        A = alpha * g
        t_ramp = 0.2
        Fext = A * (t < t_ramp ? t / t_ramp : 1.0) * sin(2 * π * res_freq * t)
        push!(Fext_array, Fext)
        
        # Calculate spring force
        Fs = AnalyticalModel.spring(z1, p_new.k1, p_new.k3, p_new.gss, p_new.kss)
        push!(Fs_array, Fs)
        
        # Calculate collision force
        Fc, collision_state = AnalyticalModel.collision(z1, z2, p_new.m2, p_new.ke, p_new.gp)
        push!(Fc_array, Fc)
        push!(CollisionState_array, collision_state)
        
        # Calculate damping force
        if collision_state == "translational"
            Fd = AnalyticalModel.damping(z1, z1dot, z2, z2dot, p_new.a, p_new.c, p_new.gp, p_new.Leff, p_new.Tf, p_new.eta, p_new.lambda, p_new.sigmap)
        else
            Fd = AnalyticalModel.damping(z1, z1dot, z2, z1dot, p_new.a, p_new.c, p_new.gp, p_new.Leff, p_new.Tf, p_new.eta, p_new.lambda, p_new.sigmap)
        end
        push!(Fd_array, Fd)
        
        # Calculate electrostatic force
        Ctotal, Fe = AnalyticalModel.electrostatic(
            z1, z2, qvar, 
            p_new.g0, p_new.gp, p_new.a, 
            p_new.e, p_new.ep, p_new.cp, 
            p_new.wt, p_new.wb, p_new.ke, 
            p_new.E, p_new.I, p_new.Leff, 
            p_new.Tf, p_new.Tp, p_new.N
        )
        push!(Fe_array, Fe)
    end
    
    # Store the forces for plotting
    resonance_forces["Fs"] = Fs_array
    resonance_forces["Fc"] = Fc_array
    resonance_forces["Fd"] = Fd_array
    resonance_forces["Fe"] = Fe_array
    resonance_forces["Fext"] = Fext_array
    resonance_forces["CollisionState"] = CollisionState_array
end

# Define a high-quality theme for journal publications
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
        margin=8mm,                    # Margins around the plot
        grid=false,                    # No grid by default
        framestyle=:box,               # Box-style frame
        palette=:default               # Color palette
    )
end

# Apply the journal theme
set_journal_theme()

# 1. Frequency Response Plots
# First, create the RMS voltage vs frequency plot
p_voltage = plot(
    freq_values, 
    rms_vout_values .* 1000,       # Convert to mV
    seriestype=:line,
    linewidth=2,
    color=:royalblue3,
    xlabel=L"Frequency $f$ (Hz)",
    ylabel=L"RMS Voltage $V_{rms}$ (mV)",
    title="",                      # Journal figures often don't have titles
    label="",                      # No label needed for single series
    xticks=0:20:160,               # Custom x-ticks
    yticks=0:2:20,                 # Custom y-ticks (adjust based on your data)
    xlims=(freq_min, freq_max),
    ylims=(0, nothing),            # Start y-axis at 0
    minorgrid=false,
    legend=false,                  # No legend for single series
    framestyle=:box,
    foreground_color_axis=:black,
    tick_direction=:out,           # Ticks pointing outward
    guidefontsize=12
)

# Add data points (markers) on top of the line
scatter!(
    p_voltage,
    freq_values, 
    rms_vout_values .* 1000,
    markersize=5,
    markerstrokewidth=0.5,
    markerstrokecolor=:black,
    markercolor=:royalblue3,
    markershape=:circle,
    label=""
)

# Create the displacement plot with the same styling
p_displacement = plot(
    freq_values, 
    max_z1_values .* 1e6,          # Convert to μm
    seriestype=:line,
    linewidth=2,
    color=:firebrick3,
    xlabel=L"Frequency $f$ (Hz)",
    ylabel=L"Maximum Displacement $x_{max}$ ($\mu$m)",
    title="",
    label="",
    xticks=0:20:160,
    yticks=0:5:30,                 # Adjust based on your data
    xlims=(freq_min, freq_max),
    ylims=(0, nothing),            # Start y-axis at 0
    minorgrid=false,
    legend=false,
    framestyle=:box,
    foreground_color_axis=:black,
    tick_direction=:out,
    guidefontsize=12
)

# Add data points to displacement plot
scatter!(
    p_displacement,
    freq_values, 
    max_z1_values .* 1e6,
    markersize=5,
    markerstrokewidth=0.5,
    markerstrokecolor=:black,
    markercolor=:firebrick3,
    markershape=:circle,
    label=""
)

# Create a multi-panel figure for frequency response
p_freq_combined = plot(
    p_voltage, 
    p_displacement, 
    layout=(2,1),
    size=(800, 800),
    dpi=600,
    margin=10mm,
    link=:x                        # Link x-axes
)

# Save frequency response plots
savefig(p_voltage, "voltage_response.pdf")
savefig(p_voltage, "voltage_response.png")
savefig(p_displacement, "displacement_response.pdf")
savefig(p_displacement, "displacement_response.png")
savefig(p_freq_combined, "frequency_response_combined.pdf")
savefig(p_freq_combined, "frequency_response_combined.png")

# 2. State Variable Plots at Resonance
if !isempty(resonance_forces)
    # Time for state plots - use section around steady state
    res_time = resonance_forces["time"]
    start_time = 0.25  # Start after initial transients
    end_time = 0.35    # Just show a small window for clarity
    
    time_indices = findall(t -> start_time <= t <= end_time, res_time)
    
    # Extract state variables in the time window
    plot_time = res_time[time_indices]
    plot_z1 = resonance_forces["z1"][time_indices]
    plot_z1dot = resonance_forces["z1dot"][time_indices]
    plot_z2 = resonance_forces["z2"][time_indices]
    plot_z2dot = resonance_forces["z2dot"][time_indices]
    plot_vout = resonance_forces["vout"][time_indices]
    
    # Get resonance frequency
    res_freq = resonance_forces["frequency"]
    
    # Create a multi-panel state plot
    p_states = plot(
        layout=(3,2),
        size=(900, 1000),
        dpi=600,
        left_margin=10mm,
        bottom_margin=10mm,
        title=["Shuttle Displacement" "Shuttle Velocity" "Electrode Displacement" "Electrode Velocity" "Output Voltage" "Relative Displacement"],
        titlelocation=:center,
        titlefontsize=10
    )
    
    # Plot each state variable
    # 1. Shuttle displacement
    plot!(
        p_states[1],
        plot_time, 
        plot_z1 .* 1e6,  # Convert to μm
        linewidth=2,
        color=:royalblue3,
        xlabel="",
        ylabel=L"x_1 ($\mu$m)",
        legend=false,
        framestyle=:box,
        xlims=(start_time, end_time),
        xticks=start_time:0.02:end_time
    )
    
    # 2. Shuttle velocity
    plot!(
        p_states[2],
        plot_time, 
        plot_z1dot .* 1e3,  # Convert to mm/s
        linewidth=2,
        color=:firebrick3,
        xlabel="",
        ylabel=L"\dot{x}_1 (mm/s)",
        legend=false,
        framestyle=:box,
        xlims=(start_time, end_time),
        xticks=start_time:0.02:end_time
    )
    
    # 3. Electrode displacement
    plot!(
        p_states[3],
        plot_time, 
        plot_z2 .* 1e6,  # Convert to μm
        linewidth=2,
        color=:green4,
        xlabel="",
        ylabel=L"x_2 ($\mu$m)",
        legend=false,
        framestyle=:box,
        xlims=(start_time, end_time),
        xticks=start_time:0.02:end_time
    )
    
    # 4. Electrode velocity
    plot!(
        p_states[4],
        plot_time, 
        plot_z2dot .* 1e3,  # Convert to mm/s
        linewidth=2,
        color=:purple3,
        xlabel="",
        ylabel=L"\dot{x}_2 (mm/s)",
        legend=false,
        framestyle=:box,
        xlims=(start_time, end_time),
        xticks=start_time:0.02:end_time
    )
    
    # 5. Output voltage
    plot!(
        p_states[5],
        plot_time, 
        plot_vout .* 1e3,  # Convert to mV
        linewidth=2,
        color=:darkorange2,
        xlabel=L"Time $t$ (s)",
        ylabel=L"V_{out} (mV)",
        legend=false,
        framestyle=:box,
        xlims=(start_time, end_time),
        xticks=start_time:0.02:end_time
    )
    
    # 6. Relative displacement (x1-x2)
    plot!(
        p_states[6],
        plot_time, 
        (plot_z1 - plot_z2) .* 1e6,  # Convert to μm
        linewidth=2,
        color=:teal,
        xlabel=L"Time $t$ (s)",
        ylabel=L"x_1 - x_2 ($\mu$m)",
        legend=false,
        framestyle=:box,
        xlims=(start_time, end_time),
        xticks=start_time:0.02:end_time
    )
    
    # Add a common title
    plot!(
        p_states,
        title=L"State Variables at Resonance Frequency $f = %$(round(Int, res_freq))$ Hz",
        titlelocation=:center,
        titlefontsize=14
    )
    
    # Save the state variables plot
    savefig(p_states, "state_variables_resonance.pdf")
    savefig(p_states, "state_variables_resonance.png")
    
    # 3. Force Plots at Resonance
    # Extract forces in the time window
    plot_Fs = resonance_forces["Fs"][time_indices]
    plot_Fc = resonance_forces["Fc"][time_indices]
    plot_Fd = resonance_forces["Fd"][time_indices]
    plot_Fe = resonance_forces["Fe"][time_indices]
    plot_Fext = resonance_forces["Fext"][time_indices]
    plot_CollisionState = resonance_forces["CollisionState"][time_indices]
    
    # Create a multi-panel force plot
    p_forces = plot(
        layout=(3,2),
        size=(900, 1000),
        dpi=600,
        left_margin=10mm,
        bottom_margin=10mm,
        title=["Spring Force" "Collision Force" "Damping Force" "Electrostatic Force" "External Force" "Collision State"],
        titlelocation=:center,
        titlefontsize=10
    )
    
    # Plot each force
    # 1. Spring force
    plot!(
        p_forces[1],
        plot_time, 
        plot_Fs .* 1e6,  # Convert to μN
        linewidth=2,
        color=:royalblue3,
        xlabel="",
        ylabel=L"F_s ($\mu$N)",
        legend=false,
        framestyle=:box,
        xlims=(start_time, end_time),
        xticks=start_time:0.02:end_time
    )
    
    # 2. Collision force
    plot!(
        p_forces[2],
        plot_time, 
        plot_Fc .* 1e6,  # Convert to μN
        linewidth=2,
        color=:firebrick3,
        xlabel="",
        ylabel=L"F_c ($\mu$N)",
        legend=false,
        framestyle=:box,
        xlims=(start_time, end_time),
        xticks=start_time:0.02:end_time
    )
    
    # 3. Damping force
    plot!(
        p_forces[3],
        plot_time, 
        plot_Fd .* 1e6,  # Convert to μN
        linewidth=2,
        color=:green4,
        xlabel="",
        ylabel=L"F_d ($\mu$N)",
        legend=false,
        framestyle=:box,
        xlims=(start_time, end_time),
        xticks=start_time:0.02:end_time
    )
    
    # 4. Electrostatic force
    plot!(
        p_forces[4],
        plot_time, 
        plot_Fe .* 1e6,  # Convert to μN
        linewidth=2,
        color=:purple3,
        xlabel="",
        ylabel=L"F_e ($\mu$N)",
        legend=false,
        framestyle=:box,
        xlims=(start_time, end_time),
        xticks=start_time:0.02:end_time
    )
    
    # 5. External force
    plot!(
        p_forces[5],
        plot_time, 
        plot_Fext .* 1e6,  # Convert to μN
        linewidth=2,
        color=:darkorange2,
        xlabel=L"Time $t$ (s)",
        ylabel=L"F_{ext} ($\mu$N)",
        legend=false,
        framestyle=:box,
        xlims=(start_time, end_time),
        xticks=start_time:0.02:end_time
    )
    
    # 6. Collision state transitions
    plot!(
        p_forces[6],
        plot_time, 
        [s == "translational" ? 0 : 1 for s in plot_CollisionState],
        linewidth=2,
        color=:black,
        xlabel=L"Time $t$ (s)",
        ylabel="State",
        yticks=([0, 1], ["Trans.", "Rot."]),
        legend=false,
        framestyle=:box,
        xlims=(start_time, end_time),
        xticks=start_time:0.02:end_time
    )
    
    # Add a common title
    plot!(
        p_forces,
        title=L"Forces at Resonance Frequency $f = %$(round(Int, res_freq))$ Hz",
        titlelocation=:center,
        titlefontsize=14
    )
    
    # Save the force plot
    savefig(p_forces, "forces_resonance.pdf")
    savefig(p_forces, "forces_resonance.png")
end

# Calculate some statistical highlights from the data
peak_voltage_idx = argmax(rms_vout_values)
peak_voltage_freq = freq_values[peak_voltage_idx]
peak_voltage_value = rms_vout_values[peak_voltage_idx] * 1000  # Convert to mV

peak_disp_idx = argmax(max_z1_values)
peak_disp_freq = freq_values[peak_disp_idx]
peak_disp_value = max_z1_values[peak_disp_idx] * 1e6  # Convert to μm

println("Statistical highlights:")
println("Peak RMS voltage: $(round(peak_voltage_value, digits=2)) mV at $(round(peak_voltage_freq, digits=1)) Hz")
println("Peak displacement: $(round(peak_disp_value, digits=2)) μm at $(round(peak_disp_freq, digits=1)) Hz")
