# ------------------------- Libraries -------------------------
using Sundials, DifferentialEquations, Plots
using Parameters, ForwardDiff, SpecialFunctions, Random

# ------------------------- Parameter Struct -------------------------
# Parameters struct
@with_kw mutable struct MEMSParams{T<:Real}
    # Fundamental geometric parameters
    g0::T = 14e-6  # Initial gap
    Tp::T = 120e-9  # Parylene-C thickness
    Tf::T = 25e-6  # Electrode thickness
    wt::T = 9e-6  # Electrode width, top
    wb::T = 30e-6  # Electrode width, bottom
    ws::T = 14.7e-6  # Suspension spring width
    wss::T = 14e-6  # Soft-stopper width
    Leff::T = 400e-6  # Effective electrode length
    Lff::T = 450e-6  # Full electrode length
    Lsp::T = 1400e-6  # Suspension spring length
    Lss::T = 1000e-6  # Soft-stopper length
    gss::T = 14e-6  # Soft-stopper position
    
    # Mass and material properties
    m1::T = 2.0933e-6  # Shuttle mass
    rho::T = 2330.0  # Density of silicon
    E::T = 170e9  # Young's modulus
    eps::T = 8.85e-12  # Permittivity of free space
    eps_p::T = 3.2  # Relative permittivity of Parylene-C
    eta::T = 1.849e-5  # Viscosity of air
    c::T = 0.015  # Damping coefficient
    
    # Electrical parameters
    N::T = 160.0  # Number of electrodes
    Cp::T = 5e-12  # Capacitance of Parylene-C
    Vbias::T = 3.0  # Bias voltage
    Rload::T = 0.42e6  # Load resistance
    
    # Dependent parameters
    gp::T = :($(g0 - 2 * Tp))  # Initial electrode gap
    a::T = :($((wb - wt) / Leff))  # Taper ratio
    m2::T = 0.0  # Modal mass of electrode
    I::T = 0.0  # Electrode moment of inertia
    ke::T = 0.0  # Electrode spring constant
    k1::T = 0.0  # Linear spring constant
    k3::T = 0.0  # Cubic spring constant
    kss::T = 0.0  # Soft-stopper spring constant  
end

# Function to create dependent parameters for struct
function create_params(; verbose = false, kwargs...)
    # Create initial params with default values and simple calculations
    params = MEMSParams(; kwargs...)

    # Calculate electrode moment of inertia
    params.I = (1/48) * params.Lff * params.Tf * (params.wb + params.wt) * (params.wb^2 + params.wt^2)
    
    # Calculate and set modal mass, m2
    modalcoeff = 0.236 + 0.045 * (1.0 - params.wt / params.wb)
    m2Physical = 0.5 * params.Lff^2 * params.rho * params.Tf * (params.wb + params.wt)
    params.m2 = modalcoeff * m2Physical
    
    # Calculate electrode spring constant, ke
    F = 1 # Normalized force used in ke derivation
    num = (params.E * params.Tf * params.wt^2 * ((params.wb - params.wt)^3))
    dem = 6 * F * params.Lff^3 * ((params.wb - 3 * params.wt) * (params.wb - params.wt) + 2 * params.wt^2 * (log(params.Lff * params.wb) - log(params.Lff * params.wt)))
    params.ke = num / dem
    
    # Calculate spring constants, k1/k3/kss
    params.k1 = (4.0/6.0) * ((params.E * params.Tf * (params.ws^3)) / (params.Lsp^3))
    params.k3 = (18.0/25.0) * ((params.E * params.Tf * params.ws) / (params.Lsp^3))
    params.kss = (params.E * params.Tf * (params.wss^3)) / (params.Lss^3)

    # Calculation checks
    if verbose
        println("\n--- Modal Mass Calculation ---")
        println("Modal coefficient terms:  ", modalcoeff)
        println("Physical mass terms:  ", m2Physical)
        println("Modal mass, m2:  ", params.m2)
        
        println("\n--- Electrode Spring Constant Calculation ---")
        println("Moment of inertia, I:  ", params.I)
        println("Electrode spring constant, ke:  ", params.ke)
        
        println("\n--- Spring Constants ---")
        println("Linear suspension spring constant, k1:  ", params.k1)
        println("Cubic suspension spring constant, k3:  ", params.k3)
        println("Soft-stopper spring constant, kss:  ", params.kss)
    end
    
    return params
end

# Create the full set of model parameters 
params = create_params();

# ------------------------- Force Expressions -------------------------
# Spring force, Fs = Fsp + Fss
function spring(x1, params)
    # Suspension spring component
    Fsp = params.k1 * x1  # Linear
    # Fsp = params.k1 * x1 + params.k3 * x1^3  # Nonlinear

    # Soft-stopper spring component
    if abs(x1) < params.gp
        Fss = 0  # No soft-stopper force
    else
        Fss = params.kss * (abs(x1) - params.gss) * sign(x1)  # soft-stopper engaged
    end

    # Total spring force
    Fs = Fsp + Fss

    return Fs
end

# Collision force, Fc = Fcc + Fnc
function collision(x1, x2, params)
    # YOUR CODE HERE: Implement collision detection and forces
    # Return effective mass and force
end

# Damping force, Fd
function damping(x2, x2dot, params)
    # YOUR CODE HERE: Implement damping force equations
    # Should handle both pre-collision and post-collision regimes
end

# Variable capacitance, Cvar
function varcapacitance(x2, params)
    # YOUR CODE HERE: Implement capacitance calculation
    # Should handle both pre-collision and post-collision regimes
end

# Total capacitance, Ctot
function totcapacitance(x2, params)
    # YOUR CODE HERE: Implement capacitance calculation
    # Should handle both pre-collision and post-collision regimes
end

# Electrostatic force, Fe
function electrostatic(x2, params)
    # YOUR CODE HERE: Implement capacitance calculation
    # Should handle both pre-collision and post-collision regimes
end

# ------------------------- Model Dynamics -------------------------
function CoupledDynamics!(du, u, p, t)
    params, external_force = p
    
    # Unpack state variables
    u1, u2, u3, u4, u5, u6 = u
    
    # Get external force
    Fext = external_force(t)
    
    # Calculate forces
    Fs = spring(u1, params)
    Fc = collision(u1, u3, params)
    Fd = damping(u3, u4, params)
    Cvar = varcapacitance(u2, params)
    Ctot = totcapacitance(u2, params)
    Fe = electrostatic(u1, u3, u5, params)
    
    # State derivatives
    du[1] = u2              
    du[2] = (Fs + (params.N / 2) * Fc) / params.m1 - Fext               
    du[3] = u4               
    du[4] = (-Fc + Fd + Fe) / params.m2 - Fext               
    du[5] = (params.Vbias - (u5 / Ctotal)) / params.Rload               
    du[6] = (params.Vbias - u5 / Ctotal - u6) / (params.Rload * Ctotal)             
end

# ------------------------- External Force -------------------------
# External force parameters
f = 20.0  # Frequency (Hz)
alpha = 3.0  # Applied acceleration constant
g = 9.81  # Gravitational constant 
A = alpha * g  # Sinusoidal amplitude
t_ramp = 0.2  # Ramp-up duration (s)
noise_amplitude = 1.0 * g  # White noise amplitude 
noise_seed = 12345 # Seed for reproducible noise

# Sinusoidal force function
function sine_force(t; A, f, use_ramp, t_ramp)
    if use_ramp
        ramp = t < t_ramp ? t / t_ramp : 1.0
        return A * ramp * sin(2π * f * t)
    else
        return A * sin(2π * f * t)
    end
end

# White noise force function
function white_noise_force(t; amplitude, seed, dt = 0.001)
    # Initialize random number generator with seed
    rng = Random.MersenneTwister(seed)
    
    # Discretize time to generate consistent noise
    discrete_t = floor(Int, t/dt)
    
    # Reset RNG to seed and advance to the current time point
    Random.seed!(rng, seed)
    for _ in 1:discrete_t
        rand(rng)
    end
    
    # Generate white noise
    return amplitude * (2*rand(rng) - 1)
end

# No force function
zero_force(t) = 0.0

# Force type selection (1: sine, 2: white noise, 3: none)
force_type = 1
use_ramp = true
function create_external_force(force_type)
    if force_type == 1  # Sinusoidal
        return t -> sine_force(t, A = A, f = f, use_ramp = use_ramp, t_ramp = t_ramp)
    elseif force_type == 2  # White noise
        return t -> white_noise_force(t, amplitude = noise_amplitude, seed = noise_seed)
    else  # No external force
        return t -> zero_force(t)
    end
end

# Create the external force function
Fext_input = create_external_force(force_type)

# ------------------------- Initial Conditions -------------------------
# Initial conditions
x1_0 = 0.0  # Initial displacement, shuttle
v1_0 = 0.0  # Initial velocity, shuttle
x2_0 = 0.0  # Initial displacement, electrode
v2_0 = 0.0  # Initial velocity, electrode

# Compute initial capacitance and charge
Cinitial = capacitance(x1_0, x2_0, 0.0, params)
q_0 = params.Vbias * (Cinitial + params.cp)  # Initial charge
V_0 = params.Vbias - (q_0 / (Cinitial + params.cp))  # Initial voltage
u0 = [x1_0, v1_0, x2_0, v2_0, q_0, V_0]

# ------------------------- Solve Analytical Model -------------------------
# Simulation parameters
tspan = (0.0, 0.5) # simulation length
# teval = () # evaluation steps
abstol = 1e-9 # absolute solver tolerance
reltol = 1e-6 # relative solver tolerance

# Define and solve the ODE problem
eqn = ODEProblem(CoupledDynamics!, z0, tspan, p_new)
# eqn = ODEProblem(AnalyticalModel.CoupledSystem!, z0, tspan, p_new)

# Solve the system using Rosenbrock23 solver
sol = solve(eqn, Rosenbrock23(); abstol=abstol, reltol=reltol, maxiters=1e7)
# If the system is too stiff, use CVODE_BDF from Sundials
# sol = solve(eqn, CVODE_BDF(), abstol=abstol, reltol=reltol, maxiters=Int(1e9))

# Verify the solution structure
println("Type of sol.u: ", typeof(sol.u))
println("Size of sol.u: ", size(sol.u))
println("Solver status: ", sol.retcode)

