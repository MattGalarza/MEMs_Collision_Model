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
        println("Tilt angle, a:  ", params.a)
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
    Fsp = -params.k1 * x1  # Linear
    # Fsp = -params.k1 * x1 - params.k3 * x1^3  # Nonlinear

    # Soft-stopper spring component
    if abs(x1) < params.gp
        Fss = 0  # No soft-stopper force
    else
        Fss = -params.kss * (abs(x1) - params.gss) * sign(x1)  # soft-stopper engaged
    end

    # Total spring force
    Fs = Fsp + Fss

    return Fs
end

# Collision force, Fc = Fcc + Fnc
function collision(x1, x2, params)
    if abs(x2) < params.gp # Independent electrode restoring force
        Fc = -params.ke * (x2 - x1)
    else # Electrode restoring force during contact
        Fc = params.ke * (abs(x1) - params.gp) * sign(x2)
    end
    return Fc
end

# Damping force, Fd
function damping(x2, x2dot, params)
    # Dynamic states
    va = x2dot  # Translational velocity of moving electrode
    vb = 0.0  # Translational velocity of fixed electrode (usually 0)

    if abs(x2) < params.gp  # translation (non-collision)
        # RHS damping force
        h = params.gp - x2  # RHS gap expression
        term1 = 12 * params.Tf * params.eta * (va - vb)
        term2 = (2 * params.Leff * params.a + (2 * h + params.Leff * params.a) * log(h / (h + params.Leff * params.a)))
        term3 = params.a^3 * (2 * h + params.Leff * params.a)
        Fr = (term1 * term2) / term3
        
        # LHS damping force
        h = params.gp + x2  # LHS gap expression
        term1 = 12 * params.Tf * params.eta * (va - vb)
        term2 = (2 * params.Leff * params.a + (2 * h + params.Leff * params.a) * log(h / (h + params.Leff * params.a)))
        term3 = params.a^3 * (2 * h + params.Leff * params.a)
        Fl = (term1 * term2) / term3   
        
        # Total force
        Fd = params.c * (Fr + Fl)
    else # rotation (collision)
    end
    return Fd
end

# Variable capacitance, Cvar
function capacitance(x2, params)
    # Parylene Layer Capacitance
    Crl = (params.e * params.ep * params.Leff * params.Tf) / params.Tp
    if abs(x2) < params.gp  # translation (non-collision)
        # RHS capacitance
        Cairr = ((params.e * params.Tf) / params.a) * log((params.gp - x2 + params.a * params.Leff) / (params.gp - x2))
        Cvarr = 1 / ((2 / Crl) + (1 / Cairr))
        # LHS capacitance
        Cairl = ((params.e * params.Tf) / params.a) * log((params.gp + x2 + params.a * params.Leff) / (params.gp + x2))
        Cvarl = 1 / ((2 / Crl) + (1 / Cairl))
        # Total variable capacitance
        Cvar = (params.N / 2) * (Cvarr + Cvarl)
    else  # rotation (collision)
        # RHS capacitance (collision)
        k = (2 * params.Tp) / (params.a * params.Leff)
        C0 = 7.105299639935359*10^-14
        Cinf =  9.95200974248769*10^-11
        Cairc = C0 + (Cinf - C0) * (log(1 + k * (x2 - params.gp)) / log(1 + k * (params.a * params.Leff)))
        Cvarc = 1 / ((2 / Crl) + (1 / Cairc))
        # LHS capacitance (non-collision)
        Cairnc = ((params.e * params.Tf) / params.a) * log((params.gp + x2 + params.a * params.Leff) / (params.gp + x2))
        Cvarnc = 1 / ((2 / Crl) + (1 / Cairnc))
        # Total variable capacitance
        Cvar = (params.N / 2) * (Cvarc + Cvarnc)
    end
    return Cvar        
end

# Electrostatic force, Fe
function electrostatic(x2, q, Cvar, params)
    # Calculate total capacitance
    Ctot = Cvar + params.Cp
    
    # Calculate capacitance derivative
    h = 1e-10  # Small step for numerical differentiation
    C1 = capacitance(x2 - h, params)
    C2 = capacitance(x2 + h, params)
    dCdx = (C2 - C1) / (2 * h)
    dC = ForwardDiff.derivative(capacitance, x2)
    
    # Calculate electrostatic force
    Fe = (-0.5 * (q^2 / Ctot^2) * dC) / (params.N / 2)
    return Ctot, Fe
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
    Cvar = capacitance(u2, params)
    Ctot, Fe = electrostatic(u3, u5, Cvar, params)
    
    # State derivatives
    du[1] = u2              
    du[2] = (Fs + (params.N / 2) * Fc) / params.m1 - Fext               
    du[3] = u4               
    du[4] = (-Fc + Fd + Fe) / params.m2 - Fext               
    du[5] = (params.Vbias - (u5 / Ctot)) / params.Rload               
    du[6] = (params.Vbias - u5 / Ctot - u6) / (params.Rload * Ctot)             
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

# ------------------------- Event Handling for System -------------------------
# Condition function to detect soft-stopper engagement, abs(x1) = gss
function stopper_condition(u, t, integrator)
    # Extract shuttle position
    x1 = u[1]
    params = integrator.p[1]
    # Root finding
    return abs(x1) - params.gss
end

# Condition function to detect electrode collision, abs(x2) = gp 
function collision_condition(u, t, integrator)
    # Extract electrode position
    x2 = u[3]
    params = integrator.p[1]
    # Root finding
    return abs(x2) - params.gp
end

# Affect function for soft-stopper engagement
function soft_stopper_affect!(integrator)
    ss_time = integrator.t
    ss_state = copy(integrator.u)
    
    # Store event information
    if !haskey(integrator.opts.callback.continuous_callbacks[2].affect!.kwargs, :ss_events)
        integrator.opts.callback.continuous_callbacks[2].affect!.kwargs[:ss_events] = [(ss_time, ss_state)]
    else
        push!(integrator.opts.callback.continuous_callbacks[2].affect!.kwargs[:ss_events], 
              (ss_time, ss_state))
    end
    
    # Tighten tolerances near discontinuity
    set_proposed_dt!(integrator, integrator.dt/2)
end

# Affect function for electrode collision
function collision_affect!(integrator)
    collision_time = integrator.t
    collision_state = copy(integrator.u)
    
    # Store event information
    if !haskey(integrator.opts.callback.continuous_callbacks[1].affect!.kwargs, :collision_events)
        integrator.opts.callback.continuous_callbacks[1].affect!.kwargs[:collision_events] = [(collision_time, collision_state)]
    else
        push!(integrator.opts.callback.continuous_callbacks[1].affect!.kwargs[:collision_events], 
              (collision_time, collision_state))
    end
    
    # Tighten tolerances near discontinuity
    set_proposed_dt!(integrator, integrator.dt/2)
end

# ------------------------- Initial Conditions -------------------------
# Initial conditions
x1_0 = 0.0  # Initial displacement, shuttle
v1_0 = 0.0  # Initial velocity, shuttle
x2_0 = 0.0  # Initial displacement, electrode
v2_0 = 0.0  # Initial velocity, electrode

# Compute initial capacitance and charge
Cinitial = capacitance(x2_0, params)
q_0 = params.Vbias * (Cinitial + params.cp)  # Initial charge
V_0 = params.Vbias - (q_0 / (Cinitial + params.cp))  # Initial voltage
u0 = [x1_0, v1_0, x2_0, v2_0, q_0, V_0]

# ------------------------- Solve Analytical Model -------------------------
# Simulation parameters
tspan = (0.0, 0.5) # simulation length
# teval = () # evaluation steps
abstol = 1e-9 # absolute solver tolerance
reltol = 1e-6 # relative solver tolerance

# Create callbacks for event detection
cb_collision = ContinuousCallback(collision_condition, collision_affect!, nothing; 
                                 rootfind=true, save_positions=(true,false))
cb_soft_stopper = ContinuousCallback(soft_stopper_condition, soft_stopper_affect!, nothing; 
                                    rootfind=true, save_positions=(true,false))
cb_set = CallbackSet(cb_collision, cb_soft_stopper)

# Problem parameters
problem_params = (params, Fext_input)

# Define the ODE problem with callbacks
prob = ODEProblem(CoupledDynamics!, u0, tspan, problem_params)

# Solve the system using an appropriate solver for stiff systems with discontinuities
sol = solve(prob, Rosenbrock23(), callback=cb_set, abstol=abstol, reltol=reltol, maxiters=1e7, dtmin=1e-15, force_dtmin=true, save_everystep=true, saveat=0.0:0.0001:0.5)

# If the system is too stiff, use CVODE_BDF from Sundials
# sol = solve(prob, CVODE_BDF(), callback=cb_set, abstol=abstol, reltol=reltol, maxiters=Int(1e9))

# Verify the solution structure
println("Type of sol.u: ", typeof(sol.u))
println("Size of sol.u: ", size(sol.u))
println("Solver status: ", sol.retcode)
println("Solution time span: ", sol.t[1], " to ", sol.t[end])
println("Number of timesteps: ", length(sol.t))

# Retrieve event information
collision_events = get(sol.prob.f.f.callbacks.continuous_callbacks[1].affect!.kwargs, :collision_events, [])
ss_events = get(sol.prob.f.f.callbacks.continuous_callbacks[2].affect!.kwargs, :ss_events, [])

println("\nCollision events: ", length(collision_events))
if !isempty(collision_events)
    for (i, (t, _)) in enumerate(collision_events)
        println("  Event $i at t = $t")
    end
end

println("\nSoft-stopper events: ", length(ss_events))
if !isempty(ss_events)
    for (i, (t, _)) in enumerate(ss_events)
        println("  Event $i at t = $t")
    end
end

# ------------------------- Model Results and Plotting -------------------------
# Extract solution components
times = sol.t
x1 = [u[1] for u in sol.u]
x1dot = [u[2] for u in sol.u]
x2 = [u[3] for u in sol.u]
x2dot = [u[4] for u in sol.u]
q = [u[5] for u in sol.u]
V = [u[6] for u in sol.u]

# Flag points where collision events occurred
ss_flags = [abs(x1_val) >= params.gss for x1_val in x1]
collision_flags = [abs(x2_val) >= params.gp for x2_val in x2]

# State components 
p1 = plot(sol.t, x1, xlabel = "Time (s)", ylabel = "x1 (m)", title = "Shuttle Mass Displacement (x1)")
display(p1)
p2 = plot(sol.t, x1dot, xlabel = "Time (s)", ylabel = "x1dot (m/s)", title = "Shuttle Mass Velocity (x1dot)")
display(p2)
p3 = plot(sol.t, x2, xlabel = "Time (s)", ylabel = "x2 (m)", title = "Mobile Electrode Displacement (x2)")
display(p3)
p4 = plot(sol.t, x2dot, xlabel = "Time (s)", ylabel = "x2dot (m/s)", title = "Mobile Electrode Velocity (x2)")
display(p4)
p5 = plot(sol.t, q, xlabel = "Time (s)", ylabel = "q (C)", title = "Charge (q)")
display(p5)
p6 = plot(sol.t, V, xlabel = "Time (s)", ylabel = "V (V)", title = "Output Voltage (V)")
display(p6)



