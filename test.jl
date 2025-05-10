# ------------------------- Libraries -------------------------
using Sundials, DifferentialEquations, Plots
using Parameters, ForwardDiff, SpecialFunctions, NaNMath, Random

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
    e::T = 8.85e-12  # Permittivity of free space
    ep::T = 3.2  # Relative permittivity of Parylene-C
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
    if abs(x1) < params.gss
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
        collision_side = sign(x2)
        boundary_position = collision_side * params.gp
        Fc = params.ke * (boundary_position - x1)
    end
    return Fc
end

# Damping force, Fd
function damping(x2, x2dot, params)
    # Dynamic states
    va = x2dot  # Translational velocity of mobile electrode
    vb = 0.0  # Translational velocity of anchored electrode (usually 0)

    if abs(x2) < params.gp  # translation (non-collision)
        # RHS damping force
        h = params.gp - x2  # RHS gap expression
        term1 = 12 * params.Tf * params.eta * (va - vb)
        term2 = (2 * params.Leff * params.a + (2 * h + params.Leff * params.a) * NaNMath.log(h / (h + params.Leff * params.a)))
        term3 = params.a^3 * (2 * h + params.Leff * params.a)
        Fr = (term1 * term2) / term3
        
        # LHS damping force
        h = params.gp + x2  # LHS gap expression
        term1 = 12 * params.Tf * params.eta * (va - vb)
        term2 = (2 * params.Leff * params.a + (2 * h + params.Leff * params.a) * NaNMath.log(h / (h + params.Leff * params.a)))
        term3 = params.a^3 * (2 * h + params.Leff * params.a)
        Fl = (term1 * term2) / term3   
        
        # Total force
        Fd = params.c * (Fr + Fl)
    else # rotation (collision)
        # Total force
        Fd = 0.0
    end
    return Fd
end

# Variable capacitance, Cvar
function capacitance(x2, params)
    # Parylene Layer Capacitance
    Crl = (params.e * params.ep * params.Leff * params.Tf) / params.Tp
    if abs(x2) < params.gp  # translation (non-collision)
        # RHS capacitance
        Cairr = ((params.e * params.Tf) / params.a) * NaNMath.log((params.gp - x2 + params.a * params.Leff) / (params.gp - x2))
        Cvarr = 1 / ((2 / Crl) + (1 / Cairr))
        # LHS capacitance
        Cairl = ((params.e * params.Tf) / params.a) * NaNMath.log((params.gp + x2 + params.a * params.Leff) / (params.gp + x2))
        Cvarl = 1 / ((2 / Crl) + (1 / Cairl))
        # Total variable capacitance
        Cvar = (params.N / 2) * (Cvarr + Cvarl)
    else  # rotation (collision)
        # RHS capacitance (collision)
        k = (2 * params.Tp) / (params.a * params.Leff)
        C0 = 7.105299639935359*10^-14
        Cinf =  9.95200974248769*10^-11
        Cairc = C0 + (Cinf - C0) * (NaNMath.log(1 + k * (abs(x2) - params.gp)) / NaNMath.log(1 + k * (params.a * params.Leff)))
        Cvarc = 1 / ((2 / Crl) + (1 / Cairc))
        # LHS capacitance (non-collision)
        Cairnc = ((params.e * params.Tf) / params.a) * NaNMath.log((params.gp + abs(x2) + params.a * params.Leff) / (params.gp + abs(x2)))
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
    
    # Calculate capacitance derivative using ForwardDiff
    dC = ForwardDiff.derivative(x -> capacitance(x, params), x2)
    
    # Calculate electrostatic force
    Fe = (-0.5 * (q^2 / Ctot^2) * dC) / (params.N / 2)
    return Ctot, Fe
end

# ------------------------- Filippov System Implementation -------------------------
# Define a regularized Filippov dynamics function that smoothly transitions between regimes
function CoupledDynamicsFilippov!(du, u, p, t)
    params, external_force = p
    
    # Unpack state variables
    u1, u2, u3, u4, u5, u6 = u
    
    # Get external force
    Fext = external_force(t)
    
    # Calculate the switching function (boundary condition)
    phi = abs(u3) - params.gp  # Positive in collision regime, negative in non-collision
    
    # Regularization parameter (smaller = closer to ideal Filippov system)
    epsilon = 1e-9  # 1 nanometer smoothing region
    
    # Compute the smooth transition function (Filippov's convex combination parameter)
    # alpha = 0 in non-collision regime, alpha = 1 in collision regime
    alpha = 0.5 * (1 + tanh(phi/epsilon))
    
    # Calculate forces in non-collision regime (alpha ≈ 0)
    # Spring force (same in both regimes)
    Fs = spring(u1, params)
    
    # Non-collision regime forces
    Fc_nc = -params.ke * (u3 - u1)  # Collision force pre-contact
    Fd_nc = damping(u3, u4, params)  # Damping force pre-contact
    Cvar_nc = capacitance(u3, params)  # Capacitance pre-contact
    Ctot_nc, Fe_nc = electrostatic(u3, u5, Cvar_nc, params)  # Electrostatic force pre-contact
    
    # Collision regime forces
    collision_side = sign(u3)
    boundary_position = collision_side * params.gp
    Fc_c = params.ke * (boundary_position - u1)  # Collision force during contact
    
    # In collision regime, damping is higher (to prevent bouncing)
    c_contact = 0.1 * sqrt(params.ke * params.m2)  # Critical damping coefficient
    Fd_c = -c_contact * u4  # Damping force during contact
    
    # Use the same electrostatic model in both regimes
    Cvar_c = Cvar_nc
    Ctot_c = Ctot_nc
    Fe_c = Fe_nc
    
    # Filippov convex combination of forces
    Fc = (1-alpha) * Fc_nc + alpha * Fc_c  # Blend collision forces
    Fd = (1-alpha) * Fd_nc + alpha * Fd_c  # Blend damping forces
    Ctot = Ctot_nc  # Capacitance doesn't change at discontinuity
    Fe = Fe_nc  # Electrostatic force doesn't change at discontinuity
    
    # State derivatives - a smooth transition between regimes
    du[1] = u2  # Shuttle velocity
    du[2] = (Fs + (params.N / 2) * Fc) / params.m1 - Fext  # Shuttle acceleration
    
    # For position and velocity of the electrode, we transition to constrained motion
    # in the collision regime:
    
    # In full collision (alpha ≈ 1), the electrode position is constrained to the boundary
    # and velocity approaches zero
    du[3] = (1-alpha) * u4  # Electrode velocity (approaches 0 in collision)
    
    # Electrode acceleration - continuous blend between regimes
    du[4] = (1-alpha) * ((-Fc + Fd + Fe) / params.m2 - Fext) - 
            alpha * (u4 / epsilon)  # Add strong damping in collision to drive velocity to zero
            
    # Electrical dynamics (same in both regimes)
    du[5] = (params.Vbias - (u5 / Ctot)) / params.Rload  # Charge dynamics
    du[6] = (params.Vbias - u5 / Ctot - u6) / (params.Rload * Ctot)  # Voltage dynamics
end

# ------------------------- External Force -------------------------
# External force parameters
f = 20.0  # Frequency (Hz)
alpha = 3.0  # Applied acceleration constant (increased from 2.0 to ensure collisions occur)
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
Cinitial = capacitance(x2_0, params)
q_0 = params.Vbias * (Cinitial + params.Cp)  # Initial charge
V_0 = params.Vbias - (q_0 / (Cinitial + params.Cp))  # Initial voltage
u0 = [x1_0, v1_0, x2_0, v2_0, q_0, V_0]

# ------------------------- Solve Using Filippov Approach -------------------------
# Simulation parameters
tspan = (0.0, 0.5)  # simulation length
abstol = 1e-8  # absolute solver tolerance
reltol = 1e-6  # relative solver tolerance

# Problem parameters
problem_params = (params, Fext_input)

# Define the ODE problem with Filippov dynamics
prob = ODEProblem(CoupledDynamicsFilippov!, u0, tspan, problem_params)

# Solve the system using an appropriate solver for stiff systems
sol = solve(prob, CVODE_BDF(), abstol=abstol, reltol=reltol, 
            maxiters=1e7, dtmin=1e-12, force_dtmin=true, 
            save_everystep=false, saveat=0.0:0.001:0.5)

# ------------------------- Solution Analysis -------------------------
# Verify the solution structure
println("Type of sol.u: ", typeof(sol.u))
println("Size of sol.u: ", size(sol.u))
println("Solver status: ", sol.retcode)
println("Solution time span: ", sol.t[1], " to ", sol.t[end])
println("Number of timesteps: ", length(sol.t))

# Extract solution components
times = sol.t
x1 = [u[1] for u in sol.u]
x1dot = [u[2] for u in sol.u]
x2 = [u[3] for u in sol.u]
x2dot = [u[4] for u in sol.u]
q = [u[5] for u in sol.u]
V = [u[6] for u in sol.u]

# Identify collision regions
collision_flags = [abs(x2_val) >= params.gp for x2_val in x2]
collision_indices = findall(collision_flags)
collision_times = times[collision_indices]
println("\nCollision regions detected: ", length(collision_indices) > 0 ? "Yes" : "No")
println("Number of timepoints in collision: ", length(collision_indices))

# ------------------------- Force Capture and Plotting -------------------------
# Initialize arrays to store forces and the blending parameter
Fs_array = Float64[] # Suspension spring force
Fc_array = Float64[] # Collision force
Fd_array = Float64[] # Damping force
Fe_array = Float64[] # Electrostatic force
Fext_array = Float64[] # External force
alpha_array = Float64[] # Filippov blending parameter

# Iterate over each solution point to compute forces
for (i, t) in enumerate(sol.t)
    # Extract state variables at time t
    u_state = sol.u[i]
    u1, u2, u3, u4, u5, u6 = u_state
    
    # Calculate the blending parameter
    phi = abs(u3) - params.gp
    epsilon = 1e-9
    alpha = 0.5 * (1 + tanh(phi/epsilon))
    push!(alpha_array, alpha)
    
    # Compute forces using the model functions
    # Fs - Spring force
    Fs = spring(u1, params)
    push!(Fs_array, Fs)
    
    # Fc - Collision force (blended)
    Fc_nc = -params.ke * (u3 - u1)
    collision_side = sign(u3)
    boundary_position = collision_side * params.gp
    Fc_c = params.ke * (boundary_position - u1)
    Fc = (1-alpha) * Fc_nc + alpha * Fc_c
    push!(Fc_array, Fc)
    
    # Fd - Damping force (blended)
    Fd_nc = damping(u3, u4, params)
    c_contact = 0.1 * sqrt(params.ke * params.m2)
    Fd_c = -c_contact * u4
    Fd = (1-alpha) * Fd_nc + alpha * Fd_c
    push!(Fd_array, Fd)
    
    # Fe - Electrostatic force
    Cvar = capacitance(u3, params)
    Ctot, Fe = electrostatic(u3, u5, Cvar, params)
    push!(Fe_array, Fe)
    
    # Fext - External force
    Fext = Fext_input(t)
    push!(Fext_array, Fext)
end

# State components 
p1 = plot(sol.t, x1, xlabel = "Time (s)", ylabel = "x1 (m)", title = "Shuttle Mass Displacement (x1)")
display(p1)
p2 = plot(sol.t, x1dot, xlabel = "Time (s)", ylabel = "x1dot (m/s)", title = "Shuttle Mass Velocity (x1dot)")
display(p2)
p3 = plot(sol.t, x2, xlabel = "Time (s)", ylabel = "x2 (m)", title = "Mobile Electrode Displacement (x2)")
display(p3)
p4 = plot(sol.t, x2dot, xlabel = "Time (s)", ylabel = "x2dot (m/s)", title = "Mobile Electrode Velocity (x2dot)")
display(p4)
p5 = plot(sol.t, q, xlabel = "Time (s)", ylabel = "q (C)", title = "Charge (q)")
display(p5)
p6 = plot(sol.t, V, xlabel = "Time (s)", ylabel = "V (V)", title = "Output Voltage (V)")
display(p6)

# Plotting the forces
p7 = plot(sol.t, Fs_array, xlabel = "Time (s)", ylabel = "Fs (N)", title = "Suspension Spring Force")
display(p7)
p8 = plot(sol.t, Fc_array, xlabel = "Time (s)", ylabel = "Fc (N)", title = "Collision Force")
display(p8)
p9 = plot(sol.t, Fd_array, xlabel = "Time (s)", ylabel = "Fd (N)", title = "Damping Force")
display(p9)
p10 = plot(sol.t, Fe_array, xlabel = "Time (s)", ylabel = "Fe (N)", title = "Electrostatic Force")
display(p10)
p11 = plot(sol.t, Fext_array, xlabel = "Time (s)", ylabel = "Fext (N)", title = "External Force")
display(p11)

# Plot the Filippov blending parameter
p12 = plot(sol.t, alpha_array, xlabel = "Time (s)", ylabel = "alpha", title = "Filippov Blending Parameter (alpha)",
          yticks=0:0.2:1, linewidth=2)
display(p12)

# Plot collision status for easier visualization
p13 = plot(sol.t, collision_flags, xlabel = "Time (s)", ylabel = "In Collision", title = "Collision Status", 
          seriestype=:step, yticks=[0, 1], linewidth=2)
display(p13)

# Add combined plots for additional insights
# Plot displacement with collision boundary
p14 = plot(sol.t, x2, label="Electrode Position (x2)", linewidth=2)
plot!(sol.t, fill(params.gp, length(sol.t)), label="Collision Boundary (+gp)", linestyle=:dash, linewidth=1)
plot!(sol.t, fill(-params.gp, length(sol.t)), label="Collision Boundary (-gp)", linestyle=:dash, linewidth=1)
xlabel!("Time (s)")
ylabel!("Position (m)")
title!("Electrode Position with Collision Boundaries")
display(p14)

# Plot total forces on electrode
total_force_electrode = [(-Fc_array[i] + Fd_array[i] + Fe_array[i] - Fext_array[i]) for i in 1:length(sol.t)]
p15 = plot(sol.t, total_force_electrode, xlabel = "Time (s)", ylabel = "Force (N)", 
          title = "Total Force on Electrode", linewidth=2)
display(p15)

# Compare displacement difference and collision force
x1_minus_x2 = [x1[i] - x2[i] for i in 1:length(sol.t)]
p16 = plot(sol.t, x1_minus_x2, label="x1 - x2", xlabel = "Time (s)", ylabel = "Distance (m)", 
          title = "Displacement Difference (x1 - x2)")
display(p16)

# Phase plot of electrode (position vs velocity)
p17 = plot(x2, x2dot, xlabel = "x2 (m)", ylabel = "x2dot (m/s)", 
          title = "Phase Portrait of Electrode", seriestype=:scatter,
          markersize=2, markerstrokewidth=0, alpha=0.5)
display(p17)

# Scatter plot of collision force vs displacement difference
p18 = plot(x1_minus_x2, Fc_array, seriestype=:scatter, markersize=2, markerstrokewidth=0,
          xlabel = "x1 - x2 (m)", ylabel = "Fc (N)", title = "Collision Force vs. Displacement Difference",
          alpha=0.5)
display(p18)