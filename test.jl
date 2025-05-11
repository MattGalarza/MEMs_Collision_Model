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
    Fs = Fsp + 0.01*Fss

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

# ------------------------- Elastic Collision Dynamics -------------------------
# Improved model with proper elastic collisions
function CoupledDynamicsFilippov!(du, u, p, t)
    params, external_force = p
    
    # Unpack state variables
    u1, u2, u3, u4, u5, u6 = u
    
    # Get external force
    Fext = external_force(t)
    
    # Check if we're at a boundary
    at_boundary = abs(u3) >= params.gp
    
    # Spring force calculation (common to both regimes)
    Fs = spring(u1, params)
    
    # Calculate forces based on position
    if !at_boundary  # Away from boundary - normal dynamics
        # Standard forces
        Fc = 0.0  # Connection spring between shuttle and electrode
        
        try
            Fd = damping(u3, u4, params)  # Normal damping
        catch
            Fd = -params.c * u4  # Simple fallback damping if calculation fails
        end
        
        # No electrostatic force for mechanical testing
        Fe = 0.0
        Ctot = params.Cp
        
        # Normal dynamics
        du[1] = u2  # Shuttle velocity
        du[2] = (Fs + (params.N / 2) * Fc) / params.m1 - Fext  # Shuttle acceleration
        du[3] = u4  # Electrode velocity
        du[4] = (-Fc + Fd + Fe) / params.m2 - Fext  # Electrode acceleration
    else  # At boundary - handle elastic collision
        # Determine which boundary is hit
        boundary_side = sign(u3)
        boundary_pos = boundary_side * params.gp
        
        # Calculate forces
        Fc = 0.0  # Using standard spring connection
        Fd = 0.0  # No damping during instantaneous collision
        Fe = 0.0  # No electrostatic force
        Ctot = params.Cp
        
        # Check if electrode is moving toward or away from boundary
        moving_toward_boundary = (boundary_side * u4 > 0)
        
        # Normal dynamics for shuttle (unchanged)
        du[1] = u2
        du[2] = (Fs + (params.N / 2) * Fc) / params.m1 - Fext
        
        if moving_toward_boundary  # Electrode is moving toward the boundary
            # Apply elastic collision with coefficient of restitution
            # For perfectly elastic collision, cor = 1.0
            cor = 0.95  # 95% energy conservation (slight damping)
            
            # Reflection direction should be opposite of current velocity
            reflected_velocity = -boundary_side * abs(u4) * cor
            
            # Set derivatives to reflect velocity and move away from boundary
            du[3] = reflected_velocity
            
            # Additional acceleration to move away from boundary
            repulsion_accel = -boundary_side * 1e-6  # Small repulsion to ensure separation
            du[4] = (-Fc) / params.m2 - Fext + repulsion_accel
        else  # Already moving away from boundary
            # Normal dynamics
            du[3] = u4
            du[4] = (-Fc) / params.m2 - Fext
        end
    end
    
    # Electrical dynamics (inactive for mechanical testing)
    du[5] = 0.0
    du[6] = 0.0
end

# ------------------------- External Force -------------------------
# External force parameters
f = 20.0  # Frequency (Hz)
alpha = 2.5  # Applied acceleration constant
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

# ------------------------- Solve System -------------------------
# Simulation parameters
tspan = (0.0, 0.5) # simulation length
abstol = 1e-8  # absolute solver tolerance
reltol = 1e-6  # relative solver tolerance

# Problem parameters
problem_params = (params, Fext_input)

# Define the ODE problem
prob = ODEProblem(CoupledDynamicsFilippov!, u0, tspan, problem_params)

# Add adaptive time stepping refinement near threshold
function affect_precision!(integrator)
    x1 = integrator.u[1]
    
    # Check if near soft-stopper threshold
    if abs(abs(x1) - params.gss) < 1e-7  # Within 0.1 μm of threshold
        # Reduce step size for better resolution near threshold
        set_proposed_dt!(integrator, integrator.dt/2)
    end
end

# Create discrete callback for precision control
precision_cb = DiscreteCallback(
    (u, t, integrator) -> true,  # Always active (checked at each step)
    affect_precision!
)

# Combine with existing callbacks
cb_set = CallbackSet(collision_cb, precision_cb)

# Use the combined callback set in the solver
sol = solve(prob, Rodas5(), 
            callback=cb_set,
            abstol=1e-8, reltol=1e-6, 
            maxiters=Int(5e8),
            dtmin=1e-14,
            force_dtmin=true,
            save_everystep=false, 
            saveat=0.0:0.001:0.5)

# Alternative solver options if needed
# sol = solve(prob, TRBDF2(), abstol=abstol, reltol=reltol, maxiters=Int(1e7))
# sol = solve(prob, Rodas5(), abstol=abstol, reltol=reltol, maxiters=Int(1e7))

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

# Detect collision events (not collision states)
collision_events = Float64[]
for i in 2:length(sol.t)
    prev_x2 = abs(sol.u[i-1][3])
    curr_x2 = abs(sol.u[i][3])
    
    # Detect crossing of boundary in either direction
    if (prev_x2 < params.gp && curr_x2 >= params.gp) || 
       (prev_x2 >= params.gp && curr_x2 < params.gp)
        push!(collision_events, sol.t[i])
    end
end

println("\nNumber of collision events: ", length(collision_events))
if !isempty(collision_events)
    for (i, t) in enumerate(collision_events)
        println("  Event $i at t = $t")
    end
end

# Also track when at boundary (for plotting)
at_boundary = [abs(x2_val) >= params.gp for x2_val in x2]
boundary_indices = findall(at_boundary)

# ------------------------- Force Calculation and Plotting -------------------------
# Initialize arrays to store forces
Fs_array = Float64[] # Suspension spring force
Fc_array = Float64[] # Collision force
Fd_array = Float64[] # Damping force
Fe_array = Float64[] # Electrostatic force
Fext_array = Float64[] # External force

# Iterate over each solution point to compute forces
for (i, t) in enumerate(sol.t)
    # Extract state variables at time t
    u_state = sol.u[i]
    u1, u2, u3, u4, u5, u6 = u_state
    
    # Compute forces
    # Fs - Spring force
    Fs = spring(u1, params)
    push!(Fs_array, Fs)
    
    # Fc - Collision force
    Fc = -params.ke * (u3 - u1)  # Use consistent expression for plotting
    push!(Fc_array, Fc)
    
    # Fd - Damping force
    if abs(u3) < params.gp
        try
            Fd = damping(u3, u4, params)
        catch
            Fd = -params.c * u4  # Simple fallback
        end
    else
        Fd = 0.0  # No damping at boundary
    end
    push!(Fd_array, Fd)
    
    # Fe - Electrostatic force (zero for mechanical test)
    push!(Fe_array, 0.0)
    
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

# Highlight collision events on position plot
p5 = plot(sol.t, x2, label="Electrode Position", linewidth=2)
# Add vertical lines at collision events
for t in collision_events
    plot!([t, t], [minimum(x2), maximum(x2)], linestyle=:dash, color=:red, label=nothing, linewidth=1)
end
# Add boundary lines
plot!(sol.t, fill(params.gp, length(sol.t)), label="Boundary (+gp)", linestyle=:dash, color=:green)
plot!(sol.t, fill(-params.gp, length(sol.t)), label="Boundary (-gp)", linestyle=:dash, color=:green)
xlabel!("Time (s)")
ylabel!("Position (m)")
title!("Electrode Position with Collision Events")
display(p5)

# Force plots
p6 = plot(sol.t, Fs_array, xlabel = "Time (s)", ylabel = "Fs (N)", title = "Suspension Spring Force")
display(p6)
p7 = plot(sol.t, Fc_array, xlabel = "Time (s)", ylabel = "Fc (N)", title = "Collision Force")
display(p7)
p8 = plot(sol.t, Fd_array, xlabel = "Time (s)", ylabel = "Fd (N)", title = "Damping Force")
display(p8)
p9 = plot(sol.t, Fe_array, xlabel = "Time (s)", ylabel = "Fe (N)", title = "Electrostatic Force (Zero)")
display(p9)
p10 = plot(sol.t, Fext_array, xlabel = "Time (s)", ylabel = "Fext (N)", title = "External Force")
display(p10)

# Combined forces plot
p11 = plot(sol.t, [Fs_array Fc_array Fe_array Fd_array], 
          label=["Spring" "Collision" "Electrostatic" "Damping"],
          xlabel = "Time (s)", ylabel = "Force (N)", title = "Force Comparison")
display(p11)

# Phase portrait of electrode
p12 = plot(x2, x2dot, xlabel = "x2 (m)", ylabel = "x2dot (m/s)", 
          title = "Phase Portrait of Electrode", seriestype=:scatter,
          markersize=2, markerstrokewidth=0, alpha=0.5)
display(p12)

# Energy analysis
# Calculate kinetic and potential energies
KE_shuttle = [0.5 * params.m1 * v^2 for v in x1dot]
KE_electrode = [0.5 * params.m2 * v^2 for v in x2dot]
PE_spring = [0.5 * params.k1 * x^2 for x in x1]
PE_conn = [0.5 * params.ke * (x1[i] - x2[i])^2 for i in 1:length(x1)]
total_energy = KE_shuttle + KE_electrode + PE_spring + PE_conn

# Energy plot
p13 = plot(sol.t, [KE_shuttle KE_electrode PE_spring PE_conn total_energy], 
          label=["KE Shuttle" "KE Electrode" "PE Spring" "PE Connection" "Total Energy"],
          xlabel = "Time (s)", ylabel = "Energy (J)", title = "Energy Components")
display(p13)

# Plot the sign of x1 vs time
p_sign = plot(sol.t, sign.(x1), ylabel="sign(x1)", xlabel="Time (s)")
display(p_sign)

# Plot both the absolute value of x1 and the soft-stopper threshold
p_abs = plot(sol.t, abs.(x1), label="|x1|")
plot!(sol.t, fill(params.gss, length(sol.t)), label="gss", linestyle=:dash)
ylabel!("Position (m)")
xlabel!("Time (s)")
display(p_abs)

# Plot the external force
p_ext = plot(sol.t, [Fext_input(t) for t in sol.t], ylabel="Fext (N)", xlabel="Time (s)")
display(p_ext)

# Plot soft-stopper force component separately
Fss_array = Float64[]
for x_val in x1
    if abs(x_val) < params.gss
        push!(Fss_array, 0.0)
    else
        push!(Fss_array, -params.kss * (abs(x_val) - params.gss) * sign(x_val))
    end
end
p_fss = plot(sol.t, Fss_array, ylabel="Fss (N)", xlabel="Time (s)")
display(p_fss)