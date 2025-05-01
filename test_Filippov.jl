using DifferentialEquations, Plots

# Parameters
m₁ = 1.0e-7  # kg, shuttle mass
m₂ = 5.0e-8  # kg, mobile electrode mass
k₁ = 1.0     # N/m, primary spring constant
kₛₛ = 10.0   # N/m, suspension spring constant
kₑ = 5.0     # N/m, collision spring constant
gₛₛ = 0.1e-6 # m, suspension spring gap
gₚ = 0.5e-6  # m, collision threshold
R = 1.0e6    # Ω, resistance
VDC = 5.0    # V, applied voltage
Fapplied = 0.0 # N, applied force
N = 2        # Number of electrodes

# Capacitance parameters
a = 0.2e-6   # m, distance parameter
ε = 8.85e-12 # F/m, permittivity of free space
εₚ = 1.0     # relative permittivity
Leff = 50e-6 # m, effective length
Tf = 2e-6    # m, finger thickness
Tp = 2e-6    # m, plate thickness
Cp = 1e-15   # F, parasitic capacitance
η = 1.8e-5   # Pa·s, air viscosity
c = 0.8      # damping coefficient

# Calculate suspension spring force
function calculate_suspension_force(x₁, k₁, kₛₛ, gₛₛ)
    Fₛₚ = -k₁ * x₁
    if abs(x₁) < gₛₛ
        Fₛₛ = 0.0
    else
        Fₛₛ = -kₛₛ * (abs(x₁) - gₛₛ) * sign(x₁)
    end
    return Fₛₚ + Fₛₛ
end

# Calculate collision force
function calculate_collision_force(x₁, x₂, kₑ, gₚ)
    if abs(x₂) < gₚ
        return -kₑ * (x₁ - x₂)
    else
        return -kₑ * (x₁ - x₂) + kₑ * (abs(x₂) - gₚ) * sign(x₂)
    end
end

# Calculate capacitance for pre-collision mode (|x₂| < gₚ)
function calculate_capacitance_pre_collision(x₂, a, ε, εₚ, gₚ, Leff, Tf, Tp, Cp)
    # Capacitance between rigid insulator and electrode
    Cᵣᵢ = (ε * εₚ * Leff * Tf) / Tp
    
    # Right side air capacitance
    Cₐᵢᵣ₋ᵣ = (ε * Tf / (2 * a)) * log((gₚ - x₂ + 2 * a * Leff) / (gₚ - x₂))
    
    # Left side air capacitance
    Cₐᵢᵣ₋ₗ = (ε * Tf / (2 * a)) * log((gₚ + x₂ + 2 * a * Leff) / (gₚ + x₂))
    
    # Variable capacitances
    Cᵥₐᵣ₋ᵣ = 1.0 / (2.0 / Cᵣᵢ + 1.0 / Cₐᵢᵣ₋ᵣ)
    Cᵥₐᵣ₋ₗ = 1.0 / (2.0 / Cᵣᵢ + 1.0 / Cₐᵢᵣ₋ₗ)
    
    # Total variable capacitance
    Cᵥₐᵣ = (N / 2) * (Cᵥₐᵣ₋ᵣ + Cᵥₐᵣ₋ₗ)
    
    return Cᵥₐᵣ
end

# Calculate capacitance for post-collision mode (|x₂| >= gₚ)
function calculate_capacitance_post_collision(x₂, a, ε, εₚ, gₚ, Leff, Tf, Tp, Cp)
    # Capacitance between rigid insulator and electrode
    Cᵣᵢ = (ε * εₚ * Leff * Tf) / Tp
    
    # Contact side air capacitance
    Cₐᵢᵣ₋ₖ = (ε * Tf * Leff / abs(x₂)) * log((2 * Tp + abs(x₂)) / (2 * Tp))
    
    # Non-contact side air capacitance (same as left side in pre-collision)
    Cₐᵢᵣ₋ₙₖ = (ε * Tf / (2 * a)) * log((gₚ + x₂ + 2 * a * Leff) / (gₚ + x₂))
    
    # Variable capacitances
    Cᵥₐᵣ₋ₖ = 1.0 / (2.0 / Cᵣᵢ + 1.0 / Cₐᵢᵣ₋ₖ)
    Cᵥₐᵣ₋ₙₖ = 1.0 / (2.0 / Cᵣᵢ + 1.0 / Cₐᵢᵣ₋ₙₖ)
    
    # Total variable capacitance
    Cᵥₐᵣ = (N / 2) * (Cᵥₐᵣ₋ₖ + Cᵥₐᵣ₋ₙₖ)
    
    return Cᵥₐᵣ
end

# Calculate electrostatic force
function calculate_electrostatic_force(x₂, q, Cᵥₐᵣ, Cp)
    Ctotal = Cᵥₐᵣ + Cp
    
    # Calculate derivative of Cᵥₐᵣ with respect to x₂ using finite difference
    δ = 1e-9
    Cᵥₐᵣ_plus = calculate_capacitance_pre_collision(x₂ + δ, a, ε, εₚ, gₚ, Leff, Tf, Tp, Cp)
    Cᵥₐᵣ_minus = calculate_capacitance_pre_collision(x₂ - δ, a, ε, εₚ, gₚ, Leff, Tf, Tp, Cp)
    dCᵥₐᵣ_dx₂ = (Cᵥₐᵣ_plus - Cᵥₐᵣ_minus) / (2 * δ)
    
    # Electrostatic force
    Fe = -(2 / N) * (q^2 / (2 * Ctotal^2)) * dCᵥₐᵣ_dx₂
    
    return Fe
end

# Calculate damping force for pre-collision (|x₂| < gₚ)
function calculate_damping_force_pre_collision(x₂, ẋ₂)
    # Left side
    A_L1 = ((gₚ + x₂) * Leff / (gₚ + x₂ + 2 * a)) * ẋ₂
    A_L2 = (12 * η * a^(-2) / (2 * (gₚ + x₂) + a * Leff)) * ẋ₂
    
    F_dL = Tf * (Leff * A_L2 - (6 * η * Leff / (a * (gₚ + x₂) * (gₚ + x₂ + a * Leff))) * A_L1 +
                6 * η * Leff * a^(-4) * (gₚ + x₂ + a * Leff) * ẋ₂ +
                12 * η * a^(-3) * ẋ₂ * log(abs((gₚ + x₂) / (gₚ + x₂ + a * Leff))))
    
    # Right side
    A_R1 = -((gₚ - x₂) * Leff / (gₚ - x₂ + 2 * a)) * ẋ₂
    A_R2 = (12 * η * a^(-2) / (2 * (gₚ - x₂) + a * Leff)) * ẋ₂
    
    F_dR = Tf * (Leff * A_R2 - (6 * η * Leff / (a * (gₚ - x₂) * (gₚ - x₂ + a * Leff))) * A_R1 +
                6 * η * Leff * a^(-4) * (gₚ - x₂ + a * Leff) * ẋ₂ +
                12 * η * a^(-3) * ẋ₂ * log(abs((gₚ - x₂) / (gₚ - x₂ + a * Leff))))
    
    return -c * (F_dL + F_dR)
end

# Calculate damping force for post-collision (|x₂| >= gₚ)
function calculate_damping_force_post_collision(x₂, ẋ₂)
    # Use pre-collision damping as an approximation with adjustment
    base_damping = calculate_damping_force_pre_collision(sign(x₂) * (2*gₚ - abs(x₂)), ẋ₂)
    
    # Add additional damping due to deformation
    deformation_factor = 1.5 * (abs(x₂) / gₚ - 1.0)
    
    return base_damping * (1.0 + deformation_factor)
end

#===============================
    Filippov Implementation
===============================
# This approach handles the discontinuity using Filippov's convex combination
# of vector fields when exactly on the switching surface
=#
function filippov_system!(du, u, p, t)
    # Unpack state variables: [x₁, ẋ₁, x₂, ẋ₂, q]
    z₁, z₂, z₃, z₄, z₅ = u
    
    # Define the switching condition: |x₂| - gₚ
    h = abs(z₃) - gₚ
    
    # Calculate suspension spring force
    Fs = calculate_suspension_force(z₁, k₁, kₛₛ, gₛₛ)
    
    # Define vector fields for both modes
    function pre_collision_dynamics(u)
        z₁, z₂, z₃, z₄, z₅ = u
        
        # Calculate forces for pre-collision mode
        Fc = -kₑ * (z₁ - z₃)
        Cᵥₐᵣ = calculate_capacitance_pre_collision(z₃, a, ε, εₚ, gₚ, Leff, Tf, Tp, Cp)
        Fe = calculate_electrostatic_force(z₃, z₅, Cᵥₐᵣ, Cp)
        Fd = calculate_damping_force_pre_collision(z₃, z₄)
        Ctotal = Cᵥₐᵣ + Cp
        
        # State derivatives
        dz₁ = z₂
        dz₂ = (1/m₁) * (Fs + N*Fc - m₁*Fapplied)
        dz₃ = z₄
        dz₄ = (1/m₂) * (Fc + Fd + Fe - m₂*Fapplied)
        dz₅ = (1/R) * (VDC - z₅/Ctotal)
        
        return [dz₁, dz₂, dz₃, dz₄, dz₅]
    end
    
    function post_collision_dynamics(u)
        z₁, z₂, z₃, z₄, z₅ = u
        
        # Calculate forces for post-collision mode
        Fc = -kₑ * (z₁ - z₃) + kₑ * (abs(z₃) - gₚ) * sign(z₃)
        Cᵥₐᵣ = calculate_capacitance_post_collision(z₃, a, ε, εₚ, gₚ, Leff, Tf, Tp, Cp)
        Fe = calculate_electrostatic_force(z₃, z₅, Cᵥₐᵣ, Cp)
        Fd = calculate_damping_force_post_collision(z₃, z₄)
        Ctotal = Cᵥₐᵣ + Cp
        
        # State derivatives with modified mass (rotation mode)
        dz₁ = z₂
        dz₂ = (1/m₁) * (Fs + N*Fc - m₁*Fapplied)
        dz₃ = z₄
        dz₄ = (1/(2*m₂)) * (Fc + Fd + Fe - m₂*Fapplied)  # Note the modified mass
        dz₅ = (1/R) * (VDC - z₅/Ctotal)
        
        return [dz₁, dz₂, dz₃, dz₄, dz₅]
    end
    
    # Filippov's method for handling the discontinuity
    if abs(h) < 1e-10  # If almost exactly on the switching surface
        # Calculate normal vector to switching surface
        grad_h = [0.0, 0.0, sign(z₃), 0.0, 0.0]
        
        # Calculate vector fields at the current point
        f1 = pre_collision_dynamics(u)
        f2 = post_collision_dynamics(u)
        
        # Calculate Lie derivatives (rate of change of h along f1 and f2)
        Lfh1 = dot(grad_h, f1)
        Lfh2 = dot(grad_h, f2)
        
        # Check for sliding mode
        if Lfh1 > 0 && Lfh2 < 0
            # Use convex combination for sliding mode
            α = Lfh2 / (Lfh2 - Lfh1)
            du .= α * f1 + (1-α) * f2
            
            # Store information about sliding mode for visualization
            # (In a full implementation, you'd use a global variable or integrator cache)
            # sliding_mode[end+1] = (t, true)
        else
            # Not a sliding mode, use the appropriate vector field
            if Lfh1 >= 0
                du .= f2  # Moving towards post-collision
            else
                du .= f1  # Moving towards pre-collision
            end
        end
    elseif h < 0
        # Pre-collision dynamics
        du .= pre_collision_dynamics(u)
    else
        # Post-collision dynamics
        du .= post_collision_dynamics(u)
    end
end

#===============================
    Hybrid Automaton Implementation
===============================
# This approach explicitly models different modes and transitions
=#

# Pre-collision dynamics
function pre_collision!(du, u, p, t)
    # Unpack state variables: [x₁, ẋ₁, x₂, ẋ₂, q]
    z₁, z₂, z₃, z₄, z₅ = u
    
    # Calculate forces for pre-collision mode
    Fs = calculate_suspension_force(z₁, k₁, kₛₛ, gₛₛ)
    Fc = -kₑ * (z₁ - z₃)
    Cᵥₐᵣ = calculate_capacitance_pre_collision(z₃, a, ε, εₚ, gₚ, Leff, Tf, Tp, Cp)
    Fe = calculate_electrostatic_force(z₃, z₅, Cᵥₐᵣ, Cp)
    Fd = calculate_damping_force_pre_collision(z₃, z₄)
    Ctotal = Cᵥₐᵣ + Cp
    
    # State derivatives
    du[1] = z₂
    du[2] = (1/m₁) * (Fs + N*Fc - m₁*Fapplied)
    du[3] = z₄
    du[4] = (1/m₂) * (Fc + Fd + Fe - m₂*Fapplied)
    du[5] = (1/R) * (VDC - z₅/Ctotal)
end

# Post-collision dynamics
function post_collision!(du, u, p, t)
    # Unpack state variables: [x₁, ẋ₁, x₂, ẋ₂, q]
    z₁, z₂, z₃, z₄, z₅ = u
    
    # Calculate forces for post-collision mode
    Fs = calculate_suspension_force(z₁, k₁, kₛₛ, gₛₛ)
    Fc = -kₑ * (z₁ - z₃) + kₑ * (abs(z₃) - gₚ) * sign(z₃)
    Cᵥₐᵣ = calculate_capacitance_post_collision(z₃, a, ε, εₚ, gₚ, Leff, Tf, Tp, Cp)
    Fe = calculate_electrostatic_force(z₃, z₅, Cᵥₐᵣ, Cp)
    Fd = calculate_damping_force_post_collision(z₃, z₄)
    Ctotal = Cᵥₐᵣ + Cp
    
    # State derivatives with modified mass
    du[1] = z₂
    du[2] = (1/m₁) * (Fs + N*Fc - m₁*Fapplied)
    du[3] = z₄
    du[4] = (1/(2*m₂)) * (Fc + Fd + Fe - m₂*Fapplied)  # Modified mass
    du[5] = (1/R) * (VDC - z₅/Ctotal)
end

# Event detection for transition to post-collision
function to_post_collision(u, t, integrator)
    z₃ = u[3]  # x₂
    return abs(z₃) - gₚ  # Becomes positive when |x₂| >= gₚ
end

# Event detection for transition to pre-collision
function to_pre_collision(u, t, integrator)
    z₃ = u[3]  # x₂
    return gₚ - abs(z₃)  # Becomes positive when |x₂| < gₚ
end

# Event handlers
function enter_post_collision!(integrator)
    integrator.f = post_collision!
    
    # For visualization, store the transition time (not shown in this simple example)
    # push!(transition_times, integrator.t)
end

function enter_pre_collision!(integrator)
    integrator.f = pre_collision!
    
    # For visualization, store the transition time (not shown in this simple example)
    # push!(transition_times, integrator.t)
end

#===============================
    Simulation and Visualization
===============================
# This section simulates both models and creates plots to visualize the states
=#

# Initial conditions
x₁₀ = 0.0       # Initial shuttle position
ẋ₁₀ = 0.0       # Initial shuttle velocity
x₂₀ = 0.2e-6    # Initial mobile electrode position (within pre-collision region)
ẋ₂₀ = 0.0       # Initial mobile electrode velocity
q₀ = 1.0e-12    # Initial charge

u₀ = [x₁₀, ẋ₁₀, x₂₀, ẋ₂₀, q₀]
tspan = (0.0, 2.0e-3)  # 2 ms simulation

# Solve using Filippov approach
prob_filippov = ODEProblem(filippov_system!, u₀, tspan)
sol_filippov = solve(prob_filippov, Tsit5(), reltol=1e-8, abstol=1e-10)

# Solve using Hybrid Automaton approach
# Create callbacks for mode transitions
cb_to_post = ContinuousCallback(to_post_collision, enter_post_collision!)
cb_to_pre = ContinuousCallback(to_pre_collision, enter_pre_collision!)
cb_set = CallbackSet(cb_to_post, cb_to_pre)

# Determine initial mode based on initial conditions
initial_mode = abs(u₀[3]) < gₚ ? pre_collision! : post_collision!

prob_hybrid = ODEProblem(initial_mode, u₀, tspan)
sol_hybrid = solve(prob_hybrid, Tsit5(), callback=cb_set, reltol=1e-8, abstol=1e-10)

# Create a function to identify when the system is in post-collision mode
function is_post_collision(x₂)
    return abs(x₂) >= gₚ
end

# Create plots
# 1. Plot states over time for Filippov approach
filippov_plt = plot(layout=(3,2), size=(900, 600), title="Filippov Approach")

# Add vertical lines where mode transitions occur (Filippov)
mode_transitions_filippov = [sol_filippov.t[i] for i in 2:length(sol_filippov.t) 
                            if is_post_collision(sol_filippov[3, i-1]) != is_post_collision(sol_filippov[3, i])]

# Shuttle position (x₁)
plot!(filippov_plt[1], sol_filippov.t, sol_filippov[1, :] .* 1e6, 
      xlabel="Time (s)", ylabel="x₁ (μm)", label="Shuttle Position",
      subplot=1)
vline!(filippov_plt[1], mode_transitions_filippov, linestyle=:dash, color=:red, alpha=0.3, label="Mode Transition")

# Shuttle velocity (ẋ₁)
plot!(filippov_plt[2], sol_filippov.t, sol_filippov[2, :] .* 1e3, 
      xlabel="Time (s)", ylabel="ẋ₁ (mm/s)", label="Shuttle Velocity",
      subplot=2)
vline!(filippov_plt[2], mode_transitions_filippov, linestyle=:dash, color=:red, alpha=0.3, label=false)

# Mobile electrode position (x₂)
plot!(filippov_plt[3], sol_filippov.t, sol_filippov[3, :] .* 1e6, 
      xlabel="Time (s)", ylabel="x₂ (μm)", label="Mobile Electrode Position",
      subplot=3)
hline!(filippov_plt[3], [gₚ, -gₚ] .* 1e6, linestyle=:dash, color=:orange, label="Collision Threshold")
vline!(filippov_plt[3], mode_transitions_filippov, linestyle=:dash, color=:red, alpha=0.3, label=false)

# Mobile electrode velocity (ẋ₂)
plot!(filippov_plt[4], sol_filippov.t, sol_filippov[4, :] .* 1e3, 
      xlabel="Time (s)", ylabel="ẋ₂ (mm/s)", label="Mobile Electrode Velocity",
      subplot=4)
vline!(filippov_plt[4], mode_transitions_filippov, linestyle=:dash, color=:red, alpha=0.3, label=false)

# Charge (q)
plot!(filippov_plt[5], sol_filippov.t, sol_filippov[5, :] .* 1e12, 
      xlabel="Time (s)", ylabel="q (pC)", label="Charge",
      subplot=5)
vline!(filippov_plt[5], mode_transitions_filippov, linestyle=:dash, color=:red, alpha=0.3, label=false)

# Mode indicator
mode_indicator = map(x -> is_post_collision(x) ? 1.0 : 0.0, sol_filippov[3, :])
plot!(filippov_plt[6], sol_filippov.t, mode_indicator, 
      xlabel="Time (s)", ylabel="Mode", label="System Mode",
      ylim=(-0.1, 1.1), yticks=[0, 1], yformatter=y->["Pre-collision", "Post-collision"][Int(y)+1],
      subplot=6)

# 2. Plot states over time for Hybrid Automaton approach
hybrid_plt = plot(layout=(3,2), size=(900, 600), title="Hybrid Automaton Approach")

# Add vertical lines where mode transitions occur (Hybrid)
mode_transitions_hybrid = [sol_hybrid.t[i] for i in 2:length(sol_hybrid.t) 
                           if is_post_collision(sol_hybrid[3, i-1]) != is_post_collision(sol_hybrid[3, i])]

# Shuttle position (x₁)
plot!(hybrid_plt[1], sol_hybrid.t, sol_hybrid[1, :] .* 1e6, 
      xlabel="Time (s)", ylabel="x₁ (μm)", label="Shuttle Position",
      subplot=1)
vline!(hybrid_plt[1], mode_transitions_hybrid, linestyle=:dash, color=:red, alpha=0.3, label="Mode Transition")

# Shuttle velocity (ẋ₁)
plot!(hybrid_plt[2], sol_hybrid.t, sol_hybrid[2, :] .* 1e3, 
      xlabel="Time (s)", ylabel="ẋ₁ (mm/s)", label="Shuttle Velocity",
      subplot=2)
vline!(hybrid_plt[2], mode_transitions_hybrid, linestyle=:dash, color=:red, alpha=0.3, label=false)

# Mobile electrode position (x₂)
plot!(hybrid_plt[3], sol_hybrid.t, sol_hybrid[3, :] .* 1e6, 
      xlabel="Time (s)", ylabel="x₂ (μm)", label="Mobile Electrode Position",
      subplot=3)
hline!(hybrid_plt[3], [gₚ, -gₚ] .* 1e6, linestyle=:dash, color=:orange, label="Collision Threshold")
vline!(hybrid_plt[3], mode_transitions_hybrid, linestyle=:dash, color=:red, alpha=0.3, label=false)

# Mobile electrode velocity (ẋ₂)
plot!(hybrid_plt[4], sol_hybrid.t, sol_hybrid[4, :] .* 1e3, 
      xlabel="Time (s)", ylabel="ẋ₂ (mm/s)", label="Mobile Electrode Velocity",
      subplot=4)
vline!(hybrid_plt[4], mode_transitions_hybrid, linestyle=:dash, color=:red, alpha=0.3, label=false)

# Charge (q)
plot!(hybrid_plt[5], sol_hybrid.t, sol_hybrid[5, :] .* 1e12, 
      xlabel="Time (s)", ylabel="q (pC)", label="Charge",
      subplot=5)
vline!(hybrid_plt[5], mode_transitions_hybrid, linestyle=:dash, color=:red, alpha=0.3, label=false)

# Mode indicator
mode_indicator = map(x -> is_post_collision(x) ? 1.0 : 0.0, sol_hybrid[3, :])
plot!(hybrid_plt[6], sol_hybrid.t, mode_indicator, 
      xlabel="Time (s)", ylabel="Mode", label="System Mode",
      ylim=(-0.1, 1.1), yticks=[0, 1], yformatter=y->["Pre-collision", "Post-collision"][Int(y)+1],
      subplot=6)

# 3. Compare x₂ between the approaches
compare_plt = plot(size=(800, 400), title="Mobile Electrode Position Comparison")
plot!(compare_plt, sol_filippov.t, sol_filippov[3, :] .* 1e6, 
      label="Filippov Approach", xlabel="Time (s)", ylabel="x₂ (μm)")
plot!(compare_plt, sol_hybrid.t, sol_hybrid[3, :] .* 1e6, 
      label="Hybrid Automaton Approach", linestyle=:dash)
hline!(compare_plt, [gₚ, -gₚ] .* 1e6, linestyle=:dash, color=:orange, label="Collision Threshold")

# 4. Phase plane plot (x₂ vs ẋ₂)
phase_plt = plot(size=(600, 500), title="Phase Plane (x₂ vs ẋ₂)")
plot!(phase_plt, sol_filippov[3, :] .* 1e6, sol_filippov[4, :] .* 1e3, 
      label="Filippov Approach", xlabel="x₂ (μm)", ylabel="ẋ₂ (mm/s)")
plot!(phase_plt, sol_hybrid[3, :] .* 1e6, sol_hybrid[4, :] .* 1e3, 
      label="Hybrid Automaton Approach", linestyle=:dash)
vline!(phase_plt, [gₚ, -gₚ] .* 1e6, linestyle=:dash, color=:orange, label="Collision Threshold")

# Display all plots
display(filippov_plt)
display(hybrid_plt)
display(compare_plt)
display(phase_plt)

# Export plots (if needed)
savefig(filippov_plt, "filippov_approach