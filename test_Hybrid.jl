using DifferentialEquations, Plots

# Constants for MEMS system
# Physical Parameters
m₁ = 1.0e-7    # kg, shuttle mass
m₂ = 5.0e-8    # kg, mobile electrode mass
k₁ = 1.0       # N/m, primary spring constant
kₛₛ = 10.0     # N/m, suspension spring constant
kₑ = 5.0       # N/m, collision spring constant
gₛₛ = 0.1e-6   # m, suspension spring gap
gₚ = 0.5e-6    # m, collision threshold
R = 1.0e6      # Ω, resistance
VDC = 5.0      # V, applied voltage
N = 2          # Number of electrodes

# Capacitance parameters
a = 0.2e-6     # m, distance parameter
ε = 8.85e-12   # F/m, permittivity of free space
εₚ = 1.0       # relative permittivity
Leff = 50e-6   # m, effective length
Tf = 2e-6      # m, finger thickness
Tp = 2e-6      # m, plate thickness
Cp = 1e-15     # F, parasitic capacitance
η = 1.8e-5     # Pa·s, air viscosity
c = 0.8        # damping coefficient

# Applied force parameters
F₀ = 1.0e-6    # N, amplitude of sinusoidal force
ω = 2π * 1000  # rad/s, angular frequency (1000 Hz)

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

# Calculate capacitance for pre-collision mode
function calculate_capacitance_pre_collision(x₂)
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

# Calculate capacitance for post-collision mode
function calculate_capacitance_post_collision(x₂)
    # Capacitance between rigid insulator and electrode
    Cᵣᵢ = (ε * εₚ * Leff * Tf) / Tp
    
    # Contact side air capacitance
    Cₐᵢᵣ₋ₖ = (ε * Tf * Leff / abs(x₂)) * log((2 * Tp + abs(x₂)) / (2 * Tp))
    
    # Non-contact side air capacitance
    Cₐᵢᵣ₋ₙₖ = (ε * Tf / (2 * a)) * log((gₚ + x₂ + 2 * a * Leff) / (gₚ + x₂))
    
    # Variable capacitances
    Cᵥₐᵣ₋ₖ = 1.0 / (2.0 / Cᵣᵢ + 1.0 / Cₐᵢᵣ₋ₖ)
    Cᵥₐᵣ₋ₙₖ = 1.0 / (2.0 / Cᵣᵢ + 1.0 / Cₐᵢᵣ₋ₙₖ)
    
    # Total variable capacitance
    Cᵥₐᵣ = (N / 2) * (Cᵥₐᵣ₋ₖ + Cᵥₐᵣ₋ₙₖ)
    
    return Cᵥₐᵣ
end

# Calculate electrostatic force
function calculate_electrostatic_force(x₂, q, mode)
    # Calculate capacitance based on mode
    if mode == :pre_collision
        Cᵥₐᵣ = calculate_capacitance_pre_collision(x₂)
    else
        Cᵥₐᵣ = calculate_capacitance_post_collision(x₂)
    end
    
    Ctotal = Cᵥₐᵣ + Cp
    
    # Calculate derivative of Cᵥₐᵣ with respect to x₂ using finite difference
    δ = 1e-10  # Small delta for derivative calculation
    
    if mode == :pre_collision
        Cᵥₐᵣ_plus = calculate_capacitance_pre_collision(x₂ + δ)
        Cᵥₐᵣ_minus = calculate_capacitance_pre_collision(x₂ - δ)
    else
        Cᵥₐᵣ_plus = calculate_capacitance_post_collision(x₂ + δ)
        Cᵥₐᵣ_minus = calculate_capacitance_post_collision(x₂ - δ)
    end
    
    dCᵥₐᵣ_dx₂ = (Cᵥₐᵣ_plus - Cᵥₐᵣ_minus) / (2 * δ)
    
    # Electrostatic force
    Fe = -(2 / N) * (q^2 / (2 * (Cᵥₐᵣ + Cp)^2)) * dCᵥₐᵣ_dx₂
    
    return Fe, Ctotal
end

# Calculate damping force
function calculate_damping_force(x₂, ẋ₂)
    # Left side
    A_L1 = ((gₚ + x₂) * Leff / (gₚ + x₂ + 2 * a)) * ẋ₂
    A_L2 = (12 * η * a^(-2) / (2 * (gₚ + x₂) + a * Leff)) * ẋ₂
    
    F_dL = Tf * (Leff * A_L2 - (6 * η * Leff / (a * (gₚ + x₂) * (gₚ + x₂ + a * Leff))) * A_L1 +
                6 * η * Leff * a^(-4) * (gₚ + x₂ + a * Leff) * ẋ₂ +
                12 * η * a^(-3) * ẋ₂ * log(abs((gₚ + x₂) / (gₚ + x₂ + a * Leff))))
    
    # Right side (simplified for stability)
    if abs(gₚ - x₂) < 1e-12  # Avoid division by zero
        F_dR = F_dL  # Approximate with left side
    else
        A_R1 = -((gₚ - x₂) * Leff / (gₚ - x₂ + 2 * a)) * ẋ₂
        A_R2 = (12 * η * a^(-2) / (2 * (gₚ - x₂) + a * Leff)) * ẋ₂
        
        F_dR = Tf * (Leff * A_R2 - (6 * η * Leff / (a * (gₚ - x₂) * (gₚ - x₂ + a * Leff))) * A_R1 +
                    6 * η * Leff * a^(-4) * (gₚ - x₂ + a * Leff) * ẋ₂ +
                    12 * η * a^(-3) * ẋ₂ * log(abs((gₚ - x₂) / (gₚ - x₂ + a * Leff))))
    end
    
    # Apply damping coefficient
    return -c * (F_dL + F_dR)
end

# Modified damping force for post-collision (adds extra damping due to deformation)
function calculate_damping_force_post_collision(x₂, ẋ₂)
    # Use pre-collision damping as a base
    base_damping = calculate_damping_force(x₂, ẋ₂)
    
    # Add additional damping due to deformation
    deformation_factor = 1.5 * (abs(x₂) / gₚ - 1.0)
    if deformation_factor < 0.0
        deformation_factor = 0.0
    end
    
    return base_damping * (1.0 + deformation_factor)
end

# Define the pre-collision mode dynamics
function pre_collision!(du, u, p, t)
    # Unpack state variables: [x₁, ẋ₁, x₂, ẋ₂, q]
    z₁, z₂, z₃, z₄, z₅ = u
    
    # Calculate the sinusoidal applied force
    Fapplied = F₀ * sin(ω * t)
    
    # Calculate forces for pre-collision mode
    Fs = calculate_suspension_force(z₁, k₁, kₛₛ, gₛₛ)
    Fc = -kₑ * (z₁ - z₃)
    Fe, Ctotal = calculate_electrostatic_force(z₃, z₅, :pre_collision)
    Fd = calculate_damping_force(z₃, z₄)
    
    # State derivatives
    du[1] = z₂
    du[2] = (1/m₁) * (Fs + N*Fc - m₁*Fapplied)
    du[3] = z₄
    du[4] = (1/m₂) * (Fc + Fd + Fe - m₂*Fapplied)
    du[5] = (1/R) * (VDC - z₅/Ctotal)
end

# Define the post-collision mode dynamics
function post_collision!(du, u, p, t)
    # Unpack state variables: [x₁, ẋ₁, x₂, ẋ₂, q]
    z₁, z₂, z₃, z₄, z₅ = u
    
    # Calculate the sinusoidal applied force
    Fapplied = F₀ * sin(ω * t)
    
    # Calculate forces for post-collision mode
    Fs = calculate_suspension_force(z₁, k₁, kₛₛ, gₛₛ)
    Fc = -kₑ * (z₁ - z₃) + kₑ * (abs(z₃) - gₚ) * sign(z₃)
    Fe, Ctotal = calculate_electrostatic_force(z₃, z₅, :post_collision)
    Fd = calculate_damping_force_post_collision(z₃, z₄)
    
    # State derivatives with modified mass (rotation mode)
    du[1] = z₂
    du[2] = (1/m₁) * (Fs + N*Fc - m₁*Fapplied)
    du[3] = z₄
    du[4] = (1/(2*m₂)) * (Fc + Fd + Fe - m₂*Fapplied)  # Note the modified mass for rotation
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

# Event handlers with mode transition recording
function enter_post_collision!(integrator)
    # Switch to post-collision dynamics
    integrator.f = post_collision!
    
    # Store mode transition info in problem data
    if !haskey(integrator.opts.callback.discrete_callbacks[1].affect!.kwargs, :mode_transitions)
        integrator.opts.callback.discrete_callbacks[1].affect!.kwargs[:mode_transitions] = [integrator.t]
    else
        push!(integrator.opts.callback.discrete_callbacks[1].affect!.kwargs[:mode_transitions], integrator.t)
    end
end

function enter_pre_collision!(integrator)
    # Switch to pre-collision dynamics
    integrator.f = pre_collision!
    
    # Store mode transition info in problem data
    if !haskey(integrator.opts.callback.discrete_callbacks[2].affect!.kwargs, :mode_transitions)
        integrator.opts.callback.discrete_callbacks[2].affect!.kwargs[:mode_transitions] = [integrator.t]
    else
        push!(integrator.opts.callback.discrete_callbacks[2].affect!.kwargs[:mode_transitions], integrator.t)
    end
end

# Function to check if a state is in post-collision mode
function is_post_collision(x₂)
    return abs(x₂) >= gₚ
end

# Simulation settings
tspan = (0.0, 0.005)  # 5 ms simulation time

# Initial conditions: for demonstrating collision
x₁₀ = 0.0             # Initial shuttle position
ẋ₁₀ = 0.0             # Initial shuttle velocity
x₂₀ = 0.3e-6          # Initial mobile electrode position (below threshold)
ẋ₂₀ = 0.0             # Initial mobile electrode velocity
q₀ = 1.0e-12          # Initial charge

u₀ = [x₁₀, ẋ₁₀, x₂₀, ẋ₂₀, q₀]

# Create callbacks for mode transitions
cb_to_post = DiscreteCallback(to_post_collision, enter_post_collision!)
cb_to_pre = DiscreteCallback(to_pre_collision, enter_pre_collision!)
cb_set = CallbackSet(cb_to_post, cb_to_pre)

# Determine initial mode based on initial conditions
initial_mode = abs(u₀[3]) < gₚ ? pre_collision! : post_collision!

# Solve using Hybrid Automaton approach with appropriate solver for stiff system
prob = ODEProblem(initial_mode, u₀, tspan)
sol = solve(prob, Rodas5(), callback=cb_set, reltol=1e-6, abstol=1e-8)

# Collect mode transitions from callback data
mode_transitions = []
if haskey(cb_to_post.affect!.kwargs, :mode_transitions)
    append!(mode_transitions, cb_to_post.affect!.kwargs[:mode_transitions])
end
if haskey(cb_to_pre.affect!.kwargs, :mode_transitions)
    append!(mode_transitions, cb_to_pre.affect!.kwargs[:mode_transitions])
end
sort!(mode_transitions)

# Create plots
plt = plot(layout=(3,2), size=(900, 600), title="Hybrid Automaton Approach for MEMS Collision")

# Shuttle position (x₁)
plot!(plt[1], sol.t, sol[1, :] .* 1e6, 
      xlabel="Time (s)", ylabel="x₁ (μm)", label="Shuttle Position",
      subplot=1)
for t in mode_transitions
    vline!(plt[1], [t], linestyle=:dash, color=:red, alpha=0.5, label=false)
end

# Shuttle velocity (ẋ₁)
plot!(plt[2], sol.t, sol[2, :] .* 1e3, 
      xlabel="Time (s)", ylabel="ẋ₁ (mm/s)", label="Shuttle Velocity",
      subplot=2)
for t in mode_transitions
    vline!(plt[2], [t], linestyle=:dash, color=:red, alpha=0.5, label=false)
end

# Mobile electrode position (x₂)
plot!(plt[3], sol.t, sol[3, :] .* 1e6, 
      xlabel="Time (s)", ylabel="x₂ (μm)", label="Mobile Electrode Position",
      subplot=3)
hline!(plt[3], [gₚ, -gₚ] .* 1e6, linestyle=:dash, color=:orange, label="Collision Threshold")
for t in mode_transitions
    vline!(plt[3], [t], linestyle=:dash, color=:red, alpha=0.5, label=false)
end

# Mobile electrode velocity (ẋ₂)
plot!(plt[4], sol.t, sol[4, :] .* 1e3, 
      xlabel="Time (s)", ylabel="ẋ₂ (mm/s)", label="Mobile Electrode Velocity",
      subplot=4)
for t in mode_transitions
    vline!(plt[4], [t], linestyle=:dash, color=:red, alpha=0.5, label=false)
end

# Charge (q)
plot!(plt[5], sol.t, sol[5, :] .* 1e12, 
      xlabel="Time (s)", ylabel="q (pC)", label="Charge",
      subplot=5)
for t in mode_transitions
    vline!(plt[5], [t], linestyle=:dash, color=:red, alpha=0.5, label=false)
end

# Applied Force
applied_force = [F₀ * sin(ω * t) * 1e6 for t in sol.t]
plot!(plt[6], sol.t, applied_force, 
      xlabel="Time (s)", ylabel="Force (μN)", label="Applied Force",
      subplot=6)
for t in mode_transitions
    vline!(plt[6], [t], linestyle=:dash, color=:red, alpha=0.5, label=false)
end

# Phase plane plot (x₂ vs ẋ₂)
phase_plt = plot(size=(600, 400), title="Phase Plane: x₂ vs ẋ₂")
plot!(phase_plt, sol[3, :] .* 1e6, sol[4, :] .* 1e3, 
      label="Trajectory", xlabel="x₂ (μm)", ylabel="ẋ₂ (mm/s)")
vline!(phase_plt, [gₚ, -gₚ] .* 1e6, linestyle=:dash, color=:orange, label="Collision Threshold")

# Mode indicator plot
mode_indicator = [is_post_collision(x₂) ? 1.0 : 0.0 for x₂ in sol[3, :]]
mode_plt = plot(size=(800, 300), title="System Mode")
plot!(mode_plt, sol.t, mode_indicator, 
      xlabel="Time (s)", ylabel="Mode", label="System Mode",
      ylim=(-0.1, 1.1), yticks=[0, 1], yformatter=y->["Pre-collision", "Post-collision"][Int(y)+1])
for t in mode_transitions
    vline!(mode_plt, [t], linestyle=:dash, color=:red, alpha=0.5, label=false)
end

# Display plots
display(plt)
display(phase_plt)
display(mode_plt)

# Save plots
savefig(plt, "hybrid_states.png")
savefig(phase_plt, "hybrid_phase.png")
savefig(mode_plt, "hybrid_mode.png")