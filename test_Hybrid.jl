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
        A_R1 = -((gₚ - x₂) * Leff / (gₚ - x₂ + 2 * a)) * ẋ