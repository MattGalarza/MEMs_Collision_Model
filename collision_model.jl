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
    E::T = 180e9             # Young's modulus (Pa)
    eta::T = 1.849e-5        # Dynamic viscosity of air (Pa·s)
    c::T = 0.015             # Damping scaling factor
    g0::T = 14e-6            # Electrode gap (m)
    Tp::T = 120e-9           # Thickness of parylene layer (m)
    Tf::T = 25e-6            # Device thickness (m)
    gss::T = 14e-6           # Soft-stopper initial gap (m)
    rho::T = 2333.0          # Density of silicon (kg/m³)
    cp::T = 5e-12            # Capacitance of parylene (F)
    wt::T = 9e-6             # Electrode width, top (m)
    wb::T = 30e-6            # Electrode width, bottom (m)
    ws::T = 14.7e-6          # Suspension spring width (m)
    Lss::T = 1400e-6         # Suspension spring length (m)
    Lff::T = 450e-6          # Electrode length (m)
    Leff::T = 400e-6         # Electrode effective overlap length (m)
    e::T = 8.85e-12          # Permittivity of free space (F/m)
    ep::T = 3.2              # Permittivity of parylene
    Vbias::T = 3.0           # Bias voltage (V)
    Rload::T = 0.42e6        # Load resistance (Ω)
    N::Int = 160             # Number of electrodes
    kss::T = 6.0             # Soft-stopper spring force (N/m)

    # Dependent parameters
    gp::T = :($(g0 - 2 * Tp))                      # Initial gap including parylene layer (m)
    wavg::T = :($((wt + wb) / 2))                  # Average electrode width (m)
    a::T = :($((wb - wt) / Leff))                  # Tilt factor
    k1::T = :($((4 / 6) * ((E * Tf * (ws^3)) / (Lss^3))))  # Spring constant (N/m)
    k3::T = :($((18 / 25) * ((E * Tf * ws) / (Lss^3))))    # Cubic spring constant (N/m³)
    I::T = :($((1 / 12) * Tf * (wavg^3)))          # Electrode moment of inertia (m⁴)
    m2::T = :($((33 / 140) * rho * Tf * Lff * wavg))  # Mass of electrode (kg)
    ke::T = :($((1 / 4) * E / Lff^3 * I))          # Electrode spring constant (N/m)
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
    # Fsp = - k1 * x1 - k3 * (x1^3) # Suspension beam force
    Fsp = -k1 * x1 # Suspension beam force
    if abs(x1) < gss
        Fss = 0.0
    else
        Fss = -kss * (abs(x1) - gss) * sign(x1) # Soft-stopper force
    end
    Fs = Fsp + Fss
    return Fs
end

# Electrode collision force, Fc
function collision(x1, x2, m2, ke, gp)
    if abs(x2) < gp
        Fc = -ke * (x1 - x2)
    else
        Fc = -ke * (x1 - x2) + ke * (abs(x2) - gp) * sign(x2)
    end
    return m2, Fc
end

# Viscous damping, Fd
function damping2(x2, x2dot, a, c, gp, Leff, Tf, eta)  
    # Damping LHS
    ul = gp + x2
    ll = ul + a * Leff
    A1_l = -ul * Leff / (ul + 2 * a) * x2dot
    A2_l = 12 * eta * a^-2 / (2 * ul + a * Leff) * x2dot
    Fd_l = Tf * (Leff * A2_l - 6 * eta * Leff / a / ul / ll * A1_l + 
                6 * eta * Leff * a^-4 / ll * x2dot + 
                12 * eta * a^-3 * x2dot * log(abs(ul / ll)))
    
    # Damping RHS
    ur = gp - x2
    lr = ur + a * Leff
    A1_r = -ur * Leff / (ur + 2 * a) * x2dot
    A2_r = 12 * eta * a^-2 / (2 * ur + a * Leff) * x2dot
    Fd_r = Tf * (Leff * A2_r - 6 * eta * Leff / a / ur / lr * A1_r + 
                6 * eta * Leff * a^-4 / lr * x2dot - 
                12 * eta * a^-3 * x2dot * log(abs(ur / lr)))
    
    # Total damping
    Fd = -c * (Fd_l + Fd_r)
    return Fd
end

function damping(x2, x2dot, a, c, gp, Leff, Tf, eta)  
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
end

# Electrostatic coupling, Fe
function electrostatic(x1, x2, Qvar, g0, gp, a, e, ep, cp, wt, wb, ke, E, I, Leff, Tf, Tp, N)
    # Function for variable capacitance, Cvar
    function Cvar_func(x2)
        if abs(x2) < gp # electrodes are NOT in contact
            crl = (e * ep * Leff * Tf) / Tp
            # RHS                           
            Cair_r = (e * Tf / a) * log((gp - x2 + 2 * a * Leff) / (gp - x2))
            Cvar_r = 1 / (1/crl + 1/Cair_r + 1/crl)
            # LHS
            Cair_l = (e * Tf / a) * log((gp + x2 + 2 * a * Leff) / (gp + x2))
            Cvar_l = 1 / (1/crl + 1/Cair_l + 1/crl)
            # Total variable capacitance, Cvar
            Cvar_value = (N / 2) * (Cvar_r + Cvar_l)
        elseif abs(x2) >= gp # electrodes are in contact
            crl = (e * ep * Leff * Tf) / Tp
            # d = ((wb - wt) + (abs(x1) - gp)) / Leff
            # d = (wb - wt) - (((ke * (abs(x2) - gp)) * (Leff^3)) / (3 * E * I))
            # Colliding electrode
            Cair_c = ((e * Tf * Leff)/abs(x2)) * (log((2 * Tp + abs(x2))/(2 * Tp)))
            Cvar_c = 1 / (1/crl + 1/Cair_c + 1/crl)
            # Non-colliding electrode
            Cair_nc = (e * Tf / a) * log((gp + abs(x2) + 2 * a * Leff) / (gp + abs(x2)))
            Cvar_nc = 1 / (1/crl + 1/Cair_nc + 1/crl)
            # Total variable capacitance, Cvar
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
    NFe = -((Qvar^2) / (2 * Ctotal^2)) * dC
    Fe = NFe / (N / 2)
    return Ctotal, Fe
end

function CoupledSystem!(dz, z, p, t, current_acceleration)
    # Unpack state variables
    z1, z2, z3, z4, z5, Vout = z

    # Compute forces
    Fs = spring(z1, p.k1, p.k3, p.gss, p.kss)
    m2, Fc = collision(z1, z3, p.m2, p.ke, p.gp)
    Fd = damping(z3, z4, p.a, p.c, p.gp, p.Leff, p.Tf, p.eta)
    Ctotal, Fe = electrostatic(z1, z3, z5, p.g0, p.gp, p.a, p.e, p.ep, p.cp, p.wt, p.wb, p.ke, p.E, p.I, p.Leff, p.Tf, p.Tp, p.N)

    # Use current_acceleration as the external force
    Fext = current_acceleration

    # Compute derivatives
    dz[1] = z2
    dz[2] = (Fs + (p.N / 2) * Fc) / p.m1 - Fext
    dz[3] = z4
    dz[4] = (-Fc + Fd + Fe) / m2 - Fext
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
alpha = 3.0 # Applied acceleration constant
g = 9.81 # Gravitational constant 
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
# eqn = ODEProblem(AnalyticalModel.CoupledSystem!, z0, tspan, p_new)

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

# Iterate over each solution point to compute forces
for (i, t) in enumerate(sol.t)
    # Extract state variables at time t
    z = sol.u[i]
    z1, z2, z3, z4, z5, Vout = z
    
    # Compute Fs (Suspension spring force)
    Fs = AnalyticalModel.spring(z1, p_new.k1, p_new.k3, p_new.gss, p_new.kss)
    push!(Fs_array, Fs)
    
    # Compute Fc (Collision force)
    m2_val, Fc = AnalyticalModel.collision(z1, z3, p_new.m2, p_new.ke, p_new.gp)
    push!(Fc_array, Fc)
    
    # Compute Fd (Damping force)
    Fd = AnalyticalModel.damping(z3, z4, p_new.a, p_new.c, p_new.gp, p_new.Leff, p_new.Tf, p_new.eta)
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
