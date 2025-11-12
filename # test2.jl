# hybrid_mems_model.jl (Standalone Simulation Script)

using DifferentialEquations
using Parameters
using ForwardDiff
using Plots

# -------------------------- System Parameters and Utilities --------------------------
@with_kw mutable struct Params{T<:Real}
    g0::T = 14e-6
    Tp::T = 120e-9
    Tf::T = 25e-6
    wt::T = 9e-6
    wb::T = 30e-6
    ws::T = 14.7e-6
    wss::T = 14e-6
    Leff::T = 400e-6
    Lf::T = 450e-6
    Lsp::T = 1400e-6
    Lss::T = 1000e-6
    gss::T = 14e-6
    m1::T = 2.0933e-6
    rho::T = 2330.0
    E::T = 170e9
    e::T = 8.85e-12
    ep::T = 3.2
    eta::T = 1.849e-5
    c::T = 0.015
    lambda::T = 70e-9
    sigmap::T = 1.016
    N::Int = 160
    cp::T = 5e-12
    Vbias::T = 3.0
    Rload::T = 0.42e6
    gp::T = g0 - 2 * Tp
    a::T = (wb - wt) / Leff
    m2::T = 0.0
    I::T = 0.0
    ke::T = 0.0
    k1::T = 0.0
    k3::T = 0.0
    kss::T = 6.0
end

function create_params(p::Params{T}; verbose=true) where T<:Real
    p.I = (1/48) * p.Lf * p.Tf * (p.wb + p.wt) * (p.wb^2 + p.wt^2)
    modalcoeff = 0.236 + 0.045 * (1.0 - p.wt / p.wb)
    m2Physical = 0.5 * p.Lf^2 * p.rho * p.Tf * (p.wb + p.wt)
    p.m2 = modalcoeff * m2Physical
    F = 1
    num = (p.E * p.Tf * p.wt^2 * ((p.wb - p.wt)^3))
    dem = 6 * F * p.Lf^3 * ((p.wb - 3 * p.wt) * (p.wb - p.wt) + 2 * p.wt^2 * (log(p.Lf * p.wb) - log(p.Lf * p.wt)))
    p.ke = num / dem
    p.k1 = (4.0/6.0) * ((p.E * p.Tf * (p.ws^3)) / (p.Lsp^3))
    p.k3 = (18.0/25.0) * ((p.E * p.Tf * p.ws) / (p.Lsp^3))
    p.kss = (p.E * p.Tf * (p.wss^3)) / (p.Lss^3)
    return p
end

function spring(x1, k1, k3, gss, kss)
    Fsp = -k1 * x1 - 0.25 * k3 * (x1^3)
    Fss = abs(x1) <= gss ? 0.0 : -kss * (abs(x1) - gss) * sign(x1)
    return Fsp + Fss
end

function collision(x1, x2, m2, ke, gp)
    if abs(x2) < gp
        return -ke * (x1 - x2), "translational"
    else
        return -ke * (gp * sign(x2) - x1), "rotational"
    end
end

function electrostatic(x1, x2, Qvar, g0, gp, a, e, ep, cp, wt, wb, ke, E, I, Leff, Tf, Tp, N)
    crl = (e * ep * Leff * Tf) / Tp
    Cair_r = (e * Tf / a) * log((gp - x2 + a * Leff) / (gp - x2))
    Cair_l = (e * Tf / a) * log((gp + x2 + a * Leff) / (gp + x2))
    Cvar_r = 1 / (2 / crl + 1 / Cair_r)
    Cvar_l = 1 / (2 / crl + 1 / Cair_l)
    Cvar_value = (N / 2) * (Cvar_r + Cvar_l)
    dC = ForwardDiff.derivative(z -> electrostatic(x1, z, Qvar, g0, gp, a, e, ep, cp, wt, wb, ke, E, I, Leff, Tf, Tp, N)[1], x2)
    Ctotal = Cvar_value + cp
    Fe = -((Qvar^2 / (2 * Ctotal^2)) * dC) / (N / 2)
    return Ctotal, Fe
end

# Sine input force
t_ramp = 0.2
f = 20.0
alpha = 2.1
g = 9.81
A = alpha * g
ramp(t) = t < t_ramp ? t / t_ramp : 1.0
Fext_input(t) = A * ramp(t) * sin(2 * pi * f * t)

# ---------------------------- Hybrid Simulation Code ----------------------------
mutable struct SimulationMode
    contact_active::Bool
end
sim_mode = SimulationMode(false)

function contact_condition(u, t, integrator)
    x2 = u[3]
    return abs(x2) - integrator.p.gp
end

function contact_affect!(integrator)
    u = integrator.u
    sim_mode.contact_active = true
    u[3] = sign(u[3]) * integrator.p.gp
    u[4] = u[2]
end

function release_condition(u, t, integrator)
    return abs(u[3]) - integrator.p.gp - 1e-9
end

function release_affect!(integrator)
    sim_mode.contact_active = false
end

contact_cb = ContinuousCallback(contact_condition, contact_affect!)
release_cb = ContinuousCallback(release_condition, release_affect!)
cbset = CallbackSet(contact_cb, release_cb)

function CoupledSystemHybrid!(dz, z, p, t)
    z1, z2, z3, z4, z5, Vout = z
    Fext = Fext_input(t)
    if sim_mode.contact_active
        z3 = sign(z3) * p.gp
        z4 = z2
    end
    Fs = spring(z1, p.k1, p.k3, p.gss, p.kss)
    Fc, _ = collision(z1, z3, p.m2, p.ke, p.gp)
    Fd = -p.c * (z2 - z4)
    Ctotal, Fe = electrostatic(z1, z3, z5, p.g0, p.gp, p.a, p.e, p.ep, p.cp, p.wt, p.wb, p.ke, p.E, p.I, p.Leff, p.Tf, p.Tp, p.N)
    dz[1] = z2
    dz[2] = (Fs + (p.N / 2) * Fc) / p.m1 - Fext
    dz[3] = z4
    dz[4] = (-Fc + Fd + Fe) / p.m2 - Fext
    dz[5] = (p.Vbias - (z5 / Ctotal)) / p.Rload
    dz[6] = (p.Vbias - z5 / Ctotal - Vout) / (p.Rload * Ctotal)
end

p = create_params(Params{Float64}(); verbose=false)
x10 = 0.0; x10dot = 0.0; x20 = 0.0; x20dot = 0.0
C0, _ = electrostatic(x10, x20, 0.0, p.g0, p.gp, p.a, p.e, p.ep, p.cp, p.wt, p.wb, p.ke, p.E, p.I, p.Leff, p.Tf, p.Tp, p.N)
Q0 = p.Vbias * C0; V0 = p.Vbias - Q0 / C0
z0 = [x10, x10dot, x20, x20dot, Q0, V0]; tspan = (0.0, 0.5)
prob = ODEProblem(CoupledSystemHybrid!, z0, tspan, p)
sol = solve(prob, Rosenbrock23(); callback=cbset, abstol=1e-9, reltol=1e-6)

plot(sol.t, [u[1] for u in sol.u], label="x1 (Shuttle)")
plot!(sol.t, [u[3] for u in sol.u], label="x2 (Electrode)", xlabel="Time (s)", ylabel="Displacement (m)", title="MEMS Displacements")
