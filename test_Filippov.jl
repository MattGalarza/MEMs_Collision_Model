
# ------------------------- Libraries -------------------------
using Sundials, DifferentialEquations, Plots
using Parameters, ForwardDiff, SpecialFunctions, NaNMath, Random

# ------------------------- Parameter Struct -------------------------
@with_kw mutable struct MEMSParams{T<:Real}
    g0::T=14e-6; Tp::T=120e-9; Tf::T=25e-6
    wt::T=9e-6; wb::T=30e-6; ws::T=14.7e-6; wss::T=14e-6
    Leff::T=400e-6; Lff::T=450e-6; Lsp::T=1400e-6; Lss::T=1000e-6; gss::T=14e-6
    m1::T=2.0933e-6; rho::T=2330.0; E::T=170e9; e::T=8.85e-12; ep::T=3.2
    eta::T=1.849e-5; c::T=0.015
    N::T=160.0; Cp::T=5e-12; Vbias::T=3.0; Rload::T=0.42e6
    gp::T=0.0; a::T=0.0; m2::T=0.0; I::T=0.0; ke::T=0.0; k1::T=0.0; k3::T=0.0; kss::T=0.0
end

# ------------------------- Create Dependent Parameters -------------------------
function create_params(; kwargs...)
    p=MEMSParams(;kwargs...)
    p.gp=p.g0-2p.Tp; p.a=(p.wb-p.wt)/p.Leff
    p.I=(1/48)*p.Lff*p.Tf*(p.wb+p.wt)*(p.wb^2+p.wt^2)
    modal=0.236+0.045*(1-p.wt/p.wb)
    phys=0.5*p.Lff^2*p.rho*p.Tf*(p.wb+p.wt); p.m2=modal*phys
    F=1; num=p.E*p.Tf*p.wt^2*(p.wb-p.wt)^3
    den=6F*p.Lff^3*((p.wb-3p.wt)*(p.wb-p.wt)+2p.wt^2*(log(p.Lff*p.wb)-log(p.Lff*p.wt)))
    p.ke=num/den
    p.k1=(4/6)*(p.E*p.Tf*p.ws^3)/(p.Lsp^3)
    p.k3=(18/25)*(p.E*p.Tf*p.ws)/(p.Lsp^3)
    p.kss=(p.E*p.Tf*p.wss^3)/(p.Lss^3)
    return p
end
params=create_params();

# ------------------------- Filippov Regularization -------------------------
const ε=1e-4  # widened smoothing
smoothH(x;ε=ε)=0.5*(1+tanh(x/ε))
smoothS(x;ε=ε)=tanh(x/ε)

# ------------------------- Force Expressions -------------------------
function spring_fss(x1,p)
    Fsp=-p.k1*x1
    Δ=abs(x1)-p.gss
    Fss=-p.kss*Δ*smoothS(x1;ε=ε)
    α=smoothH(Δ;ε=ε)
    return Fsp*(1-α)+Fss*α
end

function collision_f(p,x1,x2)
    Δ=abs(x2)-p.gp
    Fc0=-p.ke*(x2-x1)
    Fc1= p.ke*(smoothS(x2;ε=ε)*p.gp - x1)
    β=smoothH(Δ;ε=ε)
    return Fc0*(1-β)+Fc1*β
end

function damping_f(p,x2,v)
    Δ=abs(x2)-p.gp; Fd0=0.0
    if abs(x2)<p.gp
        h=p.gp-x2; t1=12*p.Tf*p.eta*v
        t2=(2*p.Leff*p.a+(2h+p.Leff*p.a)*log(h/(h+p.Leff*p.a)))
        t3=p.a^3*(2h+p.Leff*p.a); Fr=(t1*t2)/t3
        h2=p.gp+x2; t2=(2*p.Leff*p.a+(2h2+p.Leff*p.a)*log(h2/(h2+p.Leff*p.a)))
        Fl=(t1*t2)/t3; Fd0=p.c*(Fr+Fl)
    end
    Fd1=-1e3*v   # small viscous bleed
    β=smoothH(Δ;ε=ε)
    return Fd0*(1-β)+Fd1*β
end

function capacitance(x2,p)
    Crl=(p.e*p.ep*p.Leff*p.Tf)/p.Tp
    if abs(x2)<p.gp
        Cr=((p.e*p.Tf)/p.a)*NaNMath.log((p.gp-x2+p.a*p.Leff)/(p.gp-x2))
        Cl=((p.e*p.Tf)/p.a)*NaNMath.log((p.gp+x2+p.a*p.Leff)/(p.gp+x2))
        C1=1/((2/p.Cp)+(1/Cr)); C2=1/((2/p.Cp)+(1/Cl))
        return (p.N/2)*(C1+C2)
    else
        k=(2*p.Tp)/(p.a*p.Leff); C0=7.1053e-14; Cinf=9.952e-11
        Cr=C0+(Cinf-C0)*(log(1+k*(abs(x2)-p.gp))/log(1+k*p.a*p.Leff))
        C1=1/((2/p.Cp)+(1/Cr)); C2=C1
        return (p.N/2)*(C1+C2)
    end
end

function electrostatic(x2,q,Cvar,p)
    Ctot=Cvar+p.Cp
    h=1e-10; C1=capacitance(x2-h,p); C2=capacitance(x2+h,p)
    dC=(C2-C1)/(2h); Fe=(-0.5*(q^2/Ctot^2)*dC)/(p.N/2)
    return Ctot,Fe
end

# ------------------------- Filippov Dynamics -------------------------
function CoupledDynamics_Filippov!(du,u,p_,t)
    p,Fext=p_
    x1,v1,x2,v2,q,V=u
    Fs=spring_fss(x1,p); Fc=collision_f(p,x1,x2)
    Fd=damping_f(p,x2,v2); Cvar=capacitance(x2,p)
    Ctot,Fe=electrostatic(x2,q,Cvar,p)
    du[1]=v1; du[2]=(Fs+(p.N/2)*Fc)/p.m1 - Fext(t)
    du[3]=v2; du[4]=(-Fc+Fd+Fe)/p.m2 - Fext(t)
    du[5]=(p.Vbias - q/Ctot)/p.Rload
    du[6]=(p.Vbias - q/Ctot - V)/(p.Rload*Ctot)
end

# ------------------------- External Force -------------------------
f=20.0; A=2*9.81; t_ramp=0.2; use_ramp=true
function sine_force(t;A=A,f=f,use_ramp=use_ramp,t_ramp=t_ramp)
    ramp=use_ramp&&t<t_ramp?t/t_ramp:1.0; return A*ramp*sin(2π*f*t)
end
function white_noise_force(t;amplitude=9.81,seed=12345,dt=0.001)
    rng=Random.MersenneTwister(seed); k=floor(Int,t/dt)
    Random.seed!(rng,seed); for _ in 1:k rand(rng) end
    return amplitude*(2*rand(rng)-1)
end
Fext_input=t->sine_force(t)

# ------------------------- Initial Conditions -------------------------
x1_0,v1_0=0.0,0.0; x2_0,v2_0=0.0,0.0
C0=capacitance(x2_0,params); q0=params.Vbias*(C0+params.Cp)
V0=params.Vbias - q0/(C0+params.Cp)
u0=[x1_0,v1_0,x2_0,v2_0,q0,V0]

# ------------------------- Solve -------------------------
tspan=(0.0,0.5)
prob=ODEProblem(CoupledDynamics_Filippov!,u0,tspan,(params,Fext_input))
sol=solve(prob,Rosenbrock23();abstol=1e-10,reltol=1e-8,dtmin=1e-15,dtmax=1e-4,force_dtmin=true,saveat=0.0001)
println("Solved with ",length(sol.t)," steps.")

# ------------------------- Plotting -------------------------
plot(sol.t,getindex.(sol.u,1),xlabel="Time (s)",ylabel="x1 (m)",title="Shuttle Displacement")
plot(sol.t,getindex.(sol.u,2),xlabel="Time (s)",ylabel="v1 (m/s)",title="Shuttle Velocity")
plot(sol.t,getindex.(sol.u,3),xlabel="Time (s)",ylabel="x2 (m)",title="Electrode Displacement")
plot(sol.t,getindex.(sol.u,4),xlabel="Time (s)",ylabel="v2 (m/s)",title="Electrode Velocity")
plot(sol.t,getindex.(sol.u,5),xlabel="Time (s)",ylabel="q (C)",title="Charge")
plot(sol.t,getindex.(sol.u,6),xlabel="Time (s)",ylabel="V (V)",title="Output Voltage")

