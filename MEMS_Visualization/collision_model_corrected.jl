#!/usr/bin/env julia
# Corrected trapezoidal-electrode MEMS model, SI units throughout.
# Julia-only implementation and numerical verification. See the companion TeX.
# Run: julia collision_model_corrected.jl --verify
# ODE runs: julia --project=. collision_model_corrected.jl --convergence
module CorrectedMEMS
using LinearAlgebra, Printf, Test
export Params, Model, constitutive, rhs!, energy, verify, simulate, convergence

Base.@kwdef struct Params
    g0::Float64 = 14e-6
    Tp::Float64 = 120e-9
    Tf::Float64 = 25e-6
    wt::Float64 = 9e-6             # narrow clamped root
    wb::Float64 = 30e-6            # wide free tip
    ws::Float64 = 14.7e-6
    wss::Float64 = 14e-6
    Lf::Float64 = 450e-6
    Leff::Float64 = 400e-6
    Lsp::Float64 = 1400e-6
    Lss::Float64 = 1000e-6
    gss::Float64 = 14e-6
    n_beams::Int = 80             # 160 gap branches; confirm device counting
    n_parallel::Int = 4
    n_series::Int = 6
    n_stop_beams::Int = 2
    gamma3::Float64 = 1.0         # geometry correction, not fitted here
    m_shuttle::Float64 = 2.0933e-6 # MUST exclude the explicit mobile beam mass
    rho::Float64 = 2330.0
    E::Float64 = 170e9
    eps0::Float64 = 8.85e-12
    epsr::Float64 = 3.2
    eta::Float64 = 1.849e-5
    mean_free_path::Float64 = 70e-9
    slip_coefficient::Float64 = 1.016
    film_scale::Float64 = 1.0
    c1::Float64 = 0.0             # uncalibrated reference, not measured zero
    ce::Float64 = 0.0             # damping of relative beam motion, per beam
    h_eff::Float64 = 50e-9
    eps_gap::Float64 = 2e-9
    eps_wall::Float64 = 0.5e-9
    eps_stop::Float64 = 1e-9
    seal_width::Float64 = 25e-9
    seal_at_contact::Bool = true  # candidate Robin-boundary closure
    gap_slope::Float64 = NaN      # default (wb-wt)/Lf; override with metrology
    kw::Float64 = 1e6
    pw::Float64 = 1.5
    cw::Float64 = 50.0
    ghs::Float64 = Inf            # disabled until the hard-stop gap is measured
    khs::Float64 = 1e9
    phs::Float64 = 1.5
    cp::Float64 = 5e-12
    Vbias::Float64 = 3.0
    Rload::Float64 = 0.42e6
end

# Stable positive-part regularization and its exact derivative.
softpos(z,e) = z >= 0 ? (z+hypot(z,e))/2 : e^2/(2*(hypot(z,e)-z))
dsoftpos(z,e) = z >= 0 ? (1+z/hypot(z,e))/2 :
    e^2/(2*hypot(z,e)*(hypot(z,e)-z))
smootherstep(z) = z <= 0 ? 0.0 : z >= 1 ? 1.0 : z^3*(10-15*z+6*z^2)

function gausslegendre(n)
    F=eigen(SymTridiagonal(zeros(n),[j/sqrt(4*j*j-1) for j in 1:n-1]))
    F.values, 2 .* F.vectors[1,:].^2
end
const GLX, GLW = gausslegendre(96)
quad(f,a,b) = (b-a)/2*sum(GLW[j]*f((a+b)/2+(b-a)/2*GLX[j]) for j in eachindex(GLX))
width(p,s) = p.wt+(p.wb-p.wt)*s/p.Lf
EI(p,s) = p.E*p.Tf*width(p,s)^3/12
shape(p,ke,s) = s == 0 ? 0.0 : ke*quad(z->(s-z)*(p.Lf-z)/EI(p,z),0.0,s)

struct Model
    p::Params
    ke::Float64
    k1::Float64
    k3::Float64
    kss::Float64
    gc::Float64
    alpha::Float64
    hd::Float64
    kp::Float64
    M::Matrix{Float64}
    Minv::Matrix{Float64}
    beta::Vector{Float64}
    Cstruct::Matrix{Float64}
    y::Vector{Float64}
    weights::Vector{Float64}
    B1::Vector{Float64}
    B2::Vector{Float64}
    panels::Int
end

function Model(p=Params();panels=512)
    @assert p.n_beams>0 && p.n_parallel>0 && p.n_series>0 && p.n_stop_beams>0
    @assert 0<p.Leff<=p.Lf && p.g0>2*p.Tp+p.h_eff && panels>=16
    @assert all(>(0),(p.Tf,p.wt,p.wb,p.ws,p.wss,p.Lsp,p.Lss,p.m_shuttle,
        p.rho,p.E,p.eps0,p.epsr,p.eta,p.h_eff,p.eps_gap,p.eps_wall,p.eps_stop,
        p.seal_width,p.Rload,p.kw,p.pw,p.phs,p.gss))
    @assert all(>=(0),(p.c1,p.ce,p.cw,p.film_scale,p.gamma3,p.cp,p.Tp,
        p.mean_free_path,p.slip_coefficient,p.khs)) && p.ghs>=p.gss
    alpha=isnan(p.gap_slope) ? (p.wb-p.wt)/p.Lf : p.gap_slope
    @assert isfinite(alpha) && alpha>=0
    ke=1/quad(s->(p.Lf-s)^2/EI(p,s),0.0,p.Lf)
    k1=p.n_parallel/p.n_series*p.E*p.Tf*p.ws^3/p.Lsp^3
    k3=p.gamma3*p.n_parallel/p.n_series^3*0.72*p.E*p.Tf*p.ws/p.Lsp^3
    kss=p.n_stop_beams*p.E*p.Tf*p.wss^3/(4*p.Lss^3)
    mass(i,j)=p.n_beams*quad(s->begin
        phi=shape(p,ke,s); B=(1-phi,phi)
        p.rho*p.Tf*width(p,s)*B[i]*B[j]
    end,0.0,p.Lf)
    M=[p.m_shuttle+mass(1,1) mass(1,2); mass(1,2) mass(2,2)]
    @assert isposdef(Symmetric(M))
    # Graded panel endpoints; Simpson midpoints are arithmetic in physical y.
    z=range(0.0,log1p(p.Leff/(p.h_eff/max(alpha,0.001)));length=panels+1)
    edges=p.h_eff/max(alpha,0.001).*expm1.(z)
    edges[end]=p.Leff
    y=zeros(2*panels+1); weights=zero(y)
    for j in 1:panels
        i=2*j-1; a,b=edges[j],edges[j+1]; d=b-a
        y[i]=a; y[i+1]=(a+b)/2; y[i+2]=b
        weights[i]+=d/6; weights[i+1]+=2*d/3; weights[i+2]+=d/6
    end
    B2=[shape(p,ke,p.Lf-v) for v in y]; B1=1 .- B2
    Cstruct=[p.c1 0.0;0.0 0.0]+p.n_beams*p.ce*[1.0 -1.0;-1.0 1.0]
    Model(p,ke,k1,k3,kss,p.g0-2*p.Tp-p.h_eff,alpha,2*p.Tp/p.epsr,
        6*p.slip_coefficient*p.mean_free_path,M,inv(M),M*ones(2),Cstruct,
        y,weights,B1,B2,panels)
end

function cumulative_simpson(m,f)
    H=zeros(length(f))
    for j in 1:m.panels
        i=2*j-1; d=m.y[i+2]-m.y[i]
        H[i+1]=H[i]+d*(5*f[i]+8*f[i+1]-f[i+2])/24
        H[i+2]=H[i]+d*(f[i]+4*f[i+1]+f[i+2])/6
    end
    H
end

"""Total capacitance, exact coordinate gradient, and passive generalized film matrix."""
function constitutive(m,x1,x2;film=true)
    p=m.p; C=p.cp; c1=0.0; c2=0.0; D=zeros(2,2)
    for r in (-1.0,1.0)
        d=m.gc .+ m.alpha.*m.y .- r.*(m.B1.*x1 .+ m.B2.*x2)
        h=p.h_eff .+ softpos.(d,p.eps_gap)
        dh=dsoftpos.(d,p.eps_gap)
        h1=-r.*dh.*m.B1; h2=-r.*dh.*m.B2
        C+=p.n_beams*p.eps0*p.Tf*dot(m.weights,1 ./ (h .+ m.hd))
        c1-=p.n_beams*p.eps0*p.Tf*dot(m.weights,h1./(h .+ m.hd).^2)
        c2-=p.n_beams*p.eps0*p.Tf*dot(m.weights,h2./(h .+ m.hd).^2)
        if film
            H1=cumulative_simpson(m,h1); H2=cumulative_simpson(m,h2)
            W=m.weights./(h.^2 .*(h .+ m.kp)); I0=sum(W)
            mu1=dot(W,H1)/I0; mu2=dot(W,H2)/I0
            Z1=H1 .- mu1; Z2=H2 .- mu2
            chi=p.seal_at_contact ? smootherstep((r*x2-m.gc+p.seal_width)/(2*p.seal_width)) : 0.0
            fac=12*p.eta*p.Tf*p.n_beams*p.film_scale
            D[1,1]+=fac*(dot(W,Z1.^2)+chi*I0*mu1^2)
            D[1,2]+=fac*(dot(W,Z1.*Z2)+chi*I0*mu1*mu2)
            D[2,2]+=fac*(dot(W,Z2.^2)+chi*I0*mu2^2)
        end
    end
    D[2,1]=D[1,2]
    (;C,grad=[c1,c2],D)
end

function springs(m,x1,x2)
    p=m.p; rel=x2-x1
    U=m.k1*x1^2/2+m.k3*x1^4/4+p.n_beams*m.ke*rel^2/2
    F=[-m.k1*x1-m.k3*x1^3+p.n_beams*m.ke*rel,-p.n_beams*m.ke*rel]
    for r in (-1.0,1.0)
        s=softpos(r*x1-p.gss,p.eps_stop)
        U+=m.kss*s^2/2
        F[1]-=r*m.kss*s*dsoftpos(r*x1-p.gss,p.eps_stop)
        if isfinite(p.ghs)
            z=r*x1-p.ghs; s=softpos(z,p.eps_stop)
            U+=p.khs*s^(p.phs+1)/(p.phs+1)
            F[1]-=r*p.khs*s^p.phs*dsoftpos(z,p.eps_stop)
        end
    end
    (;U,F)
end

function wall(m,x2,v2)
    p=m.p; U=0.0; F=0.0; loss=0.0
    for r in (-1.0,1.0)
        z=r*x2-m.gc; vr=r*v2; s=softpos(z,p.eps_wall)
        A=p.kw*s^p.pw*dsoftpos(z,p.eps_wall)
        gate=max(1+p.cw*vr,0.0)
        U+=p.n_beams*p.kw*s^(p.pw+1)/(p.pw+1)
        F-=p.n_beams*r*A*gate
        # Includes the unloading cutoff: loss=-A*vr when gate=0.
        loss+=p.n_beams*A*vr*(gate-1)
    end
    (;U,F,loss)
end

"""State [x1,x2,v1,v2,Vout,Wbase,Wbias,ER,Dfilm,Dstruct,Dwall,throughput]."""
function rhs!(du,u,context,t)
    m,acceleration=context; p=m.p
    x1,x2,v1,v2,Vout=u[1:5]; v=[v1,v2]; Vc=p.Vbias-Vout
    c=constitutive(m,x1,x2); s=springs(m,x1,x2); w=wall(m,x2,v2)
    a=acceleration(t)
    force=s.F+[0.0,w.F]+0.5*Vc^2*c.grad-(m.Cstruct+c.D)*v-m.beta*a
    dv=m.Minv*force
    du[1]=v1;du[2]=v2;du[3]=dv[1];du[4]=dv[2]
    du[5]=-Vout/(p.Rload*c.C)+Vc/c.C*dot(c.grad,v)
    if length(du)>=12
        Pb=-a*dot(m.beta,v); Pe=p.Vbias*Vout/p.Rload
        PR=Vout^2/p.Rload; Pf=dot(v,c.D*v); Ps=dot(v,m.Cstruct*v)
        du[6]=Pb;du[7]=Pe;du[8]=PR;du[9]=Pf;du[10]=Ps;du[11]=w.loss
        du[12]=abs(Pb)+abs(Pe)+PR+Pf+Ps+w.loss
    end
    nothing
end

function energy(m,u)
    c=constitutive(m,u[1],u[2];film=false)
    dot(u[3:4],m.M*u[3:4])/2+springs(m,u[1],u[2]).U+
        wall(m,u[2],u[4]).U+c.C*(m.p.Vbias-u[5])^2/2
end
ledger(m,u,E0)=energy(m,u)-E0-u[6]-u[7]+sum(u[8:11])
function writecsv(path,header,rows)
    open(path,"w") do io
        println(io,join(header,","))
        for row in rows; println(io,join(row,",")); end
    end
end

# Independent audit of the old open-open wedge formula. These use the old
# alpha=(wb-wt)/Leff ONLY to reproduce the previously identified discrepancy.
function wedge_reference(p,h)
    a=(p.wb-p.wt)/p.Leff; kp=6*p.slip_coefficient*p.mean_free_path
    I(k)=quad(z->begin y=h*expm1(z)/a; hh=h*exp(z)
        y^k/(a*hh*(hh+kp)) end,0.0,log1p(a*p.Leff/h))
    12*p.eta*p.Tf*(I(2)-I(1)^2/I(0))
end
function wedge_legacy(p,h;n=8001)
    a=(p.wb-p.wt)/p.Leff; kp=6*p.slip_coefficient*p.mean_free_path
    y=range(0.0,p.Leff;length=n); dy=step(y)
    G=(h .+ a.*y).^2 .*(h .+ a.*y .+ kp)
    cc=12*p.eta*sum(y./G)/sum(1 ./ G)
    dp=(cc .-12*p.eta.*y)./G; pressure=0.0; result=0.0
    for j in 2:n
        pressure+=(dp[j-1]+dp[j])*dy/2; result+=pressure*dy
    end
    abs(p.Tf*result)
end
function restitution(vin,cw)
    cw==0 && return 1.0
    a=cw*vin; target=a-log1p(a); lo=0.0; hi=1-eps()
    for j in 1:100
        b=(lo+hi)/2
        if -b-log1p(-b)>target;hi=b;else;lo=b;end
    end
    (lo+hi)/(2*a)
end

function verify(;outdir=joinpath(@__DIR__,"results"))
    mkpath(outdir); m=Model(); p=m.p; mr=Model(p;panels=1024)
    massphys=p.rho*p.Tf*p.Lf*(p.wb+p.wt)/2
    frac=m.M[2,2]/p.n_beams/massphys
    caperr=0.0; filmerr=0.0; energyerr=0.0
    states=[(0.0,0.0),(m.gc-1e-6,m.gc-1e-6),
        (m.gc+0.1e-6,m.gc-2e-9),(m.gc+0.2e-6,m.gc+5e-9),
        (-m.gc-0.1e-6,-m.gc-2e-9)]
    ts=@testset "Corrected MEMS constitutive and energy checks" begin
        @test isapprox(m.ke,22.5969569824;rtol=1e-9)
        @test isapprox(m.k1,3.2799375;rtol=1e-12)
        @test isapprox(m.k3,3.035714285714e8;rtol=1e-11)
        @test isapprox(m.kss,5.831;rtol=1e-12)
        @test isapprox(m.gc,13.71e-6;rtol=1e-12)
        @test isapprox(frac,0.369820966146;rtol=1e-9)
        @test isposdef(Symmetric(m.M))
        @test isapprox(sum(m.M),p.m_shuttle+p.n_beams*massphys;rtol=1e-12)
        @test isapprox(shape(p,m.ke,p.Lf),1.0;atol=1e-12)
        @test shape(p,m.ke,0.0)==0.0
        @test softpos(-1.0,1e-12)>0
        @test dsoftpos(-1.0,1e-12)>0
        @test isapprox(sum(m.weights),p.Leff;rtol=1e-13)
        @test maximum(abs.(cumulative_simpson(m,ones(length(m.y)))-m.y))<1e-15
        for (x1,x2) in states
            c=constitutive(m,x1,x2); cr=constitutive(mr,x1,x2)
            @test c.C>p.cp
            @test eigmin(Symmetric(c.D))>=-1e-14*norm(c.D)
            de=norm(c.D-cr.D)/max(norm(cr.D),1e-30); filmerr=max(filmerr,de)
            @test de<1e-5
            dx=1e-12
            fd=[(constitutive(m,x1+dx,x2;film=false).C-constitutive(m,x1-dx,x2;film=false).C)/(2*dx),
                (constitutive(m,x1,x2+dx;film=false).C-constitutive(m,x1,x2-dx;film=false).C)/(2*dx)]
            er=norm(fd-c.grad)/max(norm(c.grad),1e-18); caperr=max(caperr,er)
            @test er<2e-5
            u=[x1,x2,0.003,-0.005,0.2,zeros(7)...]; du=zero(u)
            rhs!(du,u,(m,t->2.0),0.0)
            scales=[m.gc,m.gc,0.02,0.02,3.0]; g=zeros(5)
            for j in 1:5
                h=scales[j]*1e-6; up=copy(u);um=copy(u);up[j]+=h;um[j]-=h
                g[j]=(energy(m,up)-energy(m,um))/(2*h)
            end
            exact=du[6]+du[7]-sum(du[8:11])
            er=abs(dot(g,du[1:5])-exact)/max(du[12],1e-25);energyerr=max(energyerr,er)
            @test er<5e-5
        end
        @test wall(m,m.gc+10e-9,-0.03).loss>0
        @test wall(m,m.gc+10e-9,0.03).loss>0
        @test isapprox(wedge_reference(p,50e-9),1.15476443946e-4;rtol=1e-8)
        @test isapprox(wedge_legacy(p,50e-9),2.84470582048e-5;rtol=1e-8)
        @test isapprox(restitution(0.005,50.0),0.8568510175;rtol=1e-8)
        @test isapprox(constitutive(m,1e-6,2e-6).C,constitutive(m,-1e-6,-2e-6).C;rtol=1e-14)
        @test norm(constitutive(m,1e-6,2e-6).grad+constitutive(m,-1e-6,-2e-6).grad)<1e-20
        @test norm(constitutive(m,1e-6,2e-6).D-constitutive(m,-1e-6,-2e-6).D)<1e-15
        @test norm(constitutive(m,0.0,0.0).grad)<1e-20
    end
    open(joinpath(outdir,"verification_summary.txt"),"w") do io
        println(io,"Julia ",VERSION,"; all native verification assertions passed.")
        for (name,value) in ("ke_N_per_m"=>m.ke,"k1_N_per_m"=>m.k1,"k3_N_per_m3"=>m.k3,
            "kss_N_per_m"=>m.kss,"contact_travel_m"=>m.gc,"shape_mass_fraction"=>frac,
            "max_cap_gradient_relative_error"=>caperr,"max_film_refinement_relative_error"=>filmerr,
            "max_point_energy_relative_error"=>energyerr)
            println(io,name," = ",value)
        end
        println(io,"Numerical consistency is not experimental validation.")
    end
    rows=[]
    for h in exp.(range(log(1e-9),log(20e-6);length=70))
        a=(p.wb-p.wt)/p.Leff; K=p.eps0*p.Tf; cair=K/a*log1p(a*p.Leff/h)
        Cfilm=p.eps0*p.epsr*p.Leff*p.Tf/p.Tp
        lump=1/(2/Cfilm+1/cair)
        gradold=(lump/cair)^2*K/a*(1/h-1/(h+a*p.Leff))
        gradnew=K/a*(1/(h+m.hd)-1/(h+m.hd+a*p.Leff))
        push!(rows,(h*1e9,wedge_legacy(p,h),wedge_reference(p,h),gradold/gradnew))
    end
    writecsv(joinpath(outdir,"constitutive_audit.csv"),["gap_nm","legacy_film","corrected_film","force_ratio"],rows)
    rows=[]
    for z in range(-1e-6,1e-6;length=201)
        x1=m.gc+z; x2=z<=0 ? x1 : m.gc
        c=constitutive(m,x1,x2); v=z<=0 ? [1.0,1.0] : [1.0,0.0]
        push!(rows,(z*1e6,c.C*1e12,0.5*p.Vbias^2*c.grad[1]*1e6,
            0.5*p.Vbias^2*c.grad[2]*1e6,dot(v,c.D*v)))
    end
    writecsv(joinpath(outdir,"pinned_tip_path.csv"),["travel_um","C_pF","Q1_uN","Q2_uN","D_path"],rows)
    println("Verification outputs written to ",outdir)
    (;caperr,filmerr,energyerr)
end

# Solvers are loaded only for simulations; --verify needs Julia stdlibs only.
function simulate(;kwargs...)
    @eval import SciMLBase, OrdinaryDiffEqRosenbrock, ADTypes
    Base.invokelatest(_simulate;kwargs...)
end
function _simulate(;p=Params(),panels=512,kind=:probe,cycles=8,freq=20.0,
        acceleration=4.95*9.80665,reltol=1e-7,abstol=1e-10,
        outdir=joinpath(@__DIR__,"results"),tag=string(kind))
    @assert kind in (:probe,:drive) && cycles>4 && freq>0
    m=Model(p;panels); mkpath(outdir)
    xs=m.gc; vs=xs*sqrt(m.k1/sum(m.M)); Es=m.k1*xs^2
    scale=[xs,xs,vs,vs,max(abs(p.Vbias),1.0),fill(Es,7)...]
    u0=zeros(12)
    if kind==:probe
        u0[1:5]=[m.gc+0.15e-6,m.gc-30e-9,0.004,0.004,0.0]
    end
    ramp(t)=0.5*(1-cos(pi*min(t*freq/4,1.0)))
    accel=kind==:probe ? (t->0.0) : (t->acceleration*ramp(t)*sin(2*pi*freq*t))
    tend=kind==:probe ? 600e-6 : cycles/freq
    function scaled_rhs!(dz,z,context,t)
        rhs!(dz,z.*scale,context,t); dz ./= scale; nothing
    end
    prob=SciMLBase.ODEProblem(scaled_rhs!,u0./scale,(0.0,tend),(m,accel))
    sol=SciMLBase.solve(prob,OrdinaryDiffEqRosenbrock.Rodas5P(autodiff=ADTypes.AutoFiniteDiff());
        reltol,abstol,dtmax=kind==:probe ? 1e-6 : 2e-5,maxiters=10^7,
        save_everystep=true,dense=true)
    @assert SciMLBase.successful_retcode(sol) "Solver failed: $(sol.retcode)"
    @assert isapprox(sol.t[end],tend;rtol=1e-12) "Incomplete integration"
    physical(t)=sol(t).*scale
    E0=energy(m,u0); last=physical(tend)
    maxres=maximum(abs(ledger(m,z.*scale,E0)) for z in sol.u)
    # Independent Gauss integration of resistor power over each accepted step.
    gx,gw=gausslegendre(8); ERquad=0.0
    for j in 2:length(sol.t)
        a,b=sol.t[j-1],sol.t[j]
        ERquad+=(b-a)/2*sum(gw[k]*physical((a+b)/2+(b-a)/2*gx[k])[5]^2/p.Rload for k in eachindex(gx))
    end
    throughput=last[12]; relledger=maxres/max(throughput,eps()*Es)
    erquad=abs(ERquad-last[8])/max(abs(last[8]),eps()*Es)
    @assert relledger<1e-5 "Energy balance failed"
    @assert erquad<5e-4 "Resistor quadrature failed"
    rows=[]
    times=range(0.0,tend;length=kind==:probe ? 2401 : 8001)
    for t in times
        u=physical(t)
        push!(rows,(t,t*1e6,u[1]*1e6,u[2]*1e6,u[5]*1e3,
            u[8],ledger(m,u,E0),max(u[2]-m.gc,-u[2]-m.gc)*1e9))
    end
    writecsv(joinpath(outdir,tag*"_timeseries.csv"),
        ["t_s","t_us","x1_um","x2_um","Vout_mV","ER_J","energy_residual_J","tip_overlap_nm"],rows)
    metrics=Dict("ER_J"=>last[8],"ER_quadrature_J"=>ERquad,"ER_quadrature_relative_error"=>erquad,
        "max_energy_residual_J"=>maxres,"energy_residual_over_throughput"=>relledger,
        "saved_steps"=>length(sol.t),"duration_s"=>tend,"film_panels"=>panels,
        "reltol"=>reltol,"scaled_abstol"=>abstol)
    if kind==:drive
        u=physical(tend); prev=physical(tend-1/freq)
        metrics["last_cycle_scaled_recurrence_error"]=norm((u[1:5]-prev[1:5])./scale[1:5])
        metrics["last_four_cycles_mean_power_W"]=(u[8]-physical(tend-4/freq)[8])/(4/freq)
        # Fine window around the last sampled sign change, if a contact occurs.
        us=[z.*scale for z in sol.u]
        k=findlast(j->(abs(us[j][2])-m.gc)*(abs(us[j-1][2])-m.gc)<0,2:length(us))
        if k!==nothing
            tcontact=sol.t[k+1]; win=range(max(0,tcontact-100e-6),min(tend,tcontact+300e-6);length=2001)
            writecsv(joinpath(outdir,tag*"_contact_window.csv"),["t_s","x1_um","x2_um","Vout_mV"],
                [(t,physical(t)[1]*1e6,physical(t)[2]*1e6,physical(t)[5]*1e3) for t in win])
        end
    end
    open(joinpath(outdir,tag*"_summary.txt"),"w") do io
        println(io,"Julia ",VERSION,"; solver retcode = ",sol.retcode)
        for k in sort(collect(keys(metrics)));println(io,k," = ",metrics[k]);end
        println(io,"Finite transient; neither a periodic-attractor proof nor experimental validation.")
    end
    println(tag,": ",sol.retcode,"; energy residual / throughput = ",relledger)
    (;metrics,sol,scale,model=m)
end
function convergence(;outdir=joinpath(@__DIR__,"results"))
    a=simulate(;outdir,tag="probe")
    b=simulate(;outdir,tag="probe_refined",reltol=2e-8,abstol=2e-11)
    c=simulate(;outdir,tag="probe_grid_refined",reltol=2e-8,abstol=2e-11,panels=1024)
    d1=abs(a.metrics["ER_J"]/b.metrics["ER_J"]-1)
    d2=abs(b.metrics["ER_J"]/c.metrics["ER_J"]-1)
    @assert max(d1,d2)<0.01 "Probe load energy is not converged to 1 percent"
    open(joinpath(outdir,"probe_convergence.txt"),"w") do io
        println(io,"ER relative change, tighter tolerances = ",d1)
        println(io,"ER relative change, doubled film panels = ",d2)
    end
    latexnum(v)=begin s=@sprintf("%.3e",v); a,b=split(s,"e");a*"\\times10^{"*string(parse(Int,b))*"}" end
    open(joinpath(outdir,"probe_summary.tex"),"w") do io
        println(io,"\\begin{tabular}{lr}\\toprule Quantity & Value ","\\\\","\\midrule")
        for (label,v) in [("Integrated load energy (J)",a.metrics["ER_J"]),
            ("Energy residual / power throughput",a.metrics["energy_residual_over_throughput"]),
            ("Independent resistor quadrature error",a.metrics["ER_quadrature_relative_error"]),
            ("Load-energy change: solver refinement",d1),("Load-energy change: film refinement",d2)]
            println(io,label," & \$",latexnum(v),"\$ ","\\\\")
        end
        println(io,"\\bottomrule\\end{tabular}")
    end
    (;tolerance_change=d1,grid_change=d2)
end

function main(args=ARGS)
    mode=isempty(args) ? "--verify" : args[1]
    if mode=="--verify";verify()
    elseif mode=="--probe";simulate()
    elseif mode=="--convergence";convergence()
    elseif mode=="--drive";simulate(;kind=:drive,cycles=30)
    elseif mode=="--help"
        println("--verify | --probe | --convergence | --drive\nSee README.md for installation and parameter changes.")
    else;error("Unknown mode: $mode")
    end
end
end # module
if abspath(PROGRAM_FILE)==@__FILE__
    CorrectedMEMS.main()
end
