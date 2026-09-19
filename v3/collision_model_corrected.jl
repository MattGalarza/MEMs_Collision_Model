#!/usr/bin/env julia
# Corrected trapezoidal-electrode MEMS model, SI units throughout.
# Julia-only implementation and numerical verification. See the companion TeX.
# Run: julia collision_model_corrected.jl --verify
# Default: julia --project=. collision_model_corrected.jl (10-cycle drive + plots)
# More: --probe | --drive --cycles 30 | --convergence | --help 
module CorrectedMEMS
using LinearAlgebra, Printf, Test, Dates, TOML
export Params, Model, PlotOptions, constitutive, rhs!, energy, verify, simulate, convergence, makeplots

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
    @eval import SciMLBase, OrdinaryDiffEqRosenbrock, ADTypes, CairoMakie
    Base.invokelatest(_simulate;kwargs...)
end
function _simulate(;p=Params(),panels=512,kind=:probe,cycles=10,freq=20.0,
        acceleration=1.9*9.80665,reltol=1e-7,abstol=1e-10,
        outdir=joinpath(@__DIR__,"results"),tag=string(kind),plot_options=PlotOptions())
    @assert kind in (:probe,:drive) && cycles>4 && freq>0
    @assert occursin(r"^[A-Za-z0-9_-]+$",tag) "Use letters, digits, underscores or hyphens in tag"
    # A new directory on EVERY simulation prevents earlier runs being overwritten.
    root=abspath(outdir); mkpath(root)
    runname=tag*"__"*Dates.format(now(),"yyyymmdd_HHMMSS_sss")
    outdir=joinpath(root,runname); suffix=1
    while ispath(outdir)
        outdir=joinpath(root,runname*"_"*string(suffix)); suffix+=1
    end
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
    println("Simulating ",tag," for ",tend," s; output: ",outdir);flush(stdout)
    sol=SciMLBase.solve(prob,OrdinaryDiffEqRosenbrock.Rodas5P(autodiff=ADTypes.AutoFiniteDiff());
        reltol,abstol,dtmax=kind==:probe ? 1e-6 : 2e-5,maxiters=10^7,
        save_everystep=true,dense=true)
    @assert SciMLBase.successful_retcode(sol) "Solver failed: $(sol.retcode)"
    @assert isapprox(sol.t[end],tend;rtol=1e-12) "Incomplete integration"
    println("Integration complete: ",length(sol.t)," accepted states. Checking energy and generating plots.");flush(stdout)
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
    run=(;metrics,sol,scale,model=m,kind,freq,acceleration=accel,
        acceleration_amplitude=kind==:drive ? acceleration : 0.0,tag,outdir)
    println("Generating automatic review plots in ",outdir)
    plots=makeplots(run;options=plot_options)
    (;run...,plots)
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

# ---------------------------------------------------------------------------
# Automatic review figures. Constitutive/RHS equations above are unchanged.
# All numerical calculations and rendering in this deliverable use Julia.
# ---------------------------------------------------------------------------
Base.@kwdef struct PlotOptions
    full_points::Int = 6001
    tail_points::Int = 6001
    render_bins::Int = 1400
    png_scale::Float64 = 2.0
    recurrence_tol::Float64 = 1e-4
    sweep_points::Int = 601
    surface_positions::Int = 181
    surface_velocities::Int = 61
    sweep_limit::Union{Nothing,Float64} = nothing # metres; symmetric translation
    sweep_voltage::Union{Nothing,Float64} = nothing # fixed capacitor voltage, V
    sweep_speed::Union{Nothing,Float64} = nothing # signed-velocity range +/- v, m/s
    bending_limit::Union{Nothing,Float64} = nothing # relative displacement, m
end

const FORCE_NAMES = ["Suspension: linear", "Suspension: cubic", "Electrode bending",
    "Soft stopper", "Hard stopper", "Contact", "Electrostatic", "Fluid",
    "Structural damping", "Base excitation", "Net applied", "Force-balance residual"]
const FORCE_KEYS = ["suspension_linear", "suspension_cubic", "electrode_bending",
    "soft_stop", "hard_stop", "contact", "electrostatic", "fluid",
    "structural_damping", "base", "net", "balance_residual"]
const DIAG_KEYS = ["C_F", "Vc_V", "charge_C", "current_A", "kinetic_J",
    "spring_J", "wall_J", "electric_J", "total_J", "residual_J",
    "Pbase_W", "Pbias_W", "Presistor_W", "Pfluid_W", "Pstruct_W", "Pwall_W",
    "hmin_plus_m", "hmin_minus_m", "overlap_plus_m", "overlap_minus_m",
    "chi_plus", "chi_minus", "D11_Ns_m", "D12_Ns_m", "D22_Ns_m",
    "Dmineig_Ns_m", "Kn_max", "root_accel_m_s2", "tip_accel_m_s2",
    "wall_elastic_N", "wall_rate_N", "tip_not_min_flag"]
const REVIEW_COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#333333"]

"""Exact decomposition of the forces used by rhs!, including both coordinates.
Columns 1:10 sum to column 11. Column 12 checks against M*acceleration from rhs!.
Damping components can transfer energy between coordinates: only v'Dv is passive.
"""
function review_snapshot(m,u,a,E0)
    p=m.p; x1,x2,v1,v2,Vo=u[1:5]; v=[v1,v2]; rel=x2-x1
    c=constitutive(m,x1,x2); s=springs(m,x1,x2); w=wall(m,x2,v2)
    we=wall(m,x2,0.0); Vc=p.Vbias-Vo
    F=zeros(2,12)
    F[1,1]=-m.k1*x1; F[1,2]=-m.k3*x1^3
    F[:,3]=p.n_beams*m.ke*rel.*[1.0,-1.0]
    for r in (-1.0,1.0)
        z=r*x1-p.gss; sp=softpos(z,p.eps_stop)
        F[1,4]-=r*m.kss*sp*dsoftpos(z,p.eps_stop)
        if isfinite(p.ghs)
            z=r*x1-p.ghs; sp=softpos(z,p.eps_stop)
            F[1,5]-=r*p.khs*sp^p.phs*dsoftpos(z,p.eps_stop)
        end
    end
    F[2,6]=w.F; F[:,7]=0.5*Vc^2*c.grad
    F[:,8]=-c.D*v; F[:,9]=-m.Cstruct*v; F[:,10]=-m.beta*a
    F[:,11]=vec(sum(F[:,1:10];dims=2))
    # Independent grouping identical to RHS, avoiding a second expensive film evaluation.
    rhsforce=s.F+[0.0,w.F]+0.5*Vc^2*c.grad-(m.Cstruct+c.D)*v-m.beta*a
    acc=m.Minv*rhsforce
    F[:,12]=F[:,11]-m.M*acc
    hmin=zeros(2); non_tip=false
    for (i,r) in enumerate((1.0,-1.0))
        raw=m.gc .+m.alpha.*m.y .-r.*(m.B1.*x1.+m.B2.*x2)
        hmin[i]=p.h_eff+softpos(minimum(raw),p.eps_gap)
        non_tip |= minimum(raw)<raw[1]-max(1e-12,0.01*p.eps_gap)
    end
    chi(r)=p.seal_at_contact ? smootherstep((r*x2-m.gc+p.seal_width)/(2*p.seal_width)) : 0.0
    K=dot(v,m.M*v)/2; Ee=c.C*Vc^2/2; Et=K+s.U+w.U+Ee
    residual=length(u)>=12 ? Et-E0-u[6]-u[7]+sum(u[8:11]) : NaN
    diag=[c.C,Vc,c.C*Vc,Vo/p.Rload,K,s.U,w.U,Ee,Et,residual,
        -a*dot(m.beta,v),p.Vbias*Vo/p.Rload,Vo^2/p.Rload,dot(v,c.D*v),
        dot(v,m.Cstruct*v),w.loss,hmin[1],hmin[2],x2-m.gc,-x2-m.gc,
        chi(1.0),chi(-1.0),c.D[1,1],c.D[1,2],c.D[2,2],
        eigmin(Symmetric(c.D)),p.mean_free_path/minimum(hmin),acc[1],acc[2],
        we.F,w.F-we.F,Float64(non_tip)]
    (;F,diag)
end

# Keep accepted solver steps as well as a uniform grid. Fast contact events must
# not disappear simply because the plotted uniform grid is coarser than contact.
function review_trace(run,times)
    t=sort!(unique!(Float64.(times))); N=length(t)
    u=Matrix{Float64}(undef,12,N); F=Array{Float64}(undef,2,12,N)
    d=Matrix{Float64}(undef,length(DIAG_KEYS),N)
    E0=energy(run.model,run.sol(run.sol.t[1]).*run.scale)
    for j in eachindex(t)
        u[:,j]=run.sol(t[j]).*run.scale
        q=review_snapshot(run.model,view(u,:,j),run.acceleration(t[j]),E0)
        F[:,:,j]=q.F; d[:,j]=q.diag
    end
    (;t,u,F,d)
end

# Extrema-preserving display reduction, applied separately to each curve.
# Physical diagnostics and cycle statistics use the unreduced trace/solution.
function review_indices(t,ys;bins=1400)
    N=length(t); N<=4*bins && return collect(1:N)
    keep=falses(N); keep[1]=keep[N]=true
    edges=range(t[1],t[end];length=bins+1)
    for b in 1:bins
        a=searchsortedfirst(t,edges[b]); z=min(N,searchsortedlast(t,edges[b+1]))
        a>z && continue
        keep[a]=keep[z]=true
        for y in ys
            lo=hi=a
            for j in a+1:z
                y[j]<y[lo] && (lo=j)
                y[j]>y[hi] && (hi=j)
            end
            keep[lo]=keep[hi]=true
        end
    end
    findall(keep)
end

function review_lines!(ax,x,ys,labels;options=PlotOptions(),scale=1.0,colors=REVIEW_COLORS,legend_position=:rt)
    for (k,y) in enumerate(ys)
        idx=review_indices(x,[y];bins=options.render_bins)
        CairoMakie.lines!(ax,x[idx],y[idx].*scale;color=colors[mod1(k,length(colors))],
            linewidth=1.65,label=labels[k])
    end
    length(labels)>1 && CairoMakie.axislegend(ax;position=legend_position,labelsize=11,framevisible=false)
    ax
end

function review_figure(title,subtitle;size=(1200,920))
    f=CairoMakie.Figure(;size,fontsize=15,figure_padding=22)
    CairoMakie.Label(f[0,1:2],title*"\n"*subtitle;fontsize=18,halign=:left,tellwidth=false)
    f
end
review_axis(f,row,col,title,xlabel,ylabel;yscale=identity)=CairoMakie.Axis(f[row,col];
    title,xlabel,ylabel,yscale,titlesize=15,xlabelsize=13,ylabelsize=13,
    xticklabelsize=11,yticklabelsize=11,xgridcolor=(:gray,0.12),ygridcolor=(:gray,0.12))
function review_save(f,dir,name,caption,items,options)
    # Finish native graphics writes in a temporary directory, then copy complete
    # bytes. This also avoids incomplete exports on synchronized filesystems.
    mktempdir() do staging
        pdf=joinpath(staging,"figure.pdf");png=joinpath(staging,"figure.png")
        CairoMakie.save(pdf,f)
        CairoMakie.save(png,f;px_per_unit=options.png_scale)
        for (source,ext) in ((pdf,".pdf"),(png,".png"))
            bytes=read(source)
            isempty(bytes) && error("Empty graphics export: $name$ext")
            write(joinpath(dir,name*ext),bytes)
        end
    end
    push!(items,(name=name,caption=caption))
    if length(items)%8==0
        println("  Saved ",length(items)," figure pairs so far.");flush(stdout)
    end
    nothing
end

function review_window_plots(run,tr,idx,label,suffix,dir,items,options)
    t=tr.t[idx]; offset=t[1]; short_window=run.kind==:probe || t[end]-t[1]<1e-3
    unit=short_window ? 1e6 : 1.0
    tx=(t.-offset).*unit; xlabel=short_window ? "Time from window start (μs)" : "Time from window start (s)"
    subtitle=@sprintf("%s | actual t = %.6g to %.6g s",label,t[1],t[end])
    u=tr.u[:,idx]; d=tr.d[:,idx]; F=tr.F[:,:,idx]
    f=review_figure("State variables",subtitle)
    stateys=[u[1,:].*1e6,u[2,:].*1e6,u[3,:].*1e3,u[4,:].*1e3,u[5,:].*1e3,(u[2,:]-u[1,:]).*1e9]
    titles=["Shuttle displacement x₁","Electrode-tip displacement x₂","Shuttle velocity v₁",
        "Electrode-tip velocity v₂","Output voltage Vout","Relative bending x₂ − x₁"]
    units=["Displacement (μm)","Displacement (μm)","Velocity (mm/s)","Velocity (mm/s)","Voltage (mV)","Deflection (nm)"]
    for k in 1:6
        ax=review_axis(f,cld(k,2),mod1(k,2),titles[k],xlabel,units[k])
        review_lines!(ax,tx,[stateys[k]],[""];options)
        k==2 && CairoMakie.hlines!(ax,[-run.model.gc,run.model.gc].*1e6;color=(:black,0.5),linestyle=:dash,yautolimits=false)
    end
    review_save(f,dir,"states_"*suffix,"All five physical states and relative bending. Dashed tip lines mark nominal contact.",items,options)
    for coordinate in 1:2
        f=review_figure("Forces conjugate to x"*string(coordinate),subtitle*" | signed generalized forces; full device";size=(1200,1570))
        for k in 1:12
            ax=review_axis(f,cld(k,2),mod1(k,2),FORCE_NAMES[k],xlabel,k==12 ? "Residual (pN)" : "Force (μN)")
            review_lines!(ax,tx,[vec(F[coordinate,k,:])],[""];options,scale=k==12 ? 1e12 : 1e6)
        end
        review_save(f,dir,"forces_x$(coordinate)_"*suffix,"Each physical force, their sum, and the residual against M times acceleration. Zero panels identify inactive components.",items,options)
    end
    f=review_figure("Energy and instantaneous power",subtitle)
    specs=[("Stored mechanical energy",[5,6,7],["Kinetic","Springs","Contact"],1e12,"Energy (pJ)"),
        ("Stored electrical energy",[8],["Electrical"],1e12,"Energy (pJ)"),
        ("Input powers",[11,12],["Base","Bias source"],1e9,"Power (nW)"),
        ("Dissipated powers",[13,14,15,16],["Resistor","Fluid","Structural","Contact"],1e9,"Power (nW)"),
        ("Load power",[13],["Resistor"],1e12,"Power (pW)"),
        ("Energy-balance residual",[10],["Residual"],1e18,"Residual (aJ)")]
    for (k,(title,cols,labs,fac,ylab)) in enumerate(specs)
        ax=review_axis(f,cld(k,2),mod1(k,2),title,xlabel,ylab)
        review_lines!(ax,tx,[d[j,:] for j in cols],labs;options,scale=fac)
    end
    review_save(f,dir,"energy_power_"*suffix,"Mechanical/electrical storage, source power, every dissipation channel, and energy residual.",items,options)
    f=review_figure("Electrical variables and accelerations",subtitle)
    specs=[("Total capacitance",[1],["C"],1e12,"Capacitance (pF)"),
        ("Capacitor voltage",[2],["Vc"],1.0,"Voltage (V)"),
        ("Stored charge",[3],["Charge"],1e12,"Charge (pC)"),
        ("Load current",[4],["Current"],1e9,"Current (nA)"),
        ("Load power",[13],["Resistor"],1e12,"Power (pW)"),
        ("Mechanical accelerations",[28,29],["Shuttle","Tip"],1.0,"Acceleration (m/s²)")]
    for (k,(title,cols,labs,fac,ylab)) in enumerate(specs)
        ax=review_axis(f,cld(k,2),mod1(k,2),title,xlabel,ylab)
        review_lines!(ax,tx,[d[j,:] for j in cols],labs;options,scale=fac)
    end
    review_save(f,dir,"electrical_"*suffix,"Capacitance, capacitor voltage, charge, resistor current and power, and both mechanical accelerations.",items,options)
    f=review_figure("Integrated work and loss states",subtitle*" | cumulative values since simulation start";size=(1200,1110))
    titles=["Base work Wbase","Bias-source work Wbias","Load energy ER","Fluid loss Dfilm",
        "Structural loss Dstruct","Contact loss Dwall","Energy throughput","Total stored energy"]
    for k in 1:8
        ax=review_axis(f,cld(k,2),mod1(k,2),titles[k],xlabel,"Energy (pJ)")
        y=k<=7 ? u[k+5,:] : d[9,:]
        review_lines!(ax,tx,[y],[""];options,scale=1e12)
    end
    review_save(f,dir,"integrated_states_"*suffix,"All seven auxiliary integration states are included; these need not be periodic even when physical states recur.",items,options)
    f=review_figure("Contact and gas-film diagnostics",subtitle;size=(1200,1120))
    specs=[("Positive nominal tip overlap",[19,20],["+ side","− side"],1e9,"Contact deformation (nm)"),
        ("Minimum effective gap along electrode",[17,18],["+ side","− side"],1e9,"Gap (nm)"),
        ("Hydraulic sealing parameter",[21,22],["+ side","− side"],1.0,"χ (0 open; 1 sealed)"),
        ("Contact-force decomposition on x₂",[30,31],["Elastic","Rate contribution"],1e6,"Force (μN)"),
        ("Fluid matrix entries",[23,24,25],["D₁₁","D₁₂","D₂₂"],1e3,"Damping (mN s/m)"),
        ("Minimum fluid eigenvalue",[26],["Minimum eigenvalue"],1e3,"Damping (mN s/m)"),
        ("Maximum local Knudsen number",[27],["λ / min(h)"],1.0,"Kn"),
        ("Tip is not the minimum-gap point",[32],["Flag"],1.0,"Flag (0 no; 1 yes)")]
    for (k,(title,cols,labs,fac,ylab)) in enumerate(specs)
        ax=review_axis(f,cld(k,2),mod1(k,2),title,xlabel,ylab;yscale=k==2 ? log10 : identity)
        ys=k==1 ? [max.(d[j,:],0.0) for j in cols] : [d[j,:] for j in cols]
        review_lines!(ax,tx,ys,labs;options,scale=fac)
        k in (1,6) && CairoMakie.hlines!(ax,[0.0];color=(:black,0.4),linestyle=:dash)
    end
    review_save(f,dir,"contact_fluid_"*suffix,"Positive overlap is effective contact deformation. Kn and the minimum-gap flag are model-validity diagnostics, not validation results.",items,options)
    f=review_figure("Phase portraits and electrical response",subtitle;size=(1200,800))
    data=[(u[1,:].*1e6,u[3,:].*1e3,"Shuttle phase portrait","x₁ (μm)","v₁ (mm/s)"),
        (u[2,:].*1e6,u[4,:].*1e3,"Electrode phase portrait","x₂ (μm)","v₂ (mm/s)"),
        ((u[2,:]-u[1,:]).*1e9,vec(F[2,6,:]).*1e6,"Contact force versus relative bending","x₂ − x₁ (nm)","Qcontact,₂ (μN)"),
        (d[1,:].*1e12,d[3,:].*1e12,"Charge–capacitance trajectory","Capacitance (pF)","Charge (pC)")]
    for (k,(xx,yy,title,xlab,ylab)) in enumerate(data)
        ax=review_axis(f,cld(k,2),mod1(k,2),title,xlab,ylab)
        ii=review_indices(t,[xx,yy];bins=options.render_bins)
        CairoMakie.lines!(ax,xx[ii],yy[ii];color=REVIEW_COLORS[k],linewidth=1.2)
        CairoMakie.scatter!(ax,[xx[1],xx[end]],[yy[1],yy[end]];color=["#009E73","#D55E00"],markersize=8)
    end
    review_save(f,dir,"phase_"*suffix,"Green and orange markers identify the start and end of each trajectory window.",items,options)
end

function review_cycles(run,options,dir,items)
    run.kind==:drive || return (label="Final 20% of unforced transient; no drive cycles",period=0,errors=Float64[])
    T=1/run.freq; tend=run.sol.t[end]; n=floor(Int,tend/T+1e-8)
    # Include accepted steps within the comparison cycle to resolve rapid impacts.
    function mismatch(k,lag)
        a=(k-1)*T; b=k*T
        times=sort!(unique!(vcat(collect(range(a,b;length=1001)),run.sol.t[(run.sol.t.>=a).&(run.sol.t.<=b)])))
        maximum(maximum(abs.((run.sol(t)[1:5]-run.sol(t-lag*T)[1:5]))) for t in times)
    end
    e1=[k>=2 ? mismatch(k,1) : NaN for k in 1:n]
    e2=[k>=3 ? mismatch(k,2) : NaN for k in 1:n]
    # Three consecutive period-1 comparisons after the four-cycle forcing ramp.
    period=n>=8 && maximum(e1[end-2:end])<=options.recurrence_tol ? 1 :
        n>=8 && maximum(e2[end-1:end])<=options.recurrence_tol ? 2 : 0
    label=period==1 ? "Final two drive cycles | period-1 recurrence passed" :
        period==2 ? "Final two drive cycles | period-2 recurrence passed" :
        "Final two drive cycles | recurrence not established"
    avgP=[(run.sol(k*T)[8]-run.sol((k-1)*T)[8])*run.scale[8]/T for k in 1:n]
    rows=[(k,k*T,e1[k],e2[k],avgP[k]) for k in 1:n]
    writecsv(joinpath(dir,"cycle_metrics.csv"),["cycle","end_time_s","period1_scaled_error","period2_scaled_error","mean_load_power_W"],rows)
    f=review_figure("Cycle convergence and load output","Comparisons use the five scaled physical states throughout each cycle; threshold is a numerical criterion";size=(1200,800))
    ax=review_axis(f,1,1,"Period-1 and period-2 recurrence","Drive cycle","Maximum scaled difference";yscale=log10)
    for (e,lab,col) in [(e1,"One-cycle lag",REVIEW_COLORS[1]),(e2,"Two-cycle lag",REVIEW_COLORS[2])]
        kk=findall(isfinite,e)
        CairoMakie.scatterlines!(ax,kk,max.(e[kk],1e-16);label=lab,color=col,markersize=5)
    end
    CairoMakie.hlines!(ax,[options.recurrence_tol];color=:black,linestyle=:dash)
    CairoMakie.axislegend(ax;labelsize=11,framevisible=false)
    ax=review_axis(f,1,2,"Mean resistor power per drive cycle","Drive cycle","Power (pW)")
    CairoMakie.scatterlines!(ax,1:n,avgP.*1e12;color=REVIEW_COLORS[3],markersize=5)
    ax=review_axis(f,2,1,"Last two displacement waveforms","Phase within drive cycle","x₂ (μm)")
    ax2=review_axis(f,2,2,"Last two voltage waveforms","Phase within drive cycle","Vout (mV)")
    for (j,k) in enumerate(n-1:n)
        a=(k-1)*T; b=k*T
        ts=sort!(unique!(vcat(collect(range(a,b;length=options.tail_points)),run.sol.t[(run.sol.t.>=a).&(run.sol.t.<=b)])))
        uu=[run.sol(t).*run.scale for t in ts]; phase=(ts.-a)./T
        review_lines!(ax,phase,[[v[2]*1e6 for v in uu]],["Cycle $k"];options,colors=[REVIEW_COLORS[j]])
        review_lines!(ax2,phase,[[v[5]*1e3 for v in uu]],["Cycle $k"];options,colors=[REVIEW_COLORS[j]])
    end
    CairoMakie.axislegend(ax;labelsize=11,framevisible=false); CairoMakie.axislegend(ax2;labelsize=11,framevisible=false)
    review_save(f,dir,"cycle_convergence","Period-1/2 recurrence is checked after the forcing ramp. Passing this finite numerical check is not a Floquet-stability proof; higher-period responses are not automatically classified.",items,options)
    (;label,period,errors=e1)
end

function review_events(run,times,dir)
    m=run.model; events=Tuple{Float64,Int,String,Float64}[]
    for r in (-1,1)
        gap(t)=m.gc-r*run.sol(t)[2]*run.scale[2]
        for j in 2:length(times)
            a=times[j-1]; b=times[j]; ga=gap(a); gb=gap(b)
            ((ga*gb<0) || (gb==0 && ga!=0)) || continue
            for _ in 1:60
                b-a<=max(1e-13,8*eps(max(abs(a),abs(b)))) && break
                mid=(a+b)/2; gm=gap(mid)
                if signbit(gm)==signbit(ga); a=mid;ga=gm;else;b=mid;end
            end
            t=(a+b)/2; vr=r*run.sol(t)[4]*run.scale[4]
            label=abs(vr)<=1e-9*run.scale[4] ? "tangent_candidate" : vr>0 ? "entry" : "exit"
            push!(events,(t,r,label,vr))
        end
    end
    sort!(events;by=first)
    writecsv(joinpath(dir,"contact_events.csv"),["time_s","side","event","normal_velocity_m_s"],events)
    events
end

# Symmetric grid enriched near nominal contact so nm transitions remain resolved
# even when the complete travel is tens of micrometres.
function review_sweep_grid(m,limit,n)
    p=m.p
    extra=vcat([c+q for c in (-m.gc,m.gc) for q in
        vcat(collect(range(-4*p.seal_width,4*p.seal_width;length=101)),
             collect(range(-8*p.eps_gap,8*p.eps_gap;length=81)),
             collect(range(-8*p.eps_wall,8*p.eps_wall;length=41)))],[-limit,0.0,limit])
    sort!(unique!(vcat(collect(range(-limit,limit;length=n)),filter(x->abs(x)<=limit,extra))))
end

function review_force_maps(run,tr,dir,items,options)
    m=run.model;p=m.p
    observed=max(maximum(abs,tr.u[1,:]),maximum(abs,tr.u[2,:]))
    limit=isnothing(options.sweep_limit) ? max(observed,m.gc+max(4*p.seal_width,20*p.eps_wall)) : options.sweep_limit
    V=isnothing(options.sweep_voltage) ? abs(p.Vbias) : options.sweep_voltage
    speed=isnothing(options.sweep_speed) ? max(maximum(abs,tr.u[3:4,:]),1e-3) : options.sweep_speed
    bend=isnothing(options.bending_limit) ? max(maximum(abs,tr.u[2,:]-tr.u[1,:]),10*p.eps_gap) : options.bending_limit
    @assert limit>0 && speed>0 && bend>0 && isfinite(V)
    xx=review_sweep_grid(m,limit,options.sweep_points); X=xx.*1e6
    cc=[constitutive(m,x,x) for x in xx]
    es=reduce(hcat,[0.5*V^2*c.grad for c in cc]); C=[c.C for c in cc]
    f=review_figure("Electrostatic force across the travel range",@sprintf("Common translation: x₁ = x₂ = x; fixed capacitor voltage Vc = %.5g V; full device",V);size=(1200,980))
    ax=review_axis(f,1,1,"Signed generalized force","Displacement x (μm)","Force (μN)")
    review_lines!(ax,X,[es[1,:],es[2,:],vec(sum(es;dims=1))],["Q₁","Q₂","Q₁ + Q₂"];options,scale=1e6,legend_position=:lt)
    ax=review_axis(f,1,2,"Force magnitude","Displacement x (μm)","Magnitude (μN)")
    review_lines!(ax,X,[abs.(es[1,:]),abs.(es[2,:]),abs.(vec(sum(es;dims=1)))],["|Q₁|","|Q₂|","|Q₁ + Q₂|"];options,scale=1e6,legend_position=:ct)
    ax=review_axis(f,2,1,"Total capacitance","Displacement x (μm)","Capacitance (pF)")
    review_lines!(ax,X,[C],[""];options,scale=1e12)
    z=findall(x->abs(x-m.gc)<=4*p.seal_width,xx)
    ax=review_axis(f,2,2,"Positive-contact detail","x − gc (nm)","Force (μN)")
    if length(z)>1
        review_lines!(ax,(xx[z].-m.gc).*1e9,[es[1,z],es[2,z],vec(sum(es[:,z];dims=1))],["Q₁","Q₂","Q₁ + Q₂"];options,scale=1e6,legend_position=:lt)
        CairoMakie.vlines!(ax,[0.0];color=(:black,0.5),linestyle=:dash)
    else
        CairoMakie.text!(ax,0.5,0.5;text="Selected range excludes contact",space=:relative,align=(:center,:center))
    end
    review_save(f,dir,"map_electrostatic","Fixed-voltage constitutive slice, not the driven trajectory. Sum Q₁+Q₂ is conjugate to common translation; generalized force components are not separate electrode-side forces.",items,options)
    writecsv(joinpath(dir,"map_electrostatic.csv"),["x1_m","x2_m","Vc_V","C_F","Q1_N","Q2_N","Qtranslation_N"],
        [(x,x,V,C[j],es[1,j],es[2,j],sum(es[:,j])) for (j,x) in enumerate(xx)])
    # Suspension and stoppers act on x1; use enough travel to show soft-stop onset.
    ls=max(limit,p.gss+max(4*p.seal_width,0.02*p.gss)); sx=review_sweep_grid(m,ls,options.sweep_points)
    sy=zeros(4,length(sx))
    for (j,x) in enumerate(sx)
        sy[1,j]=-m.k1*x;sy[2,j]=-m.k3*x^3
        for r in (-1.0,1.0)
            z=r*x-p.gss; sy[3,j]-=r*m.kss*softpos(z,p.eps_stop)*dsoftpos(z,p.eps_stop)
            if isfinite(p.ghs)
                z=r*x-p.ghs;sy[4,j]-=r*p.khs*softpos(z,p.eps_stop)^p.phs*dsoftpos(z,p.eps_stop)
            end
        end
    end
    f=review_figure("Suspension and stopper force laws","Signed forces conjugate to x₁; sweep includes the soft-stop threshold";size=(1200,800))
    for (k,title) in enumerate(["Linear suspension","Cubic suspension","Soft stopper",isfinite(p.ghs) ? "Hard stopper" : "Hard stopper: disabled in this model"])
        ax=review_axis(f,cld(k,2),mod1(k,2),title,"x₁ (μm)","Force (μN)")
        review_lines!(ax,sx.*1e6,[sy[k,:]],[""];options,scale=1e6)
        k==3 && CairoMakie.vlines!(ax,[-p.gss,p.gss].*1e6;color=(:black,0.5),linestyle=:dash)
    end
    review_save(f,dir,"map_suspension_stoppers","Mechanical force curves use their actual coordinate. Any disabled hard stopper appears explicitly as zero.",items,options)
    writecsv(joinpath(dir,"map_suspension_stoppers.csv"),["x1_m","linear_N","cubic_N","soft_N","hard_N"],[(sx[j],sy[:,j]...) for j in eachindex(sx)])
    bx=collect(range(-bend,bend;length=options.sweep_points)); bv=-p.n_beams*m.ke.*bx
    f=review_figure("Electrode bending and structural damping","Baseline electrode bending is linear; force on the other coordinate is equal and opposite";size=(1200,800))
    ax=review_axis(f,1,1,"Electrode bending force","x₂ − x₁ (nm)","Force (μN)")
    review_lines!(ax,bx.*1e9,[-bv,bv],["Q₁","Q₂"];options,scale=1e6)
    ax=review_axis(f,1,2,"Stored electrode bending energy","x₂ − x₁ (nm)","Energy (pJ)")
    review_lines!(ax,bx.*1e9,[0.5*p.n_beams*m.ke.*bx.^2],[""];options,scale=1e12)
    vv=collect(range(-speed,speed;length=options.surface_velocities))
    for (k,bvec,title) in [(1,[1.0,1.0],"Structural damping: v₁ = v₂ = v"),(2,[0.0,1.0],"Structural damping: v₁ = 0, v₂ = v")]
        ff=-(m.Cstruct*bvec)*transpose(vv)
        ax=review_axis(f,2,k,title,"Signed velocity v (mm/s)","Force (μN)")
        review_lines!(ax,vv.*1e3,[ff[1,:],ff[2,:]],["Q₁","Q₂"];options,scale=1e6)
    end
    review_save(f,dir,"map_bending_structural_damping","Relative bending uses the observed deflection range (or a user override). Structural damping may be zero with the uncalibrated baseline coefficients.",items,options)
    writecsv(joinpath(dir,"map_bending.csv"),["relative_displacement_m","Q1_N","Q2_N"],[(bx[j],-bv[j],bv[j]) for j in eachindex(bx)])
    writecsv(joinpath(dir,"map_structural_damping.csv"),["velocity_m_s","Qtranslation1_N","Qtranslation2_N","Qrelative1_N","Qrelative2_N"],
        [(v,(-m.Cstruct*[v,v])...,(-m.Cstruct*[0.0,v])...) for v in vv])
    f=review_figure("Contact force across both travel limits",@sprintf("Full-device force on x₂; curves hold signed v₂ fixed at 0 and ±%.4g mm/s",speed*1e3);size=(1200,800))
    vals=[[wall(m,x,v).F for x in xx] for v in (0.0,-speed,speed)]
    ax=review_axis(f,1,1,"Signed contact force","x₂ (μm)","Force (μN)")
    review_lines!(ax,X,vals,["v₂ = 0","v₂ = −vref","v₂ = +vref"];options,scale=1e6)
    ax=review_axis(f,1,2,"Contact-force magnitude","x₂ (μm)","Magnitude (μN)")
    review_lines!(ax,X,[abs.(v) for v in vals],["v₂ = 0","v₂ = −vref","v₂ = +vref"];options,scale=1e6)
    for (col,r) in [(1,-1.0),(2,1.0)]
        jj=findall(x->abs(r*x-m.gc)<=4*p.seal_width,xx)
        ax=review_axis(f,2,col,r<0 ? "Negative-contact detail" : "Positive-contact detail","x₂ (μm)","Force (μN)")
        !isempty(jj) && review_lines!(ax,X[jj],[v[jj] for v in vals],["v₂ = 0","v₂ = −vref","v₂ = +vref"];options,scale=1e6)
    end
    review_save(f,dir,"map_contact","Loading and unloading depend on the normal velocity r*v₂. These are constitutive slices; their areas are not automatically a realized contact-cycle loss.",items,options)
    writecsv(joinpath(dir,"map_contact.csv"),["x2_m","Q_v0_N","Q_vminus_N","Q_vplus_N"],[(xx[j],(v[j] for v in vals)...) for j in eachindex(xx)])
    for (name,positions,bvec,xcoords,title) in [
        ("translation",review_sweep_grid(m,limit,options.surface_positions),[1.0,1.0],x->(x,x),"x₁ = x₂ = x; v₁ = v₂ = v"),
        ("bending",collect(range(-bend,bend;length=options.surface_positions)),[0.0,1.0],x->(0.0,x),"x₁ = 0; x₂ = Δ; v₁ = 0; v₂ = v")]
        coeff=reduce(hcat,[constitutive(m,xcoords(x)...).D*bvec for x in positions])
        surfaces=[-coeff[i,:]*transpose(vv) for i in 1:2]
        f=review_figure("Fluid-force surface: "*name,title*" | signed generalized force; dissipative model";size=(1400,750))
        xlabel=name=="translation" ? "Position x (μm)" : "Relative tip position Δ (μm)"
        for coordinate in 1:2
            ax=CairoMakie.Axis3(f[1,coordinate];title="Fluid force Q$(coordinate)",xlabel,ylabel="Signed velocity (mm/s)",
                zlabel="Force (μN)",azimuth=1.18pi,elevation=0.20pi,perspectiveness=0.35,
                xlabelsize=14,ylabelsize=14,zlabelsize=14,xticklabelsize=11,yticklabelsize=11,zticklabelsize=11)
            surf=surfaces[coordinate].*1e6; cr=max(maximum(abs,surf),1e-12)
            CairoMakie.surface!(ax,positions.*1e6,vv.*1e3,surf;colormap=:balance,colorrange=(-cr,cr),rasterize=2)
        end
        review_save(f,dir,"map_fluid_3d_"*name,"Velocity is signed. The two surfaces are coordinate slices of Qfluid = −D(x1,x2)*[v1,v2], not a unique law of x₂ and v₂ alone.",items,options)
        if name=="translation"
            f=review_figure("Total fluid force during common translation","x₁ = x₂ = x; v₁ = v₂ = v; Qtranslation = Q₁ + Q₂";size=(1100,830))
            ax=CairoMakie.Axis3(f[1,1:2];xlabel="Position x (μm)",ylabel="Signed velocity v (mm/s)",
                zlabel="Fluid force (μN)",azimuth=1.18pi,elevation=0.20pi,perspectiveness=0.35,
                xlabelsize=15,ylabelsize=15,zlabelsize=15,xticklabelsize=12,yticklabelsize=12,zticklabelsize=12)
            z=(surfaces[1]+surfaces[2]).*1e6;cr=max(maximum(abs,z),1e-12)
            sp=CairoMakie.surface!(ax,positions.*1e6,vv.*1e3,z;colormap=:balance,colorrange=(-cr,cr),rasterize=2)
            CairoMakie.Colorbar(f[2,1:2],sp;vertical=false,label="Signed fluid force (μN)",height=14)
            review_save(f,dir,"map_fluid_3d_translation_total","Total force conjugate to common translation. Its sign opposes the common velocity; the surface is a fixed-kinematic slice, not a trajectory.",items,options)
        end
        writecsv(joinpath(dir,"map_fluid_3d_"*name*".csv"),["position_parameter_m","velocity_parameter_m_s","x1_m","x2_m","v1_m_s","v2_m_s","Q1_N","Q2_N","dissipation_W"],
            [(x,v,xcoords(x)...,(bvec.*v)...,surfaces[1][i,j],surfaces[2][i,j],v^2*dot(bvec,coeff[:,i])) for (i,x) in enumerate(positions) for (j,v) in enumerate(vv)])
        # A top-down map makes the narrow force peaks easier to inspect than 3D alone.
        f=review_figure("Fluid-force map: "*name,title*" | same data as the 3D surface";size=(1200,760))
        for coordinate in 1:2
            ax=review_axis(f,1,coordinate,"Q$(coordinate) (μN)",xlabel,"Signed velocity (mm/s)")
            surf=surfaces[coordinate].*1e6; cr=max(maximum(abs,surf),1e-12)
            hm=CairoMakie.heatmap!(ax,positions.*1e6,vv.*1e3,surf;colormap=:balance,colorrange=(-cr,cr),rasterize=2)
            CairoMakie.Colorbar(f[2,coordinate],hm;vertical=false,label="Signed force (μN)",height=14)
        end
        review_save(f,dir,"map_fluid_heatmap_"*name,"Linear signed colour scale; narrow contact peaks are retained on a locally refined position grid.",items,options)
    end
    (;limit,voltage=V,speed,bending_limit=bend)
end

function makeplots(run;options=PlotOptions())
    @assert options.full_points>=101 && options.tail_points>=101 && options.render_bins>=50
    @assert options.sweep_points>=51 && options.surface_positions>=31 && options.surface_velocities>=11
    @assert options.png_scale>0 && options.recurrence_tol>0
    CairoMakie.activate!()
    dir=joinpath(run.outdir,"plots");mkpath(dir)
    items=NamedTuple{(:name,:caption),Tuple{String,String}}[]
    t0=run.sol.t[1];tf=run.sol.t[end]
    tail=run.kind==:drive ? max(t0,tf-2/run.freq) : t0+0.8*(tf-t0)
    times=vcat(run.sol.t,collect(range(t0,tf;length=options.full_points)),
        collect(range(tail,tf;length=options.tail_points)))
    tr=review_trace(run,times)
    cycles=review_cycles(run,options,dir,items)
    events=review_events(run,tr.t,dir)
    review_window_plots(run,tr,eachindex(tr.t),"Full simulation","full",dir,items,options)
    suffix=run.kind==:drive ? "last_two_cycles" : "final_transient"
    review_window_plots(run,tr,findall(>=(tail),tr.t),cycles.label,suffix,dir,items,options)
    # One high-resolution event window complements full-run / two-cycle views.
    entries=filter(e->e[3]=="entry",events)
    if !isempty(entries)
        te=entries[end][1];a=max(t0,te-30e-6);b=min(tf,te+150e-6)
        etimes=vcat(collect(range(a,b;length=1601)),run.sol.t[(run.sol.t.>=a).&(run.sol.t.<=b)])
        etr=review_trace(run,etimes)
        review_window_plots(run,etr,eachindex(etr.t),"Last detected contact entry: detailed window","contact_detail",dir,items,options)
    end
    maps=review_force_maps(run,tr,dir,items,options)
    # Common CSV grid retains extrema of EVERY exported variable in each time bin.
    values=vcat(tr.u,reshape(tr.F,24,length(tr.t)),tr.d)
    ii=review_indices(tr.t,[view(values,k,:) for k in axes(values,1)];bins=options.render_bins)
    headers=vcat(["time_s"],["x1_m","x2_m","v1_m_s","v2_m_s","Vout_V","Wbase_J","Wbias_J","ER_J","Dfilm_J","Dstruct_J","Dwall_J","throughput_J"],
        ["Q$(i)_$(k)_N" for k in FORCE_KEYS for i in 1:2],DIAG_KEYS)
    writecsv(joinpath(dir,"diagnostics_full.csv"),headers,[(tr.t[j],values[:,j]...) for j in ii])
    ti=findall(>=(tail),tr.t)
    jt=review_indices(tr.t[ti],[view(values,k,ti) for k in axes(values,1)];bins=options.render_bins)
    writecsv(joinpath(dir,"diagnostics_"*suffix*".csv"),headers,[(tr.t[j],values[:,j]...) for j in ti[jt]])
    diag=Dict{String,Any}(
        "kind"=>string(run.kind),"frequency_Hz"=>run.freq,"time_start_s"=>t0,"time_end_s"=>tf,
        "tail_start_s"=>tail,"recurrence_label"=>cycles.label,"detected_recurrence_period"=>cycles.period,
        "recurrence_tolerance"=>options.recurrence_tol,"contact_events_detected"=>length(events),
        "diagnostic_samples_before_render_reduction"=>length(tr.t),"csv_full_samples"=>length(ii),
        "minimum_effective_gap_m"=>minimum(tr.d[17:18,:]),"maximum_Kn"=>maximum(tr.d[27,:]),
        "maximum_positive_nominal_overlap_m"=>max(0.0,maximum(tr.d[19:20,:])),
        "minimum_fluid_matrix_eigenvalue_Ns_m"=>minimum(tr.d[26,:]),
        "minimum_fluid_dissipation_W"=>minimum(tr.d[14,:]),"minimum_contact_dissipation_W"=>minimum(tr.d[16,:]),
        "maximum_force_balance_residual_N"=>maximum(abs,tr.F[:,12,:]),
        "tip_not_minimum_gap_at_any_sample"=>any(>(0),tr.d[32,:]),
        "sweep_translation_limit_m"=>maps.limit,"sweep_fixed_Vc_V"=>maps.voltage,
        "sweep_signed_speed_limit_m_s"=>maps.speed,"sweep_bending_limit_m"=>maps.bending_limit,
        "figures"=>length(items),"PNG_pixels_per_figure_unit"=>options.png_scale,
        "PDF_3D_surface_note"=>"Surface layers are rasterized; text and axes remain vector.")
    open(joinpath(dir,"review_summary.toml"),"w") do io;TOML.print(io,diag);end
    metadata=Dict("parameters"=>Dict(string(k)=>getfield(run.model.p,k) for k in fieldnames(Params)),
        "run"=>Dict("tag"=>run.tag,"kind"=>string(run.kind),"frequency_Hz"=>run.freq,"film_panels"=>run.model.panels,
            "acceleration_amplitude_m_s2"=>run.acceleration_amplitude,"julia_version"=>string(VERSION),
            "plotting_package_version"=>string(Base.pkgversion(CairoMakie))),
        "derived_model"=>Dict("nominal_contact_travel_m"=>run.model.gc,"gap_slope"=>run.model.alpha,
            "electrode_stiffness_N_m"=>run.model.ke,"M11_kg"=>run.model.M[1,1],
            "M12_kg"=>run.model.M[1,2],"M22_kg"=>run.model.M[2,2]),
        "solver_metrics"=>run.metrics,"plot_options"=>Dict(string(k)=>getfield(options,k) for k in fieldnames(PlotOptions) if !isnothing(getfield(options,k))))
    open(joinpath(run.outdir,"run_configuration.toml"),"w") do io;TOML.print(io,metadata);end
    escapehtml(s)=replace(string(s),'&'=>"&amp;",'<'=>"&lt;",'>'=>"&gt;",'"'=>"&quot;")
    open(joinpath(dir,"index.html"),"w") do io
        println(io,"<!doctype html><html><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>MEMS run review</title><style>body{font:16px/1.55 system-ui;margin:32px auto;max-width:1200px;padding:0 20px;color:#17293a;background:#f6f8fa}h1,h2{line-height:1.2}article{background:white;padding:22px;margin:24px 0;border:1px solid #d9e1e8;border-radius:8px}img{width:100%;height:auto}a{color:#006ba0}code{background:#edf1f4;padding:2px 5px}.note{padding:18px;background:#e9f2f7}nav{columns:2}li{margin:6px 0}</style></head><body>")
        println(io,"<h1>MEMS simulation review</h1><p>",escapehtml(run.tag)," · ",escapehtml(cycles.label),"</p><div class='note'>All physical states, force components, integrated work/loss states, and constitutive maps are generated from the Julia model. Final-cycle recurrence is a numerical check, not experimental validation or a stability proof. Map captions specify the held coordinates and voltage.</div>")
        println(io,"<p><a href='review_summary.toml'>Diagnostic summary</a> · <a href='diagnostics_full.csv'>Full-run data</a> · <a href='diagnostics_",suffix,".csv'>Final-window data</a> · <a href='contact_events.csv'>Contact events</a> · <a href='../run_configuration.toml'>Run parameters</a></p><nav><ol>")
        for it in items;println(io,"<li><a href='#",it.name,"'>",replace(it.name,'_'=>' '),"</a></li>");end
        println(io,"</ol></nav>")
        for it in items
            println(io,"<article id='",it.name,"'><h2>",replace(it.name,'_'=>' '),"</h2><p>",escapehtml(it.caption),"</p><a href='",it.name,".pdf'>Download PDF</a> · <a href='",it.name,".png'>Open PNG</a><img loading='lazy' src='",it.name,".png' alt='",escapehtml(it.caption),"'></article>")
        end
        println(io,"</body></html>")
    end
    println("Saved ",length(items)," figures (PDF + PNG). Open: ",joinpath(dir,"index.html"))
    (;directory=dir,index=joinpath(dir,"index.html"),summary=diag)
end

function main(args=ARGS)
    mode=isempty(args) ? "--drive" : args[1]
    opts=Dict{String,String}();j=2
    while j<=length(args)
        key=args[j]
        key in ("--cycles","--freq","--tag","--outdir") || error("Unknown option: $key")
        j<length(args) || error("Missing value for $key")
        opts[key]=args[j+1];j+=2
    end
    cycles=parse(Float64,get(opts,"--cycles","10"))
    freq=parse(Float64,get(opts,"--freq","20"))
    outdir=get(opts,"--outdir",joinpath(@__DIR__,"results"))
    if mode=="--verify";verify()
    elseif mode=="--probe";simulate(;outdir,tag=get(opts,"--tag","probe"))
    elseif mode=="--convergence";convergence(;outdir)
    elseif mode=="--drive";simulate(;kind=:drive,cycles,freq,outdir,tag=get(opts,"--tag","drive"))
    elseif mode=="--help"
        println("""
        julia --project=. collision_model_corrected.jl [MODE] [OPTIONS]
        No arguments: a ten-cycle driven simulation with automatic PDF/PNG plots.
        Modes: --drive | --probe | --convergence | --verify | --help
        Drive options: --cycles 30 --freq 20 --tag my_run --outdir results
        --verify performs algebraic checks only; it does not fabricate a trajectory.
        Every simulation writes a new timestamped folder and plots/index.html.
        See README.md for VS Code setup and PlotOptions for force-map controls.
        """)
    else;error("Unknown mode: $mode")
    end
end
end # module
if abspath(PROGRAM_FILE)==@__FILE__
    CorrectedMEMS.main()
end
