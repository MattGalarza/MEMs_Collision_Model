# Included inside module test. SI units throughout.
Base.@kwdef struct Case
    name::String = "baseline"
    contact::Symbol = :penalty       # :penalty or :constraint
    electric::Symbol = :legacy      # :legacy gap floor or :geometric air gap
    contact_gap::Float64 = 50e-9     # assumed physical asperity/contact onset
    fluid_gap::Float64 = 50e-9       # independent hydraulic effective gap
    bending_scale::Float64 = 1.0
    bias::Float64 = 3.0
    alpha::Float64 = 4.95            # peak acceleration in g, NOT m/s^2
    slope_scale::Float64 = 1.0
    ce::Float64 = 0.0               # Ns/m PER BEAM, illustrative unless measured
    c1::Float64 = 0.0
    kw_scale::Float64 = 1.0
    cw_scale::Float64 = 1.0
    eps_gap::Float64 = 2e-9
    eps_wall::Float64 = 0.5e-9
    seal_width::Float64 = 25e-9
    seal::Bool = true
end
struct Plant
    b::BM.Model
    case::Case
    gc::Float64
    kb::Float64
end
function plant(c=Case();panels=128)
    @assert c.contact in (:penalty,:constraint) && c.electric in (:legacy,:geometric)
    @assert c.contact_gap>=0 && c.fluid_gap>0 && c.bending_scale>0 && c.slope_scale>=0
    @assert c.bias>=0 && c.alpha>=0 && c.ce>=0 && c.c1>=0
    p=BM.Params(h_eff=c.fluid_gap,Vbias=c.bias,ce=c.ce,c1=c.c1,
        gap_slope=(30e-6-9e-6)/450e-6*c.slope_scale,kw=1e6*c.kw_scale,
        cw=50.0*c.cw_scale,eps_gap=c.eps_gap,eps_wall=c.eps_wall,
        seal_width=c.seal_width,seal_at_contact=c.seal)
    b=BM.Model(p;panels)
    Plant(b,c,p.g0-2p.Tp-c.contact_gap,p.n_beams*b.ke*c.bending_scale)
end

# Smooth energy-consistent continuation used ONLY if an implicit solver's trial
# stage strays outside admissible geometry. Accepted states are checked separately.
# For every physical air gap H>=0, recip and its derivative are EXACTLY 1/(H+hd).
function reciprocal_trial(z,zmin)
    if z>=zmin
        return 1/z,-1/z^2
    else
        d=z-zmin
        return 1/zmin-d/zmin^2+d^2/zmin^3,-1/zmin^2+2d/zmin^3
    end
end
function constitutive(m::Plant,x1,x2;film=true)
    b=m.b;p=b.p;C=p.cp;g=zeros(2);D=zeros(2,2);hraw=Inf
    for r in (-1.,1.)
        H=p.g0-2p.Tp .+b.alpha.*b.y.-r.*(b.B1.*x1.+b.B2.*x2)
        hraw=min(hraw,minimum(H))
        d=H.-m.case.fluid_gap
        h=m.case.fluid_gap.+BM.softpos.(d,p.eps_gap)
        dh=BM.dsoftpos.(d,p.eps_gap)
        h1=-r.*dh.*b.B1;h2=-r.*dh.*b.B2
        fac=p.n_beams*p.eps0*p.Tf
        if m.case.electric==:legacy
            C+=fac*dot(b.weights,1 ./(h.+b.hd))
            g[1]-=fac*dot(b.weights,h1./(h.+b.hd).^2)
            g[2]-=fac*dot(b.weights,h2./(h.+b.hd).^2)
        else
            for j in eachindex(H)
                val,der=reciprocal_trial(H[j]+b.hd,0.1b.hd)
                C+=fac*b.weights[j]*val
                g[1]-=fac*b.weights[j]*r*b.B1[j]*der
                g[2]-=fac*b.weights[j]*r*b.B2[j]*der
            end
        end
        if film
            H1=BM.cumulative_simpson(b,h1);H2=BM.cumulative_simpson(b,h2)
            W=b.weights./(h.^2 .* (h.+b.kp));I0=sum(W)
            mu1=dot(W,H1)/I0;mu2=dot(W,H2)/I0
            Z1=H1.-mu1;Z2=H2.-mu2
            chi=p.seal_at_contact ? BM.smootherstep((r*x2-m.gc+p.seal_width)/(2p.seal_width)) : 0.
            f=12p.eta*p.Tf*p.n_beams*p.film_scale
            D[1,1]+=f*(dot(W,Z1.^2)+chi*I0*mu1^2)
            D[1,2]+=f*(dot(W,Z1.*Z2)+chi*I0*mu1*mu2)
            D[2,2]+=f*(dot(W,Z2.^2)+chi*I0*mu2^2)
        end
    end
    D[2,1]=D[1,2]
    (;C,grad=g,D,hraw)
end
function springs(m::Plant,x1,x2)
    s=BM.springs(m.b,x1,x2);rel=x2-x1
    dk=m.kb-m.b.p.n_beams*m.b.ke
    (;U=s.U+0.5dk*rel^2,F=s.F+dk*rel.*[1.,-1.])
end
function wall(m::Plant,x2,v2)
    m.case.contact==:constraint && return (;U=0.,F=0.,loss=0.)
    p=m.b.p;U=0.;F=0.;loss=0.
    for r in (-1.,1.)
        z=r*x2-m.gc;vr=r*v2;s=BM.softpos(z,p.eps_wall)
        A=p.kw*s^p.pw*BM.dsoftpos(z,p.eps_wall);gate=max(1+p.cw*vr,0.)
        U+=p.n_beams*p.kw*s^(p.pw+1)/(p.pw+1)
        F-=p.n_beams*r*A*gate;loss+=p.n_beams*A*vr*(gate-1)
    end
    (;U,F,loss)
end
function snapshot(m,u,a,phase=0)
    b=m.b;p=b.p;x1,x2,v1,v2,Vo=u[1:5];v=[v1,v2]
    c=constitutive(m,x1,x2);s=springs(m,x1,x2);w=wall(m,x2,v2)
    es=0.5*(p.Vbias-Vo)^2*c.grad;fluid=-c.D*v;structural=-b.Cstruct*v;base=-b.beta*a
    F=s.F+es+fluid+structural+base
    lambda=phase==0 ? 0. : phase*(F[2]-b.M[2,1]/b.M[1,1]*F[1])
    fc=phase==0 ? w.F : -phase*lambda
    dd=phase==0 ? b.Minv*(F+[0.,fc]) : [F[1]/b.M[1,1],0.]
    (;c,s,w,es,fluid,structural,base,lambda,fc,dd,F)
end
# 13 states: positions, velocities, output voltage, then integrated work/losses.
const HEAD=["x1_m","x2_m","v1_m_s","v2_m_s","Vout_V","Wbase_J","Wbias_J",
    "ER_J","Dfilm_J","Dstruct_J","Dwall_J","Dimpact_J","throughput_J"]
function rhs!(du,u,ctx,t)
    m,accel,phase=ctx;p=m.b.p;s=snapshot(m,u,accel(t),phase);v=u[3:4];Vo=u[5];Vc=p.Vbias-Vo
    du[1]=u[3];du[2]=phase==0 ? u[4] : 0.;du[3:4]=s.dd
    du[5]=-Vo/(p.Rload*s.c.C)+Vc/s.c.C*dot(s.c.grad,v)
    du[6]=-accel(t)*dot(m.b.beta,v);du[7]=p.Vbias*Vo/p.Rload
    du[8]=Vo^2/p.Rload;du[9]=dot(v,s.c.D*v);du[10]=dot(v,m.b.Cstruct*v)
    du[11]=s.w.loss;du[12]=0.
    du[13]=abs(du[6])+abs(du[7])+sum(du[8:11])
    nothing
end
energy(m,u)=0.5dot(u[3:4],m.b.M*u[3:4])+springs(m,u[1],u[2]).U+
    wall(m,u[2],u[4]).U+0.5constitutive(m,u[1],u[2];film=false).C*(m.b.p.Vbias-u[5])^2
ledger(m,u,E0)=energy(m,u)-E0-u[6]-u[7]+sum(u[8:12])
acceleration(c,freq)=t->c.alpha*9.80665*0.5*(1-cos(pi*min(t*freq/4,1.)))*sin(2pi*freq*t)

mutable struct Run
    model::Plant
    segments::Vector{Any}
    events::Vector{Any}
    scale::Vector{Float64}
    accel::Any
    E0::Float64
    t0::Float64
    tf::Float64
    rtol::Float64
    dtmax::Float64
    projection_error_J::Float64
end
function at(run,t)
    j=searchsortedlast([s.a for s in run.segments],t);j=clamp(j,1,length(run.segments))
    s=run.segments[j]
    s.sol(clamp(t,s.a,s.b)).*run.scale,s.phase
end
function simulate(m::Plant;freq=20.,cycles=8.,rtol=2e-7,atol=2e-10,dtmax=2e-5,
                  u0=zeros(13),tspan=(0.,cycles/freq),accel=acceleration(m.case,freq))
    xs=m.gc;vs=xs*sqrt(m.b.k1/sum(m.b.M));Es=m.b.k1*xs^2
    scale=[xs,xs,vs,vs,max(m.case.bias,1.),fill(Es,8)...]
    u=Float64.(copy(u0));u[6:13].=0.;E0=energy(m,u)
    segments=Any[];events=Any[];phase=0;t=Float64(tspan[1]);tf=Float64(tspan[2]);projerr=0.
    abs(u[2])<=m.gc+1e-13 || m.case.contact==:penalty || error("Initial tip outside contact constraint")
    while t<tf-1e-13
        hit=Ref(0)
        function f!(dz,z,ctx,tt);rhs!(dz,z.*scale,ctx,tt);dz./=scale;nothing;end
        ctx=(m,accel,phase)
        prob=SciMLBase.ODEProblem(f!,u./scale,(t,tf),ctx)
        cb=nothing
        if m.case.contact==:constraint
            if phase==0
                cp=(z,tt,integ)->(z[2]*scale[2]-m.gc)/xs
                cm=(z,tt,integ)->(-z[2]*scale[2]-m.gc)/xs
                ap=integ->(hit[]=1;SciMLBase.terminate!(integ))
                am=integ->(hit[]=-1;SciMLBase.terminate!(integ))
                cb=SciMLBase.CallbackSet(
                    SciMLBase.ContinuousCallback(cp,ap,nothing;abstol=1e-11,reltol=0.,interp_points=20),
                    SciMLBase.ContinuousCallback(cm,am,nothing;abstol=1e-11,reltol=0.,interp_points=20))
            else
                condition=(z,tt,integ)->snapshot(m,z.*scale,accel(tt),phase).lambda/1e-4
                release=integ->(hit[]=2;SciMLBase.terminate!(integ))
                cb=SciMLBase.ContinuousCallback(condition,nothing,release;
                    abstol=1e-11,reltol=0.,interp_points=20,rootfind=SciMLBase.RightRootFind)
            end
        end
        sol=SciMLBase.solve(prob,OrdinaryDiffEqRosenbrock.Rodas5P(autodiff=ADTypes.AutoFiniteDiff());
            reltol=rtol,abstol=atol,dtmax,maxiters=10^7,save_everystep=true,dense=true,callback=cb)
        SciMLBase.successful_retcode(sol) || error("$(m.case.name): solver $(sol.retcode)")
        push!(segments,(a=t,b=sol.t[end],phase=phase,sol=sol))
        u=sol.u[end].*scale;t=sol.t[end]
        hit[]==0 && (t>=tf-1e-12 ? break : error("Integration stopped early"))
        if phase==0
            r=hit[];normal=[0.,Float64(r)];v=copy(u[3:4]);vn=dot(normal,v)
            vnew=v-m.b.Minv*normal*(vn/dot(normal,m.b.Minv*normal)) # e=0 reference
            loss=0.5dot(v,m.b.M*v)-0.5dot(vnew,m.b.M*vnew)
            before=energy(m,u);u[2]=r*m.gc
            projerr+=energy(m,u)-before
            u[3:4]=vnew;u[4]=0.;u[12]+=max(loss,0.);u[13]+=max(loss,0.)
            reaction=snapshot(m,u,accel(t),r).lambda
            phase=reaction>=0 ? r : 0
            push!(events,(time=t,side=r,event=phase==0 ? "impact_no_hold" : "capture",
                speed=vn,loss=loss,reaction=reaction))
        else
            r=phase;push!(events,(time=t,side=r,event="release",speed=0.,loss=0.,reaction=snapshot(m,u,accel(t),r).lambda))
            phase=0
        end
        length(events)>1000 && error("Too many contact events; inspect grazing/capture assumptions")
        if length(segments)>3 && t-segments[end-2].a<1e-12
            error("Contact event made no time progress; inspect grazing state")
        end
    end
    Run(m,segments,events,scale,accel,E0,tspan[1],tf,rtol,dtmax,projerr)
end

function bisectroot(f,a,b;tol=1e-13)
    fa=f(a);fb=f(b);fa==0 && return a;fb==0 && return b
    fa*fb<=0 || error("Root is not bracketed")
    for j in 1:70
        c=(a+b)/2;fc=f(c)
        (b-a<tol || fc==0) && return c
        if signbit(fc)==signbit(fa);a=c;fa=fc;else;b=c;end
    end
    (a+b)/2
end
function crossings(run;nsub=4)
    run.model.case.contact==:constraint && return run.events
    events=Any[];gc=run.model.gc
    for seg in run.segments
        sol=seg.sol
        for j in 1:length(sol.t)-1
            knots=collect(range(sol.t[j],sol.t[j+1];length=nsub+1))
            # Add velocity extrema so paired position crossings are not hidden.
            for k in 1:nsub
                a=knots[k];b=knots[k+1];va=sol(a)[4];vb=sol(b)[4]
                va*vb<0 && push!(knots,bisectroot(t->sol(t)[4],a,b))
            end
            sort!(knots);unique!(knots)
            for r in (-1,1),k in 1:length(knots)-1
                a=knots[k];b=knots[k+1];f(t)=r*sol(t)[2]*run.scale[2]-gc
                fa=f(a);fb=f(b)
                if fa*fb<0
                    t=bisectroot(f,a,b);vn=r*sol(t)[4]*run.scale[4]
                    push!(events,(time=t,side=r,event=fb>fa ? "entry" : "exit",speed=vn,loss=0.,reaction=0.))
                end
            end
        end
    end
    sort!(events;by=e->e.time)
    events
end

function stiffness(m,x;V=m.case.bias)
    F(q)=springs(m,q...).F+[0.,wall(m,q[2],0.).F]+0.5V^2*constitutive(m,q...;film=false).grad
    h=1e-11
    K=hcat([-(F(x+[j==i ? h : 0. for j in 1:2])-F(x-[j==i ? h : 0. for j in 1:2]))/(2h) for i in 1:2]...)
    Symmetric((K+K')/2)
end
function local_modes(m,x)
    K=stiffness(m,x);E=eigen(K,Symmetric(m.b.M));D=constitutive(m,x...).D+m.b.Cstruct
    dv=1e-7;D[2,2]+=(wall(m,x[2],-dv).F-wall(m,x[2],dv).F)/(2dv)
    f=[v>0 ? sqrt(v)/(2pi) : NaN for v in E.values]
    z=[E.values[j]>0 ? dot(E.vectors[:,j],D*E.vectors[:,j])/(2sqrt(E.values[j])) : NaN for j in 1:2]
    (;f,z)
end

# Independent Euler-Bernoulli FEM benchmark for tapered FREE electrode; this is
# structural mesh convergence, not a distributed post-contact transient model.
function beam_fem(p,ne)
    nd=2(ne+1);K=zeros(nd,nd);M=zeros(nd,nd);l=p.Lf/ne
    gx,gw=BM.gausslegendre(5)
    for e in 1:ne
        kk=zeros(4,4);mm=zeros(4,4)
        for k in eachindex(gx)
            z=(gx[k]+1)/2;s=(e-1+z)*l
            N=[1-3z^2+2z^3,l*(z-2z^2+z^3),3z^2-2z^3,l*(-z^2+z^3)]
            B=[(-6+12z)/l^2,(-4+6z)/l,(6-12z)/l^2,(-2+6z)/l]
            kk+=gw[k]*l/2*BM.EI(p,s)*(B*B')
            mm+=gw[k]*l/2*p.rho*p.Tf*BM.width(p,s)*(N*N')
        end
        ids=(2*e-1):(2*e+2);K[ids,ids]+=kk;M[ids,ids]+=mm
    end
    K=K[3:end,3:end];M=M[3:end,3:end];f=zeros(size(K,1));f[end-1]=1
    disp=K\f;freq=sqrt.(eigvals(Symmetric(K),Symmetric(M)))./(2pi)
    (;ke=1/disp[end-1],freq=freq[1:3])
end
