function change(c::Case;kwargs...)
    Case(;merge((;(k=>getfield(c,k) for k in fieldnames(Case))...),(;kwargs...))...)
end
function static_force(m,q,V;tilt=1e-8)
    springs(m,q...).F+0.5V^2*constitutive(m,q...;film=false).grad+m.b.beta/sum(m.b.M)*tilt
end
function static_jac(m,q,V)
    h=1e-11
    hcat([(static_force(m,q+[j==i ? h : 0. for j in 1:2],V)-static_force(m,q-[j==i ? h : 0. for j in 1:2],V))/(2h) for i in 1:2]...)
end
function free_equilibrium(m,V,seed)
    q=copy(seed)
    for it in 1:45
        f=static_force(m,q,V)
        if norm(f)<1e-11
            K=Symmetric(-static_jac(m,q,V))
            ok=abs(q[2])<m.gc && eigmin(K)>0 && constitutive(m,q...;film=false).hraw>=0
            return ok ? (q=q,phase="free",reaction=0.,residual=norm(f),kmin=eigmin(K)) : nothing
        end
        J=static_jac(m,q,V);abs(det(J))<1e-15 && return nothing
        step=J\f;eta=1.;accepted=false
        for j in 1:24
            z=q-eta*step
            if abs(z[2])<m.gc && abs(z[1])<2m.gc && constitutive(m,z...;film=false).hraw>=0 && norm(static_force(m,z,V))<norm(f)
                q=z;accepted=true;break
            end
            eta/=2
        end
        accepted || return nothing
    end
    nothing
end
function contact_equilibrium(m,V)
    f(x)=static_force(m,[x,m.gc],V)[1]
    a=m.gc-3e-6;b=m.gc+3e-6
    f(a)*f(b)<0 || return nothing
    x=bisectroot(f,a,b;tol=1e-15);q=[x,m.gc];F=static_force(m,q,V)
    J=static_jac(m,q,V)
    F[2]>=0 && -J[1,1]>0 && constitutive(m,q...;film=false).hraw>=-1e-13 || return nothing
    (q=q,phase="contact",reaction=F[2],residual=abs(F[1]),kmin=-J[1,1])
end
function dc_sweep(m;vmax=60.,dv=0.5)
    rows=[];transitions=[];q=zeros(2);oldphase="free";previousV=0.
    for (direction,voltages) in [("up",collect(0:dv:vmax)),("down",collect(vmax:-dv:0))]
        for V in voltages
            ce=contact_equilibrium(m,V)
            candidates=Any[]
            # Follow the existing stable branch; a failed Newton solve alone does
            # not establish pull-in, so retry independent free-branch seeds.
            for seed in (q,zeros(2),fill(0.3m.gc,2),fill(0.75m.gc,2))
                eq=free_equilibrium(m,V,seed);eq===nothing || push!(candidates,eq)
            end
            eq=nothing
            if oldphase=="contact" && ce!==nothing
                eq=ce
            elseif !isempty(candidates)
                eq=candidates[argmin([norm(c.q-q) for c in candidates])]
            elseif ce!==nothing
                eq=ce
            end
            eq===nothing && error("No admissible stable equilibrium found at V=$V; enlarge/refine static search")
            q=eq.q
            if eq.phase!=oldphase
                push!(transitions,(m.case.name,direction,oldphase,eq.phase,min(previousV,V),max(previousV,V),dv))
            end
            C=constitutive(m,q...;film=false).C
            Pi=springs(m,q...).U-0.5V^2*C-dot(m.b.beta/sum(m.b.M)*1e-8,q)
            push!(rows,(m.case.name,direction,V,q[1],q[2],eq.phase,eq.reaction,eq.residual,eq.kmin,C,Pi))
            oldphase=eq.phase;previousV=V
        end
    end
    (;rows,transitions)
end

function invariant_checks(panels)
    checks=[]
    add(name,value,tolerance,ok)=push!(checks,(name,value,tolerance,ok ? "PASS" : "FAIL"))
    m=plant(Case();panels);p=m.b.p
    for (j,x) in enumerate(([0.,0.],[m.gc-20e-9,m.gc-30e-9],[m.gc+30e-9,m.gc+10e-9]))
        c=constitutive(m,x...);old=BM.constitutive(m.b,x...)
        e=norm(c.grad-old.grad)/max(norm(old.grad),1e-20)
        add("baseline_gradient_$j",e,1e-9,e<1e-9)
        e=norm(c.D-old.D)/norm(old.D);add("baseline_fluid_$j",e,1e-9,e<1e-9)
        ev=eigmin(Symmetric(c.D));add("passive_fluid_$j",ev,0.,ev>=-1e-12)
    end
    m=plant(Case(contact=:constraint,electric=:geometric);panels)
    x=[m.gc-30e-9,m.gc-20e-9];h=1e-11;c=constitutive(m,x...;film=false)
    numeric=[(constitutive(m,(x+[j==i ? h : 0. for j in 1:2])...;film=false).C-
        constitutive(m,(x-[j==i ? h : 0. for j in 1:2])...;film=false).C)/(2h) for i in 1:2]
    err=norm(numeric-c.grad)/norm(c.grad);add("geometric_capacitance_gradient",err,1e-5,err<1e-5)
    v=[0.003,0.002];n=[0.,1.];vp=v-m.b.Minv*n*dot(n,v)/dot(n,m.b.Minv*n)
    loss=0.5dot(v,m.b.M*v)-0.5dot(vp,m.b.M*vp)
    add("impact_normal_velocity",abs(vp[2]),1e-14,abs(vp[2])<1e-14)
    add("impact_energy_loss",loss,0.,loss>=0)
    err=abs((m.b.M*(vp-v))[1]);add("impact_tangential_momentum",err,1e-18,err<1e-18)
    # The limiting geometric-gap capacitance must not depend on the hydraulic floor.
    n=plant(change(m.case;fluid_gap=20e-9);panels)
    a=constitutive(m,x...;film=false);b=constitutive(n,x...;film=false)
    err=abs(a.C-b.C)/a.C
    # Graded quadrature nodes change with fluid_gap; compare within integration error.
    add("electric_hydraulic_separation",err,1e-5,err<1e-5)
    checks
end

function static_tables(panels)
    tables=Dict{String,Any}();rows=[]
    for hc in [0.,5.,20.,50.]*1e-9,Vs in [1.,3.,6.],bs in [0.5,1.,2.]
        c=Case(contact=:constraint,electric=:geometric,contact_gap=hc,bias=Vs,bending_scale=bs)
        m=plant(c;panels);es=0.5Vs^2*constitutive(m,m.gc,m.gc;film=false).grad[2]
        push!(rows,(hc,Vs,bs,m.kb,es,es/m.kb))
    end
    tables["HoldingScales"]=(header=["contact_gap_m","bias_V","bending_scale","array_stiffness_N_m","tip_ES_N","approx_hold_bend_m"],rows=rows)
    rows=[]
    for scale in [0.,0.25,0.5,1.,2.]
        m=plant(Case(contact=:constraint,electric=:geometric,slope_scale=scale);panels)
        F=4.5constitutive(m,m.gc,m.gc;film=false).grad[2]
        push!(rows,(scale,m.b.alpha,F,F/m.kb))
    end
    tables["GeometrySensitivity"]=(header=["gap_slope_scale","gap_slope","tip_ES_N","approx_hold_bend_m"],rows=rows)
    rows=[]
    for ne in [8,16,32,64]
        b=BM.Model();f=beam_fem(b.p,ne)
        push!(rows,(ne,f.ke,b.ke,abs(f.ke/b.ke-1),f.freq...,
            sqrt(b.p.n_beams*b.ke/b.M[2,2])/(2pi)))
    end
    tables["BeamFEM"]=(header=["elements","FEM_ke_N_m","analytical_ke_N_m","relative_stiffness_error","f1_Hz","f2_Hz","f3_Hz","single_shape_rayleigh_Hz"],rows=rows)
    rows=[]
    for electric in [:legacy,:geometric]
        m=plant(Case(electric=electric);panels)
        for dn in range(-200.,40.;length=601)
            x=[m.gc+30e-9,m.gc+dn*1e-9];c=constitutive(m,x...);w=wall(m,x[2],0.)
            push!(rows,(string(electric),dn*1e-9,4.5c.grad...,c.C,c.hraw,-w.F,
                c.D[1,1],c.D[1,2],c.D[2,2]))
        end
    end
    tables["ForceGapCurves"]=(header=["electric_model","overlap_m","Q1_ES_N","Q2_ES_N","capacitance_F","minimum_raw_gap_m","wall_force_magnitude_N","D11_Ns_m","D12_Ns_m","D22_Ns_m"],rows=rows)
    rows=[]
    for cw in [0.,25.,50.,100.,250.,500.],v in range(0.0001,0.01;length=100)
        push!(rows,(cw,v,BM.restitution(v,cw)))
    end
    tables["RestitutionLaw"]=(header=["cw_s_m","incident_speed_m_s","isolated_contact_restitution"],rows=rows)
    rows=[]
    for bs in [0.5,1.,2.],ce in [0.,2e-5,2e-4]
        m=plant(Case(bending_scale=bs,ce=ce);panels)
        md=local_modes(m,[m.gc+35e-9,m.gc+9e-9])
        push!(rows,(bs,ce,md.f...,md.z...))
    end
    tables["ContactModes"]=(header=["bending_scale","ce_per_beam_Ns_m","low_frequency_Hz","high_frequency_Hz","low_modal_damping_ratio","high_modal_damping_ratio"],rows=rows)
    tables
end
