# Plotting, metrics, workbook export, and suite orchestration (inside module test).

mutable struct Review
    directory::String
    tables::Dict{String,Any}
    figures::Vector{Any}
    checks::Vector{Any}
    errors::Vector{Any}
end
function table!(r,name,header,rows)
    r.tables[name]=(header=String.(header),rows=collect(rows))
end
function check!(r,name,value,limit,ok)
    push!(r.checks,(name,value,limit,ok ? "PASS" : "FAIL"))
end
function run_directory(root,alpha,bias,freq)
    number(x)=@sprintf("%.8g",x)
    base=joinpath(root,"RUN_$(number(alpha))g_$(number(bias))V_$(number(freq))Hz_TESTS")
    dir=base;j=2
    while ispath(dir);dir=base*@sprintf("__%02d",j);j+=1;end
    mkpath(joinpath(dir,"plots"));dir
end
const AUXHEAD=["base_acceleration_m_s2","overlap_plus_m","overlap_minus_m","x2_minus_x1_m",
    "capacitance_F","minimum_raw_gap_m","Q1_spring_N","Q2_spring_N","Q1_electric_N","Q2_electric_N",
    "Q1_fluid_N","Q2_fluid_N","Q1_structural_damping_N","Q2_structural_damping_N",
    "Q1_base_N","Q2_base_N","Q2_contact_N","compressive_reaction_N",
    "total_stored_energy_J","energy_balance_residual_J","kinetic_energy_J","elastic_energy_J","electrical_energy_J",
    "Q1_beam_bending_N","Q1_linear_suspension_N","Q1_cubic_suspension_N","Q1_secondary_stop_N"]
function trace(run,times)
    t=sort!(unique!(Float64.(times)));u=zeros(13,length(t));phase=zeros(Int,length(t));a=zeros(length(AUXHEAD),length(t))
    m=run.model
    for j in eachindex(t)
        v,ph=at(run,t[j]);u[:,j]=v;phase[j]=ph;s=snapshot(m,v,run.accel(t[j]),ph)
        ek=0.5dot(v[3:4],m.b.M*v[3:4]);us=s.s.U+s.w.U;ee=0.5s.c.C*(m.case.bias-v[5])^2
        a[:,j]=[run.accel(t[j]),v[2]-m.gc,-v[2]-m.gc,v[2]-v[1],s.c.C,s.c.hraw,
            s.s.F...,s.es...,s.fluid...,s.structural...,s.base...,s.fc,s.lambda,
            ek+us+ee,ek+us+ee-run.E0-v[6]-v[7]+sum(v[8:12]),ek,us,ee,
            m.kb*(v[2]-v[1]),-m.b.k1*v[1],-m.b.k3*v[1]^3,
            s.s.F[1]-m.kb*(v[2]-v[1])+m.b.k1*v[1]+m.b.k3*v[1]^3]
    end
    (;t,u,phase,aux=a)
end
function trace_table!(r,name,tr)
    table!(r,name,["time_s","contact_side",HEAD...,AUXHEAD...],
        ([tr.t[j],tr.phase[j],tr.u[:,j]...,tr.aux[:,j]...] for j in eachindex(tr.t)))
end
function last_sequence(events,side)
    groups=Vector{Any}[]
    for e in events
        if isempty(groups) || groups[end][1].side!=e.side;push!(groups,Any[]);end
        push!(groups[end],e)
    end
    g=filter(x->x[1].side==side,groups)
    isempty(g) ? Any[] : g[end]
end
function event_table!(r,name,events)
    table!(r,name,["time_s","side","event","normal_incident_speed_m_s","impact_loss_J","reaction_at_event_N"],
        ((e.time,e.side,e.event,e.speed,e.loss,e.reaction) for e in events))
end
function detailed_event_table!(r,name,run,events)
    rows=Any[]
    for e in events
        u,ph=at(run,e.time);s=snapshot(run.model,u,run.accel(e.time),ph)
        push!(rows,(e.time,e.side,e.event,e.speed,e.loss,e.reaction,u[1:5]...,
            e.side*(u[2]-u[1]),e.side*s.es[2],e.side*s.s.F[2],e.side*s.fluid[2],
            e.side*s.fc,s.c.C,s.c.hraw))
    end
    table!(r,name,["time_s","side","event","normal_incident_speed_m_s","impact_loss_J","reaction_at_event_N",
        HEAD[1:5]...,"opening_bend_m","outward_tip_ES_N","outward_tip_spring_N","outward_tip_fluid_N","outward_tip_contact_N","capacitance_F","minimum_raw_gap_m"],rows)
end
function chatter_table!(r,name,run,events)
    rows=Any[]
    for side in (-1,1)
        seq=last_sequence(events,side)
        for j in eachindex(seq)
            e=seq[j]
            if e.event in ("exit","release") && j<length(seq) && seq[j+1].event in ("entry","capture","impact_no_hold")
                b=seq[j+1].time
                values=[run.model.gc-side*at(run,t)[1][2] for t in range(e.time,b;length=201)]
                gap=max(maximum(values),0.)
                push!(rows,(side,e.time,b,b-e.time,gap,gap>0.1e-9,gap>1e-9,gap>10e-9))
            end
        end
    end
    table!(r,name,["side","exit_or_release_s","reentry_or_capture_s","separated_duration_s","maximum_sampled_separation_m",
        "resolved_above_0p1_nm","resolved_above_1_nm","resolved_above_10_nm"],rows)
end
function episode_times(run,events;side=1)
    seq=last_sequence(events,side);isempty(seq) && return Float64[]
    a=max(run.t0,seq[1].time-100e-6);b=min(run.tf,seq[end].time+200e-6)
    # A dense approach window and a dense release window supplement the whole
    # encounter trace. These are plot samples, not the event-counting algorithm.
    sort!(unique!(vcat(collect(range(a,b;length=2001)),
        collect(range(a,min(b,seq[1].time+400e-6);length=2501)),
        collect(range(max(a,seq[end].time-200e-6),b;length=2001)),[e.time for e in seq])))
end
function metrics(run,events;freq=20.,scope="drive")
    m=run.model;maxres=0.;mingap=Inf;maxover=-Inf;minreaction=Inf;maxconstraint=0.;n=0;held=0.
    maxspring=0.;maxes=0.;maxv1=0.;maxv2=0.
    for seg in run.segments
        seg.phase!=0 && (held+=seg.b-seg.a)
        for z in seg.sol.u
            u=z.*run.scale;n+=1
            H=m.b.p.g0-2m.b.p.Tp .+m.b.alpha.*m.b.y.-abs.(m.b.B1.*u[1].+m.b.B2.*u[2])
            mingap=min(mingap,minimum(H));maxover=max(maxover,abs(u[2])-m.gc)
            maxres=max(maxres,abs(ledger(m,u,run.E0)))
            maxv1=max(maxv1,abs(u[3]));maxv2=max(maxv2,abs(u[4]))
            maxspring=max(maxspring,abs(springs(m,u[1],u[2]).F[1]))
            if seg.phase!=0
                maxconstraint=max(maxconstraint,abs(u[2]-seg.phase*m.gc),abs(u[4])*1e-6)
            end
        end
        if seg.phase!=0
            for tt in range(seg.a,seg.b;length=101)
                uu=seg.sol(tt).*run.scale
                minreaction=min(minreaction,snapshot(m,uu,run.accel(tt),seg.phase).lambda)
            end
        end
    end
    uf=at(run,run.tf)[1];den=max(uf[13],1e-30)
    recurrence=NaN
    if scope=="drive" && run.tf-run.t0>=2/freq
        ts=range(run.tf-1/freq,run.tf;length=2001)
        aa=hcat([at(run,t)[1][1:5] for t in ts]...)
        bb=hcat([at(run,t-1/freq)[1][1:5] for t in ts]...)
        scales=run.scale[1:5]
        recurrence=maximum(abs.((aa-bb)./scales))
    end
    ts=range(max(run.t0,run.tf-2/freq),run.tf;length=4001)
    voltage=[at(run,t)[1][5] for t in ts]
    groups=[last_sequence(events,s) for s in (-1,1)]
    maxcross=isempty(events) ? 0 : maximum(length.(groups))
    entering=count(e->e.event in ("entry","capture","impact_no_hold"),events)
    leaving=count(e->e.event in ("exit","release"),events)
    (;case=m.case.name,scope,contact=string(m.case.contact),electric=string(m.case.electric),
      alpha_g=m.case.alpha,bias_V=m.case.bias,accepted_states=n,events=length(events),
      entry_or_impact_count=entering,exit_or_release_count=leaving,last_sequence_max_events=maxcross,
      max_overlap_m=maxover,min_raw_gap_m=mingap,held_time_s=held,
      min_sampled_held_reaction_N=isfinite(minreaction) ? minreaction : NaN,
      held_constraint_error_m=maxconstraint,energy_residual_J=maxres,energy_residual_fraction=maxres/den,
      projection_energy_error_J=run.projection_error_J,impact_loss_J=uf[12],
      cycle_change_scaled=recurrence,Vrms_V=sqrt(mean(abs2,voltage)),Vpeak_abs_V=maximum(abs.(voltage)),
      mean_load_power_W=mean(abs2,voltage)/m.b.p.Rload,max_v1_m_s=maxv1,max_v2_m_s=maxv2,
      max_spring_on_shuttle_N=maxspring)
end
function numerical_checks!(r,run,met)
    tag=run.model.case.name
    check!(r,"$tag: normalized energy balance",met.energy_residual_fraction,1e-4,met.energy_residual_fraction<1e-4)
    check!(r,"$tag: physical air gap >= 0",met.min_raw_gap_m,-1e-12,met.min_raw_gap_m>=-1e-12)
    if run.model.case.contact==:constraint
        check!(r,"$tag: rigid tip constraint",met.held_constraint_error_m,1e-11,met.held_constraint_error_m<1e-11)
        check!(r,"$tag: admissible sampled reaction",met.min_sampled_held_reaction_N,-1e-9,
            isnan(met.min_sampled_held_reaction_N) || met.min_sampled_held_reaction_N>=-1e-9)
    end
end
excelvalue(v)=v isa AbstractFloat && !isfinite(v) ? string(v) : v isa Symbol ? string(v) : v
function save_workbook(r,path)
    # Public API supported by pinned XLSX 0.10.4; formatting is intentionally
    # optional and cannot abort the export (no freezePanes/newxlsx assumptions).
    table!(r,"Checks",["check","measured_value","comparison_threshold","status"],r.checks)
    table!(r,"Errors",["stage","error"],r.errors)
    keys_sorted=sort!(collect(keys(r.tables)))
    index=[(name,length(r.tables[name].rows),join(r.tables[name].header,"; ")) for name in keys_sorted]
    table!(r,"SheetIndex",["sheet","rows","columns_with_units"],index)
    XLSX.openxlsx(path,mode="w") do xf
        firstsheet=true
        for name in vcat(["SheetIndex"],keys_sorted)
            tab=r.tables[name];N=length(tab.rows)
            for (part,lo) in enumerate(1:1000000:max(N,1))
                nm=first(name,min(length(name),part==1 ? 31 : 27))
                nm=part==1 ? nm : nm*"_p$part"
                if firstsheet
                    sheet=xf[1];XLSX.rename!(sheet,nm);firstsheet=false
                else
                    sheet=XLSX.addsheet!(xf,nm)
                end
                hi=min(N,lo+999999)
                cols=[Any[excelvalue(tab.rows[i][j]) for i in lo:hi] for j in eachindex(tab.header)]
                XLSX.writetable!(sheet,cols,tab.header)
            end
        end
    end
end
function figsave!(r,fig,name,title,caption="")
    p=joinpath(r.directory,"plots",name)
    # Verify native exports before publishing the pair. A renderer may return
    # without throwing despite leaving an empty image under memory pressure.
    for ext in (".pdf",".png")
        tmp=tempname()*ext;valid=false
        signature=ext==".png" ? UInt8[0x89,0x50,0x4e,0x47,0x0d,0x0a,0x1a,0x0a] : collect(codeunits("%PDF-"))
        for attempt in 1:2
            if ext==".png";CairoMakie.save(tmp,fig;px_per_unit=1.4);else;CairoMakie.save(tmp,fig);end
            valid=isfile(tmp) && filesize(tmp)>128 && open(io->read(io,length(signature))==signature,tmp)
            valid && break
            GC.gc()
        end
        valid || error("Renderer produced an invalid $ext export for $name after two attempts")
        mv(tmp,p*ext;force=true)
    end
    push!(r.figures,(name=name,title=title,caption=caption))
end
escapehtml(s)=replace(string(s),'&'=>"&amp;",'<'=>"&lt;",'>'=>"&gt;",'"'=>"&quot;")
function gallery!(r)
    open(joinpath(r.directory,"plots","index.html"),"w") do io
        print(io,"<!doctype html><meta charset='utf-8'><title>MEMS contact diagnostics</title><style>body{font:16px system-ui;max-width:1200px;margin:30px auto;padding:0 20px;background:#f6f7fa;color:#182334}article{background:white;padding:22px;margin:26px 0;border-radius:10px}img{width:100%;height:auto}p{line-height:1.5}a{color:#126080}</style><h1>MEMS contact diagnostics</h1><p>Numerical verification and uncalibrated model comparisons. A constrained, perfectly inelastic impact is a reference assumption, not experimental evidence of sticking. All chart data, metrics, parameters and errors are in <a href='../results.xlsx'>results.xlsx</a>.</p>")
        for f in r.figures
            print(io,"<article><h2>",escapehtml(f.title),"</h2><p>",escapehtml(f.caption),"</p><a href='",f.name,".pdf'>PDF</a><img loading='lazy' src='",f.name,".png' alt='",escapehtml(f.title),"'></article>")
        end
    end
end

const COLORS=[:steelblue,:darkorange,:seagreen,:firebrick,:mediumpurple,:sienna,:deeppink,:gray40]
function drive_figure!(r,run,tr,suffix,title)
    fig=Figure(size=(1250,960));Label(fig[0,:],title;fontsize=23)
    tt=tr.t.*1000
    ax=Axis(fig[1,1];xlabel="Time (ms)",ylabel="Displacement (um)")
    lines!(ax,tt,tr.u[1,:].*1e6;label="Shuttle x1",color=COLORS[1])
    lines!(ax,tt,tr.u[2,:].*1e6;label="Tip x2",color=COLORS[2])
    hlines!(ax,[-run.model.gc,run.model.gc].*1e6;color=:gray60,linestyle=:dash);axislegend(ax;position=:lb)
    ax=Axis(fig[1,2];xlabel="Time (ms)",ylabel="Velocity (mm/s)")
    for j in 3:4;lines!(ax,tt,tr.u[j,:].*1e3;label=j==3 ? "Shuttle" : "Tip",color=COLORS[j-2]);end
    axislegend(ax;position=:lb)
    ax=Axis(fig[2,1];xlabel="Time (ms)",ylabel="Output voltage (V)");lines!(ax,tt,tr.u[5,:];color=COLORS[3])
    ax=Axis(fig[2,2];xlabel="Time (ms)",ylabel="Energy balance residual (pJ)")
    lines!(ax,tt,tr.aux[20,:].*1e12;color=COLORS[4])
    for (col,indices,labels) in [(1,[7,9,11,13,15],["Spring","Electric","Fluid","Structural damping","Base"]),
                                 (2,[8,10,12,14,16,17],["Spring","Electric","Fluid","Structural damping","Base","Contact"])]
        ax=Axis(fig[3,col];xlabel="Time (ms)",ylabel=col==1 ? "Forces on shuttle (uN)" : "Forces on tip (uN)")
        for (j,i) in enumerate(indices);lines!(ax,tt,tr.aux[i,:].*1e6;label=labels[j],color=COLORS[j]);end
        axislegend(ax;labelsize=11,nbanks=2,position=:lb)
    end
    figsave!(r,fig,run.model.case.name*"_"*suffix,title,
        "Overview samples; use the contact-detail plots for fast ringing. Late cycles are shown without assuming steady state; cycle-to-cycle change is in Metrics.")
end
function ledger_figure!(r,run,tr;suffix="energy_states")
    fig=Figure(size=(1200,730));Label(fig[0,:],run.model.case.name*" — all accumulated work and loss states";fontsize=22)
    ax=Axis(fig[1,1];xlabel="Time (ms)",ylabel="Work / loss (nJ)")
    for j in 6:12;lines!(ax,tr.t.*1000,tr.u[j,:].*1e9;label=HEAD[j],color=COLORS[j-5]);end
    axislegend(ax;position=:lt,labelsize=12)
    ax=Axis(fig[1,2];xlabel="Time (ms)",ylabel="Energy (nJ)")
    lines!(ax,tr.t.*1000,tr.u[13,:].*1e9;label="Integrated throughput",color=COLORS[1])
    for (i,l,c) in [(21,"Kinetic",2),(22,"Elastic",3),(23,"Electrical",4)]
        lines!(ax,tr.t.*1000,tr.aux[i,:].*1e9;label=l,color=COLORS[c])
    end
    axislegend(ax;position=:lt)
    figsave!(r,fig,run.model.case.name*"_"*suffix,"Work, losses, and stored energies — "*run.model.case.name,
        "Includes every integrated ledger state. Perfectly inelastic impacts add an explicit positive Dimpact jump. Position-projection roundoff is reported separately.")
end
function contact_figure!(r,run,tr,events;side=1)
    isempty(tr.t) && return
    seq=last_sequence(events,side);tref=seq[1].time;tt=(tr.t.-tref).*1e6;m=run.model
    d=(side.*tr.u[2,:].-m.gc).*1e9;delta=(side.*(tr.u[2,:].-tr.u[1,:])).*1e9
    fig=Figure(size=(1250,1050));Label(fig[0,:],run.model.case.name*" — full contact sequence, side $(side)";fontsize=23)
    ax=Axis(fig[1,1];xlabel="Time from first entry / capture (us)",ylabel="Signed displacement from boundary (nm)")
    lines!(ax,tt,(side.*tr.u[1,:].-m.gc).*1e9;label="Shuttle x1",color=COLORS[1])
    lines!(ax,tt,d;label="Tip x2",color=COLORS[2]);hlines!(ax,[0.];color=:gray50,linestyle=:dash);axislegend(ax;position=:lb)
    ax=Axis(fig[1,2];xlabel="Time from first entry / capture (us)",ylabel="Opening bend: side*(x2-x1) (nm)")
    lines!(ax,tt,delta;color=COLORS[3]);hlines!(ax,[0.];color=:gray60,linestyle=:dash)
    ax=Axis(fig[2,1];xlabel="Time from first entry / capture (us)",ylabel="Tip overlap (nm)",title="Approach and initial rebounds")
    k=findall(x->-100<=x<=400,tt);lines!(ax,tt[k],d[k];color=COLORS[2]);hlines!(ax,[0.];color=:gray50,linestyle=:dash)
    for e in seq
        te=(e.time-tref)*1e6
        -100<=te<=400 && vlines!(ax,[te];color=(:gray50,0.4),linestyle=:dot)
    end
    ax=Axis(fig[2,2];xlabel="Tip overlap (nm)",ylabel="Outward tip speed (mm/s)",title="Boundary phase portrait")
    lines!(ax,d[k],side.*tr.u[4,k].*1e3;color=COLORS[1]);vlines!(ax,[0.];color=:gray50,linestyle=:dash)
    ax=Axis(fig[3,1];xlabel="Time from first entry / capture (us)",ylabel="Signed tip force (uN)")
    for (i,l,c) in [(8,"Spring",1),(10,"Electric",2),(12,"Fluid",3),(17,"Contact",4)]
        lines!(ax,tt,side.*tr.aux[i,:].*1e6;label=l,color=COLORS[c])
    end
    axislegend(ax;position=:lb,nbanks=2)
    ax=Axis(fig[3,2];xlabel="Time relative to last exit / release (us)",ylabel="Displacement from boundary (nm)",title="Final departure")
    exit_t=(seq[end].time-tref)*1e6;k=findall(x->abs(x-exit_t)<=200,tt)
    lines!(ax,tt[k].-exit_t,d[k];label="Tip",color=COLORS[2])
    lines!(ax,tt[k].-exit_t,(side.*tr.u[1,k].-m.gc).*1e9;label="Shuttle",color=COLORS[1])
    hlines!(ax,[0.];color=:gray50,linestyle=:dash);vlines!(ax,[0.];color=:gray60,linestyle=:dot);axislegend(ax;position=:lb)
    figsave!(r,fig,run.model.case.name*"_contact_"*(side>0 ? "plus" : "minus"),
        "Collision, bending, and departure — "*run.model.case.name*" ($(side))",
        "Zero overlap is the assumed physical contact onset; positive overlap is penalty penetration. The rigid reference remains at zero during capture. Event counts come from dense solver interpolation or constraint callbacks, not these plot samples.")
end
function static_figures!(r)
    data=r.tables["ForceGapCurves"].rows
    fig=Figure(size=(1200,760));Label(fig[0,:],"Electrical cutoff and contact-force audit";fontsize=24)
    ax=Axis(fig[1,1];xlabel="Tip overlap (nm)",ylabel="Tip electrostatic force (uN)",title="Same x1 and V = 3 V")
    for (j,model) in enumerate(["legacy","geometric"])
        rows=filter(x->x[1]==model,data);lines!(ax,[x[2]*1e9 for x in rows],[x[4]*1e6 for x in rows];label=model,color=COLORS[j])
    end
    vlines!(ax,[0.];color=:gray60,linestyle=:dash);axislegend(ax;position=:lt)
    ax=Axis(fig[1,2];xlabel="Tip overlap (nm)",ylabel="Contact force (uN)")
    rows=filter(x->x[1]=="legacy",data);lines!(ax,[x[2]*1e9 for x in rows],[x[7]*1e6 for x in rows]);vlines!(ax,[0.];color=:gray60,linestyle=:dash)
    ax=Axis(fig[2,1];xlabel="Relative stiffness scale",ylabel="Approximate holding bend (nm)",title="Geometric field; 50 nm physical contact gap")
    data=r.tables["HoldingScales"].rows
    for (j,V) in enumerate([1.,3.,6.])
        rows=filter(x->isapprox(x[1],50e-9;rtol=1e-12) && x[2]==V,data)
        lines!(ax,[x[3] for x in rows],[x[6]*1e9 for x in rows];label="$(V) V",color=COLORS[j]);scatter!(ax,[x[3] for x in rows],[x[6]*1e9 for x in rows];color=COLORS[j])
    end
    axislegend(ax;position=:rt)
    ax=Axis(fig[2,2];xlabel="Gap slope / original slope",ylabel="Approximate holding bend (nm)",title="Geometry sensitivity; 3 V, same area and ke")
    rows=r.tables["GeometrySensitivity"].rows
    lines!(ax,[x[1] for x in rows],[x[4]*1e9 for x in rows];color=COLORS[4]);scatter!(ax,[x[1] for x in rows],[x[4]*1e9 for x in rows];color=COLORS[4])
    figsave!(r,fig,"electrical_and_holding","Electrical cutoff, contact stiffness, and holding scales",
        "FES/(N ke) is a force scale, not a coupled pull-off prediction. Slope changes are counterfactual gap geometry with unchanged beam stiffness. The zero-slope case is a static illustration only; distributed contact is not simulated.")
    fig=Figure(size=(1200,760));Label(fig[0,:],"Independent beam and dissipation benchmarks";fontsize=24)
    ax=Axis(fig[1,1];xlabel="Beam finite elements",ylabel="Relative static-stiffness error",yscale=log10)
    rows=r.tables["BeamFEM"].rows;scatterlines!(ax,[x[1] for x in rows],[max(x[4],1e-14) for x in rows])
    ax=Axis(fig[1,2];xlabel="Beam finite elements",ylabel="Frequency (kHz)")
    for (j,col) in enumerate(5:7);scatterlines!(ax,[x[1] for x in rows],[x[col]/1000 for x in rows];label="FEM mode $j",color=COLORS[j]);end
    hlines!(ax,[rows[end][8]/1000];label="Single-shape Rayleigh",color=:gray40,linestyle=:dash);axislegend(ax;position=:rt)
    ax=Axis(fig[2,1];xlabel="Incident speed (mm/s)",ylabel="Isolated contact restitution")
    rows=r.tables["RestitutionLaw"].rows
    for (j,cw) in enumerate([0.,50.,250.,500.])
        rr=filter(x->x[1]==cw,rows);lines!(ax,[x[2]*1000 for x in rr],[x[3] for x in rr];label="cw = $(Int(cw)) s/m",color=COLORS[j])
    end
    axislegend(ax;position=:rb)
    ax=Axis(fig[2,2];xlabel="Relative stiffness scale",ylabel="Low-mode damping ratio",title="Frozen tangent at a representative contact state")
    rows=r.tables["ContactModes"].rows
    for (j,ce) in enumerate([0.,2e-5,2e-4])
        rr=filter(x->x[2]==ce,rows);scatterlines!(ax,[x[1] for x in rr],[x[5] for x in rr];label="ce = $(ce) Ns/m/beam",color=COLORS[j])
    end
    axislegend(ax;position=:lt,labelsize=11)
    figsave!(r,fig,"beam_and_damping","Beam stiffness, free modes, and damping",
        "The FEM benchmark is an independent free-electrode structural calculation. It is not a multi-mode contact simulation. Frozen tangent damping and isolated contact restitution are local diagnostics, not measured device parameters.")
end
function dc_figure!(r)
    rows=r.tables["DCBranches"].rows;isempty(rows) && return
    fig=Figure(size=(1200,720));Label(fig[0,:],"Quasistatic DC branch following — no sinusoidal drive";fontsize=24)
    ax=Axis(fig[1,1];xlabel="Bias voltage (V)",ylabel="Tip displacement (um)")
    ax2=Axis(fig[1,2];xlabel="Bias voltage (V)",ylabel="Compressive contact reaction (uN)")
    for (j,case) in enumerate(unique([x[1] for x in rows])),direction in ["up","down"]
        rr=filter(x->x[1]==case && x[2]==direction,rows);style=direction=="up" ? :solid : :dash
        lines!(ax,[x[3] for x in rr],[x[5]*1e6 for x in rr];label="$case $direction",color=COLORS[j],linestyle=style)
        lines!(ax2,[x[3] for x in rr],[x[7]*1e6 for x in rr];color=COLORS[j],linestyle=style)
    end
    axislegend(ax;position=:lt,labelsize=12)
    figsave!(r,fig,"dc_branch_sweep","DC pull-in / release search",
        "Stable-branch search on a 0.5 V grid with a constant 10 nN symmetry-breaking load. Transition brackets and equilibrium residuals are exported. These are approximate static branch changes, not driven pull-in voltages or a global bifurcation proof.")
end
function comparison_figure!(r,runs,evmap,freq)
    length(runs)>=2 || return
    fig=Figure(size=(1250,840));Label(fig[0,:],"Contact model comparison — identical forcing within each row";fontsize=23)
    for (row,pair) in enumerate([runs[1:2],length(runs)>=4 ? runs[3:4] : Run[]])
        isempty(pair) && continue
        ax=Axis(fig[row,1];xlabel="Time in final drive cycle (ms)",ylabel="Tip displacement (um)",title="$(pair[1].model.case.alpha) g peak")
        ax2=Axis(fig[row,2];xlabel="Time in final drive cycle (ms)",ylabel="Output voltage (V)")
        for (j,run) in enumerate(pair)
            t=collect(range(run.tf-1/freq,run.tf;length=5001));tr=trace(run,t)
            label=run.model.case.contact==:penalty ? "Penalty + original electrical floor" : "Inelastic constraint + geometric field"
            lines!(ax,(t.-t[1]).*1000,tr.u[2,:].*1e6;label,color=COLORS[j])
            lines!(ax2,(t.-t[1]).*1000,tr.u[5,:];color=COLORS[j])
            table!(r,"Compare_"*run.model.case.name,["time_s","x2_m","Vout_V"],
                ((t[k],tr.u[2,k],tr.u[5,k]) for k in eachindex(t)))
        end
        axislegend(ax;position=:lb,labelsize=11)
    end
    figsave!(r,fig,"model_comparison","Baseline versus constrained reference",
        "This headline comparison changes both contact and electrostatics. The factorial replay tests separately isolate those changes. The lower acceleration comparison is 1 g peak, not a calibrated experiment or an assumed RMS conversion.")
end

function replay_figures!(r,replays,evmap,mets)
    isempty(replays) && return
    cases=[m.case for m in mets if m.scope=="replay"]
    fig=Figure(size=(1250,1000));Label(fig[0,:],"One-encounter sensitivity tests — same initial physical state";fontsize=23)
    ax=Axis(fig[1,1];xlabel="Entry / capture events in replay window",yticks=(1:length(cases),cases))
    rows=filter(m->m.scope=="replay",mets)
    barplot!(ax,1:length(rows),[x.entry_or_impact_count for x in rows];direction=:x,color=:steelblue)
    ax=Axis(fig[1,2];xlabel="Maximum nominal penetration (nm)",yticks=(1:length(cases),cases))
    barplot!(ax,1:length(rows),[max(x.max_overlap_m,0.)*1e9 for x in rows];direction=:x,color=:darkorange)
    ax=Axis(fig[2,1];xlabel="Time from common start (us)",ylabel="Tip overlap (nm)",title="Factorial contact / electrical comparison")
    ax2=Axis(fig[2,2];xlabel="Time from common start (us)",ylabel="Tip outward speed (mm/s)",title="Initial approach and rebound")
    for (j,run) in enumerate(replays[1:min(4,length(replays))])
        ts=collect(range(run.t0,min(run.tf,run.t0+600e-6);length=6001));tr=trace(run,ts)
        lines!(ax,(ts.-ts[1]).*1e6,(tr.u[2,:].-run.model.gc).*1e9;label=run.model.case.name,color=COLORS[j])
        lines!(ax2,(ts.-ts[1]).*1e6,tr.u[4,:].*1000;color=COLORS[j])
        trace_table!(r,"Zoom_"*string(j),tr)
    end
    hlines!(ax,[0.];color=:gray50,linestyle=:dash);axislegend(ax;position=:lb,labelsize=11)
    figsave!(r,fig,"encounter_sensitivities","What changes the chatter?",
        "Counts combine entry for penalty cases and capture/impact for constrained cases, with separate event types in Excel. These are finite-window replays, not independently settled cycles. A missing final release is marked as censored in ReplayWindows.")
end

function convergence_figure!(r,conv)
    labels=[replace(replace(replace(x[1],"penalty"=>"Penalty","constraint"=>"Constraint"),"_"=>"\n"),"space"=>"mesh") for x in conv]
    fig=Figure(size=(1200,540));Label(fig[0,:],"Numerical refinement of the same contact encounter";fontsize=23)
    ax=Axis(fig[1,1];ylabel="Max position difference (nm)",xticks=(1:length(conv),labels))
    barplot!(ax,1:length(conv),[x[5]*1e9 for x in conv];color=:steelblue)
    ax=Axis(fig[1,2];ylabel="Max voltage difference (mV)",xticks=(1:length(conv),labels))
    barplot!(ax,1:length(conv),[x[6]*1000 for x in conv];color=:darkorange)
    figsave!(r,fig,"numerical_refinement","Time and spatial refinement",
        "Common initial state and common comparison times. Counts, event-time differences, absolute waveform differences and acceptance thresholds are in the workbook. This checks a contact episode, not every possible orbit or parameter combination.")
end

function run_suite(;alpha=4.95,bias=3.,freq=20.,cycles=8.,panels=128,
                   rtol=2e-7,atol=2e-10,dtmax=2e-5,
                   output_root=joinpath(@__DIR__,"results"),
                   save_accepted=false,strict=false,checkpoint_dir=nothing)
    @assert alpha>0 && bias>=0 && freq>0 && cycles>=6 && panels>=32
    @assert rtol>0 && atol>0 && dtmax>0
    dir=run_directory(output_root,alpha,bias,freq)
    r=Review(dir,Dict{String,Any}(),Any[],Any[],Any[])
    CairoMakie.set_theme!(fontsize=15,linewidth=1.7)
    println("MEMS contact diagnostic suite: $dir");flush(stdout)
    metadata=[("suite_version","1.0.0"),("created_UTC",string(Dates.now(Dates.UTC))),
        ("Julia_version",string(VERSION)),("alpha_peak_g",string(alpha)),("bias_V",string(bias)),
        ("forcing_frequency_Hz",string(freq)),("drive_cycles",string(cycles)),("quadrature_panels",string(panels)),
        ("relative_tolerance",string(rtol)),("scaled_absolute_tolerance",string(atol)),("maximum_step_s",string(dtmax)),
        ("solver","Rodas5P; finite-difference Jacobian; dense output"),
        ("physical_validation","NOT CALIBRATED: no raw experimental data supplied"),
        ("contact_reference","Unilateral rigid tip contact; perfectly inelastic mass-consistent capture; release at zero compressive reaction"),
        ("electrical_reference","Geometric air gap plus finite dielectric thickness; hydraulic cutoff independent"),
        ("capture_caveat","e=0 suppresses rebound by assumption; must be compared to measured impact/release, not selected merely for a smooth trace"),
        ("bias_connection","Ideal DC source with production Rload and common variable-capacitance circuit"),
        ("beam_limitation","Original two generalized coordinates and coherent N-beam motion retained; no distributed contact, contact-length evolution, or array disorder"),
        ("steady_state","Last two cycles are late-time windows; inspect cycle_change_scaled, do not assume steady state"),
        ("count_method","Penalty roots from dense accepted intervals with velocity-extremum subdivision; constraint events from root callbacks"),
        ("replay_scope","Same baseline pre-contact physical state and absolute drive time; integrated ledgers reset; not steady-state sensitivity"),
        ("static_bias_grid","0 to 60 to 0 V in 0.5 V increments; positive constant 10 nN symmetry-breaking load"),
        ("static_force_maps","Illustrative 3 V reference maps; HoldingScales includes 1, 3, 6 V"),
        ("damping_scales","ce per beam; kw per beam; cw inverse velocity; illustrative, not fitted"),
        ("numerical_thresholds","Energy residual/throughput < 1e-4; gap >= -1 pm; constraint error < 10 pm; sampled reaction >= -1 nN"),
        ("all_units","SI in workbook headers; plots use displayed engineering units"),
        ("accepted_states_exported",string(save_accepted))]
    for file in ["test.jl","test_physics.jl","test_static.jl","test_reports.jl","baseline_reference.jl","Project.toml","Manifest.toml"]
        isfile(joinpath(@__DIR__,file)) && push!(metadata,("sha256:"*file,bytes2hex(SHA.sha256(read(joinpath(@__DIR__,file))))))
    end
    table!(r,"Metadata",["item","value"],metadata)
    table!(r,"Dependencies",["name","version"],[(v.name,string(v.version)) for v in values(Pkg.dependencies()) if v.version!==nothing])
    table!(r,"StateDefinitions",["state","meaning"],zip(HEAD,["Shuttle displacement","Collective electrode-tip displacement","Shuttle velocity","Collective tip velocity","Load/output voltage",
        "Accumulated mechanical work from base","Accumulated electrical work from bias source","Resistor loss","Fluid loss","Structural damping loss","Compliant contact loss","Discrete inelastic impact loss","Absolute input work plus dissipated energy"]))
    p=BM.Params();table!(r,"BaselineParameters",["parameter","value"],[(string(k),getfield(p,k)) for k in fieldnames(BM.Params)])
    paramrows=Any[];mets=Any[];fullruns=Run[];replays=Run[];evmap=Dict{String,Any}();replaywindows=Any[]
    function guarded(stage,fn)
        println("  $stage");flush(stdout)
        try
            fn()
        catch err
            msg=sprint(showerror,err,catch_backtrace());push!(r.errors,(stage,msg))
            println(stderr,"  FAILED $stage: ",sprint(showerror,err));flush(stderr)
            nothing
        end
    end
    function compute(c;u0=zeros(13),tspan=(0.,cycles/freq),np=panels,tol=rtol,step=dtmax)
        for k in fieldnames(Case);push!(paramrows,(c.name,string(k),getfield(c,k)));end
        # Optional development checkpoint lives outside result folders. Default
        # runs write all numeric deliverables only into results.xlsx.
        key=bytes2hex(SHA.sha256(string(c,u0,tspan,np,tol,step,freq,atol)))
        path=checkpoint_dir===nothing ? nothing : joinpath(checkpoint_dir,key*".jls")
        if path!==nothing && isfile(path)
            return Serialization.deserialize(path)
        end
        run=simulate(plant(c;panels=np);freq,cycles,rtol=tol,atol,dtmax=step,u0,tspan)
        if path!==nothing;mkpath(dirname(path));Serialization.serialize(path,run);end
        run
    end
    function record(run,scope;prefix)
        m=run.model
        for (name,value) in [("actual_panels",m.b.panels),("array_bending_stiffness_N_m",m.kb),
            ("single_beam_nominal_stiffness_N_m",m.b.ke),("M11_kg",m.b.M[1,1]),("M12_kg",m.b.M[1,2]),
            ("M22_kg",m.b.M[2,2]),("actual_gap_slope",m.b.alpha),("contact_boundary_m",m.gc),
            ("dielectric_equivalent_thickness_m",m.b.hd),("actual_rtol",run.rtol),("actual_dtmax_s",run.dtmax)]
            push!(paramrows,(m.case.name,name,value))
        end
        events=crossings(run);evmap[run.model.case.name]=events
        met=metrics(run,events;freq,scope);push!(mets,met);numerical_checks!(r,run,met)
        detailed_event_table!(r,prefix*"_Events",run,events)
        chatter_table!(r,prefix*"_Recontacts",run,events)
        tr=trace(run,range(run.t0,run.tf;length=scope=="drive" ? 4001 : 2001));trace_table!(r,prefix*"_Full",tr)
        if save_accepted
            table!(r,prefix*"_Accepted",["time_s","contact_side",HEAD...],
                ([seg.sol.t[j],seg.phase,(seg.sol.u[j].*run.scale)...] for seg in run.segments for j in eachindex(seg.sol.t)))
        end
        if scope=="drive"
            guarded("Plot full states/forces: $(run.model.case.name)",()->drive_figure!(r,run,tr,"full",run.model.case.name*" — full drive"))
            guarded("Plot energy states: $(run.model.case.name)",()->ledger_figure!(r,run,tr))
            late=trace(run,range(max(run.t0,run.tf-2/freq),run.tf;length=12001));trace_table!(r,prefix*"_Late",late)
            guarded("Plot late states/forces: $(run.model.case.name)",()->drive_figure!(r,run,late,"last_two",run.model.case.name*" — final two drive cycles"))
            guarded("Plot late energy states: $(run.model.case.name)",()->ledger_figure!(r,run,late;suffix="energy_last_two"))
            for side in (-1,1)
                times=episode_times(run,events;side)
                if !isempty(times)
                    ep=trace(run,times);trace_table!(r,prefix*(side<0 ? "_ContactMinus" : "_ContactPlus"),ep)
                    guarded("Plot contact detail: $(run.model.case.name), side $side",()->contact_figure!(r,run,ep,events;side))
                end
            end
        end
        println("    $(run.model.case.name): $(met.accepted_states) accepted states, $(met.events) events; energy fraction $(met.energy_residual_fraction)");flush(stdout)
        events
    end
    guarded("Check conservation identities and original-model agreement",()->append!(r.checks,invariant_checks(panels)))
    guarded("Independent beam, force-gap, holding, and damping tests",()->merge!(r.tables,static_tables(panels)))
    if haskey(r.tables,"BeamFEM")
        fem=r.tables["BeamFEM"].rows
        check!(r,"Independent 64-element beam static stiffness",fem[end][4],1e-4,fem[end][4]<1e-4)
        df=abs(fem[end][5]/fem[end-1][5]-1)
        check!(r,"Beam first frequency: 32 to 64 elements",df,1e-4,df<1e-4)
    end
    guarded("Static diagnostic figures",()->static_figures!(r))
    dcrows=Any[];dctransitions=Any[]
    for el in [:legacy,:geometric]
        guarded("DC branch search: $el",()->begin
            s=dc_sweep(plant(Case(name="dc_"*string(el),contact=:constraint,electric=el);panels))
            append!(dcrows,s.rows);append!(dctransitions,s.transitions)
        end)
    end
    table!(r,"DCBranches",["case","direction","bias_V","x1_m","x2_m","branch","reaction_N","force_residual_N","minimum_free_or_held_stiffness_N_m","capacitance_F","effective_potential_J"],dcrows)
    table!(r,"DCTransitionBrackets",["case","direction","from","to","lower_V","upper_V","voltage_grid_V"],dctransitions)
    guarded("Plot DC branches",()->dc_figure!(r))
    base=Case(name="baseline",alpha=Float64(alpha),bias=Float64(bias))
    ref=change(base;name="constrained_reference",contact=:constraint,electric=:geometric)
    drivecases=[base,ref,change(base;name="baseline_1g",alpha=1.),change(ref;name="reference_1g",alpha=1.)]
    for (i,c) in enumerate(drivecases)
        guarded("Full drive $(i)/4: $(c.name)",()->begin
            run=compute(c);push!(fullruns,run);record(run,"drive";prefix="D$i")
        end)
    end
    if length(fullruns)==4
        guarded("Compare full-drive results",()->comparison_figure!(r,fullruns,evmap,freq))
    end
    bidx=findfirst(x->x.model.case.name=="baseline",fullruns)
    if bidx!==nothing && haskey(evmap,"baseline")
        br=fullruns[bidx];seq=last_sequence(evmap["baseline"],1)
        if !isempty(seq)
            t0=max(br.t0,seq[1].time-150e-6);tf=min(br.tf,seq[end].time+750e-6)
            u0=at(br,t0)[1]
            cases=[change(base;name="replay_baseline"),change(base;name="electric_only",electric=:geometric),
                change(base;name="constraint_only",contact=:constraint),change(ref;name="constraint_and_electric"),
                change(base;name="beam_half",bending_scale=.5),change(base;name="beam_double",bending_scale=2.),
                change(base;name="structural_damping",ce=2e-4),change(base;name="wall_10x",kw_scale=10.),
                change(base;name="impact_damping_5x",cw_scale=5.),change(base;name="wall_width_quarter",eps_wall=.125e-9),
                change(base;name="gap_width_quarter",eps_gap=.5e-9),change(base;name="seal_width_quarter",seal_width=6.25e-9),
                change(base;name="fluid_seal_off",seal=false),change(ref;name="physical_gap_5nm",contact_gap=5e-9),
                change(ref;name="reference_bias_double",bias=2Float64(bias)),change(ref;name="reference_slope_half",slope_scale=.5)]
            for (i,c) in enumerate(cases)
                guarded("Contact replay $(i)/$(length(cases)): $(c.name)",()->begin
                    run=compute(c;u0,tspan=(t0,tf));push!(replays,run);ev=record(run,"replay";prefix="P$i")
                    u,ph=at(run,tf);censored=ph!=0 || abs(u[2])>=run.model.gc
                    push!(replaywindows,(c.name,t0,tf,censored ? "CONTACT AT WINDOW END" : "FREE AT WINDOW END",u0[1:5]...))
                end)
            end
            guarded("Plot encounter sensitivities",()->replay_figures!(r,replays,evmap,mets))
            if !isempty(replays) && replays[1].model.case.name=="replay_baseline"
                rr=replays[1];reference_events=evmap[rr.model.case.name]
                guarded("Double crossing-search subdivisions",()->begin
                    refined=crossings(rr;nsub=8)
                    check!(r,"Replay crossing counts: 4 vs 8 subdivisions",length(refined)-length(reference_events),0,length(refined)==length(reference_events))
                    event_table!(r,"Crossings8Subdivisions",refined)
                end)
                conv=Any[];ts=collect(range(t0,tf;length=6001))
                refinements=[("penalty_time",1,panels,rtol/5,dtmax/2),("penalty_space",1,2panels,rtol,dtmax),
                    ("constraint_time",4,panels,rtol/5,dtmax/2),("constraint_space",4,2panels,rtol,dtmax)]
                for (label,idx,np,tol,step) in refinements
                    guarded("Numerical convergence: $label",()->begin
                        source=replays[idx]
                        reference_trace=trace(source,ts);reference_events=evmap[source.model.case.name]
                        run=compute(change(source.model.case;name=label);u0,tspan=(t0,tf),np,tol,step)
                        ev=record(run,"convergence";prefix=label)
                        fine=trace(run,ts);dx=maximum(abs.(fine.u[1:2,:]-reference_trace.u[1:2,:]));dV=maximum(abs.(fine.u[5,:]-reference_trace.u[5,:]))
                        eventerr=length(ev)==length(reference_events) && !isempty(ev) ? maximum(abs.([ev[j].time-reference_events[j].time for j in eachindex(ev)])) : NaN
                        push!(conv,(label,np,tol,step,dx,dV,length(ev),eventerr))
                        check!(r,"$label: event count",length(ev)-length(reference_events),0,length(ev)==length(reference_events))
                        check!(r,"$label: max position difference",dx,2e-9,dx<2e-9)
                        check!(r,"$label: max voltage difference",dV,0.002,dV<0.002)
                    end)
                end
                table!(r,"Convergence",["case","panels","rtol","dtmax_s","max_position_difference_m","max_voltage_difference_V","event_count","max_matched_event_time_difference_s"],conv)
                guarded("Plot numerical convergence",()->convergence_figure!(r,conv))
            end
        else
            push!(r.errors,("Contact replays","No positive contact sequence found at the selected drive; contact tests are inapplicable for this run."))
        end
    end
    if !isempty(mets);table!(r,"Metrics",string.(collect(keys(mets[1]))),[Tuple(x) for x in mets]);end
    table!(r,"CaseParameters",["case","parameter","value"],paramrows)
    table!(r,"ReplayWindows",["case","start_s","end_s","window_end_status","initial_x1_m","initial_x2_m","initial_v1_m_s","initial_v2_m_s","initial_Vout_V"],replaywindows)
    failed=count(x->x[4]=="FAIL",r.checks)
    status=isempty(r.errors) && failed==0 ? "NUMERICAL_CHECKS_PASSED" : "REQUIRES_REVIEW"
    table!(r,"Summary",["item","value"],[("status",status),("physical_validation","NOT CALIBRATED"),("successful_simulations",length(mets)),
        ("figures_PNG_and_PDF",length(r.figures)),("failed_numerical_checks",failed),("failed_stages",length(r.errors)),
        ("next_step","Compare measured contact dwell, release displacement, waveform, circuit voltage and pull-in/release with this same forcing definition.")])
    table!(r,"PlotIndex",["basename","title","caption"],[(f.name,f.title,f.caption) for f in r.figures])
    println("Writing consolidated workbook ($(length(r.tables)) sheets before index/checks)...");flush(stdout)
    workbook=joinpath(dir,"results.xlsx");save_workbook(r,workbook);gallery!(r)
    println("Finished: $status; $(length(r.figures)) figure pairs, $(length(mets)) simulations.")
    println("Workbook: $workbook\nGallery: $(joinpath(dir,"plots","index.html"))")
    strict && status!="NUMERICAL_CHECKS_PASSED" && error("Tests require review; inspect Checks and Errors in $workbook")
    (;directory=dir,workbook,gallery=joinpath(dir,"plots","index.html"),status,metrics=mets,checks=r.checks,errors=r.errors)
end
