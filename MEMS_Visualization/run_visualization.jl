using Printf, LinearAlgebra
using SciMLBase, OrdinaryDiffEqRosenbrock, ADTypes
include(joinpath(@__DIR__,"collision_model_corrected.jl"))
using .CorrectedMEMS
const OUT=joinpath(@__DIR__,"results")
mkpath(OUT)
function main()
println("Running reference model: 8 cycles, 20 Hz, 4.95 g, 3 V.");flush(stdout)
runmodel=simulate(kind=:drive,cycles=8,outdir=OUT,tag="visualization")
sol=runmodel.sol; scale=runmodel.scale; m=runmodel.model; p=m.p
physical(t)=sol(t).*scale
T=1/20; tf=sol.t[end]; E0=energy(m,physical(0.0))
header=["t_s","t_ms","tau_us","x1_um","x2_um","tip_overlap_nm","bend_nm",
    "v1_mm_s","v2_mm_s","Vout_mV","C_pF","PR_pW","Qe2_uN","Qwall2_uN",
    "Qbend2_uN","Qfilm2_uN","ER_pJ","Wbase_nJ","Wbias_nJ","Dfilm_nJ",
    "Dwall_nJ","Dstruct_nJ","deltaE_nJ","residual_aJ","a_g"]
function sample_row(t,tref)
    u=physical(t); x1,x2,v1,v2,vo=u[1:5]
    c=constitutive(m,x1,x2); w=CorrectedMEMS.wall(m,x2,v2)
    ramp=0.5*(1-cos(pi*min(t*20/4,1.0)))
    (t,t*1e3,(t-tref)*1e6,x1*1e6,x2*1e6,(abs(x2)-m.gc)*1e9,
     (x2-x1)*1e9,v1*1e3,v2*1e3,vo*1e3,c.C*1e12,vo^2/p.Rload*1e12,
     0.5*(p.Vbias-vo)^2*c.grad[2]*1e6,w.F*1e6,
     -p.n_beams*m.ke*(x2-x1)*1e6,-(c.D*[v1,v2])[2]*1e6,
     u[8]*1e12,u[6]*1e9,u[7]*1e9,u[9]*1e9,u[11]*1e9,u[10]*1e9,
     (energy(m,u)-E0)*1e9,CorrectedMEMS.ledger(m,u,E0)*1e18,
     4.95*ramp*sin(2*pi*20*t))
end
# Locate the first positive-wall entry in the final drive cycle.
entries=[j for j in 2:length(sol.t) if sol.t[j]>=tf-T &&
    sol.u[j-1][2]*scale[2]<m.gc && sol.u[j][2]*scale[2]>=m.gc]
@assert !isempty(entries) "No positive contact entry in final cycle"
j=first(entries); lo=sol.t[j-1];hi=sol.t[j]
for k in 1:50
    mid=(lo+hi)/2
    if physical(mid)[2]>m.gc;hi=mid;else;lo=mid;end
end
tc=(lo+hi)/2
println("Contact entry at t = ",tc," s. Exporting dense diagnostics.");flush(stdout)
for (name,times,ref) in [
    ("overview",collect(range(0.0,tf;length=8001)),0.0),
    ("last_cycle",sort(unique(vcat(collect(range(tf-T,tf;length=12001)),
        [t for t in sol.t if t>=tf-T]))),tf-T),
    ("contact",collect(range(tc-80e-6,tc+320e-6;length=4001)),tc)]
    CorrectedMEMS.writecsv(joinpath(OUT,name*".csv"),header,
        (sample_row(t,ref) for t in times))
end
# Extrema use accepted solver states; exact maxima would require optimization.
us=[z.*scale for z in sol.u]; lastidx=findall(>=(tf-T),sol.t)
lastus=us[lastidx]; Eload=(us[end][8]-physical(tf-4*T)[8])/(4*T)
open(joinpath(OUT,"visualization_metrics.txt"),"w") do io
    println(io,"Julia ",VERSION,"; successful native simulation; unchanged corrected model.")
    for (key,val) in ["contact_entry_time_s"=>tc,
        "last_four_cycles_mean_load_power_W"=>Eload,
        "last_cycle_sampled_max_abs_voltage_V"=>maximum(abs(u[5]) for u in lastus),
        "last_cycle_sampled_max_tip_overlap_m"=>maximum(abs(u[2])-m.gc for u in lastus),
        "last_cycle_sampled_max_abs_beam_deflection_m"=>maximum(abs(u[2]-u[1]) for u in lastus),
        "last_cycle_sampled_max_abs_shuttle_m"=>maximum(abs(u[1]) for u in lastus),
        "total_film_loss_J"=>us[end][9],"total_contact_loss_J"=>us[end][11],
        "total_resistor_energy_J"=>us[end][8],
        "energy_residual_over_throughput"=>runmodel.metrics["energy_residual_over_throughput"]]
        println(io,key," = ",val)
    end
    println(io,"Peak quantities are sampled maxima, not independently optimized extrema.")
    println(io,"Finite reference run, not experimental validation or a stability proof.")
end
println("Running the short contact probe.");flush(stdout)
simulate(kind=:probe,outdir=OUT,tag="probe")
println("Data complete. Compile figure TeX files with pdflatex from this directory.")

end
main()
