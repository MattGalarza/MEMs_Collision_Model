# Focused checks of crossing semantics, episode grouping, and run preservation.
include(joinpath(@__DIR__,"collision_model_corrected.jl"))
using .CorrectedMEMS, Test
const CM=CorrectedMEMS
struct TestSolution
    t::Vector{Float64}
    f::Function
end
(s::TestSolution)(t)=s.f(t)
@testset "Strict boundary crossings" begin
    t=collect(0.:0.1:1.)
    e=CM.level_crossings(t,t.-0.5,0.,x->x-0.5)
    @test length(e)==1
    @test e[1][1] ≈ 0.5 atol=1e-12
    @test e[1][2]==1
    @test isempty(CM.level_crossings(t,(t.-0.5).^2,0.,x->(x-0.5)^2))
    @test isempty(CM.level_crossings(t,t,0.,identity)) # endpoint is not a crossing
    @test isempty(CM.level_crossings(t,zero(t),0.,x->0.))
    @test length(CM.level_crossings([0.,0.4,0.5,0.6,1.],[-1.,0.,0.,0.,1.],0.,x->x<0.4 ? -1. : x>0.6 ? 1. : 0.))==1
end
@testset "Contact oscillation versus repeated contact loss" begin
    m=Model(;panels=32);o=PlotOptions(contact_subdivisions=8,contact_band=5e-9)
    T=100e-6
    function synthetic(offset,amplitude)
        f(t)=begin
            d=offset+amplitude*cos(2pi*t/T)
            u=zeros(12);u[1]=m.gc;u[2]=m.gc+d;u[4]=-amplitude*2pi/T*sin(2pi*t/T);u
        end
        sol=TestSolution(collect(range(0.,3T;length=31)),f)
        (;sol,scale=ones(12),model=m,kind=:probe,freq=0.,book=CM.ReviewBook())
    end
    ringing=synthetic(3e-9,1e-9)
    scan=CM.scan_contact(ringing,o)
    @test isempty(filter(e->e.level==0,scan.events))
    @test count(e->e.delta>0,scan.turning)>=5
    episodes=CM.contact_episodes(ringing,scan,o)
    @test length(episodes)==1
    @test episodes[1].left_truncated && episodes[1].right_truncated
    @test episodes[1].contact_time ≈ 3T
    recross=synthetic(0.,2e-9)
    scan=CM.scan_contact(recross,o)
    raw=filter(e->e.level==0,scan.events)
    @test length(raw)==6
    @test count(e->e.direction==1,raw)==3
    @test count(e->e.direction==-1,raw)==3
    @test length(CM.deadband_events(recross,scan,0.1e-9))==6
    tiny=synthetic(0.,0.02e-9)
    tiny_scan=CM.scan_contact(tiny,o)
    @test length(filter(e->e.level==0,tiny_scan.events))==6
    @test isempty(CM.deadband_events(tiny,tiny_scan,0.1e-9))
    double=CM.scan_contact(recross,o;nsub=16)
    @test CM.event_agreement(raw,filter(e->e.level==0,double.events)).same
    @test !CM.event_agreement(raw,raw[2:end]).same
end
@testset "Run names preserve prior outputs" begin
    mktempdir() do root
        a=CM.run_directory(root,4.95,3.,20.)
        b=CM.run_directory(root,4.95,3.,20.)
        c=CM.run_directory(root,4.95,3.,20.)
        @test basename(a)=="RUN_4.95g_3V_20Hz"
        @test basename(b)=="RUN_4.95g_3V_20Hz__02"
        @test basename(c)=="RUN_4.95g_3V_20Hz__03"
        @test length(readdir(root))==3
        @test basename(CM.run_directory(root,0.,3.,0.;kind=:probe))=="RUN_0g_3V_0Hz_PROBE"
    end
end
