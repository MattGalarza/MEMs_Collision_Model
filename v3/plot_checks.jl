# Checks that plotted quantities reproduce the simulation's equations.
# Run: julia --project=. plot_checks.jl
include(joinpath(@__DIR__, "collision_model_corrected.jl"))
using .CorrectedMEMS, LinearAlgebra, Test
const CM = CorrectedMEMS
@testset "Plot quantities agree with governing equations" begin
    for p in (Params(), Params(c1=1e-5,ce=2e-7,ghs=14.3e-6), Params(wt=30e-6,wb=30e-6))
        m=Model(p;panels=128)
        for (x1,x2,v1,v2,Vo) in [(0.,0.,0.,0.,0.),
            (m.gc-1e-6,m.gc-0.1e-6,0.001,0.002,0.2),
            (m.gc+100e-9,m.gc+20e-9,0.004,-0.03,-0.1),
            (-m.gc-100e-9,-m.gc-20e-9,-0.004,0.03,0.1),
            (14.5e-6,14.4e-6,0.01,0.008,1.0)]
            u=vcat([x1,x2,v1,v2,Vo],zeros(7));du=zeros(12);a=9.5
            CM.rhs!(du,u,(m,t->a),0.0)
            d=CM.review_snapshot(m,u,a,CM.energy(m,u))
            @test d.F[:,11] ≈ m.M*du[3:4] rtol=1e-12 atol=1e-16
            @test vec(sum(d.F[:,1:10];dims=2)) ≈ m.M*du[3:4] rtol=1e-12 atol=1e-16
            @test maximum(abs,d.F[:,12]) < 1e-15
            @test d.diag[9] ≈ CM.energy(m,u) rtol=1e-14
            @test d.diag[11:16] ≈ du[6:11] rtol=1e-13 atol=1e-20
            @test d.diag[14]>=-1e-20 && d.diag[16]>=-1e-20
            @test d.diag[30]+d.diag[31] ≈ d.F[2,6] atol=1e-16
        end
        for x in (-m.gc-10e-9,0.0,m.gc+10e-9)
            c=constitutive(m,x,x);b=[1.,1.];v=0.004
            q=-c.D*b*v
            @test -(b*v)'q >= -1e-20
            @test q ≈ -(-c.D*b*(-v)) atol=1e-16
            @test c.C ≈ constitutive(m,-x,-x).C rtol=1e-12
        end
    end
    # Display reduction must retain short force pulses and their signs.
    t=collect(range(0.0,1.0;length=10001)); y=zeros(length(t)); y[4951]=3.;y[4952]=-4.
    idx=CM.review_indices(t,[y];bins=100)
    @test 4951 in idx && 4952 in idx
    @test issorted(idx) && first(idx)==1 && last(idx)==length(t)
end
