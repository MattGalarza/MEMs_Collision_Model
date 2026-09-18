#!/usr/bin/env julia
# Checks for the worked derivation; Julia standard libraries only.
include(joinpath(@__DIR__, "collision_model_corrected.jl"))
using .CorrectedMEMS, LinearAlgebra, Test, Printf
const CM=CorrectedMEMS
function run_checks()
    m=Model(); p=m.p; d=p.wb-p.wt
    ke=p.E*p.Tf*p.wt^2*d^3/(6*p.Lf^3*((p.wb-3*p.wt)*d+2*p.wt^2*log(p.wb/p.wt)))
    function phi_closed(s)
        b=p.wt+d*s/p.Lf
        J=1+p.wb/(2*b)+b*p.wb/(2*p.wt^2)-(b+p.wb)/p.wt+log(b/p.wt)
        12*m.ke*p.Lf^3/(p.E*p.Tf*d^3)*J
    end
    quad(f,a=0.0,b=1.0)=CM.quad(f,a,b)
    @testset "Worked derivation: independent formulas" begin
        @test isapprox(ke,m.ke;rtol=1e-12) 
        for s in range(0,p.Lf;length=21)
            @test isapprox(phi_closed(s),CM.shape(p,m.ke,s);atol=2e-14)
        end
        phi(u)=1.5*u^2-0.5*u^3
        @test isapprox(quad(phi),3/8;atol=1e-14)
        @test isapprox(quad(u->phi(u)^2),33/140;atol=1e-14)
        @test isapprox(quad(u->phi(u)*(1-phi(u))),39/280;atol=1e-14)
        @test isapprox(quad(u->(1-phi(u))^2),17/35;atol=1e-14)
        @test isapprox(quad(u->(6-12*u)^2),12;atol=1e-12)
        @test isapprox(quad(u->(6*u-6*u^2)^2),6/5;atol=1e-14)
        @test isapprox(4/8*(6/5)^2,0.72;atol=1e-14)
        basis=(t->2*t^2-3*t+1,t->4*t-4*t^2,t->2*t^2-t)
        for (f,whole,half) in zip(basis,(1/6,2/3,1/6),(5/24,1/3,-1/24))
            @test isapprox(quad(f),whole;atol=1e-14)
            @test isapprox(quad(f,0.0,0.5),half;atol=1e-14)
        end
        h=200e-9; b=p.slip_coefficient*p.mean_free_path; py=-1.2e8
        vel(z)=py/(2*p.eta)*(z^2-h*z-b*h)
        dv(z)=py/(2*p.eta)*(2*z-h)
        @test isapprox(vel(0),b*dv(0);rtol=1e-14)
        @test isapprox(vel(h),-b*dv(h);rtol=1e-14)
        @test isapprox(quad(vel,0.0,h),-h^2*(h+6*b)*py/(12*p.eta);rtol=1e-13)
        L=p.Leff; G=h^2*(h+m.kp)
        I0=L/G; I1=L^2/(2*G); I2=L^3/(3*G)
        @test isapprox(12*p.eta*p.Tf*(I2-I1^2/I0),p.eta*p.Tf*L^3/G;rtol=1e-14)
        @test isapprox((6/L)*(L^4/4)/((3/L^2)*(L^5/5)),2.5;rtol=1e-14)
        e=CM.restitution(0.005,50.0); a=0.25; bb=a*e
        @test isapprox(a-log1p(a),-bb-log1p(-bb);atol=1e-15)
        @test isapprox(1/m.Minv[2,2],m.M[2,2]-m.M[1,2]^2/m.M[1,1];rtol=1e-13)
        nodes,weights=CM.gausslegendre(8)
        for degree in 0:15
            exact=isodd(degree) ? 0.0 : 2/(degree+1)
            @test isapprox(sum(weights.*nodes.^degree),exact;atol=1e-14)
        end
    end
    CM.verify(outdir=joinpath(@__DIR__,"check_outputs"))
    println("\nNumerical substitutions (SI):")
    for (name,value) in ("ke"=>m.ke,"k1"=>m.k1,"k3"=>m.k3,"kss"=>m.kss,
        "contact_travel"=>m.gc,"gap_slope"=>m.alpha,"dielectric_gap"=>m.hd,
        "slip_gap"=>m.kp,"total_mass"=>sum(m.M),"M11"=>m.M[1,1],
        "M12"=>m.M[1,2],"M22"=>m.M[2,2],"restitution"=>CM.restitution(0.005,50.0))
        @printf("%-20s %.14g\n",name,value)
    end
    println("\nAll checks passed. No new long-time simulation was run.")
end
open(joinpath(@__DIR__,"derivation_check_results.txt"),"w") do io
    redirect_stdout(io) do
        println("Julia ",VERSION)
        run_checks()
    end
end
print(read(joinpath(@__DIR__,"derivation_check_results.txt"),String))
