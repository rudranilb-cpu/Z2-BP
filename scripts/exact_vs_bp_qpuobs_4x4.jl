using Pkg
Pkg.activate(".")
using Z2GaugeBenchmarks
using Printf

const Z2 = Z2GaugeBenchmarks

function main()
    common = (
        nx=4, ny=4, K=2.2, Gamma=1.0,
        dt=0.05, steps=10,
        trotter=:lie_x_zz,
        boundary=:open_dual,
        defect_x=2, defect_y=2,
        maxdim=16, cutoff=1e-10,
        normalize_tensors=true,
        bp_tolerance=1e-6, bp_maxiter=20,
        bp_fail_on_nonconvergence=false,
    )

    cb = SimulationConfig(; common..., backend=:bp)
    ce = SimulationConfig(; common..., backend=:exact)

    r0 = resolved_defect(cb)
    targets = Tuple{String,Int,Tuple{Int,Int}}[]
    for (dir,dx,dy) in [("R",1,0),("L",-1,0),("U",0,1),("D",0,-1)]
        for d in 1:2
            v=(r0[1]+d*dx,r0[2]+d*dy)
            1<=v[1]<=cb.nx && 1<=v[2]<=cb.ny && push!(targets,(dir,d,v))
        end
    end

    pb,_ = Z2.initialize_bp(cb)
    pe = Z2.initialize_exact(ce)
    circuit = build_trotter_layer(cb,pb.graph)

    f=open("results/exact_vs_bp_qpuobs_4x4.csv","w")
    println(f,"step,time,id,observable,exact_delta,bp_delta,abs_error")

    for step in 0:cb.steps
        t=step*cb.dt

        if iseven(step)
            for (dir,d,v) in targets
                # Local X at target
                ex = Z2._expect_x(pe.loop,(v,),ce)-Z2._expect_x(pe.scv,(v,),ce)
                bx = real(only(Z2.TNQS.expect(pb.loop,[("X",[v])]))) -
                     real(only(Z2.TNQS.expect(pb.scv,[("X",[v])])))
                @printf(f,"%d,%.8f,X_%s%d,X,%.12f,%.12f,%.12e\n",
                        step,t,dir,d,ex,bx,abs(bx-ex))

                # Electric-string endpoint ZZ
                ez = Z2._expect_z(pe.loop,(r0,v),ce)-Z2._expect_z(pe.scv,(r0,v),ce)
                bz = real(only(Z2.TNQS.expect(pb.loop,[("ZZ",[r0,v])]))) -
                     real(only(Z2.TNQS.expect(pb.scv,[("ZZ",[r0,v])])))
                @printf(f,"%d,%.8f,ZZ_%s%d,ZZ,%.12f,%.12f,%.12e\n",
                        step,t,dir,d,ez,bz,abs(bz-ez))

                # Connected XX
                function cexact(p)
                    Z2._expect_x(p,(r0,v),ce) -
                    Z2._expect_x(p,(r0,),ce)*Z2._expect_x(p,(v,),ce)
                end

                valsL=real.(Z2.TNQS.expect(pb.loop,[("X",[r0]),("X",[v]),("XX",[r0,v])]))
                valsS=real.(Z2.TNQS.expect(pb.scv, [("X",[r0]),("X",[v]),("XX",[r0,v])]))
                ec=cexact(pe.loop)-cexact(pe.scv)
                bc=(valsL[3]-valsL[1]*valsL[2])-(valsS[3]-valsS[1]*valsS[2])

                @printf(f,"%d,%.8f,CXX_%s%d,CXX,%.12f,%.12f,%.12e\n",
                        step,t,dir,d,ec,bc,abs(bc-ec))
            end
            flush(f)
            @printf("validated observables at t=%.2f\n",t)
        end

        step==cb.steps && break
        pe,_ = Z2.evolve_exact!(pe,ce)
        pb,_ = Z2.evolve_bp!(pb,circuit,cb)
    end

    close(f)
    println("DONE: results/exact_vs_bp_qpuobs_4x4.csv")
end

main()
