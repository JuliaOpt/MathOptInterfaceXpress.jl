#push!(Base.LOAD_PATH,joinpath(dirname(@__FILE__),"..",".."))

using Xpress, Base.Test, MathOptInterface, MathOptInterface.Test, MathOptInterfaceXpress

const MOI = MathOptInterface
const MOIT = MathOptInterface.Test
const MOIXPR = MathOptInterfaceXpress

@testset "MathOptInterfaceXpress" begin
    @testset "Unit Tests" begin
        config = MOIT.TestConfig()
        solver = XpressOptimizer()
        MOIT.basic_constraint_tests(solver, config;
            exclude = [
            ]
        )
        MOIT.unittest(solver, config, [
            "solve_affine_interval"
        ])
    end
    @testset "Linear tests" begin
        linconfig = MOIT.TestConfig()
        solver = XpressOptimizer()
        MOIT.contlineartest(solver , linconfig, ["linear12","linear8a","linear8b","linear8c"])
        
        solver_nopresolve = XpressOptimizer(PRESOLVE = 0)
        MOIT.contlineartest(solver_nopresolve , linconfig, ["linear12"])

        linconfig_nocertificate = MOIT.TestConfig(infeas_certificates=false)
        MOIT.linear12test(solver, linconfig_nocertificate)
    end

    @testset "Quadratic tests" begin
        quadconfig = MOIT.TestConfig(atol=1e-5, rtol=1e-5, duals=false)
        solver = XpressOptimizer()
        MOIT.contquadratictest(solver, quadconfig)
    end

    @testset "Linear Conic tests" begin
        linconfig = MOIT.TestConfig()
        solver = XpressOptimizer()
        MOIT.lintest(solver, linconfig, ["lin3","lin4"])

        solver_nopresolve = XpressOptimizer(PRESOLVE=0)
        MOIT.lintest(solver_nopresolve, linconfig)
    end

    @testset "Integer Linear tests" begin
        intconfig = MOIT.TestConfig()
        solver = XpressOptimizer()
        MOIT.intlineartest(solver, intconfig)
    end

    @testset "ModelLike tests" begin
        intconfig = MOIT.TestConfig()
        solver = XpressOptimizer()
        MOIT.nametest(solver)
        @testset "validtest" begin
            MOIT.validtest(solver)
        end
        @testset "emptytest" begin
            MOIT.emptytest(solver)
        end
        @testset "orderedindicestest" begin
            MOIT.orderedindicestest(solver)
        end
        @testset "canaddconstrainttest" begin
            MOIT.canaddconstrainttest(solver, Float64, Complex{Float64})
        end
        @testset "copytest" begin
            solver2 = XpressOptimizer()
            MOIT.copytest(solver,solver2)
        end
    end
end
;
