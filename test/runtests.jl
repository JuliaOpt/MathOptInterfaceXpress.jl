using Xpress, Base.Test, MathOptInterface, MathOptInterface.Test, MathOptInterfaceXpress

const MOIT = MathOptInterface.Test
const MOIXPR = MathOptInterfaceXpress

@testset "MathOptInterfaceXpress" begin
    @testset "Linear tests" begin
        linconfig = MOIT.TestConfig()
        solver = XpressOptimizer()
        MOIT.contlineartest(solver , linconfig, ["linear10","linear12","linear8a","linear8b","linear8c"])
        
        solver_nopresolve = XpressOptimizer(PRESOLVE = 0)
        MOIT.contlineartest(solver_nopresolve , linconfig, ["linear10","linear12"])

        linconfig_nocertificate = MOIT.TestConfig(infeas_certificates=false)
        MOIT.linear12test(solver, linconfig_nocertificate)

        # Intervals
        linconfig_noquery = MOIT.TestConfig(query=false)
        MOIT.linear10test(solver, linconfig_noquery)
    end

    @testset "Quadratic tests" begin
        quadconfig = MOIT.TestConfig(atol=1e-5, rtol=1e-5, duals=false, query=false)
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
        MOIT.validtest(solver)
        MOIT.emptytest(solver)
        solver2 = XpressOptimizer()
        MOIT.copytest(solver,solver2)
    end
end
;
