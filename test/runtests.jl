using Xpress, Base.Test, MathOptInterface, MathOptInterfaceTests, MathOptInterfaceXpress

const MOIT = MathOptInterfaceTests
const MOIXPR = MathOptInterfaceXpress

@testset "MathOptInterfaceXpress" begin
    @testset "Linear tests" begin
        linconfig = MOIT.TestConfig()
        solverf() = MOIXPR.XpressSolverInstance()
        MOIT.contlineartest(solverf , linconfig, ["linear12","linear8a","linear8b","linear8c"])
        
        solverf_nopresolve() = MOIXPR.XpressSolverInstance(PRESOLVE = 0)
        MOIT.contlineartest(solverf_nopresolve , linconfig, ["linear12"])

        linconfig_nocertificate = MOIT.TestConfig(infeas_certificates = false)
        MOIT.linear12test(solverf, linconfig_nocertificate)
    end

    @testset "Quadratic tests" begin
        quadconfig = MOIT.TestConfig(atol = 1e-5, rtol = 1e-5, duals = false)
        solverf() = MOIXPR.XpressSolverInstance()
        MOIT.contquadratictest(solverf, quadconfig)
    end

    @testset "Linear Conic tests" begin
        linconfig = MOIT.TestConfig()
        solverf() = MOIXPR.XpressSolverInstance()
        MOIT.lintest(solverf, linconfig, ["lin3","lin4"])

        solverf_nopresolve() = MOIXPR.XpressSolverInstance(PRESOLVE = 0)
        MOIT.lintest(solverf_nopresolve, linconfig)
    end

    @testset "Integer Linear tests" begin
        intconfig = MOIT.TestConfig()
        solverf() = MOIXPR.XpressSolverInstance()
        MOIT.intlineartest(solverf, intconfig)
    end
end
;
