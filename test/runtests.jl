using Xpress, Base.Test, MathOptInterface, MathOptInterfaceTests, MathOptInterfaceXpress

const MOIT = MathOptInterfaceTests
const MOIXPR = MathOptInterfaceXpress

@testset "MathOptInterfaceXpress" begin
    @testset "Linear tests" begin
        linconfig = MOIT.TestConfig(1e-8,1e-8,true,true,true)
        solver = MOIXPR.XpressSolver()
        MOIT.contlineartest(solver , linconfig, ["linear12","linear8a","linear8b","linear8c"])
        
        solver_nopresolve = MOIXPR.XpressSolver(PRESOLVE = 0)
        MOIT.contlineartest(solver_nopresolve , linconfig, ["linear12"])

        linconfig_nocertificate = MOIT.TestConfig(1e-8,1e-8,true,true,false)
        MOIT.linear12test(solver, linconfig_nocertificate)
    end

    @testset "Quadratic tests" begin
        quadconfig = MOIT.TestConfig(1e-5,1e-5,false,true,true)
        solver = MOIXPR.XpressSolver()
        MOIT.contquadratictest(solver, quadconfig)
    end

    @testset "Linear Conic tests" begin
        linconfig = MOIT.TestConfig(1e-8,1e-8,true,true,true)
        solver = MOIXPR.XpressSolver()
        MOIT.lintest(solver, linconfig, ["lin3","lin4"])

        solver_nopresolve = MOIXPR.XpressSolver(PRESOLVE = 0)
        MOIT.lintest(solver_nopresolve, linconfig)
    end

    @testset "Integer Linear tests" begin
        intconfig = MOIT.TestConfig(1e-8,1e-8,true,true,true)
        solver = MOIXPR.XpressSolver()
        MOIT.intlineartest(solver, intconfig)
    end
end
;
