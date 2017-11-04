using Xpress, Base.Test, MathOptInterface, MathOptInterfaceXpress
include(joinpath(Pkg.dir("MathOptInterface"), "test", "contlinear.jl"))
include(joinpath(Pkg.dir("MathOptInterface"), "test", "intlinear.jl"))
include(joinpath(Pkg.dir("MathOptInterface"), "test", "contconic.jl"))
include(joinpath(Pkg.dir("MathOptInterface"), "test", "contquadratic.jl"))


# contlinear
linear1test(MathOptInterfaceXpress.XpressSolver())
linear2test(MathOptInterfaceXpress.XpressSolver())
linear3test(MathOptInterfaceXpress.XpressSolver())
linear4test(MathOptInterfaceXpress.XpressSolver())
linear5test(MathOptInterfaceXpress.XpressSolver())
linear6test(MathOptInterfaceXpress.XpressSolver())
linear7test(MathOptInterfaceXpress.XpressSolver())
linear8test(MathOptInterfaceXpress.XpressSolver(PRESOLVE = 0)) # infeasible/unbounded
linear9test(MathOptInterfaceXpress.XpressSolver())
linear10test(MathOptInterfaceXpress.XpressSolver())
linear11test(MathOptInterfaceXpress.XpressSolver())

# intlinear
knapsacktest(MathOptInterfaceXpress.XpressSolver())
int1test(MathOptInterfaceXpress.XpressSolver())
int2test(MathOptInterfaceXpress.XpressSolver()) # SOS
int3test(MathOptInterfaceXpress.XpressSolver())

# contconic
lin1tests(MathOptInterfaceXpress.XpressSolver())
lin2tests(MathOptInterfaceXpress.XpressSolver())
lin3test(MathOptInterfaceXpress.XpressSolver(PRESOLVE = 0)) # infeasible
lin4test(MathOptInterfaceXpress.XpressSolver(PRESOLVE = 0)) # infeasible

# contquadratic
qp1test(MathOptInterfaceXpress.XpressSolver(), atol = 1e-5)
qp2test(MathOptInterfaceXpress.XpressSolver(), atol = 1e-5)
qp3test(MathOptInterfaceXpress.XpressSolver(), atol = 1e-5)
qcp1test(MathOptInterfaceXpress.XpressSolver(), atol = 1e-5)
qcp2test(MathOptInterfaceXpress.XpressSolver(), atol = 1e-5)
qcp3test(MathOptInterfaceXpress.XpressSolver(), atol = 1e-5)
socp1test(MathOptInterfaceXpress.XpressSolver(), atol = 1e-5)
