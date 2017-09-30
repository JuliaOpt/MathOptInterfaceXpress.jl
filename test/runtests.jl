using Xpress, Base.Test, MathOptInterface, MathOptInterfaceXpress
include(joinpath(Pkg.dir("MathOptInterface"), "test", "contlinear.jl"))
include(joinpath(Pkg.dir("MathOptInterface"), "test", "intlinear.jl"))
include(joinpath(Pkg.dir("MathOptInterface"), "test", "contconic.jl"))
include(joinpath(Pkg.dir("MathOptInterface"), "test", "contquadratic.jl"))


# contlinear
linear1test(XpressSolver())
linear2test(XpressSolver())
linear3test(XpressSolver())
linear4test(XpressSolver())
linear5test(XpressSolver())
linear6test(XpressSolver())
linear7test(XpressSolver())
linear8test(XpressSolver(PRESOLVE = 0)) # infeasible/unbounded
linear9test(XpressSolver())
linear10test(XpressSolver())
linear11test(XpressSolver())

# intlinear
knapsacktest(XpressSolver())
int1test(XpressSolver())
int2test(XpressSolver()) # SOS
int3test(XpressSolver())

# contconic
lin1tests(XpressSolver())
lin2tests(XpressSolver())
lin3test(XpressSolver(PRESOLVE = 0)) # infeasible
lin4test(XpressSolver(PRESOLVE = 0)) # infeasible

# contquadratic
qp1test(XpressSolver(), atol = 1e-5)
qp2test(XpressSolver(), atol = 1e-5)
qp3test(XpressSolver(), atol = 1e-5)
qcp1test(XpressSolver(), atol = 1e-5)
qcp2test(XpressSolver(), atol = 1e-5)
qcp3test(XpressSolver(), atol = 1e-5)
socp1test(XpressSolver(), atol = 1e-5)
