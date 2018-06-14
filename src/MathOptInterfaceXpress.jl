__precompile__()
module MathOptInterfaceXpress

export XpressOptimizer

using Xpress
const XPR = Xpress
using MathOptInterface
const MOI = MathOptInterface
using LinQuadOptInterface
const LQOI = LinQuadOptInterface

const SUPPORTED_OBJECTIVES = [
    LQOI.Linear,
    LQOI.Quad
]

const SUPPORTED_CONSTRAINTS = [
    (LQOI.Linear, LQOI.EQ),
    (LQOI.Linear, LQOI.LE),
    (LQOI.Linear, LQOI.GE),
    (LQOI.Linear, LQOI.IV),
    (LQOI.Quad, LQOI.EQ),
    (LQOI.Quad, LQOI.LE),
    (LQOI.Quad, LQOI.GE),
    (LQOI.SinVar, LQOI.EQ),
    (LQOI.SinVar, LQOI.LE),
    (LQOI.SinVar, LQOI.GE),
    (LQOI.SinVar, LQOI.IV),
    (LQOI.SinVar, MOI.ZeroOne),
    (LQOI.SinVar, MOI.Integer),
    (LQOI.VecVar, LQOI.SOS1),
    (LQOI.VecVar, LQOI.SOS2),
    (LQOI.VecVar, MOI.Nonnegatives),
    (LQOI.VecVar, MOI.Nonpositives),
    (LQOI.VecVar, MOI.Zeros),
    (LQOI.VecLin, MOI.Nonnegatives),
    (LQOI.VecLin, MOI.Nonpositives),
    (LQOI.VecLin, MOI.Zeros)
]

mutable struct XpressOptimizer <: LQOI.LinQuadOptimizer
    LQOI.@LinQuadOptimizerBase
    params::Dict{Any,Any}
    XpressOptimizer(::Void) = new()
end

LQOI.LinearQuadraticModel(::Type{XpressOptimizer}, env) = XPR.Model()

function XpressOptimizer(; kwargs...)

    env = nothing
    m = XpressOptimizer(nothing)
    m.params = Dict{Any,Any}()
    MOI.empty!(m)
    for (name,value) in kwargs
        m.params[name] = value
        XPR.setparam!(m.inner, XPR.XPRS_CONTROLS_DICT[name], value)
    end
    return m
end

function MOI.empty!(m::XpressOptimizer) 
    MOI.empty!(m,nothing)
    for (name,value) in m.params
        XPR.setparam!(m.inner, XPR.XPRS_CONTROLS_DICT[name], value)
    end
end

LQOI.supported_constraints(s::XpressOptimizer) = SUPPORTED_CONSTRAINTS
LQOI.supported_objectives(s::XpressOptimizer) = SUPPORTED_OBJECTIVES

backend_type(m::XpressOptimizer, ::MOI.GreaterThan{T}) where T = Cchar('G')
backend_type(m::XpressOptimizer, ::MOI.LessThan{T}) where T    = Cchar('L')
backend_type(m::XpressOptimizer, ::MOI.EqualTo{T}) where T     = Cchar('E')
# Implemented separately
# backend_type(m::XpressOptimizer, ::MOI.Interval{T}) where T    = Cchar('R')

backend_type(m::XpressOptimizer, ::MOI.Zeros)        = Cchar('E')
backend_type(m::XpressOptimizer, ::MOI.Nonpositives) = Cchar('L')
backend_type(m::XpressOptimizer, ::MOI.Nonnegatives) = Cchar('G')



#=
not in LinQuad
=#

setparam!(instance::XpressOptimizer, name, val) = XPR.setparam!(instance.inner, XPR.XPRS_CONTROLS_DICT[name], val)

setlogfile!(instance::XpressOptimizer, path) = XPR.setlogfile(instance.inner, path::String)

cintvec(v::Vector) = convert(Vector{Int32}, v)

#=
    inner wrapper
=#

#=
Constraints
=#

MOI.canset(::XpressOptimizer, ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}) = true
MOI.canset(::XpressOptimizer, ::MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}) = true
MOI.canget(::XpressOptimizer, ::MOI.ConstraintSet, ::Type{LQOI.LCI{LQOI.IV}}) = true

LQOI.change_variable_bounds!(instance::XpressOptimizer, colvec, valvec, sensevec) = XPR.chgbounds!(instance.inner, cintvec(colvec), sensevec, valvec)

LQOI.get_variable_lowerbound(instance::XpressOptimizer, col) = XPR.get_lb(instance.inner, col, col)[1]
LQOI.get_variable_upperbound(instance::XpressOptimizer, col) = XPR.get_ub(instance.inner, col, col)[1]

LQOI.get_number_linear_constraints(instance::XpressOptimizer) = XPR.num_linconstrs(instance.inner)

LQOI.add_linear_constraints!(instance::XpressOptimizer, A::LQOI.CSRMatrix{Float64}, sensevec, rhsvec) = XPR.add_constrs!(instance.inner, A.row_pointers, A.columns, A.coefficients, sensevec, rhsvec)

function LQOI.add_ranged_constraints!(instance::XpressOptimizer, A::LQOI.CSRMatrix{Float64}, lowerbound, upperbound)
    newrows = length(lowerbound)
    rows = XPR.num_linconstrs(instance.inner)
    addedrows = collect((rows+1):(rows+newrows))

    sensevec = fill(Cchar('E'),newrows)
    XPR.add_constrs!(instance.inner, A.row_pointers, A.columns, A.coefficients, sensevec, upperbound)

    XPR.chg_rhsrange!(instance.inner, cintvec(addedrows), +upperbound-lowerbound)
end

function LQOI.modify_ranged_constraints!(instance::XpressOptimizer, rows::Vector{Int}, lowerbound::Vector{Float64}, upperbound::Vector{Float64})
    XPR.set_rhs!(instance.inner, rows, upperbound)
    XPR.chg_rhsrange!(instance.inner, cintvec(rows), +upperbound-lowerbound)
end

LQOI.get_rhs(instance::XpressOptimizer, row) = XPR.get_rhs(instance.inner, row, row)[1]

function LQOI.get_range(instance::XpressOptimizer, row)
    ub = XPR.get_rhs(instance.inner, row, row)[1]
    r = XPR.get_rhsrange(instance.inner, row, row)[1]
    return ub-r, ub
end

# TODO improve
function LQOI.get_linear_constraint(instance::XpressOptimizer, idx)
    A = XPR.get_rows(instance.inner, idx, idx)'
    return A.rowval-1, A.nzval
end

function LQOI.get_quadratic_constraint(instance::XpressOptimizer, idx)
    A = XPR.get_rows(instance.inner, idx, idx)'

    Q = XPR.get_qrowmatrix_triplets(instance.inner, idx)

    return A.rowval-1, A.nzval, Q
end

# notin LQOI
# TODO improve
function getcoef(instance::XpressOptimizer, row, col)
    A = XPR.get_rows(instance.inner, row, row)'
    cols = A.rowval
    vals = A.nzval

    pos = findfirst(cols, col)
    if pos > 0
        return vals[pos]
    else
        return 0.0
    end
end

LQOI.change_matrix_coefficient!(instance::XpressOptimizer, row, col, coef) = XPR.chg_coeffs!(instance.inner, row, col, coef)

LQOI.change_objective_coefficient!(instance::XpressOptimizer, col, coef) = XPR.set_objcoeffs!(instance.inner, Int32[col], Float64[coef])

LQOI.change_rhs_coefficient!(instance::XpressOptimizer, row, coef) = XPR.set_rhs!(instance.inner, Int32[row], Float64[coef])

LQOI.delete_linear_constraints!(instance::XpressOptimizer, rowbeg, rowend) = XPR.del_constrs!(instance.inner, cintvec(collect(rowbeg:rowend)))

LQOI.delete_quadratic_constraints!(instance::XpressOptimizer, rowbeg, rowend) = XPR.del_constrs!(instance.inner, cintvec(collect(rowbeg:rowend)))

LQOI.change_variable_types!(instance::XpressOptimizer, colvec, typevec) = XPR.chgcoltypes!(instance.inner, colvec, typevec)

LQOI.change_linear_constraint_sense!(instance::XpressOptimizer, rowvec, sensevec) = XPR.set_rowtype!(instance.inner, rowvec, sensevec)

LQOI.add_sos_constraint!(instance::XpressOptimizer, colvec, valvec, typ) = XPR.add_sos!(instance.inner, typ, colvec, valvec)

LQOI.delete_sos!(instance::XpressOptimizer, idx1, idx2) = XPR.del_sos!(instance.inner, cintvec(collect(idx1:idx2)))

# TODO improve getting processes
function LQOI.get_sos_constraint(instance::XpressOptimizer, idx)
    A, types = XPR.get_sos_matrix(instance.inner)
    line = A[idx,:] #sparse vec
    cols = line.nzind
    vals = line.nzval
    typ = types[idx] == Cchar('1') ? :SOS1 : :SOS2
    return cols, vals, typ
end

LQOI.get_number_quadratic_constraints(instance::XpressOptimizer) = XPR.num_qconstrs(instance.inner)


function scalediagonal!(V, I, J, scale)
    #  LQOI assumes 0.5 x' Q x, but Gurobi requires the list of terms, e.g.,
    #  2x^2 + xy + y^2, so we multiply the diagonal of V by 0.5. We don't
    #  multiply the off-diagonal terms since we assume they are symmetric and we
    #  only need to give one.
    #
    #  We also need to make sure that after adding the constraint we un-scale
    #  the vector because we can't modify user-data.
    for i in 1:length(I)
        if I[i] == J[i]
            V[i] *= scale
        end
    end
end

function LQOI.add_quadratic_constraint!(instance::XpressOptimizer, cols, coefs, rhs, sense, I, J, V)
    @assert length(I) == length(J) == length(V)
    scalediagonal!(V, I, J, 0.5)
    XPR.add_qconstr!(instance.inner, cols, coefs, I, J, V, sense, rhs)
    scalediagonal!(V, I, J, 2.0)
end
#=
    Objective
=#

function LQOI.set_quadratic_objective!(instance::XpressOptimizer, I, J, V)
    @assert length(I) == length(J) == length(V)
    XPR.delq!(instance.inner)
    # scalediagonal!(V, I, J, 0.5)
    XPR.add_qpterms!(instance.inner, I, J, V)
    # scalediagonal!(V, I, J, 2.0)
    return nothing
end

function LQOI.set_linear_objective!(instance::XpressOptimizer, colvec, coefvec) 
    nvars = XPR.num_vars(instance.inner)
    obj = zeros(Float64, nvars)

    for i in eachindex(colvec)
        obj[colvec[i]] = coefvec[i]
    end

    XPR.set_obj!(instance.inner, obj)
    nothing
end

function LQOI.change_objective_sense!(instance::XpressOptimizer, symbol)
    if symbol == :min
        XPR.set_sense!(instance.inner, :minimize)
    else
        XPR.set_sense!(instance.inner, :maximize)
    end
end

function LQOI.get_linear_objective!(instance::XpressOptimizer, x::Vector{Float64})
    obj = XPR.get_obj(instance.inner)
    @assert length(x) == length(obj)
    for i in 1:length(obj)
        x[i] = obj[i]
    end
end

function LQOI.get_quadratic_terms_objective(instance::XpressOptimizer)
    I, J, V = XPR.getq(instance.inner)
    return I, J, V
end

function LQOI.get_objectivesense(instance::XpressOptimizer)
    s = XPR.model_sense(instance.inner)
    if s == :maximize
        return MOI.MaxSense
    else
        return MOI.MinSense
    end
end

#=
    Variables
=#

LQOI.get_number_variables(instance::XpressOptimizer) = XPR.num_vars(instance.inner)

LQOI.add_variables!(instance::XpressOptimizer, int) = XPR.add_cvars!(instance.inner, zeros(int))

LQOI.delete_variables!(instance::XpressOptimizer, col, col2) = XPR.del_vars!(instance.inner, col)

# LQOI.lqs_addmipstarts(m, colvec, valvec)
function LQOI.add_mip_starts!(instance::XpressOptimizer, colvec, valvec) 
    x = zeros(XPR.num_vars(instance.inner))
    for i in eachindex(colvec)
        x[colvec[i]] = valvec[i]
    end
    XPR.loadbasis(instance.inner, x)
end

#=
    Solve
=#

LQOI.solve_mip_problem!(instance::XpressOptimizer) = XPR.mipoptimize(instance.inner)

LQOI.solve_quadratic_problem!(instance::XpressOptimizer) = ( writeproblem(instance, "db", "l");LQOI.solve_linear_problem!(instance) )

LQOI.solve_linear_problem!(instance::XpressOptimizer) = XPR.lpoptimize(instance.inner)

function LQOI.get_termination_status(instance::XpressOptimizer)
    stat_lp = XPR.get_lp_status2(instance.inner)
    if XPR.is_mip(instance.inner)
        stat_mip = XPR.get_mip_status2(instance.inner)
        if stat_mip == XPR.MIP_NotLoaded
            return MOI.OtherError
        elseif stat_mip == XPR.MIP_LPNotOptimal
            # MIP search incomplete but there is no linear sol
            # return MOI.OtherError
            return MOI.InfeasibleOrUnbounded
        elseif stat_mip == XPR.MIP_NoSolFound
            # MIP search incomplete but there is no integer sol
            other = xprsmoi_stopstatus(instance.inner)
            if other == MOI.OtherError
                return MOI.SlowProgress#OtherLimit
            else 
                return other
            end

        elseif stat_mip == XPR.MIP_Solution
            # MIP search incomplete but there is a solution
            other = xprsmoi_stopstatus(instance.inner)
            if other == MOI.OtherError
                return MOI.OtherLimit
            else 
                return other
            end

        elseif stat_mip == XPR.MIP_Infeasible
            if XPR.hasdualray(instance.inner)
                return MOI.Success
            else
                return MOI.InfeasibleNoResult
            end
        elseif stat_mip == XPR.MIP_Optimal
            return MOI.Success
        elseif stat_mip == XPR.MIP_Unbounded
            if XPR.hasprimalray(instance.inner)
                return MOI.Success
            else
                return MOI.UnboundedNoResult
            end
        end
        return MOI.OtherError
    else
        if stat_lp == XPR.LP_Unstarted
            return MOI.OtherError
        elseif stat_lp == XPR.LP_Optimal
            return MOI.Success
        elseif stat_lp == XPR.LP_Infeasible
            if XPR.hasdualray(instance.inner)
                return MOI.Success
            else
                return MOI.InfeasibleNoResult
            end
        elseif stat_lp == XPR.LP_CutOff
            return MOI.ObjectiveLimit
        elseif stat_lp == XPR.LP_Unfinished
            return xprsmoi_stopstatus(instance.inner)
        elseif stat_lp == XPR.LP_Unbounded
            if XPR.hasprimalray(instance.inner)
                return MOI.Success
            else
                return MOI.UnboundedNoResult
            end
        elseif stat_lp == XPR.LP_CutOffInDual
            return MOI.ObjectiveLimit
        elseif stat_lp == XPR.LP_Unsolved
            return MOI.OtherError
        elseif stat_lp == XPR.LP_NonConvex
            return MOI.InvalidModel
        end
        return MOI.OtherError
    end
end

function xprsmoi_stopstatus(instance::XpressOptimizer)
    ss = XPR.get_stopstatus(instance.inner)
    if ss == XPR.StopTimeLimit
        return MOI.TimeLimit
    elseif ss == XPR.StopControlC
        return MOI.Interrupted
    elseif ss == XPR.StopNodeLimit
        # should not be here
        warn("should not be here")
        return MOI.NodeLimit
    elseif ss == XPR.StopIterLimit
        return MOI.IterationLimit
    elseif ss == XPR.StopMIPGap
        return MOI.ObjectiveLimit
    elseif ss == XPR.StopSolLimit
        return MOI.SolutionLimit
    elseif ss == XPR.StopUser
        return MOI.Interrupted
    end
    return MOI.OtherError
end

function LQOI.get_primal_status(instance::XpressOptimizer)
    if XPR.is_mip(instance.inner)
        stat_mip = XPR.get_mip_status2(instance.inner)
        if stat_mip in [XPR.MIP_Solution, XPR.MIP_Optimal]
            return MOI.FeasiblePoint
        elseif XPR.MIP_Infeasible && XPR.hasdualray(instance.inner)
            return MOI.InfeasibilityCertificate
        elseif XPR.MIP_Unbounded && XPR.hasprimalray(instance.inner)
            return MOI.InfeasibilityCertificate
        elseif stat_mip in [XPR.MIP_LPOptimal, XPR.MIP_NoSolFound]
            return MOI.InfeasiblePoint
        end
        return MOI.UnknownResultStatus
    else
        stat_lp = XPR.get_lp_status2(instance.inner)
        if stat_lp == XPR.LP_Optimal
            return MOI.FeasiblePoint
        elseif stat_lp == XPR.LP_Unbounded && XPR.hasprimalray(instance.inner)
            return MOI.InfeasibilityCertificate
        # elseif stat_lp == LP_Infeasible
        #     return MOI.InfeasiblePoint - xpress wont return
        # elseif cutoff//cutoffindual ???
        else
            return MOI.UnknownResultStatus
        end
    end
end

function LQOI.get_dual_status(instance::XpressOptimizer) 
    if XPR.is_mip(instance.inner)
        return MOI.UnknownResultStatus
    else
        stat_lp = XPR.get_lp_status2(instance.inner)
        if stat_lp == XPR.LP_Optimal
            return MOI.FeasiblePoint
        elseif stat_lp == XPR.LP_Infeasible && XPR.hasdualray(instance.inner)
            return MOI.InfeasibilityCertificate
        # elseif stat_lp == LP_Unbounded
        #     return MOI.InfeasiblePoint - xpress wont return
        # elseif cutoff//cutoffindual ???
        else
            return MOI.UnknownResultStatus
        end
    end
end

LQOI.get_variable_primal_solution!(instance::XpressOptimizer, place) = XPR.get_solution!(instance.inner, place)

function LQOI.get_linear_primal_solution!(instance::XpressOptimizer, place)
    if num_qconstrs(instance.inner) == 0
        XPR.get_slack_lin!(instance.inner, place)
        rhs = XPR.get_rhs(instance.inner)
        for i in eachindex(place)
            place[i] = -place[i]+rhs[i]
        end
    else
        XPR.get_slack_lin!(instance.inner, place)
        rhs = XPR.get_rhs(instance.inner)
        lrows = XPR.get_lrows(instance.inner)
        for (i,v) in enumerate(lrows)
            place[i] = -place[i]+rhs[v]
        end
    end
    nothing
end

function moi_lrows(m::Xpress.Model)
    tt_rows = collect(1:num_constrs(m))
    if num_qconstrs(m) == 0
        return tt_rows
    else
        return setdiff(tt_rows, XPR.get_qrows(m))
    end
end

function LQOI.get_quadratic_primal_solution!(instance::XpressOptimizer, place)
    if num_qconstrs(instance.inner) == 0
        return nothing
    else
        qrows = XPR.get_qrows(instance.inner)
        newplace = zeros(num_constrs(instance.inner))
        XPR.get_slack!(instance.inner, newplace)
        rhs = XPR.get_rhs(instance.inner)
        for (i,v) in enumerate(qrows)
            place[i] = -newplace[v]+rhs[v]
        end
    end
    nothing
end

LQOI.get_variable_dual_solution!(instance::XpressOptimizer, place) = XPR.get_reducedcost!(instance.inner, place)

LQOI.get_linear_dual_solution!(instance::XpressOptimizer, place) = XPR.get_dual_lin!(instance.inner, place)

function LQOI.get_quadratic_dual_solution!(instance::XpressOptimizer, place)
    if num_qconstrs(instance.inner) == 0
        return nothing
    else
        qrows = XPR.get_qrows(instance.inner)
        newplace = zeros(num_constrs(instance.inner))
        XPR.get_dual!(instance.inner, newplace)
        for (i,v) in enumerate(qrows)
            place[i] = newplace[v]
        end
    end
end

LQOI.get_objective_value(instance::XpressOptimizer) = XPR.get_objval(instance.inner)

LQOI.get_objective_bound(instance::XpressOptimizer) = XPR.get_bestbound(instance.inner)

function LQOI.get_relative_mip_gap(instance::XpressOptimizer)
    best_feasible_solution = XPR.get_mip_objval(instance.inner)
    best_prossible_solution = XPR.get_bestbound(instance.inner)
    return abs((Ubest_prossible_solution-best_feasible_solution)/best_prossible_solution)
end

LQOI.get_iteration_count(instance::XpressOptimizer)  = XPR.get_simplex_iter_count(instance.inner)

LQOI.get_barrier_iterations(instance::XpressOptimizer) = XPR.get_barrier_iter_count(instance.inner)

LQOI.get_node_count(instance::XpressOptimizer) = XPR.get_node_count(instance.inner)

LQOI.get_farkas_dual!(instance::XpressOptimizer, place) = XPR.getdualray!(instance.inner, place)

LQOI.get_unbounded_ray!(instance::XpressOptimizer, place) = XPR.getprimalray!(instance.inner, place)


MOI.free!(m::XpressOptimizer) = XPR.free_model(m.onner)

"""
    writeproblem(m: :MOI.AbstractOptimizer, filename::String)
Writes the current problem data to the given file.
Supported file types are solver-dependent.
"""
writeproblem(instance::XpressOptimizer, filename::String, flags::String="") = XPR.write_model(instance.inner, filename, flags)

end # module