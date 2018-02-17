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

LQOI.LinQuadModel(::Type{XpressOptimizer},env) = XPR.Model()

function XpressOptimizer(;kwargs...)

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

LQOI.lqs_supported_constraints(s::XpressOptimizer) = SUPPORTED_CONSTRAINTS
LQOI.lqs_supported_objectives(s::XpressOptimizer) = SUPPORTED_OBJECTIVES

#=
    inner wrapper
=#

#=
    Main
=#

# LinQuadSolver # Abstract type
# done above

# LQOI.lqs_setparam!(env, name, val)
# TODO fix this one
LQOI.lqs_setparam!(instance::XpressOptimizer, name, val) = XPR.setparam!(instance.inner, XPR.XPRS_CONTROLS_DICT[name], val)

# LQOI.lqs_setlogfile!(env, path)
# TODO fix this one
LQOI.lqs_setlogfile!(instance::XpressOptimizer, path) = XPR.setlogfile(instance.inner, path::String)

# LQOI.lqs_getprobtype(m)
# TODO - consider removing, apparently useless

#=
    Constraints
=#

cintvec(v::Vector) = convert(Vector{Int32}, v)

MOI.canset(::XpressOptimizer, ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}) = true
MOI.canset(::XpressOptimizer, ::MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}) = true


# LQOI.lqs_chgbds!(m, colvec, valvec, sensevec)
LQOI.lqs_chgbds!(instance::XpressOptimizer, colvec, valvec, sensevec) = XPR.chgbounds!(instance.inner, cintvec(colvec), sensevec, valvec)

# LQOI.lqs_getlb(m, col)
LQOI.lqs_getlb(instance::XpressOptimizer, col) = XPR.get_lb(instance.inner, col, col)[1]
# LQOI.lqs_getub(m, col)
LQOI.lqs_getub(instance::XpressOptimizer, col) = XPR.get_ub(instance.inner, col, col)[1]

# LQOI.lqs_getnumrows(m)
LQOI.lqs_getnumrows(instance::XpressOptimizer) = XPR.num_linconstrs(instance.inner)

# LQOI.lqs_addrows!(m, rowvec, colvec, coefvec, sensevec, rhsvec)
LQOI.lqs_addrows!(instance::XpressOptimizer, rowvec, colvec, coefvec, sensevec, rhsvec) = XPR.add_constrs!(instance.inner, rowvec, colvec, coefvec, sensevec, rhsvec)

# LQOI.lqs_getrhs(m, rowvec)
LQOI.lqs_getrhs(instance::XpressOptimizer, row) = XPR.get_rhs(instance.inner, row, row)[1]  

# colvec, coef = LQOI.lqs_getrows(m, rowvec)
# TODO improve
function LQOI.lqs_getrows(instance::XpressOptimizer, idx)
    A = XPR.get_rows(instance.inner, idx, idx)'
    return A.rowval-1, A.nzval
end

# LQOI.lqs_getcoef(m, row, col) #??
# TODO improve
function LQOI.lqs_getcoef(instance::XpressOptimizer, row, col) #??
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

# LQOI.lqs_chgcoef!(m, row, col, coef)
# TODO SPLIT THIS ONE
function LQOI.lqs_chgcoef!(instance::XpressOptimizer, row, col, coef) 
    if row == 0
        XPR.set_objcoeffs!(instance.inner, Int32[col], Float64[coef])
    elseif col == 0
        XPR.set_rhs!(instance.inner, Int32[row], Float64[coef])
    else
        XPR.chg_coeffs!(instance.inner, row, col, coef)
    end
end

# LQOI.lqs_delrows!(m, row, row)
LQOI.lqs_delrows!(instance::XpressOptimizer, rowbeg, rowend) = XPR.del_constrs!(instance.inner, cintvec(collect(rowbeg:rowend))) 

# LQOI.lqs_chgctype!(m, colvec, typevec)
# TODO fix types
LQOI.lqs_chgctype!(instance::XpressOptimizer, colvec, typevec) = XPR.chgcoltypes!(instance.inner, colvec, typevec)

# LQOI.lqs_chgsense!(m, rowvec, sensevec)
# TODO fix types
LQOI.lqs_chgsense!(instance::XpressOptimizer, rowvec, sensevec) = XPR.set_rowtype!(instance.inner, rowvec, sensevec)

const VAR_TYPE_MAP = Dict{Symbol,Cchar}(
    :CONTINUOUS => Cchar('C'),
    :INTEGER => Cchar('I'),
    :BINARY => Cchar('B')
)
LQOI.lqs_vartype_map(m::XpressOptimizer) = VAR_TYPE_MAP

# LQOI.lqs_addsos(m, colvec, valvec, typ)
LQOI.lqs_addsos!(instance::XpressOptimizer, colvec, valvec, typ) = XPR.add_sos!(instance.inner, typ, colvec, valvec)
# LQOI.lqs_delsos(m, idx, idx)
LQOI.lqs_delsos!(instance::XpressOptimizer, idx1, idx2) = XPR.del_sos!(instance.inner, cintvec(collect(idx1:idx2)))

const SOS_TYPE_MAP = Dict{Symbol,Symbol}(
    :SOS1 => :SOS1,#Cchar('1'),
    :SOS2 => :SOS2#Cchar('2')
)
LQOI.lqs_sertype_map(m::XpressOptimizer) = SOS_TYPE_MAP

# LQOI.lqs_getsos(m, idx)
# TODO improve getting processes
function LQOI.lqs_getsos(instance::XpressOptimizer, idx)
    A, types = XPR.get_sos_matrix(instance.inner)
    line = A[idx,:] #sparse vec
    cols = line.nzind
    vals = line.nzval
    typ = types[idx] == Cchar('1') ? :SOS1 : :SOS2
    return cols, vals, typ
end

# LQOI.lqs_getnumqconstrs(m)
LQOI.lqs_getnumqconstrs(instance::XpressOptimizer) = XPR.num_qconstrs(instance.inner)

# LQOI.lqs_addqconstr(m, cols,coefs,rhs,sense, I,J,V)
LQOI.lqs_addqconstr!(instance::XpressOptimizer, cols,coefs,rhs,sense, I,J,V) = XPR.add_qconstr!(instance.inner, cols, coefs, I, J, V, sense, rhs)

# LQOI.lqs_chgrngval
LQOI.lqs_chgrngval!(instance::XpressOptimizer, rows, vals) = XPR.chg_rhsrange!(instance.inner, cintvec(rows), -vals)

const CTR_TYPE_MAP = Dict{Symbol,Cchar}(
    :RANGE => Cchar('R'),
    :LOWER => Cchar('L'),
    :UPPER => Cchar('U'),
    :EQUALITY => Cchar('E')
)
LQOI.lqs_ctrtype_map(m::XpressOptimizer) = CTR_TYPE_MAP

#=
    Objective
=#

# LQOI.lqs_copyquad(m, intvec,intvec, floatvec) #?
function LQOI.lqs_copyquad!(instance::XpressOptimizer, I, J, V)
    XPR.delq!(instance.inner)
    XPR.add_qpterms!(instance.inner, I, J, V)
    return nothing
end

# LQOI.lqs_chgobj(m, colvec,coefvec)
function LQOI.lqs_chgobj!(instance::XpressOptimizer, colvec, coefvec) 
    nvars = XPR.num_vars(instance.inner)
    obj = zeros(Float64, nvars)

    for i in eachindex(colvec)
        obj[colvec[i]] = coefvec[i]
    end

    XPR.set_obj!(instance.inner, obj)
    nothing
end

# LQOI.lqs_chgobjsen(m, symbol)
# TODO improve min max names
function LQOI.lqs_chgobjsen!(instance::XpressOptimizer, symbol)
    if symbol == :Min
        XPR.set_sense!(instance.inner, :minimize)
    else
        XPR.set_sense!(instance.inner, :maximize)
    end
end
    

# LQOI.lqs_getobj(m)
LQOI.lqs_getobj(instance::XpressOptimizer) = XPR.get_obj(instance.inner)

# lqs_getobjsen(m)
function LQOI.lqs_getobjsen(instance::XpressOptimizer)
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

# LQOI.lqs_getnumcols(m)
LQOI.lqs_getnumcols(instance::XpressOptimizer) = XPR.num_vars(instance.inner)

# LQOI.lqs_newcols!(m, int)
LQOI.lqs_newcols!(instance::XpressOptimizer, int) = XPR.add_cvars!(instance.inner, zeros(int))

# LQOI.lqs_delcols!(m, col, col)
LQOI.lqs_delcols!(instance::XpressOptimizer, col, col2) = XPR.del_vars!(instance.inner, col)

# LQOI.lqs_addmipstarts(m, colvec, valvec)
function LQOI.lqs_addmipstarts!(instance::XpressOptimizer, colvec, valvec) 
    x = zeros(XPR.num_vars(instance.inner))
    for i in eachindex(colvec)
        x[colvec[i]] = valvec[i]
    end
    XPR.loadbasis(instance.inner, x)
end
#=
    Solve
=#

# LQOI.lqs_mipopt!(m)
LQOI.lqs_mipopt!(instance::XpressOptimizer) = XPR.mipoptimize(instance.inner)

# LQOI.lqs_qpopt!(m)
LQOI.lqs_qpopt!(instance::XpressOptimizer) = LQOI.lqs_lpopt!(instance)

# LQOI.lqs_lpopt!(m)
LQOI.lqs_lpopt!(instance::XpressOptimizer) = XPR.lpoptimize(instance.inner)

# LQOI.lqs_terminationstatus(m)
function LQOI.lqs_terminationstatus(instance::XpressOptimizer)
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
            return MOI.InvalidInstance
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

function LQOI.lqs_primalstatus(instance::XpressOptimizer)
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
function LQOI.lqs_dualstatus(instance::XpressOptimizer) 
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


# LQOI.lqs_getx!(m, place)
LQOI.lqs_getx!(instance::XpressOptimizer, place) = XPR.get_solution!(instance.inner, place)

# LQOI.lqs_getax!(m, place)
function LQOI.lqs_getax!(instance::XpressOptimizer, place)
    XPR.get_slack_lin!(instance.inner, place)
    rhs = XPR.get_rhs(instance.inner)
    for i in eachindex(place)
        place[i] = -place[i]+rhs[i]
    end
    nothing
end
# LQOI.lqs_getdj!(m, place)
LQOI.lqs_getdj!(instance::XpressOptimizer, place) = XPR.get_reducedcost!(instance.inner, place)

# LQOI.lqs_getpi!(m, place)
LQOI.lqs_getpi!(instance::XpressOptimizer, place) = XPR.get_dual_lin!(instance.inner, place)

# LQOI.lqs_getobjval(m)
LQOI.lqs_getobjval(instance::XpressOptimizer) = XPR.get_objval(instance.inner)

# LQOI.lqs_getbestobjval(m)
LQOI.lqs_getbestobjval(instance::XpressOptimizer) = XPR.get_mip_objval(instance.inner)

# LQOI.lqs_getmiprelgap(m)
function LQOI.lqs_getmiprelgap(instance::XpressOptimizer)
    L = XPR.get_mip_objval(instance.inner)
    U = XPR.get_bestbound(instance.inner)
    return abs(U-L)/U
end

# LQOI.lqs_getitcnt(m)
LQOI.lqs_getitcnt(instance::XpressOptimizer)  = XPR.get_simplex_iter_count(instance.inner)

# LQOI.lqs_getbaritcnt(m)
LQOI.lqs_getbaritcnt(instance::XpressOptimizer) = XPR.get_barrier_iter_count(instance.inner)

# LQOI.lqs_getnodecnt(m)
LQOI.lqs_getnodecnt(instance::XpressOptimizer) = XPR.get_node_count(instance.inner)

# LQOI.lqs_dualfarkas(m, place)
LQOI.lqs_dualfarkas!(instance::XpressOptimizer, place) = XPR.getdualray!(instance.inner, place)

# LQOI.lqs_getray(m, place)
LQOI.lqs_getray!(instance::XpressOptimizer, place) = XPR.getprimalray!(instance.inner, place)


MOI.free!(m::XpressOptimizer) = XPR.free_model(m.onner)

"""
    writeproblem(m: :MOI.AbstractOptimizer, filename::String)
Writes the current problem data to the given file.
Supported file types are solver-dependent.
"""
writeproblem(instance::XpressOptimizer, filename::String, flags::String="") = XPR.write_model(instance.inner, filename, flags)

end # module