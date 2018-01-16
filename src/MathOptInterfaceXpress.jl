module MathOptInterfaceXpress

import Base.show, Base.copy

# Standard LP interface
# importall MathProgBase.SolverInterface

using Xpress
const XPR = Xpress
using MathOptInterface
const MOI = MathOptInterface
using LinQuadOptInterface
const LQOI = LinQuadOptInterface

import Xpress.Model
Xpress.Model(::Void) = Xpress.Model()
mutable struct XpressSolverInstance <: LQOI.LinQuadSolverInstance
    LQOI.@LinQuadSolverInstanceBase
end
function XpressSolverInstance(;kwargs...)

    env = nothing
    instance = XpressSolverInstance(
        (LQOI.@LinQuadSolverInstanceBaseInit)...,
    )
    for (name,value) in kwargs
        XPR.setparam!(instance.inner, XPR.XPRS_CONTROLS_DICT[name], value)
    end
    # csi.inner.mipstart_effort = s.mipstart_effortlevel
    # if s.logfile != ""
    #     LQOI.lqs_setlogfile!(env, s.logfile)
    # end
    return instance
end

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
    (LQOI.VecVar, MOI.SOS1),
    (LQOI.VecVar, MOI.SOS2),
    (LQOI.VecVar, MOI.Nonnegatives),
    (LQOI.VecVar, MOI.Nonpositives),
    (LQOI.VecVar, MOI.Zeros),
    (LQOI.VecLin, MOI.Nonnegatives),
    (LQOI.VecLin, MOI.Nonpositives),
    (LQOI.VecLin, MOI.Zeros)
]

const SUPPORTED_OBJECTIVES = [
    LQOI.Linear,
    LQOI.Quad
]

LQOI.lqs_supported_constraints(s::XpressSolverInstance) = SUPPORTED_CONSTRAINTS
LQOI.lqs_supported_objectives(s::XpressSolverInstance) = SUPPORTED_OBJECTIVES

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
LQOI.lqs_setparam!(instance::XpressSolverInstance, name, val) = XPR.setparam!(instance.inner, XPR.XPRS_CONTROLS_DICT[name], val)

# LQOI.lqs_setlogfile!(env, path)
# TODO fix this one
LQOI.lqs_setlogfile!(instance::XpressSolverInstance, path) = XPR.setlogfile(instance.inner, path::String)

# LQOI.lqs_getprobtype(m)
# TODO - consider removing, apparently useless

#=
    Constraints
=#

cintvec(v::Vector) = convert(Vector{Int32}, v)
cdoublevec(v::Vector) = convert(Vector{Float64}, v)

# LQOI.lqs_chgbds!(m, colvec, valvec, sensevec)
LQOI.lqs_chgbds!(instance::XpressSolverInstance, colvec, valvec, sensevec) = XPR.chgbounds!(instance.inner, cintvec(colvec), sensevec, valvec)

# LQOI.lqs_getlb(m, col)
LQOI.lqs_getlb(instance::XpressSolverInstance, col) = XPR.get_lb(instance.inner, col, col)[1]
# LQOI.lqs_getub(m, col)
LQOI.lqs_getub(instance::XpressSolverInstance, col) = XPR.get_ub(instance.inner, col, col)[1]

# LQOI.lqs_getnumrows(m)
LQOI.lqs_getnumrows(instance::XpressSolverInstance) = XPR.num_linconstrs(instance.inner)

# LQOI.lqs_addrows!(m, rowvec, colvec, coefvec, sensevec, rhsvec)
LQOI.lqs_addrows!(instance::XpressSolverInstance, rowvec, colvec, coefvec, sensevec, rhsvec) = XPR.add_constrs!(instance.inner, rowvec, colvec, coefvec, sensevec, rhsvec)

# LQOI.lqs_getrhs(m, rowvec)
LQOI.lqs_getrhs(instance::XpressSolverInstance, row) = XPR.get_rhs(instance.inner, row, row)[1]  

# colvec, coef = LQOI.lqs_getrows(m, rowvec)
# TODO improve
function LQOI.lqs_getrows(instance::XpressSolverInstance, idx)
    A = XPR.get_rows(instance.inner, idx, idx)'
    return A.rowval-1, A.nzval
end

# LQOI.lqs_getcoef(m, row, col) #??
# TODO improve
function LQOI.lqs_getcoef(instance::XpressSolverInstance, row, col) #??
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
function LQOI.lqs_chgcoef!(instance::XpressSolverInstance, row, col, coef) 
    if row == 0
        XPR.set_objcoeffs!(instance.inner, Int32[col], Float64[coef])
    elseif col == 0
        XPR.set_rhs!(instance.inner, Int32[row], Float64[coef])
    else
        XPR.chg_coeffs!(instance.inner, row, col, coef)
    end
end

# LQOI.lqs_delrows!(m, row, row)
LQOI.lqs_delrows!(instance::XpressSolverInstance, rowbeg, rowend) = XPR.del_constrs!(instance.inner, cintvec(collect(rowbeg:rowend))) 

# LQOI.lqs_chgctype!(m, colvec, typevec)
# TODO fix types
LQOI.lqs_chgctype!(instance::XpressSolverInstance, colvec, typevec) = XPR.chgcoltypes!(instance.inner, colvec, typevec)

# LQOI.lqs_chgsense!(m, rowvec, sensevec)
# TODO fix types
LQOI.lqs_chgsense!(instance::XpressSolverInstance, rowvec, sensevec) = XPR.set_rowtype!(instance.inner, rowvec, sensevec)

const VAR_TYPE_MAP = Dict{Symbol,Cchar}(
    :CONTINUOUS => Cchar('C'),
    :INTEGER => Cchar('I'),
    :BINARY => Cchar('B')
)
LQOI.lqs_vartype_map(m::XpressSolverInstance) = VAR_TYPE_MAP

# LQOI.lqs_addsos(m, colvec, valvec, typ)
LQOI.lqs_addsos!(instance::XpressSolverInstance, colvec, valvec, typ) = XPR.add_sos!(instance.inner, typ, colvec, valvec)
# LQOI.lqs_delsos(m, idx, idx)
LQOI.lqs_delsos!(instance::XpressSolverInstance, idx1, idx2) = XPR.del_sos!(instance.inner, cintvec(collect(idx1:idx2)))

const SOS_TYPE_MAP = Dict{Symbol,Symbol}(
    :SOS1 => :SOS1,#Cchar('1'),
    :SOS2 => :SOS2#Cchar('2')
)
LQOI.lqs_sertype_map(m::XpressSolverInstance) = SOS_TYPE_MAP

# LQOI.lqs_getsos(m, idx)
# TODO improve getting processes
function LQOI.lqs_getsos(instance::XpressSolverInstance, idx)
    A, types = XPR.get_sos_matrix(instance.inner)
    line = A[idx,:] #sparse vec
    cols = line.nzind
    vals = line.nzval
    typ = types[idx] == Cchar('1') ? :SOS1 : :SOS2
    return cols, vals, typ
end

# LQOI.lqs_getnumqconstrs(m)
LQOI.lqs_getnumqconstrs(instance::XpressSolverInstance) = XPR.num_qconstrs(instance.inner)

# LQOI.lqs_addqconstr(m, cols,coefs,rhs,sense, I,J,V)
LQOI.lqs_addqconstr!(instance::XpressSolverInstance, cols,coefs,rhs,sense, I,J,V) = XPR.add_qconstr!(instance.inner, cols, coefs, I, J, V, sense, rhs)

# LQOI.lqs_chgrngval
LQOI.lqs_chgrngval!(instance::XpressSolverInstance, rows, vals) = XPR.chg_rhsrange!(instance.inner, cintvec(rows), -vals)

const CTR_TYPE_MAP = Dict{Symbol,Cchar}(
    :RANGE => Cchar('R'),
    :LOWER => Cchar('L'),
    :UPPER => Cchar('U'),
    :EQUALITY => Cchar('E')
)
LQOI.lqs_ctrtype_map(m::XpressSolverInstance) = CTR_TYPE_MAP

#=
    Objective
=#

# LQOI.lqs_copyquad(m, intvec,intvec, floatvec) #?
function LQOI.lqs_copyquad!(instance::XpressSolverInstance, I, J, V)
    XPR.delq!(instance.inner)
    XPR.add_qpterms!(instance.inner, I, J, V)
    return nothing
end

# LQOI.lqs_chgobj(m, colvec,coefvec)
function LQOI.lqs_chgobj!(instance::XpressSolverInstance, colvec, coefvec) 
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
function LQOI.lqs_chgobjsen!(instance::XpressSolverInstance, symbol)
    if symbol == :Min
        XPR.set_sense!(instance.inner, :minimize)
    else
        XPR.set_sense!(instance.inner, :maximize)
    end
end
    

# LQOI.lqs_getobj(m)
LQOI.lqs_getobj(instance::XpressSolverInstance) = XPR.get_obj(instance.inner)

# lqs_getobjsen(m)
function LQOI.lqs_getobjsen(instance::XpressSolverInstance)
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
LQOI.lqs_getnumcols(instance::XpressSolverInstance) = XPR.num_vars(instance.inner)

# LQOI.lqs_newcols!(m, int)
LQOI.lqs_newcols!(instance::XpressSolverInstance, int) = XPR.add_cvars!(instance.inner, zeros(int))

# LQOI.lqs_delcols!(m, col, col)
LQOI.lqs_delcols!(instance::XpressSolverInstance, col, col2) = XPR.del_vars!(instance.inner, col)

# LQOI.lqs_addmipstarts(m, colvec, valvec)
function LQOI.lqs_addmipstarts!(instance::XpressSolverInstance, colvec, valvec) 
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
LQOI.lqs_mipopt!(instance::XpressSolverInstance) = XPR.mipoptimize(instance.inner)

# LQOI.lqs_qpopt!(m)
LQOI.lqs_qpopt!(instance::XpressSolverInstance) = LQOI.lqs_lpopt!(instance)

# LQOI.lqs_lpopt!(m)
LQOI.lqs_lpopt!(instance::XpressSolverInstance) = XPR.lpoptimize(instance.inner)

# LQOI.lqs_terminationstatus(m)
function LQOI.lqs_terminationstatus(instance::XpressSolverInstance)
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

function xprsmoi_stopstatus(instance::XpressSolverInstance)
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

function LQOI.lqs_primalstatus(instance::XpressSolverInstance)
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
function LQOI.lqs_dualstatus(instance::XpressSolverInstance) 
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
LQOI.lqs_getx!(instance::XpressSolverInstance, place) = XPR.get_solution!(instance.inner, place)

# LQOI.lqs_getax!(m, place)
function LQOI.lqs_getax!(instance::XpressSolverInstance, place)
    XPR.get_slack_lin!(instance.inner, place)
    rhs = XPR.get_rhs(instance.inner)
    for i in eachindex(place)
        place[i] = -place[i]+rhs[i]
    end
    nothing
end
# LQOI.lqs_getdj!(m, place)
LQOI.lqs_getdj!(instance::XpressSolverInstance, place) = XPR.get_reducedcost!(instance.inner, place)

# LQOI.lqs_getpi!(m, place)
LQOI.lqs_getpi!(instance::XpressSolverInstance, place) = XPR.get_dual_lin!(instance.inner, place)

# LQOI.lqs_getobjval(m)
LQOI.lqs_getobjval(instance::XpressSolverInstance) = XPR.get_objval(instance.inner)

# LQOI.lqs_getbestobjval(m)
LQOI.lqs_getbestobjval(instance::XpressSolverInstance) = XPR.get_mip_objval(instance.inner)

# LQOI.lqs_getmiprelgap(m)
function LQOI.lqs_getmiprelgap(instance::XpressSolverInstance)
    L = XPR.get_mip_objval(instance.inner)
    U = XPR.get_bestbound(instance.inner)
    return abs(U-L)/U
end

# LQOI.lqs_getitcnt(m)
LQOI.lqs_getitcnt(instance::XpressSolverInstance)  = XPR.get_simplex_iter_count(instance.inner)

# LQOI.lqs_getbaritcnt(m)
LQOI.lqs_getbaritcnt(instance::XpressSolverInstance) = XPR.get_barrier_iter_count(instance.inner)

# LQOI.lqs_getnodecnt(m)
LQOI.lqs_getnodecnt(instance::XpressSolverInstance) = XPR.get_node_count(instance.inner)

# LQOI.lqs_dualfarkas(m, place)
LQOI.lqs_dualfarkas!(instance::XpressSolverInstance, place) = XPR.getdualray!(instance.inner, place)

# LQOI.lqs_getray(m, place)
LQOI.lqs_getray!(instance::XpressSolverInstance, place) = XPR.getprimalray!(instance.inner, place)


MOI.free!(m::XpressSolverInstance) = XPR.free_model(m.onner)

"""
    writeproblem(m::AbstractSolverInstance, filename::String)
Writes the current problem data to the given file.
Supported file types are solver-dependent.
"""
writeproblem(instance::XpressSolverInstance, filename::String, flags::String="") = XPR.write_model(instance.inner, filename, flags)

end # module