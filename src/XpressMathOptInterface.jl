module XpressMathOptInterface

import Base.show, Base.copy

# Standard LP interface
# importall MathProgBase.SolverInterface

using Xpress
const XPR = Xpress
using MathOptInterface
const MOI = MathOptInterface
using LinQuadOptInterface
const LQOI = LinQuadOptInterface


export XpressSolver
struct XpressSolver <: LQOI.LinQuadSolver
    options
end
function XpressSolver(;kwargs...)
    return XpressSolver(kwargs)
end

import Xpress.Model
mutable struct XpressSolverInstance <: LQOI.LinQuadSolverInstance

    LQOI.@LinQuadSolverInstanceBase
    
end

function MOI.SolverInstance(s::XpressSolver)

    env = nothing
    m = XpressSolverInstance(
        (LQOI.@LinQuadSolverInstanceBaseInit)...,
    )
    for (name,value) in s.options
        XPR.setparam!(m.inner, XPR.XPRS_CONTROLS_DICT[name], value)
    end
    # csi.inner.mipstart_effort = s.mipstart_effortlevel
    # if s.logfile != ""
    #     LQOI.lqs_setlogfile!(env, s.logfile)
    # end
    return m
end

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
LQOI.lqs_setparam!(m::XpressSolverInstance, name, val) = XPR.setparam!(m.inner, XPR.XPRS_CONTROLS_DICT[name], val)

# LQOI.lqs_setlogfile!(env, path)
# TODO fix this one
LQOI.lqs_setlogfile!(m::XpressSolverInstance, path) = XPR.setlogfile(m.inner, path::String)

# LQOI.lqs_getprobtype(m)
# TODO - consider removing, apparently useless

#=
    Constraints
=#

cintvec(v::Vector) = convert(Vector{Int32}, v)
cdoublevec(v::Vector) = convert(Vector{Float64}, v)

# LQOI.lqs_chgbds!(m, colvec, valvec, sensevec)
LQOI.lqs_chgbds!(m::XPR.Model, colvec, valvec, sensevec) = XPR.chgbounds!(m, cintvec(colvec), sensevec, valvec)

# LQOI.lqs_getlb(m, col)
LQOI.lqs_getlb(m::XPR.Model, col) = XPR.get_lb(m, col, col)[1]
# LQOI.lqs_getub(m, col)
LQOI.lqs_getub(m::XPR.Model, col) = XPR.get_ub(m, col, col)[1]

# LQOI.lqs_getnumrows(m)
LQOI.lqs_getnumrows(m::XPR.Model) = XPR.num_linconstrs(m)

# LQOI.lqs_addrows!(m, rowvec, colvec, coefvec, sensevec, rhsvec)
LQOI.lqs_addrows!(m::XPR.Model, rowvec, colvec, coefvec, sensevec, rhsvec) = XPR.add_constrs!(m::XPR.Model, rowvec, colvec, coefvec, sensevec, rhsvec)

# LQOI.lqs_getrhs(m, rowvec)
LQOI.lqs_getrhs(m::XPR.Model, row) = XPR.get_rhs(m, row, row)[1]  

# colvec, coef = LQOI.lqs_getrows(m, rowvec)
# TODO improve
function LQOI.lqs_getrows(m::XPR.Model, idx)
    A = XPR.get_rows(m, idx, idx)'
    return A.rowval-1, A.nzval
end

# LQOI.lqs_getcoef(m, row, col) #??
# TODO improve
function LQOI.lqs_getcoef(m::XPR.Model, row, col) #??
    A = XPR.get_rows(m, row, row)'
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
function LQOI.lqs_chgcoef!(m::XPR.Model, row, col, coef) 
    if row == 0
        XPR.set_objcoeffs!(m, Int32[col], Float64[coef])
    elseif col == 0
        XPR.set_rhs!(m, Int32[row], Float64[coef])
    else
        XPR.chg_coeffs!(m, row, col, coef)
    end
end

# LQOI.lqs_delrows!(m, row, row)
LQOI.lqs_delrows!(m::XPR.Model, rowbeg, rowend) = XPR.del_constrs!(m, cintvec(collect(rowbeg:rowend))) 

# LQOI.lqs_chgctype!(m, colvec, typevec)
# TODO fix types
LQOI.lqs_chgctype!(m::XPR.Model, colvec, typevec) = XPR.chgcoltypes!(m, colvec, typevec)

# LQOI.lqs_chgsense!(m, rowvec, sensevec)
# TODO fix types
LQOI.lqs_chgsense!(m::XPR.Model, rowvec, sensevec) = XPR.set_rowtype!(m, rowvec, sensevec)

const VAR_TYPE_MAP = Dict{Symbol,Cchar}(
    :CONTINUOUS => Cchar('C'),
    :INTEGER => Cchar('I'),
    :BINARY => Cchar('B')
)
LQOI.lqs_vartype_map(m::XpressSolverInstance) = VAR_TYPE_MAP

# LQOI.lqs_addsos(m, colvec, valvec, typ)
LQOI.lqs_addsos!(m::XPR.Model, colvec, valvec, typ) = XPR.add_sos!(m, typ, colvec, valvec)
# LQOI.lqs_delsos(m, idx, idx)
LQOI.lqs_delsos!(m::XPR.Model, idx1, idx2) = XPR.del_sos!(m, cintvec(collect(idx1:idx2)))

const SOS_TYPE_MAP = Dict{Symbol,Symbol}(
    :SOS1 => :SOS1,#Cchar('1'),
    :SOS2 => :SOS2#Cchar('2')
)
LQOI.lqs_sertype_map(m::XpressSolverInstance) = SOS_TYPE_MAP

# LQOI.lqs_getsos(m, idx)
# TODO improve getting processes
function LQOI.lqs_getsos(m::XPR.Model, idx)
    A, types = XPR.get_sos_matrix(m::XPR.Model)
    line = A[idx,:] #sparse vec
    cols = line.nzind
    vals = line.nzval
    typ = types[idx] == Cchar('1') ? :SOS1 : :SOS2
    return cols, vals, typ
end

# LQOI.lqs_getnumqconstrs(m)
LQOI.lqs_getnumqconstrs(m::XPR.Model) = XPR.num_qconstrs(m)

# LQOI.lqs_addqconstr(m, cols,coefs,rhs,sense, I,J,V)
LQOI.lqs_addqconstr!(m::XPR.Model, cols,coefs,rhs,sense, I,J,V) = XPR.add_qconstr!(m, cols, coefs, I, J, V, sense, rhs)

# LQOI.lqs_chgrngval
LQOI.lqs_chgrngval!(m::XPR.Model, rows, vals) = XPR.chg_rhsrange!(m, cintvec(rows), -vals)

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
function LQOI.lqs_copyquad!(m::XPR.Model, I, J, V)
    XPR.delq!(m)
    XPR.add_qpterms!(m, I, J, V)
    return nothing
end

# LQOI.lqs_chgobj(m, colvec,coefvec)
function LQOI.lqs_chgobj!(m::XPR.Model, colvec, coefvec) 
    nvars = XPR.num_vars(m)
    obj = zeros(Float64, nvars)

    for i in eachindex(colvec)
        obj[colvec[i]] = coefvec[i]
    end

    XPR.set_obj!(m, obj)
    nothing
end

# LQOI.lqs_chgobjsen(m, symbol)
# TODO improve min max names
function LQOI.lqs_chgobjsen!(m::XPR.Model, symbol)
    if symbol == :Min
        XPR.set_sense!(m, :minimize)
    else
        XPR.set_sense!(m, :maximize)
    end
end
    

# LQOI.lqs_getobj(m)
LQOI.lqs_getobj(m::XPR.Model) = XPR.get_obj(m)

# lqs_getobjsen(m)
function LQOI.lqs_getobjsen(m::XPR.Model)
    s = XPR.model_sense(m)
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
LQOI.lqs_getnumcols(m::XPR.Model) = XPR.num_vars(m)

# LQOI.lqs_newcols!(m, int)
LQOI.lqs_newcols!(m::XPR.Model, int) = XPR.add_cvars!(m, zeros(int))

# LQOI.lqs_delcols!(m, col, col)
LQOI.lqs_delcols!(m::XPR.Model, col, col2) = XPR.del_vars!(m, col)

# LQOI.lqs_addmipstarts(m, colvec, valvec)
function LQOI.lqs_addmipstarts!(m::XPR.Model, colvec, valvec) 
    x = zeros(XPR.num_vars(m))
    for i in eachindex(colvec)
        x[colvec[i]] = valvec[i]
    end
    XPR.loadbasis(m, x)
end
#=
    Solve
=#

# LQOI.lqs_mipopt!(m)
LQOI.lqs_mipopt!(m::XPR.Model) = XPR.mipoptimize(m)

# LQOI.lqs_qpopt!(m)
LQOI.lqs_qpopt!(m::XPR.Model) = LQOI.lqs_lpopt!(m)

# LQOI.lqs_lpopt!(m)
LQOI.lqs_lpopt!(m::XPR.Model) = XPR.lpoptimize(m)

# LQOI.lqs_terminationstatus(m)
function LQOI.lqs_terminationstatus(model::XpressSolverInstance)
    m = model.inner 
    stat_lp = XPR.get_lp_status2(m)
    if XPR.is_mip(m)
        stat_mip = XPR.get_mip_status2(m)
        if stat_mip == XPR.MIP_NotLoaded
            return MOI.OtherError
        elseif stat_mip == XPR.MIP_LPNotOptimal
            # MIP search incomplete but there is no linear sol
            # return MOI.OtherError
            return MOI.InfeasibleOrUnbounded
        elseif stat_mip == XPR.MIP_NoSolFound
            # MIP search incomplete but there is no integer sol
            other = xprsmoi_stopstatus(m)
            if other == MOI.OtherError
                return MOI.SlowProgress#OtherLimit
            else 
                return other
            end

        elseif stat_mip == XPR.MIP_Solution
            # MIP search incomplete but there is a solution
            other = xprsmoi_stopstatus(m)
            if other == MOI.OtherError
                return MOI.OtherLimit
            else 
                return other
            end

        elseif stat_mip == XPR.MIP_Infeasible
            if XPR.hasdualray(m)
                return MOI.Success
            else
                return MOI.InfeasibleNoResult
            end
        elseif stat_mip == XPR.MIP_Optimal
            return MOI.Success
        elseif stat_mip == XPR.MIP_Unbounded
            if XPR.hasprimalray(m)
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
            if XPR.hasdualray(m)
                return MOI.Success
            else
                return MOI.InfeasibleNoResult
            end
        elseif stat_lp == XPR.LP_CutOff
            return MOI.ObjectiveLimit
        elseif stat_lp == XPR.LP_Unfinished
            return xprsmoi_stopstatus(m)
        elseif stat_lp == XPR.LP_Unbounded
            if XPR.hasprimalray(m)
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

function xprsmoi_stopstatus(m::XPR.Model)
    ss = XPR.get_stopstatus(m)
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

function LQOI.lqs_primalstatus(model::XpressSolverInstance)
    m = model.inner
    if XPR.is_mip(m)
        stat_mip = XPR.get_mip_status2(m)
        if stat_mip in [XPR.MIP_Solution, XPR.MIP_Optimal]
            return MOI.FeasiblePoint
        elseif XPR.MIP_Infeasible && XPR.hasdualray(m)
            return MOI.InfeasibilityCertificate
        elseif XPR.MIP_Unbounded && XPR.hasprimalray(m)
            return MOI.InfeasibilityCertificate
        elseif stat_mip in [XPR.MIP_LPOptimal, XPR.MIP_NoSolFound]
            return MOI.InfeasiblePoint
        end
        return MOI.UnknownResultStatus
    else
        stat_lp = XPR.get_lp_status2(m)
        if stat_lp == XPR.LP_Optimal
            return MOI.FeasiblePoint
        elseif stat_lp == XPR.LP_Unbounded && XPR.hasprimalray(m)
            return MOI.InfeasibilityCertificate
        # elseif stat_lp == LP_Infeasible
        #     return MOI.InfeasiblePoint - xpress wont return
        # elseif cutoff//cutoffindual ???
        else
            return MOI.UnknownResultStatus
        end
    end
end
function LQOI.lqs_dualstatus(model::XpressSolverInstance)
    m = model.inner    
    if XPR.is_mip(m)
        return MOI.UnknownResultStatus
    else
        stat_lp = XPR.get_lp_status2(m)
        if stat_lp == XPR.LP_Optimal
            return MOI.FeasiblePoint
        elseif stat_lp == XPR.LP_Infeasible && XPR.hasdualray(m)
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
LQOI.lqs_getx!(m::XPR.Model, place) = XPR.get_solution!(m, place)

# LQOI.lqs_getax!(m, place)
function LQOI.lqs_getax!(m::XPR.Model, place)
    XPR.get_slack_lin!(m, place)
    rhs = XPR.get_rhs(m)
    for i in eachindex(place)
        place[i] = -place[i]+rhs[i]
    end
    nothing
end
# LQOI.lqs_getdj!(m, place)
LQOI.lqs_getdj!(m::XPR.Model, place) = XPR.get_reducedcost!(m, place)

# LQOI.lqs_getpi!(m, place)
LQOI.lqs_getpi!(m::XPR.Model, place) = XPR.get_dual_lin!(m, place)

# LQOI.lqs_getobjval(m)
LQOI.lqs_getobjval(m::XPR.Model) = XPR.get_objval(m)

# LQOI.lqs_getbestobjval(m)
LQOI.lqs_getbestobjval(m::XPR.Model) = XPR.get_mip_objval(m)

# LQOI.lqs_getmiprelgap(m)
function LQOI.lqs_getmiprelgap(m::XPR.Model)
    L = XPR.get_mip_objval(m)
    U = XPR.get_bestbound(m)
    return abs(U-L)/U
end

# LQOI.lqs_getitcnt(m)
LQOI.lqs_getitcnt(m::XPR.Model)  = XPR.get_simplex_iter_count(m)

# LQOI.lqs_getbaritcnt(m)
LQOI.lqs_getbaritcnt(m::XPR.Model) = XPR.get_barrier_iter_count(m)

# LQOI.lqs_getnodecnt(m)
LQOI.lqs_getnodecnt(m::XPR.Model) = XPR.get_node_count(m)

# LQOI.lqs_dualfarkas(m, place)
LQOI.lqs_dualfarkas!(m::XPR.Model, place) = XPR.getdualray!(m, place)

# LQOI.lqs_getray(m, place)
LQOI.lqs_getray!(m::XPR.Model, place) = XPR.getprimalray!(m, place)


MOI.free!(m::XpressSolverInstance) = XPR.free_model(m.onner)

"""
    writeproblem(m::AbstractSolverInstance, filename::String)
Writes the current problem data to the given file.
Supported file types are solver-dependent.
"""
MOI.writeproblem(m::XpressSolverInstance, filename::String, flags::String="") = XPR.write_model(m.inner, filename, flags)

end # module