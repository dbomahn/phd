using DataStructures,DataFrames,DelimitedFiles,JuMP,CPLEX,LinearAlgebra,StatsBase,MathProgBase,MathOptInterface
#,CSV,CPUTime,JLD2,
const MPB = MathProgBase;
mutable struct Data
    lpfile::String; m::Int; n::Int; C::Array{}; B::Array{}; RHS::Dict{}; signs::Array{}; vub::Array{}
    function Data(lpfile::String)
        lpmodel=buildlp([-1,0],[2 1],'<',1.5, CplexSolver(CPX_PARAM_SCRIND=0))
        # lpmodel = CPLEX.CplexMathProgModel();
        MPB.loadproblem!(lpmodel,lpfile)
        Bmtx = MPB.getconstrmatrix(lpmodel);B = Bmtx[3:end,:]
        C = Bmtx[1:2,:]
        m,n=size(B)
        vub = MPB.getvarUB(lpmodel)
        lb = MPB.getconstrLB(lpmodel)[3:end]; ub = MPB.getconstrUB(lpmodel)[3:end]
        RHS = Dict()
        for i=1:m
            if ub[i]==Inf
                RHS[i] = lb[i]
            else
                RHS[i] = ub[i]
            end
        end
        signs = []
        for i=1:m
            if ub[i] == Inf
                push!(signs,"l")
            elseif lb[i] == -Inf
                push!(signs,"u")
            else
                push!(signs, "s")
            end
        end
        new(lpfile,m,n,C,B,RHS,signs,vub)
    end
end

dt = Data("./lp/Test4S1.lp") #calling lp file with cplex.jl 0.6
st = findall(i->i!=1,dt.vub)[1]

m1 = Model(CPLEX.Optimizer)
MOI.set(m1, MOI.RawParameter("CPX_PARAM_SCRIND"), false )
MOI.set(m1, MOI.RawParameter("CPX_PARAM_THREADS"),1  )
@variables(m1, begin
    yu[1:st-1], Bin
    0<= xh[st:dt.n]
end)
@variable(m1, ex[1:dt.n] )
@variable(m1, 0 <= ep)
@variable(m1, 0<= ep′ )
@constraint(m1, [k=1:st-1], ex[k] == yu[k] )
@constraint(m1, [k=st:dt.n], ex[k] == xh[k] )
@constraint(m1, epcon1, dot(ex,dt.C[2,:]) <= ep)

for k=1:dt.m
    if dt.signs[k] == "l"
        @constraint(m1, dot(dt.B[k,1:st-1],yu)+dot(dt.B[k,st:end],xh) >= dt.RHS[k])
    elseif dt.signs[k] == "u"
        @constraint(m1, dot(dt.B[k,1:st-1],yu)+dot(dt.B[k,st:end],xh) <= dt.RHS[k])
    else
        @constraint(m1, dot(dt.B[k,1:st-1],yu)+dot(dt.B[k,st:end],xh) == dt.RHS[k])
    end
end
@objective(m1, Min, dot(ex,dt.C[1,:]) )


function dominated(y,P)
    st = false
    for k=1:length(P)
        if all( y .>= P[k])# && any(x > P[k])
            st = true; break
        else
            continue
        end
    end
    return st
end
function opt(ϵ,C)
    JuMP.fix(ep, ϵ; force = true);
    optimize!(m1)
    if termination_status(m1) == MOI.OPTIMAL
        return JuMP.value.(ex), [objective_value(m1),dot(JuMP.value.(ex),C[2,:])]
    else
        return nothing,nothing
    end
end

function epsilon(C)
    P = []; Y = []; ϵ = 2*10^(6); δ =2*10^(5); lb = 11^(6); fval = [0,ϵ]
    while fval[2] >= lb
        s,fval = opt(ϵ,C)
        println(fval)
        if s == nothing
            break
        end
        if dominated(fval,Y)==false
            push!(P,s); push!(Y,fval);
        end
        ϵ = ϵ-δ
    end
    return P,Y
end

ex,ey = epsilon(dt.C)
domFilter(ex,ey)
