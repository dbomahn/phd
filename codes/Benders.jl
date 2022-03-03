using JuMP,CPLEX,LinearAlgebra,CPLEX
using MathProgBase
const MPB = MathProgBase

https://github.com/matbesancon/SimpleBenders.jl/blob/master/test/runtests.jl
https://matbesancon.xyz/post/2019-05-08-simple-benders/
https://co-at-work.zib.de/slides/Donnerstag_24.9/Benders_decomposition-Fundamentals.pdf

mutable struct CallModel
    lpfile::String; m::Int; n::Int; C::Array{}; B::Array{}; RHS::Dict{}; signs::Array{}; vub::Array{}
    function CallModel(lpfile::String)
        lpmodel=buildlp([-1,0],[2 1],'<',1.5, CplexSolver(CPX_PARAM_SCRIND=0))
        # lpmodel = CPLEX.CplexMathProgModel();
        MPB.loadproblem!(lpmodel,lpfile)
        Bmtx = MPB.getconstrmatrix(lpmodel);
        B = Bmtx[3:end,:]; C = Bmtx[1:2,:]
        # cut = find(i-> varub[i]==1 &&varub[i+1]!=1, 1:length(varub))[end]
        # vub = varub[1:cut]; B = Bmtx[3:end,1:cut]; C = Bmtx[1:2,1:cut]
        m,n=size(B)
        vub = MPB.getvarUB(lpmodel)
        lb = MPB.getconstrLB(lpmodel)[3:end]
        ub = MPB.getconstrUB(lpmodel)[3:end]
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
mt = CallModel("F:/scnd/test01S2.lp")
# mt = CallModel("/home/ak121396/Desktop/instances/SCND/test01S2.lp")

struct SubProblemData
    b::Vector{Float64}
    D::Matrix{Float64}
    A::Matrix{Float64}
    c::Vector{Float64}
end

struct DualSubProblem
    data::SubProblemData
    α::Vector{VariableRef}
    m::Model
end

function DualSubProblem(d::SubProblemData, m::Model)
    α = @variable(m, α[i = 1:size(d.A, 1)] >= 0)
    @constraint(m, dot(d.A, α) .<= d.c)
    return DualSubProblem(d, α, m)
end

function JuMP.optimize!(sp::DualSubProblem, yh)
    obj = sp.data.b .- sp.data.D * yh
    @objective(sp.m, Max, dot(obj, sp.α))
    optimize!(sp.m)
    st = termination_status(sp.m)
    if st == MOI.OPTIMAL
        α = JuMP.value.(sp.α)
        return (:OptimalityCut, α)
    elseif st == MOI.DUAL_INFEASIBLE
        return (:FeasibilityCut, α)
    else
        error("DualSubProblem error: status $status")
    end
end

function benders_optimize!(m::Model, y::Vector{VariableRef}, sd::SubProblemData, f::Union{Function,Type}; eta_bound::Real = -1000.0) #sp_optimizer,
    subproblem = Model(CPLEX.Optimizer)
    MOI.set(subproblem, MOI.RawParameter("CPX_PARAM_SCRIND"), false );
    MOI.set(subproblem, MOI.RawParameter("CPX_PARAM_THREADS"),1  );
    dsp = DualSubProblem(sd, subproblem)
    @variable(m, η >= eta_bound)
    @objective(m, Min, f(y) + η)
    optimize!(m)
    st = MOI.get(m, MOI.TerminationStatus())
    # restricted master has a solution or is unbounded
    nopt_cons, nfeas_cons = (0, 0)
    @info "Initial status $st"
    cuts = Tuple{Symbol, Vector{Float64}}[]
    while (st == MOI.DUAL_INFEASIBLE) || (st == MOI.OPTIMAL)
        optimize!(m)
        st = MOI.get(m, MOI.TerminationStatus())
        ŷ = JuMP.value.(y)
        η0 = JuMP.value(η)
        (res, α) = optimize!(dsp, ŷ)
        if res == :OptimalityCut
            @info "Optimality cut found"
            if η0 ≥ dot(α, (dsp.data.b - dsp.data.D * ŷ))
                break
            else
                nopt_cons += 1
                @constraint(m, η ≥ dot(α, (dsp.data.b - dsp.data.D * y)))
            end
        else
            @info "Feasibility cut found"
            nfeas_cons += 1
            @constraint(m, 0 ≥ dot(α, (dsp.data.b - dsp.data.D * y)))
        end
        push!(cuts, (res, α))
    end
    return (m, y, cuts, nopt_cons, nfeas_cons)
end
c1b = [mt.C[1,i] for i in bvar]; c1r = [mt.C[1,i] for i in rvar];
c2b = [mt.C[2,i] for i in bvar]; c2r = [mt.C[2,i] for i in rvar];
c1 = [c1b;c1r]; c2 = [c2b;c2r]
hcat(c1,c2)'

mt.B[:,cut1+1:cut2-1]


bvar
function SetData(mt)
    bvar = findall(i->i==1,mt.vub); rvar = findall(i->i!=1,mt.vub);

    c1r = [mt.C[1,i] for i in rvar];
    c2r = [mt.C[2,i] for i in rvar];
    c = c1r; #c2r

    cut1 = findall(i-> mt.vub[i]==1 && mt.vub[i+1]!=1, 1:length(mt.vub)-1)[end]
    cut2 = findall(i-> mt.vub[i]!=1 && mt.vub[i+1]==1, 1:length(mt.vub)-1)[end]+1
    A = hcat(mt.B[:,1:cut1],mt.B[:,cut2:end])
    D = zeros(2, 1) .+ [1, 3]
    b = [3, 4]
    return SubProblemData(b, D, A, c)
end

data =  test_data()
# objective function on y
c1b = [mt.C[1,i] for i in bvar]; c2b = [mt.C[2,i] for i in bvar]; #c1 = [c1b;c1r]; c2 = [c2b;c2r]
f(v) = dot(c1b,v[1:length(c1b)])
# initialize the problem
m = Model(CPLEX.Optimizer)
@variable(m, y[j=1:1] >= 0)
# solve and voilà
(m, y, cuts, nopt_cons, nfeas_cons) = benders_optimize!(m, y, data, f)
