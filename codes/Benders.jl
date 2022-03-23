using JuMP,CPLEX,LinearAlgebra,MathProgBase
const MPB = MathProgBase

"If we use weighted sum for BD, obj values must be normalised"
"We can use the weighted sum code by Gandibluex (Vopt?) and just solve sub problems with BD"

https://github.com/matbesancon/SimpleBenders.jl/blob/master/test/runtests.jl
https://matbesancon.xyz/post/2019-05-08-simple-benders/
https://co-at-work.zib.de/slides/Donnerstag_24.9/Benders_decomposition-Fundamentals.pdf

mutable struct SCNDModel
    lpfile::String; m::Int; n::Int; C::Array{}; B::Array{}; RHS::Array{}; signs::Array{}; vub::Array{}
    function SCNDModel(lpfile::String)
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

        RHS = []
        for i=1:m
            if ub[i]==Inf
                push!(RHS,lb[i])
            else
                push!(RHS,ub[i])
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
mt = SCNDModel("F:/scnd/Test1S2.lp")
# mt = CallModel("/home/ak121396/Desktop/instances/SCND/test01S2.lp")

mutable struct FLPModel
    lpfile::String; m::Int; n::Int; C::Array{}; B::Array{}; RHS::Array{}; signs::Array{}; vub::Array{}
    function FLPModel(datafile::String)
        datafile = readdlm("F:paper_instances/Khombole-10.txt", '#')
        data = filter(x->x!="", datafile)
        nt = findall(x->data[x][1]==' ', 1:length(data))[1]
        a = split(data[nt])
        b = parse.(Int, split(data[1]," "))
        par = Dict(zip(a,b))



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

        RHS = []
        for i=1:m
            if ub[i]==Inf
                push!(RHS,lb[i])
            else
                push!(RHS,ub[i])
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
    @constraints(m, begin
        [i=1:size(d.A,2)], dot(d.A'[i,:], α) .<= d.c[i]
    end)
    return DualSubProblem(d, α, m)
end
function JuMP.optimize!(sp::DualSubProblem, yh)
    Dy = [dot(sp.data.D[i,:],yh) for i=1:size(sp.data.D,1)]
    obj = sp.data.b .- Dy #sp.data.D * yh
    @objective(sp.m, Max, dot(obj, sp.α))
    optimize!(sp.m)
    @show termination_status(sp.m)
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

function benders_optimize!(m::Model, y::Vector{VariableRef}, sd::SubProblemData, sp_optimizer, f::Union{Function,Type}; eta_bound::Real = -1000.0) #sp_optimizer,
    subproblem = Model(sp_optimizer)
    # MOI.set(subproblem, MOI.RawParameter("CPX_PARAM_SCRIND"), false );
    # MOI.set(subproblem, MOI.RawParameter("CPX_PARAM_THREADS"),1  );
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
        Dy = [dot(dsp.data.D[i,:],yh) for i=1:size(dsp.data.D,1)]
        if res == :OptimalityCut
            @info "Optimality cut found"
            if η0 ≥ dot(α, (dsp.data.b - Dy))
                break
            else
                nopt_cons += 1
                @constraint(m, η ≥ dot(α, (dsp.data.b - Dy)))
            end
        else
            @info "Feasibility cut found"
            nfeas_cons += 1
            @constraint(m, 0 ≥ dot(α, (dsp.data.b - Dy)))
        end
        push!(cuts, (res, α))
    end
    return (m, y, cuts, nopt_cons, nfeas_cons)
end
function SetData(mt)
    rvar = findall(i->i!=1,mt.vub);
    c1r = [mt.C[1,i] for i in rvar];  #1st obj
    c2r = [mt.C[2,i] for i in rvar];  #2nd obj
    c = c1r; #c2r
    range1 = findall(i-> mt.vub[i]==1 && mt.vub[i+1]!=1, 1:length(mt.vub)-1)[end]
    range2 = findall(i-> mt.vub[i]!=1 && mt.vub[i+1]==1, 1:length(mt.vub)-1)[end]
    D = hcat(mt.B[:,1:range1],mt.B[:,range2+1:end])
    A = mt.B[:,range1+1:range2]
    b = mt.RHS
    return SubProblemData(b, D, A, c)
end

data = SetData(mt)
# objective function on y
# c1b = [mt.C[1,i] for i in bvar];
# f(v) = dot(c1b,v[1:length(c1b)])
# initialize the problem
# sp_optimizer = CPLEX.Optimizer
m = Model(sp_optimizer)
@variable(m, y[j=1:length(c1b)], Bin)
# solve and voilà
(m, y, cuts, nopt_cons, nfeas_cons) = benders_optimize!(m, y, data, CPLEX.Optimizer ,f)


bvar = findall(i->i==1,mt.vub);

c1r = [mt.C[1,i] for i in rvar];
c2b = [mt.C[2,i] for i in bvar]; c2r = [mt.C[2,i] for i in rvar];
c1 = [c1b;c1r]; c2 = [c2b;c2r]
hcat(c1,c2)'

rvar = findall(i->i!=1,mt.vub);
c1r = [mt.C[1,i] for i in rvar];  #1st obj
c2r = [mt.C[2,i] for i in rvar];  #2nd obj
c = c1r; #c2r
range1 = findall(i-> mt.vub[i]==1 && mt.vub[i+1]!=1, 1:length(mt.vub)-1)[end]
range2 = findall(i-> mt.vub[i]!=1 && mt.vub[i+1]==1, 1:length(mt.vub)-1)[end]
D = hcat(mt.B[:,1:range1],mt.B[:,range2+1:end])
A = mt.B[:,range1+1:range2]
b = mt.RHS
