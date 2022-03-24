using JuMP,CPLEX,LinearAlgebra,MathProgBase,CPUTime
const MPB = MathProgBase
# https://github.com/matbesancon/SimpleBenders.jl/blob/master/test/runtests.jl
# https://matbesancon.xyz/post/2019-05-08-simple-benders/
# https://co-at-work.zib.de/slides/Donnerstag_24.9/Benders_decomposition-Fundamentals.pdf

mutable struct CallModel
    lpfile::String; m::Int; n::Int; C::Array{}; B::Array{}; RHS::Array{}; signs::Array{}; vub::Array{}
    function CallModel(lpfile::String)
        lpmodel=buildlp([-1,0],[2 1],'<',1.5, CplexSolver(CPX_PARAM_SCRIND=0))
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
    signs::Array{}
end

"""
    DualSubProblem
Benders dual subproblem:
max (b - Dŷ)ᵀ α
s.t. Aᵀα ⩽ c
     α ⩾ 0
"""
struct DualSubProblem
    data::SubProblemData
    α::Vector{VariableRef}
    m::Model
end

function DualSubProblem(d::SubProblemData, m::Model)
    # α = @variable(m, α[i = 1:size(d.A, 1)] >= 0)
    eqpoint = findall(i-> d.signs[i]=="l", 1:length(d.signs))
    @variable(m, α[i = 1:size(d.A, 1)] )
    for i in eqpoint
      @constraint(m, α[i] >= 0)
    end
    @constraint(m, d.A' * α .<= d.c)
    return DualSubProblem(d, α, m)
end

function JuMP.optimize!(sp::DualSubProblem, ŷ)
    obj = sp.data.b .- sp.data.D * ŷ
    @objective(sp.m, Max, obj' * sp.α)
    optimize!(sp.m)
    st = termination_status(sp.m)
    if st == MOI.OPTIMAL
        α = JuMP.value.(sp.α)
        return (:OptimalityCut, α)
    elseif st == MOI.DUAL_INFEASIBLE
        return (:FeasibilityCut, α)
    else
        error("DualSubProblem error: status $st")
    end
end

function SubData(mt,range1,range2,firstcon,w)
    rvar = findall(i->i!=1,mt.vub);
    c1r = [mt.C[1,i] for i in rvar];  #1st obj
    c2r = [mt.C[2,i] for i in rvar];  #2nd obj
    c = c1r*w + c2r*(1-w)
    D = hcat(mt.B[firstcon[end]+1:end,1:range1],mt.B[firstcon[end]+1:end,range2+1:end])
    A = mt.B[firstcon[end]+1:end,range1+1:range2]
    b = mt.RHS[firstcon[end]+1:end]
    signs = mt.signs[firstcon[end]+1:end]
    # D0 = vcat(mt.B[1:firstcon[1]-1,1:range1],mt.B[firstcon[end]+1:end,1:range1])
    # D1 = vcat(mt.B[1:firstcon[1]-1,range2+1:end],mt.B[firstcon[end]+1:end,range2+1:end])
    # D = hcat(D0,D1)
    # A = vcat(mt.B[1:firstcon[1]-1,range1+1:range2],mt.B[firstcon[end]+1:end,range1+1:range2])
    # b = vcat(mt.RHS[1:firstcon[1]-1],mt.RHS[firstcon[end]+1:end])
    # signs = vcat(mt.signs[1:firstcon[1]-1],mt.signs[firstcon[end]+1:end])
    return SubProblemData(b, D, A, c, signs)
end

function benders_optimize!(m::Model, y::Vector{VariableRef}, sd::SubProblemData, sp_optimizer, f::Union{Function,Type}; eta_bound::Real = -1000.0)
    subproblem = Model(sp_optimizer)
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
            if η0 ≥ α' * (dsp.data.b - dsp.data.D * ŷ)
                break
            else
                nopt_cons += 1
                @constraint(m, η ≥ α' * (dsp.data.b - dsp.data.D * y))
            end
        else
            @info "Feasibility cut found"
            nfeas_cons += 1
            @constraint(m, 0 ≥ α' * (dsp.data.b - dsp.data.D * y))
        end
        push!(cuts, (res, α))
    end
    return (m, y, cuts, nopt_cons, nfeas_cons)
end
######################       Importing Data         #####################
# mt = CallModel("F:/scnd/Test1S2.lp")
# mt = CallModel("/home/k2g00/k2g3475/scnd/lp/test01S2.lp")
mt = CallModel("/home/ak121396/Desktop/instances/SCND/test01S2.lp")
bvar = findall(i->i==1,mt.vub); rvar = findall(i->i!=1,mt.vub);
c1b = [mt.C[1,i] for i in bvar]; #c2b is zeros
# c1r = [mt.C[1,i] for i in rvar]; c2r = [mt.C[2,i] for i in rvar];
range1 = findall(i-> mt.vub[i]==1 && mt.vub[i+1]!=1, 1:length(mt.vub)-1)[1]
range2 = findall(i-> mt.vub[i]!=1 && mt.vub[i+1]==1, 1:length(mt.vub)-1)[1]
firstcon = findall(i->sum(mt.B[i,range1+1:range2])==0, 1:length(mt.RHS))
mas0 = mt.B[firstcon[1]:firstcon[end],:]
mascon = hcat(mas0[:,1:range1],mas0[:,range2+1:end])
mRHS = mt.RHS[firstcon[1]:firstcon[end]]
msigns = mt.signs[firstcon[1]:firstcon[end]]
# W = [0.3,0.6,0.9]

# initialize the problem
w=1
data = SubData(mt,range1,range2,firstcon,w)
sp_optimizer = CPLEX.Optimizer
function defineMP(c1b,w,msigns,mascon,mRHS)
    m = Model(CPLEX.Optimizer); set_silent(m);
    MOI.set(m, MOI.RawParameter("CPX_PARAM_SCRIND"), false )
    MOI.set(m, MOI.RawParameter("CPX_PARAM_THREADS"),1  )
    c_y = c1b*w
    @variable(m, y[j=1:length(c_y)], Bin); # @variable(m, 0<=y[j=1:length(c_y)]<=1);
    for k=1:length(msigns)
        if msigns[k] == "l"
            @constraint( m, dot(mascon[k,:],y) >= mRHS[k] );
        elseif msigns[k] == "u"
            @constraint( m, dot(mascon[k,:],y) <= mRHS[k] );
        else
            @constraint( m, dot(mascon[k,:],y) == mRHS[k] );
        end
    end
    #obj function for the first stage variables
    f(y) = dot(c_y, y);
    return m,y,f
end
m,y,fy = defineMP(c1b,w,msigns,mascon,mRHS)
# solve and voilà
runtime = @CPUelapsed (m, y, cuts, nopt_cons, nfeas_cons) = benders_optimize!(m, y, data, sp_optimizer, fy)

objective_value(m)
value.(y)
findall(i->i!=0, value.(η))
cuts
nfeas_cons
# file = open("/home/k2g00/k2g3475/scnd/lp/testBD.txt","w")
# output = open("/home/ak121396/Desktop/testBD.txt")
# writedlm(file, runtime, objective_value(m), )
# close(file)

sol_w = []
for w ∈ W
    data = SubData(mt,range1,range2,firstcon,w)
    m = Model(sp_optimizer); set_silent(m);
    MOI.set(m, MOI.RawParameter("CPX_PARAM_SCRIND"), false )
    # MOI.set(m, MOI.RawParameter("CPX_PARAM_THREADS"),1  )
    @variable(m, y[j=1:length(c_y)], Bin);  # @variable(m, 0<=y[j=1:length(c_y)]<=0);
    for k=1:length(msigns)
        if msigns[k] == "l"
            @constraint( m, dot(mascon[k,:],y) >= mRHS[k] )
        elseif msigns[k] == "u"
            @constraint( m, dot(mascon[k,:],y) <= mRHS[k] )
        else
            @constraint( m, dot(mascon[k,:],y) == mRHS[k] )
        end
    end
    #obj function for the first stage variables
    c_y = c1b*w
    f(y) = dot(c_y, y[1:length(c_y)])
    runtime = @CPUelapsed (m, y, cuts, nopt_cons, nfeas_cons) = benders_optimize!(m, y, data, sp_optimizer, f)
    push!(sol_w,[runtime, m, y, cuts, nopt_cons, nfeas_cons])
end
function test_data()
    c = [2., 3.]
    A = [1 2;2 -1]
    D = zeros(2, 1) .+ [1, 3]
    b = [3, 4]
    return SubProblemData(b, D, A, c)
end

data = test_data()
# objective function on y
f(v) = 2v[1]
# initialize the problem
m = Model(CPLEX.Optimizer)
@variable(mas, y[j=1:1] >= 0)
# solve and voilà
(mas, y, cuts, nopt_cons, nfeas_cons) = benders_optimize!(mas, y, data, sp_optimizer, f)

objective_value(m)
JuMP.value.(y)

using JuMP,CPLEX, LinearAlgebra
mt = CallModel("/home/ak121396/Desktop/instances/SCND/test01S2.lp")
bvar = findall(i->i==1,mt.vub); rvar = findall(i->i!=1,mt.vub);
c1b = [mt.C[1,i] for i in bvar]; #c2b is zeros
# c1r = [mt.C[1,i] for i in rvar]; c2r = [mt.C[2,i] for i in rvar];
range1 = findall(i-> mt.vub[i]==1 && mt.vub[i+1]!=1, 1:length(mt.vub)-1)[1]
range2 = findall(i-> mt.vub[i]!=1 && mt.vub[i+1]==1, 1:length(mt.vub)-1)[1]
firstcon = findall(i->sum(mt.B[i,range1+1:range2])==0, 1:length(mt.RHS))
mas0 = mt.B[firstcon[1]:firstcon[end],:]
mascon = hcat(mas0[:,1:range1],mas0[:,range2+1:end])
mRHS = mt.RHS[firstcon[1]:firstcon[end]]
msigns = mt.signs[firstcon[1]:firstcon[end]]
w = 1
####### Master problem
mas = Model(CPLEX.Optimizer); set_silent(mas);
MOI.set(mas, MOI.RawParameter("CPX_PARAM_SCRIND"), false )
MOI.set(mas, MOI.RawParameter("CPX_PARAM_THREADS"),1  )
c_y = c1b#*w
@variable(mas, y[j=1:length(c_y)], Bin); # @variable(mas, 0<=y[j=1:length(c_y)]<=1);
@variable(mas, η); # @variable(mas, 0<=y[j=1:length(c_y)]<=1);
@objective(mas, Min, dot(c_y,y)+η)
for k=1:length(msigns)
    if msigns[k] == "l"
        @constraint( mas, dot(mascon[k,:],y) >= mRHS[k] );
    elseif msigns[k] == "u"
        @constraint( mas, dot(mascon[k,:],y) <= mRHS[k] );
    else
        @constraint( mas, dot(mascon[k,:],y) == mRHS[k] );
    end
end
function solveMP(γ,sd::SubProblemData)
    @constraint(mas, η ≥ γ * (sd.b - (sd.D * y)))
    optimize!(mas)
    return objective_value(mas)
end

########### Sub problem
function solve_sub(ybar,sd::SubProblemData)
    eqpoint = findall(i-> sd.signs[i]=="l" && sd.signs[i+1]=="s", 1:length(sd.signs)-1)[1]
    obj = sd.b .- sd.D*ybar
    sub = Model(CPLEX.Optimizer)
    @variable(sub, α[1:eqpoint] >= 0)
    @variable(sub, β[1:length(sd.signs)-eqpoint])
    @variable(sub, γ[1:length(sd.signs)]>=0)
    for i=1:eqpoint
        @constraint(sub,γ[i]==α[i])
    end
    @objective(sub, Max, dot(obj,γ) );
    @constraint(sub, sd.A' * γ .<= sd.c)
    optimize!(sub)
    return objective_value(sub),value.(γ)
end

let
    UB = Inf
    LB = -Inf
    Delta = 0
    ybar = zeros(length(y))
    iter = 1
    while UB-LB>Delta
        sub_obj,γ = solve_sub(ybar,data);
        LB = max(LB,sub_obj+dot(c_y,ybar));
        mas_obj = solveMP(γ,data);
        ybar = value.(y);
        UB = mas_obj
        println("$it UB: $UB...LB: $LB...Sub: $sub_obj")
        iter+=1
    end
end

# struct SubProblemData
#     b::Vector{Float64}
#     D::Matrix{Float64}
#     A::Matrix{Float64}
#     c::Vector{Float64}
#     signs::Array{}
# end
# struct DualSubProblem
#     data::SubProblemData
#     α1::Vector{VariableRef}
#     α2::Vector{VariableRef}
#     m::Model
# end
# function DualSubProblem(d::SubProblemData, m::Model)
#     eqpoint = findall(i-> d.signs[i]=="l" && d.signs[i+1]=="s", 1:length(d.signs)-1)[1]
#     # α = @variable(m, α[i = 1:length(d.signs)])# >= 0)
#     α1 = @variable(m, α1[1:eqpoint] >=0)
#     α2 = @variable(m, α2[1:length(d.signs)-eqpoint])
#     for i=1:length(d.signs)
#         if sum(d.A'[i,:])!=0
#             @constraint(m, dot(d.A'[i,1:eqpoint],α1) + dot(d.A'[i,1+eqpoint:length(d.signs)],α2) <= d.c[i])
#         else
#             println(i,"th column is empty")
#         end
#     end
#     # @constraint(m, d.A' * α .<= d.c)
#     return DualSubProblem(d, α1, α2, m)
# end
# function JuMP.optimize!(sp::DualSubProblem, ŷ)
#     obj = sp.data.b .- sp.data.D * ŷ
#     @objective(sp.m, Max, dot(obj[1:length(sp.α1)]',sp.α1) + dot(obj[1+length(sp.α1):end]',sp.α2))
#     optimize!(sp.m)
#     st = termination_status(sp.m)
#     α1 = JuMP.value.(sp.α1); α2 = JuMP.value.(sp.α2)
#     α = [α1;α2]
#     if st == MOI.OPTIMAL
#         return (:OptimalityCut, α)
#     elseif st == MOI.DUAL_INFEASIBLE
#         return (:FeasibilityCut, α)
#     else
#         error("DualSubProblem error: status $st")
#     end
# end
# function benders_optimize!(m::Model, y::Vector{VariableRef}, sd::SubProblemData, sp_optimizer, f::Union{Function,Type}; eta_bound::Real = -1000.0)
#     subproblem = Model(sp_optimizer)
#     dsp = DualSubProblem(sd, subproblem)
#     @variable(m, η >= eta_bound)
#     @objective(m, Min, f(y) + η)
#     optimize!(m)
#     st = termination_status(m) #MOI.get(m, MOI.TerminationStatus())
#     # restricted master has a solution or is unbounded
#     nopt_cons, nfeas_cons = (0, 0)
#     # @info "Initial status $st"
#     cuts = Tuple{Symbol, Vector{Float64}}[]
#     # t0 = time()
#     while (st == MOI.DUAL_INFEASIBLE) || (st == MOI.OPTIMAL) # && time()-t0<TL
#         optimize!(m)
#         st = termination_status(m) # MOI.get(m, MOI.TerminationStatus())
#         ŷ = JuMP.value.(y)
#         η0 = JuMP.value(η)
#         (res, α) = optimize!(dsp, ŷ)
#         if res == :OptimalityCut
#             # @info "Optimality cut found"
#             if η0 ≥ α' * (dsp.data.b - dsp.data.D * ŷ)
#                 break
#             else
#                 nopt_cons += 1
#                 @constraint(m, η ≥ α' * (dsp.data.b - dsp.data.D * y))
#             end
#         else
#             # @info "Feasibility cut found"
#             nfeas_cons += 1
#             @constraint(m, 0 ≥ α' * (dsp.data.b - dsp.data.D * y))
#         end
#         push!(cuts, (res, α))
#     end
#     return (m, y, cuts, nopt_cons, nfeas_cons)
# end
