using JuMP,CPLEX,LinearAlgebra,MathProgBase,CPUTime,DelimitedFiles,StatsBase
const MPB = MathProgBase

mutable struct FLPdata
    datafile::String; d::Matrix{}; W::Matrix{}; par::Dict{};
    opencost::Array{}; popsize::Array{}; capa::Array{}; V::Int
    function FLPdata(datafile::String)
        dt = readdlm(datafile,'\t') #
        # dt = readdlm("F:paper_instances/Thienaba-10.txt",'\n')
        hash = findall(i->'#' in dt[i], 1:length(dt))
        sc = dt[hash[end]+1:end]
        nodes = length(split(sc[1]," "))
        W = zeros(Int,length(sc),nodes)
        for i=1:length(sc)
            for j=1:nodes
                scn = parse.(Int,split(sc[i]," "))
                W[i,j] = scn[j]
            end
        end
        stca = split.(dt[hash[end-1]+1:hash[end-1]+nodes]," ")
        capa = [parse(Int,stca[i][2]) for i=1:nodes]
        pop = split.(dt[hash[end-2]+1:hash[end-2]+nodes]," ")
        popsize = [parse(Int,pop[i][2]) for i=1:nodes]
        opening = split.(dt[hash[end-3]+1:hash[end-3]+nodes], " ")
        opencost= [parse(Int,opening[i][2]) for i=1:nodes]
        v = length(dt[hash[end-5]+1:hash[end-4]-1])
        dist = dt[hash[end-4]+1:hash[end-4]+v]
        d = zeros(Int,v,v)
        for i=1:v
            for j=1:v
                distm = parse.(Int,split(dist[i]," "))
                d[i,j] = distm[j]
            end
        end
        a = split(dt[hash[1]]," ")[2:end]
        b = parse.(Int,split(dt[hash[2]+1]," "))
        par = Dict(zip(a,b))
        V = v-1

        new(datafile,d,W,par,opencost,popsize,capa,V)
    end
end
dpath = "/home/ak121396/Desktop/instances/ORspectrum_instances/Thienaba-10.txt"
dt = FLPdata(dpath)

####### Master problem
mas = Model(CPLEX.Optimizer); set_silent(mas);
MOI.set(mas, MOI.RawParameter("CPX_PARAM_SCRIND"), false )
MOI.set(mas, MOI.RawParameter("CPX_PARAM_THREADS"),1  )
@variable(mas, z[1:dt.V], Bin)
@variable(mas, θ[1:dt.par["nSamples"]])
@objective(mas, Min, dot(dt.opencost,z) + 1/dt.par["nSamples"]*sum(θ))

function solveMP(m::Model,W,N,zbar,𝜋b,σb,δb)
    @constraint(m, [v=1:N], -sum(W[v,:]) <= θ[v] <= 0)
    # @constraint(m, [l in Sv,v=1:N], θ[v] >= Q-sum(𝜋b[j]*dt.capa[j]*z[j] for j=1:dt.V)-sum(σb[i,j]*W[v,i]*z[j] for i=1:dt.V for j=1:dt.V)-sum(W[i]*δb[i] for i=1:dt.V) )
    # @constraint(m, [v=1:N], -sum(z[j]*capa[j]) <= θ[v])
    optimize!(m)
    return objective_value(m),value.(z)
end
function solveSP(zbar,capa,V,W,s)
    sub = Model(CPLEX.Optimizer)
    @variable(sub, λ[1:dt.V] >= 0)
    @variable(sub, 𝜋[1:dt.V] >= 0)
    @variable(sub, σ[1:dt.V,1:dt.V]>=0)
    @variable(sub, δ[1:dt.V]>=0)
    @constraint(sub, [j=1:dt.V], λ[j]+𝜋[j] >= 1)
    @constraint(sub, [i=1:dt.V,j=1:dt.V], λ[j]-σ[i,j]-δ[i] <= 0)
    @objective(sub, Max, -sum(𝜋[j]*capa[j]*zbar[j] for j=1:dt.V)-sum(σ[i,j]*W[s,i]*zbar[j] for i=1:dt.V for j=1:dt.V)-sum(W[s,i]*δ[i] for i=1:dt.V) );
    optimize!(sub)
    st = termination_status(sub)
    if st == MOI.OPTIMAL
        return (:OptimalityCut, objective_value(sub),value.(λ),value.(𝜋),value.(σ),value.(δ))
    elseif st == MOI.DUAL_INFEASIBLE
        return (:FeasibilityCut, objective_value(sub),value.(λ),value.(𝜋),value.(σ),value.(δ))
    else
        error("DualSubProblem error: status $st")
    end
    # return objective_value(sub),value.(λ),value.(𝜋),value.(σ),value.(δ)
end


########### Sub problem
# function benders_optimize!(zbar,W,V)
#     st = termination_status(sub)
#     if st == MOI.OPTIMAL
#         return (:OptimalityCut,value.(λ),value.(𝜋),value.(σ),value.(δ))
#     elseif st == MOI.DUAL_INFEASIBLE
#         return (:FeasibilityCut,value.(λ),value.(𝜋),value.(σ),value.(δ))
#     else
#         error("DualSubProblem error: status $st")
#     end
# end
let
    # m = defineMP(dt.V,dt.par["nSamples"],dt.opencost,[0.5,0.5])
    S = length(dt.par["nSamples"])
    UB = Inf #ones(S)*sum(dt.opencost)  #UB of the whole objfunction
    LB = -Inf #LB
    sub_obj = zeros(S)
    ϵ = 1e-6
    nopt_cons, nfeas_cons = (0, 0)
    (zbar,λ,𝜋,σ,δ) = zeros(dt.V),zeros(S,dt.V),zeros(S,dt.V),zeros(S,dt.V),zeros(S,dt.V)
    iter = 1 #iteration
    cuts = Tuple{Symbol, Vector{Float64}}[]

    while UB-LB ≥ ϵ
        for s=1:S
            @show res,sub_obj[s],λ,𝜋,σ,δ = solveSP(zbar,dt.capa,dt.V,dt.W,s)
            if res == :OptimalityCut
                @info "Optimality cut found"
                if sub_obj[s] ≥ UB
                    println("BREAK::LB ≥ UB")
                    break
                else
                    nopt_cons += 1
                    @constraint(mas, θ[s] ≥ sub_obj[s] -sum(𝜋[j]*dt.capa[j]*(zbar[j]-z[j]) for j=1:dt.V)-sum(σ[i,j]*dt.W[s,i]*zbar[j] for i=1:dt.V for j=1:dt.V)-sum(dt.W[s,i]*δ[i] for i=1:dt.V) )
                end
            else
                @info "Feasibility cut found"
                nfeas_cons += 1
                @constraint(mas, 0 ≥ sub_obj[s] -sum(𝜋[j]*dt.capa[j]*(zbar[j]-z[j]) for j=1:dt.V)-sum(σ[i,j]*dt.W[s,i]*zbar[j] for i=1:dt.V for j=1:dt.V)-sum(dt.W[s,i]*δ[i] for i=1:dt.V) )
            end
            push!(cuts, (res,sub_obj[s]))
        end
        @show LB = mean(sub_obj)
        mas_obj,zbar = solveMP(mas,dt.W,S,zbar,𝜋,σ,δ);
        UB = mas_obj;
        println("$iter UB: $UB...LB: $LB")# "...Sub: $sub_obj")
    end
    return z, nopt_cons, nfeas_cons
end
value.(z)
##############################  BOSFLP Mathematical model  ######################################
flp = Model(CPLEX.Optimizer)
MOI.set(flp, MOI.RawParameter("CPX_PARAM_SCRIND"), false )
MOI.set(flp, MOI.RawParameter("CPX_PARAM_THREADS"),1  )
@variable(flp, z[1:dt.V], Bin)
@variable(flp, y[1:dt.V,1:dt.V,v=1:dt.par["nSamples"]] >=0 )
@variable(flp, u[1:dt.V,1:dt.par["nSamples"]] >=0)

#objectives
@constraint(flp, obj1, dot(dt.opencost,z)<=0)
@constraint(flp, obj2, -1/dt.par["nSamples"]*(sum(u[j,v] for j=1:dt.V for v=1:dt.par["nSamples"])) <= 0 )
#con13
@constraint( flp, [j=1:dt.V,v=1:dt.par["nSamples"]], -u[j,v] >= -sum(y[i,j,v] for i=1:dt.V)  )
#con14
@constraint( flp, [j=1:dt.V,v=1:dt.par["nSamples"]], -u[j,v] >= -dt.capa[j]*z[j] )
#con15
@constraint( flp, [i=1:dt.V,j=1:dt.V,v=1:dt.par["nSamples"]], -y[i,j,v] >= -dt.W[v,i]*z[j])
#con16
@constraint( flp, [i=1:dt.V,v=1:dt.par["nSamples"]], -sum(y[i,j,v] for j=1:dt.V) >= -dt.W[v,i])

write_to_file(flp, dpath*".lp")



function loadlp(filename,solver=CplexSolver(CPX_PARAM_SCRIND=0))
    model=buildlp([-1,0],[2 1],'<',1.5, solver) # create dummy model with correct solver
    MPB.loadproblem!(model,filename) # load what we actually want
    return model
end
lpmodel = loadlp(dpath*".lp")
Bmtx = MPB.getconstrmatrix(lpmodel);
B = Bmtx[3:end,:];C = Bmtx[1:2,:]; vub = MPB.getvarUB(lpmodel)
m,n=size(B)
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

#############################    Julia pkg    ##################################
using StochasticPrograms
dpath = "/home/ak121396/Desktop/instances/ORspectrum_instances/Thienaba-10.txt"
dt = FLPdata(dpath)

@stochastic_model bosflp begin
    @stage 1 begin
        @decision(bosflp, z[1:dt.V], Bin)
        @objective(bosflp, Min, dot(dt.opencost,z) )
    end
    @stage 2 begin
        @uncertain W
        @recourse(bosflp, 0 <= y[1:dt.V,1:dt.V])
        @recourse(bosflp, 0 <= u[1:dt.V])
        @objective(bosflp, Max, sum(u))
        #con5
        @constraint( bosflp, [j=1:dt.V], u[j] <= sum(y[i,j] for i=1:dt.V)  )
        #con6
        @constraint( bosflp, [j=1:dt.V], u[j] <= dt.capa[j]*z[j] )
        #con7
        @constraint( bosflp, [i=1:dt.V,j=1:dt.V], y[i,j] <= dt.W[i]*z[j] )
        #con8
        @constraint( bosflp, [i=1:dt.V], sum(y[i,j] for j=1:dt.V) <= dt.W[i])
    end
end

Ξ = []
for i=1:dt.par["nSamples"]
    ξ = @scenario W = dt.W[i,:] probability = (1/dt.par["nSamples"])
    push!(Ξ, ξ)
end

# sp = instantiate(bosflp, [Ξ[i] for i=1:dt.par["nSamples"]], optimizer = CPLEX.Optimizer)
# optimize!(sp)
# println("Objective value in scenario 2: $(objective_value(sp, 2))")
# println("Optimal decision: $(optimal_decision(sp))")
# objective_value(sp)
# optimal_decision(sp)


flp_lshaped = instantiate(bosflp, [Ξ[i] for i=1:dt.par["nSamples"]], optimizer = LShaped.Optimizer)
# print(flp_lshaped)
set_optimizer_attribute(flp_lshaped, MasterOptimizer(), CPLEX.Optimizer)
set_optimizer_attribute(flp_lshaped, SubProblemOptimizer(), CPLEX.Optimizer)
set_suboptimizer_attribute(flp_lshaped, MOI.RawParameter("CPX_PARAM_SCRIND"), false) # Silence CPLEX
set_suboptimizer_attribute(flp_lshaped, MOI.RawParameter("CPX_PARAM_THREADS"), 1) #
# "CPXXbtrans"
optimize!(flp_lshaped; cache = true)
# println("Objective value in scenario 2: $(objective_value(flp_lshaped, 2))")
objective_value(flp_lshaped)
optimal_decision(flp_lshaped)
termination_status(flp_lshaped)


lazy_model = Model(CPLEX.Optimizer)
@variable(lazy_model, z[1:dt.V] , Bin)
@variable(lazy_model, θ >= -1000)
@objective(lazy_model, Min, θ)
print(lazy_model)
k=0
function my_callback(cb_data)
    global k += 1
    x_k = callback_value.(cb_data, x)
    θ_k = callback_value(cb_data, θ)
    lower_bound = c_1' * x_k + θ_k
    ret = solve_subproblem(x_k)
    upper_bound = c_1' * x_k + c_2' * ret.y
    gap = (upper_bound - lower_bound) / upper_bound
    print_iteration(k, lower_bound, upper_bound, gap)
    if gap < ABSOLUTE_OPTIMALITY_GAP
        println("Terminating with the optimal solution")
        return
    end
    cut = @build_constraint(θ >= ret.obj + -ret.π' * A_1 * (x .- x_k))
    MOI.submit(model, MOI.LazyConstraint(cb_data), cut)
    return
end
MOI.set(lazy_model, MOI.LazyConstraintCallback(), my_callback)

1
