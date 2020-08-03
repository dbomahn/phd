using MultiJuMP, JuMP, GLPK, GLPKMathProgInterface,vOptGeneric
# mm = multi_model(solver = CplexSolver(), linear=true)
# @variable(mm, x[1:4,1:4], Bin)
# @constraint( mm, [i=1:4], sum{ x[i,j], j=1:4}== 1)
# @constraint( mm, [j=1:4], sum{ x[i,j], i=1:4}== 1 )
#
# @expression( mm, ex1, sum{C_1[i,j]*x[i,j], i=1:4, j=1:4} )
# @expression( mm, ex2, sum{C_2[i,j]*x[i,j], i=1:4, j=1:4} )
# @expression( mm, ex3, sum{C_3[i,j]*x[i,j], i=1:4, j=1:4} )
# const obj1 = SingleObjective( ex1 )
# const obj2 = SingleObjective( ex2 )
# const obj3 = SingleObjective( ex3 )
#
# multi = get_multidata(mm)
# multi.objectives= [obj1, obj2, obj3]
# solve(mm, method = WeightedSum())
#
# plot(multi)
# title!("Extrema of Pareto front")
# print(getvalue(a))

#VOptGeneric model
m = vModel(solver = CplexSolver())
C_1 =[2 5 4 7; 3 3 5 7; 3 8 4 2; 6 5 2 5]
C_2 = [3 3 6 2; 5 3 7 3; 5 2 7 4; 4 6 3 5]
C_3 = [4 2 5 3; 5 3 4 3; 4 3 5 2; 6 4 7 3]

@variable(m, x[1:size], Bin)
@addobjective(m, Max, dot(x, p1))
@addobjective(m, Max, dot(x, p2))
@constraint(m, dot(x, w) <= c)



###########################################################################
include("Lexico.jl")
####################################Biobjective##########################3
using MultiJuMP,JuMP, LinearAlgebra
#INITINALISE
nobj = 3
u =range(1,length(S))
W = Dict()
A = Dict()
P = Dict()
for i=1:length(S)
    A[i] = []
    P[i] = filter(x-> x!=i && i<x, u)
    W[i] =[]
end

k=1
kprime=3

b = multi_model(solver = CplexSolver(), linear=true)
@variable(b, λ[1:nobj] >= 0)

@expression( b, e1, sum{λ[i]*S[k][i],i=1:nobj} )
@expression( b, e2, sum{λ[i]*S[kprime][i],i=1:nobj} )
@constraint( b, sum{λ[i], i=1:nobj} == 1 )
@constraint( b, e1 <= e2)

obj1 = SingleObjective( -e1 )
obj2 = SingleObjective( e2 )

mo = get_multidata(b)
mo.objectives = [obj1,obj2]
# print(b)
solve(b, method = WeightedSum())

    # push!(new_X, XE)
# new_λ = unique(mo.paretofront)
λ = unique!(mo.paretovarvalues)
v = filter(x-> x!=k && x!=kprime ,u)
# for i=1:length(λ)
#     λ[i][v[1]]=0
# end
push!(W[k],λ)

print(λ)
