# THIS FILE IS COMPETIBLE WITH JULIA 1.1 && JuMP 0.19 VERSION
include("newLexi.jl")
#############################Biobjective##########################3
using LinearAlgebra

# for i=1:n
#     j=1
#     while j>i
#         # biobj prob
#         j=j+1
#     end
# end
dF = Model( with_optimizer(GLPK.Optimizer) )
@variable(dF, lam[1:n] >= 0)
@expression( dF, ex1, dot(Y[1],lam) )
@expression( dF, ex2,  dot(Y[2],lam) )
@constraint(dF, sum(lam[i] for i=1:n)== 1)
@constraint( dF, ex1==ex2 )
@objective(dF, Min,ex1-ex2 )
optimize!(dF)

lambd = value.(lam)
for i=1:length(lambd)
    lambd[i] = round(lambd[i]; digits=4)
end
@show termination_status(dF)
termination_status(dF)
while termination_status(dF) == :OPTIMAL
    print(value.(lam))
end



#Solve biobj prob again with a new coefficients
nC1 = lambd[1]*C_1 + lambd[2]*C_2 + lambd[3]*C_3

bi = Model( with_optimizer(GLPK.Optimizer) )
@variable(bi, lam[1:n] >= 0)
@constraint(dF, sum(lam[i] for i=1:n)== 1)
@expression( dF, ex1, dot(Y[1],lam) )
@expression( dF, ex2, dot(Y[3],lam) )

@objective(dF, Min ,ex1-ex2)
optimize!(dF)

 @sprintf "%.20f" x



#INITINALISE
m = 3
u =range(1,length(S))
W = Dict()
A = Dict()
P = Dict()
for i=1:length(S)
    A[i] = []
    P[i] = filter(x-> x!=i && i<x, u)
    W[i] =[]
end

# Define W[y_i] for ∀ i
for k=1:length(S)-1
    # Find W[y_i]
    for kprime=k+1:(length(S))
        b = multi_model( with_optimizer(GLPK.Optimizer) )
        @variable(b, λ[1:3] >= 0)

        @expression( b, e1, sum{λ[i]*S[k][i],i=1:3} )
        @expression( b, e2, sum{λ[i]*S[kprime][i],i=1:3} )
        @constraint( b, sum{λ[i], i=1:3} == 1 )
        @constraint( b, e1 <= e2)

        @expression(b, ob1, )
        @expression
        obj1 = SingleObjective( -e1 )
        obj2 = SingleObjective( e2 )

        mo = get_multidata(b)
        mo.objectives = [obj1,obj2]
        # print(b)
        solve(b, method = WeightedSum())

            # push!(new_X, XE)
        new_X = mo.paretofront
        println(C_1*new_X, C_2*new_X,C_3*new_X)
        λ = unique!(mo.paretovarvalues)
        v = filter(x-> x!=k && x!=kprime ,u)
        for i=1:length(λ)
            λ[i][v[1]]=0
        end
        push!(W[k],λ)
    end
    # UPDATE A,P
    # CASE1
    for j=1:length(W[k])
        if round(dot(λ[j], S[k])) == round(dot(λ[j],S[kprime]))
            push(kprint, A[k])

        end

k=1
kprime=k+1
case = zeros(1,3)
for j=1:length(W[1])

    if round(dot(λ[j], S[kprime])) == round(dot(λ[j],S[k]))
        case[1]+=1
    elseif round(dot(λ[j], S[kprime])) < round(dot(λ[j],S[k]))
        case[2]+=1
    else


new_X = mo.paretofront
if (new_X in values(XE)) == false
    print(new_X)
end




c1 = Model(solver = CplexSolver())
# c1 = Model(solver=IpoptSolver())
@variable(c1, λ[1:3] >=0)

for j=1:(length(S)-1)
    @constraint( c1, sum{λ[i]*S[j][i], i=1:3} <= sum{λ[i]*S[j+1][i], i=1:3} )
end
@constraint( c1, sum{λ[i], i=1:3} == 1 )
@objective(c1, Min, sum{λ[i]*S[1][i], i=1:3})


print(c1)
status = solve(c1)
getobjectivevalue(c1)
println(getvalue(λ[1]))







# mc1 = Model(solver = GurobiSolver())
# @variable(mc1, x[1:4,1:4], Bin)
# @variable(mc1, w[1:2] >=0 )
# @constraint( mc1, [i=1:4], sum{ x[i,j], j=1:4}== 1)
# @constraint( mc1, [j=1:4], sum{ x[i,j], i=1:4}== 1 )
# @constraint( mc1, sum{w[i], i=1:2} == 1 )
# @constraint( mc1, w[2] == 0)
# @expression( mc1, bi1, sum{nC_1[i,j]*x[i,j], i=1:4, j=1:4} )
# @expression( mc1, bi2, sum{nC_2[i,j]*x[i,j], i=1:4, j=1:4} )
# # @expression( mc1, ex3, (sum{C_3[i,j]*x[i,j], i=1:4, j=1:4}) );
# @objective( mc1, Min, w[1]*bi1 + w[2]*bi2 );
# # print(mc1)
# status = solve(mc1);
# println("Objective value: ", getobjectivevalue(mc1))
# println("x = ", getvalue(x))
# println("w = ", getvalue(w))
#
# newx = getvalue(x)
# newy = C_1*newx
# # solution_x=getvalue(x);
# # oo = zeros(1,3)
# # oo[1] = sum(C_1[i,j]*solution_x[i,j] for i=1:4 for j=1:4)
# # oo[2] = sum(C_2[i,j]*solution_x[i,j] for i=1:4 for j=1:4)
# # ob[3] = sum(C_3[i,j]*solution_x[i,j] for i=1:4 for j=1:4)
# # println(oo[1],"\n", oo[2], "\n")
#
#
# c2 = Model(solver = GurobiSolver() )
# @variable(c2, x[1:4,1:4], Bin)
# @variable(c2, λ[1:2] >= 0 )
#
# @constraint( c2, [i=1:4], sum{ x[i,j], j=1:4}== 1)
# @constraint( c2, [j=1:4], sum{ x[i,j], i=1:4}== 1 )
# @constraint( c2, sum{λ[i], i=1:2} == 1 )
# @constraint(c2, λ[2] == 1)
# @expression( c2, ex1, (sum{C_1[i,j]*x[i,j], i=1:4, j=1:4}) )
# @expression( c2, ex2, (sum{C_2[i,j]*x[i,j], i=1:4, j=1:4}) )
# @expression( c2, ex3, (sum{C_3[i,j]*x[i,j], i=1:4, j=1:4}) );
#
# @objective( c2, Min, λ[1]*(ex1*(1/4)+ex3*(3/4)) + λ[2]*( (3/10)*ex2 +(7/10)*ex3));
#
# print(c2)
# status = solve(c2);
#
# println("Objective value: ", getobjectivevalue(c2))
# println("x = ", getvalue(x))
# println("λ = ", getvalue(λ))
# solution_x=getvalue(x);
# ob = zeros(1,3)
# ob[1] = sum(C_1[i,j]*solution_x[i,j] for i=1:4 for j=1:4)
# ob[2] = sum(C_2[i,j]*solution_x[i,j] for i=1:4 for j=1:4)
# ob[3] = sum(C_3[i,j]*solution_x[i,j] for i=1:4 for j=1:4)
# println(ob[1],"\n", ob[2], "\n", ob[3])
