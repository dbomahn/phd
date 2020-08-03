using DelimitedFiles, CPLEX, LinearAlgebra, JuMP
path = "C:\\cygwin64\\home\\AK121396\\multiobjective\\instances\\moolibrary\\KP\\"
fl = readdlm(path*"KP_p-3_n-10_ins-1.dat", '\t', String, '\n')
# files=readdir(path*"/JKU/kp/KP")
# fl = readdlm(path*"/JKU/kp/KP/"*files[1], '\t', String, '\n')
#global
obj=parse(Int,fl[1])
n=parse(Int,fl[2])
ub=parse(Int,fl[3])

################### objective function coefficient (P) matrix ################
b = fl[4:length(fl)-1]
C= ones(length(b),n)
C = round.(Int,C)
global ct=0;
for x=1:length(b)
    a = split(b[x],('[',']',','))
    aa=filter!(e->e ∉ [ "" ,"[","," ,"]"] ,a)
    for y=1:length(aa)
        c=parse(Int64,aa[y])
        C[x,y] = c
        if c==0
            global ct = ct+1;
        end
    end
end

######################## coefficient matrix (B) #############################
weight =ones(1,n)
weight = round.(Int,weight)
item = fl[length(fl)]
w1 = split(item, ('[',']',','))
w2 = filter!(e->e ∉ ["" ,"[", "]"] ,w1)
for i=1:n
    weight[i] = parse(Int64,w2[i])
end

########################## MIP model ####################################
# JuMP
function solveKP(weight,C,ub,λ_val)
    kp_m = Model(with_optimizer(CPLEX.Optimizer))
    @variable(kp_m, 0<=x[1:n]<=1 )
    @constraint( kp_m, sum(weight[i]*x[i] for i=1:n) <= ub )
    @expression( kp_m, fl[i=1:obj], sum(C[i,j]*x[j] for j=1:n) )
    @objective( kp_m, Max, sum(fl[i]*λ_val[i] for i=1:obj) )
    # objective_function(kp_m, AffExpr)
    # xval = value.(x)
    optimize!(kp_m)
    rval = Float64[]
    for i=1:obj
        P = C[((i-1)*n)+1:i*n]
        push!( rval, dot(P,value.(x)) )
    end
    # print(termination_status(kp_m),rval)
    return rval
end

##########################  Calculate λ  ############################
function calλ(R)
    lam = Model(with_optimizer(CPLEX.Optimizer))
    @variable(lam, λ[1:3])
    @constraint( lam, sum(λ[1:3]) == 1)  #Pkg add https://github.com/VMLS-book/VMLS.jl
    @constraint( lam, dot(λ,R[1]) == dot(λ,R[2]) )
    @constraint( lam, dot(λ,R[2]) == dot(λ,R[3])     )
    @constraint( lam, dot(λ,R[1]) == dot(λ,R[3]) )
    # @objective( lam, Max, -sum(λ[i]*R[i] for i=1:obj) )
    optimize!(lam)
    # obj_value = JuMP.objective_value(lam)
    return value.(λ)
end

############################  Find a new λ #############################
# function findnewλ(fy,R)
#     lam_m = Model(with_optimizer(CPLEX.Optimizer))
#     @variable(lam_m, nλ[1:obj] )
#
#     @constraint( lam_m, dot(nλ,fy) <= dot(nλ,R[1]) )
#     @constraint( lam_m, dot(nλ,fy) <= dot(nλ,R[2]) )
#     @constraint( lam_m, dot(nλ,fy) <= dot(nλ,R[3]) )
#     @expression( lam_m, ex1, dot(nλ,fy) - dot(nλ,R[1]) )
#     @expression( lam_m, ex2, dot(nλ,fy) - dot(nλ,R[2]) )
#     @expression( lam_m, ex3, dot(nλ,fy) - dot(nλ,R[3]) )
#
#     @objective( lam_m, Min, ex1+ex2+ex3 )
#     optimize!(lam_m)
#     if termination_status(lam_m) == MOI.OPTIMAL && [value.(nλ)] ∈ R
#         return newλ = [value.(nλ)]
#     end
#     # print(termination_status(lam_m))
# end


# MultiJuMP
# kp_m = multi_model(solver = CplexSolver(), linear = true)
# multi_model( solver = CplexSolver(), linear=true)
# @variable(kp_m, 0 <= x[1:n] <= 1)
# @constraint( kp_m, sum{ weight[i]*x[i], i=1:n } <= ub )
#
# @expression( kp_m, fl1, sum{C[i,j]*x[j], i=1, j=1:n} )
# @expression( kp_m, fl2, sum{C[i,j]*x[j], i=2, j=1:n} )
# @expression( kp_m, fl3, sum{C[i,j]*x[j], i=3, j=1:n} )
# const obj1 = SingleObjective( fl1, sense = :Max )
# const obj2 = SingleObjective( fl2, sense = :Max )
# const obj3 = SingleObjective( fl3, sense = :Max )
#
# multi = get_multidata(kp_m)
# multi.objectives= [obj1, obj2, obj3]
#
# solve(kp_m, method = WeightedSum())
#
# getobjectivevalue(kp_m)
# println(getvalue(x))
# getvalue(obj3)

# function solveMIP(obj,λ_val,R)
#     kp_m = Model(with_optimizer(CPLEX.Optimizer))
#     @variable(kp_m, r[1:obj]>=0 )
#
#     @constraint( kp_m, dot(λ_val,R[1]) >= dot(λ_val,r[1:obj])  )
#     @constraint( kp_m, dot(λ_val,R[2]) >= dot(λ_val,r[1:obj]) )
#     @constraint( kp_m, dot(λ_val,R[3]) >= dot(λ_val,r[1:obj]) )
#
#     @expression( kp_m, ex1, dot(λ_val,R[1]) - dot(λ_val,r[1:obj])  )
#     @expression( kp_m, ex2, dot(λ_val,R[2]) - dot(λ_val,r[1:obj]) )
#     @expression( kp_m, ex3, dot(λ_val,R[3]) - dot(λ_val,r[1:obj]) )
#
#     @objective( kp_m, Min, ex1+ex2+ex3 )
#     # objective_function(kp_m, AffExpr)
#     # objective_value(kp_m)
#     optimize!(kp_m)
#     if termination_status(kp_m) == MOI.OPTIMAL
#         return r_val = [objective_value(kp_m)]
#     end
# end
