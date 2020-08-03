using DelimitedFiles, CPLEX, LinearAlgebra, JuMP
path = "C:\\cygwin64\\home\\AK121396\\multiobjective\\instances\\moolibrary\\AP\\"
f = readdlm(path*"AP_p-3_n-10_ins-1.dat", '\t', String, '\n')
obj=parse(Int,f[1])
n=parse(Int,f[2])
################### objective function coefficient (P) matrix ################
P=[]
for i=3:length(f)
    a = split(f[i], ('[',']',','))
    a = filter!(e->e ∉ [ "" ,"[","," ,"]"] ,a)
    if length(a)>2
        push!(P,a)
    end
end

P2 = ones(obj,n*n)
P2 = round.(Int,P2)

global ct=0;
for x=1:length(P)
    for y=1:n
        p=parse(Int64,P[x][y])
        idx = Int(floor((x-1)/n))
        ind = ((x-1)*n+y)%(n*n)
        if ind != 0
            P2[idx+1,((x-1)*n+y)%(n*n)] = p
        else
            P2[idx+1,(n*n)] = p
        end
        if p==0
            ct =ct+1;
        end
    end
end


######################## coefficient matrix (B) #############################
A1 =zeros(n,n*n)
A2 =zeros(n,n*n)
A1  = round.(Int,A1)
A2 = round.(Int,A2)

for i=1:n
    for j=1:n
        A1[i,((i-1)*n)+j] = 1
            if A1[i,((i-1)*n)+j]==1
            end
    end
end

for i=1:n
    for j=1:n
        A2[i,((j-1)*n)+i] = 1
        if A2[i, ((j-1)*n)+i]==1
        end
    end
end
B=[A1;A2]

########################## MIP model ####################################
# JuMP
function solveAP(B,P2,λ_val)
    ap_m = Model(with_optimizer(CPLEX.Optimizer))
    # λ_val =[1/3,1/3,1/3]
    # coeff = [P2[i,:]*λ_val[i] for i=1:3]
    @variable( ap_m, 0<=x[1:n*n]<=1 )
    @constraint( ap_m, cons[i=1:2*n], sum(B[i,j]*x[j] for j=1:n*n) == 1 )
    @expression( ap_m, fl[i=1:obj], sum(P2[i,j]*x[j] for j=1:n*n) )
    @objective( ap_m, Min, sum((fl[i]*λ_val[i]) for i=1:obj) )
    # objective_function(ap_m, AffExpr)
        # xval = value.(x)
    optimize!(ap_m)
    rval = []
    for i=1:obj
        cost = P2[i,:]
        push!( rval, dot(cost,value.(x)) )
    end
    # print(termination_status(ap_m),rval)
    return rval
end

##########################  Calculate λ  ############################
function calλ(R)
    lam = Model(with_optimizer(CPLEX.Optimizer))
    @variable(lam, λ[1:3])
    @constraint( lam, sum(λ[1:3]) == 1)  #Pkg add https://github.com/VMLS-book/VMLS.jl
    @constraint( lam, dot(λ,nth(R,1)) == dot(λ,nth(R,2)) )
    @constraint( lam, dot(λ,nth(R,2)) == dot(λ,nth(R,3)) )
    # @constraint( lam, dot(λ,R[1]) == dot(λ,R[3]) )
        # @objective( lam, Max, -sum(λ[i]*R[i] for i=1:obj) )
    optimize!(lam)
        # obj_value = JuMP.objective_value(lam)
    return value.(λ)
end


############################  Find a new λ #############################
# function findnewλ(fy,R)
#     lam_m = Model(with_optimizer(CPLEX.Optimizer))
#     @variable(lam_m, nλ[1:obj]>=0 )
#     @constraint( lam_m, sum(nλ[1:obj]) == 1)
#     @constraint( lam_m, cons[i=1:obj] ,dot(nλ,fy) <= dot(nλ,R[i]) )
#     # @objective( lam_m, Min, sum(ex) )
#     optimize!(lam_m)
#     newλ = [value.(nλ)]
#
#     if termination_status(lam_m) == MOI.OPTIMAL
#         return newλ = [value.(nλ)]
#     else
#         print("there's no value ============")
#     end
# end
