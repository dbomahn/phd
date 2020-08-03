using DelimitedFiles, CPLEX, LinearAlgebra, JuMP,TickTock, IterTools,CSV, DataFrames
##########################  Read input files #########################
f=readdlm(ARGS[1], '\t', String, '\n')
# path = "C:\\cygwin64\\home\\AK121396\\multiobjective\\instances\\moolibrary\\AP\\"
# f = readdlm(path*"AP_p-3_n-5_ins-6.dat", '\t', String, '\n')
obj=parse(Int,f[1])
n=parse(Int,f[2])
################## objective function coefficient (P) matrix ################
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
######################### MIP model ####################################
function solveAP(B,P2,lam_val)
    ap_m = Model(with_optimizer(CPLEX.Optimizer))
    # lam_val =[1/3,1/3,1/3]
    # coeff = [P2[i,:]*lam_val[i] for i=1:3]
    @variable( ap_m, 0<=x[1:n*n]<=1 )
    @constraint( ap_m, cons[i=1:2*n], sum(B[i,j]*x[j] for j=1:n*n) == 1 )
    @expression( ap_m, fl[i=1:obj], sum(P2[i,j]*x[j] for j=1:n*n) )
    @objective( ap_m, Min, sum((fl[i]*lam_val[i]) for i=1:obj) )
    # objective_function(ap_m, AffExpr)
        # xval = value.(x)
    optimize!(ap_m)
    rval = []
    for i=1:obj
        cost = P2[i,:]
        push!( rval, dot(cost,value.(x)) )
    end
    # print(termination_status(ap_m),rval)
    return rval,value.(x)
end
#########################  Calculate lam  ############################
function callam(R)
    lamb = Model(with_optimizer(CPLEX.Optimizer))
    @variable(lamb, lam[1:3])
    @constraint( lamb, sum(lam[1:3]) == 1)  #Pkg add https://github.com/VMLS-book/VMLS.jl
    @constraint( lamb, dot(lam,nth(R,1)) == dot(lam,nth(R,2)) )
    @constraint( lamb, dot(lam,nth(R,2)) == dot(lam,nth(R,3)) )
    # @constraint( lam, dot(lam,R[1]) == dot(lam,R[3]) )
        # @objective( lam, Max, -sum(lam[i]*R[i] for i=1:obj) )
    optimize!(lamb)
        # obj_value = JuMP.objective_value(lam)
    return value.(lam)
end
############################ ExA Algorithm ###############################
# Step1
countLP = 0
V = Set()
F = Set()
itval = 20
BigM = (itval*n*obj)*(itval*n)
m1 = Float64[BigM,0,0]
m2 = Float64[0,BigM,0]
m3 = Float64[0,0,BigM]
Y_em = [m1,m2,m3]
Y_e = []
X = []
k=obj
y = Dict((1=>m1,2=>m2,3=>m3))
L = Set([Set([m1,m2,m3])])
# Step2

 #measure CPUtime
start = time()
while isempty(L)==false
    global L,V,Y_em,Y_e,F,y,X,countLP
    R = Set((nth(L,1)))
    V = union(V,Set([R]))
    #Step3 calculate lam (or just use function to get normal)
    lam_val = callam(R)
    #Step4
    if all(x->x>=0, lam_val)==true
        r_st, xval = solveAP(B,P2,lam_val) #Step4.1 (#attach "[1]" for AP instances)
        countLP+=1
        if r_st ∈ R                 #Step4.2
            F = union(Set([F]),Set([R]))
        else                        #Step4.3
            subR = [nth(R,1),nth(R,2),nth(R,3)]
            for i ∈ subsets(subR,2)   #Step4.3.1
                itv = push!(i,r_st)
                sitv = Set(itv)
                L = push!(L, sitv) #Set is not duplicated by using push!
            end
            interLV = intersect(L,V)
            L = setdiff(L,interLV)
            Y_e = Y_em[4:end]
            if r_st ∉ Y_e   #Step4.3.2
                global k = k+1
                y2 = Dict([k => r_st])
                temp = merge!(y,y2)
                y = temp
                Y_em = union(Y_em, [y[k]])
                if xval ∉ X
                    push!(X,xval)
                end
            end
        end
    else                    #Step5:  any lam_val < 0
        for i in k:-1:1
            if (y[i]∈Y_em && y[i]∉ R)
                r_st = y[i]
                break
            end
        end

        subR = [nth(R,1),nth(R,2),nth(R,3)]
        for i ∈ subsets(subR,2)   #Go to Step4.3 Again
            itv = push!(i,r_st)
            sitv = Set(itv)
            L = push!(L, sitv) #Set is not duplicated by using push!
        end
        interLV = intersect(L,V)
        L = setdiff(L,interLV)
        Y_e = Y_em[4:end]
        if r_st ∉ Y_e
            k = k+1
            y2 = Dict([k => r_st])
            temp = merge!(y,y2)
            y = temp
            Y_em = union(Y_em, [y[k]])
            if xval ∉ X
                push!(X,xval)
            end
        end
    end
    L = setdiff(L,Set([R]))  #Step6
    # print("\n length of interval array: ", length(L))
end
elapsed = time() - start
timeLP = DataFrame( CPUtime=Float64(elapsed), solvedLP=countLP )
solvar = DataFrame( sol=Y_em[4:end], var=X )
CSV.write("Dicho_AP_time_LP.csv", timeLP, append=true, writeheader=true, delim=';' )
CSV.write("Dicho_AP_sol_var.csv", solvar,  append=true, writeheader=true, delim=';' )
# open(ARGS[1]*".log", "a") do io
#     # writedlm(io, "CPU time,solved LP,var_val,PF\n")
#     writedlm(io, [Float64(elapsed),countLP, Y_e, X])
# end

# io = open(ARGS[1]*".log", "w")
# write(io, "CPU time: ",string(elapsed),"\n")
# write(io, "# of solved LP: ", countLP,"\n")
# write(io, "varialb values: \n", X, "\n")
# write(io,"Pareto Frontier: ", string(Y_e))
# close(io)
# print("CPU time: ", string(elapsed),"\n")
# print("# of solved LP: ", countLP,"\n")
# print("varialb values: \n", string(X), "\n")
# print("Pareto Frontier: ", string(Y_e))
