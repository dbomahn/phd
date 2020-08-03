using DelimitedFiles, CPLEX, LinearAlgebra, JuMP,TickTock, IterTools
# using GLPK (faster than CPLEX)
##########################  Read input files #########################
path = "C:\\cygwin64\\home\\AK121396\\multiobjective\\instances\\moolibrary\\ILP\\"
f = readdlm(path*"ILP_p-3_n-10_m-5_ins-1.dat", '\t', String, '\n')

global obj=parse(Int, f[1])
global n=parse(Int,f[2])
global m=parse(Int,f[3])
########################    MIP instances    #####################
##################objective function coefficient (P) matrix ################
c = f[4:obj+3]
P= ones(obj,n)
P=round.(Int,P)

global zp=0;
for x=1:length(c)
    a = split(c[x],('[',']',','))
    aa = filter!(e->e ∉ [ "" ,"[","," ,"]"] ,a)
    for y = 1:length(aa)
        p = parse(Int64,aa[y])
        P[x,y] = p
        if p==0
            global zp = zp+1;
        end
    end
end
####### technical coefficient (a) ###########
TC = f[obj+4:length(f)-1]
global za=0;
a= ones(m,n)
a=round.(Int,a)
for x=1:length(TC)
    t = split(TC[x],('[',']',','))
    tt = filter!(e->e ∉ [ "" ,"[","," ,"]"] ,t)
    for y=1:length(tt)
        tc=parse(Int64,tt[y])
        a[x,y] = tc
        if tc==0
            global za = za+1;
        end
    end
end
####### RHS values (b) ##########
b = ones(1,m)
b = round.(Int,b)
r = f[length(f)]
r1 = split(r, ('[',']',','))
r2 = filter!(e->e ∉ ["" ,"[", "]"] ,r1)
for i=1:m
    b[i] = parse(Int64,r2[i])
end

################################       MIP MODEL   ###############################
function solveMIP(a,P,b,lam_val)
    mip = Model(with_optimizer(CPLEX.Optimizer))
    @variable(mip, x[1:n]>=0, Int)
    @constraint(mip, con[i=1:n], sum(a[i][j]*x[j] for j=1:m) <= b[i])
    @expression(mip, ex[i=1:obj], sum(P[i][j]*x[j] for j=1:m))
    @objective(mip, Max, sum(ex[i]*lam_val[i] for i=1:obj))

    optimize!(mip)
    rval = Float64[]
    for i=1:obj
        C = a[((i-1)*n)+1:i*n]
        push!( rval, dot(C,value.(x)) )
    end
    return rval
end
##############################################################################
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
V = Set()
F = Set()
maxUB = 100
BigM = (maxUB*n*obj)*(maxUB*n)
m1 = Float64[BigM,0,0]
m2 = Float64[0,BigM,0]
m3 = Float64[0,0,BigM]
Y_em = [m1,m2,m3]
Y_e = []
k=obj
y = Dict((1=>m1,2=>m2,3=>m3))
L = Set([Set([m1,m2,m3])])
# Step2
# count =1
 #measure CPUtime

start = time()
while isempty(L)==false
# while count<10
    # print("\n",count,"th Loop " , "\n")
    global L,V,Y_em,Y_e,F,y
    R = Set((nth(L,1)))
    # print("\n ======show R======== ", R)
    V = union(V,Set([R]))
    #Step3 calculate lam (or just use function to get normal)
    lam_val = callam(R)
    # print("\n lam value =", lam_val)
    #Step4
    if all(x->x>=0, lam_val)==true
        r_st = solveMIP(a,P,b,lam_val) #Step4.1 (#attach "[1]" for AP instances)
        if r_st ∈ R                 #Step4.2
            F = union(Set([F]),Set([R]))
            # print("\n positive case, r_st ∈ R  " )
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
                # print("\n positive case:: r_st ∈ R but  r_st ∉ Y_e ")
                # print("\n positive case:: y[k]= ",k, y[k])

            # else
                # print("\n positive case::::r_st is already in Y_e : ", r_st)
            end
        end
    else                    #Step5:  any lam_val < 0
        for i in k:-1:1
            if (y[i]∈Y_em && y[i]∉ R)
                r_st = y[i]
                # print("\n Here is fy from y[k] :\n", fy)
                break
            end
        end
        # print("\n negative case:: r_st",r_st)
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
            # print("\n negative case:: new r_st: ",r_st)
        end
    end
    L = setdiff(L,Set([R]))  #Step6
    print("\n length of interval array: ", length(L))
    # global count=count+1
end
elapsed = time() - start
# Step7
Y_e = Y_em[4:end] #remove dummy points from Y_em
# return Y_e
print("\n computational time: ", elapsed)
print("\n Y_e is : " ,Y_e, "\n length is  ", length(Y_e))

# io = open(ARGS[1]*".log", "w")
# write(io,string(elapsed),"\n")
# write(io,"Pareto Frontier: ", string(Y_e))
# close(io)


# for i=1:length(Y_e)
#     println(io,Y_e[i][1]," ",Y_e[i][2]," ",Y_e[i][3])
# end
