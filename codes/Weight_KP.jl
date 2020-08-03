using Revise, DelimitedFiles, GLPK,GLPKMathProgInterface, MathProgBase, LinearAlgebra, JuMP,CPUTime, IterTools, CSV, DataFrames
##########################  Read input files #########################
# fl=readdlm(ARGS[1], '\t', String, '\n')
# path = "C:\\Users\\AK121396\\Desktop\\Dropbox\\JKU\\Julia\\KP\\"
path = pwd()*"/Dropbox/JKU/Julia/KP/"
# path ="/home/ak121396/multiobjective/instances/moolibrary/KP"
f = readdir(path)
# f =readdlm(path*"KP_p-3_n-30_ins-"*".dat", '\t', String, '\n')
for i=1:10
    fl = readdlm(path*f[i], '\t', String, '\n')

    obj=parse(Int,fl[1])
    n=parse(Int,fl[2])
    ub=parse(Int,fl[3])
    ########################    KP     #####################
    b = fl[4:length(fl)-1]
    C= ones(length(b),n)
    C = round.(Int,C)
    global ct=0;
    for x=1:length(b)
        a = split(b[x],('[',']',','))
        aa = filter!(e->!(e in [ "" ,"[","," ,"]"]) ,a)
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
    w2 = filter!(e->!(e in ["" ,"[", "]"]) ,w1)
    for i=1:n
        weight[i] = parse(Int64,w2[i])
    end
    ########################## MIP model ####################################
    function solveKP(weight,C,ub,lam_val)
        kp_m = Model(with_optimizer(GLPK.Optimizer))#,CPX_PARAM_THREADS=1,CPX_PARAM_SCRIND=0))
        @variable(kp_m, 0<=x[1:n]<=1 )
        @constraint( kp_m, sum(weight[i]*x[i] for i=1:n) <= ub )
        @expression( kp_m, ex[i=1:obj], dot(C[i,:],x[:]) )
        @objective( kp_m, Max, sum(ex[i]*lam_val[i] for i=1:obj) )
        optimize!(kp_m)
        time = MOI.get(kp_m, MOI.SolveTime())
        rval = Float64[]
        for i=1:obj
            # P = C[i,:]
            # push!( rval, dot(P,value.(x)) )
            push!(rval, round(dot(C[i,:],value.(x))))
        end
        return rval,value.(x),time
    end

    function callam(R)
        lamb = Model(with_optimizer(GLPK.Optimizer)) #,CPX_PARAM_THREADS=1,CPX_PARAM_SCRIND=0))
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

    #small example for precompiling
    solveKP(weight,C,ub,[1/3,1/3,1/3])

    ############################ ExA Algorithm ###############################
    # Step1
    global countLP = 0
    global V = Set()
    global F = Set()
    global itval = 1000
    global BigM = (itval*n*obj)*(itval*n)
    global m1 = Float64[BigM,0,0]
    global m2 = Float64[0,BigM,0]
    global m3 = Float64[0,0,BigM]
    global Y_em = [m1,m2,m3]
    global X = []
    global Y_e = []
    global k=obj
    global y = Dict((1=>m1,2=>m2,3=>m3))
    global L = Set([Set([m1,m2,m3])])
    global MIPtime = []
    global ItvTime = []
    # Step2

    #measure CPUtime
    CPUtic()
    while isempty(L)==false
        global L,V,Y_em,Y_e,F,y,X,countLP,MIPtime,ItvTime
        R = Set((nth(L,1)))
        V = union(V,Set([R]))
        #Step3 calculate lam (or just use function to get normal)
        lam_val = callam(R)
        #Step4
        ############Record time of making intervals########
        CPUtic()
        if all(x->x>=0, lam_val)==true
            global r_st, xval,Ptime = solveKP(weight,C,ub,lam_val) #Step4.1 (#attach "[1]" for AP instances)
            countLP+=1
            push!(MIPtime,Ptime)    #Record solving time of MIPs

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
        L = setdiff(L,Set([R]))  #Step6
        # Record total time of making Intervals
        push!(ItvTime,CPUtoq())
    end
    fin = CPUtoq()

    # make output.csv file to the folder to ./home
    timeLP = DataFrame( CPUtime=fin, solvedLP=countLP)
    solvar = DataFrame( sol=Y_em[4:end], var=X )
    CSV.write("WD_KP_n10_time_#LP.csv", timeLP; append=true, writeheader=false, delim=';' )
    # CSV.write("WD_KP_n10_sol_var_LP.csv", solvar;  append=true, writeheader=true, delim=';' )
    MIPTime = DataFrame(mip_time = MIPtime)
    IntervalTime = DataFrame(interval_time = ItvTime)
    # CSV.write("WD_KP_n10_Interval_Time.csv", IntervalTime; append=true, writeheader=true)
    # CSV.write("WD_KP_n10_MIP_Time.csv", MIPTime; append=true, writeheader=true)
end
print("done")
