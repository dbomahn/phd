using DataStructures,DelimitedFiles,DataFrames,JuMP,CPLEX,LinearAlgebra,CPUTime,CSV,Statistics,MathProgBase,MathOptInterface,MathOptFormat,JLD2,SparseArrays,StatsBase
const MPB = MathProgBase; const MOF = MathOptFormat; const MOI = MathOptInterface;
mutable struct Data
    mpsfile::String; vlpfile::String
    m::Int; n::Int; C::Array{}; B::Array{}; RHS::Dict{}; signs::Array{};
    function Data(mpsfile::String, vlpfile::String)
        mpsmodel = CPLEX.CplexMathProgModel(); MPB.loadproblem!(mpsmodel,mpsfile)
        B = MPB.getconstrmatrix(mpsmodel); m,n=size(B)
        lb = MPB.getconstrLB(mpsmodel);ub = MPB.getconstrUB(mpsmodel);
        RHS = Dict()
        for i=1:m
            if ub[i]==Inf
                RHS[i] = lb[i]
            else
                RHS[i] = ub[i]
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
        vlp = readdlm(vlpfile); os = findall(i->vlp[i,1]=="o", 1:length(vlp[:,1]))
        inio = os[1]; endo = os[end]; Pmtx = vlp[inio:endo,2:4]; C = zeros(3,n)
        for k=1:length(Pmtx[:,1]) #(len,id) enumerate
            x = Pmtx[k,:][1]; y = Pmtx[k,:][2];
            C[x,y] = Pmtx[k,:][3]
        end
        new(mpsfile,vlpfile,m,n,C,B,RHS,signs)
    end
end
mutable struct Valu
    x::String; y::String; dvar::Array{}; LB::Array{}; LBmtx::Array{}; PF::Array{}
    function Valu(x,y)
        JLD2.@load x dv
        dv1 = round.(dv; digits=4)
        objs = round.(digits=4, readdlm(y))
        ind = findall(i-> 0 in objs[i,:]  , 1:size(objs)[1])
        dv2 = dv1[setdiff(1:end, ind), :];
        LBmtx = objs[setdiff(1:end, ind), 2:end];
        dvar = [Vector(dv2[i, :]) for i = 1:size(LBmtx)[1]]
        LB = [Vector(LBmtx[i, :]) for i = 1:size(LBmtx)[1]]
        PF = [];
        for i=1:length(LB)
            if count(x->0<x<1,dvar[i])==0 #if all variable values are binary
                push!(PF,LB[i])
            end
        end
        new(x,y,dvar,LB,LBmtx,PF)
    end
end
################################   Functions  ##################################
function flip(x_h,j,e)
    if x_h[e[j]]==1
        x_h[e[j]] = 0
    else
        x_h[e[j]] = 1
    end
    return x_h
end
function flipoper(Tabu,x_t,x_r)
    e = sortperm(abs.(x_t-x_r),rev=true)
    xi = []
    x_h = copy(x_r)
    j = 1
    M=length(x_t) #
    while j<=M && xi==[]
        x_h = flip(x_h,j,e)
        if x_h ∉ Tabu
            xi=x_h
        else
            j+=1
        end
    end
    if xi==[]
        while j<=M
            x_h=copy(x_r)
            Num = Int64(rand(ceil(length(x_r)/2):length(x_r)-1))
            R = sample(1:M,Num, replace=false)
            for i in R
                x_h = flip(x_h,r)
                if x_h ∉ Tabu
                    xi = x_h
                end
            end
            j+=1
        end
    end
    return xi
end
function FBcheck(xx,n)
    for k=1:n
        JuMP.fix(mx[k],xx[k]; force=true)
    end
    optimize!(mip)
    if termination_status(mip) == MOI.OPTIMAL
        return true
    else
        return false
    end
end
function fbsearch(x_r) #solveLP
    idx0 = findall(k->k==0, x_r)
    idx1 = findall(k->k==1, x_r)
    @objective( dist, Min, sum(x[i] for i in idx0) + sum(1-(x[j]) for j in idx1) )
    optimize!(dist)
    if termination_status(dist) == MOI.OPTIMAL
        return JuMP.value.(x)
    else
        return 0;
    end
end
function getobjval(x,C)
    return [dot(x,C[1,:]),dot(x,C[2,:]),dot(x,C[3,:])]
end
function dominated(x,P)
    st = false
    for k=1:length(P)
        if all( x .>= P[k]) && any(x > P[k])
            st=true; break
        else
            continue
        end
    end

    return st
end
function domFilter(sol,obj)
    copysol = Dict(); copyobj = Dict();
    for i=1:length(obj)
        copysol[i] = sol[i]
        copyobj[i] = obj[i]
    end
    for i=1:length(obj)-1
        for j=i+1:length(obj)
            if all(obj[i] .>= obj[j]) == true #dominated by PF[j]
                copyobj[i]=nothing; copysol[i]=nothing;break
            elseif all(obj[j] .>= obj[i]) == true
                copyobj[j]=nothing; copysol[j]=nothing;
            end
        end
    end
    finalsol = filter!(a->a!=nothing, collect(values(copysol)))
    finalobj = filter!(a->a!=nothing, collect(values(copyobj)))
    return finalsol,finalobj
end

############################ Feasibility Pump ver2  ##########################


#############################    Feasibility Pump ver1  ############################
function FP1(candX,n,C,TL)
    X= []; Y = []; Tabu = []; newsol=0; LPcount = 0; t0=time();
    for k=1:length(candX)
        if time()-t0 >= TL
            break
        end
        # print("========= Feasi Pump",k," th candidate sol =========\n")
        x_t = candX[k];
        SearchDone = false
        itr = 1
        Max_itr = 10 #max(count(x->0<x<1,x_t),1)
        while itr<Max_itr && SearchDone==false
            x_r = round.(Int,x_t)
            fx = getobjval(x_r,C)
            if ( (FBcheck(x_r,n) == true) && (dominated(fx,Y)==false) ) #checking feasibility and dominance   #z[1]<=fmax[1]-1 && z[2]<=fmax[2]-1 && z[3]<=fmax[3]-1)
                push!(X,x_r); push!(Y,fx) #add new solval to Y
                newsol+=1; #print("ROUNING => newsol added \n")
                SearchDone = true
            else
                if x_r ∈ Tabu
                    x_r = flipoper(Tabu,x_t,x_r)
                    if x_r==[]
		    	#print("FLIP didn't work \n")
                        SearchDone = true
                    else
                        fx = getobjval(x_r,C)
                        if ( (FBcheck(x_r,n) == true) && (dominated(fx,Y)==false) )  #z[1]<=fmax[1]-1 && z[2]<=fmax[2]-1 && z[3]<=fmax[3]-1 )
                            push!(X,x_r); push!(Y,fx) #add new solval to Y
                            newsol+=1; #print("FLIP=> newsol added \n")
                            SearchDone = true
                        end
                    end
                end
                if time()-t0 >= TL
                    break
                end
                if SearchDone == false
                    push!(Tabu,x_r) #when break
                    x_t = fbsearch(x_r)
                    LPcount+=1
                    if x_t == 0 #when there's no new feasible lp sol
                        # print("no lp sol's found");
                        break
                    end
                end
            end
            itr+=1;
        end
    end
    return X,Y,newsol,LPcount
end
################################  Data  ####################################
data = Data("/home/ak121396/Desktop/instances/MIPLIB(official)/mps/eil33-2.mps","/home/ak121396/Desktop/instances/MIPLIB(official)/vlp/eil33-2.vlp")
pr = Valu("/home/ak121396/Desktop/solvers/Bensolve/MIPLIB/official/benout//X/irp_X.jld2","/home/ak121396/Desktop/solvers/Bensolve/MIPLIB/official/benout/Y/irp_img_p.sol")

######################### Mathematical Model #############################
mip = Model(with_optimizer(CPLEX.Optimizer))
MOI.set(mip, MOI.RawParameter("CPX_PARAM_SCRIND"), false)
MOI.set(mip, MOI.RawParameter("CPX_PARAM_THREADS"),1  )
@variable(mip, mx[1:data.n], Bin)
for k=1:data.m
    if data.signs[k] == "l"
        @constraint(mip, dot(data.B[k,:],mx) >= data.RHS[k])
    elseif data.signs[k] == "u"
        @constraint(mip, dot(data.B[k,:],mx) <= data.RHS[k])
    else
        @constraint(mip, dot(data.B[k,:],mx) == data.RHS[k])
    end
end
optimize!(mip)
#####################  Feasibility Search Mathematical Model  ##################
dist = Model(with_optimizer(CPLEX.Optimizer))
MOI.set(dist, MOI.RawParameter("CPX_PARAM_SCRIND"), false)
MOI.set(dist, MOI.RawParameter("CPX_PARAM_THREADS"),1  )
@variables(dist, begin
    0 <= x[1:data.n] <=1
end)
for k=1:data.m
    if data.signs[k] == "l"
        @constraint(dist, dot(data.B[k,:],x) >= data.RHS[k])
    elseif data.signs[k] == "u"
        @constraint(dist, dot(data.B[k,:],x) <= data.RHS[k])
    else
        @constraint(dist, dot(data.B[k,:],x) == data.RHS[k])
    end
end
optimize!(dist)
#Time Feasibility Pump
FPtime = @CPUelapsed X2,PF2,newsol,Tabu,LPcount = FP(pr.dvar,pr.PF,data.n,data.C,5)
#Filter dominated solutions
fpX, fpY = domFilter(X2,PF2)
####################### Running Path Relinking  ##########################
function createNB(SI,C,dif,exploredSI)
    neibour = []; neiobj = [];
    for i in dif
        cpSI = copy(SI)
        if cpSI[i] == 1
            cpSI[i] = 0
        else
            cpSI[i] = 1
        end
        push!(neibour, cpSI);
        push!(neiobj, getobjval(cpSI,C))
    end
    idx = findall(i-> neibour[i] in exploredSI, 1:length(neibour))
    neibour = setdiff(neibour,exploredSI)
    neiobj = neiobj[setdiff(1:end, idx),:]

    return neibour,neiobj
end
function nextSI(neibour,neiobj,C,SI)
    SIobj = getobjval(SI,C)
    for i=1:length(neiobj)
        if length(neiobj) == 1  #if there is one candiate sol
            return neibour[1]#,neiobj[1]
        else length(neiobj) > 1 # if there are multiple candiates, check the improved ratio
            ratiotb = zeros(length(neiobj),3)
            for i=1:length(neiobj)
                ratiotb[i,:] = neiobj[i]./SIobj
            end
            ranktb = zeros(length(neiobj),3)
            for i=1:3
                ranktb[:,i] = tiedrank(ratiotb[:,i])
            end
            ranksum = [sum(ranktb[i,:]) for i=1:length(neiobj)]
            mostimp = findall(x-> x == maximum(ranksum), ranksum)
            k = rand(mostimp)[1]
            return neibour[k]#, neiobj[k]
        end
    end
end
function GPR(C,n,dvar,LB,TL)
    IGPair=[]; exploredSI = []; newsol=0; t0=time();
    for i=1:n*10
    	if time()-t0 >= TL
            break
        end
        I,G = sample(1:length(dvar), 2, replace=false)
        SI = dvar[I]; SG = dvar[G]; Maxiter = length(dvar)*10; infeas=0; #if the algorithm finds infeasible solutions in a row, move to the next iteration
        while all.(SI != SG) && [I,G]∉IGPair && infeas<Maxiter && (time()-t0<TL)
            dif = findall(i-> SI[i]!=SG[i], 1:n)
            neibour,neiobj = createNB(SI,C,dif,exploredSI)
            if length(neiobj)==0
                break
            else
                for l=1:length(neiobj)
                    feasi = FBcheck(neibour[l],n)
                    if feasi == true && neibour[l]∉ dvar
                        push!(dvar, neibour[l]); push!(LB, neiobj[l]);
                        newsol+=1; print("newsol added \n")
                    else
                        infeas+=1; print("infeasible \n")
                    end
                end
                SI = nextSI(neibour,neiobj,C,SI)
                if SI∉dvar
                    push!(exploredSI,SI);
                end
            end
        end
        push!(IGPair,[I,G])
    end
    return dvar,LB,newsol
end
GPR(data.C,data.n,fpX,fpY,10)
GPRtime = @CPUelapsed candset,candobj,newsol2 = GPR(data.C,data.n,fpX,fpY,3600-FPtime-Bentime)

#########################  Record outputs  ############################
ins = ARGS[1][end-12:end-7] #CPUtime recorded, naming after the instance
colname = ARGS[1][end-12:end-4]
record1 = DataFrame( Initial_IPsol = length(Xf1), Initial_LPsol = length(candX),LP=FPLPcount, sol=length(Pz), removed=length(P)-length(Pz))
insertcols!(record1,6, Symbol("$colname")=>GroupingTime+FPTime)
GFPX=DataFrame(matP); GFPY=DataFrame(dfP);

CSV.write("/home/k2g00/k2g3475/multiobjective/solvers/ep+FP/GFPresults/GFPlex_"*"$ins"*"_record.csv",record1, append=true, writeheader=false )#, delim=',' )
#CSV.write(ARGS[1]*"_GFP_X.csv",GFPX, append=true, writeheader=false)
CSV.write(ARGS[1]*"_GFP_Y_.csv",GFPY, writeheader=false, delim=' ' )
print(colname," GroupFP Done!")
