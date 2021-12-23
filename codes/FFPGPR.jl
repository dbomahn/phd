using DataStructures,DelimitedFiles,DataFrames,JuMP,CPLEX,GLPK,LinearAlgebra,CPUTime,CSV,Statistics,StatsBase,MathProgBase,MathOptInterface,MathOptFormat,JLD2,SparseArrays
const MPB = MathProgBase; const MOF = MathOptFormat; const MOI = MathOptInterface;
@show ARGS[1]; @show ARGS[2]; @show ARGS[3]; @show ARGS[4]; @show ARGS[5];
TL = 3600
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
    x::String; y::String; dvar::Array{}; LB::Array{}; LBmtx::Array{};
    function Valu(x,y)
        JLD2.@load x dv;
        dv0 = Array(dv);
        dv1 = round.(dv0; digits=4);
        objs = round.(digits=4, readdlm(y));
        ind = findall(i-> 0 in objs[i,:]  , 1:size(objs)[1]);
        dv2 = dv1[setdiff(1:end, ind), :];
        LBmtx = objs[setdiff(1:end, ind), 2:end];
        dvar = [Vector(dv2[i, :]) for i = 1:size(LBmtx)[1]];
        LB = [Vector(LBmtx[i, :]) for i = 1:size(LBmtx)[1]];
        new(x,y,dvar,LB,LBmtx)
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
        if all( x .>= P[k])# && any(x > P[k])
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
function Postpro(initX,initY,dvar,LB,newsol)
    P = union(initX,dvar[end-newsol+1:end]); Pobj = union(initY,LB[end-newsol+1:end])
    #Filter dominated solutions
    copysol = Dict(); copyobj = Dict();
    for i=1:length(Pobj)
        copysol[i] = P[i]; copyobj[i] = Pobj[i]
    end
    for i=1:length(Pobj)-1
        for j=i+1:length(Pobj)
            if all(Pobj[i] .>= Pobj[j]) == true #dominated by PF[j]
                copyobj[i]=nothing; copysol[i]=nothing; break
            elseif all(Pobj[j] .>= Pobj[i]) == true
                copyobj[j]=nothing; copysol[j]=nothing;
            end
        end
    end
    finalsol = filter!(a->a!=nothing, collect(values(copysol)))
    finalobj = filter!(a->a!=nothing, collect(values(copyobj)))
    return finalsol,finalobj
end
#############################    Feasibility Pump   ############################
function fractionalFP(dvar,n,C,TL)
    X = copy(dvar); IGPair=[]; Tabu = []; newsol=0; t0=time(); #Y = copy(LB);
    while time()-t0 < TL
        I,G = sample(1:length(X), 2, replace=false)
        x1 = X[I]; x2 = X[G];
        λ = round(rand(Float64, 1)[1]; digits=1)
        x_t = x1*λ + x2*(1-λ); SearchDone = false; iter=0; Max_iter = 10
        while [I,G]∉IGPair && iter<Max_iter && SearchDone == false
            x_r = round.(Int,x_t)
            if ( (FBcheck(x_r,n) == true) && (x_r∉X) )
                push!(X,x_r); #Y = [Y; getobjval(x_r,C)'];
                newsol+=1; SearchDone = true
            else
                if x_r ∈ Tabu
                    x_r = flipoper(Tabu,x_t,x_r)
                    if x_r==[]
                        SearchDone = true
                    else
                        if ( (FBcheck(x_r,n) == true) && (x_r∉X) )
                            push!(X,x_r); #Y = [Y; getobjval(x_r,C)'];
                            newsol+=1; SearchDone = true;
                        end
                    end
                end
                if time()-t0 >= TL
                    break
                end
                if SearchDone == false
                    push!(Tabu,x_r)
                    x_t = fbsearch(x_r)
                    if x_t == 0 #when there's no new feasible lp sol
                        SearchDone = true
                    end
                end
            end
			iter+=1
        end
        push!(IGPair,[I,G])
    end
    return X,newsol
end
################################  Data  ####################################
data = Data(ARGS[1],ARGS[2]); pre = Valu(ARGS[3],ARGS[4]);
######################### Mathematical Model #############################
mip = Model(GLPK.Optimizer)
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
optimize!(mip);
#####################  Feasibility Search Mathematical Model  ##################
dist = Model(GLPK.Optimizer)
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
optimize!(dist);
######################## Running Feasibility Pump ##########################
Bentime = readdlm(ARGS[5])[1]; FPTL = (TL-Bentime)/2;
weightFP(pre.dvar,data.n,data.C,5)
runtime1 = @CPUelapsed fpX,newsol = fractionalFP(pre.dvar,data.n,data.C,FPTL)

# fpX,fpY = domFilter(X2,PF2);
# clistY = []
# for i=1:length(candlist)
# 	push!(clistY,getobjval(candlist[1],data.C))
# end
# cX = [fpX;candlist]; cY = [fpY;clistY]
tx = copy(fpX); ty = copy(fpY);

####################### Running Path Relinking  ##########################
function createNB(SI,C,dif,exploredSI)
    neibour = []; neiobj = [];
    for i in dif
        cpSI = copy(SI)
        if cpSI[i] == round(cpSI[i]) #if vari is int value
            if cpSI[i] == 1
                cpSI[i] = 0
            else
                cpSI[i] = 1
            end
            push!(neibour, cpSI);
            push!(neiobj, getobjval(cpSI,C))
        else #if vari is fractional
            cpSI[i] = 1; push!(neibour, cpSI); push!( neiobj, getobjval(cpSI,C) )
            cpSI = copy(SI)
            cpSI[i] = 0; push!(neibour, cpSI); push!( neiobj, getobjval(cpSI,C) )
        end
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

function GPR(candX,candY,C,n,TL)
    IGPair=[]; exploredSI = []; newsol=0; t0=time();
    while time()-t0 < TL
	    I,G = sample(1:length(candX), 2, replace=false)
        SI = candX[I]; SG = candX[G]; infeas=0; #Maxiter = length(candX)*10;  #if the algorithm finds infeasible solutions in a row, move to the next iteration
        while all.(SI != SG) && [I,G]∉IGPair && iter<20 && (time()-t0<TL)
            dif = findall(i-> SI[i]!=SG[i], 1:n)
            neibour,neiobj = createNB(SI,C,dif,exploredSI)
            if ( (length(neiobj)==0) || (time()-t0 >= TL) )
                break
            else
                for l=1:length(neiobj)
                    if ( (FBcheck(neibour[l],n)==true) && neibour[l]∉candX ) #dominance check (dominated(neiobj[l],LB)==false)
                        push!(candX, neibour[l]); push!(candY, neiobj[l]);
                        newsol+=1;
                    end
                end
            end
            SI = nextSI(neibour,neiobj,C,SI)
            if SI∉candX
                push!(exploredSI,SI);
            end
            iter+=1
        end
        push!(IGPair,[I,G])
    end
    return candX,candY,newsol
end

GPR(data.C,data.n,tx,ty,5) #compiling
runtime1 = @CPUelapsed X,nsol = fractionalFP(pre.dvar,data.n,data.C,TL)
runtime2 = @CPUelapsed candset,candobj,newsol2 = GPR(data.C,data.n,fpX,fpY,TL-FPtime)
totaltime = runtime1+Bentime+runtime2
finalX,finalY = domFilter(candset,candobj);
otable = zeros(length(finalY),3)
for i=1:length(finalY)
    for j=1:3
        otable[i,j] = finalY[i][j]
    end
end
newsol3 = length(setdiff(finalY,fpY))

ins = ARGS[2][1:end-4];
print("$ins", totaltime)
CSV.write(ins*"_1Y.log",DataFrame(otable, :auto),append=false, header=false, delim=' ' )
record1 = DataFrame(Filename=ins[48:end],FPsol=length(fpY),PRsol=newsol2, totalsol=length(finalY)) #,Bentime=Bentime,FPtime=FPtime,PRtime=GPRtime, totaltime=totaltime )
CSV.write("/home/k2g00/k2g3475/clusterhome/multiobjective/generalPR/record.csv", record1,append=true, header=false )#, delim=','i )
# matriX = zeros(Int,length(finalX),data.n)
# for i=1:length(finalX)
#     for j=1:data.n
#         matriX[i,j] = finalX[i][j]
#     end
# end
function FPGPR(dvar,LB,n,C,TL)
    candX = copy(dvar); candY = copy(LB); IGPair=[];  newsol=0;
	gMax = 20; giter=0;	exploredSI = []; Tabu = []; t0=time();
    while time()-t0 < TL #&& newsol<=length(candY)*0.5
        I,G = sample(1:length(candY), 2, replace=false) #BASIC VER
		SI = candX[I]; SG = candX[G];
        while all.(SI != SG) && [I,G]∉IGPair && giter<gMax
			dif = findall(i-> SI[i]!=SG[i], 1:n)
			neibour,neiobj = createNB(SI,C,dif,exploredSI)
			if length(neibour)==0
				break
			else
				for l=1:length(neibour)
					SearchDone = false; fiter=0; fMax = 10;
					while fiter<fMax && SearchDone == false && time()-t0 <TL
						x_r = round.(Int,neibour[l])
						fx = getobjval(x_r,C)
			            if ( (FBcheck(x_r,n) == true) && (fx ∉ candY) ) #checking feasibility and dominance   #z[1]<=fmax[1]-1 && z[2]<=fmax[2]-1 && z[3]<=fmax[3]-1)
			                push!(candX,x_r); push!(candY,fx)#Y = [Y; getobjval(x_r,C)'];
			                newsol+=1; SearchDone = true
			            else
			                if x_r ∈ Tabu
			                    x_r = flipoper(Tabu,neibour[l],x_r)
			                    if x_r==[]
			                        SearchDone = true
			                    else
			                        fx = getobjval(x_r,C)
			                        if ( (FBcheck(x_r,n) == true) && (fx ∉ candY) )  #z[1]<=fmax[1]-1 && z[2]<=fmax[2]-1 && z[3]<=fmax[3]-1 )
			                            push!(candX,x_r); push!(candY,fx)#Y = [Y; getobjval(x_r,C)'];
			                            newsol+=1; SearchDone = true;# push!(sucess,[I,G])
			                        end
			                    end
			                end
			                if SearchDone == false
			                    push!(Tabu,x_r)
			                    x_t = fbsearch(x_r)
			                    if x_t == 0 #when there's no new feasible lp sol
			                        SearchDone = true
			                    end
			                end
			            end
						fiter+=1
        			end
				end
			end
			SI = nextSI(neibour,neiobj,C,SI)
			if SI∉candX
				push!(exploredSI,SI);
			end
			giter+=1
		end
		push!(IGPair,[I,G])
    end

	while time()-t0 < TL
		I,G = sample(1:length(candY), 2, replace=false) #BASIC VER
		SI = candX[I]; SG = candX[G]; giter=0

        while all.(SI != SG) && [I,G]∉IGPair && giter<gMax
			dif = findall(i-> SI[i]!=SG[i], 1:n)
			neibour,neiobj = createNB(SI,C,dif,exploredSI)
			if length(neibour)==0
				break
			else
				for l=1:length(neibour)
					if FBcheck(neibour[l],n) && neibour[l]∉ candX
						push!(candX, neibour[l]); candY = [candY; neiobj[l]']; #push!(candY, neiobj[l]);
						newsol+=1; gn+=1;
					end
				end
			end
			SI = nextSI(neibour,neiobj,C,SI)
			if SI∉candX
				push!(exploredSI,SI);
			end
			giter+=1
		end
		push!(IGPair,[I,G])
	end

    return candX,candY,newsol
end
