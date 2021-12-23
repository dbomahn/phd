using DataStructures,DelimitedFiles,DataFrames,JuMP,CPLEX,GLPK,LinearAlgebra,CPUTime,CSV,StatsBase,MathProgBase,MathOptInterface,MathOptFormat,JLD2,SparseArrays,Clustering
const MPB = MathProgBase; const MOF = MathOptFormat; const MOI = MathOptInterface;
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
mutable struct Val
    x::String; y::String; dvar::Array{}; LB::Array{}; L::Array{}; k::Int
    function Val(x,y)
        JLD2.@load x dv;
        dv0 = Array(dv);
        dv1 = round.(dv0; digits=4);
        objs = round.(digits=4, readdlm(y));
        ind = findall(i-> 0 in objs[i,:]  , 1:size(objs)[1]);
        dv2 = dv1[setdiff(1:end, ind), :];
        L = objs[setdiff(1:end, ind), 2:end];
        dvar = [Vector(dv2[i, :]) for i = 1:size(L)[1]];
        LB = [Vector(L[i, :]) for i = 1:size(L)[1]];
        if length(LB)<1000
            k = 20
        elseif length(LB)<3000
            k = 60
        else
            k = 200
        end
        new(x,y,dvar,LB,L,k)
    end
end
function getobjval(x,C)
    return [dot(x,C[1,:]),dot(x,C[2,:]),dot(x,C[3,:])]
end
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
    @objective( dist, Min, sum(dx[i] for i in idx0) + sum(1-(dx[j]) for j in idx1) )
    optimize!(dist)
    if termination_status(dist) == MOI.OPTIMAL
        return JuMP.value.(dx)
    else
        return 0;
    end
end
function dominated(x,P)
    st = false
    for k=1:length(P)
        if all( x .>= P[k]) #&& any(x > P[k])
            st=true; break
        else
            continue
        end
    end

    return st
end
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
function Postpro(candX,candY,newsol)
    #Filter fractional solutions from LB
    initdv = candX[1:end-newsol]
    Y = [Vector(candY[i,:]) for i=1:size(candY)[1]]
    initLB = Y[1:end-newsol]
    frac = findall(j-> trunc.(initdv[j])!= initdv[j], 1:length(initdv))
    dv2 = initdv[setdiff(1:end,frac)]; LB2 = initLB[setdiff(1:end,frac)]
    P = union(dv2,candX[end-newsol+1:end])
    Pobj = union(LB2,Y[end-newsol+1:end])

    #Filter dominated solutions
    copysol = Dict(); copyobj = Dict();
    for i=1:length(Pobj)
        copysol[i] = P[i]
        copyobj[i] = Pobj[i]
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
function weightFP(candX,candY,x1,x2,n,C,Tabu,newsol,wn)
    λ = round(rand(Float64, 1)[1]; digits=1)
    x_t = x1*λ + x2*(1-λ); SearchDone = false; iter=0; Max_iter = 10
    while iter<Max_iter && SearchDone == false
        x_r = round.(Int,x_t); fx = getobjval(x_r,C)
        if ( (FBcheck(x_r,n) == true) && (x_r∉candX) )
            push!(candX,x_r); candY = [candY; getobjval(x_r,C)'];
            newsol+=1; wn+=1; SearchDone = true
        else
            if x_r ∈ Tabu
                x_r = flipoper(Tabu,x_t,x_r)
                if x_r==[]
                    SearchDone = true
                else
                    if ( (FBcheck(x_r,n) == true) && (x_r∉candX) )  #z[1]<=fmax[1]-1 && z[2]<=fmax[2]-1 && z[3]<=fmax[3]-1 )
                        push!(candX,x_r); candY = [candY; getobjval(x_r,C)'];
                        newsol+=1; wn+=1; SearchDone = true
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
		iter+=1
    end
    return candX,candY,newsol,wn,Tabu
end
function GPR(candX,candY,SI,SG,n,C,exploredSI,newsol,gn)
    iter=0;#t0=time();
    while all.(SI != SG) && iter<20
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
        iter+=1
    end
    return candX,candY,newsol,gn,exploredSI
end
function FP(candX,candY,x_t,n,C,Tabu,newsol,fn,t0,TL)
    SearchDone = false; iter=0; Max_iter = 10
    while iter<Max_iter && SearchDone == false
        x_r = round.(Int,x_t)
        if ( (FBcheck(x_r,n) == true) && x_r∉candX )
            push!(candX,x_r); candY = [candY; getobjval(x_r,C)']; #add new solval to Y
            newsol+=1; fn+=1;
            SearchDone = true
        else
            if x_r ∈ Tabu
                x_r = flipoper(Tabu,x_t,x_r)
                if x_r==[]
                    SearchDone = true
                else
                    if ( (FBcheck(x_r,n) == true) && x_r∉candX )
                        push!(candX,x_r); candY = [candY; getobjval(x_r,C)'];
                        newsol+=1; fn+=1;
                        SearchDone = true
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
		iter+=1
    end

    return candX,candY,newsol,fn,Tabu
end
function clumath(dvar,L,n,C,k,TL)
    SolPair = Set(); exploredX = [];  t0=time();
    IGPair=[]; exploredSI = []; Tabu = []; candX = copy(dvar); candY = copy(L);
    newsol=0; wn = 0; gn=0; fn=0;
    while time()-t0 < TL
        # @show time()-t0

        clu = kmeans(candY',max(k,2)); nc = nclusters(clu); sc = counts(clu);
        targec = findall(i->i in nsmallest(2, sc),sc)
        candx1 = findall(i->i==targec[1],clu.assignments); candx2 = findall(i->i==targec[2],clu.assignments)
        id1 = sample(candx1, 1, replace=false)[1]; id2 = sample(candx2, 1, replace=false)[1]
        x1 = candX[id1]; x2 = candX[id2];
        frt = (count(i->0<i<1,x1)+count(i->0<i<1,x2))/(n*2);

        # if frt<=0.02 # 2% fractionality
            # while Set([id1,id2])∉SolPair
            #     candX,candY,newsol,wn,Tabu = weightFP(candX,candY,x1,x2,n,C,Tabu,newsol,wn)
            #     push!(SolPair,Set([id1,id2]))
            # end
        # elseif  frt<=0.1 # 10% fractionality
        while [id1,id2]∉IGPair &&  time()-t0<TL
            candX,candY,newsol,gn,exploredSI = GPR(candX,candY,x1,x2,n,C,exploredSI,newsol,gn)
            push!(IGPair,[id1,id2])
        end
        # else # higher fractionality
        #     targex = findall(i->i in nsmallest(1, sc),sc);
        #     candxt = findall(i->i==targex[1],clu.assignments); xtid = sample(candxt, 1, replace=false)[1];
        #     x_t = candX[xtid];
        #     while xtid∉exploredX &&  time()-t0<TL
        #         candX,candY,newsol,fn,Tabu = FP(candX,candY,x_t,n,C,Tabu,newsol,fn,t0,TL)
        #         push!(exploredX,xtid)
        #     end
        # end
    end
    return candX,candY,newsol,wn,gn,fn
end
################################  Data  ####################################
data = Data(ARGS[1],ARGS[2]); pre = Val(ARGS[3],ARGS[4]);
######################### Mathematical Model #############################
mip = Model(with_optimizer(GLPK.Optimizer)); set_silent(mip);
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
dist = Model(with_optimizer(GLPK.Optimizer));set_silent(dist);
@variables(dist, begin
    0 <= dx[1:data.n] <=1
end)
for k=1:data.m
    if data.signs[k] == "l"
        @constraint(dist, dot(data.B[k,:],dx) >= data.RHS[k])
    elseif data.signs[k] == "u"
        @constraint(dist, dot(data.B[k,:],dx) <= data.RHS[k])
    else
        @constraint(dist, dot(data.B[k,:],dx) == data.RHS[k])
    end
end
optimize!(dist)
##################################################
Bentime = readdlm(ARGS[5])[1];
weightFP(pre.dvar,pre.L,pre.dvar[1],pre.dvar[2],data.n,data.C,[],0,0,0,5)
GPR(pre.dvar,pre.L,pre.dvar[1],pre.dvar[2],data.n,data.C,[],0,0,0,5)
FP(pre.dvar,pre.L,pre.dvar[1],data.n,data.C,[],0,0,0,5)
clumath(pre.dvar,pre.L,data.n,data.C,pre.k,5)# compiling
clutime = @CPUelapsed fpX,fpY,newsol,wn,gn,fn = clumath(pre.dvar,pre.L,data.n,data.C,pre.k,10)#TL-Bentime)

finalX,finalY = Postpro(fpX,fpY,newsol);
otable = zeros(length(finalY),3)
for i=1:length(finalY)
    for j=1:3
        otable[i,j] = finalY[i][j]
    end
end
ins = ARGS[2][45:end-4];
print("miptime: ", clutime+Bentime)
CSV.write("/home/k2g00/kg23475/generic/instances/"*"$ins"*"_Y.log",DataFrame(otable, :auto),append=false, header=false, delim=' ' )
record1 = DataFrame(Filename=ins,wnsol=wn, gnsol=gn, fpsol=fn, totalsol = length(finalY))
CSV.write("./miprecord.csv", record1, append=true, header=false )
