using DataStructures,DelimitedFiles,DataFrames,JuMP,CPLEX,LinearAlgebra,CSV,Clustering,StatsBase
mutable struct Data
    dtfile::String; n::Int; C::Array{}; ub::Int; weight::Array{};#P::Array{}
    function Data(dtfile::String)
        d = readdlm(dtfile)
        data = readdlm(dtfile,'\t', String, '\n')
        b = data[4:length(data)-1]
        n=parse(Int,data[2])
        C= ones(length(b),n)
        C = round.(Int,C)
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
        ub=parse(Int,data[3])
        weight =ones(1,n)
        weight = round.(Int,weight)
        item = data[length(data)]
        w1 = split(item, ('[',']',','))
        w2 = filter!(e->!(e in ["" ,"[", "]"]) ,w1)
        for i=1:n
            weight[i] = parse(Int64,w2[i])
        end
        new(dtfile,n,C,ub,weight)#,P)
    end
end

mutable struct Val
    x::String; y::String; n::Int; dvar::Array{}; LB::Array{}; k::Int
    function Val(x,y,n)
        dv = round.(digits=4, readdlm(x))
        objs = round.(digits=2, readdlm(y))#[:,1:end]
        ind = findall(i-> 0 in objs[i,:]  , 1:size(objs)[1])
        dv2 = dv[setdiff(1:end, ind), :];
        LB = -objs[setdiff(1:end, ind), 2:end];
        dvar = [Vector(dv2[i, :]) for i = 1:size(LB)[1]]
        k = ceil(Int,size(LB)[1]/(n/2));
        new(x,y,n,dvar,LB,k)
    end
end

function getobjval(x,C)
    return [-dot(x,C[1,:]),-dot(x,C[2,:]),-dot(x,C[3,:])]
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
function FBcheck(xx,n)
    for k=1:n
        JuMP.fix(x[k],xx[k])
    end
    optimize!(kp_m)
    if termination_status(kp_m) == MOI.OPTIMAL
        return true
    else
        return false
    end
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
     #t0=time();
    λ = round(rand(Float64, 1)[1]; digits=1)
    x_t = x1*λ + x2*(1-λ); SearchDone = false; iter=0; Max_iter = 10 #max( round(Int,count(x->0<x<1,x_t)/5), 1 )
    while iter<Max_iter && SearchDone == false #time()-t0 < TL &&
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
            # if time()-t0 >= TL-1
            #     break
            # end
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
    while all.(SI != SG)  && iter<20 #&& (time()-t0<TL)
        dif = findall(i-> SI[i]!=SG[i], 1:n)
	    neibour = createNB(SI,C,dif,exploredSI)
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
        SI = nextSI(neibour,C,SI)
        if SI∉dvar
            push!(exploredSI,SI);
        end
        iter+=1
    end
    return candX,candY,newsol,gn
end

function FP(candX,candY,x_t,n,C,Tabu,newsol,fn)
     #t0=time();
    SearchDone = false; iter=0; Max_iter = 10 #max( round(Int,count(x->0<x<1,x_t)/5), 1 )
    while iter<Max_iter && SearchDone == false #time()-t0 < TL &&
        x_r = round.(Int,x_t)
        if ( (FBcheck(x_r,n) == true) && x_r∉candX ) #checking feasibility and incumbent sols
            push!(candX,x_r); candY = [candY; getobjval(x_r,C)']; #add new solval to Y
            newsol+=1; fn=+1;
            SearchDone = true
        else
            if x_r ∈ Tabu
                x_r = flipoper(Tabu,x_t,x_r)
                if x_r==[]
                    SearchDone = true
                else
                    if ( (FBcheck(x_r,n) == true) && x_r∉candX )
                        push!(candX,x_r); candY = [candY; getobjval(x_r,C)'];
                        newsol+=1; fn=+1;
                        SearchDone = true
                    end
                end
            end
            # if time()-t0 >= TL
            #     break
            # end
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

dt = Data("/home/ak121396//multiobjective/instances/KP/gpr/dat/KP_p-3_n-030_ins-1.dat")
pr = Val("/home/ak121396//multiobjective/instances/KP/gpr/X/KP_p-3_n-030_ins-1.x.sol","/home/ak121396//multiobjective/instances/KP/gpr/Y/KP_p-3_n-030_ins-1.y.sol")
kp_m = Model(CPLEX.Optimizer);
MOI.set(kp_m, MOI.RawParameter("CPX_PARAM_SCRIND"), false );
MOI.set(kp_m, MOI.RawParameter("CPX_PARAM_THREADS"),1  );
@variable(kp_m, x[1:dt.n], Bin );
@constraint( kp_m, -dot(dt.weight,x) >= -(dt.ub) );
optimize!(kp_m);

dist = Model(CPLEX.Optimizer);
MOI.set(dist, MOI.RawParameter("CPX_PARAM_SCRIND"), false );
MOI.set(dist, MOI.RawParameter("CPX_PARAM_THREADS"),1  );
@variable(dist, dx[1:dt.n], Bin );
@constraint( dist, -dot(dt.weight,dx) >= -(dt.ub) );
optimize!(dist);

function clumath(dvar,LB,n,C,k,TL)
    SolPair = Set(); exploredX = [];  t0=time();
    IGPair=[]; exploredSI = []; Tabu = []; candX = copy(dvar); candY = copy(LB);
    newsol=0; wn = 0; gn=0; fn=0;
    while time()-t0 < TL
        clu = kmeans(candY',max(k,2)); nc = nclusters(clu); sc = counts(clu);
        targec = findall(i->i in nsmallest(2, sc),sc)
        candx1 = findall(i->i==targec[1],clu.assignments); candx2 = findall(i->i==targec[2],clu.assignments)
        id1 = sample(candx1, 1, replace=false)[1]; id2 = sample(candx2, 1, replace=false)[1]
        x1 = candX[id1]; x2 = candX[id2];

        if x1 == round.(x1) && x2 == round.(x2) #both sols are int => wFP
            # println("weighted FP used")
            while Set([id1,id2])∉SolPair &&  time()-t0<TL
                candX,candY,newsol,wn,Tabu = weightFP(candX,candY,x1,x2,n,C,Tabu,newsol,wn)
                push!(SolPair,Set([id1,id2]))
            end
        elseif  x1 != round.(x1) && x2 != round.(x2) #both sols are frac => PR
            while [I,G]∉IGPair &&  time()-t0<TL
                candX,candY,newsol,gn = GPR(candX,candY,x1,x2,n,C,exploredSI,newsol,gn)
                push!(IGPair,[I,G])
            end
        else #one sol is frac/Int =>FP
            targex = findall(i->i in nsmallest(1, sc),sc);
            candxt = findall(i->i==targex[1],clu.assignments); xtid = sample(candxt, 1, replace=false)[1];
            x_t = candX[xtid];
            while xtid∉exploredX &&  time()-t0<TL
                candX,candY,newsol,fn,Tabu = FP(candX,candY,x_t,n,C,Tabu,newsol,fn)
                push!(exploredX,xtid)
            end
        end
    end
    return candX,candY,newsol,wn,gn,fn
end

X,Y,nsol,wn,gn,fn = clumath(pr.dvar,pr.LB,dt.n,dt.C,pr.k,10)
finalx,finaly = Postpro(X,Y,nsol)
candset,candobj,nbtime,SItime,FBtime,newsol = GPR(dt.C,dt.n,pr.dvar,pr.LB,10);
finalX,finalY = Postpro(candset,candobj,newsol)

data = Data("/home/k2g00/k2g3475/multiobjective/instances/KP/gpr/dat/KP_p-3_n-10_ins-1.dat")
pre = Val("/home/k2g00/k2g3475/multiobjective/instances/KP/gpr/X/KP_p-3_n-10_ins-1.x.sol","/home/k2g00/k2g3475/multiobjective/instances/KP/gpr/Y/KP_p-3_n-10_ins-1.y.sol",dt.n)
kp_m = Model(CPLEX.Optimizer);
MOI.set(kp_m, MOI.RawParameter("CPX_PARAM_SCRIND"), false );
MOI.set(kp_m, MOI.RawParameter("CPX_PARAM_THREADS"),1  );
@variable(kp_m, x[1:data.n], Bin );
@constraint( kp_m, -dot(data.weight,x) >= -(data.ub) );
optimize!(kp_m);

if data.n<=60
    TL = 10;
elseif data.n<=80
    TL = 30;
else
    TL = 160;
end
runtime1 = @CPUelapsed candset,candobj,nbtime,SItime,FBtime,newsol = GPR(data.C,data.n,pre.dvar,pre.LB,TL)
runtime2 = @CPUelapsed finalX,finalY = Postpro(candset,candobj,newsol)
push!(finalY,[0,0,0])
totaltime = runtime1+runtime2
otable = zeros(Int, length(finalY),3)
for i=1:length(finalY)
    for j=1:3
        otable[i,j] = finalY[i][j]
    end
end
matriX = zeros(Int,length(finalX),data.n)
for i=1:length(finalX)
    for j=1:data.n
        matriX[i,j] = finalX[i][j]
    end
end
@show ins = ARGS[1][end-15:end-4]
class = ins[1:end-5]

record1 = DataFrame(totalsol=length(finalY), numnewsol = newsol, createnb = nbtime, nextSI = SItime, feastime = FBtime)
CSV.write("/home/k2g00/k2g3475/multiobjective/solvers/generalPR/goutputs/records/"*"$class"*".csv",record1,append=true, header=false )#, delim=',' )
io = open("/home/k2g00/k2g3475/multiobjective/solvers/generalPR/goutputs/time/"*"$ins"*".txt", "a")
println(io,"$totaltime"); close(io) # header=false )#, delim=',' )
CSV.write("/home/k2g00/k2g3475/multiobjective/solvers/generalPR/goutputs/Y/"*"$ins"*"_Y.log",DataFrame(otable, :auto), header=false, delim=' ' )
sparX = sparse(matriX); JLD2.@save "/home/k2g00/k2g3475/multiobjective/solvers/generalPR/goutputs/X/"*"$ins"*"_X.jdl2" sparX
