using DataStructures,DelimitedFiles,DataFrames,JuMP,CPLEX,LinearAlgebra,CSV,StatsBase,Clustering
# using CPUtime

mutable struct Data
    input::String; n::Int; C::Array{}; B::Array{};
    function Data(input::String)
        f=readdlm(input, '\t', String, '\n')
        obj=parse(Int,f[1])
        n=parse(Int,f[2])
        P=[]
        for i=3:length(f)
            a = split(f[i], ('[',']',','))
            a = filter!(e->e ∉ [ "" ,"[","," ,"]"] ,a)
            if length(a)>2
                push!(P,a)
            end
        end
        C = ones(obj,n*n); C = round.(Int,C);
        for x=1:length(P)
            for y=1:n
                p=parse(Int64,P[x][y])
                idx = Int(floor((x-1)/n))
                ind = ((x-1)*n+y)%(n*n)
                if ind != 0
                    C[idx+1,((x-1)*n+y)%(n*n)] = p
                else
                    C[idx+1,(n*n)] = p
                end
            end
        end
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
        new(input,n*n,C,B)
    end
end
mutable struct Val
    x::String; y::String; n::Int; dvar::Array{}; LB::Array{}; k::Int
    function Val(x,y,n)
        dv = round.(digits=4, readdlm(x))
        objs = round.(digits=4, readdlm(y))
        ind = findall(i-> 0 in objs[i,:]  , 1:size(objs)[1])
        dv2 = dv[setdiff(1:end, ind), :];
        LB = objs[setdiff(1:end, ind), 2:end];
        dvar = [Vector(dv2[i, :]) for i = 1:size(LB)[1]]
        # LB = [Vector(LBmtx[i, :]) for i = 1:size(LBmtx)[1]]
        k = ceil(Int,size(LB)[1]/Int(sqrt(n)));
        new(x,y,n,dvar,LB,k)
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
        JuMP.fix(x[k],xx[k]; force=true)
    end
    optimize!(ap_m)
    if termination_status(ap_m) == MOI.OPTIMAL
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
# function wPR(candX,candY,n,C,TL)
#     dvar = copy(candX); LB = copy(candY)
#     IGPair=Set(); exploredSI = []; newsol=0; t0=time();
#     while time()-t0 < TL
#         I,G = sample(1:length(dvar), 2, replace=false)
#         SI = dvar[I]; SG = dvar[G]; iter=0;
#         λ = round(rand(Float64, 1)[1]; digits=1); x_t = round.(SI*λ + SG*(1-λ); digits=1);
#         if  Set([I,G])∉IGPair
#             while all.(SI != x_t) && iter<10 && (time()-t0<TL)
#                 dif = findall(i-> SI[i]!=x_t[i], 1:n)
#                 neibour,neiobj = createNB(SI,C,dif,exploredSI)
#                 if length(neiobj)==0
#                     break
#                 else
#                     for l=1:length(neiobj)
#                         feasi = FBcheck(neibour[l],n)
#                         if feasi == true && neibour[l]∉ dvar
#                             push!(dvar, neibour[l]); push!(LB, neiobj[l]);
#                             newsol+=1
#                         end
#                     end
#                     SI = nextSI(neibour,neiobj,C,SI)
#                     if SI∉dvar
#                         push!(exploredSI,SI);
#                     end
#                 end
#                 iter+=1
#             end
#             while all.(x_t!=SG) && iter<20 && (time()-t0<TL)
#                 dif = findall(i-> x_t[i]!=SG[i], 1:n)
#                 neibour,neiobj = createNB(SI,C,dif,exploredSI)
#                 if length(neiobj)==0
#                     break
#                 else
#                     for l=1:length(neiobj)
#                         feasi = FBcheck(neibour[l],n)
#                         if feasi == true && neibour[l]∉ dvar
#                             push!(dvar, neibour[l]); push!(LB, neiobj[l]);
#                             newsol+=1
#                         end
#                     end
#                     SI = nextSI(neibour,neiobj,C,SI)
#                     if SI∉dvar
#                         push!(exploredSI,SI);
#                     end
#                 end
#                 iter+=1
#             end
#             push!(IGPair,Set([I,G]))
#         end
#     end
#     return dvar,LB,newsol
# end

# dt = Data("/home/k2g00/k2g3475/multiobjective/instances/AP/gpr/dat/AP_p-3_n-05_ins-01.dat")
# pr = Val("/home/k2g00/k2g3475/multiobjective/instances/AP/gpr/X/AP_p-3_n-5_ins-01_X.sol","/home/ak121396//multiobjective/instances/KP/gpr/Y/AP_p-3_n-5_ins-01_Y.sol")
dt = Data("/home/ak121396/Desktop/instances/AP/dat/AP_p-3_n-05_ins-1.dat")
pr = Val("/home/ak121396/Desktop/solvers/Bensolve/APoutputs/X/AP_p-3_n-05_ins-01_X.sol","/home/ak121396/Desktop/solvers/Bensolve/APoutputs/Y/AP_p-3_n-05_ins-01_Y.sol",dt.n)
ap_m = Model(CPLEX.Optimizer);
MOI.set(ap_m, MOI.RawParameter("CPX_PARAM_SCRIND"), false );
MOI.set(ap_m, MOI.RawParameter("CPX_PARAM_THREADS"),1  );
@variable( ap_m, x[1:dt.n], Bin );
for i=1:Int(sqrt(dt.n))*2
    @constraint( ap_m, dot(dt.B[i,:],x) == 1 );
end
optimize!(ap_m);

runtime1 = @CPUelapsed candset,candobj,newsol = GPR(data.C,data.n,pre.dvar,pre.LB,TL)
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
class = ins[1:end-6]

####################Lab AP#########
# pr = Val("/home/ak121396/Desktop/solvers/Bensolve/APoutputs/X/AP_p-3_n-05_ins-01_X.sol","/home/ak121396/Desktop/solvers/Bensolve/APoutputs/Y/AP_p-3_n-05_ins-01_Y.sol")
# dt = Data("/home/ak121396/Desktop/instances/AP/dat/AP_p-3_n-05_ins-1.dat")
dt = Data("E:/bensolve/ap/dat/AP_p-3_n-5_ins-1.dat"); pr = Val("E:/bensolve/ap/X/AP_p-3_n-05_ins-01_X.sol","E:/bensolve/ap/Y/AP_p-3_n-05_ins-01_Y.sol")
ap_m = Model(with_optimizer(CPLEX.Optimizer));
@variable( ap_m, x[1:dt.n],Bin);
for i=1:Int(sqrt(dt.n))*2
    @constraint( ap_m, dot(dt.B[i,:],x) == 1 );
end
optimize!(ap_m);
dist = Model(with_optimizer(CPLEX.Optimizer));
@variable(dist, 0<=dx[1:dt.n]<=1)
for i=1:Int(sqrt(dt.n))*2
    @constraint( dist, dot(dt.B[i,:],dx) == 1 );
end
optimize!(dist);
MOI.set(ap_m, MOI.RawParameter("CPX_PARAM_SCRIND"), false );
MOI.set(ap_m, MOI.RawParameter("CPX_PARAM_THREADS"),1  );
MOI.set(dist, MOI.RawParameter("CPX_PARAM_SCRIND"), false );
MOI.set(dist, MOI.RawParameter("CPX_PARAM_THREADS"),1  );

fx,fy,fn = fracFP(pr.dvar,pr.LB,dt.n,dt.C,120)

wfx,wfy = Postpro(fx,fy,fn)
gx,gy,gn = wPR(pr.dvar,pr.LB,dt.n*dt.n,dt.C,2)
wpx,wpy = Postpro(gx,gy,gn)


##############kp
pr = Val("/home/ak121396/Desktop/solvers/Bensolve/KPoutputs/X/KP_p-3_n-20_ins-2.x.sol","/home/ak121396/Desktop/solvers/Bensolve/KPoutputs/Y/KP_p-3_n-20_ins-2.y.sol")
dt = Data("/home/ak121396/Desktop/instances/KP/dat/KP_p-3_n-20_ins-2.dat")
kp_m = Model(with_optimizer(CPLEX.Optimizer));
@variable(kp_m, x[1:dt.n], Bin );
@constraint( kp_m, -dot(dt.weight,x) >= -(dt.ub) );
optimize!(kp_m);
dist = Model(with_optimizer(CPLEX.Optimizer));
@variable(dist, 0<=dx[1:dt.n]<=1)
@constraint( dist, -dot(dt.weight,dx) >= -(dt.ub) );
optimize!(dist);
MOI.set(kp_m, MOI.RawParameter("CPX_PARAM_SCRIND"), false );
MOI.set(kp_m, MOI.RawParameter("CPX_PARAM_THREADS"),1  );
MOI.set(dist, MOI.RawParameter("CPX_PARAM_SCRIND"), false );
MOI.set(dist, MOI.RawParameter("CPX_PARAM_THREADS"),1  );

############flp
flp = Model(with_optimizer(CPLEX.Optimizer));
@variable(flp, x[1:dt.n] ,Bin);
@constraint(flp, con1[b=1:dt.j], sum(x[dt.i+dt.j+dt.j*(a-1)+b] for a in 1:dt.i) == x[dt.i+b]);
@constraint(flp, con2[a=1:dt.i,b=1:dt.j], x[dt.i+dt.j+dt.j*(a-1)+b] <= x[a]);
optimize!(flp)
dist = Model(with_optimizer(CPLEX.Optimizer));
@variable(dist, 0<=dx[1:dt.n]<=1)
@constraint(dist, con1[b=1:dt.j], sum(dx[dt.i+dt.j+dt.j*(a-1)+b] for a in 1:dt.i) == dx[dt.i+b]);
@constraint(dist, con2[a=1:dt.i,b=1:dt.j], dx[dt.i+dt.j+dt.j*(a-1)+b] <= dx[a]);
optimize!(dist);
MOI.set(flp, MOI.RawParameter("CPX_PARAM_SCRIND"), false );
MOI.set(flp, MOI.RawParameter("CPX_PARAM_THREADS"),1  );
MOI.set(dist, MOI.RawParameter("CPX_PARAM_SCRIND"), false );
MOI.set(dist, MOI.RawParameter("CPX_PARAM_THREADS"),1  );

xx,yy,nn=fracFP(pr.dvar,pr.LB,dt.n,dt.C,5,dt.i,dt.j)
