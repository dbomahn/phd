using DataStructures,DelimitedFiles,DataFrames,JuMP,GLPK,LinearAlgebra,CSV,StatsBase,Clustering,CPUTime
mutable struct Data
    dtfile::String; n::Int; C::Array{}; ub::Int; weight::Array{};
    function Data(dtfile::String)
        d = readdlm(dtfile)
        dt = readdlm(dtfile,'\t', String, '\n')
        b = dt[4:length(dt)-1]
        n=parse(Int,dt[2])
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
        ub=parse(Int,dt[3])
        weight =ones(1,n)
        weight = round.(Int,weight)
        item = dt[length(dt)]
        w1 = split(item, ('[',']',','))
        w2 = filter!(e->!(e in ["" ,"[", "]"]) ,w1)
        for i=1:n
            weight[i] = parse(Int64,w2[i])
        end
        new(dtfile,n,C,ub,weight)#,P)
    end
end
mutable struct Val
    x::String; y::String; n::Int; dvar::Array{}; LB::Array{};
    function Val(x,y,n)
        dv = round.(digits=4, readdlm(x))
        objs = round.(digits=2, readdlm(y))#[:,1:end]
        ind = findall(i-> 0 in objs[i,:]  , 1:size(objs)[1])
        dv2 = dv[setdiff(1:end, ind), :];
        L = -objs[setdiff(1:end, ind), 2:end];
        dvar = [Vector(dv2[i, :]) for i=1:size(L)[1]]
		LB = [Vector(L[i,:]) for i=1:size(L)[1]]
        new(x,y,n,dvar,LB)
    end
end
function getobjval(x,C)
    return [-dot(x,C[1,:]),-dot(x,C[2,:]),-dot(x,C[3,:])]
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
        JuMP.fix(x[k],xx[k]; force=true)
    end
    optimize!(kp_m)
    if termination_status(kp_m) == MOI.OPTIMAL
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


ti = [[1],[2],[3]]
tl = copy(ti)

while tl != []
	k = rand(1:length(tl))
	println(k)
	deleteat!(tl,k);
	@show tl
end
1
#############################    Feasibility Pump   ############################
function FP(candX,n,C,TL)
	X = []; Y = [];	Tabu = []; newsol=0; candlist = copy(candX);k=1; t0=time();
	# k = rand(1:length(candlist))
    while candlist != [] &&  time()-t0 < TL && k < length(candX)+1
        x_t = candlist[k]; SearchDone = false
		iter=0; Max_iter = 10
        while iter<Max_iter && SearchDone == false
            x_r = round.(Int,x_t)
            fx = getobjval(x_r,C)
            if ( (FBcheck(x_r,n) == true) &&  x_r∉[candX;X])
				push!(X,x_r); push!(Y,getobjval(x_r,C))
                newsol+=1; setdiff!(candlist,candX[k,:]); SearchDone = true
				#deleteat!(candlist,k)
            else
                if x_r ∈ Tabu
                    x_r = flipoper(Tabu,x_t,x_r)
                    if x_r==[]
                        SearchDone = true
                    else
                        fx = getobjval(x_r,C)
                        if ( (FBcheck(x_r,n) == true) &&  x_r∉[candX;X])
							push!(X,x_r); push!(Y,getobjval(x_r,C))
                            newsol+=1; setdiff!(candlist,candX[k,:]);SearchDone = true
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
        k+=1
    end
    return X,Y,candlist
end
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

function GPR(C,n,dvar,LB,TL)
    IGPair=[]; exploredSI = []; newsol=0; candX = copy(dvar); candY = copy(LB); t0=time(); #[vec(LB[i,:]) for i=1:size(LB)[1] ]
    while time()-t0 < TL
        I,G = sample(1:length(candX), 2, replace=false)
        SI = candX[I]; SG = candX[G]; iter=0; #Maxiter = length(dvar)*10;  #if the algorithm finds infeasible solutions in a row, move to the next iteration
        while all.(SI != SG) && [I,G]∉IGPair && iter<20 && (time()-t0<TL)
            dif = findall(i-> SI[i]!=SG[i], 1:n)
            neibour,neiobj = createNB(SI,C,dif,exploredSI)
            if length(neiobj)==0
                break
            else
                for l=1:length(neiobj)
                    if ( (FBcheck(neibour[l],n)==true) && neibour[l]∉candX ) #dominance check
                        push!(candX, neibour[l]); push!(candY,neiobj[l])
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
    return candX,candY
end
data = Data(ARGS[1]); pre = Val(ARGS[2],ARGS[3],data.n);
kp_m = Model(with_optimizer(GLPK.Optimizer));set_silent(kp_m)
@variable(kp_m, x[1:data.n], Bin );
@constraint( kp_m, sum(data.weight[i]*x[i] for i=1:data.n) <= data.ub );
optimize!(kp_m);
dist = Model(with_optimizer(GLPK.Optimizer));;set_silent(dist)
@variable(dist, 0<=dx[1:data.n]<=1)
@constraint( dist, sum(data.weight[i]*dx[i] for i=1:data.n) <= data.ub );
optimize!(dist);
if data.n<=10
    TL = 1;
elseif data.n<=20
    TL = 1.2;
elseif data.n<=30
    TL = 1.5;
elseif data.n<=40
    TL = 2;
elseif data.n<=50
    TL = 3.5;
elseif data.n<=60
    TL = 7;
elseif data.n<=70
    TL = 15;
elseif data.n<=80
    TL = 31;
elseif data.n<=90
    TL = 50;
else
    TL = 119;
end
FP(pre.dvar,data.n,data.C,5); GPR(data.C,data.n,pre.dvar,pre.LB,5); #compiling
FPtime = @CPUelapsed fx,fy,xu = FP(pre.dvar,data.n,data.C,TL/2)
GPRtime = @CPUelapsed candset,candobj = GPR(data.C,data.n,fx,fy,TL-FPtime)
finalX,finalY = domFilter(candset,candobj)

totaltime = FPtime+GPRtime
otable = zeros(Int, length(finalY),3)
for i=1:length(finalY)
    for j=1:3
        otable[i,j] = finalY[i][j]
    end
end

@show ins = ARGS[1][end-15:end-4]
print("kptime: ", totaltime)
# record1 = DataFrame(ins = ins, totalsol=length(finalY), time=totaltime) #fpsol = newsol, gpsol = newsol2,
# CSV.write("/home/k2g00/k2g3475/generic/kprecords.csv",record1,append=true, header=false )#, delim=',' )
CSV.write("/home/k2g00/k2g3475/generic/send/"*"$ins"*"kpY.log",DataFrame(otable, :auto), header=false, delim=' ' )
