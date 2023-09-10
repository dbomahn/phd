
using DataStructures,DelimitedFiles,DataFrames,JuMP,GLPK,LinearAlgebra,CPUTime,CSV,Statistics,Clustering,MathOptInterface,StatsBase,SparseArrays,JLD2

mutable struct Data
    input::String
    i::Int; j::Int; n::Int; C::Array{}
    function Data(input::String)
        data=readdlm(input)
        i,j = filter!(a->typeof(a)==Int, data[1,:]) # facility_i, customer_j # for FLP
        n = i+j+i*j
        fixcost,capa,cost = [],[],[]
        for k=2:i+1
            fix,cp = filter!(a->typeof(a)==Int, data[k,:])
            push!(fixcost,fix)
        end
        demand = filter!(a->typeof(a)==Int, data[2+i,:])
        a,b = size(data)
        for k=i+3:a
            ct = filter!(a->typeof(a)==Int, data[k,:])
            # push!(cost2,ct)
            for l=1:length(ct)
                push!(cost,ct[l])
            end
        end
        C = [fixcost;demand;cost]
        new(input,i,j,n,C)#fixcost,demand,cost,cost2
    end
end
mutable struct Val
    x::String; y::String; j::Int; dvar::Array{}; LB::Array{}; LBmtx::Array{}; radius::Float64
    function Val(x,y,j)
        dv = round.(digits=4, readdlm(x))
        objs = round.(digits=4, readdlm(y))#[:,1:end]
        ind = findall(i-> 0 in objs[i,:]  , 1:size(objs)[1])
        dv2 = dv[setdiff(1:end, ind), :];
        LBmtx = objs[setdiff(1:end, ind), 2:end];
        dvar = [Vector(dv2[i, :]) for i = 1:size(LBmtx)[1]]
        LB = [Vector(LBmtx[i, :]) for i = 1:size(LBmtx)[1]]
        maxobj = [maximum(LBmtx[:,i]) for i=1:3]; minobj = [minimum(LBmtx[:,i]) for i=1:3];
        steps = round.(Int,(maxobj-minobj)/j)  #length(LB[:,1])
        interval = round(Int,mean(maxobj-minobj)/j)
        # points = transpose(LB)
        α = max(round(length(LB)*0.05),1)
        radius = interval*α
        new(x,y,j,dvar,LB,LBmtx,radius)
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
function getobjval(x,C) #the order: fixcost,-demand,cost
    return [dot(x[1:dt.i],C[1:dt.i]),-dot(x[dt.i+1:dt.i+dt.j],C[dt.i+1:dt.i+dt.j]),dot(x[dt.i+dt.j+1:end],C[dt.i+dt.j+1:end])]
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
function FBcheck(xx,n)
    for k=1:n
        JuMP.fix(x[k],xx[k])
    end
    optimize!(flp)
    # print("status: ", termination_status(flp), "\n" )
    if termination_status(flp) == MOI.OPTIMAL
        return true
    else
        return false
    end
end
function Postpro(dvar, LB, newsol)
    #Filter fractional solutions from LB
    initdv = dvar[1:end-newsol]
    initLB = LB[1:end-newsol]
    frac = findall(j -> trunc.(initdv[j]) != initdv[j], 1:length(initdv))
    dv2 = initdv[setdiff(1:end, frac)]
    LB2 = initLB[setdiff(1:end, frac)]
    P = union(dv2, dvar[end-newsol+1:end])
    Pobj = union(LB2, LB[end-newsol+1:end])

    #Filter dominated solutions
    copysol = Dict()
    copyobj = Dict()
    for i = 1:length(Pobj)
        copysol[i] = P[i]
        copyobj[i] = Pobj[i]
    end
    for i = 1:length(Pobj)-1
        for j = i+1:length(Pobj)
            if all(Pobj[i] .>= Pobj[j]) == true #dominated by PF[j]
                copyobj[i] = nothing
                copysol[i] = nothing
                break
            elseif all(Pobj[j] .>= Pobj[i]) == true
                copyobj[j] = nothing
                copysol[j] = nothing
            end
        end
    end

    finalsol = filter!(a -> a != nothing, collect(values(copysol)))
    finalobj = filter!(a -> a != nothing, collect(values(copyobj)))

    return finalsol, finalobj
end


function FFP(dvar,LB,n,C,TL)
    X = copy(dvar); IGPair=[]; Tabu = []; newsol=0; t0=time(); Y = copy(LB);
    while time()-t0 < TL
        I,G = sample(1:length(X), 2, replace=false)
        x1 = X[I]; x2 = X[G];
        λ = round(rand(Float64, 1)[1]; digits=1)
        x_t = x1*λ + x2*(1-λ); SearchDone = false; iter=0; Max_iter = 10
        while [I,G]∉IGPair && iter<Max_iter && SearchDone == false
            x_r = round.(Int,x_t)
            if ( (FBcheck(x_r,n) == true) && (x_r∉X) )
                push!(X,x_r); push!(Y, getobjval(x_r,C));
                newsol+=1; SearchDone = true
            else
                if x_r ∈ Tabu
                    x_r = flipoper(Tabu,x_t,x_r)
                    if x_r==[]
                        SearchDone = true
                    else
                        if ( (FBcheck(x_r,n) == true) && (x_r∉X) )
                            push!(X,x_r); push!(Y, getobjval(x_r,C)); #Y = [Y; getobjval(x_r,C)'];
                            SearchDone = true; newsol+=1; 
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
    return X,Y,newsol
end

#dt = Data("/home/k2g00/k2g3475/pastinstance/flp/dat/05_010_04.txt")
#pre = Val("/home/k2g00/k2g3475/pastinstance/flp/varval/05_10_04_pre_img_p.sol","/home/k2g00/k2g3475/pastinstance/flp/PF/05_10_04_img_p.sol",dt.j)


dt = Data(ARGS[1]); pre = Val(ARGS[2],ARGS[3],dt.j);
#Bentime = readdlm(ARGS[4])[1];
tnum = dt.i;
if tnum == 5
    TL =2 ; Bentime = 0.1;
elseif tnum == 10
    TL =6; Bentime = 0.3;
elseif tnum == 15
    TL = 50; Bentime = 1.2;
elseif tnum == 20
    TL =270; Bentime = 3.7;
elseif tnum == 25
    TL =599; Bentime = 9.4;
elseif tnum == 30
    TL =599; Bentime = 20.1;
elseif tnum == 35
    TL =599; Bentime = 39.1;
elseif tnum == 40
    TL =599; Bentime = 62.5;
elseif tnum == 45
    TL =599; Bentime = 105.9;
elseif tnum == 50
    TL =599; Bentime = 165.3;
elseif tnum == 55
    TL =600; Bentime = 246.1;
else
    TL = 600; Bentime = 380.3;
end

FFPTL = TL-Bentime;
dist = Model(GLPK.Optimizer);
@variable(dist, 0<= dx[1:dt.n] <= 1);
@constraint(dist, cond1[b=1:dt.j], sum(dx[dt.i+dt.j+dt.j*(a-1)+b] for a in 1:dt.i) == dx[dt.i+b]);
@constraint(dist, cond2[a=1:dt.i,b=1:dt.j], dx[dt.i+dt.j+dt.j*(a-1)+b] <= dx[a]);
optimize!(dist);

flp = Model(GLPK.Optimizer);
@variable(flp, x[1:dt.n] ,Bin);
@constraint(flp, con1[b=1:dt.j], sum(x[dt.i+dt.j+dt.j*(a-1)+b] for a in 1:dt.i) == x[dt.i+b]);
@constraint(flp, con2[a=1:dt.i,b=1:dt.j], x[dt.i+dt.j+dt.j*(a-1)+b] <= x[a]);
optimize!(flp);
#COMPILING
FFP(pre.dvar,pre.LB,dt.n,dt.C,2);

runtime1 = @CPUelapsed candset,candobj,newsol = FFP(pre.dvar,pre.LB,dt.n,dt.C,FFPTL);
runtime2 = @CPUelapsed finalX,finalY = Postpro(candset,candobj,newsol);

#push!(finalY,[0,0,0])
totaltime = runtime1+runtime2+Bentime
otable = zeros(length(finalY),3)
for i=1:length(finalY)
    for j=1:3
        otable[i,j] = finalY[i][j]
    end
end
@show ins = ARGS[2][1:end-10]
print("$ins time: ", totaltime)
#io = open(ARGS[3]*"_time.txt", "a"); println(io,"$totaltime"); close(io)
CSV.write(ins*"_Y.log",DataFrame(otable, :auto), header=false, delim=' ' )

