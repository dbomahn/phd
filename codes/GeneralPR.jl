using DataStructures,DelimitedFiles,DataFrames,JuMP,CPLEX,LinearAlgebra,CPUTime,CSV,Statistics,Clustering,MathOptInterface,StatsBase

mutable struct Data
    input::String
    i::Int; j::Int; C::Array{}
    # fixcost::Array{}
    # demand::Array{}
    # cost::Array{}
    # cost2::Array{}

    function Data(input::String)
        data=readdlm(input)
        i,j = filter!(a->typeof(a)==Int, data[1,:]) # facility_i, customer_j # for FLP
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

        new(input,i,j,C)#fixcost,demand,cost,cost2
    end
end

mutable struct Val
    x::String; y::String; j::Int; dvar::Array{}; LB::Array{};radius::Float64
    function Val(x,y,j)
        dv = round.(digits=4, readdlm(x))
        objs = round.(digits=4, readdlm(y))#[:,1:end]
        ind = findall(i-> 0 in objs[i,:]  , 1:size(objs)[1])
        dvar = dv[setdiff(1:end, ind), :];
        LB = objs[setdiff(1:end, ind), 2:end];
        maxobj = [maximum(LB[:,i]) for i=1:3]; minobj = [minimum(LB[:,i]) for i=1:3];
        steps = round.(Int,(maxobj-minobj)/j)  #length(LB[:,1])
        interval = round(Int,mean(maxobj-minobj)/j)
        # points = transpose(LB)
        α = max(round(length(LB[:,1])*0.05),1)
        radius = interval*α

        new(x,y,j,dvar,LB,radius)
    end
end

function clustering(LB,radius)
    initclu = dbscan(LB', radius) #, min_neighbors = 1, min_cluster_size = 1)
    core = getproperty.(initclu, :core_indices)
    boundary = getproperty.(initclu, :boundary_indices)
    LBclu = union.(core,boundary); Icluster = Dict();
    clsize = getproperty.(initclu, :size)
    clulist = unique(sort(clsize))
    candI = findall(x-> x == clulist[1], clsize)
    candG = findall(x-> x == clulist[end], clsize)
    clumean = Dict()
    for i=1:length(LBclu)
        if core[i] != []
            clumean[i] = round.(digits=2, mean.(LB[core[i],:][:,k] for k=1:3))
        else
            clumean[i] = round.(digits=2, mean.(LB[boundary[i],:][:,k] for k=1:3))
        end
    end
    cores = collect(values(SortedDict(clumean)))
    return LBclu,clsize,clulist,cores,candI,candG
end

function createNB(di,dj,SI,C,dif,exploredSI)
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
            push!(neiobj, getobjval(cpSI,C,di,dj))
        else #if vari is fractional
            cpSI[i] = 1; push!(neibour, cpSI); push!( neiobj, getobjval(cpSI,C,di,dj) )
            cpSI = copy(SI)
            cpSI[i] = 0; push!(neibour, cpSI); push!( neiobj, getobjval(cpSI,C,di,dj) )
        end
    end
    idx = findall(i-> neibour[i] in exploredSI, 1:length(neibour))
    neibour = setdiff(neibour,exploredSI)    #removed already explored SI
    #or neibour = neibour[setdiff(1:end, idx),:]
    neiobj = neiobj[setdiff(1:end, idx),:]

    return neibour,neiobj
end
function getobjval(x,C,i,j) #the order: fixcost,-demand,cost
    return [dot(x[1:i],C[1:i]),-dot(x[i+1:i+j],C[i+1:i+j]),dot(x[i+j+1:end],C[i+j+1:end])]
end

function nextSI(di,dj,neibour,neiobj,C,SI,candI,cores,clulist,clsize,LBclu,dvar)
    closclu = []
    for k=1:length(neiobj)
        newsol = neiobj[k,:][1]
        gap = [cores[i]-newsol for i=1:length(cores)]
        dis = [sum(gap[i].^2) for i=1:length(gap)]
        clos = findall( x-> x==minimum(dis), dis )[1]
        push!(closclu,clos)
    end

    SIobj = getobjval(SI,C,di,dj)
     #check if intermediate path is in the smallest cluster
    candid = findall(x-> x in candI, closclu)
    for i in clulist
        if length(candid) == 1  #if there is one candiate sol
            return neibour[candid[1]],neiobj[candid[1]]
        elseif length(candid) > 1 # if there are multiple candiates, check the improved ratio
            ratiotb = zeros(length(candid),3)
            for (i,j) in enumerate(candid)
                ratiotb[i,:] = neiobj[j]./SIobj
            end
            ranktb = zeros(length(candid),3)
            for i=1:3
                ranktb[:,i] = tiedrank(ratiotb[:,i])
            end
            ranksum = [sum(ranktb[i,:]) for i=1:length(candid)]
            hval = maximum(ranksum)
            mostimp = findall(x-> x == hval, ranksum)
            k = candid[rand(mostimp)[1]]
            return neibour[k], neiobj[k]
        else #there is no candidate
            # print("randomly select a candidate \n")
            k = rand(1:length(neibour)) #randomly select a candidate
            return neibour[k],neiobj[k]
        end
    end
    # @show k,neibour

    # return neibour[k],neiobj[k]
end


function ConstCheck(xx,i,j)
    for k=1:i+j+(i*j)
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
flp = Model(CPLEX.Optimizer)
# set_optimizer(flp, CPLEX.Optimizer)
MOI.set(flp, MOI.RawParameter("CPX_PARAM_SCRIND"), false );
MOI.set(flp, MOI.RawParameter("CPX_PARAM_THREADS"),1  );
@variable(flp, x[1:data.i+data.j+(data.i*data.j)] );
# unregister(flp, :con1),unregister(flp, :con2)
@constraint(flp, con1[b=1:data.j], sum(x[data.i+data.j+data.j*(a-1)+b] for a in 1:data.i) == x[data.i+b]);
@constraint(flp, con2[a=1:data.i,b=1:data.j], x[data.i+data.j+data.j*(a-1)+b] <= x[a]);
optimize!(flp)
ConstCheck(pre.dvar[1,:],5,10)



function FBcheck(x,i,j)
    for k=1:i+j+(i*j)
        if (x[k] == 0 || x[k]==1)
        # if JuMP.is_binary(x[k]) == true #or use set_binary(x[k]); optimize!(flp);
        else
            return false
        end
    end
    return true
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
                copyobj[i]=nothing; copysol[i]=nothing; break
            elseif all(obj[j] .>= obj[i]) == true
                copyobj[j]=nothing; copysol[j]=nothing;
            end
        end
    end

    finalsol = filter!(a->a!=nothing, collect(values(copysol)))
    finalobj = filter!(a->a!=nothing, collect(values(copyobj)))

    return finalsol,finalobj
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

function Postpro(dvar,LB,C,candX,candobj)
    frac = findall(j-> trunc.(LB[j,:])!= LB[j,:], 1:size(LB)[1])
    dvar = dvar[setdiff(1:end,frac),:]; LB = LB[setdiff(1:end,frac),:]

    Pobj = []; P = copy(dvar)
    for i=1:size(LB)[1]
        push!(Pobj,LB[i,:])
    end
    newsol = 0

    for i=1:length(candX)
        if dominated(candobj[i],Pobj)==false && candobj[i]∉Pobj
            P = [P; transpose(candX[i][:])]
            push!(Pobj,candobj[i][:])
            newsol+=1
        end
    end
    return P,Pobj,newsol
end
function GeneralPR(di,dj,C,dvar,LB,radius)
    candset = []; candobj=[]; IGPair=[];
    clutime = 0; nbtime = 0; SItime = 0; Contime=0; FBtime = 0;
    clutime= clutime + @CPUelapsed LBclu,clsize,clulist,cores,candI,candG = clustering(LB,radius)

    for i=1:100
        I = rand(LBclu[rand(candI)]); G = rand(LBclu[rand(candG)])
        SI = dvar[I,:]; SG = dvar[G,:]; iter=1; exploredSI = [];

        while all.(SI != SG) && [I,G]∉IGPair && iter<di*dj+di+dj
            rg = range(1, length=di*dj+di+dj)
            dif = findall(i-> SI[i]!=SG[i], rg)
            nbtime = nbtime + @CPUelapsed neibour,neiobj = createNB(di,dj,SI,C,dif,exploredSI)
            # if length(neibour) == 0
            #     break
            # end

            for l=1:length(neiobj)
                Contime = Contime + @CPUelapsed constcheck = ConstCheck(neibour[l],di,dj)
                if constcheck == true && neibour[l]∉ candobj
                    dvar = [dvar; transpose(neibour[l])]
                    LB = [LB; transpose(neiobj[l])]
                end
            end
            SItime = SItime + @CPUelapsed SI,SIobj = nextSI(di,dj,neibour,neiobj,C,SI,candI,cores,clulist,clsize,LBclu,dvar)

            FBtime = FBtime + @CPUelapsed fbcheck = FBcheck(SI,di,dj)
            if fbcheck == true && SI∉candset
                push!(exploredSI,SI); push!(candset,SI); push!(candobj,SIobj);
                # print("feasible new sol added \n")
            end
            iter+=1
        end
        push!(IGPair,[I,G])
    end
    return candset,candobj,clutime,nbtime,SItime,Contime,FBtime
end




# data=Data(ARGS[1]); pre=Val(ARGS[2],ARGS[3],data.j)
# totaltime = @CPUelapsed candX,candobj,clutime,nbtime,SItime,FBtime = GeneralPR(data.i,data.j,data.C,pre.dvar,pre.LB,pre.radius)
# P,Pobj,newsol = Postpro(pre.dvar,pre.LB,data.C,candX,candobj) #runtime2
# print(" Instance number: "*"$i"*" &  newly found sol: ", newsol,"\n") #" # final sols : ",length(Pobj),
#
# otable = ones(Int, length(Pobj),3)
# for i=1:length(Pobj)
#     for j=1:3
#         otable[i,j] = Pobj[i][j]
#     end
# end
# dir = ARGS[1][1:end-12]
# ins = ARGS[1][end-12:end-3]
# record1 = DataFrame(initsol=size(pre.LB)[1], newsol = newsol, clustering = clutime, createnb = nbtime, nextSI = SItime, feasicheck=FBtime, Algotime=totaltime)#,
# CSV.write("$dir"*"/record/"*"$rname"*".csv",record1, append=true, header=false )#, delim=',' )
# CSV.write("$dir"*"/Y/"*"$ins"*"_Y.log",DataFrame(otable),header=false, delim=' ' )

paths = ("/home/ak121396/Desktop//instances/FLP/instances/","/home/ak121396/Desktop//instances/FLP/varval/","/home/ak121396/Desktop//instances/FLP/PF/")
# readdir(paths[1]);readdir(paths[2]);readdir(paths[3])
data = Data( paths[1]*readdir(paths[1])[4] )
pre = Val(paths[2]*readdir(paths[2])[4],paths[3]*readdir(paths[3])[4],data.j)



flp = Model(CPLEX.Optimizer)
MOI.set(flp, MOI.RawParameter("CPX_PARAM_SCRIND"), false );
MOI.set(flp, MOI.RawParameter("CPX_PARAM_THREADS"),1  );
@variable(flp, x[1:data.i+data.j+(data.i*data.j)] );
@constraint(flp, con1[b=1:data.j], sum(x[data.i+data.j+data.j*(a-1)+b] for a in 1:data.i) == x[data.i+b]);
@constraint(flp, con2[a=1:data.i,b=1:data.j], x[data.i+data.j+data.j*(a-1)+b] <= x[a]);
optimize!(flp)
1
for i=11:20

    data = Data( paths[1]*readdir(paths[1])[i] )
    pre = Val(paths[2]*readdir(paths[2])[i],paths[3]*readdir(paths[3])[i],data.j)

    flp = Model(CPLEX.Optimizer) #with_optimizer(CPLEX.Optimizer) for 0.20 version
    MOI.set(flp, MOI.RawParameter("CPX_PARAM_SCRIND"), false );
    MOI.set(flp, MOI.RawParameter("CPX_PARAM_THREADS"),1  );
    @variable(flp, x[1:data.i+data.j+(data.i*data.j)] );
    @constraint(flp, con1[b=1:data.j], sum(x[data.i+data.j+data.j*(a-1)+b] for a in 1:data.i) == x[data.i+b]);
    @constraint(flp, con2[a=1:data.i,b=1:data.j], x[data.i+data.j+data.j*(a-1)+b] <= x[a]);
    optimize!(flp);

    totaltime = @CPUelapsed candX,candobj,clutime,nbtime,SItime,Contime,FBtime = GeneralPR(data.i,data.j,data.C,pre.dvar,pre.LB,pre.radius)
    P,Pobj,newsol = Postpro(pre.dvar,pre.LB,data.C,candX,candobj) #runtime2
    print(" Instance number: "*"$i"*" &  newly found sol: ", newsol,"\n") #" # final sols : ",length(Pobj),

    otable = ones(Int, length(Pobj),3)
    for i=1:length(Pobj)
        for j=1:3
            otable[i,j] = Pobj[i][j]
        end
    end
    ins = readdir(paths[1])[i][1:end-4]
    class = ins[1:6]
    record1 = DataFrame(initsol=size(pre.LB)[1], newsol = newsol, clustering = clutime, createnb = nbtime,
        nextSI = SItime, Constime = Contime, feastime = FBtime, Algotime=totaltime)#,
    CSV.write("/home/ak121396/Desktop/GeneralPR1/rc2/"*"$class"*"_iter100.csv",record1, append=true, header=false )#, delim=',' )
    # CSV.write("/home/ak121396/Desktop/GeneralPR1/Y/"*"$fname"*"_Y.log",DataFrame(otable),header=false, delim=' ' )
    # # #########################  Record outputs  ############################
    # CSV.write("/home/ak121396/Desktop/GeneralPR1/Y/"*"$ins"*"_Y.log",DataFrame(otable),header=false, delim=' ' )

end

init
new
cluster
neibour
nextSI
fbcheck
total
