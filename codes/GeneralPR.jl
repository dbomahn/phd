using DataStructures,DelimitedFiles,DataFrames,JuMP,CPLEX,LinearAlgebra,CPUTime,CSV,Statistics,Clustering,MathOptInterface,StatsBase,JLD2,SparseArrays

mutable struct Data
    input::String
    i::Int; j::Int; C::Array{}

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
function clustering(obj,radius)
    initclu = dbscan(obj', radius) #, min_neighbors = 1, min_cluster_size = 1)
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
            clumean[i] = round.(digits=2, mean.(obj[core[i],:][:,k] for k=1:3))
        else
            clumean[i] = round.(digits=2, mean.(obj[boundary[i],:][:,k] for k=1:3))
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
            return neibour[candid[1]]#,neiobj[candid[1]]
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
            mostimp = findall(x-> x == maximum(ranksum), ranksum)
            k = candid[rand(mostimp)[1]]
            return neibour[k]#, neiobj[k]
        else #there is no candidate
            # print("randomly select a candidate \n")
            k = rand(1:length(neibour)) #randomly select a candidate
            return neibour[k]#,neiobj[k]
        end
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
function Postpro(dvar,LB,newsol)
    #Filter fractional solutions from LB
    initdv = dvar[1:end-newsol]
    initLB = LB[1:end-newsol]
    frac = findall(j-> trunc.(initdv[j])!= initdv[j], 1:length(initdv))
    dv2 = initdv[setdiff(1:end,frac)]; LB2 = initLB[setdiff(1:end,frac)]
    P = union(dv2,dvar[end-newsol+1:end])
    Pobj = union(LB2,LB[end-newsol+1:end])

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
function clusterPR(di,dj,C,dvar,LB,LBmtx,radius,TL,n)
    newsol=0; IGPair=[];
    clutime = 0; nbtime = 0; SItime = 0; FBtime = 0; t0=time();
    clutime= clutime + @CPUelapsed LBclu,clsize,clulist,cores,candI,candG = clustering(LBmtx,radius)

    for i=1:dj
        if time()-t0 >= TL
            break
        end
        I = rand(LBclu[rand(candI)]); G = rand(LBclu[rand(candG)])
        SI = dvar[I]; SG = dvar[G]; iter=1; exploredSI = [];

        while all.(SI != SG) && [I,G]∉IGPair && iter<20 && (time()-t0<TL)
            rg = range(1, length=n)
            dif = findall(i-> SI[i]!=SG[i], rg)
            nbtime = nbtime + @CPUelapsed neibour,neiobj = createNB(di,dj,SI,C,dif,exploredSI)

            for l=1:length(neiobj)
                FBtime = FBtime + @CPUelapsed feasi = FBcheck(neibour[l],n)
                if feasi == true && neibour[l]∉ dvar
                    push!(dvar, neibour[l]); push!(LB, neiobj[l]);
                    newsol+=1
                end
            end
            SItime = SItime + @CPUelapsed SI = nextSI(di,dj,neibour,neiobj,C,SI,candI,cores,clulist,clsize,LBclu,dvar)

            if SI∉dvar; push!(exploredSI,SI); end
            iter+=1
        end
        push!(IGPair,[I,G])
    end
    return dvar,LB,clutime,nbtime,SItime,FBtime,newsol
end

###########################      GPR without clustering   ######################################################

function nextSI2(di,dj,neibour,neiobj,C,SI)
    SIobj = getobjval(SI,C,di,dj)
    for i=1:length(neiobj)
        if length(neiobj) == 1  #if there is one candiate sol
            return neibour[1]#,neiobj[1]
        elseif length(neiobj) > 1 # if there are multiple candiates, check the improved ratio
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
        else #there is no candidate
            # k = rand(1:length(neibour)) #randomly select a candiate
            # return neibour[k],neiobj[k]
        end
    end
end

function GPR(di,dj,C,n,dvar,LB,TL)
    IGPair=[]; exploredSI = []; newsol=0;
    nbtime = 0; SItime = 0; FBtime = 0; t0=time();

    for i=1:dj
    	if time()-t0 >= TL
            break
        end
        I,G = sample(1:length(dvar), 2, replace=false)
        SI = dvar[I]; SG = dvar[G]; iter=0;

        while all.(SI != SG) && [I,G]∉IGPair && iter<20 && (time()-t0<TL)
            # rg = range(1, length=n)
            dif = findall(i-> SI[i]!=SG[i], 1:n); #1:n
            nbtime = nbtime + @CPUelapsed neibour,neiobj = createNB(di,dj,SI,C,dif,exploredSI)
            # print("neighbours: ",length(neiobj),"\n")

            if length(neiobj)==0
                break
            else
                for l=1:length(neiobj)
                    FBtime = FBtime + @CPUelapsed feasi = FBcheck(neibour[l],n)
                    if feasi == true && neibour[l]∉ dvar
                        push!(dvar, neibour[l]); push!(LB, neiobj[l]);
                        newsol+=1
                    end
                end
                SItime = SItime + @CPUelapsed SI = nextSI2(di,dj,neibour,neiobj,C,SI)
                if SI∉dvar
                    push!(exploredSI,SI);
                end
            end
            iter+=1
        end
        push!(IGPair,[I,G])
    end
    return dvar,LB,nbtime,SItime,FBtime,newsol
end
# FLP
paths = ("/home/ak121396/Desktop//instances/FLP/instances/","/home/ak121396/Desktop//instances/FLP/varval/","/home/ak121396/Desktop//instances/FLP/PF/")
data = Data( paths[1]*readdir(paths[1])[3] )
pre = Val(paths[2]*readdir(paths[2])[3],paths[3]*readdir(paths[3])[3],data.j)
flp = direct_model(CPLEX.Optimizer())
MOI.set(flp, MOI.RawParameter("CPX_PARAM_SCRIND"), false );
MOI.set(flp, MOI.RawParameter("CPX_PARAM_THREADS"),1  );
@variable(flp, x[1:data.i+data.j+(data.i*data.j)],Bin );
@constraint(flp, con1[b=1:data.j], sum(x[data.i+data.j+data.j*(a-1)+b] for a in 1:data.i) == x[data.i+b]);
@constraint(flp, con2[a=1:data.i,b=1:data.j], x[data.i+data.j+data.j*(a-1)+b] <= x[a]);
optimize!(flp)
if data.i<=15
    TL = 60;
elseif data.i<=25
    TL = 180;
else
    TL = 1800;
end

# GPR(data.i,data.j,data.C,pre.dvar,pre.LB,data.i+data.j+(data.i*data.j),TL)
runtime1 = @CPUelapsed candset,candobj,clutime,nbtime,SItime,FBtime,newsol = clusterPR(data.i,data.j,data.C,pre.dvar,pre.LB,pre.LBmtx,pre.radius,TL,n)
runtime2 = @CPUelapsed finalX,finalY = Postpro(candset,candobj,newsol)
push!(finalY,[0,0,0])
# totaltime = runtime1 + runtime2
# finalX = push!(finalX, [zeros(Int,n)'])

for i=1:length(finalX)
    if FBcheck(finalX[i],65) == true
    # all(i-> (i==0||i==1), finalX[i] )
        print("y \n")
    end
end


otable = zeros(Int, length(finalY),3)
for i=1:length(finalY)
    for j=1:3
        otable[i,j] = finalY[i][j]
    end
end
matriX = zeros(Int,length(finalX),n)
for i=1:length(finalX)
    for j=1:n
        matriX[i,j] = finalX[i][j]
    end
end
sparX = sparse(matriX); #JLD2.@save "/home/k2g00/k2g3475/multiobjective/solvers/generalPR/outputs/X/"*"$ins"*"_X.jdl2" sparX
JLD2.@save "/home/ak121396/Desktop/PP.jdl2" sparX
JLD2.@load "/home/ak121396/Desktop/PP.jdl2" sparX





########################################## Call dv values ###################
dir2 = "/home/ak121396/Desktop/GeneralPR/cluoutputs/X/"
dvar = readdir(dir2)
for i=1:length(dvar)
    JLD2.@load dir2*dvar[i] sparX
    dvmatrix = Array(sparX)

    if all(i-> (i==0||i==1), dvmatrix ) == false
        @show i
    end
end

if any(dvmatrix[i,:])
    @show i
else
end


all(i-> (i==0||i==1), dvmatrix )
################################# runs on the local PC #####################
for i=1:1#20
    data = Data( paths[1]*readdir(paths[1])[i] )
    pre = Val(paths[2]*readdir(paths[2])[i],paths[3]*readdir(paths[3])[i],data.j)

    # flp = Model(CPLEX.Optimizer)
    flp = direct_model(CPLEX.Optimizer())
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


# function oldnextSI2(di,dj,neibour,neiobj,C,SI,dvar)
#     SIobj = getobjval(SI,C,di,dj)
#      #check if intermediate path is in the smallest cluster
#     for i=1:length(neiobj)
#         if length(neiobj) == 1  #if there is one candiate sol
#             return neibour[1],neiobj[1]
#         elseif length(neiobj) > 1 # if there are multiple candiates, check the improved ratio
#             ratiotb = zeros(length(neiobj),3)
#             for i=1:length(neiobj)
#                 ratiotb[i,:] = neiobj[i]./SIobj
#             end
#             ranktb = zeros(length(neiobj),3)
#             for i=1:3
#                 ranktb[:,i] = tiedrank(ratiotb[:,i])
#             end
#             ranksum = [sum(ranktb[i,:]) for i=1:length(neiobj)]
#             mostimp = findall(x-> x == maximum(ranksum), ranksum)
#             k = rand(mostimp)[1]
#             return neibour[k], neiobj[k]
#         else #there is no candidate
#             # k = rand(1:length(neibour)) #randomly select a candiate
#             # return neibour[k],neiobj[k]
#         end
#     end
# end
# function oldGPR(di,dj,C,dvar,LB,n,TL)
#     candX = []; candY=[]; IGPair=[]; cpX = [vec(pre.dvar[i,:]) for i=1:size(pre.dvar)[1]]; exploredSI = [];
#     nbtime = 0; SItime = 0; Contime=0; FBtime = 0; t0=time();
#
#     for i=1:dj
#         if time()-t0 >= TL
#             break
#         end
#         I,G = sample(1:length(cpX), 2, replace=false)
#         SI = cpX[I]; SG = cpX[G]; iter=1;
#
#         while all.(SI != SG) && [I,G]∉IGPair && iter<20 && time()-t0<TL
#
#             rg = range(1, length=n)
#             dif = findall(i-> SI[i]!=SG[i], rg);
#             nbtime = nbtime + @CPUelapsed neibour,neiobj = createNB(di,dj,SI,C,dif,exploredSI)
#             # print("neighbours: ",length(neiobj),"\n")
#
#             if length(neiobj)==0
#                 break
#             else
#                 for l=1:length(neiobj)
#                     Contime = Contime + @CPUelapsed constcheck = ConstCheck(neibour[l],di,dj)
#                     if constcheck == true && neibour[l]∉ candY
#                         dvar = [dvar; transpose(neibour[l])]
#                         LB = [LB; transpose(neiobj[l])]
#                     end
#                 end
#                 SItime = SItime + @CPUelapsed SI,SIobj = nextSI2(di,dj,neibour,neiobj,C,SI,dvar)
#                 FBtime = FBtime + @CPUelapsed fbcheck = FBcheck(SI,di,dj)
#                 if fbcheck == true && SI∉candX
#                     push!(exploredSI,SI); push!(candX,SI); push!(candY,SIobj);
#                 end
#             end
#             iter+=1
#         end
#         push!(IGPair,[I,G])
#     end
#     return candX,candY,nbtime,SItime,Contime,FBtime
# end
