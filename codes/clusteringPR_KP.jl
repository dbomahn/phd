using DataStructures,DelimitedFiles,DataFrames,LinearAlgebra,Clustering,StatsBase,Random,CPUTime,CSV
mutable struct Data
    xval::String; yval::String; dtfile::String;
    dvar::Array{}; LB::Array{}
    n::Int; C::Array{}; ub::Int; weight::Array{}
    radius::Float64
    function Data(xval::String,yval::String,dtfile::String)
        dv = round.(digits=4, readdlm(xval))
        objs = readdlm(yval)[:,2:end]
        ind = findall(i-> -1 in objs[i,:]  , 1:size(objs)[1])
        dvar = dv[setdiff(1:end, ind), :]; LB = objs[setdiff(1:end, ind), :];
        maxobj = [maximum(LB[:,i]) for i=1:3]; minobj = [minimum(LB[:,i]) for i=1:3];
        steps = round.(Int,(maxobj-minobj)/length(LB[:,1]))
        interval = round(Int,mean(maxobj-minobj)/length(LB[:,1]))
        # points = transpose(LB)
        α = max(round(length(LB[:,1])*0.05),1)
        radius = interval*α

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
        new(xval,yval,dtfile,dvar,LB,n,C,ub,weight,radius)
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
            clumean[i] = mean.(LB[core[i],:][:,k] for k=1:3)
        else
            clumean[i] = mean.(LB[boundary[i],:][:,k] for k=1:3)
        end
    end
    cores = collect(values(SortedDict(clumean)))
    return LBclu,clsize,clulist,cores,candI,candG
end

function KPfbcheck(x,n,weight,ub)
    result = false
    if all(i->round(i)==i,x) && dot(weight,x) <= ub
        # print("total weight: ", dot(weight,x),"outof ",ub, "\n")
        result = true
    else
        result = false
    end
    return result
end

function createNB(SI,C,dif,exploredSI,explSIobj)
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
            push!(neiobj,[dot(C[k,:],cpSI) for k=1:3])
        else #if vari is fractional
            cpSI[i] = 1; push!(neibour, cpSI); push!(neiobj,[dot(C[k,:],cpSI) for k=1:3])
            cpSI = copy(SI)
            cpSI[i] = 0; push!(neibour, cpSI); push!(neiobj,[dot(C[k,:],cpSI) for k=1:3])
        end
    end
    neibour = setdiff(neibour,exploredSI) #removed already explored SI
    neiobj = setdiff(neiobj,explSIobj)
    return neibour,neiobj
end
function domFilter(sol,obj)
    copysol = Dict(); copyobj = Dict();
    for i=1:length(obj)
        copysol[i] = sol[i]
        copyobj[i] = obj[i]
    end

    for i=1:length(obj)-1
        for j=i+1:length(obj)
            if all(obj[i] .<= obj[j]) == true #dominated by PF[j]
                copyobj[i]=nothing; copysol[i]=nothing;break
            elseif all(obj[j] .<= obj[i]) == true
                copyobj[j]=nothing; copysol[j]=nothing;
            end
        end
    end

    finalsol = filter!(a->a!=nothing, collect(values(copysol)))
    finalobj = filter!(a->a!=nothing, collect(values(copyobj)))

    return finalsol,finalobj
end
function nextSI(neibour,neiobj,C,SI,candI,cores,clulist,clsize,LBclu,dvar)
    closclu = []
    for k=1:length(neiobj)
        newsol = neiobj[k,:][1]
        gap = [cores[i]-newsol for i=1:length(cores)]
        dis = [sum(gap[i].^2) for i=1:length(gap)]
        clos = findall( x-> x==minimum(dis), dis )[1]
        push!(closclu,clos)
    end

    SIobj = [dot(C[k,:],SI) for k=1:3]
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
        else #there is no candid
            k = rand(1:length(neibour)) #randomly select one candidate

        end
    end
    return neibour[k],neiobj[k]
end
######################## Path Relinking  ###################################
# LB = kp.LB; dvar=kp.dvar; n=kp.n; weight = kp.weight; ub=kp.ub; C=kp.C; n=kp.n; radius = kp.radius

# runtime1 = @CPUelapsed candset,candobj = clusterPR(n,C,weight,ub,LB,dvar)

# h = 0
# for i=1:3
#     h = h+@CPUelapsed 1+1
#     # return h
# end
# h
# 1
function clusterPR(n,C,weight,ub,LB,dvar,radius)
    candset = []; candobj=[]; IGPair=[];
    clutime = 0; nbtime = 0; SItime = 0;
    for i=1:200#size(LB)[1]*10)
        clutime = clutime + @CPUelapsed LBclu,clsize,clulist,cores,candI,candG = clustering(LB,radius)
        I = rand(LBclu[rand(candI)]); G = rand(LBclu[rand(candG)])
        SI = dvar[I,:]; SG = dvar[G,:]; iter=1; exploredSI = []; explSIobj = [];

        while all.(SI != SG) && [I,G]∉IGPair && iter<n*10
            # @show iter
            rg = range(1, length=n)
            dif = findall(i-> SI[i]!=SG[i], rg)
            nbtime = nbtime + @CPUelapsed neibour,neiobj = createNB(SI,C,dif,exploredSI,explSIobj)
            if length(neibour) == 0
                break
            end
            # print("# neighbours:",length(unique(neibour)),"\n")

            # find out feasible sol in neighbour and add to LB set
            for i=1:length(neibour)
                if KPfbcheck(neibour[i],n,weight,ub)==true && neibour[i]∉candset
                    # print("feasi sol in neighbours \n")
                    dvar = [dvar; transpose(neibour[i])]
                    LB = [LB; transpose(neiobj[i])]
                end
            end
            SItime = SItime + @CPUelapsed SI,SIobj = nextSI(neibour,neiobj,C,SI,candI,cores,clulist,clsize,LBclu,dvar)
            push!(exploredSI,SI); push!(explSIobj,SIobj)
            # Feasibilitycheck
            if KPfbcheck(SI,n,weight,ub)==true && SI∉candset
                push!(candset,SI)
                push!(candobj,SIobj)

                # print("feasible new sol added \n")
            end
            iter+=1
        end
        push!(IGPair,[I,G])
    end
    return candset,candobj,clutime,nbtime,SItime
end

function PostProc(P,C,cand,candobj,ub,n,weight)
    Pobj = []
    for i=1:size(P)[1]
        push!(Pobj,[dot(C[k,:],P[i,:]) for k=1:3])
    end
    newsol = 0
    for i=1:length(cand)
        if dominated(candobj[i],Pobj)==false && candobj[i]∉Pobj
            P = [P; transpose(cand[i][:])]
            push!(Pobj,candobj[i][:])
            newsol+=1
        end
    end
    return P,Pobj,newsol
end

paths = ("/home/ak121396/Desktop/Bensolve_KP/X/", "/home/ak121396/Desktop/Bensolve_KP/Y/", "/home/ak121396/Desktop/instances/KP_for_Bensolve/")
xfile = readdir(paths[1]); yfile = readdir(paths[2]); ins = readdir(paths[3])
for i=1:10
    kp = Data(paths[1]*xfile[i],paths[2]*yfile[i],paths[3]*ins[i])
    runtime1 = @CPUelapsed candX,candobj,clutime,nbtime,SItime = clusterPR(kp.n,kp.C,kp.weight,kp.ub,kp.LB,kp.dvar,kp.radius)
    print("#sol Befor filtering : ", length(candobj),"\n")
    runtime2 = @CPUelapsed P,Pobj,newsol = PostProc(kp.P,kp.C,cand,candobj,kp.ub,kp.n,kp.weight)
    for i=1:length(Pobj)
        for j=1:3
            otable[i,j] = Pobj[i][j]
        end
    end
    # runtime2 = @CPUelapsed effsol,ndpoints, = domFilter(candX,candobj); print("# final sols : ",length(effsol),"\n")
    # otable = ones(Int, length(effsol),3);
    # for i=1:length(effsol)
    #     for j=1:3
    #         otable[i,j] = -ndpoints[i][j]
    #     end
    # end
    fname = xfile[i][1:end-6]
    # record1 = DataFrame(sol=length(ndpoints),CPUtime=runtime1+runtime2)#newsol=effsol,
    record1 = DataFrame(clustering = clutime, createnb = nbtime,nextSI = SItime,Algo=runtime1)
    CSV.write("/home/ak121396/Desktop/clusterPR/record/iter200_functime.csv",record1, append=true, header=false )#, delim=',' )
    # CSV.write("/home/ak121396/Desktop/clusterPR/F500Wn/"*"$fname"*"_Y.log",DataFrame(otable),header=false, delim=' ' )
    # # #########################  Record outputs  ############################
end

############################  for Visualisation  #############################
points = DataFrame(obj1=LB[:,1], obj2=LB[:,2], obj3=LB[:,3])
# IB = readdlm("/home/ak121396/Desktop/Bensolve_KP/RoundDown/Y/"*readdir("/home/ak121396/Desktop/Bensolve_KP/RoundDown/Y/")[i])
# points = DataFrame(obj1=-IB[:,1], obj2=-IB[:,2], obj3=-IB[:,3])

clusort = Dict()
for i=1:length(LBclu)
    for j=1:length(LBclu[i])
        clusort[LBclu[i][j]] = i
    end
end

clusid = []
for i=1:length(clusort)
    push!(clusid,clusort[i])
end
insertcols!(points, 4, :cluster => clusid)
clustering(points)

######################## clustering PR&FP ###########################
initclu = dbscan(LB', radius) #, min_neighbors = 1, min_cluster_size = 1)
core = getproperty.(initclu, :core_indices)
boundary = getproperty.(initclu, :boundary_indices)
LBclu = union.(core,boundary); Icluster = Dict();
clsize = getproperty.(initclu, :size)
clulist = unique(sort(clsize))
