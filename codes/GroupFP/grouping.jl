using DelimitedFiles,DataFrames,JuMP,CPLEX,LinearAlgebra,CPUTime,CSV,Statistics #,PlotlyJS,RDatasets,Colors
global TimeLimit = 300
################################  Data  ####################################
struct path{S<:String}
    dir1::S; dir2::S; dir3::S
end
pathes=path("/home/ak121396/Desktop/instances/AP/dat/",
    "/home/ak121396/Desktop/solvers/Bensolve/APoutputs/X/",
    "/home/ak121396/Desktop/solvers/Bensolve/APoutputs/Y/")
dir1 = pathes.dir1; dir2 = pathes.dir2; dir3 = pathes.dir3
apf=readdir(dir1); apx = readdir(dir2);apy = readdir(dir3)
dt = Data(dir1*apf[14]);x = dir2*apx[14];y = dir3*apy[14]
####################### current Val part #########
dv = round.(digits=4, readdlm(x))
objs = round.(digits=4, readdlm(y))
ind = findall(i-> 0 in objs[i,:]  , 1:size(objs)[1])
dv2 = dv[setdiff(1:end, ind), :];
L = objs[setdiff(1:end, ind), 2:end];
candX = [Vector(dv2[i, :]) for i = 1:size(LBmtx)[1]]
LB = [Vector(LBmtx[i, :]) for i = 1:size(LBmtx)[1]]


##################################clustering####################################

# (clu.centers)' # get the cluster centers # assignments(clu) # get the assignments of points to clusters

k = round(Int,length(LB)/n); clu = kmeans(L',max(k,2));
targec = indexin(nsmallest(2,sc),sc)

candx1 = findall(i->i==targec[1],clu.assignments); candx2 = findall(i->i==targec[2],clu.assignments)
un = Set()
push!(un,Set([4,5]))
Set([4,3]) == Set([3,4])
Set([3,9]) ∉ un
function matheuristic(candX,LB,n,C,TL,clu)
    nc = nclusters(clu);  SolPair = Set();
    while t0-<TL

    end

    sc = counts(clu);
    targec = findall(i->i in nsmallest(2, sc),sc)
    candx1 = findall(i->i==targec[1],clu.assignments); candx2 = findall(i->i==targec[2],clu.assignments)
    id1 = sample(candx1, 1, replace=false)[1]; id2 = sample(candx2, 1, replace=false)[1]
    x1 = candX[id1]; x2 = candX[id2];
    push!(SolPair,Set([id1,id2]))

    while Set([id1,id2]) ∉ SolPair
        if x1 == round.(x1) && x2 == round.(x2) #both sols are int
            fx,fy,fns = weightFP(candX,LB,x1,x2,SolPair,n,C,TL)


#################### Groouping ###############
fmin = [minimum(L2[:,i]) for i=1:3]
fmax = [maximum(L2[:,i]) for i=1:3]
steps = [round.(Int,abs(fmax[i]-fmin[i])/length(L)*n) for i=1:3] #determine steps according to #customers(j)
cube = Dict();

for iter=1:length(L)
    loca = [round.(Int,((L2[iter,k]-fmin[k])/steps[k])+1) for k=1:3]
    if !haskey(cube,loca)
        cube[loca] = [iter]
    else
        push!(cube[loca], iter)
    end
end
################################   Functions  ##################################
########################        for Feasibility Pump      ######################
function flip(hxi,j,e)
    if hxi[e[j]]==1
        hxi[e[j]] = 0
    else
        hxi[e[j]] = 1
    end
    return hxi
end
function flipoper(Tabu,tx,txi)
    e = sortperm(abs.(tx-txi),rev=true)
    xi = []
    hxi = copy(txi)
    j = 1
    M=length(tx) #
    while j<=M && xi==[]
        hxi = flip(hxi,j,e)
        if hxi ∉ Tabu
            xi=hxi
        else
            j+=1
        end
    end

    if xi==[]
        while j<=M
            hxi=copy(txi)
            Num = Int64(rand(ceil(length(txi)/2):length(txi)-1))
            R = sample(1:M,Num, replace=false)
            for i in R
                hxi = flip(hxi,r)
                if hxi ∉ Tabu
                    xi = hxi
                end
            end
            j+=1
        end
    end
    return xi
end
function fbcheck(txi)
    fb = Model(CPLEX.Optimizer)
    MOI.set(fb, MOI.RawParameter("CPX_PARAM_SCRIND"), false)
    MOI.set(fb, MOI.RawParameter("CPX_PARAM_THREADS"),1  )
    @variable(fb, x[1:n*n], Bin)
    @objective(fb, Min, 1)
    for i=1:2*n
        @constraint( fb,  sum(B[i,j]*x[j] for j=1:n*n) == 1 )
    end
    optimize!(fb)
    if termination_status(fb) == MOI.OPTIMAL
        return true
    else
        return false
    end
end
function fbsearch(idx,steps,fmin,txi,fmax) #solveLP
    yub = [idx[j]*steps[j]+fmin[j] for j=1:3]
    ylb = [(idx-[1,1,1])[j]*steps[j]+fmin[j] for j=1:3]
    idx0 = findall(x->x==0, txi)
    idx1 = findall(x->x==1, txi)
    ap_m = Model(CPLEX.Optimizer)
    MOI.set(ap_m, MOI.RawParameter("CPX_PARAM_SCRIND"), false)
    MOI.set(ap_m, MOI.RawParameter("CPX_PARAM_THREADS"),1  )
    @variables(ap_m, begin
        0 <= x[1:n*n] <=1
    end)
    @objective( ap_m, Min, sum(x[n] for n in idx0) + sum(1-(x[m]) for m in idx1) )
    for i=1:2*n
        @constraint( ap_m,  sum(B[i,j]*x[j] for j=1:n*n) == 1 )
    end
    @constraint( ap_m, c1, sum(P2[1,j]*x[j] for j=1:n*n) <= yub[1]-1 )#ylb[1]+1 <=
    @constraint( ap_m, c2, sum(P2[2,j]*x[j] for j=1:n*n) <= yub[2]-1 ) #ylb[2]+1 <=
    @constraint( ap_m, c3, sum(P2[3,j]*x[j] for j=1:n*n) <= yub[3]-1 ) #ylb[3]+1 <=
    optimize!(ap_m)
    if termination_status(ap_m) == MOI.OPTIMAL
        return JuMP.value.(x)
    else
        return 0; #print("no new lp solution found");
    end
end
function mindist(sameval,groups,groupkeys,L2,l,que) #find a point which has the minimum distance to the avg value of points in the same cuboid
    # yub = [groupkeys[l][j]*steps[j]+fmin[j] for j=1:3]
    # ylb = [(groupkeys[l]-[1,1,1])[j]*steps[j]+fmin[j] for j=1:3]

    dis = Model(CPLEX.Optimizer)
    MOI.set(dis, MOI.RawParameter("CPX_PARAM_SCRIND"), false)
    @variable(dis, 0<=x[1:n*n]<=1)
    @constraint(dis, con[k in sameval], x[k] == que[k] )
    cubemean = [mean(L2[groups[l],:][:,i]) for i=1:3]
    @constraint(dis, sum(P2[1,j]*x[j] for j=1:n*n)  == cubemean[1])
    @constraint(dis, sum(P2[2,j]*x[j] for j=1:n*n) == cubemean[2])
    @constraint(dis, sum(P2[3,j]*x[j] for j=1:n*n) == cubemean[3])
    for i=1:2*n
        @constraint( dis,  sum(B[i,j]*x[j] for j=1:n*n) == 1 )
    end
    @objective(dis, Min, 1)
    optimize!(dis)
    if JuMP.termination_status(dis) ==MOI.OPTIMAL
        newlpsol = JuMP.value.(x)
        return newlpsol
    else
        return nothing
    end #solveLP
end
function Grouping(groups,groupkeys,lpX,steps,fmin,L2,n) #find common elements(variable values) from solutions in the same cuboid
    LPcount1 = 0;
    Xf= Dict(); candX = Dict(); loca_check = [];
    for l=1:length(groups)
        glength = length(groups[l])
        if glength==1 #if there is one point in a group
            if (all.(x->trunc(x)==x, lpX[groups[l]])[1] == 1) #all varval are integer => insert to Xf
                Xf[groupkeys[l]] = lpX[groups[l][1]]
                # print("lp==ip \n")
            else
                candX[groupkeys[l]] = lpX[groups[l]][1] #candX is a set of solutions that should be processed to be an integer solution
                # print("only lp sol in cuboid \n")
            end

        else
            cogroup = lpX[groups[l]]
            que = round.(hcat(cogroup...))
            sameval = findall(x->x==true,[all(x->x==que[l][1], que[l,:]) for l=1:n*n])
            newlpsol = mindist(sameval,groups,groupkeys,L2,l,que)
            LPcount1+=1

            if newlpsol == nothing
                # print("there is no LPsolution \n")
            else #when newlpsol is integer, check if all values are integer
                if (all.(x->trunc(x)==x, newlpsol) == 1) #all varval is integer => insert to Xf
                    Xf[groupkeys[l]] = newlpsol
                    # print("newlpsol==ip \n")
                else
                    push!(loca_check, newlpsol);
                    # print("newlpsol in cuboid is found\n")
                    candX[groupkeys[l]] = newlpsol;
                end
            end

        end
    end
    return Xf,candX,loca_check,LPcount1
end
function getlocation(x,fmin,steps)
    L = giveobjval(x)
    # cubes = []
    # for x=1:length(L)
    loca = [round.(Int,((L[k]-fmin[k])/steps[k])+1) for k=1:3]
    #     push!(cubes,loca)
    # end
    return loca
end
## for tri- epsilon
function giveobjval(x)
    return [sum(P2[1,j]*x[j] for j=1:n*n),sum(P2[2,j]*x[j] for j=1:n*n),sum(P2[3,j]*x[j] for j=1:n*n)]
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

#############################    Group Feasibility Pump   ############################
groups = collect(values(cube)); groupkeys = collect(keys(cube))
########################### TEST ########################
Grouping(groups[1:min(length(groups),5)],groupkeys[1:min(length(groups),5)],lpX,steps,fmin,L2,n*n)
########################################################
GroupingTime = @CPUelapsed Xf,cand,loca_check,LPcount = Grouping(groups,groupkeys,lpX,steps,fmin,L2,n*n)
Xf1 = copy(Xf)
# candX = filter(x->(last(x)!==nothing), candX)# candkeys = collect(keys(candX)); candvals = collect(values(candX));
candX = vcat(collect(values(cand)),loca_check)
PFset = giveobjval.(collect(values(Xf1)))
function GroupFP(Xf,PFset,candX,LPcount,steps,fmin,fmax)
    Tabu = []; elaps=0;
    algo_start = CPUtime_us()
    for k=1:length(candX)
        # print("===============Feasi Pump",k," th candidate sol ==============\n")
        x_t = candX[k];
        SearchDone = false
        itr = 1
        Max_itr = n*10 #max(count(x->0<x<1,x_t),1)  #maximum number of attempts => How to set
        while itr<Max_itr && SearchDone==false
            xi_t = round.(Int,x_t)
            fxi = giveobjval(xi_t)
            if ( (fbcheck(xi_t) == true) && (dominated(fxi,PFset)==false) ) #checking feasibility and dominance   #z[1]<=fmax[1]-1 && z[2]<=fmax[2]-1 && z[3]<=fmax[3]-1)
                idx = getlocation(xi_t,fmin,steps)
                push!(PFset,fxi) #add new solval to PFset
                if !haskey(Xf,idx)
                    # print("ROUNING => another solution added \n")
                    Xf[idx] = xi_t
                else
                    # print("ROUNING =>new intsol to a cuboid \n")
                    push!([Xf[idx]], xi_t)
                end
                SearchDone = true
            else
                if xi_t ∈ Tabu
                    xi_t = flipoper(Tabu,x_t,xi_t)
                    if xi_t==[]
                        # print("FLIP didn't work \n")
                        SearchDone = true
                    else
                        fxi = giveobjval(xi_t)
                        if ( (fbcheck(xi_t) == true) && (dominated(fxi,PFset)==false) )  #z[1]<=fmax[1]-1 && z[2]<=fmax[2]-1 && z[3]<=fmax[3]-1 )
                            idx = getlocation(xi_t,fmin,steps)
                            push!(PFset,fxi) #add new solval to PFset
                            if !haskey(Xf,idx)
                                # print("FLIP=> another solution's added \n")
                                Xf[idx] = xi_t
                            else
                                # print("FLIP=> new intsol to a cuboid \n")
                                push!([Xf[idx]], xi_t)
                            end
                            SearchDone = true
                        end
                    end
                end
                if SearchDone == false
                    push!(Tabu,xi_t) #when break
                    idx = getlocation(xi_t,fmin,steps)
                    x_t = fbsearch(idx,steps,fmin,xi_t,fmax)
                    LPcount+=1
                    if x_t == 0 #when there's no new feasible lp sol
                        # print("no lp sol's found");
                        break
                    # else
                        # print("\n New lp obj: ",giveobjval(x_t), "\n")
                    end
                end
            end
            itr+=1; #@show itr
            if CPUtime_us()-algo_start < TimeLimit*1e6
                continue
            else
                return Xf,PFset,Tabu,LPcount
            end
        end
    end
    return Xf,PFset,Tabu,LPcount
end
#Time Feasibility Pump
FPTime = @CPUelapsed Xf2,PFset2,Tabu,FPLPcount = GroupFP(Xf1,PFset,candX,LPcount,steps,fmin,fmax)
P = convert.(Array{Int,1},collect(values(Xf2)))
Pcoordi = convert.(Array{Int,1},collect(keys(Xf2)))
Pobj = convert.(Array{Int,1}, giveobjval.(P))
#Filter dominated solutions
FPsol, FPPF = domFilter(P,Pobj)
Pz,Py,Px=[hcat(FPPF...)[i,:] for i=1:3]; dfP = Px,Py,Pz
matP = zeros(length(FPsol),n)
for i=1:length(FPsol)
    matP[i,:] = FPsol[i]
end

#########################  Record outputs  ############################
ins = ARGS[1][end-12:end-7] #CPUtime recorded, naming after the instance
colname = ARGS[1][end-12:end-4]
record1 = DataFrame( LP=FPLPcount, sol=length(Pz), CPUtime=GroupingTime+FPTime )
# insertcols!(record1,3, Symbol("$colname")=>GroupingTime+FPTime)
GFPX=DataFrame(matP); GFPY=DataFrame(dfP);
CSV.write("/home/ak121396/Desktop/phd/GFPresults/KP_"*"$ins"*"_record.csv",record1, append=true, writeheader=false )#, delim=',' )
#CSV.write(ARGS[1]*"_GFP_X.csv",GFPX, append=true, writeheader=false)
CSV.write("/home/ak121396/Desktop/phd/GFPresults/Y/"*ARGS[end-12:end-4]*_"KP_Y_.txt",GFPY, writeheader=false, delim=' ' )
print(colname," GroupFP Done!")
