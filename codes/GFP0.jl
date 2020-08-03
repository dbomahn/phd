using DelimitedFiles,DataFrames,JuMP,CPLEX,LinearAlgebra,CPUTime,CSV,Statistics #,PlotlyJS,RDatasets,Colors
global TimeLimit = 7200
################################  Data  ####################################
data=readdlm(ARGS[1])#, '\t', String, '\n')
presol=readdlm(ARGS[2])#, '\t', String, '\n')
pf=readdlm(ARGS[3])[:,2:end] #, '\t', String, '\n')

print("Importing files...\n")
var = round.(presol, digits=2)
objs = round.(Int,pf) #pfs[10*(fi-1)+sf]

# facility_i, customer_j # for FLP
global i,j = filter!(a->typeof(a)==Int, data[1,:])
(global fixcost,capa,cost,cost2 = [],[],[],[] )
for k=2:i+1
    fix,cp = filter!(a->typeof(a)==Int, data[k,:])
    push!(fixcost,fix)
    push!(capa,cp)
end
global demand = filter!(a->typeof(a)==Int, data[2+i,:])
a,b = size(data)
for k=i+3:a
    ct = filter!(a->typeof(a)==Int, data[k,:])
    push!(cost2,ct)
    for l=1:length(ct)
        push!(cost,ct[l])
    end
end
s,ss = size(var)
indx= findall(x->x!=0,[sum(var[i,:]) for i=1:s]) #NZ sol
lpX = [var[i,:] for i in indx]

pf = hcat(objs[:,2],objs[:,1],objs[:,3])
L = reshape([pf[i,:] for i in indx],length(lpX),1)
L2 = transpose(hcat(L...))
fmin = Array{Int}([minimum(L2[:,i]) for i=1:3])
fmax = Array{Int}([maximum(L2[:,i]) for i=1:3])
steps = [round.(Int,abs(fmax[i]-fmin[i])/length(L)*j) for i=1:3] #determine steps according to #customers(j)
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

    result = true
    if all(x-> 0==x, txi)
        result = false
    else
        #1st constraint checking
        for b=i+1:j+i
            if sum(txi[(a)*j+b] for a=1:i) == txi[b]
                continue
            else
                result = false
                return result
            end
        end

        #2nd constraint checking
        # if result == true
        for a=1:i
            for b=1:j
                if txi[i+j+j*(a-1)+b]<=txi[a]
                    continue
                else
                    result = false
                    return result
                end
            end
        end
        # end
    end
    return result
end

function fbsearch(idx,steps,fmin,txi,fmax) #solveLP
    yub = [idx[j]*steps[j]+fmin[j] for j=1:3]
    ylb = [(idx-[1,1,1])[j]*steps[j]+fmin[j] for j=1:3]
    idx0 = findall(x->x==0, txi)
    idx1 = findall(x->x==1, txi)
    lp = Model(with_optimizer(CPLEX.Optimizer))
    MOI.set(lp, MOI.RawParameter("CPX_PARAM_SCRIND"), false)
    MOI.set(lp, MOI.RawParameter("CPX_PARAM_THREADS"),1  )
    @variables(lp, begin
        0 <= x[1:i+j+i*j] <=1
    end)
    @objective( lp, Min, sum(x[n] for n in idx0) + sum(1-(x[m]) for m in idx1) )
    @constraint(lp, con1[b=1:j], sum(x[i+j+j*(a-1)+b] for a in 1:i) == x[i+b])
    @constraint(lp, con2[a=1:i,b=1:j], x[i+j+j*(a-1)+b]<=x[a])
    @constraint( lp, c1,  ylb[2]+1 <=sum(Array{Float64,1}(fixcost)[a]*x[a] for a=1:i) <= yub[2]-1 ) #  fmax[2])#
    @constraint( lp, c2,  ylb[1]+1 <=-dot(transpose(Array{Float64,1}(demand)),x[i+1:j+i]) <= yub[1]-1 ) #  fmax[1])#
    @constraint( lp, c3,  ylb[3]+1 <=dot(transpose(Array{Float64,1}(cost)),x[i+j+1:end]) <= yub[3]-1 ) #fmax[3]) #
    optimize!(lp)
    if termination_status(lp) == MOI.OPTIMAL
        return JuMP.value.(x)
    else
        return 0; #print("no new lp solution found");
    end
end

function mindist(sameval,groups,groupkeys,L2,l,que) #find a point which has the minimum distance to the avg value of points in the same cuboid
    # yub = [groupkeys[l][j]*steps[j]+fmin[j] for j=1:3]
    # ylb = [(groupkeys[l]-[1,1,1])[j]*steps[j]+fmin[j] for j=1:3]

    dis = Model(with_optimizer(CPLEX.Optimizer))
    MOI.set(dis, MOI.RawParameter("CPX_PARAM_SCRIND"), false)
    @variable(dis, 0<=x[1:i+j+i*j]<=1)
    @constraint(dis, con[k in sameval], x[k] == que[k] )
    cubemean = [mean(L2[groups[l],:][:,i]) for i=1:3]

    @constraint(dis, dot(x[1:i],transpose(Array{Float64}(fixcost))) == cubemean[2]) #yub[2]
    @constraint(dis, dot(x[i+1:j+i],transpose(Array{Float64}(demand))) == -cubemean[1]) #transpose(Array{Float64}
    @constraint(dis, dot(x[i+j+1:end],transpose(Array{Float64}(cost)))== cubemean[3]) #transpose(Array{Float64}
    @constraint(dis, con1[b=1:j], sum(x[i+j+j*(a-1)+b] for a in 1:i) == x[i+b])
    @constraint(dis, con2[a=1:i,b=1:j], x[i+j+j*(a-1)+b]<=x[a])

    @objective(dis, Min, 1)
    optimize!(dis)
    if JuMP.termination_status(dis) ==MOI.OPTIMAL
        newlpsol = JuMP.value.(x)
        return newlpsol
    else
        return nothing
    end #solveLP
end

function Grouping(groups,groupkeys,lpX,steps,fmin,L2) #find common elements(variable values) from solutions in the same cuboid
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
            sameval = findall(x->x==true,[all(x->x==que[l][1], que[l,:]) for l=1:length(lpX)])
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
    return [dot(x[j+i+1:end],cost[:]),dot(x[1:i],fixcost[:]),-dot(x[i+1:j+i],demand[:]),]
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
Grouping(groups[1:min(length(groups),5)],groupkeys[1:min(length(groups),5)],lpX,steps,fmin,L2)
########################################################
GroupingTime = @CPUelapsed Xf,cand,loca_check,LPcount = Grouping(groups,groupkeys,lpX,steps,fmin,L2)
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
        Max_itr = j+i #max(count(x->0<x<1,x_t),1)  #maximum number of attempts => How to set
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
push!(Pobj,[0,0,0])
push!(P,zeros(Int,i*j+i+j))
#Filter dominated solutions
FPsol, FPPF = domFilter(P,Pobj)
############## remove the location of dominated solutions #########
cpcoordi = Dict();
for i=1:length(Pcoordi)
    cpcoordi[i] = Pcoordi[i]
end
for i=1:length(Pobj)-1
    for j=i+1:length(Pobj)
        if all(Pobj[i] .>= Pobj[j]) == true #dominated by PF[j]
            cpcoordi[i]=nothing;break
        elseif all(Pobj[j] .>= Pobj[i]) == true
            cpcoordi[j]=nothing;
        end
    end
end
FPcoordi = filter!(a->a!=nothing, collect(values(cpcoordi)))

detectlex = [count(i->i==0,FPPF[k]) for k=1:length(FPPF)]
lexid = findall(x->x>=2,detectlex)
for i=1:length(lexid)
    deleteat!(FPsol,lexid[i]); deleteat!(FPPF,lexid[i]); deleteat!(FPcoordi,lexid[i])
end

Pz,Py,Px=[hcat(FPPF...)[i,:] for i=1:3]; dfP = Px,Py,Pz
matP = zeros(length(FPsol),i+j+i*j)
for i=1:length(FPsol)
    matP[i,:] = FPsol[i]
end

#########################  Record outputs  ############################
ins = ARGS[1][end-12:end-7] #CPUtime recorded, naming after the instance
colname = ARGS[1][end-12:end-4]
record1 = DataFrame( Initial_IPsol = length(Xf1), Initial_LPsol = length(candX),LP=FPLPcount, sol=length(Pz), removed=length(P)-length(Pz))
insertcols!(record1,6, Symbol("$colname")=>GroupingTime+FPTime)
GFPX=DataFrame(matP); GFPY=DataFrame(dfP);

CSV.write("/home/k2g00/k2g3475/multiobjective/solvers/ep+FP/GFPresults/GFPlex_"*"$ins"*"_record.csv",record1, append=true, writeheader=false )#, delim=',' )
#CSV.write(ARGS[1]*"_GFP_X.csv",GFPX, append=true, writeheader=false)
CSV.write(ARGS[1]*"_GFP_Y_.csv",GFPY, writeheader=false, delim=' ' )
print(colname," GroupFP Done!")
