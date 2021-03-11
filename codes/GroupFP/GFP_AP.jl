using DelimitedFiles,DataFrames,JuMP,CPLEX,LinearAlgebra,CPUTime,CSV,Statistics #,PlotlyJS,RDatasets,Colors
global TimeLimit = 300
################################  Data  ####################################
struct path{S<:String}
    dir1::S; dir2::S; dir3::S
end
pathes=path("/home/ak121396/Desktop/AP/data/", "/home/ak121396/Desktop/AP/X/","/home/ak121396/Desktop/AP/Y/")
dir1 = pathes.dir1; dir2 = pathes.dir2; dir3 = pathes.dir3
data = readdir(dir1)
x = readdir(dir2)
y = readdir(dir3)

f=readdlm(dir1*data[i], '\t', String, '\n')
var = round.(readdlm(dir2*x[i]), digits=4) #round numbers till 2 digits
objs = round.(readdlm(dir3*y[i])[:,2:end], digits=4)

# fl = readdlm(ARGS[1], '\t', String, '\n')
# var = round.(readdlm(ARGS[2]), digits=4)
# objs = round.(readdlm(ARGS[3])[:,2:end] , digits=4)

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

P2 = ones(obj,n*n)
P2 = round.(Int,P2)

global ct=0;
for x=1:length(P)
    for y=1:n
        p=parse(Int64,P[x][y])
        idx = Int(floor((x-1)/n))
        ind = ((x-1)*n+y)%(n*n)
        if ind != 0
            P2[idx+1,((x-1)*n+y)%(n*n)] = p
        else
            P2[idx+1,(n*n)] = p
        end
        if p==0
            ct =ct+1;
        end
    end
end
######################## coefficient matrix (B) #############################
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
#############################
s,ss = size(var)
indx= findall(x->x!=0,[sum(var[i,:]) for i=1:s]) #NZ sol
lpX = [var[i,:] for i in indx]

pf = hcat(objs[:,1],objs[:,2],objs[:,3])
L = reshape([pf[i,:] for i in indx],length(lpX),1)
L2 = transpose(hcat(L...))
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
