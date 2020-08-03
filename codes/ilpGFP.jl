using DelimitedFiles, CPLEX, LinearAlgebra, JuMP, CSV, DataFrames, Statistics, CPUTime
# using GLPK (faster than CPLEX)
##########################  Read input files #########################
# path = "/home/ak121396/Desktop/ILP/"
# fff = readdir(path*"data/")
# f = readdlm(path*"data/ILP_p-3_n-30_m-15_ins-2.dat", '\t', String, '\n')
f=readdlm(ARGS[1], '\t', String, '\n')
global nobj=parse(Int, f[1])
global n=parse(Int,f[2])
global m=parse(Int,f[3])
presol = readdlm(ARGS[2])
vari = round.(presol, digits=2)
objs = readdlm(ARGS[3])
########################    fbsearch instances    #####################
##################objective function coefficient (P) matrix ################
c = f[4:nobj+3]
P= ones(nobj,n)
P=round.(Int,P)

global zp=0;
for x=1:length(c)
    a = split(c[x],('[',']',','))
    aa = filter!(e->e ∉ [ "" ,"[","," ,"]"] ,a)
    for y = 1:length(aa)
        p = parse(Int64,aa[y])
        P[x,y] = p
        if p==0
            global zp = zp+1;
        end
    end
end
####### technical coefficient (a) ###########
TC = f[nobj+4:length(f)-1]
global za=0;
a= ones(m,n)
a=round.(Int,a)
for x=1:length(TC)
    t = split(TC[x],('[',']',','))
    tt = filter!(e->e ∉ [ "" ,"[","," ,"]"] ,t)
    for y=1:length(tt)
        tc=parse(Int64,tt[y])
        a[x,y] = tc
        if tc==0
            global za = za+1;
        end
    end
end
####### RHS values (b) ##########
b = ones(1,m)
b = round.(Int,b)
r = f[length(f)]
r1 = split(r, ('[',']',','))
r2 = filter!(e->e ∉ ["" ,"[", "]"] ,r1)
for i=1:m
    b[i] = parse(Int64,r2[i])
end



global i,j = n,m

s,ss = size(vari)
indx= findall(x->x!=0,[sum(vari[i,:]) for i=1:s]) #NZ sol
lpX = [vari[i,:] for i in indx]

pf = hcat(objs[:,2],objs[:,3],objs[:,4])
L = reshape([pf[i,:] for i in indx],length(lpX),1)
L2 = transpose(hcat(L...))
fmin = ([minimum(L2[:,i]) for i=1:3])
fmax = ([maximum(L2[:,i]) for i=1:3])
steps = [round.(Int,abs(fmax[i]-fmin[i])/(length(L)*j)) for i=1:3] #determine steps according to #customers(j)
cube = Dict();

for iter=1:length(L)
    loca = [round.(Int,((L2[iter,k]-fmin[k])/steps[k])+1) for k=1:3]
    if !haskey(cube,loca)
        cube[loca] = [iter]
    else
        push!(cube[loca], iter)
    end
end


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


function fbsearch(txi,a,b)
    # yub = [idx[j]*steps[j]+fmin[j] for j=1:3]
    # ylb = [(idx-[1,1,1])[j]*steps[j]+fmin[j] for j=1:3]
    idx0 = findall(x->x==0, txi)
    idx1 = findall(x->x==1, txi)
    lp = Model(CPLEX.Optimizer)
    MOI.set(lp, MOI.RawParameter("CPX_PARAM_SCRIND"), false)
    MOI.set(lp, MOI.RawParameter("CPX_PARAM_THREADS"),1  )
    @variables(lp, 0 <= x[1:n])
    @objective( lp, Min, sum(x[a] for a in idx0) + sum(1-(x[b]) for b in idx1) )
    @constraint(lp, con[i=1:m], sum(a[i,j]*x[j] for j=1:n) <= b[i])

    optimize!(lp)
    if termination_status(lp) == MOI.OPTIMAL
        return JuMP.value.(x)
    else
        return 0; #print("no new lp solution found");
    end
end
1
################################       MIP MODEL   ###############################
function mindist(sameval,groups,L2,l,que,a,P,b)
    dis = Model(CPLEX.Optimizer)
    MOI.set(dis, MOI.RawParameter("CPX_PARAM_SCRIND"), false)
    MOI.set(dis, MOI.RawParameter("CPX_PARAM_THREADS"),1  )
    @variable(dis, 0<=x[1:n])
    @constraint(dis, discon[k in sameval], x[k] == que[k] )
    cubemean = round.([mean(L2[groups[l],:][:,k]) for k=1:3])

    @constraint(dis, con[i=1:m], sum(a[i,j]*x[j] for j=1:n) <= b[i])
    @constraint(dis, objcon[k=1:3], dot(P[k,:],x) >= cubemean[k])

    @objective(dis, Max, 1);
    optimize!(dis)
    if JuMP.termination_status(dis) ==MOI.OPTIMAL
        newlpsol = JuMP.value.(x)
        return newlpsol
    else
        return nothing
    end #solveLP
end

function fbcheck(txi)
    result = true
    if all(x-> 0==x, txi)
        result = false
    else
        #1st constraint checking
        for k=1:m
            if dot(txi,a[k,:]) <= b[k]
                continue
            else
                result = false
                return result
            end
        end
    end
    return result
end

function giveobjval(x)
    return [dot(P[1,:],x),dot(P[2,:],x),dot(P[3,:],x)]
end

function getlocation(x,fmin,steps)
    L = giveobjval(x)
    loca = [round.(Int,((L[k]-fmin[k])/steps[k])+1) for k=1:3]
    return loca
end



function dominated(x,P1)
    st = false
    for k=1:length(P1)
        if all( x .>= P1[k]) && any(x > P1[k])
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

function Grouping(groups,groupkeys,lpX,steps,fmin,L2) #find common elements(variable values) from solutions in the same cuboid
    LPcount1 = 0;
    Xf= Dict(); candX = Dict(); loca_check = [];
    for l=1:length(groups)
        glength = length(groups[l])
        if glength==1 #if there is one point in a group
            if (all.(x->trunc(x)==x, lpX[groups[l]])[1] == 1) #all varval are integer => insert to Xf
                Xf[groupkeys[l]] = lpX[groups[l][1]]
                print("lp==ip \n")
            else
                candX[groupkeys[l]] = lpX[groups[l]][1] #candX is a set of solutions that should be processed to be an integer solution
                print("only lp sol in cuboid \n")
            end

        else
            cogroup = lpX[groups[l]]
            que = round.(hcat(cogroup...))
            sameval = findall(x->x==true,[all(x->x==que[l][1], que[l,:]) for l=1:length(lpX)])
            newlpsol = mindist(sameval,groups,L2,l,que,a,P,b)
            LPcount1+=1

            if newlpsol == nothing
                print("there is no LPsolution \n")
            else #when newlpsol is integer, check if all values are integer
                if (all.(x->trunc(x)==x, newlpsol) == 1) #all varval is integer => insert to Xf
                    Xf[groupkeys[l]] = newlpsol
                    # print("newlpsol==ip \n")
                else
                    push!(loca_check, newlpsol);
                    print("newlpsol in cuboid is found\n")
                    candX[groupkeys[l]] = newlpsol;
                end
            end

        end
    end
    return Xf,candX,loca_check,LPcount1
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
    Tabu = []; elaps=0; TimeLimit = 3600;
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
                    x_t = fbsearch(xi_t,a,b)
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

FPTime = @CPUelapsed Xf2,PFset2,Tabu,FPLPcount = GroupFP(Xf1,PFset,candX,LPcount,steps,fmin,fmax)
PP = convert.(Array{Int,1},collect(values(Xf2)))
Pcoordi = convert.(Array{Int,1},collect(keys(Xf2)))
Pobj = convert.(Array{Int,1}, giveobjval.(PP))

push!(PP,zeros(Int,m));push!(Pobj,[0,0,0])
#Filter dominated solutions
FPsol, FPPF = domFilter(PP,Pobj)
Pz,Py,Px=[hcat(FPPF...)[i,:] for i=1:3]; dfP = Px,Py,Pz
matP = zeros(length(FPsol),n)
for i=1:length(FPsol)
    matP[i,:] = FPsol[i]
end

genew = 0;
for i=1:length(FPsol)
    if FPsol[1] ∉ lpX
        global genew+=1
    end
end

#########################  Record outputs  ############################
ins = ARGS[1][end-2:end-0] #CPUtime recorded, naming after the instance
# colname = ARGS[1][end-12:end-4]

record1 = DataFrame(InitialInt = length(Xf1)+1, CandiX = length(candX),
    commonality = length(fixednum),avgfixedvar = mean(fixednum),FoundIPs=length(P)-(length(Xf1)+1) ,
    solvedLP=FPLPcount, totalsol=length(FPPF),CPUtime=GroupingTime+FPTime, newsol=genew)
GFPX=DataFrame(matP); GFPY=DataFrame(dfP);

CSV.write(ARGS[1]*"/GFP_"*"$ins"*"_ILPrecord.csv",record1, append=true, writeheader=false )#, delim=',' )
CSV.write(ARGS[1]*"_GFP_ILP_X.csv",GFPX, append=true, writeheader=false)
CSV.write(ARGS[1]*"_GFP_ILP_Y_.csv",GFPY, writeheader=false, delim=' ' )
print(colname," GroupFP Done!")
