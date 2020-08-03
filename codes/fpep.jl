using DelimitedFiles,DataFrames,JuMP,CPLEX,LinearAlgebra,CPUTime,CSV #,PlotlyJS,RDatasets,Colors
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
    lp = Model(CPLEX.Optimizer)
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
    # cubemax = [maximum(L2[groups[l],:][:,i]) for i=1:3]
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
function getConstraints(i,e,P)
    ϵ = [0,0];
    ϵ′ = [Inf,Inf];
    for j=2:3
        d = convert(Int,i%(length(P)+1))
        i = convert(Int,(i-d)/(length(P)+1))
        ϵ[j-1] = e[j][d+1]
        ϵ′[j-1] = e[j][d+2]
    end
    return ϵ,ϵ′
end

function updateConstraints(f,e)
    for j=2:3
        i = 1
        for k=1:length(e[j])
            if e[j][k]<f[j]
                i+=1
            else
                break
            end
        end
        insert!(e[j],i,f[j])
    end
    return e
end

function giveobjval(x)
    return [-dot(x[i+1:j+i],demand[:]),dot(x[1:i],fixcost[:]),dot(x[j+i+1:end],cost[:])]
end

function dominated(x,P)
    st = false
    if x==nothing
        return true
    else
        for k=1:length(P)
            if all( x .>= P[k])
                st=true; break
            else
                continue
            end
        end
        return st
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

###################FLP model ###############
lex1 = Model( with_optimizer(CPLEX.Optimizer )) #with_optimizer(CPLEX.Optimizer)
MOI.set(lex1, MOI.RawParameter("CPX_PARAM_SCRIND"), false )
MOI.set(lex1, MOI.RawParameter("CPX_PARAM_THREADS"),1  )
@variable(lex1, x1[1:i,1:j], Bin)
@variable(lex1, y1[1:i], Bin)
@variable(lex1, z1[1:j], Bin)
@variable(lex1, ep[1:2], Int)
@variable(lex1, ep′[1:2], Int)
@constraint(lex1, con1[b=1:j], sum(x1[a,b] for a in 1:i) == z1[b])
@constraint(lex1, con2[a=1:i,b=1:j], x1[a,b]<=y1[a])
@constraint(lex1, epsilon11,  sum(fixcost[a]*y1[a] for a=1:i) <= ep′[1]-1)
@constraint(lex1, epsilon12,  ep[1] <=sum(fixcost[a]*y1[a] for a=1:i) )
@constraint(lex1, epsilon21, sum(cost2[a][b]*x1[a,b] for a=1:i for b=1:j) <= ep′[2]-1)
@constraint(lex1, epsilon22, ep[2] <= sum(cost2[a][b]*x1[a,b] for a=1:i for b=1:j) )
@objective(lex1, Min, -sum(demand[b]*z1[b] for b=1:j))
optimize!(lex1)

lex2 = Model( CPLEX.Optimizer )
MOI.set(lex2, MOI.RawParameter("CPX_PARAM_SCRIND"), false)
MOI.set(lex2, MOI.RawParameter("CPX_PARAM_THREADS"),1  )
@variable(lex2, x2[1:i,1:j], Bin)
@variable(lex2, y2[1:i], Bin)
@variable(lex2, z2[1:j], Bin)
@variable(lex2, ep2[1:2], Int)
@variable(lex2, ep2′[1:2], Int)
@variable(lex2, obj1)
@constraint(lex2, con1[b=1:j], sum(x2[a,b] for a in 1:i) == z2[b])
@constraint(lex2, con2[a=1:i,b=1:j], x2[a,b]<=y2[a])
@constraint(lex2, con3, -sum(z2[b]*demand[b] for b=1:j) <= obj1)
@constraint(lex2, epsilon11,   sum(fixcost[a]*y2[a] for a=1:i) <= ep2′[1]-1 )
@constraint(lex2, epsilon12,  ep2[1] <= sum(fixcost[a]*y2[a] for a=1:i) )
@constraint(lex2, epsilon21, sum(cost2[a][b]*x2[a,b] for a=1:i for b=1:j) <= ep2′[2]-1 )
@constraint(lex2, epsilon22,  ep2[2] <= sum(cost2[a][b]*x2[a,b] for a=1:i for b=1:j) )
@objective(lex2, Min, sum(fixcost[a]*y2[a] for a=1:i) )
optimize!(lex2)

lex3 = Model( CPLEX.Optimizer )
MOI.set(lex3, MOI.RawParameter("CPX_PARAM_SCRIND"), false)
MOI.set(lex3, MOI.RawParameter("CPX_PARAM_THREADS"),1  )
@variable(lex3, x3[1:i,1:j], Bin)
@variable(lex3, y3[1:i], Bin)
@variable(lex3, z3[1:j], Bin)
@variable(lex3, ep3[1:2], Int)
@variable(lex3, ep3′[1:2], Int)
@variable(lex3, obj2)
@variable(lex3, obj3)
@constraint(lex3, con1[b=1:j], sum(x3[a,b] for a in 1:i) == z3[b])
@constraint(lex3, con2[a=1:i,b=1:j], x3[a,b]<=y3[a])
@constraint(lex3, con3, -sum(z3[b]*demand[b] for b=1:j) <= obj2 )
@constraint(lex3, con4 , sum(fixcost[a]*y3[a] for a=1:i) <= obj3 )
@constraint(lex3, epsilon11,  sum(fixcost[a]*y3[a] for a=1:i) <= ep3′[1]-1 )
@constraint(lex3, epsilon12,  ep3[1] <= sum(fixcost[a]*y3[a] for a=1:i))
@constraint(lex3, epsilon21,  sum(cost2[a][b]*x3[a,b] for a=1:i for b=1:j) <= ep3′[2]-1 )
@constraint(lex3, epsilon22,  ep3[2] <= sum(cost2[a][b]*x3[a,b] for a=1:i for b=1:j)  )
@objective(lex3, Min, sum(cost2[a][b]*x3[a,b] for a=1:i for b=1:j))
optimize!(lex3)

function opt(ϵ,ϵ′,solvedIP)
    JuMP.fix(ep[1], ϵ[1]); JuMP.fix(ep[2], ϵ[2])
    JuMP.fix(ep′[1],ϵ′[1]); JuMP.fix(ep′[2], ϵ′[2])
    optimize!(lex1)

    if termination_status(lex1) == MOI.OPTIMAL
        solvedIP+=1
        JuMP.fix(ep2[1], ϵ[1]); JuMP.fix(ep2[2], ϵ[2])
        JuMP.fix(ep2′[1], ϵ′[1]); JuMP.fix(ep2′[2], ϵ′[2])
        JuMP.fix(obj1, objective_value(lex1))
        optimize!(lex2)

        if (termination_status(lex2) == MOI.OPTIMAL)
            solvedIP+=1
            JuMP.fix(ep3[1], ϵ[1]); JuMP.fix(ep3[2], ϵ[2])
            JuMP.fix(ep3′[1],ϵ′[1]); JuMP.fix(ep3′[2], ϵ′[2])
            JuMP.fix(obj2, objective_value(lex1)); JuMP.fix(obj3, objective_value(lex2))
            optimize!(lex3)
            if (termination_status(lex3) == MOI.OPTIMAL)
                solvedIP+=1
                lexY = round.(Int,value.(y3)); lexZ = round.(Int,value.(z3)); lexX = round.(Int,value.(x3));
                s =  append!( append!(lexY,lexZ),transpose(lexX) )
                return s,giveobjval(s),solvedIP
            else
                # print("Lex3: ",termination_status(lex3),"\n")
                return nothing,nothing,solvedIP
            end
        else
            # print("Lex2: ",termination_status(lex2),"\n")
            return nothing,nothing,solvedIP
        end
    else
        # print("Lex1: ",termination_status(lex1),"\n")
        return nothing,nothing,solvedIP
    end
end
##
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
            xi_t = round.(x_t)
            fxi = giveobjval(xi_t)
            if ( (fbcheck(xi_t) == true) && (dominated(fxi,PFset)==false) ) #checking easibility and dominance   #z[1]<=fmax[1]-1 && z[2]<=fmax[2]-1 && z[3]<=fmax[3]-1)
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
            if CPUtime_us()-algo_start < 3600*1e6
                continue
            else
                return Xf,PFset,Tabu,LPcount
            end

        end
    end
    # end
    return Xf,PFset,Tabu,LPcount
end
#Time Feasibility Pump
FPTime = @CPUelapsed Xf2,PFset2,Tabu,FPLPcount = GroupFP(Xf1,PFset,candX,LPcount,steps,fmin,fmax)
P = convert.(Array{Int,1},collect(values(Xf2)))
Pcoordi = convert.(Array{Int,1},collect(keys(Xf2)))
Pobj = convert.(Array{Int,1}, giveobjval.(P))
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
# detectlex = [count(i->i==0,FPPF[k]) for k=1:length(FPPF)]
lexid = findall(x->x>=2,detectlex)
for i=1:length(lexid)
    deleteat!(FPsol,lexid[i]); deleteat!(FPPF,lexid[i]); deleteat!(FPcoordi,lexid[i])
end

Pz,Py,Px=[hcat(FPPF...)[i,:] for i=1:3]; dfP = Px,Py,Pz
matP = zeros(length(FPsol),i+j+i*j)
for i=1:length(FPsol)
    matP[i,:] = FPsol[i]
end

Solpair = Dict(); PFpair = Dict()
for i=1:length(FPPF)
    Solpair[FPcoordi[i]] = FPsol[i];
    PFpair[FPcoordi[i]] = FPPF[i]
end
keyunion = union(FPcoordi,groupkeys)
searchdone = Dict()
for i = 1:length(keyunion)
    searchdone[keyunion[i]] = false
end
####################   Loop Feasibility Pump + ϵ-constraint    ######################
function FPep2(S,T,keyunion,allsol)
    nosol = copy(T); global solvedIP=0; algo_start = CPUtime_us();iter=1;F=[]
    while all(collect(values(searchdone)))!=true #true #iter<5*length(T)
        for k in keyunion
            print( "=============FP+ϵ ",iter,"th==+ ",k," cuboid===========\n")
            # Setting the stopping condition. either % or a ratio of total solution exceeds
            if length(T[k])>(allsol*0.3) || nosol[k]==true
                searchdone[k] = true
                @goto next_group
            else
                @show nosol[k]
                temp1 = (fmin+steps.*k)[2:3]
                temp2 = (fmin+steps.*(k+[1,1,1]))[2:3]
                global e = Dict( 2=>[temp1[1],temp2[1]], 3=>[temp1[2],temp2[2]] )
                fval = [ T[k][i] for i=1:length(T[k]) ]
                # Put all given pareto solutions in "Array e"
                for i=1:length(fval)
                    e =updateConstraints(fval[i],e)
                end
                m = 0; c = (length(T[k])+1)^2;
                print("\n======== m,c = ", m," ",c,"  ========\n" )
                # Proceed iteration until explore whole divided search region
                newsol=0
                while m<c
                    global (ϵ,ϵ′) = getConstraints(m,e,T[k])
                    # print("new epsilons= ",ϵ,ϵ′,"\n")
                    if [ϵ,ϵ′] ∉ F
                        global (s,objs,solvedIP) = opt(ϵ,ϵ′,solvedIP) #fval
                        solvedIP+=3
                        if ( s == nothing || dominated(objs, [T[k][j] for j=1:length(T[k])] )==true)
                            F = union(F,[[ϵ,ϵ′]]);

                            # print("USELESS new solution \n")
                        else
                            push!(S[k], s)
                            push!(T[k], objs)
                            F = union(F, [ [objs[2],objs[3]], [ϵ′[1],ϵ′[2]] ] )
                            e = updateConstraints(objs,e)
                            newsol=+1
                            print("found solution: ",objs , "\n")
                        end
                    end
                    m+=1;
                    # print("Next region \n")
                    if pretime+(CPUtime_us()-algo_start) < TimeLimit*1e6
                        continue
                    else
                        return S,T,solvedIP
                    end
                end
                if c==0 || newsol==0
                    # print("no solution in a cuboid\n")
                    nosol[k] = true
                end
            end

            @label next_group
            allsol = sum([length(T[k]) for k in keyunion])
        end
        iter+=1
    end
    return S,T,solvedIP
end

S = Dict();T = Dict()
for i in keyunion
    if i in FPcoordi
        S[i] = [ Solpair[i] ];T[i] = [ PFpair[i] ]
    else
        S[i] = [];T[i] = []
    end
end
# print("#cubes= ", length(FPPF));
pretime = GroupingTime+FPTime
eFP2Time = @CPUelapsed S2,T2,eFP2IPcount = FPep2(S,T,keyunion,length(FPPF))
T3=collect(values(T2)); T4 = []; S3 = collect(values(S2));S4 = [];
for k=1:length(T3)
    for l=1:length(T3[k])
        push!(T4,T3[k][l])
        push!(S4,S3[k][l])
    end
end
ndS,ndT = domFilter(S4,T4)
Tz,Ty,Tx=[hcat(ndT...)[i,:] for i=1:3]; dfT=Tx,Ty,Tz
matT = zeros(length(ndS),i+j+i*j)
for i=1:length(ndS)
    matT[i,:] = ndS[i]
end
ins = ARGS[1][end-12:end-7] #CPUtime recorded, naming after the instance
colname = ARGS[1][end-12:end-4]
#record1 = DataFrame( Initial_IPsol = length(Xf1), Initial_LPsol = length(candX),LP=FPLPcount, sol=length(Pz), removed=length(P)-length(Pz))
#insertcols!(record1,6, Symbol("$colname")=>GroupingTime+FPTime)
#GFPX=DataFrame(matP); GFPY=DataFrame(dfP);
#CSV.write("/home/k2g00/k2g3475/multiobjective/solvers/ep+FP/GFPresults/GFPlex_"*"$ins"*"_record.csv",record1, append=true, writeheader=false )#, delim=',' )
#CSV.write(ARGS[1]*"GFP_X.csv",GFPX, append=false, writeheader=false)
#CSV.write(ARGS[1]*"GFPlex_Y_.csv",GFPY, writeheader=false, delim=' ' )
#print(colname," GroupFP Done!")

record2 = DataFrame(solvedIP=eFP2IPcount,sols=length(ndT),removed=length(T4)-length(ndT))
insertcols!(record2,4,Symbol("$colname")=>eFP2Time)
eFP2X=DataFrame(matT);eFP2Y=DataFrame(dfT);
CSV.write("/home/k2g00/k2g3475/multiobjective/solvers/ep+FP/FPepresults/FPep_2hr_"*"$ins"*"_record.csv",record2, append=true, writeheader=false )#, delim=',' )
#CSV.write(ARGS[1]*"FPep_X_2hr.csv",eFP2X, append=false, writeheader=false)
CSV.write(ARGS[1]*"FPep_2hr_Y_.csv",eFP2Y, writeheader=false, delim=' ' )
print(colname," FP+ep Done!")
