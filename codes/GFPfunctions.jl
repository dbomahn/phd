using DelimitedFiles,DataFrames,JuMP,CPLEX,LinearAlgebra,CPUTime,CSV,Statistics #,PlotlyJS,RDatasets,Colors
global TimeLimit = 300
################################  Data  ####################################
struct path{S<:String}
    dir1::S; dir2::S; dir3::S
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
function KPfbcheck(txi)
    result = true
    if all(x-> 0==x, txi)
        result = false
    else
        #constraint checking
        if -sum(weight[i]*txi[i] for i=1:n) >= -ub
        else
            result = false
            return result
        end
    end
    return result
end
function fbsearch(idx,steps,fmin,txi,fmax) #solveLP
    # yub = [idx[j]*steps[j]+fmin[j] for j=1:3]
    # ylb = [(idx-[1,1,1])[j]*steps[j]+fmin[j] for j=1:3]
    idx0 = findall(x->x==0, txi)
    idx1 = findall(x->x==1, txi)
    kp_m = Model(CPLEX.Optimizer)
    MOI.set(kp_m, MOI.RawParameter("CPX_PARAM_SCRIND"), false)
    MOI.set(kp_m, MOI.RawParameter("CPX_PARAM_THREADS"),1  )
    @variables(kp_m, begin
        0 <= x[1:n] <=1
    end)
    @objective( kp_m, Min, sum(x[n] for n in idx0) + sum(1-(x[m]) for m in idx1) )
    @constraint( kp_m, -sum(weight[i]*x[i] for i=1:n) >= -ub )
    # @constraint( kp_m, c1, -dot(C[1,:],x[:]) <= yub[1]-1 )#ylb[1]+1 <=
    # @constraint( kp_m, c2, -dot(C[2,:],x[:]) <= yub[2]-1 ) #ylb[2]+1 <=
    # @constraint( kp_m, c3, -dot(C[3,:],x[:]) <= yub[3]-1 ) #ylb[3]+1 <=
    optimize!(kp_m)
    if termination_status(kp_m) == MOI.OPTIMAL
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
    @variable(dis, 0<=x[1:n]<=1)
    @variable(dis, y==1)
    @constraint(dis, con[k in sameval], x[k] == que[k] )
    cubemean = [mean(L2[groups[l],:][:,i]) for i=1:3]

    @constraint(dis,-dot(C[1,:],x[:])  == -cubemean[1])
    @constraint(dis, -dot(C[2,:],x[:]) == -cubemean[2])
    @constraint(dis, -dot(C[3,:],x[:]) == -cubemean[3])
    @constraint( dis, -sum(weight[i]*x[i] for i=1:n) >= -ub )
    @objective(dis, Min, y)
    optimize!(dis)
    if JuMP.termination_status(dis) ==MOI.OPTIMAL
        newlpsol = JuMP.value.(x)
        return newlpsol
    else
        return nothing
    end #solveLP
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
function Grouping(groups,groupkeys,lpX,steps,fmin,L2,n) #find common elements(variable values) from solutions in the same cuboid
    LPcount1 = 0;
    Xf= Dict(); candX = Dict(); loca_check = [];
    for l=1:length(groups)
        glength = length(groups[l])
        if glength==1 #if there is one point in a group
            if (all.(x->trunc(x)==x, lpX[groups[l]])[1] == 1) #all varival are integer => insert to Xf
                Xf[groupkeys[l]] = lpX[groups[l][1]]
                # print("lp==ip \n")
            else
                candX[groupkeys[l]] = lpX[groups[l]][1] #candX is a set of solutions that should be processed to be an integer solution
                # print("only lp sol in cuboid \n")
            end

        else
            cogroup = lpX[groups[l]]
            que = round.(hcat(cogroup...))
            sameval = findall(x->x==true,[all(x->x==que[l][1], que[l,:]) for l=1:n])
            newlpsol = mindist(sameval,groups,groupkeys,L2,l,que)
            LPcount1+=1

            if newlpsol == nothing
                # print("there is no LPsolution \n")
            else #when newlpsol is integer, check if all values are integer
                if (all.(x->trunc(x)==x, newlpsol) == 1) #all varival is integer => insert to Xf
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
## for tri- epsilon
function giveobjval(x)
    return [-dot(C[1,:],x),-dot(C[2,:],x),-dot(C[3,:],x)]
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
