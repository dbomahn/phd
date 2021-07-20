using DelimitedFiles,DataFrames,JuMP,CPLEX,LinearAlgebra,CPUTime,CSV
TimeLimit=Inf
############### Importing Data ############
mutable struct Data
    input::String
    i::Int; j::Int;
    fixcost::Array{}
    demand::Array{}
    cost::Array{}
    cost2::Array{}

    function Data(input::String)
        data=readdlm(input)
        i,j = filter!(a->typeof(a)==Int, data[1,:]) # facility_i, customer_j # for FLP
        fixcost,capa,cost,cost2 = [],[],[],[]
        for k=2:i+1
            fix,cp = filter!(a->typeof(a)==Int, data[k,:])
            push!(fixcost,fix)
        end
        demand = filter!(a->typeof(a)==Int, data[2+i,:])
        a,b = size(data)
        for k=i+3:a
            ct = filter!(a->typeof(a)==Int, data[k,:])
            push!(cost2,ct)
            for l=1:length(ct)
                push!(cost,ct[l])
            end
        end
        new(input,i,j,fixcost,demand,cost,cost2)
    end
end
data=Data("/home/ak121396/multiobjective/instances/triFLP/instances/05_010_03.txt") #05_010_09.txt

struct path{S<:String}
    dir1::S; dir2::S;
end
# input
pathes=path("/home/ak121396/multiobjective/instances/triFLP/", "/home/ak121396/multiobjective/instances/triFLP/PF/")
dir1 = pathes.dir1 # "/home/ak121396/Desktop/FLPInstances/"
dir2 = pathes.dir2 #"/home/ak121396/Desktop/FLPInstances/FLPvlp/PF/"
subdir = readdir(dir1*"/instances/")[3]
f = readdlm(dir1*"/instances/"*subdir) #f = readdlm(dir1*"/acc-tight4_pre_img_p.sol")
presol = readdir(dir1*"/varval/")[3] #FLPvlp
img_p = round.(readdlm(dir1*"/varval/"*presol), digits=2) #round numbers till 2 digits
pf = round.(Int,readdlm(dir2*presol[1:9]*"img_p.sol")[:,2:end])
# facility_i, customer_j # for FLP
const i,j = filter!(a->typeof(a)==Int, f[1,:])
(const fixcost,capa,cost,cost2 = [],[],[],[] )
for k=2:i+1
    fix,cp = filter!(a->typeof(a)==Int, f[k,:])
    push!(fixcost,fix)
    push!(capa,cp)
end
const demand = filter!(a->typeof(a)==Int, f[2+i,:])
a,b = size(f)
for k=i+3:a
    ct = filter!(a->typeof(a)==Int, f[k,:])
    push!(cost2,ct)
    for l=1:length(ct)
        push!(cost,ct[l])
    end
end
s,ss = size(img_p)
indx= findall(x->x!=0,[sum(img_p[i,:]) for i=1:s]) #NZ sol
lpX = [img_p[i,:] for i in indx]
pf = hcat(pf[:,3],pf[:,1],pf[:,2])
L = reshape([pf[i,:] for i in indx],length(lpX),1)
L2 = transpose(hcat(L...))
fmin = Array{Int}([minimum(L2[:,i]) for i=1:3])
fmax = Array{Int}([maximum(L2[:,i]) for i=1:3])
steps = [round.(Int,abs(fmax[i]-fmin[i])/length(L)) for i=1:3] #determine steps according to #customers(j)
cube = Dict();
for iter=1:length(L)
    loca = [round.(Int,((L2[iter,k]-fmin[k])/steps[k])+1) for k=1:3]
    if !haskey(cube,loca)
        cube[loca] = [iter]
    else
        push!(cube[loca], iter)
    end
end
1
###################FLP model ###############
# mutable struct Problem
#     model::JuMP.Model
#     x::Array{VariableRef,2}
#     y::Array{VariableRef,1}
#     z::Array{VariableRef,1}
#     ep::Array{VariableRef,1}
#     ep′::Array{VariableRef,1}
#     obj1::VariableRef
#     obj2::VariableRef
# end
#
# function FLPlex(FLP::Data)::Problem
#     flp = JuMP.Model(CPLEX.Optimizer)
#     MOI.set(flp, MOI.RawParameter("CPX_PARAM_SCRIND"), false)
#     @variable(flp, x[1:FLP.i,1:FLP.j], Bin)
#     @variable(flp, y[1:FLP.i], Bin)
#     @variable(flp, z[1:FLP.j], Bin)
#     @variable(flp, ep[1:2], Int)
#     @variable(flp, ep′[1:2], Int)
#     @variable(flp, obj1)
#     @variable(flp, obj2)
#     @constraint(flp, con1[b=1:FLP.j], sum(x[a,b] for a in 1:FLP.i) == z[b])
#     @constraint(flp, con2[a=1:FLP.i,b=1:FLP.j], x[a,b]<=y[a])
#     @constraint(flp, con3, -sum(z[b]*FLP.demand[b] for b=1:FLP.j) <= obj1 )
#     @constraint(flp, con4 , sum(FLP.fixcost[a]*y[a] for a=1:FLP.i) <= obj2 )
#     @constraint(flp, epsilon11,  sum(FLP.fixcost[a]*y[a] for a=1:FLP.i) <= ep′[1]-1)
#     @constraint(flp, epsilon12,  ep[1] <=sum(FLP.fixcost[a]*y[a] for a=1:FLP.i) )
#     @constraint(flp, epsilon21, sum(FLP.cost2[a][b]*x[a,b] for a=1:FLP.i for b=1:FLP.j) <= ep′[2]-1)
#     @constraint(flp, epsilon22, ep[2] <= sum(FLP.cost2[a][b]*x[a,b] for a=1:FLP.i for b=1:FLP.j) )
#     Problem(flp,x,y,z,ep,ep′,obj1,obj2)
# end
# function lex1(lex::Problem,FLP::Data)
#     @objective(lex.model,Min,-sum(FLP.demand[b]*lex.z[b] for b=1:FLP.j))
#     Problem(lex.model,lex.x,lex.y,lex.z,lex.ep,lex.ep′,lex.obj1,lex.obj2)
# end
# function lex2(lex::Problem,FLP::Data)
#     @objective(lex.model,Min,sum(FLP.fixcost[a]*lex.y[a] for a=1:FLP.i))
#     Problem(lex.model,lex.x,lex.y,lex.z,lex.ep,lex.ep′,lex.obj1,lex.obj2)
# end
# function lex3(lex::Problem,FLP::Data)
#     @objective(lex.model,Min, sum(FLP.cost2[a][b]*lex.x[a,b] for a=1:FLP.i for b=1:FLP.j))
#     Problem(lex.model,lex.x,lex.y,lex.z,lex.ep,lex.ep′,lex.obj1,lex.obj2)
# end
# lexi1 = lex1(FLPlex(data),data);lexi2 = lex2(FLPlex(data),data);lexi3 = lex3(FLPlex(data),data)
# optimize!(lexi1.model) #test run
#
# function opt(lex1::Problem,lex2::Problem,lex3::Problem,FLP::Data,ϵ,ϵ′,solvedIP)
#     JuMP.fix(lex1.ep[1], ϵ[1]); JuMP.fix(lex1.ep[2], ϵ[2])
#     JuMP.fix(lex1.ep′[1],ϵ′[1]); JuMP.fix(lex1.ep′[2], ϵ′[2])
#     optimize!(lex1.model)
#
#     if termination_status(lex1.model) == MOI.OPTIMAL
#         solvedIP+=1; val1 = objective_value(lex1.model)
#         JuMP.fix(lex2.ep[1], ϵ[1])
#         JuMP.fix(lex2.ep[2], ϵ[2])
#         JuMP.fix(lex2.ep′[1],ϵ′[1])
#         JuMP.fix(lex2.ep′[2], ϵ′[2])
#         JuMP.fix(lex2.obj1, val1);
#
#         optimize!(lex2.model)
#
#         if (termination_status(lex2.model) == MOI.OPTIMAL)
#             solvedIP+=1; val2 = objective_value(lex2.model)
#             JuMP.fix(lex3.ep[1], ϵ[1]); JuMP.fix(lex3.ep[2], ϵ[2])
#             JuMP.fix(lex3.ep′[1],ϵ′[1]); JuMP.fix(lex3.ep′[2], ϵ′[2])
#             JuMP.fix(lex3.obj1, val1); JuMP.fix(lex3.obj2, val2)
#             optimize!(lex3.model)
#             if (termination_status(lex3.model) == MOI.OPTIMAL)
#                 solvedIP+=1
#                 lexY = round.(Int,value.(lex3.y)); lexZ = round.(Int,value.(lex3.z)); lexX = round.(Int,value.(lex3.x));
#                 s =  append!( append!(lexY,lexZ),transpose(lexX) )
#                 return s,giveobjval(s,FLP),solvedIP
#             else
#                 print("Lex3: ",termination_status(lex3.model),"\n")
#                 return nothing,nothing,solvedIP
#             end
#         else
#             print("Lex2: ",termination_status(lex2.model),"\n")
#             return nothing,nothing,solvedIP
#         end
#     else
#         print("Lex1: ",termination_status(lex1.model),"\n")
#         return nothing,nothing,solvedIP
#     end
# end
########################        for Feasibility Pump      ##########################

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
        return result
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

function fbsearch(idx,steps,fmin,txi,fmax,FLP::Data) #solveLP
    demand = FLP.demand; fixcost = FLP.fixcost; cost = FLP.cost
    yub = [idx[j]*steps[j]+fmin[j] for j=1:3]
    ylb = [(idx-[1,1,1])[j]*steps[j]+fmin[j] for j=1:3]
    idx0 = findall(x->x==0, txi)
    idx1 = findall(x->x==1, txi)
    lp = Model(CPLEX.Optimizer)
    MOI.set(lp, MOI.RawParameter("CPX_PARAM_SCRIND"), false)
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

    dis = Model(CPLEX.Optimizer)
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
            if (all.(x->x∈[0,1], lpX[groups[l]])[1] == 1) #all varval are integer => insert to Xf
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
            if length(sameval)>0
                newlpsol = mindist(sameval,groups,groupkeys,L2,l,que,data)
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
            else #If there is no common feature among solutions. choose one solution at random
                @show candX[groupkeys[l]]  = lpX[[groups]][rand(1:glength)]
                print("pass solving LP \n")
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




######################  Loop Feasibility Pump + ϵ-constraint    ######################
function FPep2(S,T,keyunion,allsol)
    nosol = copy(T); global solvedIP=0; algo_start = CPUtime_us();iter=1;
    while true #iter<5*length(T)
        for k in keyunion
            print( "=============FP+ϵ ",iter,"th==+ ",k," cuboid===========\n")
            # Setting the stopping condition. either % or a ratio of total solution exceeds
            if length(T[k])>(allsol*0.3) || nosol[k]==true
                @goto next_group
            else
                temp1 = (fmin+steps.*k)[2:3]
                temp2 = (fmin+steps.*(k+[1,1,1]))[2:3]
                global e = Dict( 2=>[temp1[1],temp2[1]], 3=>[temp1[2],temp2[2]] )
                fval = [ T[k][i] for i=1:length(T[k]) ]
                # Put all given pareto solutions in "Array e"
                for i=1:length(fval)
                    e =updateConstraints(fval[i],e)
                end
                m = 0; c = (length(T[k])+1)^2; F=[]
                # print("\n======== m,c = ", m," ",c,"  ========\n" )
                # Proceed iteration until explore whole divided search region
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
                            # print("found solution: ",objs , "\n")
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
            end
            if c==0
                # print("no solution in a cuboid\n")
                nosol[k] = true
            end
            @label next_group
            allsol = sum([length(T[k]) for k in keyunion])
        end

        #iter+=1
    end
    return S,T,solvedIP
end

########################        for tri_epsilon        ##########################
# function getConstraints(i,e,P)
#     ϵ = [0,0];
#     ϵ′ = [Inf,Inf];
#     for j=2:3
#         d = round(Int,i%(length(P)+1))
#         i = round(Int,(i-d)/(length(P)+1))
#         ϵ[j-1] = e[j][d+1]
#         ϵ′[j-1] = e[j][d+2]
#     end
#     return ϵ,ϵ′
# end
#
# function updateConstraints(f,e)
#     for j=2:3
#         i = 1
#         for k=1:length(e[j])
#             if e[j][k]<f[j]
#                 i+=1
#             else
#                 break
#             end
#         end
#         insert!(e[j],i,f[j])
#     end
#     return e
# end
#
# function giveobjval(x,FLP::Data)
#     if x!=nothing
#         return [-dot(x[FLP.i+1:FLP.j+FLP.i],FLP.demand[:]),dot(x[1:FLP.i],FLP.fixcost[:]),dot(x[FLP.j+FLP.i+1:end],(FLP.cost[:]))]
#     else
#         return nothing
#     end
# end
#
# function dominated(x,P)
#     st = false
#     if x==nothing
#         return true
#     else
#         for k=1:length(P)
#             if all( x .>= P[k])&& any(x > P[k])
#                 st=true; break
#             else
#                 continue
#             end
#         end
#         return st
#     end
# end
#
function domFilter(sol,obj)
    copysol = Dict(); copyobj = Dict();
    for i=1:length(obj)
        copysol[i] = sol[i]
        copyobj[i] = obj[i]
    end

    for i=1:length(obj)-1
        for j=i+1:length(obj)
            if all(obj[i] .>= obj[j]) == true #dominated by PF[j]
                copyobj[i]=nothing; copysol[i]=nothing; print(obj[i],"dom by ", obj[j],"\n");break
            elseif all(obj[j] .>= obj[i]) == true
                copyobj[j]=nothing; copysol[j]=nothing; print(obj[j],"dom by ", obj[i],"\n")
            end
        end
    end

    finalsol = filter!(a->a!=nothing, collect(values(copysol)))
    finalobj = filter!(a->a!=nothing, collect(values(copyobj)))

    return finalsol,finalobj
end


#####################  Feasibility Pump + ϵ-constraint    ######################

function FPep2(S,T,keyunion,allsol)
    nosol = copy(T); global solvedIP=0; algo_start = CPUtime_us();iter=1;F=[]
    while all(collect(values(searchdone)))!=true #true #iter<5*length(T)
        for k in keyunion
            print( "=============FP+ϵ ",iter,"th==+ ",k," cuboid===========\n")
            # Setting the stopping condition. either % or a ratio of total solution exceeds
            if length(T[k])>(allsol*0.3) || nosol[k]==true
                searchdone[k] = true
                print("moving to the next")
                @goto next_group
            else
                temp1 = (fmin+steps.*k)[2:3]
                temp2 = (fmin+steps.*(k+[1,1,1]))[2:3]
                global e = Dict( 2=>[temp1[1],temp2[1]], 3=>[temp1[2],temp2[2]] )
                fval = [ T[k][i] for i=1:length(T[k]) ]
                # Put all given pareto solutions in "Array e"
                for i=1:length(fval)
                    e =updateConstraints(fval[i],e)
                end
                m = 0; c = (length(T[k])+1)^2;
                # print("\n======== m,c = ", m," ",c,"  ========\n" )
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
                            # print("found solution: ",objs , "\n")
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


print("function.jl compiled!")

# function opt2(ϵ,ϵ′)
#     JuMP.fix(ep[1], ϵ[1]); JuMP.fix(ep[2], ϵ[2])
#     JuMP.fix(ep′[1],ϵ′[1]); JuMP.fix(ep′[2], ϵ′[2])
#     optimize!(lex1)
#
#     if termination_status(lex1) == MOI.OPTIMAL
#         JuMP.fix(ep2[1], ϵ[1]); JuMP.fix(ep2[2], ϵ[2])
#         JuMP.fix(ep2′[1], ϵ′[1]); JuMP.fix(ep2′[2], ϵ′[2])
#         JuMP.fix(obj1, objective_value(lex1))
#         optimize!(lex2)
#
#         if (termination_status(lex2) == MOI.OPTIMAL)
#             JuMP.fix(ep3[1], ϵ[1]); JuMP.fix(ep3[2], ϵ[2])
#             JuMP.fix(ep3′[1],ϵ′[1]); JuMP.fix(ep3′[2], ϵ′[2])
#             JuMP.fix(obj2, objective_value(lex1)); JuMP.fix(obj3, objective_value(lex2))
#             optimize!(lex3)
#             if (termination_status(lex3) == MOI.OPTIMAL)
#                 lexY = round.(Int,value.(y3)); lexZ = round.(Int,value.(z3)); lexX = round.(Int,value.(x3));
#                 s =  append!( append!(lexY,lexZ),transpose(lexX) )
#                 return s,giveobjval(s)
#             else
#                 # print("Lex3: ",termination_status(lex3),"\n")
#                 return nothing,nothing
#             end
#         else
#             # print("Lex2: ",termination_status(lex2),"\n")
#             return nothing,nothing
#         end
#     else
#         # print("Lex1: ",termination_status(lex1),"\n")
#         return nothing,nothing
#     end
# end
# function FeasPump(Xf,candX,LPcount,steps,fmin,fmax)
#     Tabu = [];
#     # candkeys = collect(keys(candX)); candvals = collect(values(candX));
#     getlocation()
#     for idx in candkeys
#         x_t = candX[idx];
#         SearchDone = false
#         itr = 1
#         Max_itr = (i+j)*100 #max(count(x->0<x<1,x_t),1) #maximum number of attempts => How to set
#         while itr<=Max_itr && SearchDone==false
#             # if all(x->x==1,[trunc(Int,x_t[k]) == x_t[k] for k=1:i+j+i*j])
#             #     if !haskey(Xf,idx)
#             #         Xf[idx] = x_t
#             #     else
#             #         Xf = union(Xf[idx], [x_t])
#             #         # push!(Xf[i],xi_t)
#             #     end
#             #     SearchDone= true; print("LPsol == Intsol \n") #checking integer value is already done in GroupFP
#             # else
#             xi_t = round.(x_t)
#             # z = giveobjval(xi_t)
#             # loca = [round.(Int,((z[k]-fmin[k])/steps[k])+1) for k=1:3] #restrict the search region
#             if ( xi_t == fbcheck(xi_t) && dominated(xi_t,collect(values(Xf))==false) ) #checking easibility and dominance   #z[1]<=fmax[1]-1 && z[2]<=fmax[2]-1 && z[3]<=fmax[3]-1)
#                 LPcount+=1
#                 if !haskey(Xf,idx)
#                     Xf[idx] = xi_t
#                 else
#                     Xf = union(Xf[idx], [xi_t])
#                     # push!(Xf[i],xi_t)
#                 end
#                 SearchDone = true ; print("Rounding works! \n ")
#             else
#                 if xi_t ∈ Tabu
#                     xi_t = flipoper(Tabu,x_t,xi_t)
#                     if xi_t==[]
#                         SearchDone = true
#                     else
#                         if ( xi_t == fbcheck(xi_t) && dominated(xi_t,collect(values(Xf))==false) ) #z[1]<=fmax[1]-1 && z[2]<=fmax[2]-1 && z[3]<=fmax[3]-1 )
#                             LPcount+=1
#                             if !haskey(Xf,idx)
#                                 Xf[idx] = xi_t
#                             else
#                                 Xf = union(Xf[idx], [xi_t])
#                             end
#                             SearchDone = true; print("Flip operation works! \n ")
#                         end
#                     end
#                 end
#                 if SearchDone == false
#                     push!(Tabu,xi_t) #when break
#                     x_t = fbsearch(idx,steps,fmin,xi_t)
#                     LPcount+=1
#                     print("New lp found \n")
#                 end
#             end
#             itr+=1; @show itr
#         end
#     end
#     return Xf,Tabu,LPcount
# end

# function FeasiPump2019(y,lpX)
#     Tabu = [];
#     Xf= [];
#     for k = 1:length(lpX)
#         tx = lpX[k]
#         SearchDone = false
#         itr = 1
#         Max_itr = (i+j)*10 #max(count(x->0<x<1,tx),1) #maximum number of attempts
#         while itr<=Max_itr && SearchDone==false
#             if all(x->x==1,[trunc(Int,tx[x]) == tx[x] for x=1:i+j+i*j])
#                 push!(Xf,tx)
#                 SearchDone = true
#             else
#                 txi = round.(tx)
#                 z = Giveobjval(txi)
#                 # zval = dot(txi[1:i],fixcost),-dot(txi[i+1:i+j],demand),dot(txi[i+j+1:end],cost)
#                 if (txi == fbcheck(txi)  && z[1]<=y[1]-1 && z[2]<=y[2]-1 && z[3]<=y[3]-1 )
#                     push!(Xf,txi)
#                     print("Rounding works! \n ")
#                     SearchDone = true
#                 else
#                     if txi ∈ Tabu
#                         txi = flipoper(Tabu,tx,txi)
#                         if txi==[]
#                             SearchDone = true
#                         else
#                             if (txi == fbcheck(txi) &&  z[1]<=y[1]-1 && z[2]<=y[2]-1 && z[3]<=y[3]-1 )
#                                 print("Rounding works! \n ")
#                                 SearchDone = true
#                             end
#                         end
#                     end
#                     if SearchDone == false
#                         push!(Tabu,txi) #when break
#                         tx = fbsearch2(y,txi)
#                         if tx==0
#                             print("no new lp sol's found")
#                         else
#                             print("New lp found \n")
#                         end
#                     end
#                 end
#             end
#             itr+=1; #@show itr
#         end
#     end
#     return Xf,Tabu
# end
