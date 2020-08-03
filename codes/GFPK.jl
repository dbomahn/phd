using DelimitedFiles,DataFrames,JuMP,MathOptFormat,CPLEX,LinearAlgebra,CPUTime,CSV,Statistics

function writemodel(model::Model, filename::String)
     if endswith(lowercase(filename), ".lp")
         file = MathOptFormat.LP.Model()
     elseif endswith(lowercase(filename), ".mps")
         file = MathOptFormat.MPS.Model()
     else
         println("Unknown file extension to write model: ",
                 split(filename, ".")[end])
         exit(8)
     end
     MOI.copy_to(file, backend(model))
     MOI.write_to_file(file, filename)
     println("wrote model to ", filename)
end


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
mutable struct Variable
    input::String; lpX::Array{}; indx::Array{}
    function Variable(input::String,d::Data)
        i=data.i; j=data.j
        var=round.(readdlm(input), digits=2)
        s,ss = size(var)
        var2 = hcat(var[:,i+j+1:end],var[:,1:i],var[:,i+1:i+j])
        indx= findall(x->x!=0,[sum(var2[k,:]) for k=1:s]) #NZ sol
        lpX = [var2[k,:] for k in indx]
        new(input,lpX,indx)
    end
end
mutable struct Prelim
    input::String; L2::Array{}; fmin::Array{}; fmax::Array{}; steps::Array{};
    cubes::Dict{}; groups::Array{}; groupkeys::Array{};

    function Prelim(input::String,LPsol::Variable,data::Data)
        pfs = readdlm(input)[:,2:end]
        objs = round.(Int,pfs)
        pf = hcat(objs[:,3],objs[:,1],objs[:,2])
        L = reshape([pf[i,:] for i in LPsol.indx],length(LPsol.lpX),1)
        L2 = transpose(hcat(L...))
        fmin = Array{Int}([minimum(L2[:,i]) for i=1:3])
        fmax = Array{Int}([maximum(L2[:,i]) for i=1:3])
        steps = [round.(Int,abs(fmax[k]-fmin[k])/length(L)*data.i) for k=1:3] #determine steps according to #customers(j)
        cubes = Dict();
        for iter=1:length(L)
            loca = [round.(Int,((L2[iter,k]-fmin[k])/steps[k])+1) for k=1:3]
            if !haskey(cubes,loca)
                cubes[loca] = [iter]
            else
                push!(cubes[loca], iter)
            end
        end
        groups = collect(values(cubes)); groupkeys = collect(keys(cubes))

        new(input,L2,fmin,fmax,steps,cubes,groups,groupkeys)
    end
end

data=Data("/home/ak121396/multiobjective/instances/triFLP/instances/20_040_02.txt") #05_010_09.txt
vari = Variable("/home/ak121396/multiobjective/instances/triFLP/varval/20_40_02_pre_img_p.sol",data)
pre = Prelim("/home/ak121396/multiobjective/instances/triFLP/PF/20_40_02_img_p.sol",vari,data)

function fbsearch(idx,txi,steps,fmin,fmax,i,j,cost,fixcost,demand)
    # yub = [idx[k]*steps[k]+fmin[k] for k=1:3]
    # ylb = [(idx-[1,1,1])[k]*steps[k]+fmin[k] for k=1:3]
    idx0 = findall(x->x==0, txi)
    idx1 = findall(x->x==1, txi)
    lp = Model(CPLEX.Optimizer)
    MOI.set(lp, MOI.RawParameter("CPX_PARAM_SCRIND"), false)
    MOI.set(lp, MOI.RawParameter("CPX_PARAM_THREADS"),1 )
    @variables(lp, begin
        0 <= x[1:i+j+i*j] <=1
    end)
    @objective( lp, Min, sum(x[n] for n in idx0) + sum(1-(x[m]) for m in idx1) )
    @constraint(lp, con1[b=1:j], sum(x[j*(a-1)+b] for a in 1:i) == x[i*j+i+b])
    @constraint(lp, con2[a=1:i,b=1:j], x[j*(a-1)+b]<=x[i*j+a])
    optimize!(lp)
    if termination_status(lp) == MOI.OPTIMAL
        return round.(JuMP.value.(x),digits=4)
    else
        return 0; #print("no new lp solution found");
    end
end

function mindist(sameval,i,j,L2,groups,cost,fixcost,demand,l,que) #find a point which has the minimum distance to the avg value of points in the same cuboid
    dis = Model(CPLEX.Optimizer)
    MOI.set(dis, MOI.RawParameter("CPX_PARAM_SCRIND"), false)
    MOI.set(dis, MOI.RawParameter("CPX_PARAM_THREADS"),1  )

    @variable(dis, 0<=x[1:i+j+i*j]<=1)
    @constraint(dis, con[k in sameval], x[k] == que[k] )
    cubemean = [mean(L2[groups[l],:][:,k]) for k=1:3]
    @constraint(dis, dot(x[1:i*j],transpose(Array{Float64}(cost)))== cubemean[1]) #transpose(Array{Float64}
    @constraint(dis, dot(x[i*j+1:i*j+i],transpose(Array{Float64}(fixcost))) == cubemean[2]) #yub[2]
    @constraint(dis, dot(x[i*j+i+1:end],transpose(Array{Float64}(demand))) == -cubemean[3]) #transpose(Array{Float64}
    @constraint(dis, con1[b=1:j], sum(x[j*(a-1)+b] for a in 1:i) == x[i*j+i+b])
    @constraint(dis, con2[a=1:i,b=1:j], x[j*(a-1)+b]<=x[i*j+a])
    @objective(dis, Min, 1)
    optimize!(dis)
    if JuMP.termination_status(dis) ==MOI.OPTIMAL
        return round.(JuMP.value.(x), digits=4)
    else
        return nothing
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

function fbcheck(txi,i,j)
    result = true
    #1st constraint checking
    for b=1:j
        if sum(txi[(a-1)*j+b] for a=1:i) == txi[i*j+i+b]
            continue
        else
            result = false
            return result
        end
    end

    #2nd constraint checking
    for a=1:i
        for b=1:j
            if txi[j*(a-1)+b]<=txi[i*j+a]
                continue
            else
                result = false
                return result
            end
        end
    end
    return result
end

function giveobjval(x,i,j,cost,fixcost,demand)
    return [dot(x[1:j*i],cost[:]),dot(x[j*i+1:j*i+i],fixcost[:]),-dot(x[j*i+i+1:end],demand[:]),]
end

function getlocation(x,i,j,fmin,steps,cost,fixcost,demand)
    L = giveobjval(x,i,j,cost,fixcost,demand)
    loca = [floor.(Int,((L[k]-pre.fmin[k])/pre.steps[k])+1) for k=1:3]
    return loca
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

# function Grouping(data::Data, pre::Prelim,vari::Variable) #find common elements(variable values) from solutions in the same cuboid
#     lpX = vari.lpX; i=data.i; j=data.j; cost=data.cost; fixcost=data.fixcost; demand=data.demand;
#     groups = pre.groups; groupkeys = pre.groupkeys; L2 = pre.L2
#     LPcount1 = 0; Xf= Dict(); candX = Dict(); loca_check = []; fixednum = [];
#
#     for l=1:length(groups)
#         glength = length(groups[l])
#         if glength==1 #if there is one point in a group
#             if (all.(x->trunc(x)==x, lpX[groups[l]])[1] == 1) #all varval are integer => insert to Xf
#                 Xf[groupkeys[l]] = lpX[groups[l][1]]
#                 # print("lp==ip \n")
#             else
#                 candX[groupkeys[l]] = lpX[groups[l]][1] #candX is a set of solutions that should be processed to be an integer solution
#                 # print("only lp sol in cuboid \n")
#             end
#         else
#             cogroup = lpX[groups[l]]
#             que = round.(hcat(cogroup...))
#             sameval = findall(x->x==true,[all(x->x==que[l][1], que[l,:]) for l=1:length(lpX)])
#             push!(fixednum,length(sameval));
#             newlpsol = mindist(sameval,i,j,L2,groups,cost,fixcost,demand,l,que)
#             LPcount1+=1
#             if newlpsol == nothing
#                 # print("there is no LPsolution \n")
#             else #when newlpsol is integer, check if all values are integer
#                 if (all.(x->trunc(x)==x, newlpsol) == 1) #all varval is integer => insert to Xf
#                     Xf[groupkeys[l]] = newlpsol
#                     # print("newlpsol==ip \n")
#                 else
#                     push!(loca_check, newlpsol);
#                     # print("newlpsol in cuboid is found\n")
#                     candX[groupkeys[l]] = newlpsol;
#                 end
#             end
#
#         end
#     end
#     return Xf,candX,loca_check,LPcount1,fixednum
# end
GroupingTime = @CPUelapsed Xf,cand,loca_check,LPcount,fixednum = Grouping(data,pre,vari)

# mindist ver.2
GroupingTime = @CPUelapsed Xf,cand,loca_check,LPcount,fixednum = Grouping3(data,pre,vari)


Xf1 = copy(Xf)
candX = vcat(collect(values(cand)),loca_check)
PFset = [giveobjval(collect(values(Xf1))[k],data.i,data.j,data.cost,data.fixcost,data.demand) for k=1:length(Xf1)]

function GroupFP(Xf,PFset,candX,LPcount,pre::Prelim,data::Data)
    Tabu = []; elaps=0; TimeLimit = 3600;
    cost=data.cost; fixcost=data.fixcost; demand=data.demand;fmin=pre.fmin; fmax=pre.fmax;steps=pre.steps
    algo_start = CPUtime_us()
    for k=1:length(candX)
        # print("===============Feasi Pump",k," th candidate sol ==============\n")
        x_t = candX[k];
        SearchDone = false
        itr = 1
        Max_itr = 100*data.i#data.i+data.j #data.i*data.j+data.i+data.j #max(count(x->0<x<1,x_t),1)  #maximum number of attempts => How to set
        while itr<Max_itr && SearchDone==false
            i=data.i; j=data.j;
            xi_t = round.(Int,x_t)
            fxi = giveobjval(xi_t,i,j,cost,fixcost,demand)
            if ( (fbcheck(xi_t,i,j) == true) && (dominated(fxi,PFset)==false) ) #checking feasibility and dominance   #z[1]<=fmax[1]-1 && z[2]<=fmax[2]-1 && z[3]<=fmax[3]-1)
                idx = getlocation(xi_t,i,j,fmin,steps,cost,fixcost,demand)

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
                        fxi = giveobjval(xi_t,i,j,cost,fixcost,demand)
                        if ( (fbcheck(xi_t,i,j) == true) && (dominated(fxi,PFset)==false) )  #z[1]<=fmax[1]-1 && z[2]<=fmax[2]-1 && z[3]<=fmax[3]-1 )
                            idx = getlocation(xi_t,i,j,fmin,steps,cost,fixcost,demand)
                            push!(PFset,fxi) #add new solval to PFset
                            if !haskey(Xf,idx)
                                Xf[idx] = xi_t
                            else
                                push!([Xf[idx]], xi_t)
                            end
                            SearchDone = true
                        end
                    end
                end
                if SearchDone == false
                    push!(Tabu,xi_t) #when break
                    idx = getlocation(xi_t,i,j,fmin,steps,cost,fixcost,demand)
                    x_t = fbsearch(idx,xi_t,steps,fmin,fmax,i,j,cost,fixcost,demand)
                    LPcount+=1
                    if x_t == 0 #when there's no new feasible lp sol
                        # print("no lp sol's found");
                        break
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
FPTime = @CPUelapsed Xf2,PFset2,Tabu,FPLPcount = GroupFP(Xf,PFset,candX,LPcount,pre,data)
P = convert.(Array{Int,1},collect(values(Xf2)))
Pcoordi = convert.(Array{Int,1},collect(keys(Xf2)))
Pobj = [giveobjval(P[k],data.i,data.j,data.cost,data.fixcost,data.demand) for k=1:length(P)]

#Filter dominated solutions
FPsol, FPPF = domFilter(P,Pobj)
Pz,Py,Px=[hcat(FPPF...)[i,:] for i=1:3]; dfP = Px,Py,Pz
matP = zeros(length(FPsol),m)
for i=1:length(FPsol)
    matP[i,:] = FPsol[i]
end

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



#####################Count real new solutions #############
genew = 0;
for i=1:length(FPsol)
    if FPsol[1] ∉ vari.lpX
        global genew+=1
    end
end
record1 = DataFrame(InitialInt = length(Xf1)+1, CandiX = length(candX),
    commonality = length(fixednum),avgfixedvar = mean(fixednum),FoundIPs=length(P)-(length(Xf1)+1) ,
    solvedLP=FPLPcount, totalsol=length(FPPF),CPUtime=GroupingTime+FPTime, newsol=genew
)


Solpair = Dict(); PFpair = Dict()
for i=1:length(FPPF)
    Solpair[FPcoordi[i]] = FPsol[i];
    PFpair[FPcoordi[i]] = FPPF[i]
end
keyunion = union(FPcoordi,pre.groupkeys)
searchdone = Dict()
for i = 1:length(keyunion)
    searchdone[keyunion[i]] = false
end

S = Dict();T = Dict()
for i in keyunion
    if i in FPcoordi
        S[i] = [ Solpair[i] ];T[i] = [ PFpair[i] ]
    else
        S[i] = [];T[i] = []
    end
end

mutable struct Problem
    model::JuMP.Model
    x::Array{VariableRef,2}
    y::Array{VariableRef,1}
    z::Array{VariableRef,1}
    ub::Array{VariableRef,1}
    lb::Array{VariableRef,1}
    dummy::VariableRef
end
function flpmodel(data::Data)::Problem
    flp = JuMP.Model(CPLEX.Optimizer)
    MOI.set(flp, MOI.RawParameter("CPX_PARAM_SCRIND"), false)
    MOI.set(flp, MOI.RawParameter("CPX_PARAM_THREADS"),1  )
    @variable(flp, x[1:data.i,1:data.j], Bin)
    @variable(flp, y[1:data.i], Bin)
    @variable(flp, z[1:data.j], Bin)
    @variable(flp, ub[1:3], Int)
    @variable(flp, lb[1:3], Int)
    @variable(flp, dummy==1)
    @constraint(flp, o1,sum(data.cost2[a][b]*x[a,b] for a=1:data.i for b=1:data.j) <=1)
    @constraint(flp, o2, sum(data.fixcost[a]*y[a] for a=1:data.i) <=1 )
    @constraint(flp, o3, -sum(data.demand[b]*z[b] for b=1:data.j) <= 1)
    @constraint(flp, con1[b=1:data.j], sum(x[a,b] for a in 1:data.i) == z[b])
    @constraint(flp, con2[a=1:data.i,b=1:data.j], x[a,b]<=y[a])

    @constraint(flp, y11, lb[1] <= sum(data.cost2[a][b]*x[a,b] for a=1:data.i for b=1:data.j) )
    @constraint(flp, y21,  lb[2] <= sum(data.fixcost[a]*y[a] for a=1:data.i) )
    @constraint(flp, y31,  lb[3] <= -sum(data.demand[b]*z[b] for b=1:data.j) )
    @constraint(flp, y12, sum(data.cost2[a][b]*x[a,b] for a=1:data.i for b=1:data.j) <= ub[1])
    @constraint(flp, y22, sum(data.fixcost[a]*y[a] for a=1:data.i) <= ub[2])
    @constraint(flp, y32, -sum(data.demand[b]*z[b] for b=1:data.j)<= ub[3])
    # @constraint(flp, ,  lb[1] <=sum(data.fixcost[a]*y[a] for a=1:data.i) )
    # @constraint(flp, lb2, lb[2] <= -sum(data.demand[b]*z[b] for b=1:data.j) )
    @objective(flp, Min, dummy)
    Problem(flp,x,y,z,ub,lb,dummy)
end


function submodel(S,T)
    key = collect(keys(T)); val = collect(values(T));
    for l=1:length(key)
        m = flpmodel(data);
        yub = [key[l][k]*pre.steps[k]+pre.fmin[k] for k=1:3]
        ylb = [(key[l]-[1,1,1])[k]*pre.steps[k]+pre.fmin[k] for k=1:3]
        JuMP.fix(m.ub[1], yub[1]); JuMP.fix(m.ub[2], yub[2]); JuMP.fix(m.ub[3], yub[3])
        JuMP.fix(m.lb[1], ylb[1]); JuMP.fix(m.lb[2], ylb[2]); JuMP.fix(m.lb[3], ylb[3])
        if val[l] != []
            vs = S[key[l]][1]
            vs1 = reshape(vs[1:data.i*data.j],(data.i,data.j))
            vs2 = vs[data.i*data.j+1:data.i*data.j+data.i]
            vs3 = vs[data.i*data.j+data.i+1:end]
            for i=1:length(val[l])-1
                append!(val[l][1],val[l][i+1])
            end
            v1 = val[l]; l1 = length(val[l][1])
        else
            v1 = [[0,0,0]]; l1 = 0
        end

        path = "/home/ak121396/Downloads/KirlikSayin2014/sub"
        fname = path*"/"*data.input[end-12:end-4]*"_$l"*"_hyb.lp"
        f2name = path*"/"*data.input[end-12:end-4]*"_$l"*"_hyb.txt"
        writemodel(m.model, "$fname")
        dt0 = readlines("$fname")
        dt1 = lowercase.(dt0)
        dt2 = replace.(dt1, ['_','[',']'] =>"")
        writedlm("$fname", dt2); writedlm("$f2name", [l1 v1]);
    end
end
submodel(S,T)

# function makesubmodel(S,T)
#     key = collect(keys(T)); val = collect(values(T));
#     for l=1:length(key)
#         m = flpmodel(data);
#         yub = [key[l][k]*pre.steps[k]+pre.fmin[k] for k=1:3]
#         ylb = [(key[l]-[1,1,1])[k]*pre.steps[k]+pre.fmin[k] for k=1:3]
#         JuMP.fix(m.ub[1], yub[1]); JuMP.fix(m.ub[2], yub[2]); JuMP.fix(m.ub[3], yub[3])
#         JuMP.fix(m.lb[1], ylb[1]); JuMP.fix(m.lb[2], ylb[2]); JuMP.fix(m.lb[3], ylb[3])
#         if val[l] != []
#             vs = S[key[l]][1]
#             vs1 = reshape(vs[1:data.i*data.j],(data.i,data.j))
#             vs2 = vs[data.i*data.j+1:data.i*data.j+data.i]
#             vs3 = vs[data.i*data.j+data.i+1:end]
#         end
#
#         path = "/home/k2g00/k2g3475/multiobjective/solvers/"
#         fname = path*"fpkirlik/"*data.input[end-12:end-7]*"/"*data.input[end-12:end-3]*"_$l"*"_hyb.lp"
#         writemodel(m.model, "$fname")
#         dt0 = readlines("$fname")
#         dt1 = lowercase.(dt0)
#         dt2 = replace.(dt1, ['_','[',']'] =>"")
#         writedlm("$fname", dt2)
#     end
# end
# makesubmodel(S,T)

# function runKS(S,T)
#     makesubmodel(S,T); val = collect(values(T));
#     path = "/home/k2g00/k2g3475/multiobjective/solvers/"
#     files = readdir(path*"fpkirlik/"*data.input[end-12:end-7])
#
#     for l=1:length(files)
#         if val[l] != [] #if there is a nd point in a group
#             for i=1:length(val[l])-1
#                 append!(val[l][1],val[l][i+1])
#             end
#             v1 = val[l][1]; l1 = length(val[l][1])
#             lpfile = path*"fpkirlik/"*data.input[end-12:end-7]*"/"*files[i]
#             run(`$path./main ./$lpfile $v1 $l1`)
#         else #if there is no nd point in a group
#             v1 = [0,0,0]; l1 = 0
#             lpfile = path*"fpkirlik/"*data.input[end-12:end-7]*"/"*files[i]
#             run(`$path./main ./$lpfile $v1 $l1`)
#         end
#     end
# end
# runKS(S,T)





ins = ARGS[1][end-12:end-7] #CPUtime recorded, naming after the instance
colname = ARGS[1][end-12:end-4]


# insertcols!(record1,8, Symbol("$colname")=>GroupingTime+FPTime)
GFPX=DataFrame(matP); GFPY=DataFrame(dfP);
CSV.write("/home/k2g00/k2g3475/multiobjective/solvers/ep+FP/GFPresults/GFP_"*"$ins"*"_2record.csv",record1, append=true, writeheader=false )#, delim=',' )
#CSV.write(ARGS[1]*"_GFP_X.csv",GFPX, append=true, writeheader=false)
# CSV.write(ARGS[1]*"_GFP_Y_.csv",GFPY, writeheader=false, delim=' ' )
print(colname," GroupFP Done!")


# function testGrouping(groups,groupkeys,lpX,steps,fmin,L2) #find common elements(variable values) from solutions in the same cuboid
#     LPcount1 = 0; cost=data.cost; fixcost=data.fixcost; demand=data.demand
#     Xf= Dict(); candX = Dict(); loca_check = [];
#     for l=1:length(groups)
#         glength = length(groups[l])
#         if glength==1 #if there is one point in a group
#             if (all.(x->trunc(x)==x, lpX[groups[l]])[1] == 1) #all varval are integer => insert to Xf
#                 Xf[groupkeys[l]] = lpX[groups[l][1]]
#             else
#                 candX[groupkeys[l]] = lpX[groups[l]][1] #candX is a set of solutions that should be processed to be an integer solution
#             end
#         else
#             cogroup = lpX[groups[l]]
#             que = round.(hcat(cogroup...))
#             sameval = findall(x->x==true,[all(x->x==que[l][1], que[l,:]) for l=1:length(lpX)])
#             newlpsol = mindist(sameval,data.i,data.j,L2,groups,cost,fixcost,demand,l,que)
#             LPcount1+=1
#
#
#             if newlpsol == nothing
#             else #when newlpsol is integer, check if all values are integer
#                 if (all.(x->trunc(x)==x, newlpsol) == 1) #all varval is integer => insert to Xf
#                     Xf[groupkeys[l]] = newlpsol
#                 else
#                     push!(loca_check, newlpsol);
#                     candX[groupkeys[l]] = newlpsol;
#                 end
#             end
#
#         end
#     end
#     # return Xf,candX,loca_check,LPcount1
# end
#initiate function (warming up)
# @CPUelapsed testGrouping(pre.groups[1:min(length(pre.groups),5)], pre.groupkeys[1:min(length(pre.groups),5)],vari.lpX,pre.steps,pre.fmin,pre.L2)


# function mindist2(sameval,i,j,L2,groups,cost,fixcost,demand,l,que,minpoint) #find a point which has the minimum distance to the avg value of points in the same cuboid
#     dis = Model(CPLEX.Optimizer)
#     MOI.set(dis, MOI.RawParameter("CPX_PARAM_SCRIND"), false)
#     MOI.set(dis, MOI.RawParameter("CPX_PARAM_THREADS"),1  )
#
#     @variable(dis, 0<=x[1:i+j+i*j]<=1)
#     @constraint(dis, con[k in sameval], x[k] == que[k] )
#     @constraint(dis, dot(x[1:i*j],transpose(Array{Float64}(cost)))>= minpoint[1]) #transpose(Array{Float64}
#     @constraint(dis, dot(x[i*j+1:i*j+i],transpose(Array{Float64}(fixcost))) >= minpoint[2]) #yub[2]
#     @constraint(dis, dot(x[i*j+i+1:end],transpose(Array{Float64}(demand))) >= -minpoint[3]) #transpose(Array{Float64}
#     @constraint(dis, con1[b=1:j], sum(x[j*(a-1)+b] for a in 1:i) == x[i*j+i+b])
#     @constraint(dis, con2[a=1:i,b=1:j], x[j*(a-1)+b]<=x[i*j+a])
#     @objective(dis, Min, 1)
#     optimize!(dis)
#     if JuMP.termination_status(dis) ==MOI.OPTIMAL
#         return round.(JuMP.value.(x), digits=4)
#     else
#         return nothing
#     end
# end
