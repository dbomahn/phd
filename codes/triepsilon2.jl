using DelimitedFiles,DataFrames,JuMP,CPLEX,LinearAlgebra,CPUTime,CSV #,PlotlyJS,RDatasets,Colors
global TimeLimit = Inf
################################  Data  ####################################
mutable struct Daat{S<:String,I<:Int64,A<:Array{}}
    input::S
    i::I
    j::I
    # data::Array{}
    fixcost::A
    cost::A
    data.cost2::A
    function Data(input::String)
        data=readdlm(input[1])
        i,j = filter!(a->typeof(a)==Int, data[1,:]) # facility_i, customer_j # for FLP
        for k=2:i+1
            fix = filter!(a->typeof(a)==Int, data[k,:])
            push!(fixcost,fix)
        end
        data.demand = filter!(a->typeof(a)==Int, data[2+i,:])
        a,b = size(data)
        for k=i+3:a
            ct = filter!(a->typeof(a)==Int, data[k,:])
            push!(data.cost2,ct)
            for l=1:length(ct)
                push!(cost,ct[l])
            end
        end
        new{S,I,A}(input,i,j,fixcost,cost,data.cost2)
    end
end

data=Data("/home/ak121396/multiobjective/instances/triFLP/instances/05_010_01.txt") #05_010_09.txt
demand = data.demand; fixcost=data.fixcost; cost2=data.cost;cost=data.cost

# global i,j = filter!(a->typeof(a)==Int, data[1,:])
# (global fixcost,capa,cost,data.cost2 = [],[],[],[] )
# for k=2:i+1
#     fix,cp = filter!(a->typeof(a)==Int, data[k,:])
#     push!(fixcost,fix)
#     push!(capa,cp)
# end
# global data.demand = filter!(a->typeof(a)==Int, data[2+i,:])
# a,b = size(data)
# for k=i+3:a
#     ct = filter!(a->typeof(a)==Int, data[k,:])
#     push!(data.cost2,ct)
#     for l=1:length(ct)
#         push!(cost,ct[l])
#     end
# end
################################   Functions  ##################################

##############  tri- epsilon ################
function getConstraints(i,e,P)
    ϵ = [0,0];
    ϵ′ = [Inf,Inf];
    for j=2:3
        d = round(Int,i%(length(P)+1))
        i = round(Int,(i-d)/(length(P)+1))
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
    if x!=nothing
        return [-dot(x[i+1:j+i],demand[:]),dot(x[1:i],fixcost[:]),dot(x[j+i+1:end],(cost[:]))]
    else
        return nothing
    end
end

function dominated(x,P)
    st = false
    if x==nothing
        return true
    else
        for k=1:length(P)
            if all( x .>= P[k])#&& any(x > P[k])
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
                copyobj[i]=nothing; copysol[i]=nothing; print(i,"dom by ", j,"\n");break
            elseif all(obj[j] .>= obj[i]) == true
                copyobj[j]=nothing; copysol[j]=nothing; print(j,"dom by ", i,"\n")
            end
        end
    end

    finalsol = filter!(a->a!=nothing, collect(values(copysol)))
    finalobj = filter!(a->a!=nothing, collect(values(copyobj)))

    return finalsol,finalobj
end

###################FLP model ###############
lex1 = Model(CPLEX.Optimizer   ) #with_optimizer(CPLEX.Optimizer)
MOI.set(lex1, MOI.RawParameter("CPX_PARAM_SCRIND"), false )
MOI.set(lex1, MOI.RawParameter("CPX_PARAM_THREADS"),1  )
@variable(lex1, x1[1:i,1:j], Bin)
@variable(lex1, y1[1:i], Bin)
@variable(lex1, z1[1:j], Bin)
@variable(lex1, ep[1:2], Int)
@variable(lex1, ep′[1:2], Int)
@constraint(lex1, con1[b=1:j], sum(x1[a,b] for a in 1:i) == z1[b])
@constraint(lex1, con2[a=1:i,b=1:j], x1[a,b]<=y1[a])
@constraint(lex1, epsilon11,  sum(data.fixcost[a]*y1[a] for a=1:i) <= ep′[1]-1)
@constraint(lex1, epsilon12,  ep[1] <=sum(data.fixcost[a]*y1[a] for a=1:i) )
@constraint(lex1, epsilon21, sum(data.cost2[a][b]*x1[a,b] for a=1:i for b=1:j) <= ep′[2]-1)
@constraint(lex1, epsilon22, ep[2] <= sum(data.cost2[a][b]*x1[a,b] for a=1:i for b=1:j) )
@objective(lex1, Min, -sum(data.demand[b]*z1[b] for b=1:j))
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
@constraint(lex2, con3, -sum(z2[b]*data.demand[b] for b=1:j) <= obj1)
@constraint(lex2, epsilon11,   sum(data.fixcost[a]*y2[a] for a=1:i) <= ep2′[1]-1 )
@constraint(lex2, epsilon12,  ep2[1] <= sum(data.fixcost[a]*y2[a] for a=1:i) )
@constraint(lex2, epsilon21, sum(data.cost2[a][b]*x2[a,b] for a=1:i for b=1:j) <= ep2′[2]-1 )
@constraint(lex2, epsilon22,  ep2[2] <= sum(data.cost2[a][b]*x2[a,b] for a=1:i for b=1:j) )
@objective(lex2, Min, sum(data.fixcost[a]*y2[a] for a=1:i) )
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
@constraint(lex3, con3, -sum(z3[b]*data.demand[b] for b=1:j) <= obj2 )
@constraint(lex3, con4 , sum(data.fixcost[a]*y3[a] for a=1:i) <= obj3 )
@constraint(lex3, epsilon11,  sum(data.fixcost[a]*y3[a] for a=1:i) <= ep3′[1]-1 )
@constraint(lex3, epsilon12,  ep3[1] <= sum(data.fixcost[a]*y3[a] for a=1:i))
@constraint(lex3, epsilon21,  sum(data.cost2[a][b]*x3[a,b] for a=1:i for b=1:j) <= ep3′[2]-1 )
@constraint(lex3, epsilon22,  ep3[2] <= sum(data.cost2[a][b]*x3[a,b] for a=1:i for b=1:j)  )
@objective(lex3, Min, sum(data.cost2[a][b]*x3[a,b] for a=1:i for b=1:j))
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
#############################  Tri-epsilon ######################
function triepsilon()
    E = []; F = []; global e = Dict(2=>[1,Inf],3=>[1,Inf]); fval = [];
    m = 0; c = 1; global solvedIP=0;

    algo_start = CPUtime_us()
    while true
        while true
            # print("\n======== m,c = ", m," ",c,"  ========\n" )
            if m>=c || (CPUtime_us()-algo_start)>=TimeLimit*1e6
                return E,solvedIP

            else
                # @goto cont
            end
            # @label cont
            global (ϵ,ϵ′) = getConstraints(m,e,E)
            # print("new epsilons= ",ϵ,ϵ′,"\n")
            if [ϵ,ϵ′] ∉ F
                global (s,objs,solvedIP) = opt(ϵ,ϵ′,solvedIP)
                if ( s == nothing || dominated(objs,fval)==true)
                    F = union(F,[[ϵ,ϵ′]]);
                else
                    fval = union(fval,[objs]);
                    break #
                end
            end
            # print("already searched \n")
            m+=1
        end
        # print("added obj= ",objs,"\n")
        E = union(E, [s]);
        F = union(F, [ [objs[2],objs[3]], [ϵ′[1],ϵ′[2]] ] )
        updateConstraints(objs,e)
        c = (length(E)+1)^2
        m = 0
    end
    return E,solvedIP
end

epTime = @CPUelapsed epsol,epIPcount = triepsilon()
epPF=[giveobjval(epsol[j]) for j=1:length(epsol)]

# Filter dominated solutions
ϵsol,ϵPF = domFilter(epsol,epPF)
epz,epy,epx=[hcat(ϵPF...)[i,:] for i=1:3]; dfE = epx,epy,epz;
matE = zeros(length(ϵsol),i+j+i*j)
for i=1:length(ϵsol)
    matE[i,:] = ϵsol[i]
end

ins = ARGS[1][end-12:end-7] #CPUtime recorded, naming after the instance
colname = ARGS[1][end-12:end-4]
record1 = DataFrame(solvedIP=epIPcount,sols=length(epz))
insertcols!(record1,3, Symbol("$colname")=>epTime)
epX=DataFrame(matE); epY=DataFrame(dfE);

CSV.write("/home/k2g00/k2g3475/multiobjective/solvers/ep+FP/ep_results/ep_2hr_"*"$ins"*"_record.csv",record1, append=true, writeheader=false)
#CSV.write(ARGS[1]*"ep_X.csv", epX, append=true, writeheader=false)
CSV.write(ARGS[1]*"ep_2hr_Y.csv",epY, writeheader=false, delim=' ' )
print(colname,"Tri epsilon done")
