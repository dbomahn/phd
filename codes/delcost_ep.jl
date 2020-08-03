using DelimitedFiles,DataFrames,JuMP,CPLEX,LinearAlgebra,CPUTime,CSV #,PlotlyJS,RDatasets,Colors
################################  Data  ####################################
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

data=Data(ARGS[1])
################################   Functions  ##################################
function getConstraints(i,e,P)
    ϵ = [0,0];
    ϵ′ = [0,0];
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
function giveobjval(x,FLP::Data)
    if x!=nothing
        return [dot(x[FLP.j+FLP.i+1:end],(FLP.cost[:])),dot(x[1:FLP.i],FLP.fixcost[:]),-dot(x[FLP.i+1:FLP.j+FLP.i],FLP.demand[:])]
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
            if all( x .>= P[k]) #&& any(x > P[k])
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
mutable struct Problem
    model::JuMP.Model
    x::Array{VariableRef,2}
    y::Array{VariableRef,1}
    z::Array{VariableRef,1}
    ep::Array{VariableRef,1}
    ep′::Array{VariableRef,1}
    obj1::VariableRef
    obj2::VariableRef
end
function FLPlex(FLP::Data)::Problem
    flp = JuMP.Model(with_optimizer(CPLEX.Optimizer))
    MOI.set(flp, MOI.RawParameter("CPX_PARAM_SCRIND"), false)
    MOI.set(flp, MOI.RawParameter("CPX_PARAM_THREADS"),1  )
    @variable(flp, x[1:FLP.i,1:FLP.j], Bin)
    @variable(flp, y[1:FLP.i], Bin)
    @variable(flp, z[1:FLP.j], Bin)
    @variable(flp, ep[1:2], Int)
    @variable(flp, ep′[1:2], Int)
    @variable(flp, obj1)
    @variable(flp, obj2)
    @constraint(flp, con1[b=1:FLP.j], sum(x[a,b] for a in 1:FLP.i) == z[b])
    @constraint(flp, con2[a=1:FLP.i,b=1:FLP.j], x[a,b]<=y[a])
    @constraint(flp, con3, sum(FLP.cost2[a][b]*x[a,b] for a=1:FLP.i for b=1:FLP.j) <= obj1 )
    @constraint(flp, con4 , sum(FLP.fixcost[a]*y[a] for a=1:FLP.i) <= obj2 )
    @constraint(flp, epsilon11,  sum(FLP.fixcost[a]*y[a] for a=1:FLP.i) <= ep′[1]-1)
    @constraint(flp, epsilon12,  ep[1] <=sum(FLP.fixcost[a]*y[a] for a=1:FLP.i) )
    @constraint(flp, epsilon21, -sum(FLP.demand[b]*z[b] for b=1:FLP.j) <= ep′[2]-1)
    @constraint(flp, epsilon22, ep[2] <= -sum(FLP.demand[b]*z[b] for b=1:FLP.j) )
    Problem(flp,x,y,z,ep,ep′,obj1,obj2)
end
function lex1(lex::Problem,FLP::Data)
    @objective(lex.model,Min,sum(FLP.cost2[a][b]*lex.x[a,b] for a=1:FLP.i for b=1:FLP.j))
    Problem(lex.model,lex.x,lex.y,lex.z,lex.ep,lex.ep′,lex.obj1,lex.obj2)
end
function lex2(lex::Problem,FLP::Data)
    @objective(lex.model,Min,sum(FLP.fixcost[a]*lex.y[a] for a=1:FLP.i))
    Problem(lex.model,lex.x,lex.y,lex.z,lex.ep,lex.ep′,lex.obj1,lex.obj2)
end
function lex3(lex::Problem,FLP::Data)
    @objective(lex.model,Min,-sum(FLP.demand[b]*lex.z[b] for b=1:FLP.j))
    Problem(lex.model,lex.x,lex.y,lex.z,lex.ep,lex.ep′,lex.obj1,lex.obj2)
end


lexi1 = lex1(FLPlex(data),data);lexi2 = lex2(FLPlex(data),data);lexi3 = lex3(FLPlex(data),data)
optimize!(lexi1.model) #test run

function opt(lex1::Problem,lex2::Problem,lex3::Problem,FLP::Data,ϵ,ϵ′,solvedIP)
    JuMP.fix(lex1.ep[1], ϵ[1]); JuMP.fix(lex1.ep[2], ϵ[2])
    JuMP.fix(lex1.ep′[1],ϵ′[1]); JuMP.fix(lex1.ep′[2], ϵ′[2])
    optimize!(lex1.model)

    if termination_status(lex1.model) == MOI.OPTIMAL
        solvedIP+=1; val1 = objective_value(lex1.model)
        JuMP.fix(lex2.ep[1], ϵ[1]); JuMP.fix(lex2.ep[2], ϵ[2])
        JuMP.fix(lex2.ep′[1],ϵ′[1]); JuMP.fix(lex2.ep′[2], ϵ′[2])
        JuMP.fix(lex2.obj1, val1);
        optimize!(lex2.model)

        if (termination_status(lex2.model) == MOI.OPTIMAL)
            solvedIP+=1; val2 = objective_value(lex2.model)
            JuMP.fix(lex3.ep[1], ϵ[1]); JuMP.fix(lex3.ep[2], ϵ[2])
            JuMP.fix(lex3.ep′[1],ϵ′[1]); JuMP.fix(lex3.ep′[2], ϵ′[2])
            JuMP.fix(lex3.obj1, val1); JuMP.fix(lex3.obj2, val2)
            optimize!(lex3.model)
            if (termination_status(lex3.model) == MOI.OPTIMAL)
                solvedIP+=1
                lexY = round.(Int,value.(lex3.y)); lexZ = round.(Int,value.(lex3.z)); lexX = round.(Int,value.(lex3.x));
                s =  append!( append!(lexY,lexZ),transpose(lexX) )
                return s,giveobjval(s,FLP),solvedIP
            else
                # print("Lex3: ",termination_status(lex3.model),"\n")
                return nothing,nothing,solvedIP
            end
        else
            # print("Lex2: ",termination_status(lex2.model),"\n")
            return nothing,nothing,solvedIP
        end
    else
        # print("Lex1: ",termination_status(lex1.model),"\n")
        return nothing,nothing,solvedIP
    end
end

#############################  Tri-epsilon ######################

mutable struct Laumanns
    E::Array{}; F::Array{}; fval::Array{}
    e::Dict
    m::Int; c::Int; solvedIP::Int
    objs::Array{}; s::Array{}
    ϵ::Array{}; ϵ′::Array{}
    Times::Array{}
    function Laumanns()
        TimeLimit = Inf
        E = []; F = []; e = Dict(2=>[0,9999999999],3=>[-99999999999,0]); fval = [];
        m = 0; c = 1; solvedIP=0; ϵ=[]; ϵ′=[]; objs=[]; s =[]; #Times=[0.0,0.0,0.0,0.0,0.0]
        algo_start = CPUtime_us()
        while true
            while true
                # print("\n======== m,c = ", m," ",c,"  ========\n" )
                if m>=c || (CPUtime_us()-algo_start)>=TimeLimit*1e6
                    return E,solvedIP#,Times
                else
                end
                # Time_getCon = @CPUelapsed
                ϵ,ϵ′ = getConstraints(m,e,E)
                # Times[1] = Times[1]+Time_getCon
                # print("new epsilons= ",ϵ,ϵ′,"\n")
                if [ϵ,ϵ′] ∉ F
                    # Time_opt= @CPUelapsed
                    s,objs,solvedIP = opt(lexi1,lexi2,lexi3,data,ϵ,ϵ′,solvedIP)
                    # Times[2] = Times[2]+Time_opt
                    # Time_domi = @CPUelapsed
                    dominance = dominated(objs,fval)
                    # Times[3] = Times[3]+Time_domi
                    if ( s == nothing || dominance == true) #dominated(objs,fval)==true
                        F = union(F,[[ϵ,ϵ′]]);
                    else
                        fval = union(fval,[objs]);
                        break
                    end
                end
                # print("already searched \n")
                m+=1
            end
            # print("added obj= ",objs,"\n")
            CPUtic()
            E = union(E, [s]);
            F = union(F, [ [objs[2],objs[3]], [ϵ′[1],ϵ′[2]] ] )
            # Time_updateCon = @CPUelapsed dummy =
            updateConstraints(objs,e)
            # Times[4] = Times[4]+Time_updateCon
            c = (length(E)+1)^2
            m = 0
            # Times[5] = Times[5]+ CPUtoq()
        end
        return E,solvedIP#,Times
    end
end
epTime = @CPUelapsed epsol,epIPcount = Laumanns()
epPF=[giveobjval(epsol[j],data) for j=1:length(epsol)]

# Filter dominated solutions
i = data.i; j=data.j
ϵsol,ϵPF = domFilter(epsol,epPF)
epz,epy,epx=[hcat(ϵPF...)[i,:] for i=1:3]; dfE = epx,epy,epz;
matE = zeros(length(ϵsol),i+j+i*j)
for i=1:length(ϵsol)
    matE[i,:] = ϵsol[i]
end

ins = ARGS[1][end-12:end-7] #CPUtime recorded, naming after the instance
colname = ARGS[1][end-12:end-4]
record1 = DataFrame(solvedIP=epIPcount,sols=length(epz))  #, getCon=Timerecord[1], opt=Timerecord[2],dom=Timerecord[3],updatCon=Timerecord[4],union = Timerecord[5])
insertcols!(record1,3, Symbol("$colname")=>epTime)
epX=DataFrame(matE); epY=DataFrame(dfE);

CSV.write("/home/k2g00/k2g3475/multiobjective/solvers/ep+FP/ep_results/ep_inf_delcost_"*"$ins"*"_record.csv",record1, append=true, writeheader=false)
#CSV.write(ARGS[1]*"ep_del_record.csv",record1 ,writeheader=false)
#CSV.write(ARGS[1]*"ep_X.csv", epX, append=true, writeheader=false)
CSV.write(ARGS[1]*"ep_inf_del_Y.csv",epY, writeheader=false, delim=' ' )
print(colname,"Tri epsilon done")
