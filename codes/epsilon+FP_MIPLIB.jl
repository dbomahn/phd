# To provide good quality of integer solutions to tri_epsilon method, FP proceed first.

using DelimitedFiles,DataFrames,JuMP,CPLEX,MathOptFormat,LinearAlgebra,CPUTime
# dir1 = "C:\\Users\\AK121396\\Desktop\\FLP\\"
# dir2 = dir1*"PF\\"
# include("C:\\Users\\AK121396\\Desktop\\Dropbox\\JKU\\Julia\\functions.jl")
# folders = readdir(dir1)[12]
# subf = readdir(dir1*folders)# [fi])
# f = readdlm(dir1*folders[:]*"/"*subf[4])

include("/home/ak121396/Dropbox/JKU/Julia/functions.jl")
dir1 = "/home/ak121396/Desktop/phd/MIPLIB/"; dir2 = dir1*"objval/"
sols = readdlm(dir1*"/varval/"*"neos-1620770_pre_img_p.sol")
s,ss = size(sols)
indx= findall(x->x!=0,[sum(sols[i,:]) for i=1:s]) #NZ sol
lpX = [sols[i,:] for i in indx]
pf = round.(readdlm(dir2*"neos-1620770_img_p.sol")[:,2:end]; digits=2)
pf = hcat(pf[:,2],pf[:,1],pf[:,3])
L = reshape([pf[i,:] for i in indx],length(lpX),1)
L2 = transpose(hcat(L...))
fmin = Array{Int}([minimum(L2[:,i]) for i=1:3])
fmax = Array{Int}([maximum(L2[:,i]) for i=1:3])

unique(L2)

# steps = [round.(Int,abs(fmax[i]-fmin[i])/length(L)) for i=1:3] #determine steps according to #customers(j)
cube = Dict();
for iter=1:length(L)
    loca = [round.(Int,((L2[iter,k]-fmin[k])/steps[k])+1) for k=1:3]
    if !haskey(cube,loca)
        cube[loca] = [iter]
    else
        push!(cube[loca], iter)
    end
end
groups = collect(values(cube)); groupkeys = collect(keys(cube));
GroupingTime = @CPUelapsed Xf,candX,LPcount = GroupFP(groups,groupkeys,lpX,steps,fmin,L2)
Xf1 = copy(Xf)
candX = filter(x->(last(x)!==nothing), candX)

#Time Feasibility Pump
FPTime = @CPUelapsed Xf2,Tabu,LPcount2 = FeasPump(Xf1,candX,LPcount,steps,fmin,fmax)

P = collect(values(Xf2));
Pcoordi = collect(keys(Xf2))
Pobj = giveobjval.(P)



function opt(f,ϵ′,ϵ″)
    lex1 = Model( with_optimizer(CPLEX.Optimizer))
    @variable(lex1, x1[1:i,1:j], Bin)
    @variable(lex1, y1[1:i], Bin)
    @variable(lex1, z1[1:j], Bin)
    @constraint(lex1, con1[b=1:j], sum(x1[a,b] for a in 1:i) == z1[b])
    @constraint(lex1, con2[a=1:i,b=1:j], x1[a,b]<=y1[a])
    @constraint(lex1, epsilon1,  ϵ[1] <= sum(fixcost[a]*y1[a] for a=1:i) <= ϵ′[1])
    @constraint(lex1, epsilon2,  ϵ[2] <= sum(cost2[a][b]*x1[a,b] for a=1:i for b=1:j) <= ϵ′[2])
    @objective(lex1, Min, -sum(demand[b]*z1[b] for b=1:j))
    optimize!(lex1)

    if termination_status(lex1)==MOI.OPTIMAL
        lex2 = Model( with_optimizer(CPLEX.Optimizer))
        UB3 = 10;
        @variable(lex2, x2[1:i,1:j], Bin)
        @variable(lex2, y2[1:i], Bin)
        @variable(lex2, z2[1:j], Bin)
        @constraint(lex2, con1[b=1:j], sum(x2[a,b] for a in 1:i) == z2[b])
        @constraint(lex2, con2[a=1:i,b=1:j], x2[a,b]<=y2[a])
        @constraint(lex2, con3, -sum(z2[b]*demand[b] for b=1:j) <= objective_value(lex1))
        @constraint(lex2, epsilon1,  ϵ[1] <= dot(fixcost,transpose(y2)) <= ϵ′[1])
        @constraint(lex2, epsilon2,  ϵ[2] <= sum(cost2[a][b]*x2[a,b] for a=1:i for b=1:j) <= ϵ′[2])
        @objective(lex2, Min, sum(fixcost[a]*y2[a] for a=1:i)+ sum(cost2[a][b]*x2[a,b] for a=1:i for b=1:j)/UB3 )
        optimize!(lex2)
        objective_value(lex2)
        lex2y = round.(Int,value.(y2)); lex2z = round.(Int,value.(z2)); lex2x = round.(Int,value.(x2));
        # newf = JuMP.objective_value(lex2)
        if termination_status(lex2) == MOI.OPTIMAL
            if ( giveobjval(append!(append!(lex2y,lex2z), lex2x))∉f )
                lex2y = round.(Int,value.(y2)); lex2z = round.(Int,value.(z2)); lex2x = round.(Int,value.(x2));
                return append!(append!(lex2y,lex2z), lex2x)
            end
        else
            return nothing
        end
    else
        return nothing
    end
end


function triepsilon(P,Pobj,Pcoordi,limit)
    EE = [];
    for k=1:length(P)
        global E = []; F = [];
        m = 0; c = 4; solvedLP = 0;
        fval = Pobj[k]
        temp1 = (fmin+steps.*Pcoordi[k])[2:3]
        temp2 = (fmin+steps.*(Pcoordi[k]+[1,1,1]))[2:3]
        union(E,[P[k]]);
        global e = Dict( 2=>[temp1[1],temp2[1]], 3=>[temp1[2],temp2[2]] )
        @show updateConstraints(fval,e)
        while solvedLP<=limit
            while true
                if m>=c
                    return E; print("Done. \n")
                end
                global (ϵ,ϵ′) = getConstraints(m,e,E)
                if [ϵ,ϵ′] ∉ F
                    global s = opt(fval,ϵ,ϵ′)
                    if ( s == nothing || dominated(s,E)==true)
                        F = union(F,[[ϵ,ϵ′]]);
                    else
                        break
                    end
                end
                m+=1
            end
            E = union(E, [s]);
            F = union(F, [ [giveobjval(s)[2],giveobjval(s)[3]], [ϵ′[1],ϵ′[2]] ] )
            updateConstraints(giveobjval(s),e)
            c = (length(E)+1)^2
            m = 0
            solvedLP+=2
        end
        EE = union(EE,E);
    end
    return EE
end

# result=triepsilon(P,Pobj,Pcoordi,1000) #limit is set temporaryily
# PF=[giveobjval(result[j]) for j=1:length(result)]
# z1,y1,x1=[hcat(PF...)[i,:] for i=1:3]
