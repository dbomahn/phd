using DelimitedFiles,DataFrames,JuMP,CPLEX,LinearAlgebra,CPUTime,CSV,MathOptFormat,MathProgBase
const MPB=MathProgBase

path1 = "/home/ak121396/Desktop/MIPLIB/"
# data = readdlm(path1*"neos-1620770.mps") #f = readdlm(dir1*"/acc-tight4_pre_varval.sol")
varval = round.(readdlm(path1*"/varval/neos-1620770_pre_img_p.sol"), digits=2) #round numbers till 2 digits
objval = round.(readdlm(path1*"/objval/neos-1620770_img_p.sol")[:,2:end], digits=2) #pfs[10*(fi-1)+sf]
vlp = readdlm(path1*"/vlp/neos-1620770.lp.vlp")

###############   Load lp file to get obj coefficents ###
function loadfile(filename,solver=CplexSolver(CPX_PARAM_SCRIND=0))
    model=buildlp([-1,0],[2 1],'<',1.5, solver) # create dummy model with correct solver
    MPB.loadproblem!(model,filename) # load what we actually want
    return model
end
mip = loadfile(path1*"/instances/neos-1620770.mps")
B = MPB.getconstrmatrix(mip)
# lb = MPB.getconstrLB(mip)
# ub = MPB.getconstrUB(mip)
# m,n=size(B)
s,ss = size(varval)


coef = zeros(3,ss)
for k=1:length(vlp)
    if vlp[k]=="o"
        coef[ vlp[k,:][2],vlp[k,:][3] ] = vlp[k,:][4]
    end
end



indx= findall(x->x!=0,[sum(varval[i,:]) for i=1:s]) #NZ sol
lpX = [varval[i,:] for i in indx]

pf = hcat(objval[:,2],objval[:,1],objval[:,3])
L = reshape([pf[i,:] for i in indx],length(lpX),1)
L2 = transpose(hcat(L...))
fmin = Array{Float64}([minimum(L2[:,i]) for i=1:3])
fmax = Array{Float64}([maximum(L2[:,i]) for i=1:3])
steps = [abs(fmax[i]-fmin[i])/length(L) for i=1:3] #determine steps according to #customers(j)
cube = Dict();
for iter=1:length(L)
    loca = [round.(Int,((L2[iter,k]-fmin[k])/steps[k])+1) for k=1:3]
    if !haskey(cube,loca)
        cube[loca] = [iter]
    else
        push!(cube[loca], iter)
    end
end
groups = collect(values(cube)); groupkeys = collect(keys(cube))
########################### TEST ########################
Grouping(groups[1:min(length(groups),5)],groupkeys[1:min(length(groups),5)],lpX,steps,fmin,L2)
########################################################
GroupingTime = @CPUelapsed Xf,cand,loca_check,LPcount = Grouping(groups,groupkeys,lpX,steps,fmin,L2)
Xf1 = copy(Xf)
# candX = filter(x->(last(x)!==nothing), candX)# candkeys = collect(keys(candX)); candvals = collect(values(candX));
candX = vcat(collect(values(cand)),loca_check)
PFset = giveobjval.(collect(values(Xf1)))
#Time Feasibility Pump
FPTime = @CPUelapsed Xf2,PFset2,Tabu,FPLPcount = mipFP(Xf1,PFset,candX,LPcount,steps,fmin,fmax)

############ adding more to jump_model
info = VariableInfo(true, 0, false, NaN, false, NaN, false, NaN, false, false)
JuMP.add_variable(miplib, JuMP.build_variable(error, info), "x")
JuMP.add_constraint()
JuMP.add_constraint(jump_model, build_constraint(_error, func, set))
@constraint(jump_model, ep1, 2x<=1)
