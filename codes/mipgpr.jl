using DataStructures,DelimitedFiles,DataFrames,JuMP,CPLEX,LinearAlgebra,CPUTime,CSV,Statistics,MathProgBase,MathOptInterface,MathOptFormat,StatsBase,JLD2,SparseArrays
const MPB = MathProgBase; const MOF = MathOptFormat; const MOI = MathOptInterface;
# @show ARGS[1]
# function loadlp(filename,solver=CplexSolver(CPX_PARAM_SCRIND=0))
#     model=buildlp([-1,0],[2 1],'<',1.5, solver) # create dummy model with correct solver
#     MPB.loadproblem!(model,filename) # load what we actually want
#     return model
# end

mutable struct Data
    mpsfile::String; vlpfile::String
    m::Int; n::Int; C::Array{}; B::Array{}; RHS::Dict{}; signs::Array{};
    function Data(mpsfile::String, vlpfile::String)
        mpsmodel = CPLEX.CplexMathProgModel(); MPB.loadproblem!(mpsmodel,mpsfile)
        B = MPB.getconstrmatrix(mpsmodel); m,n=size(B)
        lb = MPB.getconstrLB(mpsmodel);ub = MPB.getconstrUB(mpsmodel);
        RHS = Dict()
        for i=1:m
            if ub[i]==Inf
                RHS[i] = lb[i]
            else
                RHS[i] = ub[i]
            end
        end
        signs = []
        for i=1:m
            if ub[i] == Inf
                push!(signs,"l")
            elseif lb[i] == -Inf
                push!(signs,"u")
            else
                push!(signs, "s")
            end
        end
        vlp = readdlm(vlpfile); os = findall(i->vlp[i,1]=="o", 1:length(vlp[:,1]))
        inio = os[1]; endo = os[end]; Pmtx = vlp[inio:endo,2:4]; C = zeros(3,n)
        for k=1:length(Pmtx[:,1]) #(len,id) enumerate
            x = Pmtx[k,:][1]; y = Pmtx[k,:][2];
            C[x,y] = Pmtx[k,:][3]
        end
        new(mpsfile,vlpfile,m,n,C,B,RHS,signs)
    end
end

mutable struct valu
    x::String; y::String; dvar::Array{}; LB::Array{}; LBmtx::Array{}; #j::Int;radius::Float64
    function valu(x,y)
        JLD2.@load x dv
        dv1 = round.(dv; digits=4) #readdlm(x)
        objs = round.(digits=4, readdlm(y))#[:,1:end]
        ind = findall(i-> 0 in objs[i,:]  , 1:size(objs)[1])
        dv2 = dv1[setdiff(1:end, ind), :];
        LBmtx = objs[setdiff(1:end, ind), 2:end];
        dvar = [Vector(dv2[i, :]) for i = 1:size(LBmtx)[1]]
        LB = [Vector(LBmtx[i, :]) for i = 1:size(LBmtx)[1]]
        # maxobj = [maximum(LBmtx[:,i]) for i=1:3]; minobj = [minimum(LBmtx[:,i]) for i=1:3];
        # steps = round.(Int,(maxobj-minobj)/j)  #length(LB[:,1])
        # interval = round(Int,mean(maxobj-minobj)/j)
        # points = transpose(LB)
        # α = max(round(length(LB)*0.05),1)
        # radius = interval*α
        new(x,y,dvar,LB,LBmtx)#j,radius)
    end
end
mutable struct Valu
    x::String; y::String; roundv::Array{}; LB::Array{}; C::Array{}; #j::Int;radius::Float64
    function Valu(x,y,C)
        # dv = readdlm(x)
        JLD2.@load x dv
        c1 = ceil.(dv); c2 = [Vector(c1[i, :]) for i = 1:size(dv)[1]]
        f1 = floor.(dv); f2 = [Vector(f1[i, :]) for i = 1:size(dv)[1]]
        r1 = round.(dv); r2 = [Vector(r1[i, :]) for i = 1:size(dv)[1]]
        rdv = union(c2,f2,r2); rdv2 = unique!(rdv)
        idx = []
        for j=1:length(rdv2)
            if all(i->i==0, rdv2[j])==true
                push!(idx,j)
            end
        end
        roundv = rdv2[setdiff(1:end, idx)]
        LB = []
        for i=1:length(roundv)
            roundLB = [dot(C[j,:], roundv[i]) for j=1:3]
            push!(LB,roundLB)
        end
        new(x,y,roundv,LB,C)
    end
end
function getobjval(x,C)
    return [dot(x,C[1,:]),dot(x,C[2,:]),dot(x,C[3,:])]
end
function createNB(SI,C,dif,exploredSI)
    neibour = []; neiobj = [];
    for i in dif
        cpSI = copy(SI)
        if cpSI[i] == round(cpSI[i]) #if vari is int value
            if cpSI[i] == 1
                cpSI[i] = 0
            else
                cpSI[i] = 1
            end
            push!(neibour, cpSI);
            push!(neiobj, getobjval(cpSI,C))
        else #if vari is fractional
            cpSI[i] = 1; push!(neibour, cpSI); push!( neiobj, getobjval(cpSI,C) )
            cpSI = copy(SI)
            cpSI[i] = 0; push!(neibour, cpSI); push!( neiobj, getobjval(cpSI,C) )
        end
    end
    idx = findall(i-> neibour[i] in exploredSI, 1:length(neibour))
    neibour = setdiff(neibour,exploredSI)
    neiobj = neiobj[setdiff(1:end, idx),:]
    return neibour,neiobj
end
function FBcheck(xx,n)
    for k=1:n
        JuMP.fix(x[k],xx[k])
    end
    optimize!(mip)
    # print("status: ", termination_status(mip), "\n" )
    if termination_status(mip) == MOI.OPTIMAL
        return true
    else
        return false
    end
end
function nextSI(neibour,neiobj,C,SI)
    SIobj = getobjval(SI,C)
    for i=1:length(neiobj)
        if length(neiobj) == 1  #if there is one candiate sol
            return neibour[1]#,neiobj[1]
        else length(neiobj) > 1 # if there are multiple candiates, check the improved ratio
            ratiotb = zeros(length(neiobj),3)
            for i=1:length(neiobj)
                ratiotb[i,:] = neiobj[i]./SIobj
            end
            ranktb = zeros(length(neiobj),3)
            for i=1:3
                ranktb[:,i] = tiedrank(ratiotb[:,i])
            end
            ranksum = [sum(ranktb[i,:]) for i=1:length(neiobj)]
            mostimp = findall(x-> x == minimum(ranksum), ranksum)
            k = rand(mostimp)[1]
            return neibour[k]#, neiobj[k]
        end
    end
end
function Postpro(dvar,LB,newsol)
    #Filter integer feasible sols from the initial solution set
    initX = []; initY = [];
    for i=1:length(roundv)-newsol
        if (FBcheck(roundv[i],n) == true)
            push!(initX,roundv[i]); push!(initY,LB[i]);
        end
    end
    P = union(initX,dvar[end-newsol+1:end]); Pobj = union(initY,LB[end-newsol+1:end])
    #Filter dominated solutions
    copysol = Dict(); copyobj = Dict();
    for i=1:length(Pobj)
        copysol[i] = P[i]; copyobj[i] = Pobj[i]
    end
    for i=1:length(Pobj)-1
        for j=i+1:length(Pobj)
            if all(Pobj[i] .>= Pobj[j]) == true #dominated by PF[j]
                copyobj[i]=nothing; copysol[i]=nothing; break
            elseif all(Pobj[j] .>= Pobj[i]) == true
                copyobj[j]=nothing; copysol[j]=nothing;
            end
        end
    end
    finalsol = filter!(a->a!=nothing, collect(values(copysol)))
    finalobj = filter!(a->a!=nothing, collect(values(copyobj)))
    return finalsol,finalobj
end
function GPR(C,n,dvar,LB,TL)
    IGPair=[]; exploredSI = []; newsol=0; t0=time(); # nbtime = 0; SItime = 0; FBtime = 0;
    for i=1:1000#n*10
    	if time()-t0 >= TL
            break
        end
        I,G = sample(1:length(dvar), 2, replace=false)
        SI = dvar[I]; SG = dvar[G]; infeas=0;
        while all.(SI != SG) && [I,G]∉IGPair && infeas<20 && (time()-t0<TL)
            dif = findall(i-> SI[i]!=SG[i], 1:n)
            neibour,neiobj = createNB(SI,C,dif,exploredSI) # nbtime = nbtime + @CPUelapsed
            if length(neiobj)==0
                break
            else
                for l=1:length(neiobj)
                    feasi = FBcheck(neibour[l],n) #FBtime = FBtime + @CPUelapsed
                    if feasi == true && neibour[l]∉ dvar
                        push!(dvar, neibour[l]); push!(LB, neiobj[l]);
                        newsol+=1
                        print(" new sol added \n")
                    else
                        print(" infeasible \n")
                    end
                end
                SI = nextSI(neibour,neiobj,C,SI) # SItime = SItime + @CPUelapsed                
                push!(exploredSI,SI);
        end
        push!(IGPair,[I,G])
    end
    return dvar,LB,newsol #,nbtime,SItime,FBtime
end

data = Data("/home/ak121396/Desktop/instances/MIPLIB(official)/mps/cvs16r128-89.mps","/home/ak121396/Desktop/instances/MIPLIB(official)/vlp/cvs16r128-89.vlp/")
pre = Val2("/home/ak121396/Desktop/solvers/Bensolve/MIPLIB/official/X/cvs16r128-89_pre_img_p.sol","/home/ak121396/Desktop/solvers/Bensolve/MIPLIB/official/Y/cvs16r128-89_img_p.sol", data.C)
# pre = valu("/home/ak121396/Desktop/solvers/Bensolve/MIPLIB/common/sparsX/cvs16r89-60_preX.jld2","/home/ak121396/Desktop/solvers/Bensolve/MIPLIB/common/Y/cvs16r89-60_img_p.sol")
Bentime = readdlm("/home/ak121396/Desktop/solvers/Bensolve/MIPLIB/common/bentime/cvs16r89-60.log")[1];
mip = JuMP.Model(with_optimizer(CPLEX.Optimizer));
@variable(mip, x[1:data.n], Bin)#data.n
MOI.set(mip, MOI.RawParameter("CPX_PARAM_SCRIND"), false );
MOI.set(mip, MOI.RawParameter("CPX_PARAM_THREADS"),1  );

for k=1:data.m
    if data.signs[k] == "l"
        @constraint(mip, dot(data.B[k,:],x) >= data.RHS[k])
    elseif data.signs[k] == "u"
        @constraint(mip, dot(data.B[k,:],x) <= data.RHS[k])
    else
        @constraint(mip, dot(data.B[k,:],x) == data.RHS[k])
    end
end
optimize!(mip);

runtime1 = @CPUelapsed candset,candobj,newsol = GPR(data.C,data.n,pre.dvar,pre.LB,30-Bentime);
runtime2 = @CPUelapsed finalX,finalY = Postpro(candset,candobj,newsol)
# push!(finalY,[0,0,0])
