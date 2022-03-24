using DataStructures,DataFrames,DelimitedFiles,JuMP,CPLEX,LinearAlgebra,StatsBase,MathProgBase,MathOptInterface
# using JLD2,CSV,CPUTime
const MPB = MathProgBase;

mutable struct CallModel
    lpfile::String; m::Int; n::Int; C::Array{}; B::Array{}; RHS::Dict{}; signs::Array{}; vub::Array{}
    function CallModel(lpfile::String)
        lpmodel=buildlp([-1,0],[2 1],'<',1.5, CplexSolver(CPX_PARAM_SCRIND=0))
        # lpmodel = CPLEX.CplexMathProgModel();
        MPB.loadproblem!(lpmodel,lpfile)
        Bmtx = MPB.getconstrmatrix(lpmodel);
        B = Bmtx[3:end,:]; C = Bmtx[1:2,:]
        # cut = findall(i-> varub[i]==1 &&varub[i+1]!=1, 1:length(varub))[end]
        # vub = varub[1:cut]; B = Bmtx[3:end,1:cut]; C = Bmtx[1:2,1:cut]
        m,n=size(B)
        vub = MPB.getvarUB(lpmodel)
        lb = MPB.getconstrLB(lpmodel)[3:end]
        ub = MPB.getconstrUB(lpmodel)[3:end]
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
        new(lpfile,m,n,C,B,RHS,signs,vub)
    end
end
mutable struct Valu
    x::String; y::String; dvar::Array{}; LB::Array{}; LBmtx::Array{};
    function Valu(x,y)
        # JLD2.@load x dv; dv0 = Array(dv);
        dv0 = readdlm(x)
        dv1 = round.(dv0; digits=4);
        objs = round.(readdlm(y); digits=4);
        ind = findall(i-> 0 in objs[i,:]  , 1:size(objs)[1]);
        dv2 = dv1[setdiff(1:end, ind), :];
        LBmtx = objs[setdiff(1:end, ind), 2:end];
        dvar = [Vector(dv2[i, :]) for i = 1:size(LBmtx)[1]];
        LB = [Vector(LBmtx[i, :]) for i = 1:size(LBmtx)[1]];
        new(x,y,dvar,LB,LBmtx)
    end
end
mt = CallModel("F:/scnd/test01S2.lp")
pr = Valu("F:/scnd/test01S1_pre_img_p.sol","F:/scnd/test01S2_img_p.sol")
mt = CallModel("/home/ak121396/Desktop/instances/SCND/test01S2.lp")
pr = Valu("/home/ak121396/Desktop/instances/SCND/test01S2_pre_img_p.sol","/home/ak121396/Desktop/instances/SCND/test01S2_img_p.sol")

function getobjval(x,C)
    return [dot(x,C[1,:]),dot(x,C[2,:])]
end
function flip(x_h,j,e)
    if x_h[e[j]]==1
        x_h[e[j]] = 0
    else
        x_h[e[j]] = 1
    end
    return x_h
end
function flipoper(Tabu,x_t,x_r)
    e = sortperm(abs.(x_t-x_r),rev=true)
    xi = []
    x_h = copy(x_r)
    j = 1
    M=length(x_t) #
    while j<=M && xi==[]
        x_h = flip(x_h,j,e)
        if x_h ∉ Tabu
            xi=x_h
        else
            j+=1
        end
    end

    if xi==[]
        j = 1
        while j<=M
            x_h=copy(x_r)
            Num = Int64(rand(ceil(length(x_r)/2):length(x_r)-1))
            R = sample(1:M,Num, replace=false)
            for i in R
                x_h = flip(x_h,r)
                if x_h ∉ Tabu
                    xi = x_h
                end
            end
            j+=1
        end
    end
    return xi
end

function findsol(x_r,var)
    for k in var
        JuMP.fix(x[k],x_r[k]; force=true)
    end
    optimize!(scnd)
    if termination_status(scnd) == MOI.OPTIMAL
        return JuMP.value.(x)
    else
        return []
    end
end
function fbcheck(xx,n)
    for k=1:n
        JuMP.fix(x[k],xx[k]; force=true)
    end
    optimize!(scnd)
    if termination_status(scnd) == MOI.OPTIMAL
        return true
    else
        return false
    end
end
function fbsearch(x,bvar,C,θ) #solveLP
    idx0 = findall(k->x[k]==0, bvar)
    idx1 = findall(k->x[k]==1, bvar)
    @objective( dist, Min, (1-θ)*(sum(dx[i] for i in idx0) + sum(1-(dx[j]) for j in idx1)) +
        θ*((dot(x,C[1,:])+dot(x,C[2,:]))/sqrt(norm(C[1,:])+norm(C[2,:]))) )
    optimize!(dist)
    if termination_status(dist) == MOI.OPTIMAL
        return JuMP.value.(dx)
    else
        return false;
    end
end

function fbsearch(x_r,bvar,C) #solveLP
    idx0 = findall(k->x[k]==0, bvar)
    idx1 = findall(k->x[k]==1, bvar)
    @objective( dist, Min, sum(dx[i] for i in idx0) + sum(1-(dx[j]) for j in idx1) )
    optimize!(dist)
    if termination_status(dist) == MOI.OPTIMAL
        return JuMP.value.(dx)
    else
        return false;
    end
end
function dominated(y,P)
    bvar = false
    for k=1:length(P)
        if all( y .>= P[k])# && any(x > P[k])
            bvar=true; break
        else
            continue
        end
    end

    return bvar
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
function Postpro(candX,candY,newsol)
    #Filter fractional solutions from LB
    initdv = candX[1:end-newsol]
    Y = [Vector(candY[i,:]) for i=1:size(candY)[1]]
    initLB = Y[1:end-newsol]
    frac = findall(j-> trunc.(initdv[j])!= initdv[j], 1:length(initdv))
    dv2 = initdv[setdiff(1:end,frac)]; LB2 = initLB[setdiff(1:end,frac)]
    P = union(dv2,candX[end-newsol+1:end])
    Pobj = union(LB2,Y[end-newsol+1:end])

    #Filter dominated solutions
    copysol = Dict(); copyobj = Dict();
    for i=1:length(Pobj)
        copysol[i] = P[i]
        copyobj[i] = Pobj[i]
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

# dt = CallModel(ARGS[1]); pr = Valu(ARGS[2],ARGS[3])
# Bentime = readdlm(ARGS[4])[1];
#################### SCND model #########################
scnd = Model(CPLEX.Optimizer);
MOI.set(scnd, MOI.RawParameter("CPX_PARAM_SCRIND"), false )
MOI.set(scnd, MOI.RawParameter("CPX_PARAM_THREADS"),1  )
bvar = findall(i->i==1,mt.vub); rvar = findall(i->i!=1,mt.vub);
@variable(scnd, yu[i in bvar], Bin)
@variable(scnd, xh[i in rvar] >=0 )
@variable(scnd, x[1:mt.n] )

for i=1:mt.n
    if i in bvar
        @constraint(scnd, x[i]==yu[i]  )
    else
        @constraint(scnd, x[i]==xh[i] )
    end
end
for k=1:mt.m
    if mt.signs[k] == "l"
        @constraint(scnd, dot(mt.B[k,:],x) >= mt.RHS[k])
    elseif mt.signs[k] == "u"
        @constraint(scnd, dot(mt.B[k,:],x) <= mt.RHS[k])
    else
        @constraint(scnd, dot(mt.B[k,:],x) == mt.RHS[k])
    end
end

optimize!(scnd)
##################### Feasibility Search model ######################
dist = Model(CPLEX.Optimizer);
MOI.set(dist, MOI.RawParameter("CPX_PARAM_SCRIND"), false );
MOI.set(dist, MOI.RawParameter("CPX_PARAM_THREADS"),1  );

@variable(dist, 0 <= dyu[i in bvar] <= 1)
@variable(dist, dxh[i in rvar] >=0 )
@variable(dist, dx[1:mt.n] )
for i=1:mt.n
    if i in bvar
        @constraint(dist, dx[i]==dyu[i]  )
    else
        @constraint(dist, dx[i]==dxh[i] )
    end
end
for k=1:mt.m
    if mt.signs[k] == "l"
        @constraint(dist, dot(mt.B[k,:],dx) >= mt.RHS[k])
    elseif mt.signs[k] == "u"
        @constraint(dist, dot(mt.B[k,:],dx) <= mt.RHS[k])
    else
        @constraint(dist, dot(mt.B[k,:],dx) == mt.RHS[k])
    end
end
optimize!(dist);
𝚯 = [0,ℯ/(ℯ+ℯ^2),ℯ^2/(ℯ+ℯ^2)]
function FP(candX,n,C,TL,bvar)#,𝚯)
    X= []; Y = []; Tabu = []; t0=time(); newsol=0; #LPcount=0;
    candlist = copy(candX)
    # k=1
    while candlist != [] #&&  time()-t0 < TL
        k = rand(1:length(candlist)); SearchDone = false;
        x_t = candlist[k]; iter=0; Max_iter = 5#length(findall(i-> 0<i<1,x_t))
        # for θ ∈ 𝚯
        while iter<Max_iter #&& time()-t0 < TL && SearchDone == false
            x_r = x_t
            for i in bvar
                x_r[i] = round(x_t[i])
            end
            x_n = findsol(x_r,bvar)
            if x_n∈X
                deleteat!(candlist,k); break;
            end
            if x_n!=[]
                push!(X,x_n); push!(Y,getobjval(x_n,C));
                newsol+=1; deleteat!(candlist,k); SearchDone = true;
                println("Rounding worked")
            else
                x1_r = [x_r[i] for i in bvar]
                if x1_r ∈ Tabu
                    x1_t = [x_t[i] for i in bvar]
                    x1_r = flipoper(Tabu,x1_t,x1_r)
                    if x1_r==[]
                        SearchDone = true;
                    else
                        for i=1:length(bvar)
                            x_r[bvar[i]] = x1_r[i]
                        end
                        x_n = findsol(x_r,bvar)
                        if x_n∈X
                            deleteat!(candlist,k); break
                        end
                        if x_n!=[] # && dominated(getobjval(x_n,C),Y)==false)
                            push!(X,x_n); push!(Y,getobjval(x_n,C)); # candY = [candY; getobjval(x_r,C)'];
                            newsol+=1; deleteat!(candlist,k); SearchDone = true;
                            println("Flip worked")
                        end
                    end
                end
                if SearchDone == false
                    push!(Tabu,x1_r)
                    x_t = fbsearch(x_r,bvar,C)#,θ)
                    if x_t == false #when there's no new feasible lp sol
                        deleteat!(candlist,k); SearchDone = true;
                    end
                end
            end
    		iter+=1
        end
        # end
    end
    return X,Y,candlist,newsol,Tabu
end
# FP(pr.dvar,pr.LB,mt.n,mt.C,60,bvar)# compiling
# FPTL = (TL-Bentime)
FPtime = @CPUelapsed fx,fy,candlist,fn,tabu = FP(pr.dvar,mt.n,mt.C,60,bvar)#,𝚯)
fy,length(candlist)
fpx,fpy = domFilter(fx,fy)
fpy
otable = zeros(length(fpy),2)
for i=1:length(fpy)
    for j=1:2
        otable[i,j] = fpy[i][j]
    end
end

ins = ARGS[2][1:end-4];
record1 = DataFrame(file=ins[26:end], totalsol = length(fpy), t=round(FPtime; digits=2))
CSV.write("/home/k2g00/k2g3475/scnd/record.csv", record1,append=true, header=false )#, delim=','i )
CSV.write(ins*"_Y.log",DataFrame(otable, :auto),append=false, header=false, delim=' ' )

1
# function FP(candX,candY,n,C,TL,bvar,rvar)
#     X= []; Y = []; Tabu = []; t0=time(); newsol=0; #LPcount=0;
#     candlist = copy(candX)
#     # k=1
#     while candlist != [] &&  time()-t0 < TL
#         k = rand(1:length(candlist))
#         x_t = candlist[k]; SearchDone = false; iter=0; Max_iter = 10
#         while iter<Max_iter && SearchDone == false && time()-t0 < TL
#             x_r = x_t
#             for i in bvar
#                 x_r[i] = round(x_t[i])
#             end
#
#             x_n = findsol(x_r,bvar)
#
#             if (x_n !=false && dominated(getobjval(x_n,C),Y)==false)
#                 push!(X,x_n); push!(Y,getobjval(x_n,C)); # candY = [candY; getobjval(x_r,C)'];
#                 newsol+=1;
#                 deleteat!(candlist,k); SearchDone = true;
#                 println(fx,"Rounding worked")
#             else
#                 if x_n ∈ Tabu
#                     x_r = flipoper(Tabu,x_t,x_r)
#                     if x_r==[]
#                         deleteat!(candlist,k); SearchDone = true;
#                     else
#                         x_n = findsol(x_r,bvar)
#                         if (x_n != false && dominated(getobjval(x_n,C),Y)==false)
#                             push!(X,x_n); push!(Y,getobjval(x_n,C)); # candY = [candY; getobjval(x_r,C)'];
#                             newsol+=1;
#                             deleteat!(candlist,k); SearchDone = true;
#                             println("flip worked")
#                         end
#                     end
#                 end
#                 if SearchDone == false
#                     push!(Tabu,x_n)
#                     x_t = fbsearch(x_r)
#                     if x_t == false #when there's no new feasible lp sol
#                         deleteat!(candlist,k); SearchDone = true;
#                     end
#                 end
#             end
#     		iter+=1
#         end
#     end
#
#     return X,Y,candlist,newsol
# end
