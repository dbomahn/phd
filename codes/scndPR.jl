using DataStructures,DataFrames,DelimitedFiles,JuMP,JLD2,CPLEX,LinearAlgebra,StatsBase,CPUTime,MathProgBase,MathOptInterface
# using CSV
const MPB = MathProgBase;

mutable struct CallModel
    lpfile::String; m::Int; n::Int; C::Array{}; B::Array{}; RHS::Dict{}; signs::Array{}; vub::Array{}
    function CallModel(lpfile::String)
        lpmodel=buildlp([-1,0],[2 1],'<',1.5, CplexSolver(CPX_PARAM_SCRIND=0))
        # lpmodel = CPLEX.CplexMathProgModel();
        MPB.loadproblem!(lpmodel,lpfile)
        Bmtx = MPB.getconstrmatrix(lpmodel);
        B = Bmtx[3:end,:]; C = Bmtx[1:2,:]
        # cut = find(i-> varub[i]==1 &&varub[i+1]!=1, 1:length(varub))[end]
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
        JLD2.@load x dv;
        dv0 = Array(dv);
        dv1 = round.(dv0; digits=4);
        objs = round.(digits=4, readdlm(y));
        ind = findall(i-> 0 in objs[i,:]  , 1:size(objs)[1]);
        dv2 = dv1[setdiff(1:end, ind), :];
        LBmtx = objs[setdiff(1:end, ind), 2:end];
        dvar = [Vector(dv2[i, :]) for i = 1:size(LBmtx)[1]];
        LB = [Vector(LBmtx[i, :]) for i = 1:size(LBmtx)[1]];
        new(x,y,dvar,LB,LBmtx)
    end
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

function PRfindsol(x_r,bvar)
    for (i,v) in enumerate(bvar)
        JuMP.fix(x[v], x_r[i]; force=true)
    end
    optimize!(scnd2)
    if termination_status(scnd2) == MOI.OPTIMAL
        return JuMP.value.(x)
    else
        return []
    end
end

function dominated(x,P)
    st = false
    for k=1:length(P)
        if all( x .>= P[k])# && any(x > P[k])
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
function createNB(SI,C,dif,exploredSI,bvar)
    neibour = []; neiobj = [];
    for i in dif
        cpSI = copy(SI)
        if cpSI[i] == round(cpSI[i]) #if vari is int value
            if cpSI[i] == 1
                cpSI[i] = 0
            else
                cpSI[i] = 1
            end
            push!(neibour, cpSI); push!(neiobj, binobj(cpSI,bvar))
        else #if vari is fractional
            cpSI[i] = 1; push!(neibour, cpSI); push!( neiobj, binobj(cpSI,bvar) )
            cpSI = copy(SI)
            cpSI[i] = 0; push!(neibour, cpSI); push!( neiobj, binobj(cpSI,bvar) )
        end
    end
    idx = findall(i-> neibour[i] in exploredSI, 1:length(neibour))
    neibour = setdiff(neibour,exploredSI)
    neiobj = neiobj[setdiff(1:end, idx),:]

    return neibour,neiobj
end

function nextSI(neibour,neiobj,C,SI,bvar)
    SIobj = binobj(SI,bvar)#,C,bvar)
    # neiobj = [getobjval(neibour[i],C) for i=1:length(neibour)]
    for i=1:length(neiobj)
        if length(neiobj) == 1  #if there is one candiate sol
            return neibour[1]#,neiobj[1]
        else length(neiobj) > 1 # if there are multiple candiates, check the improved ratio
            ratiotb = zeros(length(neiobj),length(neiobj[1]))
            for i=1:length(neiobj)
                ratiotb[i,:] = neiobj[i]./SIobj
            end
            ranktb = zeros(length(neiobj),length(neiobj[1]))
            for i=1:length(neiobj[1])
                ranktb[:,i] = tiedrank(ratiotb[:,i])
            end
            ranksum = [sum(ranktb[i,:]) for i=1:length(neiobj)]
            mostimp = findall(x-> x == maximum(ranksum), ranksum)
            k = rand(mostimp)[1]
            return neibour[k]#,neiobj[k]
        end
    end
end
function Postpro(candX,candY,newsol)
    #Filter fractional solutions from LB
    initdv = candX[1:end-newsol]
    initLB = candY[1:end-newsol]
    frac = findall(j-> trunc.(initdv[j])!= initdv[j], 1:length(initdv))
    dv2 = initdv[setdiff(1:end,frac)]; LB2 = initLB[setdiff(1:end,frac)]
    P = union(dv2,candX[end-newsol+1:end])
    Pobj = union(LB2,candY[end-newsol+1:end])

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


mt = CallModel("/home/ak121396/Desktop/relise/newlp.lp");
pr = Valu("/home/ak121396/Desktop/relise/test01S2_X.jld2","/home/ak121396/Desktop/relise/test01S2_img_p.sol");
# mt = Data(ARGS[1]); pr = Valu(ARGS[2],ARGS[3])
# Bentime = readdlm(ARGS[4])[1];
#################### SCND model #########################
scnd2 = Model(CPLEX.Optimizer);
MOI.set(scnd2, MOI.RawParameter("CPX_PARAM_SCRIND"), false );
# MOI.set(scnd2, MOI.RawParameter("CPX_PARAM_THREADS"),1  )
bvar = findall(i->i==1,mt.vub); # rvar = findall(i->i!=1,mt.vub);
@variable(scnd2, x[1:mt.n]>=0 );
for i=1:mt.n
    if i in bvar
        set_binary(x[i])
    # else
    #     set_lower_bound(x[i],0);
    end
end
for k=1:mt.m
    if mt.signs[k] == "l"
        @constraint(scnd2, dot(mt.B[k,:],x) >= mt.RHS[k])
    elseif mt.signs[k] == "u"
        @constraint(scnd2, dot(mt.B[k,:],x) <= mt.RHS[k])
    else
        @constraint(scnd2, dot(mt.B[k,:],x) == mt.RHS[k])
    end
end
optimize!(scnd2)

function FBcheck(xx)
    JuMP.fix.(scnd2[:x],xx; force=true)
    optimize!(scnd2)
    if termination_status(scnd2) == MOI.OPTIMAL
        return true
    else
        return false
    end
end
function getobjval(x)
    return [ dot(x,mt.C[1,:]), dot(x,mt.C[2,:]) ]
end
function binobj(x,bvar)
    return [ dot(x,[mt.C[1,:][k] for k in bvar]), dot(x,[mt.C[2,:][k] for k in bvar]) ]
end
function GPR(lpX,lpY,TL,bvar,C)
    candX = copy(lpX); candY = copy(lpY);
    IGPair=[]; exploredSI = []; t0=time();newsol=0;
    while time()-t0 < TL
	    I,G = StatsBase.sample(1:length(candX), 2, replace=false)
        SI = candX[I][bvar[1]:bvar[end]]; SG = candX[G][bvar[1]:bvar[end]]; iter=0;
        SI_r = round.(SI); SG_r = round.(SG)
        dif = findall(i-> SI_r[i]!=SG_r[i], 1:length(bvar))
        # println("dif is: ", length(dif))
        Max_iter = length(findall(i-> 0<i<1,candX[I]))
        while length(dif)>0 && [I,G]∉IGPair && iter<Max_iter && (time()-t0<TL)
            neibour,neiobj = createNB(SI,C,dif,exploredSI,bvar)
            if (length(neiobj)==0) #(time()-t0 >= TL)
                break
            else
                for l=1:length(neibour)
                    xn = PRfindsol(round.(neibour[l]),bvar)
                    if ( xn!=[] && xn∉candX )
                        push!(candX, xn); push!(candY, getobjval(xn));
                        newsol+=1;println("new sol");
                    end
                end
            end
            SI = nextSI(neibour,neiobj,C,SI,bvar)
            if SI∉candX
                push!(exploredSI,SI);
            end
            @show iter+=1
        end
        push!(IGPair,[I,G])
    end
    return candX,candY,newsol
end
GPRtime = @CPUelapsed px,py,pn = GPR(pr.dvar,pr.LB,60,bvar,mt.C)
totaltime = FPtime+Bentime+GPRtime
prx,pry = Postpro(px,py,pn)

otable = zeros(length(pry),2)
for i=1:length(pry)
    for j=1:2
        otable[i,j] = pry[i][j]
    end
end

ins = ARGS[2][1:end-4];
record1 = DataFrame(file=ins[26:end], totalsol = length(fpry), t=round(totaltime; digits=2))
CSV.write("/home/k2g00/k2g3475/scnd/fpr_record.csv", record1,append=true, header=false )#, delim=','i )
CSV.write(ins*"_Y.log",DataFrame(otable, :auto),append=false, header=false, delim=' ' )


# function GPR(candX,candY,C,n,TL,bvar)
#     IGPair=[]; exploredSI = []; t0=time();newsol=0;
#     while time()-t0 < TL
# 	    I,G = sample(1:length(candX), 2, replace=false)
#         SI = candX[I]; SG = candX[G]; iter=0; #Maxiter = length(candX)*10;  #if the algorithm finds infeasible solutions in a row, move to the next iteration
#         dif = findall(i-> SI[i]!=SG[i], bvar)
#         while all.(SI != SG) && [I,G]∉IGPair && iter<length(dif) && (time()-t0<TL)
#             neibour,neiobj = createNB(SI,C,dif,exploredSI)
#             # newnb: fix 0,1 and calculate frac var
#             if ( (length(neiobj)==0) || (time()-t0 >= TL) )
#                 break
#             else
#                 for l=1:length(neiobj)
#                     if ( (FBcheck(neibour[l],n)==true) && neibour[l]∉candX ) #dominance check (dominated(neiobj[l],LB)==false)
#                         push!(candX, neibour[l]); push!(candY, neiobj[l]);
#                         println("new sol"); newsol+=1;
#                     end
#                 end
#             end
#             SI = nextSI(neibour,neiobj,C,SI)
#             if SI∉candX
#                 push!(exploredSI,SI);
#             end
#             iter+=1
#         end
#         push!(IGPair,[I,G])
#     end
#     return candX,candY,newsol
# end

########################################
# function FP(candX,candY,n,C,TL,st)
#     X= []; Y = []; Tabu = []; newsol=0; t0=time(); #LPcount=0;
#     candlist = copy(candX)
#     # k=1
#     while candlist != [] &&  time()-t0 < TL
#         k = rand(1:length(candlist))
#         x_t = candlist[k]; SearchDone = false; iter=0; Max_iter = 10
#         while iter<Max_iter && SearchDone == false && time()-t0 < TL
#             x1 = x_t[1:st-1]; x1_r = round.(Int,x1)
#             x_n = findsol(x1_r,st-1)
#             if ( (x_n != false) && (x_n ∉ X) )
#                 push!(X,x_n); push!(Y,getobjval(x_n,C)); # candY = [candY; getobjval(x_r,C)'];
#                 newsol+=1; SearchDone = true;deleteat!(candlist,k)
#                 # println("Rounding worked")
#             else
#                 if x_n ∈ Tabu
#                     x1_r = flipoper(Tabu,x1,x1_r)
#                     if x1_r==[]
#                         SearchDone = true; deleteat!(candlist,k)
#                     else
#                         x_n = findsol(x1_r,st-1)
#                         if ( (x_n != false) && (x_n ∉ X) )
#                             push!(X,x_n); push!(Y,getobjval(x_n,C)); # candY = [candY; getobjval(x_r,C)'];
#                             newsol+=1; SearchDone = true;deleteat!(candlist,k)
#                             # println("flip worked")
#                         end
#                     end
#                 end
#                 if SearchDone == false
#                     push!(Tabu,x_n)
#                     x_t = fbsearch(x1_r)
#                     if x_t == false #when there's no new feasible lp sol
#                         SearchDone = true
#                         deleteat!(candlist,k)
#                     end
#                 end
#             end
#     		iter+=1
#         end
#     end
#
#     return X,Y,candlist,newsol
# endFP(pr.dvar,pr.LB,mt.n,mt.C,10,st)# compiling
# FPTL = (TL-Bentime)/2
# FPtime = @CPUelapsed fcanx,fcany,X_u,fnsol = FP(pr.dvar,pr.LB,mt.n,mt.C,FPTL,st)
