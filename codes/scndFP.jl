using DataStructures,DataFrames,DelimitedFiles,JuMP,LinearAlgebra,JLD2,CPLEX,LinearAlgebra,CSV,StatsBase,CPUTime,MathProgBase,MathOptInterface
const MPB = MathProgBase;
TL = 3600
mutable struct Data
    lpfile::String; m::Int; n::Int; C::Array{}; B::Array{}; RHS::Dict{}; signs::Array{}; vub::Array{}
    function Data(lpfile::String)
        lpmodel=buildlp([-1,0],[2 1],'<',1.5, CplexSolver(CPX_PARAM_SCRIND=0))
        # lpmodel = CPLEX.CplexMathProgModel();
        MPB.loadproblem!(lpmodel,lpfile)
        Bmtx = MPB.getconstrmatrix(lpmodel);B = Bmtx[3:end,:]
        C = Bmtx[1:2,:]
        m,n=size(B)
        vub = MPB.getvarUB(lpmodel)
        lb = MPB.getconstrLB(lpmodel)[3:end]; ub = MPB.getconstrUB(lpmodel)[3:end]
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

function findsol(x_r,n)
    for k=1:n
        JuMP.fix(x[k],x_r[k]; force=true)
    end
    optimize!(scnd)
    if termination_status(scnd) == MOI.OPTIMAL
        return JuMP.value.(x)
    else
        return false
    end
end
function fbsearch(x_r) #solveLP
    idx0 = findall(k->k==0, x_r)
    idx1 = findall(k->k==1, x_r)
    @objective( dist, Min, sum(dx[i] for i in idx0) + sum(1-(dx[j]) for j in idx1) )
    optimize!(dist)
    if termination_status(dist) == MOI.OPTIMAL
        return JuMP.value.(dx)
    else
        return false;
    end
end
function dominated(y,P)
    st = false
    for k=1:length(P)
        if all( y .>= P[k])# && any(x > P[k])
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

# dt = Data("/home/ak121396/Desktop/tests1.lp")
# pr = Valu("/home/ak121396/Desktop/relise/Test1S1_X.jld2","/home/ak121396/Desktop/relise/Test1S1_img_p.sol")
dt = Data(ARGS[1]); pr = Valu(ARGS[2],ARGS[3])
Bentime = readdlm(ARGS[4])[1];
#################### SCND model #########################
st = findall(i->i!=1,dt.vub)[1]
scnd = Model(CPLEX.Optimizer);
MOI.set(scnd, MOI.RawParameter("CPX_PARAM_SCRIND"), false )
MOI.set(scnd, MOI.RawParameter("CPX_PARAM_THREADS"),1  )
@variables(scnd, begin
    yu[1:st-1], Bin
    0<= xh[st:dt.n]
end)
@variable(scnd, x[1:dt.n] )
@constraint(scnd, [k=1:st-1], x[k] == yu[k] )
@constraint(scnd, [k=st:dt.n], x[k] == xh[k] )
for k=1:dt.m
    if dt.signs[k] == "l"
        @constraint(scnd, dot(dt.B[k,1:st-1],yu)+dot(dt.B[k,st:end],xh) >= dt.RHS[k])
    elseif dt.signs[k] == "u"
        @constraint(scnd, dot(dt.B[k,1:st-1],yu)+dot(dt.B[k,st:end],xh) <= dt.RHS[k])
    else
        @constraint(scnd, dot(dt.B[k,1:st-1],yu)+dot(dt.B[k,st:end],xh) == dt.RHS[k])
    end
end
optimize!(scnd);
##################### Feasibility Search model ######################
dist = Model(CPLEX.Optimizer);
MOI.set(dist, MOI.RawParameter("CPX_PARAM_SCRIND"), false );
MOI.set(dist, MOI.RawParameter("CPX_PARAM_THREADS"),1  );
@variables(dist, begin
    0<= dyu[1:st-1] <=1
    0<= dxh[st:dt.n]
end)
@variable(dist, dx[1:dt.n] )
@constraint(dist, [k=1:st-1], dx[k] == dyu[k] )
@constraint(dist, [k=st:dt.n], dx[k] == dxh[k] )
for k=1:dt.m
    if dt.signs[k] == "l"
        @constraint(dist, dot(dt.B[k,1:st-1],dyu)+dot(dt.B[k,st:end],dxh) >= dt.RHS[k])
    elseif dt.signs[k] == "u"
        @constraint(dist, dot(dt.B[k,1:st-1],dyu)+dot(dt.B[k,st:end],dxh) <= dt.RHS[k])
    else
        @constraint(dist, dot(dt.B[k,1:st-1],dyu)+dot(dt.B[k,st:end],dxh) == dt.RHS[k])
    end
end
optimize!(dist);


function FP(candX,candY,n,C,TL,st)
    X= []; Y = []; Tabu = []; newsol=0; t0=time(); #LPcount=0;
    candlist = copy(candX)
    # k=1
    while candlist != [] &&  time()-t0 < TL
        k = rand(1:length(candlist))
        x_t = candlist[k]; SearchDone = false; iter=0; Max_iter = 10
        while iter<Max_iter && SearchDone == false && time()-t0 < TL
            x1 = x_t[1:st-1]; x1_r = round.(Int,x1)
            x_n = findsol(x1_r,st-1)
            fx = getobjval(x_n,C)
            if ( (x_n != false) && (dominated(fx,Y)==false) )
                push!(X,x_n); push!(Y,fx); # candY = [candY; getobjval(x_r,C)'];
                newsol+=1; SearchDone = true;deleteat!(candlist,k)
                # println("Rounding worked")
            else
                if x_n ∈ Tabu
                    x1_r = flipoper(Tabu,x1,x1_r)
                    if x1_r==[]
                        SearchDone = true; deleteat!(candlist,k)
                    else
                        x_n = findsol(x1_r,st-1)
                        fx = getobjval(x_n,C)
                        if ( (x_n != false) && (dominated(fx,Y)==false) )
                            push!(X,x_n); push!(Y,fx); # candY = [candY; getobjval(x_r,C)'];
                            newsol+=1; SearchDone = true;deleteat!(candlist,k)
                            # println("flip worked")
                        end
                    end
                end
                if SearchDone == false
                    push!(Tabu,x_n)
                    x_t = fbsearch(x1_r)
                    if x_t == false #when there's no new feasible lp sol
                        SearchDone = true
                        deleteat!(candlist,k)
                    end
                end
            end
    		iter+=1
        end
    end

    return X,Y,candlist,newsol
end

FP(pr.dvar,pr.LB,dt.n,dt.C,10,st)# compiling
FPTL = (TL-Bentime)
FPtime = @CPUelapsed fcanx,fcany,candlist,newsol = FP(pr.dvar,pr.LB,dt.n,dt.C,FPTL,st)
fpx,fpy = domFilter(fcanx,fcany)
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

###########################################################
FPtime = @CPUelapsed fcanx,fcany,X_u,newsol = FP(pr.dvar,pr.LB,dt.n,dt.C,300,st)
fpx,fpy = domFilter(fcanx,fcany)

newsol
length(candlist)
cdr = copy(pr.dvar)
append!(dt.n,cdr[1,:],cdr[5,:])

function GPR(dvar,LB,C,n,TL,st)
    IGPair=[]; exploredSI = []; newsol=0; t0=time();
	X = copy(dvar); Y = copy(LB);
	while time()-t0 < TL
        I,G = sample(1:length(Y), 2, replace=false); infeasi=0; Maxiter = 100
        SI = round.(X[I][1:st-1]); SG = round.(X[G][1:st-1]);
		# SI = X[I][1:st-1]; SG = X[G][1:st-1];
        while all.(SI != SG) && [I,G]∉IGPair && infeasi<Maxiter && (time()-t0<TL)
            dif = findall(i-> SI[i]!=SG[i], 1:st-1)
            neibour,neiobj = createNB(SI,C,dif,exploredSI)
			if length(neiobj)==0
                break
            else
	            for l=1:length(neibour)
					x_n = findsol(neibour[l],st-1)
					if ( (x_n != false) && (dominated(x_n,X)==false) )
	                    push!(X,x_n); push!(Y,getobjval(x_n,C));
	                    newsol+=1;
						print("newsol added \n")
	                else
						infeasi+=1;
						print("infeasible \n")
	                end
	            end
	            SI = nextSI(neibour,neiobj,C,SI)
	            if SI∉X
	                push!(exploredSI,SI);
	            end
			end
        end
        push!(IGPair,[I,G])

    end
    return X,Y,newsol
end


GPRtime = @CPUelapsed pcanx,pcany,pnew = GPR(pr.dvar,pr.LB,dt.C,dt.n,300,st)

Y_u = [getobjval(X_u[i],dt.C) for i=1:length(X_u)]
GPRtime = @CPUelapsed pcanx,pcany,pnew = GPR([fcanx;X_u],[fcany;Y_u],dt.C,dt.n,300,st)
# totaltime = FPtime+Bentime+GPRtime

prx,pry = Postpro(pcanx,pcany,pnew,st)
