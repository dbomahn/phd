using DataStructures,DelimitedFiles,DataFrames,JuMP,CPLEX,LinearAlgebra,CPUTime,CSV,Statistics,Clustering,MathOptInterface,StatsBase,JLD2,SparseArrays
mutable struct Data
    dtfile::String; n::Int; C::Array{}; ub::Int; weight::Array{};#P::Array{}
    function Data(dtfile::String)
        d = readdlm(dtfile)
        data = readdlm(dtfile,'\t', String, '\n')
        b = data[4:length(data)-1]
        n=parse(Int,data[2])
        C= ones(length(b),n)
        C = round.(Int,C)
        for x=1:length(b)
            a = split(b[x],('[',']',','))
            aa = filter!(e->!(e in [ "" ,"[","," ,"]"]) ,a)
            for y=1:length(aa)
                c=parse(Int64,aa[y])
                C[x,y] = c
                if c==0
                    global ct = ct+1;
                end
            end
        end
        ub=parse(Int,data[3])
        weight =ones(1,n)
        weight = round.(Int,weight)
        item = data[length(data)]
        w1 = split(item, ('[',']',','))
        w2 = filter!(e->!(e in ["" ,"[", "]"]) ,w1)
        for i=1:n
            weight[i] = parse(Int64,w2[i])
        end
        new(dtfile,n,C,ub,weight)#,P)
    end
end

mutable struct Val
    x::String; y::String; dvar::Array{}; LB::Array{}; LBmtx::Array{};
    function Val(x,y)
        dv = round.(digits=4, readdlm(x))
        objs = round.(digits=2, readdlm(y))#[:,1:end]
        ind = findall(i-> 0 in objs[i,:]  , 1:size(objs)[1])
        dv2 = dv[setdiff(1:end, ind), :];
        LBmtx = -objs[setdiff(1:end, ind), 2:end];
        dvar = [Vector(dv2[i, :]) for i = 1:size(LBmtx)[1]]
        LB = [Vector(LBmtx[i, :]) for i = 1:size(LBmtx)[1]]
        # maxobj = [maximum(LBmtx[:,i]) for i=1:3]; minobj = [minimum(LBmtx[:,i]) for i=1:3];
        # steps = round.(Int,(maxobj-minobj)/j)  #length(LB[:,1])
        # interval = round(Int,mean(maxobj-minobj)/j)
        # points = transpose(LB)
        # α = max(round(length(LB)*0.05),1)
        # radius = interval*α
        new(x,y,dvar,LB,LBmtx)#,radius)
    end
end

function getobjval(x,C)
    return [-dot(x,C[1,:]),-dot(x,C[2,:]),-dot(x,C[3,:])]
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
    optimize!(kp_m)
    # print("status: ", termination_status(flp), "\n" )
    if termination_status(kp_m) == MOI.OPTIMAL
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
            mostimp = findall(x-> x == maximum(ranksum), ranksum)
            k = rand(mostimp)[1]
            return neibour[k]#, neiobj[k]
        end
    end
end
function Postpro(dvar,LB,newsol)
    #Filter fractional solutions from LB
    initdv = dvar[1:end-newsol]
    initLB = LB[1:end-newsol]
    frac = findall(j-> trunc.(initdv[j])!= initdv[j], 1:length(initdv))
    dv2 = initdv[setdiff(1:end,frac)]; LB2 = initLB[setdiff(1:end,frac)]
    P = union(dv2,dvar[end-newsol+1:end])
    Pobj = union(LB2,LB[end-newsol+1:end])

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
function GPR(C,n,dvar,LB,TL)
    IGPair=[]; exploredSI = []; newsol=0;
    nbtime = 0; SItime = 0; FBtime = 0; t0=time();
    for i=1:n*5
    	if time()-t0 >= TL
            break
        end
        I,G = sample(1:length(dvar), 2, replace=false)
        SI = dvar[I]; SG = dvar[G]; iter=0;
        while all.(SI != SG) && [I,G]∉IGPair && iter<20 && (time()-t0<TL)
            dif = findall(i-> SI[i]!=SG[i], 1:n)
            nbtime = nbtime + @CPUelapsed neibour,neiobj = createNB(SI,C,dif,exploredSI)
            if length(neiobj)==0
                break
            else
                for l=1:length(neiobj)
                    FBtime = FBtime + @CPUelapsed feasi = FBcheck(neibour[l],n)
                    if feasi == true && neibour[l]∉ dvar
                        push!(dvar, neibour[l]); push!(LB, neiobj[l]);
                        newsol+=1
                    end
                end
                SItime = SItime + @CPUelapsed SI = nextSI(neibour,neiobj,C,SI)
                if SI∉dvar
                    push!(exploredSI,SI);
                end
            end
            iter+=1
        end
        push!(IGPair,[I,G])
    end
    return dvar,LB,nbtime,SItime,FBtime,newsol
end

dt = Data("/home/ak121396//multiobjective/instances/KP/gpr/dat/KP_p-3_n-030_ins-1.dat")
pr = Val("/home/ak121396//multiobjective/instances/KP/gpr/X/KP_p-3_n-030_ins-1.x.sol","/home/ak121396//multiobjective/instances/KP/gpr/Y/KP_p-3_n-030_ins-1.y.sol")
kp_m = Model(CPLEX.Optimizer);
MOI.set(kp_m, MOI.RawParameter("CPX_PARAM_SCRIND"), false );
MOI.set(kp_m, MOI.RawParameter("CPX_PARAM_THREADS"),1  );
@variable(kp_m, x[1:dt.n], Bin );
@constraint( kp_m, -dot(dt.weight,x) >= -(dt.ub) );
optimize!(kp_m);
candset,candobj,nbtime,SItime,FBtime,newsol = GPR(dt.C,dt.n,pr.dvar,pr.LB,10);
finalX,finalY = Postpro(candset,candobj,newsol)

data = Data("/home/k2g00/k2g3475/multiobjective/instances/KP/gpr/dat/KP_p-3_n-10_ins-1.dat")
pre = Val("/home/k2g00/k2g3475/multiobjective/instances/KP/gpr/X/KP_p-3_n-10_ins-1.x.sol","/home/k2g00/k2g3475/multiobjective/instances/KP/gpr/Y/KP_p-3_n-10_ins-1.y.sol")
kp_m = Model(CPLEX.Optimizer);
MOI.set(kp_m, MOI.RawParameter("CPX_PARAM_SCRIND"), false );
MOI.set(kp_m, MOI.RawParameter("CPX_PARAM_THREADS"),1  );
@variable(kp_m, x[1:data.n], Bin );
@constraint( kp_m, -dot(data.weight,x) >= -(data.ub) );
optimize!(kp_m);

if data.n<=60
    TL = 10;
elseif data.n<=80
    TL = 30;
else
    TL = 160;
end
runtime1 = @CPUelapsed candset,candobj,nbtime,SItime,FBtime,newsol = GPR(data.C,data.n,pre.dvar,pre.LB,TL)
runtime2 = @CPUelapsed finalX,finalY = Postpro(candset,candobj,newsol)
push!(finalY,[0,0,0])
totaltime = runtime1+runtime2
otable = zeros(Int, length(finalY),3)
for i=1:length(finalY)
    for j=1:3
        otable[i,j] = finalY[i][j]
    end
end
matriX = zeros(Int,length(finalX),data.n)
for i=1:length(finalX)
    for j=1:data.n
        matriX[i,j] = finalX[i][j]
    end
end
@show ins = ARGS[1][end-15:end-4]
class = ins[1:end-5]

record1 = DataFrame(totalsol=length(finalY), numnewsol = newsol, createnb = nbtime, nextSI = SItime, feastime = FBtime)
CSV.write("/home/k2g00/k2g3475/multiobjective/solvers/generalPR/goutputs/records/"*"$class"*".csv",record1,append=true, header=false )#, delim=',' )
io = open("/home/k2g00/k2g3475/multiobjective/solvers/generalPR/goutputs/time/"*"$ins"*".txt", "a")
println(io,"$totaltime"); close(io) # header=false )#, delim=',' )
CSV.write("/home/k2g00/k2g3475/multiobjective/solvers/generalPR/goutputs/Y/"*"$ins"*"_Y.log",DataFrame(otable, :auto), header=false, delim=' ' )
sparX = sparse(matriX); JLD2.@save "/home/k2g00/k2g3475/multiobjective/solvers/generalPR/goutputs/X/"*"$ins"*"_X.jdl2" sparX
