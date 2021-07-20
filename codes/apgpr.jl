using DataStructures,DelimitedFiles,DataFrames,JuMP,CPLEX,LinearAlgebra,CPUTime,CSV,Statistics,Clustering,MathOptInterface,StatsBase,JLD2,SparseArrays

mutable struct Data
    input::String; n::Int; C::Array{}; B::Array{};
    function Data(input::String)
        f=readdlm(input, '\t', String, '\n')
        obj=parse(Int,f[1])
        n=parse(Int,f[2])
        P=[]
        for i=3:length(f)
            a = split(f[i], ('[',']',','))
            a = filter!(e->e ∉ [ "" ,"[","," ,"]"] ,a)
            if length(a)>2
                push!(P,a)
            end
        end
        C = ones(obj,n*n); C = round.(Int,C);
        for x=1:length(P)
            for y=1:n
                p=parse(Int64,P[x][y])
                idx = Int(floor((x-1)/n))
                ind = ((x-1)*n+y)%(n*n)
                if ind != 0
                    C[idx+1,((x-1)*n+y)%(n*n)] = p
                else
                    C[idx+1,(n*n)] = p
                end
            end
        end
        A1 =zeros(n,n*n)
        A2 =zeros(n,n*n)
        A1  = round.(Int,A1)
        A2 = round.(Int,A2)
        for i=1:n
            for j=1:n
                A1[i,((i-1)*n)+j] = 1
                    if A1[i,((i-1)*n)+j]==1
                    end
            end
        end
        for i=1:n
            for j=1:n
                A2[i,((j-1)*n)+i] = 1
                if A2[i, ((j-1)*n)+i]==1
                end
            end
        end
        B=[A1;A2]
        new(input,n,C,B)
    end
end
mutable struct Val
    x::String; y::String; dvar::Array{}; LB::Array{}; LBmtx::Array{};
    function Val(x,y)
        dv = round.(digits=4, readdlm(x))
        objs = round.(digits=4, readdlm(y))
        ind = findall(i-> 0 in objs[i,:]  , 1:size(objs)[1])
        dv2 = dv[setdiff(1:end, ind), :];
        LBmtx = objs[setdiff(1:end, ind), 2:end];
        dvar = [Vector(dv2[i, :]) for i = 1:size(LBmtx)[1]]
        LB = [Vector(LBmtx[i, :]) for i = 1:size(LBmtx)[1]]
        new(x,y,dvar,LB,LBmtx)
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
    optimize!(ap_m)
    # print("status: ", termination_status(flp), "\n" )
    if termination_status(ap_m) == MOI.OPTIMAL
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
    for i=1:n*10
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

dt = Data("/home/k2g00/k2g3475/multiobjective/instances/AP/gpr/dat/AP_p-3_n-05_ins-01.dat")
pr = Val("/home/k2g00/k2g3475/multiobjective/instances/AP/gpr/X/AP_p-3_n-5_ins-01_X.sol","/home/ak121396//multiobjective/instances/KP/gpr/Y/AP_p-3_n-5_ins-01_Y.sol")
# dt = Data("/home/ak121396/Desktop/instances/AP/dat/AP_p-3_n-35_ins-1.dat")
# pr = Val("/home/ak121396/Desktop/solvers/Bensolve/APoutputs/X/AP_p-3_n-35_ins-01_X.sol","/home/ak121396//multiobjective/instances/AP/gpr/Y/AP_p-3_n-35_ins-01_Y.sol")
ap_m = Model(CPLEX.Optimizer);
MOI.set(ap_m, MOI.RawParameter("CPX_PARAM_SCRIND"), false );
MOI.set(ap_m, MOI.RawParameter("CPX_PARAM_THREADS"),1  );
@variable( ap_m, x[1:dt.n*dt.n], Bin );
for i=1:dt.n*2
    @constraint( ap_m, dot(dt.B[i,:],x) == 1 );
end
optimize!(ap_m);
GPR(dt.C,dt.n,pr.dvar,pr.LB,10);


data = Data(ARGS[1]); pre = Val(ARGS[2],ARGS[3]);
ap_m = Model(CPLEX.Optimizer);
MOI.set(ap_m, MOI.RawParameter("CPX_PARAM_SCRIND"), false );
MOI.set(ap_m, MOI.RawParameter("CPX_PARAM_THREADS"),1  );
@variable( ap_m, x[1:data.n*data.n], Bin );
for i=1:data.n*2
    @constraint( ap_m, dot(data.B[i,:],x) == 1 );
end
optimize!(ap_m);

if data.n<=30
    TL = 40;
elseif data.n<=40
    TL = 160;
else
    TL = 460;
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
class = ins[1:end-6]

record1 = DataFrame(totalsol=length(finalY), createnb = nbtime, nextSI = SItime, feastime = FBtime)
CSV.write("/home/k2g00/k2g3475/multiobjective/solvers/generalPR/goutputs/records/"*"$class"*".csv",record1,append=true, header=false )#, delim=',' )
io = open("/home/k2g00/k2g3475/multiobjective/solvers/generalPR/goutputs/time/"*"$ins"*".txt", "a")
println(io,"$totaltime"); close(io) # header=false )#, delim=',' )
CSV.write("/home/k2g00/k2g3475/multiobjective/solvers/generalPR/goutputs/Y/"*"$ins"*"_Y.log",DataFrame(otable, :auto), header=false, delim=' ' )
sparX = sparse(matriX); JLD2.@save "/home/k2g00/k2g3475/multiobjective/solvers/generalPR/goutputs/X/"*"$ins"*"_X.jdl2" sparX
