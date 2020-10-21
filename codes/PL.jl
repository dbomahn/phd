using DataStructures,DelimitedFiles,DataFrames,CSV,StatsBase,Random,LinearAlgebra,CPUTime
#Statistics,
# struct path{S<:String}
#     dir1::S; dir2::S;
# end
mutable struct Data
    dvar::String; dtfile::String
    n::Int; C::Array{}
    ub::Int; weight::Array{}
    P::Array{}
    function Data(dvar::String,dtfile::String)
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
        P = readdlm(dvar, Int)
        new(dvar,dtfile,n,C,ub,weight,P)
    end
end
function KPfbcheck(txi,n,weight,ub)
    result = false
    if sum(weight[i]*txi[i] for i=1:n) <= ub
        result = true
    else
        result = false
    end
    return result
end

function domFilter(sol,obj)
    copysol = Dict(); copyobj = Dict();
    for i=1:length(obj)
        copysol[i] = sol[i,:]
        copyobj[i] = obj[i]
    end

    for i=1:length(obj)-1
        for j=i+1:length(obj)
            if all(obj[i] .<= obj[j]) == true #dominated by PF[j]
                copyobj[i]=nothing; copysol[i]=nothing;break
            elseif all(obj[j] .<= obj[i]) == true
                copyobj[j]=nothing; copysol[j]=nothing;
            end
        end
    end

    finalsol = filter!(a->a!=nothing, collect(values(copysol)))
    finalobj = filter!(a->a!=nothing, collect(values(copyobj)))

    return finalsol,finalobj
end
function dominated(x,P)
    st = false
    for k=1:length(P)
        if all( x .<= P[k]) && any(x < P[k])
            st=true; break
        else
            continue
        end
    end
    return st
end

function PathRelinking(P,n,C,weight,ub)
    iter=1; candset = []; candobj=[]; IGPair=[]; cpP = [vec(P[i,:]) for i=1:size(P)[1]]
    for i=1:round(Int,size(P)[1]*50)
        I,G = sample(1:length(cpP), 2, replace=false)
        SI = cpP[I]; SG = cpP[G];
        while all.(SI != SG) && [I,G]∉IGPair
            rg = range(1, length=n)
            disim = findall(i-> SI[i]!=SG[i], rg); dl = length(disim)
            numnei = dl#minimum([dl,Int(n*0.2)])
            neibour = zeros(Int,numnei,n)
            cpSI = copy(SI)
            Ptable = []
            for i=1:numnei
                if cpSI[disim[i]] == 1
                    cpSI[disim[i]] = 0
                else
                    cpSI[disim[i]] = 1
                end
                neibour[i,:] = cpSI;
                push!(Ptable,[dot(C[k,:],neibour[i,:]) for k=1:3])
            end

            if rand()<0.7
                ftneibour,ftPtable = domFilter(neibour,Ptable)
                if length(ftPtable) == 1
                    SI = ftneibour[1]
                else
                    SI = rand(ftneibour)
                end
            else
                SI = rand([neibour[i,:] for i=1:numnei])
            end

            # Feasibilitycheck
            if KPfbcheck(SI,n,weight,ub)==true && SI∉cpP
                push!(cpP,SI)
                push!(candset,SI)
                push!(candobj,[dot(C[k,:],SI) for k=1:3])
                # candset = union(candset,[SI])
                # candobj = union(candobj,[[dot(C[k,:],SI) for k=1:3]])
            end
            iter+=1
        end
        push!(IGPair,[I,G])
    end
    return candset,candobj,iter
end

function improvedND(C,SI,ftPtable)
    SIobj = [dot(C[k,:],SI) for k=1:3]
    ratiotb = zeros(length(ftPtable),3)
    for i=1:length(ftPtable)
        ratiotb[i,:] = ftPtable[i]./SIobj
    end
    ranktb = zeros(length(ftPtable),3)
    for i=1:3
        ranktb[:,i] = tiedrank(ratiotb[:,i])
    end
    ranksum = [sum(ranktb[i,:]) for i=1:length(ftPtable)]
    hval = maximum(ranksum)
    ndid = findall(x-> x == hval, ranksum)
    return rand(ndid)
end

function similarity(SI,cpP,n)
    sim = Dict()
    for k=1:length(cpP)
        if SI ≠ cpP[k]
            ct = count(SI[i] == cpP[k][i] for i=1:n)
            sim[k] =ct
        end
    end
    return sim
end

function PostProc(P,C,cand,candobj,ub,n,weight)
    Pobj = []
    for i=1:size(P)[1]
        push!(Pobj,[dot(C[k,:],P[i,:]) for k=1:3])
    end
    newsol = 0
    for i=1:length(cand)
        if dominated(candobj[i],Pobj)==false && candobj[i]∉Pobj
            P = [P; transpose(cand[i][:])]
            push!(Pobj,candobj[i][:])
            newsol+=1
        end
    end
    return P,Pobj,newsol
end
function mostdiffPR(P,n,C,weight,ub)
    iter=1; candset = []; candobj=[]; IGPair=[]; cpP = [vec(P[i,:]) for i=1:size(P)[1]]
    for i=1:round(Int,size(P)[1]*50)
        I = rand(1:length(cpP))
        SI = cpP[I];
        simil = similarity(SI,cpP,n)
        lsim = minimum(collect(values(simil)))
        mostdiff = findall(x-> simil[x] == lsim, collect(keys(simil)))
        for i = 1:length(mostdiff)
            global G = rand(mostdiff)
            if [I,G] ∉ IGPair
                global SG = cpP[G]
            end
        end
        while all.(SI != SG) && [I,G]∉IGPair
            rg = range(1, length=n)
            disim = findall(i-> SI[i]!=SG[i], rg); dl = length(disim)
            numnei = dl#minimum([dl,Int(n*0.2)])
            neibour = zeros(Int,numnei,n)
            cpSI = copy(SI)
            Ptable = []
            for i=1:numnei
                if cpSI[disim[i]] == 1
                    cpSI[disim[i]] = 0
                else
                    cpSI[disim[i]] = 1
                end
                neibour[i,:] = cpSI;
                push!(Ptable,[dot(C[k,:],neibour[i,:]) for k=1:3])
            end

            if rand()<0.7
                ftneibour,ftPtable = domFilter(neibour,Ptable)
                if length(ftPtable) == 1
                    SI = ftneibour[1]
                else
                    SI = rand(ftneibour)
                end
            else
                SI = rand([neibour[i,:] for i=1:numnei])
            end

            # Feasibilitycheck
            if KPfbcheck(SI,n,weight,ub)==true && SI∉cpP
                push!(cpP,SI)
                push!(candset,SI)
                push!(candobj,[dot(C[k,:],SI) for k=1:3])
                # candset = union(candset,[SI])
                # candobj = union(candobj,[[dot(C[k,:],SI) for k=1:3]])
            end
            iter+=1
        end
        push!(IGPair,[I,G])
    end
    return candset,candobj,iter
end
function mostsimPR(P,n,C,weight,ub)
    iter=1; candset = []; candobj=[]; IGPair=[]; cpP = [vec(P[i,:]) for i=1:size(P)[1]]
    for i=1:round(Int,size(P)[1]*50)
        I = rand(1:length(cpP))
        SI = cpP[I];
        simil = similarity(SI,cpP,n)
        hsim = maximum(collect(values(simil)))
        mostsim = findall(x-> simil[x] == hsim, collect(keys(simil)))
        for i = 1:length(mostsim)
            global G = rand(mostsim)
            if [I,G] ∉ IGPair
                global SG = cpP[G]
            end
        end
        while all.(SI != SG) && [I,G]∉IGPair
            rg = range(1, length=n)
            disim = findall(i-> SI[i]!=SG[i], rg); dl = length(disim)
            numnei = dl#minimum([dl,Int(n*0.2)])
            neibour = zeros(Int,numnei,n)
            cpSI = copy(SI)
            Ptable = []
            for i=1:numnei
                if cpSI[disim[i]] == 1
                    cpSI[disim[i]] = 0
                else
                    cpSI[disim[i]] = 1
                end
                neibour[i,:] = cpSI;
                push!(Ptable,[dot(C[k,:],neibour[i,:]) for k=1:3])
            end

            if rand()<0.7
                ftneibour,ftPtable = domFilter(neibour,Ptable)
                if length(ftPtable) == 1
                    SI = ftneibour[1]
                else
                    SI = rand(ftneibour)
                end
            else
                SI = rand([neibour[i,:] for i=1:numnei])
            end

            # Feasibilitycheck
            if KPfbcheck(SI,n,weight,ub)==true && SI∉cpP
                push!(cpP,SI)
                push!(candset,SI)
                push!(candobj,[dot(C[k,:],SI) for k=1:3])
                # candset = union(candset,[SI])
                # candobj = union(candobj,[[dot(C[k,:],SI) for k=1:3]])
            end
            iter+=1
        end
        push!(IGPair,[I,G])
    end
    return candset,candobj,iter
end
function ratioPR(P,n,C,weight,ub)
    iter=1; candset = []; candobj=[]; IGPair=[]; cpP = [vec(P[i,:]) for i=1:size(P)[1]]
    for i=1:round(Int,size(P)[1]*50)
        I,G = sample(1:length(cpP), 2, replace=false)
        SI = cpP[I]; SG = cpP[G];
        while all.(SI != SG) && [I,G]∉IGPair
            rg = range(1, length=n)
            disim = findall(i-> SI[i]!=SG[i], rg); dl = length(disim)
            numnei = dl#minimum([dl,Int(n*0.2)])
            neibour = zeros(Int,numnei,n)
            cpSI = copy(SI)
            Ptable = []
            for i=1:numnei
                if cpSI[disim[i]] == 1
                    cpSI[disim[i]] = 0
                else
                    cpSI[disim[i]] = 1
                end
                neibour[i,:] = cpSI;
                push!(Ptable,[dot(C[k,:],neibour[i,:]) for k=1:3])
            end

            if rand()<0.7
                ftneibour,ftPtable = domFilter(neibour,Ptable)
                if length(ftPtable) == 1
                    SI = ftneibour[1]
                else
                    id = improvedND(C,SI,ftPtable)
                    SI = ftneibour[id]
                end
            else
                SI = rand([neibour[i,:] for i=1:numnei])
            end

            # Feasibilitycheck
            if KPfbcheck(SI,n,weight,ub)==true && SI∉cpP
                push!(cpP,SI)
                push!(candset,SI)
                push!(candobj,[dot(C[k,:],SI) for k=1:3])
                # candset = union(candset,[SI])
                # candobj = union(candobj,[[dot(C[k,:],SI) for k=1:3]])
            end
            iter+=1
        end
        push!(IGPair,[I,G])
    end
    return candset,candobj,iter
end
function diffratioPR(P,n,C,weight,ub)
    iter=1; candset = []; candobj=[]; IGPair=[]; cpP = [vec(P[i,:]) for i=1:size(P)[1]]
    for i=1:round(Int,size(P)[1]*50)
        I = rand(1:length(cpP))
        SI = cpP[I];
        simil = similarity(SI,cpP,n)
        lsim = minimum(collect(values(simil)))
        mostdiff = findall(x-> simil[x] == lsim, collect(keys(simil)))
        for i = 1:length(mostdiff)
            global G = rand(mostdiff)
            if [I,G] ∉ IGPair
                global SG = cpP[G]
            end
        end
        while all.(SI != SG) && [I,G]∉IGPair
            rg = range(1, length=n)
            disim = findall(i-> SI[i]!=SG[i], rg); dl = length(disim)
            numnei = dl#minimum([dl,Int(n*0.2)])
            neibour = zeros(Int,numnei,n)
            cpSI = copy(SI)
            Ptable = []
            for i=1:numnei
                if cpSI[disim[i]] == 1
                    cpSI[disim[i]] = 0
                else
                    cpSI[disim[i]] = 1
                end
                neibour[i,:] = cpSI;
                push!(Ptable,[dot(C[k,:],neibour[i,:]) for k=1:3])
            end

            if rand()<0.7
                ftneibour,ftPtable = domFilter(neibour,Ptable)
                if length(ftPtable) == 1
                    SI = ftneibour[1]
                else
                    id = improvedND(C,SI,ftPtable)
                    SI = ftneibour[id]
                end
            else
                SI = rand([neibour[i,:] for i=1:numnei])
            end

            # Feasibilitycheck
            if KPfbcheck(SI,n,weight,ub)==true && SI∉cpP
                push!(cpP,SI)
                push!(candset,SI)
                push!(candobj,[dot(C[k,:],SI) for k=1:3])
                # candset = union(candset,[SI])
                # candobj = union(candobj,[[dot(C[k,:],SI) for k=1:3]])
            end
            iter+=1
        end
        push!(IGPair,[I,G])
    end
    return candset,candobj,iter
end

function simratioPR(P,n,C,weight,ub)
    iter=1; candset = []; candobj=[]; IGPair=[]; cpP = [vec(P[i,:]) for i=1:size(P)[1]]
    for i=1:round(Int,size(P)[1]*50)
        I = rand(1:length(cpP))
        SI = cpP[I];
        simil = similarity(SI,cpP,n)
        hsim = maximum(collect(values(simil)))
        mostsim = findall(x-> simil[x] == hsim, collect(keys(simil)))
        for i = 1:length(mostsim)
            global G = rand(mostsim)
            if [I,G] ∉ IGPair
                global SG = cpP[G]
            end
        end
        while all.(SI != SG) && [I,G]∉IGPair
            rg = range(1, length=n)
            disim = findall(i-> SI[i]!=SG[i], rg); dl = length(disim)
            numnei = dl#minimum([dl,Int(n*0.2)])
            neibour = zeros(Int,numnei,n)
            cpSI = copy(SI)
            Ptable = []
            for i=1:numnei
                if cpSI[disim[i]] == 1
                    cpSI[disim[i]] = 0
                else
                    cpSI[disim[i]] = 1
                end
                neibour[i,:] = cpSI;
                push!(Ptable,[dot(C[k,:],neibour[i,:]) for k=1:3])
            end

            if rand()<0.7
                ftneibour,ftPtable = domFilter(neibour,Ptable)
                if length(ftPtable) == 1
                    SI = ftneibour[1]
                else
                    id = improvedND(C,SI,ftPtable)
                    SI = ftneibour[id]
                end
            else
                SI = rand([neibour[i,:] for i=1:numnei])
            end

            # Feasibilitycheck
            if KPfbcheck(SI,n,weight,ub)==true && SI∉cpP
                push!(cpP,SI)
                push!(candset,SI)
                push!(candobj,[dot(C[k,:],SI) for k=1:3])
                # candset = union(candset,[SI])
                # candobj = union(candobj,[[dot(C[k,:],SI) for k=1:3]])
            end
            iter+=1
        end
        push!(IGPair,[I,G])
    end
    return candset,candobj,iter
end

# paths = ("E:\\Bensolve_KP\\RoundDown\\X\\", "E:\\Bensolve_KP\\data\\")
paths = ("/home/ak121396/Desktop/Bensolve_KP/RoundDown/X/", "/home/ak121396/Desktop/Bensolve_KP/data/")
file = readdir(paths[1]); ins = readdir(paths[2])


# u=1
# while u<11
for i=1:100
    kp = Data(paths[1]*file[i],paths[2]*ins[i])
    runtime1 = @CPUelapsed cand,candobj,iter = ratioPR(kp.P,kp.n,kp.C,kp.weight,kp.ub)
    print("#sol Befor filtering : ", length(candobj),"\n")
    runtime2 = @CPUelapsed P,Pobj,newsol = PostProc(kp.P,kp.C,cand,candobj,kp.ub,kp.n,kp.weight)

    otable = ones(Int, length(Pobj),3)
    for i=1:length(Pobj)
        for j=1:3
            otable[i,j] = -Pobj[i][j]
        end
    end
    fname = kp.dtfile[end-21:end-3]
    # CSV.write("C:\\Users\\AK121396\\Desktop\\simY\\"*"$fname"*"Y.log",DataFrame(otable),header=false, delim=' ')
    CSV.write("/home/ak121396/Desktop/PR_KP/"*"ratio2"*"/"*"$fname"*"Y.log",DataFrame(otable),header=false, delim=' ' )
    # #########################  Record outputs  ############################
    record1 = DataFrame(newsol=newsol, sol=length(Pobj),CPUtime=runtime1+runtime2, iter=iter)
    CSV.write("/home/ak121396/Desktop/PR_KP/"*"ratio2"*"/2ndKP_record.csv",record1, append=true, header=false )#, delim=',' )

    # CSV.write("/home/ak121396/Desktop/BENKP/"*"$u"*"/diff_KP_record.csv",record1, append=true, header=false )#, delim=',' )
    # CSV.write("C:\\Users\\AK121396\\Desktop\\simY\\"*"/sim_KP_record.csv",record1, append=true, header=false )
end

#     global u=u+1
# end



# function PickNeibour()
#     ptw = kp.C./kp.weight
#
#     ranks = [ordinalrank(ptw[i,:]) for i=1:3]
#
#     findall(x->x==kp.n,rankptw)
#
#     rankptw[4]
#
#     ptw[rankptw[6]]
#     maximum(ptw)
# end
