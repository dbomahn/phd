
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
        new(input,n*n,C,B)
    end
end

mutable struct Val
    x::String; y::String; n::Int; dvar::Array{}; L::Array{}; LB::Array{}; k::Int
    function Val(x,y,n)
        dv = round.(digits=4, readdlm(x))
        objs = round.(digits=4, readdlm(y))
        ind = findall(i-> 0 in objs[i,:]  , 1:size(objs)[1])
        dv2 = dv[setdiff(1:end, ind), :];
        L = objs[setdiff(1:end, ind), 2:end];
        dvar = [Vector(dv2[i, :]) for i = 1:size(L)[1]]
        LB = [Vector(L[i, :]) for i = 1:size(L)[1]]
        k = ceil(Int,size(LB)[1]/Int(sqrt(n)));
        new(x,y,n,dvar,L,LB,k)
    end
end
function getobjval(x,C)
    return [dot(x,C[1,:]),dot(x,C[2,:]),dot(x,C[3,:])]
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
        JuMP.fix(x[k],xx[k]; force=true)
    end
    optimize!(ap_m)
    if termination_status(ap_m) == MOI.OPTIMAL
        return true
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
        return 0;
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
function dominated(x,P)
    st = false
    for k=1:length(P)
        if all( x .>= P[k]) #&& any(x > P[k])
            st=true; break
        else
            continue
        end
    end

    return st
end

function Postpro2(candX,candY,newsol)
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
function fractionalFP(dvar,n,C,TL)
    X = copy(dvar); IGPair=[]; Tabu = []; newsol=0; t0=time(); #Y = copy(LB);
    while time()-t0 < TL
        I,G = sample(1:length(X), 2, replace=false)
        x1 = X[I]; x2 = X[G];
        λ = round(rand(Float64, 1)[1]; digits=1)
        x_t = x1*λ + x2*(1-λ); SearchDone = false; iter=0; Max_iter = 10
        while [I,G]∉IGPair && iter<Max_iter && SearchDone == false
            x_r = round.(Int,x_t)
            if ( (FBcheck(x_r,n) == true) && (x_r∉X) )
                push!(X,x_r); #Y = [Y; getobjval(x_r,C)'];
                newsol+=1; SearchDone = true
            else
                if x_r ∈ Tabu
                    x_r = flipoper(Tabu,x_t,x_r)
                    if x_r==[]
                        SearchDone = true
                    else
                        if ( (FBcheck(x_r,n) == true) && (x_r∉X) )
                            push!(X,x_r); #Y = [Y; getobjval(x_r,C)'];
                            newsol+=1; SearchDone = true;
                        end
                    end
                end
                if time()-t0 >= TL
                    break
                end
                if SearchDone == false
                    push!(Tabu,x_r)
                    x_t = fbsearch(x_r)
                    if x_t == 0 #when there's no new feasible lp sol
                        SearchDone = true
                    end
                end
            end
			iter+=1
        end
        push!(IGPair,[I,G])
    end
    return X,newsol
end
