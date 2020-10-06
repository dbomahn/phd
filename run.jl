
######################### KP ###################

for i=1:100
    fl = readdlm(dir1*data[i],'\t', String, '\n')
    var = round.(readdlm(dir2*x[i]), digits=4) #round numbers till 2 digits
    objs = round.(readdlm(dir3*y[i])[:,2:end], digits=4)
    global obj=parse(Int,fl[1])
    global n=parse(Int,fl[2])
    global ub=parse(Int,fl[3])

    ########################    KP     #####################
    global b = fl[4:length(fl)-1]
    C= ones(length(b),n)
    global C = round.(Int,C)
    global ct=0;
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
    ######################## coefficient matrix (B) #############################
    weight =ones(1,n)
    global weight = round.(Int,weight)
    item = fl[length(fl)]
    w1 = split(item, ('[',']',','))
    w2 = filter!(e->!(e in ["" ,"[", "]"]) ,w1)
    for i=1:n
        weight[i] = parse(Int64,w2[i])
    end
    #######################################
    s,ss = size(var)
    indx= findall(x->x!=0,[sum(var[i,:]) for i=1:s]) #NZ sol
    lpX = [var[i,:] for i in indx]

    pf = hcat(objs[:,1],objs[:,2],objs[:,3])
    L = reshape([pf[i,:] for i in indx],length(lpX),1)
    L2 = transpose(hcat(L...))
    fmin = [minimum(L2[:,i]) for i=1:3]
    fmax = [maximum(L2[:,i]) for i=1:3]
    steps = [round.(Int,abs(fmax[i]-fmin[i])/length(L)) for i=1:3] #determine steps according to #customers(j)
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
    Grouping(groups[1:min(length(groups),5)],groupkeys[1:min(length(groups),5)],lpX,steps,fmin,L2,n)

    GroupingTime = @CPUelapsed Xf,cand,loca_check,LPcount = Grouping(groups,groupkeys,lpX,steps,fmin,L2,n)
    Xf1 = copy(Xf)
    candX = vcat(collect(values(cand)),loca_check)
    PFset = giveobjval.(collect(values(Xf1)))

    FPTime = @CPUelapsed Xf2,PFset2,Tabu,FPLPcount = GroupFP(Xf1,PFset,candX,LPcount,steps,fmin,fmax)
    P = convert.(Array{Int,1},collect(values(Xf2)))
    Pcoordi = convert.(Array{Int,1},collect(keys(Xf2)))
    Pobj = convert.(Array{Int,1}, giveobjval.(P))
    #Filter dominated solutions
    FPsol, FPPF = domFilter(P,Pobj)
    Pz,Py,Px=[hcat(FPPF...)[i,:] for i=1:3]; dfP = Px,Py,Pz
    matP = zeros(length(FPsol),n)
    for i=1:length(FPsol)
        matP[i,:] = FPsol[i]
    end

    # ins = data[i][1:end-10]
    record1 = DataFrame( LP=FPLPcount, sol=length(Pz), CPUtime=GroupingTime+FPTime )
    GFPX=DataFrame(matP); GFPY=DataFrame(dfP);

    CSV.write("/home/ak121396/Desktop/phd/KPresults//"*data[i][1:12]*"_record.csv",record1, append=true, header=false )#, delim=',' )
    CSV.write("/home/ak121396/Desktop/phd/KPresults/Y/"*data[i][1:end-4]*".txt",GFPY, header=false, delim=' ' )
end




######################### AP ###################

for i=1:100

    f=readdlm(dir1*data[i], '\t', String, '\n')
    var = round.(readdlm(dir2*x[i]), digits=4) #round numbers till 2 digits
    objs = round.(readdlm(dir3*y[i])[:,2:end], digits=4)
    obj=parse(Int,f[1])
    global n=parse(Int,f[2])
    P=[]
    for i=3:length(f)
        a = split(f[i], ('[',']',','))
        a = filter!(e->e ∉ [ "" ,"[","," ,"]"] ,a)
        if length(a)>2
            push!(P,a)
        end
    end

    P2 = ones(obj,n*n)
    global P2 = round.(Int,P2)

    global ct=0;
    for x=1:length(P)
        for y=1:n
            p=parse(Int64,P[x][y])
            idx = Int(floor((x-1)/n))
            ind = ((x-1)*n+y)%(n*n)
            if ind != 0
                P2[idx+1,((x-1)*n+y)%(n*n)] = p
            else
                P2[idx+1,(n*n)] = p
            end
            if p==0
                ct =ct+1;
            end
        end
    end
    ######################## coefficient matrix (B) #############################
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
    global B=[A1;A2]
    s,ss = size(var)
    indx= findall(x->x!=0,[sum(var[i,:]) for i=1:s]) #NZ sol
    lpX = [var[i,:] for i in indx]

    pf = hcat(objs[:,1],objs[:,2],objs[:,3])
    L = reshape([pf[i,:] for i in indx],length(lpX),1)
    L2 = transpose(hcat(L...))
    fmin = [minimum(L2[:,i]) for i=1:3]
    fmax = [maximum(L2[:,i]) for i=1:3]
    steps = [round.(Int,abs(fmax[i]-fmin[i])/length(L)*n) for i=1:3] #determine steps according to #customers(j)
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
    # Grouping(groups[1:min(length(groups),5)],groupkeys[1:min(length(groups),5)],lpX,steps,fmin,L2,n)

    GroupingTime = @CPUelapsed Xf,cand,loca_check,LPcount = Grouping(groups,groupkeys,lpX,steps,fmin,L2,n)
    Xf1 = copy(Xf)
    candX = vcat(collect(values(cand)),loca_check)
    PFset = giveobjval.(collect(values(Xf1)))

    FPTime = @CPUelapsed Xf2,PFset2,Tabu,FPLPcount = GroupFP(Xf1,PFset,candX,LPcount,steps,fmin,fmax)
    P = convert.(Array{Int,1},collect(values(Xf2)))
    Pcoordi = convert.(Array{Int,1},collect(keys(Xf2)))
    Pobj = convert.(Array{Int,1}, giveobjval.(P))
    #Filter dominated solutions
    FPsol, FPPF = domFilter(P,Pobj)
    Pz,Py,Px=[hcat(FPPF...)[i,:] for i=1:3]; dfP = Px,Py,Pz


    # ins = data[i][1:end-10]
    record1 = DataFrame( LP=FPLPcount, sol=length(Pz), CPUtime=GroupingTime+FPTime )
    # GFPX=DataFrame(matP);
    GFPY=DataFrame(dfP);

    CSV.write("/home/ak121396/Desktop/phd/APresults/AP_"*data[i][1:12]*"_record.csv",record1, append=true, header=false )#, delim=',' )
    CSV.write("/home/ak121396/Desktop/phd/APresults/Y/"*data[i][1:end-4]*"_Y_.txt",GFPY, header=false, delim=' ' )
end
