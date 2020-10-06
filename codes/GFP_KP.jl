pathes=path("/home/ak121396/Desktop/BEN+FP/data/", "/home/ak121396/Desktop/BENoutputs//KP/X/","/home/ak121396/Desktop/BENoutputs//KP/Y/")
dir1 = pathes.dir1; dir2 = pathes.dir2; dir3 = pathes.dir3
data = readdir(dir1)
x = readdir(dir2)
y = readdir(dir3)

timerecord = []
for i=11:100
    fl = readdlm(dir1*data[i],'\t', String, '\n')
    vari= round.(readdlm(dir2*x[i]), digits=4) #round numbers till 2 digits
    objs = round.(readdlm(dir3*y[i])[:,2:end], digits=4)

    # fl = readdlm(ARGS[1], '\t', String, '\n')
    # vari = round.(readdlm(ARGS[2]), digits=4)
    # objs = round.(readdlm(ARGS[3])[:,2:end] , digits=4)

    obj=parse(Int,fl[1])
    global n=parse(Int,fl[2])
    global ub=parse(Int,fl[3])
    ########################    KP     #####################
    b = fl[4:length(fl)-1]
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
    weight = round.(Int,weight)
    item = fl[length(fl)]
    w1 = split(item, ('[',']',','))
    w2 = filter!(e->!(e in ["" ,"[", "]"]) ,w1)
    for i=1:n
        weight[i] = parse(Int64,w2[i])
    end

    s,ss = size(vari)
    indx= findall(x->x!=0,[sum(vari[i,:]) for i=1:s]) #NZ sol
    lpX = [vari[i,:] for i in indx]

    pf = hcat(objs[:,1],objs[:,2],objs[:,3])
    L = reshape([pf[i,:] for i in indx],length(lpX),1)
    L2 = transpose(hcat(L...))
    fmin = [minimum(L2[:,i]) for i=1:3]
    fmax = [maximum(L2[:,i]) for i=1:3]
    steps = [1,1,1]#[round.(Int,abs(fmax[i]-fmin[i])/length(L)) for i=1:3] #determine steps according to #customers(j)
    cube = Dict();

    for iter=1:length(L)
        loca = [round.(Int,((L2[iter,k]-fmin[k])/steps[k])+1) for k=1:3]
        if !haskey(cube,loca)
            cube[loca] = [iter]
        else
            push!(cube[loca], iter)
        end
    end

    #############################    Group Feasibility Pump   ############################
    groups = collect(values(cube)); groupkeys = collect(keys(cube))
    ########################### TEST ########################
    # Grouping(groups,groupkeys,lpX,steps,fmin,L2,n)
    ########################################################
    # GroupingTime = @CPUelapsed Xf,cand,loca_check,LPcount = Grouping(groups,groupkeys,lpX,steps,fmin,L2,n)
    # Xf1 = copy(Xf)
    # candX = vcat(collect(values(cand)),loca_check)
    Xf1 = Dict()
    candX = lpX #let's try to consider everything

    PFset = giveobjval.(collect(values(Xf1)))
    function GroupFP(Xf,PFset,candX,steps,fmin,fmax)
        Tabu = []; elaps=0;LPcount=0
        algo_start = CPUtime_us()
        for k=1:length(candX)
            # print("===============Feasi Pump",k," th candidate sol ==============\n")
            x_t = candX[k];
            SearchDone = false
            itr = 1
            Max_itr = n*5 #max(count(x->0<x<1,x_t),1)  #maximum number of attempts => How to set
            while itr<Max_itr && SearchDone==false
                xi_t = round.(Int,x_t)
                fxi = giveobjval(xi_t)
                if ( (KPfbcheck(xi_t) == true) && (dominated(fxi,PFset)==false) ) #checking feasibility and dominance   #z[1]<=fmax[1]-1 && z[2]<=fmax[2]-1 && z[3]<=fmax[3]-1)
                    idx = getlocation(xi_t,fmin,steps)
                    push!(PFset,fxi) #add new solval to PFset
                    if !haskey(Xf,idx)
                        # print("ROUNING => another solution added \n")
                        Xf[idx] = xi_t
                    else
                        # print("ROUNING =>new intsol to a cuboid \n")
                        push!([Xf[idx]], xi_t)
                    end
                    SearchDone = true
                else
                    if xi_t ∈ Tabu
                        xi_t = flipoper(Tabu,x_t,xi_t)
                        if xi_t==[]
                            # print("FLIP didn't work \n")
                            SearchDone = true
                        else
                            fxi = giveobjval(xi_t)
                            if ( (KPfbcheck(xi_t) == true) && (dominated(fxi,PFset)==false) )  #z[1]<=fmax[1]-1 && z[2]<=fmax[2]-1 && z[3]<=fmax[3]-1 )
                                idx = getlocation(xi_t,fmin,steps)
                                push!(PFset,fxi) #add new solval to PFset
                                if !haskey(Xf,idx)
                                    # print("FLIP=> new cuboid new solution \n")
                                    Xf[idx] = xi_t
                                else
                                    # print("FLIP=> new intsol to a cuboid \n")
                                    push!([Xf[idx]], xi_t)
                                end
                                SearchDone = true
                            end
                        end
                    end
                    if SearchDone == false
                        push!(Tabu,xi_t) #when break
                        idx = getlocation(xi_t,fmin,steps)
                        x_t = fbsearch(idx,steps,fmin,xi_t,fmax)
                        LPcount+=1
                        if x_t == 0 #when there's no new feasible lp sol
                            # print("no lp sol's found");
                            break
                        # else
                            # print("\n New lp obj: ",giveobjval(x_t), "\n")
                        end
                    end
                end
                itr+=1; #@show itr
                if CPUtime_us()-algo_start < TimeLimit*1e6
                    continue
                else
                    return Xf,PFset,Tabu,LPcount
                end
            end
        end
        return Xf,PFset,Tabu,LPcount
    end
    #Time Feasibility Pump
    FPTime = @CPUelapsed Xf2,PFset2,Tabu,FPLPcount = GroupFP(Xf1,PFset,candX,steps,fmin,fmax)
    push!(timerecord,FPTime)
    P = convert.(Array{Int,1},collect(values(Xf2)))
    Pcoordi = convert.(Array{Int,1},collect(keys(Xf2)))
    Pobj = convert.(Array{Int,1}, giveobjval.(P))
    # push!(Pobj,[0,0,0]);# push!(P,zeros(Int,n))
    #Filter dominated solutions
    FPsol, FPPF = domFilter(P,Pobj)

    Pz,Py,Px=[hcat(FPPF...)[i,:] for i=1:3]; dfP = Px,Py,Pz
    ndpoint = DataFrame(dfP)
    fname = data[i][1:19]
    CSV.write("/home/ak121396/Desktop/BEN+FP/Y/"*"$fname"*"Y.log",ndpoint,header=false, delim=' ' )
    matP = zeros(Int,length(FPsol),n)
    for i=1:length(FPsol)
        matP[i,:] = FPsol[i]
    end
    sols = DataFrame(matP)
    CSV.write("/home/ak121396/Desktop/BEN+FP/X/"*"$fname"*"X.log",sols,header=false, delim=' ' )
end

#########################  Record outputs  ############################
# ins = ARGS[1][end-12:end-7] #CPUtime recorded, naming after the instance
# colname = ARGS[1][end-12:end-4]
# record1 = DataFrame( LP=FPLPcount, sol=length(Pz), CPUtime=GroupingTime+FPTime )
# # insertcols!(record1,3, Symbol("$colname")=>GroupingTime+FPTime)
# GFPX=DataFrame(matP); GFPY=DataFrame(dfP);
# CSV.write("/home/ak121396/Desktop/phd/GFPresults/KP_"*"$ins"*"_record.csv",record1, append=true, writeheader=false )#, delim=',' )
# #CSV.write(ARGS[1]*"_GFP_X.csv",GFPX, append=true, writeheader=false)
# CSV.write("/home/ak121396/Desktop/phd/GFPresults/Y/"*ARGS[end-12:end-4]*_"KP_Y_.txt",GFPY, writeheader=false, delim=' ' )
# print(colname," GroupFP Done!")
