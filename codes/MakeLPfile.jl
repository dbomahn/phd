using DelimitedFiles,CPLEX,LinearAlgebra,JuMP,MathOptInterface,MathOptFormat

mutable struct Data
    input::String
    i::Int; j::Int;
    fixcost::Array{}
    demand::Array{}
    cost::Array{}
    cost2::Array{}

    function Data(input::String)
        data=readdlm(input)
        i,j = filter!(a->typeof(a)==Int, data[1,:]) # facility_i, customer_j # for FLP
        fixcost,capa,cost,cost2 = [],[],[],[]
        for k=2:i+1
            fix,cp = filter!(a->typeof(a)==Int, data[k,:])
            push!(fixcost,fix)
        end
        demand = filter!(a->typeof(a)==Int, data[2+i,:])
        a,b = size(data)
        for k=i+3:a
            ct = filter!(a->typeof(a)==Int, data[k,:])
            push!(cost2,ct)
            for l=1:length(ct)
                push!(cost,ct[l])
            end
        end
        new(input,i,j,fixcost,demand,cost,cost2)
    end
end
path = "/home/ak121396/multiobjective/instances/triFLP/instances/"
files = readdir(path)


function writemodel(model::Model, filename::String)
     if endswith(lowercase(filename), ".lp")
         file = MathOptFormat.LP.Model()
     elseif endswith(lowercase(filename), ".mps")
         file = MathOptFormat.MPS.Model()
     else
         println("Unknown file extension to write model: ",
                 split(filename, ".")[end])
         exit(8)
     end
     MOI.copy_to(file, backend(model))
     MOI.write_to_file(file, filename)
     println("wrote model to ", filename)
end


######################### Remove redundant characters from FLP LP file
for i=1:120
    data = Data(path*files[i])
    flp = JuMP.Model(CPLEX.Optimizer)
    @variable(flp, x[1:data.i,1:data.j], Bin)
    @variable(flp, y[1:data.i], Bin)
    @variable(flp, z[1:data.j], Bin)
    @objective(flp, Min, sum(data.cost2[a][b]*x[a,b] for a=1:data.i for b=1:data.j))
    # @variable(flp, dummy==1)
    # @constraint(flp, o1,sum(data.cost2[a][b]*x[a,b] for a=1:data.i for b=1:data.j) <=1)
    @constraint(flp, c[b=1:data.j], sum(x[a,b] for a in 1:data.i) == z[b])
    @constraint(flp, c1[a=1:data.i,b=1:data.j], x[a,b]<=y[a])
    @constraint(flp, o2, -sum(data.demand[b]*z[b] for b=1:data.j) == 0)
    @constraint(flp, o3, sum(data.fixcost[a]*y[a] for a=1:data.i) == 0)

    fname = files[i]
    writemodel(flp,"/home/ak121396/Desktop/FPBH/flp_lp/"*"$fname"*".lp")

    dt0 = readlines("/home/ak121396/Desktop/FPBH/flp_lp/"*"$fname"*".lp")
    dt1 = lowercase.(dt0)
    dt2 = replace.(dt1, ['_','[',']'] =>"")
    writedlm("/home/ak121396/Desktop/FPBH/flp_lp/"*"$fname"*".lp", dt2)
    # print

end


#############################       ILP MODEL      ###############################
for i=1:length(files)
    f = readdlm(path*files[i], '\t', String, '\n')
    global obj=parse(Int, f[1])
    global n=parse(Int,f[2])
    global m=parse(Int,f[3])

    c = f[4:obj+3]
    P= ones(obj,n)
    P=round.(Int,P)

    global zp=0;
    for x=1:length(c)
        a = split(c[x],('[',']',','))
        aa = filter!(e->e ∉ [ "" ,"[","," ,"]"] ,a)
        for y = 1:length(aa)
            p = parse(Int64,aa[y])
            P[x,y] = p
            if p==0
                global zp = zp+1;
            end
        end
    end
    ####### technical coefficient (a) ###########
    TC = f[obj+4:length(f)-1]
    global za=0;
    a= ones(m,n)
    a=round.(Int,a)
    for x=1:length(TC)
        t = split(TC[x],('[',']',','))
        tt = filter!(e->e ∉ [ "" ,"[","," ,"]"] ,t)
        for y=1:length(tt)
            tc=parse(Int64,tt[y])
            a[x,y] = tc
            if tc==0
                global za = za+1;
            end
        end
    end
    ####### RHS values (b) ##########
    b = ones(1,m)
    b = round.(Int,b)
    r = f[length(f)]
    r1 = split(r, ('[',']',','))
    r2 = filter!(e->e ∉ ["" ,"[", "]"] ,r1)
    for i=1:m
        b[i] = parse(Int64,r2[i])
    end

    ################################       MIP MODEL   ###############################
    mip = Model(CPLEX.Optimizer)
    @variable(mip, x[1:n]>=0)
    @variable(mip, y==1 )
    @objective(mip, Min, y)
    for i=1:m
            @constraint(mip, sum(a[i,j]*x[j] for j=1:n) <= b[i])
    end
    @constraint(mip, -sum(P[1,j]*x[j] for j=1:n) <= 1)
    @constraint(mip, -sum(P[2,j]*x[j] for j=1:n) <= 1)
    @constraint(mip, -sum(P[3,j]*x[j] for j=1:n) <= 1)
    # print(mip)
    # optimize!(mip)
    # writemodel(mip,"C:\\cygwin64\\home\\AK121396\\multiobjective\\instances\\Kirlik\\ILP\\"*files[i]*".lp)
#
    writemodel(mip,path*files[i]*".lp")

end

#############################       KP MODEL      ###############################
path = "/home/ak121396/Desktop/KP/data/"
f = readdir(path)

for i=1:length(f)
    fl = readdlm(path*f[i], '\t', String, '\n')

    obj=parse(Int,fl[1])
    n=parse(Int,fl[2])
    ub=parse(Int,fl[3])
    ########################    KP     #####################
    b = fl[4:length(fl)-1]
    C= ones(length(b),n)
    C = round.(Int,C)
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

    kp_m = Model(CPLEX.Optimizer )#,CPX_PARAM_THREADS=1,CPX_PARAM_SCRIND=0))
    @variable(kp_m, x[1:n], Bin )
    # @variable(kp_m, y==1 )
    # @objective( kp_m, Min, y )
    @objective(kp_m, Min, -dot(C[1,:],x[:]))
    # @constraint( kp_m, <=1 )
    @constraint( kp_m, -sum(weight[i]*x[i] for i=1:n) >= -ub )
    @constraint( kp_m, -dot(C[2,:],x[:]) == 0 )
    @constraint( kp_m, -dot(C[3,:],x[:]) == 0 )
    writemodel(kp_m,path*"/kp_lp/"*f[i]*".lp")
end
##########################      AP MODEL          #########################
path = "/home/ak121396/Desktop/AP/data/"
files = readdir(path)

for i=1:length(files)
    f=readdlm(path*files[i], '\t', String, '\n')
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

    P2 = ones(obj,n*n)
    P2 = round.(Int,P2)

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
    B=[A1;A2]


    ap_m = Model(CPLEX.Optimizer)
    @variable( ap_m, x[1:n*n], Bin )
    # @variable(ap_m, y==1 )
    @objective( ap_m, Min, sum(P2[1,j]*x[j] for j=1:n*n) )
    # @constraint( ap_m, cons[i=1:2*n], sum(B[i,j]*x[j] for j=1:n*n) == 1 )
    # @constraint( ap_m,  >= 0 )
    for i=1:2*n
        @constraint( ap_m,  sum(B[i,j]*x[j] for j=1:n*n) == 1 )
    end
    @constraint( ap_m, sum(P2[2,j]*x[j] for j=1:n*n) == 0 )
    @constraint( ap_m, sum(P2[3,j]*x[j] for j=1:n*n) == 0)

    writemodel(ap_m,path*"/ap_lp/"*files[i]*".lp")
end



############################# ILP model ####################################
path = "/home/ak121396/Desktop/ILP/data/"
fff = readdir(path)
for i=1:length(fff)
    f = readdlm(path*fff[i], '\t', String, '\n')
    # f=readdlm(ARGS[1], '\t', String, '\n')
    global obj=parse(Int, f[1])
    global n=parse(Int,f[2])
    global m=parse(Int,f[3])
    ########################    MIP instances    #####################
    ##################objective function coefficient (P) matrix ################
    c = f[4:obj+3]
    P= ones(obj,n)
    P=round.(Int,P)

    global zp=0;
    for x=1:length(c)
        a = split(c[x],('[',']',','))
        aa = filter!(e->e ∉ [ "" ,"[","," ,"]"] ,a)
        for y = 1:length(aa)
            p = parse(Int64,aa[y])
            P[x,y] = p
            if p==0
                global zp = zp+1;
            end
        end
    end
    ####### technical coefficient (a) ###########
    TC = f[obj+4:length(f)-1]
    global za=0;
    a= ones(m,n)
    a=round.(Int,a)
    for x=1:length(TC)
        t = split(TC[x],('[',']',','))
        tt = filter!(e->e ∉ [ "" ,"[","," ,"]"] ,t)
        for y=1:length(tt)
            tc=parse(Int64,tt[y])
            a[x,y] = tc
            if tc==0
                global za = za+1;
            end
        end
    end
    ####### RHS values (b) ##########
    b = ones(1,m)
    b = round.(Int,b)
    r = f[length(f)]
    r1 = split(r, ('[',']',','))
    r2 = filter!(e->e ∉ ["" ,"[", "]"] ,r1)
    for i=1:m
        b[i] = parse(Int64,r2[i])
    end

    ################################       MIP MODEL   ###############################
    mip = Model(CPLEX.Optimizer)
    @variable(mip, 0<=x[1:n], Int )
    @objective(mip, Min, sum(P[1,j]*x[j] for j=1:n))
    for i=1:m
        @constraint(mip, -sum(a[i,j]*x[j] for j=1:n) >= -b[i])
    end
    @constraint(mip, c1, sum(P[2,j]*x[j] for j=1:n) == 0)
    @constraint(mip, c2, sum(P[3,j]*x[j] for j=1:n) == 0)
    writemodel(mip,path*"/ilp_lp/"*fff[i]*".lp")


    # Kirlik&Sayin LP file format
    # mip = Model(CPLEX.Optimizer)
    # @variable(mip, 0<=x[1:n], Int )
    # @variable(mip, dummy == 1 )
    # @objective(mip, Min, dummy)
    # for i=1:m
    #     @constraint(mip, -sum(a[i,j]*x[j] for j=1:n) >= -b[i])
    # end
    # @constraint(mip, o1, sum(P[1,j]*x[j] for j=1:n) <= 1 )
    # @constraint(mip, o2, sum(P[2,j]*x[j] for j=1:n) <= 1)
    # @constraint(mip, o3, sum(P[3,j]*x[j] for j=1:n) <= 1)
    # writemodel(mip, "/home/ak121396/Downloads/KirlikSayin2014/ILP/"*fff[i]*".lp")
    #
    # dt0 = readlines("/home/ak121396/Downloads/KirlikSayin2014/ILP/"*fff[i]*".lp")
    # dt1 = lowercase.(dt0)
    # dt2 = replace.(dt1, ['_','[',']'] =>"")
    # writemodel(mip, "/home/ak121396/Downloads/KirlikSayin2014/ILP/"*fff[i]*".lp")

end
