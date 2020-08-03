using DelimitedFiles

# pwd()
# cd("")
# path=pwd()
path = "/home/ak121396/Desktop/ILP/"
# path2="C:/Users/okokq/Dropbox/JKU/ap/"
files=readdir(path*"/data/")

for i=1:length(files)
    ins = open(path*"/data/"*files[i]*".vlp","a")
    f = readdlm(path*"/data/"*files[i], '\t', String, '\n')
    global obj=parse(Int, f[1])
    global n=parse(Int,f[2])
    global m=parse(Int,f[3])

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

    ########################### vlp file for BENSOLVE ####################
    wholearray=[];
    arr=["p vlp max",m,n,(m*n)-za,obj,(obj*n)-zp]
    push!(wholearray,arr)

    for i=1:m
        for j=1:n
            if a[i,j]!=0
                push!(wholearray,("a", i,j,a[i,j]) )
            end
        end
    end

    # Parray=[];
    for i=1:obj
        for j=1:n
            if P[i,j]!=0
                push!(wholearray,("o",i,j,P[i,j]))
            end
        end
    end

    for i=1:m
        if b[i]!= 0
            push!(wholearray,("i",i,"u",b[i]))
        end
    end

    for i=1:m
        push!(wholearray, ("j",i,"l",0))
    end

    push!(wholearray,"e")
    # print(wholearray)
    writedlm(ins,wholearray)
    close(ins)
end
