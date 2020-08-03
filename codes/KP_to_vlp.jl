using DelimitedFiles

pwd()
# cd("")
# path=pwd()
path = "/home/ak121396/Desktop/KP/"
# path2="C:/Users/okokq/Dropbox/JKU/ap/"
files=readdir(path*"/data/")
for i=1:length(files)
    ins = open(path*"//JKU/kp/"*files[i]*".vlp","a")
    f = readdlm(path*"/JKU/kp/KP/"*files[i], '\t', String, '\n')
    global obj=parse(Int, f[1])
    global k=parse(Int,f[2])
    global ub=parse(Int,f[3])

    ##################objective function coefficient (P) matrix ################
    b = f[4:length(f)-1]
    # b=filter!(e->e≠",",profit)
    P= ones(length(b),k)
    P=round.(Int,P)
    global ct=0;
    for x=1:length(b)
        a = split(b[x],('[',']',','))
        aa=filter!(e->e ∉ [ "" ,"[","," ,"]"] ,a)
        for y=1:length(aa)
            p=parse(Int64,aa[y])
            P[x,y] = p
            if p==0
                ct = ct+1;
            end
        end
    end


    ####### coefficient matrix (B) ##########
    weight =ones(1,k)
    weight = round.(Int,weight)
    item = f[length(f)]
    w1 = split(item, ('[',']',','))
    w2 = filter!(e->e ∉ ["" ,"[", "]"] ,w1)
    for i=1:k
        weight[i] = parse(Int64,w2[i])
    end

    ########################### vlp file for BENSOLVE ####################
    wholearray=[];
    arr=["p vlp max",k,k,k,obj,(k*obj)-ct]
    push!(wholearray,arr)
    for i=1:k
        if weight[i]!=0
            push!(wholearray,("a",1,i,weight[i]))
        end
    end

    for i=1:obj
        for j=1:k
            if P[i,j]!=0
                push!(wholearray,("o",i,j,P[i,j]))
            end
        end
    end

    for j=1:k
        push!(wholearray,("i",j,"u",ub))
    end


    for i=1:k
        if weight[i]!=0
            push!(wholearray,("j",i,"d",0,1))
        end
    end
    push!(wholearray,"e")
    writedlm(ins,wholearray)
    close(ins)
end
