using DelimitedFiles

pwd()
# cd("")
# path=pwd()
path = "/home/ak121396/Desktop/AP/"
# path2="C:/Users/okokq/Dropbox/JKU/ap/"
files=readdir(path*"/data/")

for i=1:length(files)
    ins = open(path*"//JKU/ap/"*files[i]*".vlp","a")
    f = readdlm(path*"/JKU/ap/AP/"*files[i], '\t', String, '\n')
    global obj=parse(Int, f[1])
    global k=parse(Int,f[2])
    ##################objective function coefficient (P) matrix ################
    P=[]
    for i=3:length(f)
        a = split(f[i], ('[',']',','))
        a = filter!(e->e ∉ [ "" ,"[","," ,"]"] ,a)
        if length(a)>2
            push!(P,a)
        end
    end

    P2 = ones(obj,k*k)
    P2 = round.(Int,P2)

    global ct=0;
    for x=1:length(P)
        for y=1:k
            p=parse(Int64,P[x][y])
            idx = Int(floor((x-1)/k))
            ind = ((x-1)*k+y)%(k*k)
            if ind != 0
                P2[idx+1,((x-1)*k+y)%(k*k)] = p
            else
                P2[idx+1,(k*k)] = p
            end
            if p==0
                ct =ct+1;
            end
        end
    end


    ####### coefficient matrix (B) ##########
    A1 =zeros(k,k*k)
    A2 =zeros(k,k*k)
    A1  = round.(Int,A1)
    A2 = round.(Int,A2)

    for i=1:k
        for j=1:k
            A1[i,((i-1)*k)+j] = 1
                if A1[i,((i-1)*k)+j]==1
                end
        end
    end

    for i=1:k
        for j=1:k
            A2[i,((j-1)*k)+i] = 1
            if A2[i, ((j-1)*k)+i]==1
            end
        end
    end
    B=[A1;A2]



    ########################### vlp file for BENSOLVE ####################
    wholearray=[];
    arr=["p vlp min",2*k,k*k,2*k*k,obj,(k*k*obj)-ct]
    push!(wholearray,arr)
    for i=1:2*k
        for j=1:k^2
            if B[i,j]!=0
                push!(wholearray,("a",i,j,B[i,j]))
            end
        end
    end

    # Parray=[];
    for i=1:obj
        for j=1:k*k
            if P2[i,j]!=0
                push!(wholearray,("o",i,j,P2[i,j]))
            end
        end
    end

    for j=1:2*k
        push!(wholearray,("i",j,"s",1))
    end

    for i=1:k*k


        push!(wholearray,("j", i,'d',0,1))
    end
    push!(wholearray,"e")

    writedlm(ins,wholearray)
    close(ins)
end
