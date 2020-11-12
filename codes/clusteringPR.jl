using DataStructures,DelimitedFiles,DataFrames,LinearAlgebra,Clustering

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

# paths = ("/home/ak121396/Desktop/Bensolve_KP/RoundDown/X/", "/home/ak121396/Desktop/Bensolve_KP/data/")
dvfile = readdir(paths[1]);  ins = readdir(paths[3]);
i=13
pathY= "/home/ak121396/Desktop/Bensolve_KP/Y/"; objfile = readdir(pathY);
# kp = Data(paths[1]*file[i],paths[2]*ins[i])
points = transpose(readdlm(pathY*objfile[i])[:,2:end])

dbscan(points, 0.05, min_neighbors = 3, min_cluster_size = 20)
