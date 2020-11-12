using DataStructures,DelimitedFiles,DataFrames,LinearAlgebra,Clustering,StatsBase

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
# dvfile = readdir(paths[1]);  ins = readdir(paths[3]);
i=13
pathY = "E:\\Bensolve_KP\\Y\\"; objfile = readdir(pathY);
pathY= "/home/ak121396/Desktop/Bensolve_KP/Y/"; objfile = readdir(pathY);
# kp = Data(paths[1]*file[i],paths[2]*ins[i])
objs = readdlm(pathY*objfile[i])[:,2:end]
ind = findall(i-> -1 in LB[i,:]  , 1:size(LB)[1])
for i=1:3

    push!(LB, deleteat!(objs[:,i], ind[j] for j=1:length(ind)) )
end

LB = hcat(deleteat!(objs[:,1],ind[j] for j=1:length(ind)),deleteat!(objs[:,2],ind[j] for j=1:length(ind)),
    deleteat!(objs[:,3],ind[j] for j=1:length(ind)) )
maxobj = [maximum(LB[:,i]) for i=1:3]; minobj = [minimum(LB[:,i]) for i=1:3];
radius = round(Int,mean(maxobj-minobj)/length(LB[:,1]))
points = transpose(LB)
LBcluster = dbscan(points, radius, min_neighbors = 1, min_cluster_size = 1)



LBres = kmeans(points, 10;)# maxiter=200, display=:iter

scatter(,LBres.Petalwidth, marker_z=LBres.assignments,
        color=:lightrainbow, legend=true)
