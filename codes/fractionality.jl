using StatsBase,JLD2,DelimitedFiles
mutable struct Valu
    x::String; y::String; dvar::Array{}; LB::Array{}; LBmtx::Array{};
    function Valu(x,y)
        JLD2.@load x dv;
        dv0 = Array(dv);
        dv1 = round.(dv0; digits=4);
        objs = round.(digits=4, readdlm(y));
        ind = findall(i-> 0 in objs[i,:]  , 1:size(objs)[1]);
        dv2 = dv1[setdiff(1:end, ind), :];
        LBmtx = objs[setdiff(1:end, ind), 2:end];
        dvar = [Vector(dv2[i, :]) for i = 1:size(LBmtx)[1]];
        LB = [Vector(LBmtx[i, :]) for i = 1:size(LBmtx)[1]];
        new(x,y,dvar,LB,LBmtx)
    end
end


x = readdir("/home/k2g00/k2g3475/clusterhome/miplib/vlp/final/x/")
y = readdir("/home/k2g00/k2g3475/clusterhome/miplib/vlp/final/y/")

avg = []
for i=1:5
    pre = Valu("/home/k2g00/k2g3475/clusterhome/miplib/vlp/final/x/"*x[i],"/home/k2g00/k2g3475/clusterhome/miplib/vlp/final/y/"*y[i])
    # l = length(pre.dvar[1])
    l = length(pre.dvar)
    ilist = []
    for i=1:length(pre.dvar)
        # ct = count(i->(0<i<1), pre.dvar[i]) #(i!=1&&i!=0)
        # push!(ilist,(ct/l)*100)
        idx = findall(i->(i!=0), pre.dvar[i])
        push!(ilist, [idx])
    end
    # push!(avg, mean(ilist))
    uq = union(ilist)
    push!(avg,(length(uq)/l)*100)
end

mean(avg)


pre = Valu("/home/k2g00/k2g3475/clusterhome/miplib/vlp/final/x/"*x[1],"/home/k2g00/k2g3475/clusterhome/miplib/vlp/final/y/"*y[1])

D1 = Dict(); D0 = Dict(); Df = Dict()
for i=1:length(pre.dvar[1])
    D1[i] = 0; D0[i] = 0; Df[i] = 0
end
l = length(pre.dvar[1])
for i=1:length(pre.dvar)
    for j=1:length(pre.dvar[1])
        if pre.dvar[i][j]==1
            D1[j] = D1[j]+1
        elseif pre.dvar[i][j]==0
            D0[j] = D0[j]+1
        else
            Df[j]= Df[j]+1
        end
    end
end

cd1 = collect(values(D1)); cd0 = collect(values(D0)); cdf = collect(values(Df));
(count(i->(i>l*0.4),cd1)/l)*100
(count(i->(i>l*0.4),cd0)/l)*100
(count(i->(i>l*0.4),cdf)/l)*100


maximum([1,2])
union([[1,2]],[[3,3]],[[1,2]])



#######################    FLP    ########################
x = readdir("/home/k2g00/k2g3475/clusterhome/pastinstance/flp/varval/")
y = readdir("/home/k2g00/k2g3475/clusterhome/pastinstance/flp/PF/")
mutable struct Vall
    x::String; y::String; dvar::Array{}; LB::Array{}; LBmtx::Array{};
    function Vall(x,y)
        dv = round.(digits=4, readdlm(x))
        objs = round.(digits=4, readdlm(y))
        ind = findall(i-> 0 in objs[i,:]  , 1:size(objs)[1])
        dv2 = dv[setdiff(1:end, ind), :];
        LBmtx = objs[setdiff(1:end, ind), 2:end];
        dvar = [Vector(dv2[i, :]) for i = 1:size(LBmtx)[1]]
        LB = [Vector(LBmtx[i, :]) for i = 1:size(LBmtx)[1]]
        new(x,y,dvar,LB,LBmtx)
    end
end



xpath = "/home/k2g00/k2g3475/clusterhome/pastinstance/ap/gpr/X/"
ypath = "/home/k2g00/k2g3475/clusterhome/pastinstance/ap/gpr/Y/"
x = readdir(xpath)
y = readdir(ypath)

f = [] ; z=[]; o =[];
for i=1:length(x)
    pre = Vall(xpath*x[i],ypath*y[i])
    D1 = Dict(); D0 = Dict(); Df = Dict()
    for i=1:length(pre.dvar[1])
        D1[i] = 0; D0[i] = 0; Df[i] = 0
    end
    l = length(pre.dvar[1])
    for i=1:length(pre.dvar)
        for j=1:length(pre.dvar[1])
            if pre.dvar[i][j]==1
                D1[j] = D1[j]+1
            elseif pre.dvar[i][j]==0
                D0[j] = D0[j]+1
            else
                Df[j]= Df[j]+1
            end
        end
    end

    cd1 = collect(values(D1)); cd0 = collect(values(D0)); cdf = collect(values(Df));
    push!(o,(count(i->(i>l*0.4),cd1)/l)*100); push!(z,(count(i->(i>l*0.4),cd0)/l)*100); push!(f,(count(i->(i>l*0.4),cdf)/l)*100)
end
mean(f)
mean(o)
mean(z)


x = readdir("/home/k2g00/k2g3475/clusterhome/pastinstance/ap/gpr/X/")
y = readdir("/home/k2g00/k2g3475/clusterhome/pastinstance/ap/gpr/Y/")
avg = []
for i=1:100
    pre = Vall("/home/k2g00/k2g3475/clusterhome/pastinstance/ap/gpr/X/"*x[i],"/home/k2g00/k2g3475/clusterhome/pastinstance/ap/gpr/Y/"*y[i])
    l = length(pre.dvar[1])
    # l = length(pre.dvar)
    ilist = []
    for i=1:length(pre.dvar)
        ct = count(i->(0<i<1), pre.dvar[i]) #(i!=1&&i!=0)
        push!(ilist,(ct/l)*100)
        # idx = findall(i->(i!=0), pre.dvar[i])
        # push!(ilist, [idx])
    end
    push!(avg, mean(ilist))
    # uq = union(ilist)
    # push!(avg,(length(uq)/l)*100)
end

mean(avg)

rt = []
for i=1:120
    f = round.(readdlm(files[i]); digits=2)
    ct = 0
    for j=1:size(f)[1]
        if any(i->0<i<1, f[j,:])==true
            ct+=1
        end
    end
    push!(rt,ct/(size(f)[1]-3))
end

        
