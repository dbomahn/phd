using DataStructures,DelimitedFiles,DataFrames,CSV,StatsBase,Random,LinearAlgebra,CPUTime
struct path{S<:String}
    dir1::S; dir2::S;
end

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
function KPfbcheck(txi,n,weight,ub)
    result = false
    if sum(weight[i]*txi[i] for i=1:n) <= ub
        result = true
    else
        result = false
    end
    return result
end

function domFilter(sol,obj)
    copysol = Dict(); copyobj = Dict();
    for i=1:length(obj)
        copysol[i] = sol[i,:]
        copyobj[i] = obj[i]
    end

    for i=1:length(obj)-1
        for j=i+1:length(obj)
            if all(obj[i] .<= obj[j]) == true #dominated by PF[j]
                copyobj[i]=nothing; copysol[i]=nothing;break
            elseif all(obj[j] .<= obj[i]) == true
                copyobj[j]=nothing; copysol[j]=nothing;
            end
        end
    end

    finalsol = filter!(a->a!=nothing, collect(values(copysol)))
    finalobj = filter!(a->a!=nothing, collect(values(copyobj)))

    return finalsol,finalobj
end
function dominated(x,P)
    st = false
    for k=1:length(P)
        if all( x .<= P[k]) && any(x < P[k])
            st=true; break
        else
            continue
        end
    end
    return st
end
pathes = ("/home/ak121396/Desktop/BENoutputs/KP/X/", "/home/ak121396/Desktop/BENKP/data/")
file = readdir(pathes[1]); ins = readdir(pathes[2])
file[13], ins[12]
