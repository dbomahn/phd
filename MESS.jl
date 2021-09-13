using DataFrames,DelimitedFiles,JuMP,CSV

dir = "/home/ak121396/Desktop/instances/PublicInstances/"
data = readdlm(dir*ins[1], ';')

fdata = []
for i in [1,2,size(data)[1]-1]
    if occursin("=", data[i])==true
        dt = parse(Int,split(data[i],"= ")[2])
        push!(fdata, dt)
    end
end
for i in [3,4,5]
    if occursin("=", data[i])==true
        dt = split(data[i],"= ")[2]

        push!(fdata, dt)
    end
end
sc = zeros(Int,(fdata[2],fdata[1]))

for i=6:(5+fdata[2])
    tmp = split(data[i],r"[|,]")
    for j=2:1+fdata[1]
        sc[i-5,j-1] = parse(Int,tmp[j])
    end
end

mutable struct Data
    filepath::String; W::Int; S::Int; C::Array{}; F::Array{}; G::Array{};
    SC::Array{}; IC::Int; IP::Array{}

    function Data(filepath)
        data = readdlm(filepath, ';' )
        fdata = []
        for i in [1,2,size(data)[1]-1]
            if occursin("=", data[i])==true
                tmp = parse(Int,split(data[i],"= ")[2])
                push!(fdata, tmp)
            end
        end
        for i in [3,4,5]
            if occursin("=", data[i])==true
                tmp = split(data[i],"= ")[2]
                push!(fdata, tmp)
            end
        end
        W = fdata[1]; S = fdata[2]; IC = fdata[3];

        C = []; F= []
        for i=4:5
            tmp = split(fdata[i],r"[],[]")
            if i == 4
                for j=2:1+W
                    push!(C,parse(Int,tmp[j]))
                end
            else
                for j=2:1+W
                    push!(F,parse(Int,tmp[j]))
                end
            end
        end

        tp = split(fdata[6],r"[],[]")
        G = [ parse(Int,tp[x]) for x=2:S+1]
        SC = zeros(Int,(fdata[2],fdata[1]))

        for i=6:(5+S)
            tmp = split(data[i],r"[|,]")
            for j=2:(1+W)
                SC[i-5,j-1] = parse(Int,tmp[j])
            end
        end

        tm = split(data[size(data)[1]],r"[|]")
        IP = []
        for i=2:1+IC
            tm2 = parse.(Int,split(tm[i], ","))
            push!(IP,tm2)
        end

        new(filepath,W,S,C,F,G,SC,IC,IP)
    end
end

data = Data(dir*ins[1])
