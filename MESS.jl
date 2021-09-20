using DataFrames,DelimitedFiles,JuMP,LinearAlgebra

@show file = ARGS[1][1:end-4]

# dir = "/home/ak121396/Desktop/instances/PublicInstances/"
# data = readdlm(dir*ins[1], ';')
file = "C:/Users/AK121396/Desktop/phd/toy.dzn"
# data = readdlm(dir, ';')

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
        SC = reshape(transpose(SC), (1,S*W))

        tm = split(data[size(data)[1]],r"[|]")
        IP = []
        for i=2:1+IC
            tm2 = parse.(Int,split(tm[i], ","))
            push!(IP,tm2)
        end

        new(filepath,W,S,C,F,G,SC,IC,IP)
    end
end

data = Data(file)


M = 9999999999999999
###########################   Mathematical Model  #############################
# wlp = Model() #GLPK.Optimizer
# @variable( wlp, y[1:data.W] ) #,binary=true
# @variable( wlp, x[1:data.S,1:data.W])
# @variable( wlp, z[1:data.S,1:data.W])
# for w=1:data.W
#     @constraint( wlp, sum(x[s,w] for s=1:data.S) <= data.C[w] )
# end
# for s=1:data.S
#     @constraint( wlp, sum(x[s,w] for w=1:data.W) <= data.G[s])
# end
# for w=1:data.W
#     @constraint( wlp, sum(x[s,w] for s=1:data.S) <= (M*y[w]) )
# end
# for w=1:data.W
#     for i=1:data.IC
#         @constraint( wlp, x[data.IP[i][1],w] <= M*z[data.IP[i][1],w] )
#         @constraint( wlp, x[data.IP[i][2],w] <= M*(1-z[data.IP[i][1],w]) )
#     end
# end
# @objective( wlp, Min, sum(data.F[w]*y[w] for w=1:data.W) + sum(sum(data.SC[s,w]*x[s,w] for w=1:data.W) for s=1:data.S) )
################################  convert to vlp file  ###############################################

A1 = zeros(Int,data.W,data.S*data.W+data.W+data.W*data.S)
for w=1:data.W
    for s=1:data.S
        A1[w,data.W*(s-1)+w] = 1
    end
end

A2 = zeros(Int,data.S,data.S*data.W+data.W+data.W*data.S)
for s=1:data.S
    for w=1:data.W
        A2[s,w+data.W*(s-1)] = 1
    end
end

A3 = copy(A1)
for w=1:data.W
    A3[w,data.S*data.W+w] = -M
end
# A32 = Matrix{Int}(I, data.W, data.W)*(-M)

A41= zeros(Int,data.IC*data.W,data.S*data.W+data.W)
A42 = zeros(Int,data.IC*data.W,data.S*data.W)
A51 = copy(A41); A52 =copy(A42)
for i=1:data.IC
    for w=1:data.W
        k = data.W*(i-1)+w
        A41[k, data.W*(data.IP[i][1]-1)+w] = 1
        A42[k, data.W*(data.IP[i][1]-1)+w] = -M
    end
    for w=1:data.W
        k = data.W*(i-1)+w
        A51[k, data.W*(data.IP[i][2]-1)+w] = 1
        A52[k, data.W*(data.IP[i][2]-1)+w] = M
    end
end

A4 = hcat(A41, A42)
A5 = hcat(A51, A52)
P = hcat(data.SC,transpose(data.F), zeros(Int, (1,data.S*data.W)) )

B = [A1;A2;A3;A4;A5]
m,n = size(B)
nz = count(i->(i!=0),B)
objnz = count(i->(i!=0),P)
RHS = [data.C;data.G;zeros(Int,(data.W+data.W*data.IC,1));ones(Int,(data.W*data.IC,1))*M]

wholearray=[];
arr=["p vlp min",m,n,nz,1,objnz]
push!(wholearray,arr)

for i=1:m
    for j=1:n
        if B[i,j]!=0
            push!(wholearray,("a",i,j,B[i,j]))
        end
    end
end

for i=1:length(P)
    if P[i]!=0
        push!(wholearray,("o",1,i,P[i]))
    end
end

for j=1:m
    if data.W<j<data.W+data.S+1
        push!(wholearray,("i",j,"s",RHS[j]))
    else
        push!(wholearray,("i",j,"u",RHS[j]))
    end
end


push!(wholearray,"e")

ins = open("/home/k2g00/k2g3475/mess/vlp/"*file[end-9:end-4]*".vlp","a")
# ins = open("C:/Users/AK121396/Desktop/mess.vlp","a")

writedlm(ins,wholearray)
close(ins)

for i=1:length(wholearray)
    println(wholearray[i])
end
wholearray
