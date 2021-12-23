using DataFrames,DelimitedFiles,JuMP,LinearAlgebra,CPLEX,MathProgBase
const MPB = MathProgBase

function loadlp(filename,solver=CplexSolver(CPX_PARAM_SCRIND=0))
    model=buildlp([-1,0],[2 1],'<',1.5, solver) # create dummy model with correct solver
    MPB.loadproblem!(model,filename) # load what we actually want
    return model
end
@show file = ARGS[1]

# file = "/home/ak121396/Desktop/instances/PublicInstances/toy.dzn"
# data = readdlm(dir*ins[1], ';')
# file = "C:/Users/AK121396/Desktop/phd/toy.dzn"
# data = readdlm(dir, ';')

mutable struct Data
    filepath::String; W::Int; S::Int; C::Array{}; F::Array{}; G::Array{};
    TC::Array{}; SC::Array{}; IC::Int; IP::Array{}

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
        TC = copy(SC)
        SC = reshape(transpose(SC), (1,S*W))

        tm = split(data[size(data)[1]],r"[|]")
        IP = []
        for i=2:1+IC
            tm2 = parse.(Int,split(tm[i], ","))
            push!(IP,tm2)
        end

        new(filepath,W,S,C,F,G,TC,SC,IC,IP)
    end
end

data = Data(file)
M = 9999999999999999
###########################   Mathematical Model  #############################
wlp = Model(CPLEX.Optimizer) #GLPK.Optimizer
@variable( wlp, 0<=x[1:data.S,1:data.W])
@variable( wlp, 0<=y[1:data.W]<=1 ) #,binary=true
@variable( wlp, 0<=z[1:data.S,1:data.W]<=1)
for w=1:data.W
    @constraint( wlp, sum(x[s,w] for s=1:data.S) <= data.C[w] )
end
for s=1:data.S
    @constraint( wlp, sum(x[s,w] for w=1:data.W) == data.G[s])
end
for w=1:data.W
    @constraint( wlp, sum(x[s,w] for s=1:data.S) <= (M*y[w]) )
end
for w=1:data.W
    for i=1:data.IC
        @constraint( wlp, x[data.IP[i][1],w] <= M*z[data.IP[i][1],w] )
        @constraint( wlp, x[data.IP[i][2],w] <= M*(1-z[data.IP[i][1],w]) )
    end
end
@objective( wlp, Min, sum(sum(data.TC[s,w]*x[s,w] for w=1:data.W) for s=1:data.S) + sum(data.F[w]*y[w] for w=1:data.W) )
print(wlp)
write_to_file(wlp , "/home/ak121396/Desktop/wlptest.lp")

# optimize!(wlp)
# MOI.get(wlp, MOI.ObjectiveValue())
lpmodel = loadlp("/home/ak121396/Desktop/wlptest.lp")

P = MPB.getobj(lpmodel)
B = MPB.getconstrmatrix(lpmodel)
m,n=size(B)
lb = MPB.getconstrLB(lpmodel)
ub = MPB.getconstrUB(lpmodel)
RHS = Dict()
for i=1:m
    if ub[i]==Inf
        RHS[i] = lb[i]
    else
        RHS[i] = ub[i]
    end
end
signs = []
for i=1:m
    if ub[i] == Inf
        push!(signs,"l")
    elseif lb[i] == -Inf
        push!(signs,"u")
    else
        push!(signs, "s")
    end
end
nz = count(i->(i!=0),B)
objnz = count(i->(i!=0),P)
obj=1
wholearray=[];
arr=["p vlp min",m,n,nz,obj,objnz]
push!(wholearray,arr)

for i=1:m
   for j=1:n
       if (B[i,j]!=0)
           if (B[i,j]%1) == 0 #if B[i,j] is Int
               push!(wholearray,("a",i,j,Int128(B[i,j])))
           else# B[i,j] is Float
               push!(wholearray,("a",i,j,Float64(B[i,j])))
           end
       end
   end
end

for i=1:length(P)
    if P[i]!=0
        push!(wholearray,("o",1,i,P[i]))
    end
end

for i=1:m
   push!(wholearray,("i",i,signs[i],RHS[i]))
end

for j=1:n
    if j<=data.W*data.S
        push!(wholearray,("j",j,"l",0))
    else
        push!(wholearray,("j", j,'d',0,1))
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





################################  convert to vlp file  ###############################################
# RHS = [data.C;data.G;zeros(Int,(data.W+data.W*data.IC,1));ones(Int,(data.W*data.IC,1))*M]
# A1 = zeros(Int,data.W,data.S*data.W+data.W+data.W*data.S)
# for w=1:data.W
#     for s=1:data.S
#         A1[w,data.W*(s-1)+w] = 1
#     end
# end
#
# A2 = zeros(Int,data.S,data.S*data.W+data.W+data.W*data.S)
# for s=1:data.S
#     for w=1:data.W
#         A2[s,w+data.W*(s-1)] = 1
#     end
# end
#
# A3 = copy(A1)
# for w=1:data.W
#     A3[w,data.S*data.W+w] = -M
# end
# # A32 = Matrix{Int}(I, data.W, data.W)*(-M)
#
# A41= zeros(Int,data.IC*data.W,data.S*data.W+data.W)
# A42 = zeros(Int,data.IC*data.W,data.S*data.W)
# A51 = copy(A41); A52 =copy(A42)
# for i=1:data.IC
#     for w=1:data.W
#         k = data.W*(i-1)+w
#         A41[k, data.W*(data.IP[i][1]-1)+w] = 1
#         A42[k, data.W*(data.IP[i][1]-1)+w] = -M
#     end
#     for w=1:data.W
#         k = data.W*(i-1)+w
#         A51[k, data.W*(data.IP[i][2]-1)+w] = 1
#         A52[k, data.W*(data.IP[i][2]-1)+w] = M
#     end
# end
#
# A4 = hcat(A41, A42)
# A5 = hcat(A51, A52)
# P = hcat(data.SC,transpose(data.F), zeros(Int, (1,data.S*data.W)) )
#
# B = [A1;A2;A3;A4;A5]
# for i=1:m
#     if data.W<i<data.W+data.S+1
#         push!(wholearray,("i",i,"s",RHS[i]))
#     else
#         push!(wholearray,("i",i,"u",RHS[i]))
#     end
# end
#
# for j=1:n
#     if j<data.W*data.S+1
#         push!(wholearray,("j",j,"l",0))
#     else
#         push!(wholearray,("j",j,"d",0,1))
#     end
# end
