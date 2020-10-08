using DataStructures,DelimitedFiles,DataFrames,CSV,StatsBase,Random,LinearAlgebra,CPUTime,JuMP,CPLEX
struct path{S<:String}
    dir1::S; dir2::S;
end

mutable struct Data
    dvar::String; dtfile::String
    n::Int; C::Array{}
    ub::Int; weight::Array{}
    L::Array{}
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
        L = round.(readdlm(dvar),digits=4)
        new(dvar,dtfile,n,C,ub,weight,L)
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
function solveKP(minid,xfix,n) #solveLP ,fmin,fmax
    # yub = [idx[j]*steps[j]+fmin[j] for j=1:3]
    # ylb = [(idx-[1,1,1])[j]*steps[j]+fmin[j] for j=1:3]
    kp_m = Model(CPLEX.Optimizer)
    MOI.set(kp_m, MOI.RawParameter("CPX_PARAM_SCRIND"), false)
    MOI.set(kp_m, MOI.RawParameter("CPX_PARAM_THREADS"),1  )
    @variables(kp_m, begin
        0 <= x[1:n] <=1
    end)
    @variable(kp_m, y==1)
    @variable(kp_m, x[minid] == xfix)
    @objective( kp_m, Min, y )
    @constraint( kp_m, -sum(weight[i]*x[i] for i=1:n) >= -ub )
    # @constraint( kp_m, c1, -dot(C[1,:],x[:]) <= yub[1]-1 )#ylb[1]+1 <=
    # @constraint( kp_m, c2, -dot(C[2,:],x[:]) <= yub[2]-1 ) #ylb[2]+1 <=
    # @constraint( kp_m, c3, -dot(C[3,:],x[:]) <= yub[3]-1 ) #ylb[3]+1 <=
    optimize!(kp_m)
    if termination_status(kp_m) == MOI.OPTIMAL
        return JuMP.value.(x)
    else
        return 0; #print("no new lp solution found");
    end
end
pathes = ("E:\\Bensolve_KP\\X\\","E:\\Bensolve_KP\\data\\")
# pathes = ("/home/ak121396/Desktop/BENoutputs/KP/X/", "/home/ak121396/Desktop/BENKP/data/")
file = readdir(pathes[1]); ins = readdir(pathes[2])
i=13
kp = Data(pathes[1]*file[i],pathes[2]*ins[i])
function fractionality(n,L)
    f_x = []
    for i=1:10
        push!(f_x, sum( [ abs(L[:,i][k] - floor(L[:,i][k]+0.5) ) for k=1:size(L)[1]]) )
    end
    return f_x
end

frac_x = fractionality(kp.n, kp.L)
fr = filter!(x->x!=0,frac_x)
# I = findall(x->(x!=0) && (x!=1), frac_x)




archive = []

for i=1:size(kp.L)[1]
c_max = Inf; ctx = 0; iter = 0

while isempty(I)==false && ctx < c_max
iter+=1
minid = findall(x->x==minimum(fr),fr)
xfix = round(kp.L[i,minid][1])

solveKP(minid,xfix,kp.n)
cpx = copy(x[i,:])
cpx[I[1]] = floor(Int,x[i,:][I[1]])
if KPfbcheck(txi,n,weight,ub)== true
    xbar = cpx
else
    break
end

# since all coefficients are nonnegative (trivially roundable), round down
