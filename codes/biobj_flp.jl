using DelimitedFiles,DataFrames,JuMP,CPLEX,MathOptFormat,LinearAlgebra #all pkgs should be written like this
using Revise
using Pkg
cd(Pkg.devdir()*"/vOptPkg/")
using vOptPkg
# edit(pathof(vOptPkg))

dir1 = "/home/ak121396/Desktop/FLPInstances/"
dir = readdir(dir1)
dir2 = "/home/ak121396/Desktop/FLPInstances/flpLP/"
dir3 = "/home/ak121396/Desktop/FLPInstances/FLPvlp/"

# dir1 = "C:\\Users\\AK121396\\Desktop\\FLP\\size_5_10\\"
# dir2 = dir1*"varval\\PF\\"

# for s=4:15
# s = 1
subdir = dir1 * "size_20_40/"
files = readdir(subdir)
for idx=1:10
# idx = 10
    f = readdlm(subdir * files[idx])

    #facility_i, customer_j
    global i, j = filter!(a -> typeof(a) == Int, f[1, :])
    # filter!(a->typeof(a)==Int, f[2,:])
    fixcost = [];
    capa = [];
    for k = 2:i+1
        fix, cp = filter!(a -> typeof(a) == Int, f[k, :])
        push!(fixcost, fix)
        push!(capa, cp)
    end

    demand = filter!(a -> typeof(a) == Int, f[2+i, :])

    cost = [];
    cost2 = [];
    a, b = size(f)
    for k = i+3:a
        ct = filter!(a -> typeof(a) == Int, f[k, :])
        push!(cost2, ct)
        for l = 1:length(ct)
            push!(cost, ct[l])
        end
    end


    # vOpt MODEL
    biflp = vModel(solver = CplexSolver())
    @variables(biflp, begin
        x[1:i, 1:j], Bin
        y[1:i], Bin
    end)
    @constraint(biflp, con1[b = 1:j], sum(x[a, b] for a = 1:i) == 1)
    @constraint(biflp, con2[a = 1:i], sum(demand[b] * x[a, b] for b = 1:j) <= sum(capa[a] * y[a]))
    @addobjective(biflp, Min, sum(cost2[a][b] * x[a, b] for a = 1:i for b = 1:j))
    @addobjective(biflp, Min, sum(fixcost[a] * y[a] for a = 1:i))
    solve(biflp, method = :epsilon, step = 0.5) #dichotomy)
    Y_N = getY_N(biflp)
    obj1, obj2 = map(x -> x[1], Y_N), map(x -> x[2], Y_N)
    df = DataFrame(val1=obj1,val2=obj2)

    using CSV
    CSV.write("/home/ak121396/Desktop/complog.csv",df,append=true)
end



# JUMP MODEL. Cplex cannot solve tri-obj problems
triflp = Model(with_optimizer(CPLEX.Optimizer))
@variables(flp, begin
    x[1:i,1:j], Bin
    y[1:i] ,Bin
    z[1:j] ,Bin
    # k == 1
end)
@constraint(triflp, con1[b=1:j], sum(x[a,b] for a in 1:i) == z[b])
@constraint(triflp, con2[a=1:i,b=1:j], x[a,b]<=y[a])
@objective( triflp, Min, obj1, sum(cost2[a][b]*x[a,b] for a=1:i for b=1:j) )
@objective( triflp, Min, obj2, sum(fixcost[a]*y[a] for a=1:i) )
@objective( triflp, Min, obj3, -sum(demand[b]*z[b] for b=1:j) )



# printX_E(biflp)
# X = [vOptGeneric.getvalue(x,n) for n=1:length(Y_N)]
# for n = 1:length(Y_N)
#     X= vOptGeneric.getvalue(x,n)
#     print(findall(X.==1.0),"\n")
#     println("| z = ",Y_N[n])
# end
# using PyPlot
# clf()
# obj1, obj2 = map(x -> x[1], Y_N), map(x -> x[2], Y_N)
# PyPlot.title("Pareto Front")
# xlabel("obj1 to minimise");
# ylabel("obj2 to minimise");
# grid()
# plot(obj1, obj2, "o", markersize = "8",color="red", label = "\$Y_N\$")
# legend(loc = 1, fontsize = "small")
# show()
# gcf()
# savefig("/home/ak121396/Desktop/"*"fig_5")
