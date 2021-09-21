using DelimitedFiles,DataFrames,JuMP#,SparseArrays,CPLEX#CPUTime, #DataStructures,
fpath = "/home/ak121396/Desktop/instances/SCND/"
fpath = "F:/scnd/"
dt = readdlm(fpath*"Test1S1") #, '\t'

notafile = readdlm(fpath*"Notations.txt", '=')
nota = notafile[1:end,1]
N= Dict()



i=11
id1 = findall(x->x==nota[i], dt)[1][1]
id2 = findall(x->x==nota[i+1], dt)[1][1]
id2-id1

W = []
for x=id1+1:id1+(id2-id1-1)
    tmp = filter(x->x!="", dt[x,:])[1:end]
    println(tmp)
    push!(W,tmp)
end
W
filter(x->x!="", dt[id1+1:id1+(id2-id1-1),:])
dt[id1+1:id1+(id2-id1-1),:]


for i=1:length(nota)-1
    id1 = findall(x->x==nota[i], dt)[1][1]
    id2 = findall(x->x==nota[i+1], dt)[1][1]
    if id2-id1<3
        tmp = filter(x->x!="",  dt[id1+(id2-id1-1),:])
        if length(tmp)<2
            N[nota[i]] = tmp[1]
        else
            N[nota[i]] = tmp
        end
    else
        W = []
        for x=id1+1:id1+(id2-id1-1)
            tmp = filter(x->x!="", dt[x,:])
            push!(W,tmp)
        end
        # tmp = [filter(x->x!="", dt[x,:]) for x in id1+1:id1+(id2-id1-1)]
        N[nota[i]] = W
    end
end
d = N["demand"]
c = append!(N["fcp"],N["fcd"])
a = N["vcs"]
e = append!(N["vcp"],N["vcd"])
cap = append!(N["cas"],N["cap"],N["cad"])
b = reshape( N["ves"], (N["supplier"],Int(length(N["ves"])/N["supplier"])) )
q = append!(N["vep"],N["ved"])

Mij = transpose(reshape(N["ModeIJ"], (N["plant"],N["supplier"])))
Mjk = transpose(reshape(N["ModeJK"], (N["distribution"],N["plant"])))
Mkl = transpose(reshape(N["ModeKL"], (N["customer"],N["distribution"])))

sum(length(N["fixedcostModepd"][i]) for i=1:length(N["fixedcostModepd"]))

reshape(N["fixedcostModepd"], (1,sum(N["ModeJK"])) )

N["fixedcostModepd"]
N["fcd"][2][1]
sum(Mjk[j,1:k])

Mjk
N["tcd"][2]
Mjk
N["fixedcostModepd"]
[1]

narc = 1:N["supplier"]*N["plant"]+N["plant"]*N["distribution"]+N["distribution"]*N["customer"]
cou = 0
for i=1:length(N["LcapacityModepd"])
    cou = cou+length(N["LcapacityModepd"][i])
    # if length(N["LcapacityModedc"][i])== length(N["fixedcostModedc"][i])
        println(length(N["LcapacityModepd"][i])) # yes!)
    # end
end

cou
N["tcp"]
N["supplier"][1]
i=12
id1 = findall(x->x==nota[i], dt)[1][1]
id2 = findall(x->x==nota[i+1], dt)[1][1]

filter(x->x!="",  dt[id1+(id2-id1-1),:])


##########################  Mathematical model  #########################
scnd = Model()#CPLEX.Optimizer
# @variable(scnd, x[1:N["supplier"][1], 1:N["plant"][1],1:N["mode"][1]])
@variable(scnd, y[1:N["plant"]+N["distribution"],1:2] )
@variable(scnd, u[1:narc,] )
@variable(scnd, xij[1:N["supplier"]*N["plant"]+N["plant"]*N["distribution"]+N["distribution"]*N["customer"],] )
@variable(scnd, h[1:N["plant"]+N["distribution"]] )


@constraint(scnd, sum(xij[i,m]) == sum(xjk[k,m]))
@constraint(scnd, dot(data.B[k,:],x) <= dot(xij[1:N[]]))
@constraint(scnd, dot(data.B[k,:],x) == data.RHS[k])
optimize!(scnd);

6*6*12*60

N["tcp"][1]

c[j][t]*y[j,t]
i=6
id1 = findall(x->x==nota[i], dt)[1][1]
id2 = findall(x->x==nota[i+1], dt)
[1][1]
if id2-id1 <3
