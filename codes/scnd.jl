using DelimitedFiles,DataFrames,JuMP,CPLEX#CPUTime, #DataStructures,
fpath = "/home/ak121396/Desktop/instances/SCND/"
fpath = "F:/SCND/"
dt = readdlm(fpath*"Test1S1", '\t')

notafile = readdlm(fpath*"Notations.txt", '=')
nota = notafile[1:end,1]
N= Dict()
for i=1:length(nota)-1
    id1 = findall(x->x==nota[i], dt)[1][1]
    id2 = findall(x->x==nota[i+1], dt)[1][1]
    if id2-id1 <3
        tmp = filter(x->x!="",  dt[id1+(id2-id1-1),:])
        N[nota[i]] = tmp
    else
        tmp = [filter(x->x!="", dt[x,:]) for x in id1+1:id1+(id2-id1-1)]
        N[nota[i]] = tmp
    end
end

##########################  Mathematical model  #########################
# scnd = Model(with_optimizer(GLPK.Optimizer))
scnd = Model(CPLEX.Optimizer)
# @variable(scnd, x[1:N["supplier"][1], 1:N["plant"][1],1:N["mode"][1]])
@variable(scnd, y[1:N["plant"][1]+N["distribution"][1]] )
@variable(scnd, u[1:1:N["supplier"][1],1:N["plant"][1]] )
@variable(scnd, xij[1:N["supplier"][1],1:N["plant"][1]] )
@variable(scnd, xjk[1:N["plant"][1],1:N["distribution"][1]] )
@variable(scnd, xkl[1:N["distribution"][1],1:N["customer"][1]] )
@variable(scnd, h[1:N["plant"][1]+N["distribution"][1]] )


@constraint(scnd, sum(xij[i,m]) == sum(xjk[k,m]))
@constraint(scnd, dot(data.B[k,:],x) <= dot(xij[1:N[]]))
@constraint(scnd, dot(data.B[k,:],x) == data.RHS[k])
optimize!(scnd);

N["cep"][1]
[1][1]
N["customer"]
