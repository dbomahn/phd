using SparseArrays,DelimitedFiles,JuMP,MathOptInterface
cd("/Users/AK121396/Desktop/ProjectBenders/")
using MathOptInterface,CPLEX,JuMP,JLD2,DelimitedFiles
struct Valu
    x::String
    y::String
    dvar::Array{}
    LB::Array{}
    LBmtx::Array{}
    function Valu(x, y)
        JLD2.@load x dv
        dv0 = Array(dv)
        # dv0 = readdlm(x)
        dv1 = round.(dv0; digits = 4)
        objs = round.(readdlm(y); digits = 4)
        ind = findall(i -> 0 in objs[i, :], 1:size(objs)[1])
        dv2 = dv1[setdiff(1:end, ind), :]
        LBmtx = objs[setdiff(1:end, ind), 2:end]
        dvar = [Vector(dv2[i, :]) for i = 1:size(LBmtx)[1]]
        LB = [Vector(LBmtx[i, :]) for i = 1:size(LBmtx)[1]]
        new(x, y, dvar, LB, LBmtx)
    end
end
m1 = read_from_file("F:scnd/Test1S21dim.lp")
set_optimizer(m1, CPLEX.Optimizer)
allvar = all_variables(m1)
bvar = findall(i->i==1, is_binary.(allvar))
relax_integrality(m1)
print(m1)

pr = Valu("F:scnd/test01S2_X.jld2","F:scnd/test01S2_img_p.sol");
findall(i->0<i<1, pr.dvar[1])
for k in bvar
    println(pr.dvar[1][k])
    # JuMP.fix(allvar[k], pr.dvar[1][k]; force = true)
end

# py = pr.dvar[1][length(rvar)+1:length(rvar)+length(y1)]
# pui = pr.dvar[1][length(rvar)+1+length(y1):length(rvar)+length(y1)+length(uij1)]
# puj = pr.dvar[1][length(rvar)+1+length(y1)+length(uij1):length(rvar)+length(y1)+length(uij1)+length(ujk1)]
# puk = pr.dvar[1][length(rvar)+1+length(y1)+length(uij1)+length(ujk1):end]
# for k=1:length(y1)
#     JuMP.fix(y1[k], py[k]; force = true)
# end
# for k=1:length(uij1)
#     JuMP.fix(uij1[k], pui[k]; force = true)
# end
# for k=1:length(ujk1)
#     JuMP.fix(ujk1[k], puj[k]; force = true)
# end
# for k=1:length(ukl1)
#     JuMP.fix(ukl1[k], puk[k]; force = true)
# end
#
# optimize!(m1); termination_status(m1)




write_to_file(scnd,"/home/ak121396/Desktop/relise/noobj.lp")

model = read_from_file("/home/ak121396/Desktop/relise/noobj.lp")
1
# set_optimizer(model, CPLEX.Optimizer)
# optimize!(model)
# objective_value(model)
objective_function(model)
coefficient(model,h1)
obj = MOI.get(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
obj.terms
length(obj.terms)
c = zeros(length(obj.terms))
for term in obj.terms
    c[term.variable_index.value] = term.coefficient
end
@show(c)

#constraints
list_of_constraint_types(model)
con = Dict()
con_rhs = Dict()
idx_con = Dict()
con_idx = Dict()

constraint = all_constraints(model, GenericAffExpr{Float64,VariableRef}, MathOptInterface.LessThan{Float64})
for i in 1:length(constraint)
    con_name = name(constraint[i])
    con_rhs[con_name] = constraint_object(constraint[i]).set.upper
    con[con_name]= :Less
    idx_con[constraint[i].index.value] = con_name
    con_idx[con_name] = constraint[i].index.value
end
constraint = all_constraints(model, GenericAffExpr{Float64,VariableRef}, MathOptInterface.GreaterThan{Float64})
for i in 1:length(constraint)
    con_name = name(constraint[i])
    con_rhs[con_name] = constraint_object(constraint[i]).set.lower
    con[con_name]= :Greater
    idx_con[constraint[i].index.value] = con_name
    con_idx[con_name] = constraint[i].index.value
end

constraint = all_constraints(model, GenericAffExpr{Float64,VariableRef}, MathOptInterface.EqualTo{Float64})
for i in 1:length(constraint)
    con_name = name(constraint[i])
    con_rhs[con_name] = constraint_object(constraint[i]).set.value
    con[con_name]= :Equal
    idx_con[constraint[i].index.value] = con_name
    con_idx[con_name] = constraint[i].index.value
end
# variable
var = Dict()
var_idx = Dict()
idx_var = Dict()
var_lb = Dict()
var_ub = Dict()
var_ref = Dict()
for var_name in (all_variables(model))
    var[var_name] = :Con
    var_idx[var_name] = var_name.index.value
    idx_var[var_name.index.value] = var_name
    var_lb[var_name] = -Inf
    var_ub[var_name] = Inf
    var_ref[var_name] = :except
end


variable = all_constraints(model, VariableRef, MathOptInterface.Interval{Float64})
for i in (1:length(variable))
    var_name = constraint_object(variable[i]).func
    var_lb[var_name] = constraint_object(variable[i]).set.lower
    var_ub[var_name] = constraint_object(variable[i]).set.upper
    var_ref[var_name] = :Interval
end


variable = all_constraints(model, VariableRef, MathOptInterface.EqualTo{Float64})
for i in (1:length(variable))
    var_name = constraint_object(variable[i]).func
    var_lb[var_name] = constraint_object(variable[i]).set.value
    var_ub[var_name] = constraint_object(variable[i]).set.value
    var_ref[var_name] = :EqualTo
end

variable = all_constraints(model, VariableRef, MathOptInterface.Integer)
for i in (1:length(variable))
    var_name = constraint_object(variable[i]).func
    var[var_name] = :Int
end

variable = all_constraints(model, VariableRef, MathOptInterface.ZeroOne)
for i in (1:length(variable))
    var_name = constraint_object(variable[i]).func
    var[var_name] = :Bin
    var_lb[var_name] = 0
    var_ub[var_name] = 1
end

# Sparse matrix
II = Int[]
JJ = Int[]
VV = Float64[]
uu = Dict()
ll = Dict()


con_set = [k for (k,v) in con if v==:Less]
for i in (con_set)
    con_term =collect(linear_terms(constraint_object(constraint_by_name(model, i )).func))
    uu[con_idx[i]] = con_rhs[i]
    ll[con_idx[i]] = -Inf
    for j in 1:length(con_term)
        push!(II, con_idx[i])
        push!(JJ, var_idx[con_term[j][2]])
        push!(VV, con_term[j][1])
    end
end

con_set = [k for (k,v) in con if v==:Greater]
for i in (con_set)
    uu[con_idx[i]] = -(con_rhs[i])
    ll[con_idx[i]] = -Inf
    con_term =collect(linear_terms(constraint_object(constraint_by_name(model, i )).func))
    for j in 1:length(con_term)
        push!(II, con_idx[i])
        push!(JJ, var_idx[con_term[j][2]])
        push!(VV, -(con_term[j][1]))
    end
end
con_set = [k for (k,v) in con if v==:Equal]
for i in (con_set)
    uu[con_idx[i]] = con_rhs[i]
    ll[con_idx[i]] = con_rhs[i]
    con_term =collect(linear_terms(constraint_object(constraint_by_name(model, i )).func))
    for j in 1:length(con_term)
        push!(II, con_idx[i])
        push!(JJ, var_idx[con_term[j][2]])
        push!(VV, con_term[j][1])
    end
end

A = sparse(II,JJ,VV)

con # constraint with
con_rhs
idx_con
con_idx

var # variable with type
var_idx
idx_var
var_lb
var_ub
var_ref

A # sparse matrix for constraint (standard form)
uu # constraint upper bound
ll # constraint lower bound
idx_con[1]
findall(i->idx_con[i]=="obj1",1:length(idx_con))
idx_con[2]
