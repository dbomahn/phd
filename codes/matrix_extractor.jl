using SparseArrays,DelimitedFiles,JuMP,MathOptInterface

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
