using DelimitedFiles,MathProgBase,MathOptFormat,MathOptInterface,CPLEX,Random,LinearAlgebra,JuMP
const MPB = MathProgBase; const MOF = MathOptFormat; const MOI = MathOptInterface;
function writemodel(model::Model, filename::String)
     if endswith(lowercase(filename), ".lp")
         file = MathOptFormat.LP.Model()
     elseif endswith(lowercase(filename), ".mps")
         file = MathOptFormat.MPS.Model()
     else
         println("Unknown file extension to write model: ",
                 split(filename, ".")[end])
         exit(8)
     end
     MOI.copy_to(file, backend(model))
     MOI.write_to_file(file, filename)
     println("wrote model to ", filename)
end
function loadlp(filename,solver=CplexSolver(CPX_PARAM_SCRIND=0))
    model=buildlp([-1,0],[2 1],'<',1.5, solver) # create dummy model with correct solver
    MPB.loadproblem!(model,filename) # load what we actually want
    return model
end
@show ARGS[1]

##################### FPBH.lp file to Kirlik.lp file ######################

file = ARGS[1];
lpmodel = loadlp(file)
C0 = MPB.getobj(lpmodel); B0 = MPB.getconstrmatrix(lpmodel)
B = B0[1:end-2,:]; m,n=size(B);
lb = MPB.getconstrLB(lpmodel)[1:end-2]; ub = MPB.getconstrUB(lpmodel)[1:end-2]
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
        push!(signs,"s")
    end
end
C = [C0,B0[end-1,:],B0[end,:]]

ks_model = JuMP.Model();
@variable(ks_model, x[1:n], Bin)#data.n
@variable(ks_model, dummy==1);
MOI.set(ks_model, MOI.ObjectiveFunction{MOI.SingleVariable}(),MOI.SingleVariable(dummy));
@constraint(ks_model, o1, dot(C[1],x) <= 9999999999999999 );
@constraint(ks_model, o2, dot(C[2],x) <= 9999999999999999 );
@constraint(ks_model, o3, dot(C[3],x) <= 9999999999999999 );  
for k=1:m
    if signs[k] == "l"
        @constraint(ks_model, dot(B[k,:],x) >= RHS[k])
    elseif signs[k] == "u"
        @constraint(ks_model, dot(B[k,:],x) <= RHS[k])
    else
        @constraint(ks_model, dot(B[k,:],x) == RHS[k])
    end
end
writemodel(ks_model,file[1:end-6]*"_ks.lp")

# removing _ in variables : unnecessary
# dt0 = readlines(file[1:end-6]*"_ks.lp")
# dt1 = lowercase.(dt0)
# dt2 = replace.(dt1, ['_','[',']'] =>"")
# writedlm(file[1:end-6]*"_ks.lp", dt2)
