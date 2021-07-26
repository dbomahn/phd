using DelimitedFiles,MathProgBase,MathOptFormat,MathOptInterface,CPLEX,Random,LinearAlgebra,JuMP #Random
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

######################## Recall Coefficient matrix from vlpfile  ###############
vlp = readdlm(ARGS[2]);
os = findall(i->vlp[i,1]=="o", 1:length(vlp[:,1]))
inio = os[1]; endo = os[end]; Pmtx = vlp[inio:endo,2:4]; C = zeros(3,n)
for k=1:length(Pmtx[:,1]) #(len,id) enumerate
    x = Pmtx[k,:][1]; y = Pmtx[k,:][2];
    C[x,y] = Pmtx[k,:][3]
end

# CPLEX multi-objective LP file format
cplex_model = JuMP.Model();
@variable(cplex_model, x[1:n], Bin)
@objective(cplex_model, obj1, dot(C[1,:],x))
@constraint(cplex_model, obj2, dot(C[2,:],x)  );
@constraint(cplex_model, obj3, dot(C[3,:],x)  );
for k=1:m
    if signs[k] == "l"
        @constraint(cplex_model, dot(B[k,:],x) >= RHS[k])
    elseif signs[k] == "u"
        @constraint(cplex_model, dot(B[k,:],x) <= RHS[k])
    else
        @constraint(cplex_model, dot(B[k,:],x) == RHS[k])
    end
end
writemodel(cplex_model,mpsfile[1:end-4]*"ibm.lp")


# FPBH LP file format
fpbh_model = JuMP.Model();
@variable(fpbh_model, x[1:n], Bin)
@objective(fpbh_model, Min, dot(C[1,:],x))
@constraint(fpbh_model, obj2, dot(C[2,:],x) == 0 );
@constraint(fpbh_model, obj3, dot(C[3,:],x) == 0 );
for k=1:m
    if signs[k] == "l"
        @constraint(fpbh_model, dot(B[k,:],x) >= RHS[k])
    elseif signs[k] == "u"
        @constraint(fpbh_model, dot(B[k,:],x) <= RHS[k])
    else
        @constraint(fpbh_model, dot(B[k,:],x) == RHS[k])
    end
end
writemodel(fpbh_model,mpsfile[1:end-4]*"_fp.lp")

########################  Copy mps file to JuMP model
# ks_model = JuMP.Model();
# @variable(ks_model, x[1:n], Bin)
# @variable(ks_model, dummy==1);
# MOI.set(ks_model, MOI.ObjectiveFunction{MOI.SingleVariable}(),MOI.SingleVariable(dummy));
# @constraint(ks_model, obj1, dot(C[1,:],x) <= 9999999999999999 );
# @constraint(ks_model, obj2, dot(C[2,:],x) <= 9999999999999999 );
# @constraint(ks_model, obj3, dot(C[3,:],x) <= 9999999999999999 );
# for k=1:m
#     if signs[k] == "l"
#         @constraint(ks_model, dot(B[k,:],x) >= RHS[k])
#     elseif signs[k] == "u"
#         @constraint(ks_model, dot(B[k,:],x) <= RHS[k])
#     else
#         @constraint(ks_model, dot(B[k,:],x) == RHS[k])
#     end
# end
# writemodel(ks_model,mpsfile[1:end-4]*"_ks.lp")
