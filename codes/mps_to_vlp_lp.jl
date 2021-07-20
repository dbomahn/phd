using DelimitedFiles,MathProgBase,MathOptFormat,MathOptInterface,CPLEX,Random,LinearAlgebra,JuMP,Random
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
@show ARGS[1]
mpsfile = "/home/ak121396/Desktop/instances/MIPLIB(official)/mps/cvs16r128-89.mps"

mpsmodel = CPLEX.CplexMathProgModel()
MPB.loadproblem!(mpsmodel,mpsfile)
C = MPB.getobj(mpsmodel)
B = MPB.getconstrmatrix(mpsmodel)
m,n=size(B)
lb = MPB.getconstrLB(mpsmodel)
ub = MPB.getconstrUB(mpsmodel)
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

###################################### make Coefficient matrix for obj2&3
uC = unique(C)
if length(uC)==1
    if uC[1]>=0
        coef = rand(1:10,n)
        #  Coefficient Permutation for the 2nd and 3rd obj
        coef2 = shuffle(coef); coef3 = shuffle(coef);
        P =[coef,coef2,coef3]
    else #negative Coefficient
        coef = rand(-10:-1,n)
        coef2 = shuffle(coef); coef3 = shuffle(coef);
        P =[coef,coef2,coef3]
    end
else
    coef = [(C[x]) for x=1:length(C)] #Int(C)
    coef2 = shuffle(coef); coef3 = shuffle(coef);
    P =[coef,coef2,coef3]
end

# P2 = [coef;coef2;coef3]
# P3 = reshape(P2,3,n)
# sP = sparse(P3)
# JLD2.@save "/home/ak121396/Desktop/instances/MIPLIB(official)/coef.jld2" sP

obj = 3
objnz = n*obj-(count(iszero,coef)+count(iszero,coef2)+count(iszero,coef3))
NZ = length(B)-count(iszero,B)
# ######################### vlp file for BENSOLVE ####################
ins = open(mpsfile[1:end-4]*".vlp","a")
wholearray=[]; #NZ
arr=["p vlp min",m,n,NZ,obj,objnz]
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

for i=1:obj
   for j=1:n
       if P[i][j]!=0
           push!(wholearray,("o",i,j,P[i][j]))
       end
   end
end

for i=1:m
   push!(wholearray,("i",i,signs[i],RHS[i]))
end

for j=1:n
   push!(wholearray,("j", j,'d',0,1))
end
push!(wholearray,"e")
writedlm(ins,wholearray)
close(ins)

########################  Copy mps file to JuMP model
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
writemodel(ks_model,fpfile[1:end-4]*"_ks.lp")

# FPBH LP file format
# mof_model = MathOptFormat.MPS.Model();
# MOI.read_from_file(mof_model, mpsfile);
# MOI.copy_to(backend(fpbh_model), mof_model);
# allvar2 = all_variables(fpbh_model)
fpbh_model = JuMP.Model();
@variable(fpbh_model, x[1:n], Bin)#data.n
@objective(fpbh_model, Min, dot(C[1],x))
@constraint(fpbh_model, obj2, dot(C[2],x) == 0 );
@constraint(fpbh_model, obj3, dot(C[3],x) == 0 );
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


############################  coef matrix from vlpfile   ###################
vlp = readdlm("/home/ak121396/Desktop/instances/MIPLIB(official)/vlp/cvs16r128-89.vlp");
os = findall(i->vlp[i,1]=="o", 1:length(vlp[:,1]))
inio = os[1]; endo = os[end]; Pmtx = vlp[inio:endo,2:4]; C = zeros(3,n)
for k=1:length(Pmtx[:,1]) #(len,id) enumerate
    x = Pmtx[k,:][1]; y = Pmtx[k,:][2];
    C[x,y] = Pmtx[k,:][3]
end

fpbh_model = JuMP.Model();
@variable(fpbh_model, x[1:n], Bin)#data.n
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
writemodel(fpbh_model,mpsfile[1:end-20]*"_fp.lp")
