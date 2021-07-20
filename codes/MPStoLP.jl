using DelimitedFiles,MathProgBase, MathOptInterface,MathOptFormat,CPLEX,Random,LinearAlgebra,JuMP
const MPB = MathProgBase;
#const MOI=MathOptInterface; const MOF=MathOptFormat

################################ convert .mps to .lp ###########################
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




lpath = "/home/ak121396/Desktop/instances//MIPLIB(official)/LP/"
files = readdir(lpath)
mpspath = "/home/ak121396/Desktop/instances/MIPLIB(official)/mps/"
mpss = readdir(mpspath)

for i =1:length(files)
    # print("LP: ",files[i], "\n mps: ", mpss[i],"\n")

    lpmodel = loadlp(lpath*files[i])
    C = MPB.getobj(lpmodel)

    #Copy mps file to JuMP model
    ks_model = Model(); fpbh_model = Model();
    mof_model = MOF.MPS.Model();
    MOI.read_from_file(mof_model, mpspath*mpss[i]);

    # Kirlik solver LP file format
    MOI.copy_to(backend(ks_model), mof_model);
    allvar = all_variables(ks_model)
    @constraint(ks_model, o1, dot(P[1],allvar) <= 1 );
    @constraint(ks_model, o2, dot(P[2],allvar) <= 1 );
    @constraint(ks_model, o3, dot(P[3],allvar) <= 1 );
    @variable(ks_model, dummy==1);
    MOI.set(ks_model, MOI.ObjectiveFunction{MOI.SingleVariable}(),MOI.SingleVariable(dummy));
    writemodel(ks_model,"/home/ak121396/Desktop/instances/MIPLIB(official)/kirlik/"*files[i][1:end-3]*".lp")

    # FPBH LP file format
    MOI.copy_to(backend(fpbh_model), mof_model);
    @constraint(fpbh_model, obj2, dot(P[2],allvar) == 0 );
    @constraint(fpbh_model, obj3, dot(P[3],allvar) == 0 );
    writemodel(fpbh_model,"/home/ak121396/Desktop/instances/MIPLIB(official)/fpbh/"*files[i][1:end-3]*".lp")
end

i=8

mps = CPLEX.CplexMathProgModel()
MPB.loadproblem!(mps,mpspath*mpss[i])
C = MPB.getobj(mps)
uC = unique(C)
if length(uC)==1
    if uC[1]>=0
        coef = rand(1:10,n)
        #  Coefficient Permutation for the 2nd and 3rd obj
        coef2 = shuffle(coef); coef3 = shuffle(coef);
        P =[coef,coef2,coef3]
    else #negative Coefficient
        coef = rand(-10,1,n)
        coef2 = shuffle(coef); coef3 = shuffle(coef);
        P =[coef,coef2,coef3]
else
    coef = [(C[x]) for x=1:length(C)] #Int(C)
    coef2 = shuffle(coef); coef3 = shuffle(coef);
    P =[coef,coef2,coef3]
end

#Copy mps file to JuMP model
ks_model = JuMP.Model(); fpbh_model = JuMP.Model();
mof_model = MathOptFormat.MPS.Model();
MOI.read_from_file(mof_model, mpspath*mpss[i]);

# Kirlik solver LP file format
MOI.copy_to(backend(ks_model), mof_model);
allvar = all_variables(ks_model)
@constraint(ks_model, o1, dot(P[1],allvar) <= 1 );
@constraint(ks_model, o2, dot(P[2],allvar) <= 1 );
@constraint(ks_model, o3, dot(P[3],allvar) <= 1 );
@variable(ks_model, dummy==1);
MOI.set(ks_model, MOI.ObjectiveFunction{MOI.SingleVariable}(),MOI.SingleVariable(dummy));
writemodel(ks_model,"/home/ak121396/Desktop/instances/MIPLIB(official)/kirlik/"*files[i][1:end-3]*".lp")

# FPBH LP file format
MOI.copy_to(backend(fpbh_model), mof_model);
@constraint(fpbh_model, obj2, dot(P[2],allvar) == 0 );
@constraint(fpbh_model, obj3, dot(P[3],allvar) == 0 );
writemodel(fpbh_model,"/home/ak121396/Desktop/instances/MIPLIB(official)/fpbh/"*files[i][1:end-3]*".lp")
end

##################Load



jump_model = Model()
mof_model = MathOptFormat.MPS.Model()
MOI.read_from_file(mof_model, mpss*mpss[31])
MOI.copy_to(backend(jump_model), mof_model)
writemodel("/home/ak121396/Desktop/instances/MIPLIB(official)/test"*mpss[30]*".lp")
