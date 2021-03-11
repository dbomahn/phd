using JuMP,MathOptFormat,CPLEX,MathProgBase
const MPB = MathProgBase

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

cd("../../instances/MIPLIB(official)/")

flist = readdir()
for i =2:length(flist)
    jump_model = Model()
    mof_model = MathOptFormat.MPS.Model()
    # "C:\\Users\\AK121396\\Downloads\\"
    MOI.read_from_file(mof_model, flist[i])
    MOI.copy_to(backend(jump_model), mof_model)
    # optimize!(jump_model, with_optimizer=CPLEX.Optimizer)
    # print(jump_model)
    writemodel(jump_model,pwd()*"/LP/"*flist[i][1:end-3]*"lp")
end

##############################convert .lp to .vlp ################################
function loadlp(filename,solver=CplexSolver(CPX_PARAM_SCRIND=0))
    model=buildlp([-1,0],[2 1],'<',1.5, solver) # create dummy model with correct solver
    MPB.loadproblem!(model,filename) # load what we actually want
    return model
end

c = MPB.getobj(model)

files = readdir(pwd()*"/LP/")
model = loadlp(files[1])
