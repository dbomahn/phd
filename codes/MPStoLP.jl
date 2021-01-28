using JuMP,MathOptFormat,CPLEX

jump_model = Model()
mof_model = MathOptFormat.MPS.Model()
# cd("../../instances/MIPLIB/")
MOI.read_from_file(mof_model, "C:\\Users\\AK121396\\Downloads\\n3div36.mps")
MOI.copy_to(backend(jump_model), mof_model)
# optimize!(jump_model, with_optimizer=CPLEX.Optimizer)
print(jump_model)


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

writemodel(jump_model,"C:\\Users\\AK121396\\Downloads\\copymodel.lp")
