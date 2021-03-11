using DelimitedFiles,DataFrames,JuMP,CPLEX,MathOptInterface,LinearAlgebra
# MathOptFormat
dir1 = "/home/ak121396/Desktop/instances/FLP/instances/"
dir2 = "/home/ak121396/multiobjective/solvers/kirlikSayin2014/lpfiles/"
# dir = readdir(dir1)

files = readdir(dir1)[1]
1

for idx=2:120
    # subdir = dir1*(dir[s])*'/'
    files = readdir(dir1)
    f = readdlm(dir1*files[idx])
        #facility_i, customer_j
    global i,j = filter!(a->typeof(a)==Int, f[1,:])
    # filter!(a->typeof(a)==Int, f[2,:])
    fixcost = []; capa = [];
    for k=2:i+1
        fix,cp = filter!(a->typeof(a)==Int, f[k,:])
        push!(fixcost,fix)
        push!(capa,cp)
    end

    demand = filter!(a->typeof(a)==Int, f[2+i,:])

    cost = [];
    cost2 = [];
    a,b = size(f)
    for k=i+3:a
        ct = filter!(a->typeof(a)==Int, f[k,:])
        push!(cost2,ct)
        for l=1:length(ct)
            push!(cost,ct[l])
        end
    end
    # JUMP MODEL for FPBH LPfiles
    # flp = Model(CPLEX.Optimizer)
    # @variables(flp, begin
    #     x[1:i,1:j], Bin
    #     y[1:i] ,Bin
    #     z[1:j] ,Bin
    #     k == 1
    # end);
    # @constraint(flp, con1[b=1:j], sum(x[a,b] for a in 1:i) == z[b]);
    # @constraint(flp, con2[a=1:i,b=1:j], x[a,b]<=y[a]);
    # @constraint(flp, conobj2, -sum(demand[b]*z[b] for b=1:j) == 0 );
    # @constraint(flp,conobj3, sum(cost2[a][b]*x[a,b] for a=1:i for b=1:j) == 0 );
    # @objective(flp,Min, sum(fixcost[a]*y[a] for a=1:i));

    # Kirlik solver model
    flp = Model(CPLEX.Optimizer)
    @variables(flp, begin
        x[1:i,1:j], Bin
        y[1:i] ,Bin
        z[1:j] ,Bin
        dummy
    end)
    @constraint(flp, obj1, sum(fixcost[a]*y[a] for a=1:i) <= 1  );
    @constraint(flp, obj2, -sum(demand[b]*z[b] for b=1:j) <= 1 );
    @constraint(flp,obj3, sum(cost2[a][b]*x[a,b] for a=1:i for b=1:j) <= 1 );
    @constraint(flp, con1[b=1:j], sum(x[a,b] for a in 1:i) == z[b]);
    @constraint(flp, con2[a=1:i,b=1:j], x[a,b]<=y[a]);
    @objective(flp, Min, dummy);
    write_to_file(flp,dir2*files[idx][1:end-3]*"lp")
    #normal model
    # flp = Model(CPLEX.Optimizer)
    # @variables(flp, begin
    #     x[1:i,1:j], Bin
    #     y[1:i] ,Bin
    #     z[1:j] ,Bin
    #     k == 1
    # end)
    # @constraint(flp, con1[b=1:j], sum(x[a,b] for a in 1:i) == z[b]);
    # @constraint(flp, con2[a=1:i,b=1:j], x[a,b]<=y[a]);
    # @constraint(flp, conobj1, sum(fixcost[a]*y[a] for a=1:i) >= -1  );
    # @constraint(flp, conobj2, -sum(demand[b]*z[b] for b=1:j) >= -1 );
    # @constraint(flp,conobj3, sum(cost2[a][b]*x[a,b] for a=1:i for b=1:j) >= -1 );
    # @objective(flp,Min,k);
    # optimize!(flp);
    # termination_status(flp);
    # objective_value(flp)
    # xval = JuMP.value.(x)
    # JuMP.value.(y)
    # JuMP.value.(z)
    # print(flp);
    # writemodel(flp,dir1[1:end-10]*"/LP/"*files[idx][1:end-4]*".lp")
end
1
        # For Bensolve vlp format
        # for x=1:length(xval)
        #     print(xval[x],"\n")
        # end
#
#
#         #Parray
#         P = zeros(3,i+j+(i*j))
#
#         for y=1:i
#             P[1,y] = fixcost[y]
#         end
#         for y=1:j
#             P[2,i+y] = -demand[y]
#         end
#         for y=1:i*j
#             P[3,i+j+y] = cost[y]
#         end
#
#         # Barray
#         # 1st constraint
#         B = zeros(j+i*j,(i*j)+i+j)
#         for x=1:j
#             B[x,i+x]= -1
#             for y=1:i
#                 B[x,i+j+j*(y-1)+x] = 1
#            end
#         end
#         # 2nd constraint
#         for x=1:i
#             for y=1:j
#                 B[j*x+y,x] = -1
#                 B[j*x+y,i+j+j*(x-1)+y] = 1
#             end
#         end
#
#         NZ = (i+1)*j+(2*i*j)
#         objNZ = j+i+(i*j)
#         ins = open(dir3*"/"*files[idx]*".vlp","a")
#         wholearray=[];
#         arr=["p vlp min",(i*j)+j,(j*i)+i+j,NZ,3,objNZ]
#         push!(wholearray,arr)
#
#         # Parray
#         for k=1:i
#             push!(wholearray,("o",1,k,fixcost[k]))
#         end
#
#         for k=1:j
#             push!(wholearray,("o",2,i+k,-demand[k]))
#         end
#
#         for k=1:i*j
#             push!(wholearray,("o",3,i+j+k,cost[k]))
#         end
#
#         for x=1:(i*j)+j
#            for y=1:(i*j)+j+i
#                if (B[x,y]!=0)
#                    push!(wholearray,("a",x,y,Int128(B[x,y])))
#                end
#            end
#         end
#
#         for x=1:i*j+j
#             if x<=j
#                 push!(wholearray,("i",x,"s",0))
#             else
#                 push!(wholearray,("i",x,"u",0))
#             end
#         end
#
#         for y=1:(j*i)+j+i
#            push!(wholearray,("j", y,"d",0,1))
#         end
#         push!(wholearray,"e")
#
#         # println(ins,wholearray)
#         writedlm(ins,wholearray)
#         close(ins)
#         print("It's done")
#     end
# end

# old version of write_to_file by using MathOptFormat
# function writemodel(model::Model, filename::String)
#      if endswith(lowercase(filename), ".lp")
#          file = MathOptFormat.LP.Model()
#      elseif endswith(lowercase(filename), ".mps")
#          file = MathOptFormat.MPS.Model()
#      else
#          println("Unknown file extension to write model: ",
#                  split(filename, ".")[end])
#          exit(8)
#      end
#      MOI.copy_to(file, backend(model))
#      MOI.write_to_file(file, filename)
#      # println("wrote model to ", filename)
# end
