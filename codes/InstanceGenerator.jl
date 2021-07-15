using JuMP,LinearAlgebra,CPLEX,DelimitedFiles

#############################  Generating the 3rd obj  ########################
path = "C:\\Users\\AK121396\\Downloads\\tri-objBoland2015b\\"
folders = readdir(path)
# for i=1:length(folders)
#     input = readdir(path*folders[i])
#     for j=1:length(input)
#         data=readdlm(path*folders[i]*"\\"*input[j])
#         m,nc,nb = data[1:3]
#         C3 = rand(-10:10,(nc,1)); F3 = rand(-200:200,(nb,1))
#         extC = vcat(C3,F3)
#
#         open(path*folders[i]*"\\"*input[j],"a") do f
#             print(f,"\n")
#             for i in extC
#                 print(f,i," ")
#             end
#         end
#     end
# end

######################       formating the data        #####################
input = readdir(path*folders[2])
data=readdlm(path*folders[2]*"\\"*input[1])
m,nc,nb = data[1:3]
C = zeros(Int,(3,2*nc));
C[1,:] = vcat(data[4,:][1:nc],data[5,:][1:nb])
C[2,:] = vcat(data[6,:][1:nc],data[7,:][1:nb])
A = data[8:7+nc,1:nb]
Ap = data[end-2,1:nb]
B = data[end-1,1:m-1]
C[3,:] = data[end,1:nc+nb]

#####################    Mathematical Model  ################################
model = Model(CPLEX.Optimizer)
@variable(model, x[1:nc+nb], Bin)
@constraint(model, con1[j=1:nb], sum(A[i,j]*x[i] for i=1:nc) + Ap[j]*x[nc+j] <= B[j])
@constraint(model, con2[j=nb+1:m-1], sum(A[i,j]*x[i] for i=1:nc) <= B[j])

# @objective(model, dot(C1,x))
# @objective(model, dot(C2,x))
