using DelimitedFiles,MathProgBase,CPLEX,Random,LinearAlgebra
const MPB = MathProgBase; const MOF = MathOptFormat

function loadlp(filename,solver=CplexSolver(CPX_PARAM_SCRIND=0))
    model=buildlp([-1,0],[2 1],'<',1.5, solver) # create dummy model with correct solver
    MPB.loadproblem!(model,filename) # load what we actually want
    return model
end
lpath = "/home/ak121396/Desktop/instances//MIPLIB(official)/LP/"
files = readdir(lpath)
mpspath = "/home/ak121396/Desktop/instances/MIPLIB(official)/mps/"
mpss = readdir(mpspath)

for i=1:35
    # print("file_$i : ",files[i],"\n")
    model = loadlp(lpath*files[i])
    C = MPB.getobj(model)
    # print("num of coeff: ",length(unique(C)), "\n")
    if all(x->x==C[1],C)==true
        print(i,"th ",files[i], "\n")
    end
end


i=10
print("$i file: ",files[i],"\n")

lpmodel = loadlp(lpath*files[i])
C = MPB.getobj(lpmodel)
B = MPB.getconstrmatrix(lpmodel)
m,n=size(B)
lb = MPB.getconstrLB(lpmodel)
ub = MPB.getconstrUB(lpmodel)
RHS = Dict()
for i=1:m
    if ub[i]==Inf
        RHS[i] = lb[i]
    else
        RHS[i] = ub[i]
    end
end
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
obj = 3
objnz = n*obj-(count(iszero,coef)+count(iszero,coef2)+count(iszero,coef3))
NZ = length(B)-count(iszero,B)
# # mpss[i]
# # mpsf = readdlm(mpspath*mpss[i], String, header=true)
# k = findall(i->(i=="L") ||(i=="G")||(i=="E"),mpsf[1][:])[1]
# bound = mpsf[1][k:m+k-1]
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
# RHS = Dict()  #For AP/KP/FLP
# for k=1:m
#     RHS[k]= 1
# end
# for i=1:m
#     if bound[i] == "G"
#         push!(signs,"l")
#     elseif bound[i] == "L"
#         push!(signs,"u")
#     elseif bound[i] == "E"
#         push!(signs, "s")
#     end
# end

ins = open("/home/ak121396/Desktop/instances/MIPLIB(official)/vlp/test2"*files[i][1:end-3]*".vlp","a")
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


1
# mpsf2 = readdlm(path*readdir(path)[i],'\t',String,header=true)
# function rhsindex(k)
#     if k==1
#         for j=1:length(mpsf2[1])
#             if ((occursin("rhs",mpsf2[k][j])==true)||(occursin("RHS",mpsf2[k][j])==true))
#                 # print(j,"\n")
#                 return j
#                 break
#             end
#         end
#     else
#         for j=length(mpsf2[1]):-1:1
#             if ((occursin("rhs",mpsf2[1][j])==true )||(occursin("RHS",mpsf2[1][j])==true))
#                 # print(j,"\n")
#                 return j
#                 break
#             end
#         end
#     end
# end
# function rhsfin()
#     # if k==1
#     for j=1:length(mpsf2[1])
#         if occursin("BOUND",mpsf2[1][j])==true
#             # print(j,"\n")
#             return j
#             break
#         end
#     end
# end
# idx=rhsindex(1)
# finidx=rhsfin()-1
# for k = idx+1:finidx
#     tmp = mpsf2[1][k]
#     tmpp = split(tmp," ")
#     ttmp = filter!(e->e∉["","rhs","RHS","RHS1","B"],tmpp)
#     @show ttmp
#
#     # for y=1:Int(length(ttmp)/2)
#     y=1
#     sptmp = split(ttmp[2*y-1],"EQ")
#     # sptmp = split(ttmp[2*y-1],r"r_|m|c|R|C|EQ|")
#     # if sptmp[2] == "a"
#     key = parse(Int128,sptmp[2])+1
#     val = parse(Float64,ttmp[2*y])
#     RHS[key] = val
#         # if length(sptmp)>2
#         #     sptemp = split(sptmp[1],"(")
#         #     sptemp2 = split(sptemp[2],")")
#         #     key = parse(Int128,sptemp2[3])+1
#         #     val = parse(Float64,ttmp[2*y])
#         #     RHS[key] = val
#         # else
#         #     key = parse(Int128,sptmp[3])+1
#         #     val = parse(Float64,ttmp[2*y])
#         #     RHS[key] = val
#         # end
#     # end
# end
# setdiff(keys(RHS))
# count(iszero,values(RHS))
# Generation method (a)
# for j=1:length(solval)
#     if iszero(solval[j])==false
#     for i=1:length(coef)
#         lb = min(-coef[i],coef[i])
#         ub = max(-coef[i], coef[i])
#         push!( coef2, rand(lb:ub) )
#         push!( coef3, rand(lb:ub) )
#     end
#     for i=1:length(coef)
#         if iszero(solval[i])==false
#             lb = min(-coef[i],coef[i])
#             ub = max(-coef[i], coef[i])
#             if ( ((coef[i]>0) && (solval[i]!= lb)) || ((coef[i]<0) && (solval[i]!= ub)) )
#                 coef2[i] = -coef[i]
#             else
#                 coef2[i] = coef[i]
#             end
#         end
#     end
# end
