#################### Write lp file & Call it as mps format #####################
using MathProgBase
const MPB = MathProgBase

function loadlp(filename,solver=CplexSolver(CPX_PARAM_SCRIND=0))
    model=buildlp([-1,0],[2 1],'<',1.5, solver) # create dummy model with correct solver
    MPB.loadproblem!(model,filename) # load what we actually want
    return model
end


write_to_file(scnd, file*"md.lp")
lpmodel = loadlp(file*"md.lp")
Bmtx = MPB.getconstrmatrix(lpmodel);
# cut = findall(i-> varub[i]==1 && varub[i+1]!=1, 1:length(varub))[end]
# B = Bmtx[3:end,1:cut];P = Bmtx[1:2,1:cut]; vub = varub[1:cut]
B = Bmtx[3:end,:];P = Bmtx[1:2,:]; vub = MPB.getvarUB(lpmodel)
m,n=size(B)
lb = MPB.getconstrLB(lpmodel)[3:end]
ub = MPB.getconstrUB(lpmodel)[3:end]
RHS = []
for i=1:m
    if ub[i]==Inf
        push!(RHS,lb[i])
    else
        push!(RHS,ub[i])
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
########################## Make vlp file for Bensolve #######################
nz = count(i->(i!=0),B)
objnz = count(i->(i!=0),P)
obj=size(P)[1]
wholearray=[];
arr=["p vlp min",m,n,nz,obj,objnz]
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
       if P[i,j]!=0
           push!(wholearray,("o",i,j,P[i,j]))
       end
   end
end
for i=1:m
   push!(wholearray,("i",i,signs[i],RHS[i]))
end
# for j=1:n
#     if j in bvar
#         push!(wholearray,("j",j,"s",fpx[1][j])) #assign FP int var values
#     else
#         push!(wholearray,("j", j,'l',0))
#     end
# end
for j=1:n
    if vub[j]==1
        push!(wholearray,("j",j,"d",0,1))
    else
        push!(wholearray,("j", j,'l',0))
    end
end
push!(wholearray,"e")

ins = open("/home/k2g00/k2g3475/scnd/vlp/"*file[36:end]*"md.vlp","w")
# ins = open(file*".vlp","w")
writedlm(ins,wholearray)
close(ins)

# function voptmodel()
#     ##########################  Mathematical model  #########################
#     scnd = vModel(CPLEX.Optimizer); set_silent(scnd)
#     # # MOI.set(scnd, MOI.NumberOfThreads(), 1);
#     @variable(scnd, 0<= y[1:dt.N["plant"]+dt.N["distribution"],1:2] <= 1 );
#     @variable(scnd, 0<= uij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:dt.Mij[i,j]] <= 1);
#     @variable(scnd, 0<= ujk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:dt.Mjk[j,k]] <= 1);
#     @variable(scnd, 0<= ukl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:dt.Mkl[k,l]] <= 1);
#     ##############################  IP   #####################################
#     # scnd = vModel(CPLEX.Optimizer); set_silent(scnd);
#     #MOI.set(scnd, MOI.NumberOfThreads(), 1);
#     ############
#     @variable(scnd, 0<= xij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:dt.Mij[i,j],1:5] );
#     @variable(scnd, 0<= xjk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:dt.Mjk[j,k],1:5] );
#     @variable(scnd, 0<= xkl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:dt.Mkl[k,l],1:5] );
#     @variable(scnd, 0<= h[1:dt.N["plant"]+dt.N["distribution"],1:5,1:2] );
#     exg = AffExpr(0);
#     for i=1:dt.N["supplier"]
#         for j=1:dt.N["plant"]
#             add_to_expression!(exg, 10000*uij[i,j,1]);
#         end
#     end
#     for j=1:dt.N["plant"]
#         for k=1:dt.N["distribution"]
#             add_to_expression!(exg, 10000*ujk[j,k,1]);
#         end
#     end
#     for k=1:dt.N["distribution"]
#         for l=1:dt.N["customer"]
#             add_to_expression!(exg,10000*ukl[k,l,1]);
#         end
#     end

#     # 1st obj
#     @addobjective(scnd, Min, sum(dt.c[j][t]*y[j,t] for j=1:1:dt.N["plant"]+dt.N["distribution"] for t=1:2) + exg +
#         sum(dt.N["vcs"][i][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) +
#         sum(dt.vij[i][j][m][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
#         sum(dt.vjk[j][k][m][p]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
#         sum(dt.vkl[k][l][m][p]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5)+
#         sum(dt.N["vcp"][j][2*(p-1)+t]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
#         sum(dt.N["vcd"][k][2*(p-1)+t]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2) )
#     #2nd obj
#     @addobjective(scnd, Min, sum(dt.b[i,p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
#         sum(dt.N["vep"][j][2*(p-1)+t]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
#         sum(dt.N["ved"][k][2*(p-1)+t]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2)+
#         sum(dt.rij[i][j][m]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
#         sum(dt.rjk[j][k][m]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
#         sum(dt.rkl[k][l][m]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5) );
   
#     @constraint(scnd,[j=1:dt.N["plant"],p=1:5], sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j]) == sum(xjk[j,k,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]));
#     @constraint(scnd,[k=1:dt.N["distribution"],p=1:5], sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k]) == sum(xkl[k,l,m,p] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l]));
#     ########### constraint 4-6 #############
#     # @constraints(scnd, begin
#     #     [j=1:dt.N["plant"],p=1:5], sum(h[j,p,:]) == sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j])
#     #     [k=1:dt.N["distribution"],p=1:5], sum(h[k+dt.N["plant"],p,:] ) == sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k])
#     #     [l=1:dt.N["customer"],p=1:5], sum(xkl[k,l,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mkl[k,l]) >= dt.d[l][p]
#     # end );
#     @constraint(scnd,[j=1:dt.N["plant"],p=1:5], sum(h[j,p,:]) == sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j]));
#     @constraint(scnd,[k=1:dt.N["distribution"],p=1:5], sum(h[k+dt.N["plant"],p,:] ) == sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k]));
#     @constraint(scnd, [l=1:dt.N["customer"],p=1:5], sum(xkl[k,l,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mkl[k,l]) >= dt.d[l][p]);
#     ########### constraint 7-9 #############
#     @constraint(scnd,[i=1:dt.N["supplier"]], sum(xij[i,j,m,p] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) <= dt.N["cas"][i]);
#     @constraint(scnd,[j=1:dt.N["plant"]+dt.N["distribution"], t=1:2], sum(h[j,:,t]) <=[dt.N["cap"];dt.N["cad"]][j]*y[j,t]);
#     @constraint(scnd,[j=1:dt.N["plant"]+dt.N["distribution"]], sum(y[j,:]) <= 1);
#     # # ########### constraint 10 #############
#     @constraint(scnd,[i=1:dt.N["supplier"],j=1:dt.N["plant"]], sum(uij[i,j,m] for m=1:dt.Mij[i,j]) <= 1);
#     @constraint(scnd,[j=1:dt.N["plant"],k=1:dt.N["distribution"]], sum(ujk[j,k,m] for m=1:dt.Mjk[j,k]) <= 1);
#     @constraint(scnd,[k=1:dt.N["distribution"],l=1:dt.N["customer"]], sum(ukl[k,l,m] for m=1:dt.Mkl[k,l]) <= 1);
#     ########### constraint 11 #############
#     @constraint(scnd,[i=1:dt.N["supplier"],j=1:dt.N["plant"],m=1:dt.Mij[i,j]], sum(xij[i,j,m,p] for p=1:5) <= dt.bigM*uij[i,j,m]);
#     @constraint(scnd,[j=1:dt.N["plant"],k=1:dt.N["distribution"],m=1:dt.Mjk[j,k]] ,sum(xjk[j,k,m,p] for p=1:5) <= dt.bigM*ujk[j,k,m]);
#     @constraint(scnd,[k=1:dt.N["distribution"],l=1:dt.N["customer"],m=1:dt.Mkl[k,l]], sum(xkl[k,l,m,p] for p=1:5) <= dt.bigM*ukl[k,l,m]);
#     # ########### constraint 12 #############
#     @constraint(scnd,[i=1:dt.N["supplier"], j=1:dt.N["plant"], m=1:dt.Mij[i,j]], sum(xij[i,j,m,p] for p=1:5) >= dt.Vij[m]*uij[i,j,m] );
#     @constraint(scnd,[j=1:dt.N["plant"], k=1:dt.N["distribution"], m=1:dt.Mjk[j,k]], sum(xjk[j,k,m,p] for p=1:5) >= dt.Vjk[m]*ujk[j,k,m]);
#     # @constraint(scnd,[k=1:dt.N["distribution"], l=1:dt.N["customer"], m=1:dt.Mkl[k,l]], sum(xkl[k,l,m,p] for p=1:5) >= dt.Vkl[k][l][m]*ukl[k,l,m]);
#     # ########### constraint 13-14 #############
#     @constraint(scnd,sum(y[j,t] for j=1:dt.N["plant"] for t=1:2) <= dt.upl);
#     @constraint(scnd,sum(y[j,t] for j=dt.N["plant"]+1:dt.N["distribution"]+dt.N["plant"] for t=1:2) <= dt.udc);
#     return scnd
# end
# scndm = voptmodel()
# vSolve(scndm, method=:dicho, verbose=false) 

# function Defexp()
#     #a_ip*x_ijmp
#     exa = AffExpr(0);
#     for i=1:dt.N["supplier"]
#         for j=1:dt.N["plant"]
#             for m=1:dt.Mij[i,j]
#                 for p=1:5
#                     # if is_valid(scnd,xij[i,j,m,p])==true
#                     add_to_expression!(exa,sum(dt.N["vcs"][i][p]*xij[i,j,m,p]));
#                     # end
#                 end
#             end
#         end
#     end
#

#
#     #v_ijmp*x_ijmp expression
#     exv = AffExpr(0);
#     for i=1:dt.N["supplier"]
#         idx = 1;
#         for j=1:dt.N["plant"]
#             # if dt.Mij[i,j]<dt.m
#             #     for d=1:dt.m-dt.Mij[i,j]
#             #         delete(scnd,xij[i,j,dt.m-d+1,:])
#             #     end
#             # end
#             for m=1:dt.Mij[i,j]
#                 add_to_expression!(exv,sum(dot.(dt.N["tcp"][i][idx:idx+4],(xij[i,j,m,p] for p=1:5))))
#                 idx+=5
#             end
#         end
#     end
#     for j=1:dt.N["plant"]
#         idx = 1;
#         for k=1:dt.N["distribution"]
#             # if dt.Mjk[j,k]<dt.m
#             #     for d=1:dt.m-dt.Mjk[j,k]
#             #         delete(scnd,xjk[j,k,dt.m-d+1,:])
#             #     end
#             # end
#             for m=1:dt.Mjk[j,k]
#                 add_to_expression!(exv,sum(dot.(dt.N["tcd"][j][idx:idx+4],xjk[j,k,m,p] for p=1:5)))#*sqrt((dt.N["pointplant"][1][j]-dt.N["pointdistribution"][1][k])^2+(dt.N["pointplant"][2][j]-dt.N["pointdistribution"][2][k])^2)) );
#                 idx+=5
#             end
#         end
#     end
#     for k=1:dt.N["distribution"]
#         idx = 1;
#         for l=1:dt.N["customer"]
#             # if dt.Mkl[k,l]<dt.m
#             #     for d=1:dt.m-dt.Mkl[k,l]
#             #         delete(scnd,xkl[k,l,dt.m-d+1,1:5])
#             #     end
#             # end
#             for m=1:dt.Mkl[k,l]
#                 add_to_expression!(exv,sum(dot.(dt.N["tcc"][k][idx:idx+4],xkl[k,l,m,p] for p=1:5)))#*sqrt((dt.N["pointdistribution"][1][k]-dt.N["pointCustmoer"][1][l])^2+(dt.N["pointdistribution"][2][k]-dt.N["pointCustmoer"][2][l])^2)) );
#                 idx+=5
#             end
#         end
#     end
#     #b_ip*x_ijmp
#     exb = AffExpr(0);
#     for i=1:dt.N["supplier"]
#         for j=1:dt.N["plant"]
#             for m=1:dt.Mij[i,j]
#                 for p=1:5
#                     # if is_valid(scnd,xij[i,j,m,p])==true
#                     add_to_expression!(exb,sum(dt.b[i,p]*xij[i,j,m,p]) );
#                     # end
#                 end
#             end
#         end
#     end
#     exr = AffExpr(0);
#     for i=1:dt.N["supplier"]
#         idx = 1;
#         for j=1:dt.N["plant"]
#             for m=1:dt.Mij[i,j]
#                 add_to_expression!(exr,sum(dot.(dt.N["cep"][i][idx:idx+4],xij[i,j,m,p] for p=1:5)))#*sqrt((dt.N["pointsupplier"][1][i]-dt.N["pointplant"][1][j])^2+(dt.N["pointsupplier"][2][i]-dt.N["pointplant"][2][j])^2)) );
#                 idx+=5
#             end
#         end
#     end
#     for j=1:dt.N["plant"]
#         idx = 1;
#         for k=1:dt.N["distribution"]
#             for m=1:dt.Mjk[j,k]
#                 add_to_expression!(exr,sum(dot.(dt.N["ced"][j][idx:idx+4],xjk[j,k,m,p] for p=1:5)))#*sqrt((dt.N["pointplant"][1][j]-dt.N["pointdistribution"][1][k])^2+(dt.N["pointplant"][2][j]-dt.N["pointdistribution"][2][k])^2)) );
#                 idx+=5
#             end
#         end
#     end
#     for k=1:dt.N["distribution"]
#         idx = 1;
#         for l=1:dt.N["customer"]
#             for m=1:dt.Mkl[k,l]
#                 add_to_expression!(exr,sum(dot.(dt.N["cec"][k][idx:idx+4],xkl[k,l,m,p] for p=1:5)))#*sqrt((dt.N["pointdistribution"][1][k]-dt.N["pointCustmoer"][1][l])^2+(dt.N["pointdistribution"][2][k]-dt.N["pointCustmoer"][2][l])^2)) );
#                 idx+=5
#             end
#         end
#     end
#     return exa,exg,exv,exb,exr
# end
# exa,exg,exv,exb,exr = Defexp();
