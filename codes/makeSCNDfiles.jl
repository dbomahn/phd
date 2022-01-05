using DataFrames,DelimitedFiles,JuMP,LinearAlgebra,CPLEX,MathProgBase,MathOptInterface
const MPB = MathProgBase

function loadlp(filename,solver=CplexSolver(CPX_PARAM_SCRIND=0))
    model=buildlp([-1,0],[2 1],'<',1.5, solver) # create dummy model with correct solver
    MPB.loadproblem!(model,filename) # load what we actually want
    return model
end

mutable struct Data
    filepath::String; N::Dict{}; d::Array{}; m::Int; c::Array{}; e::Array{}; gij::Array{}; gjk::Array{}; gkl::Array{};
    Mij::Array{}; Mjk::Array{}; Mkl::Array{}; Vij::Array{}; Vjk::Array{}; Vkl::Array{}; b::Array{}; q::Array{};
    upl::Int; udc::Int
    # rij::Array{}; rjk::Array{}; rkl::Array{}; vij::Array{}; vjk::Array{}; vkl::Array{};
    function Data(filepath)
        dt = readdlm(filepath);
        # notafile = readdlm("/home/ak121396/Desktop/instances/SCND/Notations.txt", '=');
        # notafile = readdlm("E:/scnd/Notations.txt", '=');
        notafile = readdlm("/home/k2g00/k2g3475/scnd/Notations.txt", '=');
        nota = notafile[1:end,1];  N= Dict();
        for i=1:length(nota)-1
            id1 = findall(x->x==nota[i], dt)[1][1];
            id2 = findall(x->x==nota[i+1], dt)[1][1];
            if id2-id1<3
                tmp = filter(x->x!="",  dt[id1+(id2-id1-1),:])
                if length(tmp)<2
                    N[nota[i]] = tmp[1];
                else
                    N[nota[i]] = tmp;
                end
            else
                W = []
                for x=id1+1:id1+(id2-id1-1)
                    tmp = filter(x->x!="", dt[x,:]);
                    push!(W,tmp);
                end
                # tmp = [filter(x->x!="", dt[x,:]) for x in id1+1:id1+(id2-id1-1)]
                N[nota[i]] = W;
            end
        end
        d = N["demand"];  m = N["transportation"];
        c = append!(N["fcp"],N["fcd"]); e = append!(N["vcp"],N["vcd"]);
        # gij = N["fixedcostModesp"]; gjk = N["fixedcostModepd"];  gkl = N["fixedcostModedc"];
        gij = replace.(N["fixedcostModesp"], 0=>10^(-3));
        gjk = replace.(N["fixedcostModepd"], 0=>10^(-3));
        gkl = replace.(N["fixedcostModedc"], 0=>10^(-3));
        Mij = transpose(reshape(N["ModeIJ"], (N["plant"],N["supplier"])));
        Mjk = transpose(reshape(N["ModeJK"], (N["distribution"],N["plant"])));
        Mkl = transpose(reshape(N["ModeKL"], (N["customer"],N["distribution"])));

        Vij = [];
        for i=1:N["supplier"]
            idx = 1; push!(Vij,[]);
            for j=1:N["plant"]
                th = []
                for m=1:Mij[i,j]
                    push!(th, N["LcapacityModesp"][i][idx]);
                    idx+=1
                end
                push!(Vij[i],th);
            end
        end
        Vjk = [];
        for j=1:N["plant"]
            idx = 1; push!(Vjk,[]);
            for k=1:N["distribution"]
                th = []
                for m=1:Mjk[j,k]
                    push!(th, N["LcapacityModepd"][j][idx]);
                    idx+=1
                end
                push!(Vjk[j],th);
            end
        end
        Vkl = [];
        for k=1:N["distribution"]
            idx = 1; push!(Vkl,[]);
            for l=1:N["customer"]
                th= []
                for m=1:Mkl[k,l]
                    push!(th, N["LcapacityModedc"][k][idx]);
                    idx+=1
                end
                push!(Vkl[k],th);
            end
        end
        b = reshape( N["ves"], (N["supplier"],Int(length(N["ves"])/N["supplier"])) );
        q = append!(N["vep"],N["ved"]);
        upl = N["upperpants"]; udc = N["upperdistribution"]

        new(filepath,N,d,m,c,e,gij,gjk,gkl,Mij,Mjk,Mkl,Vij,Vjk,Vkl,b,q,upl,udc); #cap,Mij,Mjk,Mkl,
    end
end
# @show file = ARGS[1];
# file = "/home/ak121396/Desktop/instances/SCND/test04S3"
# file = "E:/scnd/Test4S3"
file = "/home/k2g00/k2g3475/scnd/instances/test04S3"
dt = Data(file);
##########################  Mathematical model  #########################
scnd = Model(CPLEX.Optimizer); set_silent(scnd)
@variable(scnd, 0<= y[1:dt.N["plant"]+dt.N["distribution"],1:2] <= 1 );
@variable(scnd, 0<= uij[1:dt.N["supplier"],1:dt.N["plant"],1:dt.m] <= 1);
@variable(scnd, 0<= ujk[1:dt.N["plant"],1:dt.N["distribution"],1:dt.m] <= 1);
@variable(scnd, 0<= ukl[1:dt.N["distribution"],1:dt.N["customer"],1:dt.m] <= 1);
##############################  IP   #####################################
# @variable(scnd, y[1:dt.N["plant"]+dt.N["distribution"],1:2], Bin  );
# MOI.set(scnd, MOI.RawParameter("CPX_PARAM_SCRIND"), false )
# MOI.set(scnd, MOI.RawParameter("CPX_PARAM_THREADS"),1  )
# @variable(scnd, uij[1:dt.N["supplier"],1:dt.N["plant"], 1:dt.m] , Bin);
# @variable(scnd, ujk[1:dt.N["plant"],1:dt.N["distribution"],1:dt.m] , Bin);
# @variable(scnd, ukl[1:dt.N["distribution"],1:dt.N["customer"],1:dt.m] , Bin);
############
@variable(scnd, 0<= xij[1:dt.N["supplier"],1:dt.N["plant"],1:dt.m,1:5] );
@variable(scnd, 0<= xjk[1:dt.N["plant"],1:dt.N["distribution"],1:dt.m,1:5] );
@variable(scnd, 0<= xkl[1:dt.N["distribution"],1:dt.N["customer"],1:dt.m,1:5] );
@variable(scnd, 0<= h[1:dt.N["plant"]+dt.N["distribution"],1:5,1:2] );

#1st obj
@constraint(scnd, obj1,
    sum(dt.c[j][t]*y[j,t] for j=1:1:dt.N["plant"]+dt.N["distribution"] for t=1:2) + exa +
    sum(dt.e[j][(p-1)*2+t]*h[j,p,t] for j=1:dt.N["plant"] + dt.N["distribution"] for p=1:5 for t=1:2)+
    exg + exv <=0);
# @objective(scnd, Min, sum(dt.c[j][t]*y[j,t] for j=1:1:dt.N["plant"]+dt.N["distribution"] for t=1:2) + exa +
#     sum(dt.e[j][(p-1)*2+t]*h[j,p,t] for j=1:(dt.N["plant"]+dt.N["distribution"]) for p=1:5 for t=1:2)+
#     exg + exv);
#2nd obj
@constraint(scnd, obj2, exb+sum(dt.q[j][(p-1)*2+t]*h[j,p,t] for j=1:dt.N["plant"]+dt.N["distribution"] for p=1:5 for t=1:2) +exr <=0);
# @objective(scnd,Min,exb+sum(dt.q[j][(p-1)*2+t]*h[j,p,t] for j=1:dt.N["plant"]+dt.N["distribution"] for p=1:5 for t=1:2) +exr );

########### constraint 3 #############
@constraints(scnd, begin
    [j=1:dt.N["plant"],p=1:5], sum(xij[:,j,:,p]) == sum(xjk[j,:,:,p])
    [k=1:dt.N["distribution"],p=1:5], sum(xjk[:,k,:,p]) == sum(xkl[k,:,:,p])
end);
########### constraint 4-6 #############
@constraints(scnd, begin
    [j=1:dt.N["plant"],p=1:5], sum(h[j,p,:]) == sum(xij[:,j,:,p])
    [k=1:dt.N["distribution"],p=1:5], sum(h[k+dt.N["plant"],p,:] ) == sum(xjk[:,k,:,p])
    [l=1:dt.N["customer"],p=1:5], sum(xkl[:,l,:,p] ) >= dt.d[l][p]
end )
########### constraint 7-9 #############
@constraint(scnd,[i=1:dt.N["supplier"]], sum(xij[i,:,:,:] ) <= dt.N["cas"][i]);
@constraint(scnd,[j=1:dt.N["plant"]+dt.N["distribution"], t=1:2], sum(h[j,:,t]) <= [dt.N["cap"];dt.N["cad"]][j]*y[j,t]);
@constraint(scnd,[j=1:dt.N["plant"]+dt.N["distribution"]], sum(y[j,:]) <= 1);
########### constraint 10 #############
@constraint(scnd,[i=1:dt.N["supplier"],j=1:dt.N["plant"]], sum(uij[i,j,:]) <= 1);
@constraint(scnd,[j=1:dt.N["plant"],k=1:dt.N["distribution"]], sum(ujk[j,k,:]) <= 1);
@constraint(scnd,[k=1:dt.N["distribution"],l=1:dt.N["customer"]], sum(ukl[k,l,:]) <= 1);
########### constraint 12 #############
@constraint(scnd,[i=1:dt.N["supplier"], j=1:dt.N["plant"], m=1:dt.Mij[i,j]], sum(xij[i,j,m,:] ) >= dt.Vij[i][j][m]*uij[i,j,m]);
@constraint(scnd,[j=1:dt.N["plant"], k=1:dt.N["distribution"], m=1:dt.Mjk[j,k]], sum(xjk[j,k,m,:]) >= dt.Vjk[j][k][m]*ujk[j,k,m]);
@constraint(scnd,[k=1:dt.N["distribution"], l=1:dt.N["customer"], m=1:dt.Mkl[k,l]], sum(xkl[k,l,m,:]) >= dt.Vkl[k][l][m]*ukl[k,l,m]);

########### constraint 13-14 #############
@constraint(scnd, sum(y[j,t] for j=1:dt.N["plant"] for t=1:dt.m) <= dt.upl);
@constraint(scnd, sum(y[j,t] for j=dt.N["plant"]+1:dt.N["distribution"]+dt.N["plant"] for t=1:dt.m) <= dt.udc);
########### products can be delivered only by chosen transportation mode #############
BigM = 10^(20)
@constraint(scnd,[i=1:dt.N["supplier"], j=1:dt.N["plant"], m=1:dt.m], sum(xij[i,j,m,:] ) <= BigM*uij[i,j,m] );
@constraint(scnd,[j=1:dt.N["plant"], k=1:dt.N["distribution"], m=1:dt.m], sum(xjk[j,k,m,:]) <= BigM*ujk[j,k,m]);
@constraint(scnd,[k=1:dt.N["distribution"], l=1:dt.N["customer"], m=1:dt.m], sum(xkl[k,l,m,:]) <= BigM*ukl[k,l,m]);
########### Suppliers Availibility  ################
for i=1:dt.N["supplier"]
    for p=1:5
        if dt.N["SuppliersAvailibility"][i][p]==0
            @constraint(scnd, sum(xij[i,:,:,p])==0)
        end
    end
end

optimize!(scnd)
objective_value(scnd)
termination_status(scnd)

sum(value.(uij))
(value.(y))[:,2]
#################### Write lp file & Call it as mps format #####################
write_to_file(scnd, file*".lp")

write_to_file(scnd , "/home/k2g00/k2g3475/scnd/lp/"*file[36:end]*".lp")
lpmodel = loadlp("/home/k2g00/k2g3475/scnd/lp/"*file[36:end]*".lp")

Bmtx = MPB.getconstrmatrix(lpmodel);
P = Bmtx[1:2,:]
B = Bmtx[3:end,:]
m,n=size(B)
varlb = MPB.getvarLB(lpmodel)
varub = MPB.getvarUB(lpmodel)
lb = MPB.getconstrLB(lpmodel)[3:end]
ub = MPB.getconstrUB(lpmodel)[3:end]
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
for j=1:n
    if varub[j]==1
        push!(wholearray,("j",j,"d",0,1))
    else
        push!(wholearray,("j", j,'l',0))
    end
end
push!(wholearray,"e")

ins = open("/home/k2g00/k2g3475/scnd/vlp/"*file[36:end]*".vlp","w")
writedlm(ins,wholearray)
close(ins)


##################   Data  ########################
# gij = [];
# for i=1:N["supplier"]
#     idx = 1; push!(gij,[]);
#     for j=1:N["plant"]
#         fc = []
#         for m=1:Mij[i,j]
#             push!(fc, N["fixedcostModesp"][i][idx]);
#             idx+=1
#         end
#         push!(gij[i],fc);
#     end
# end
# gjk = [];
# for j=1:N["plant"]
#     idx = 1; push!(gjk,[]);
#     for k=1:N["distribution"]
#         fc = []
#         for m=1:Mjk[j,k]
#             push!(fc, N["fixedcostModepd"][j][idx]);
#             idx+=1
#         end
#         push!(gjk[j],fc);
#     end
# end
# gkl = [];
# for k=1:N["distribution"]
#     idx = 1; push!(gkl,[]);
#     for l=1:N["customer"]
#         fc = []
#         for m=1:Mkl[k,l]
#             push!(fc, N["fixedcostModedc"][k][idx]);
#             idx+=1
#         end
#         push!(gkl[k],fc);
#     end
# end
# vij = [];
# for i=1:N["supplier"]
#     idx = 1; push!(vij,[]);
#     for j=1:N["plant"]
#         tc = []
#         for m=1:Mij[i,j]
#             push!(tc, N["tcp"][i][idx:idx+4]);
#             idx+=5
#         end
#         push!(vij[i],tc);
#     end
# end
# vjk = [];
# for j=1:N["plant"]
#     idx = 1; push!(vjk,[]);
#     for k=1:N["distribution"]
#         tc = []
#         for m=1:Mjk[j,k]
#             push!(tc, N["tcd"][j][idx:idx+4]);
#             idx+=5
#         end
#         push!(vjk[j],tc);
#     end
# end
# vkl = [];
# for k=1:N["distribution"]
#     idx = 1; push!(vkl,[]);
#     for l=1:N["customer"]
#         tc = []
#         for m=1:Mkl[k,l]
#             push!(tc, N["tcc"][k][idx:idx+4]);
#             idx+=5
#         end
#         push!(vkl[k],tc);
#     end
# end
# rij = [];
# for i=1:N["supplier"]
#     idx = 1; push!(rij,[]);
#     for j=1:N["plant"]
#         em = []
#         for m=1:Mij[i,j]
#             push!(em, N["cep"][i][idx:idx+4]);
#             idx+=5
#         end
#         push!(rij[i],em);
#     end
# end
# rjk = []
# for j=1:N["plant"]
#     idx = 1; push!(rjk,[]);
#     for k=1:N["distribution"]
#         em = []
#         for m=1:Mjk[j,k]
#             push!(em, N["ced"][j][idx:idx+4]);
#             idx+=5
#         end
#         push!(rjk[j],em);
#     end
# end
# rkl = []
# for k=1:N["distribution"]
#     idx = 1; push!(rkl,[]);
#     for l=1:N["customer"]
#         em = []
#         for m=1:Mkl[k,l]
#             push!(em, N["cec"][k][idx:idx+4]);
#             idx+=5
#         end
#         push!(rkl[k],em);
#     end
# end
