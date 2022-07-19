using DataFrames,DelimitedFiles,JuMP,LinearAlgebra,CPLEX,MathProgBase,MathOptInterface,CPUTime
const MPB = MathProgBase
mutable struct Data1
    file::String; N::Dict{}; d::Array{}; c::Array{}; a::Array{}; e::Array{}; gij::SparseVector{}; gjk::SparseVector{}; gkl::SparseVector{};
    vij::Array{}; vjk::Array{}; vkl::Array{}; Vij::SparseVector{}; Vjk::SparseVector{}; Vkl::SparseVector{}; Mij::Array{}; Mjk::Array{}; Mkl::Array{};
    b::Array{}; q::Array{}; rij::Array{}; rjk::Array{}; rkl::Array{}; upl::Int; udc::Int
    function Data1(file)
        dt1 = readdlm(file);
        notafile = readdlm("/home/ak121396/Desktop/instances/SCND/Notations.txt", '=');
        # notafile = readdlm("E:/scnd1/Notations.txt", '=');
        # notafile = readdlm("/home/k2g00/k2g3475/scnd1/Notations.txt", '=');
        nota = notafile[1:end,1];  N= Dict();

        for i=1:length(nota)-1
            id1 = findall(x->x==nota[i], dt1)[1][1];
            id2 = findall(x->x==nota[i+1], dt1)[1][1];
            if id2-id1<3
                tmp = filter(x->x!="",  dt1[id1+(id2-id1-1),:])
                if length(tmp)<2
                    N[nota[i]] = tmp[1];
                else
                    N[nota[i]] = tmp;
                end
            else
                W = []
                for x=id1+1:id1+(id2-id1-1)
                    tmp = filter(x->x!="", dt1[x,:]);
                    append!(W,tmp)
                    # push!(W,tmp);
                end
                # tmp = [filter(x->x!="", dt1[x,:]) for x in id1+1:id1+(id2-id1-1)]
                N[nota[i]] = W;
            end
        end
        d = reshape(N["demand"],5,N["customer"])'; c = append!(N["fcp"],N["fcd"]);
        a = reshape(N["vcs"],5,N["supplier"])';    e = append!(N["vcp"],N["vcd"]);
        gij = sparse(N["fixedcostModesp"]); gjk = sparse(N["fixedcostModepd"]); gkl = sparse(N["fixedcostModedc"]);
        vij = N["tcp"]; vjk = N["tcd"]; vkl = N["tcc"];
        Vij = sparse(N["LcapacityModesp"]); Vjk = sparse(N["LcapacityModepd"]); Vkl =  sparse(N["LcapacityModedc"]);
        Mij = transpose(reshape(N["ModeIJ"], (N["plant"],N["supplier"])));
        Mjk = transpose(reshape(N["ModeJK"], (N["distribution"],N["plant"])));
        Mkl = transpose(reshape(N["ModeKL"], (N["customer"],N["distribution"])));
        b = reshape(N["ves"],5,N["supplier"])';  q = append!(N["vep"],N["ved"]);
        rij = N["cep"]; rjk = N["ced"]; rkl = N["cec"];
        upl = N["upperpants"]; udc = N["upperdistribution"];

        new(file,N,d,c,a,e,gij,gjk,gkl,vij,vjk,vkl,Vij,Vjk,Vkl,Mij,Mjk,Mkl,b,q,rij,rjk,rkl,upl,udc);
    end
end

file = "/home/ak121396/Desktop/instances/SCND/test01S2"
dt = Data1(file);

scnd0 = Model(CPLEX.Optimizer); set_silent(scnd0)
# @variable(scnd0, 0<=y[1:dt.N["plant"]+dt.N["distribution"],1:2]<=1  );
# @variable(scnd0, 0<=uij[i=1:dt.N["supplier"],j=1:dt.N["plant"], m=1:dt.Mij[i,j]] <=1 )
# @variable(scnd0, 0<=ujk[j=1:dt.N["plant"],k=1:dt.N["distribution"],m=1:dt.Mjk[j,k]] <=1);
# @variable(scnd0, 0<=ukl[k=1:dt.N["distribution"],l=1:dt.N["customer"],m=1:dt.Mkl[k,l]] <=1);
# MOI.set(scnd0, MOI.RawParameter("CPX_PARAM_SCRIND"), false )
# MOI.set(scnd0, MOI.RawParameter("CPX_PARAM_THREADS"),1  )
##############################  IP   #####################################
@variable(scnd0, y[1:dt.N["plant"]+dt.N["distribution"],1:2], Bin  );
@variable(scnd0, uij[i=1:dt.N["supplier"],j=1:dt.N["plant"], m=1:dt.Mij[i,j]] , Bin);
@variable(scnd0, ujk[j=1:dt.N["plant"],k=1:dt.N["distribution"],m=1:dt.Mjk[j,k]] , Bin);
@variable(scnd0, ukl[k=1:dt.N["distribution"],l=1:dt.N["customer"],m=1:dt.Mkl[k,l]] , Bin);

@variable(scnd0, 0<= xij[i=1:dt.N["supplier"],j=1:dt.N["plant"],m=1:dt.Mij[i,j],1:5] );
@variable(scnd0, 0<= xjk[j=1:dt.N["plant"],k=1:dt.N["distribution"],m=1:dt.Mjk[j,k],1:5] );
@variable(scnd0, 0<= xkl[k=1:dt.N["distribution"],l=1:dt.N["customer"],m=1:dt.Mkl[k,l],1:5] );
@variable(scnd0, 0<= h[1:dt.N["plant"]+dt.N["distribution"],1:5,1:2] );

#a_ip*x_ijmp
exa = AffExpr(0);
for i=1:dt.N["supplier"]
    for j=1:dt.N["plant"]
        for m=1:dt.Mij[i,j]
            for p=1:5
                # if is_valid(scnd0,xij[i,j,m,p])==true
                add_to_expression!(exa,sum(dt.N["vcs"][i][p]*xij[i,j,m,p]) );
                # end
            end
        end
    end
end
#g_ijm*u_ijm expression
exg = AffExpr(0);
for i=1:dt.N["supplier"]
    idx = 1;
    for j=1:dt.N["plant"]
        # if dt.Mij[i,j]<dt.m
        #     for d=1:dt.m-dt.Mij[i,j]
        #         delete(scnd0,uij[i,j,dt.m-d+1])
        #     end
        # end
        for m=1:dt.Mij[i,j]
            add_to_expression!(exg, dt.gij[i][idx]*uij[i,j,m]);
            idx+=1
        end
    end
end
for j=1:dt.N["plant"]
    idx = 1;
    for k=1:dt.N["distribution"]
        # if dt.Mjk[j,k]<dt.m
        #     for d=1:dt.m-dt.Mjk[j,k]
        #         delete(scnd0,ujk[j,k,dt.m-d+1])
        #     end
        # end
        for m=1:dt.Mjk[j,k]
            add_to_expression!(exg, dt.gjk[j][idx]*ujk[j,k,m]);
            idx+=1
        end
    end
end
for k=1:dt.N["distribution"]
    idx = 1;
    for l=1:dt.N["customer"]
        # if dt.Mkl[k,l]<dt.m
        #     for d=1:dt.m-dt.Mkl[k,l]
        #         delete(scnd0,ukl[k,l,dt.m-d+1])
        #     end
        # end
        for m=1:dt.Mkl[k,l]
            add_to_expression!(exg,dt.gkl[k][idx]*ukl[k,l,m]);
            idx+=1
        end
    end
end

#v_ijmp*x_ijmp expression
exv = AffExpr(0);
for i=1:dt.N["supplier"]
    idx = 1;
    for j=1:dt.N["plant"]
        # if dt.Mij[i,j]<dt.m
        #     for d=1:dt.m-dt.Mij[i,j]
        #         delete(scnd0,xij[i,j,dt.m-d+1,:])
        #     end
        # end
        for m=1:dt.Mij[i,j]
            add_to_expression!(exv,sum(dot.(dt.N["tcp"][i][idx:idx+4],(xij[i,j,m,p] for p=1:5))))#*sqrt((dt.N["pointsupplier"][1][i]-dt.N["pointplant"][1][j])^2+(dt.N["pointsupplier"][2][i]-dt.N["pointplant"][2][j])^2)) );
            idx+=5
        end
    end
end
for j=1:dt.N["plant"]
    idx = 1;
    for k=1:dt.N["distribution"]
        # if dt.Mjk[j,k]<dt.m
        #     for d=1:dt.m-dt.Mjk[j,k]
        #         delete(scnd0,xjk[j,k,dt.m-d+1,:])
        #     end
        # end
        for m=1:dt.Mjk[j,k]
            add_to_expression!(exv,sum(dot.(dt.N["tcd"][j][idx:idx+4],xjk[j,k,m,p] for p=1:5)))#*sqrt((dt.N["pointplant"][1][j]-dt.N["pointdistribution"][1][k])^2+(dt.N["pointplant"][2][j]-dt.N["pointdistribution"][2][k])^2)) );
            idx+=5
        end
    end
end
for k=1:dt.N["distribution"]
    idx = 1;
    for l=1:dt.N["customer"]
        # if dt.Mkl[k,l]<dt.m
        #     for d=1:dt.m-dt.Mkl[k,l]
        #         delete(scnd0,xkl[k,l,dt.m-d+1,1:5])
        #     end
        # end
        for m=1:dt.Mkl[k,l]
            add_to_expression!(exv,sum(dot.(dt.N["tcc"][k][idx:idx+4],xkl[k,l,m,p] for p=1:5)))#*sqrt((dt.N["pointdistribution"][1][k]-dt.N["pointCustmoer"][1][l])^2+(dt.N["pointdistribution"][2][k]-dt.N["pointCustmoer"][2][l])^2)) );
            idx+=5
        end
    end
end
#b_ip*x_ijmp
exb = AffExpr(0);
for i=1:dt.N["supplier"]
    for j=1:dt.N["plant"]
        for m=1:dt.Mij[i,j]
            for p=1:5
                # if is_valid(scnd0,xij[i,j,m,p])==true
                add_to_expression!(exb,sum(dt.b[i,p]*xij[i,j,m,p]) );
                # end
            end
        end
    end
end

exr = AffExpr(0);
for i=1:dt.N["supplier"]
    idx = 1;
    for j=1:dt.N["plant"]
        for m=1:dt.Mij[i,j]
            add_to_expression!(exr,sum(dot.(dt.N["cep"][i][idx:idx+4],xij[i,j,m,p] for p=1:5)))#*sqrt((dt.N["pointsupplier"][1][i]-dt.N["pointplant"][1][j])^2+(dt.N["pointsupplier"][2][i]-dt.N["pointplant"][2][j])^2)) );
            idx+=5
        end
    end
end
for j=1:dt.N["plant"]
    idx = 1;
    for k=1:dt.N["distribution"]
        for m=1:dt.Mjk[j,k]
            add_to_expression!(exr,sum(dot.(dt.N["ced"][j][idx:idx+4],xjk[j,k,m,p] for p=1:5)))#*sqrt((dt.N["pointplant"][1][j]-dt.N["pointdistribution"][1][k])^2+(dt.N["pointplant"][2][j]-dt.N["pointdistribution"][2][k])^2)) );
            idx+=5
        end
    end
end
for k=1:dt.N["distribution"]
    idx = 1;
    for l=1:dt.N["customer"]
        for m=1:dt.Mkl[k,l]
            add_to_expression!(exr,sum(dot.(dt.N["cec"][k][idx:idx+4],xkl[k,l,m,p] for p=1:5)))#*sqrt((dt.N["pointdistribution"][1][k]-dt.N["pointCustmoer"][1][l])^2+(dt.N["pointdistribution"][2][k]-dt.N["pointCustmoer"][2][l])^2)) );
            idx+=5
        end
    end
end
#1st obj
# @constraint(scnd0, obj1, sum(dt.c[j][t]*y[j,t] for j=1:1:dt.N["plant"]+dt.N["distribution"] for t=1:2) + exa + sum(dt.e[j][(p-1)*2+t]*h[j,p,t] for j=1:dt.N["plant"] + dt.N["distribution"] for p=1:5 for t=1:2) + exg + exv <=0);
@objective(scnd0, Min, sum(dt.c[j][t]*y[j,t] for j=1:1:dt.N["plant"]+dt.N["distribution"] for t=1:2) + exg +exa + exv + sum(dt.e[j][(p-1)*2+t]*h[j,p,t] for j=1:(dt.N["plant"]+dt.N["distribution"]) for p=1:5 for t=1:2));
#2nd obj
# @constraint(scnd0, obj2, exb+sum(dt.q[j][(p-1)*2+t]*h[j,p,t] for j=1:dt.N["plant"]+dt.N["distribution"] for p=1:5 for t=1:2) +exr <=0);
# @objective(scnd0,Min,exb + exr + sum(dt.q[j][(p-1)*2+t]*h[j,p,t] for j=1:dt.N["plant"]+dt.N["distribution"] for p=1:5 for t=1:2) );

########### constraint 3 ###############
@constraints(scnd0, begin
    [j=1:dt.N["plant"],p=1:5], sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j]) == sum(xjk[j,k,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k])
    [k=1:dt.N["distribution"],p=1:5], sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k]) == sum(xkl[k,l,m,p] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])
end);
########### constraint 4-6 #############
@constraints(scnd0, begin
    [j=1:dt.N["plant"],p=1:5], sum(h[j,p,:]) == sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j])
    [k=1:dt.N["distribution"],p=1:5], sum(h[k+dt.N["plant"],p,:] ) == sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k])
    [l=1:dt.N["customer"],p=1:5], sum(xkl[k,l,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mkl[k,l]) >= dt.d[l][p]
end )
########### constraint 7-9 #############
@constraint(scnd0,[i=1:dt.N["supplier"]], sum(xij[i,j,m,p] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) <= dt.N["cas"][i]);
@constraint(scnd0,[j=1:dt.N["plant"]+dt.N["distribution"], t=1:2], sum(h[j,1:5,t]) <= [dt.N["cap"];dt.N["cad"]][j]*y[j,t]);
@constraint(scnd0,[j=1:dt.N["plant"]+dt.N["distribution"]], sum(y[j,:]) <= 1);
########### constraint 10 #############
@constraint(scnd0,[i=1:dt.N["supplier"],j=1:dt.N["plant"]], sum(uij[i,j,m] for m=1:dt.Mij[i,j]) <= 1);
@constraint(scnd0,[j=1:dt.N["plant"],k=1:dt.N["distribution"]], sum(ujk[j,k,m] for m=1:dt.Mjk[j,k]) <= 1);
@constraint(scnd0,[k=1:dt.N["distribution"],l=1:dt.N["customer"]], sum(ukl[k,l,m] for m=1:dt.Mkl[k,l]) <= 1);
########### constraint 12 #############
@constraint(scnd0,[i=1:dt.N["supplier"], j=1:dt.N["plant"], m=1:dt.Mij[i,j]], sum(xij[i,j,m,p] for p=1:5) >= dt.Vij[i][j][m]*uij[i,j,m] );
@constraint(scnd0,[j=1:dt.N["plant"], k=1:dt.N["distribution"], m=1:dt.Mjk[j,k]], sum(xjk[j,k,m,p] for p=1:5) >= dt.Vjk[j][k][m]*ujk[j,k,m]);
@constraint(scnd0,[k=1:dt.N["distribution"], l=1:dt.N["customer"], m=1:dt.Mkl[k,l]], sum(xkl[k,l,m,p] for p=1:5) >= dt.Vkl[k][l][m]*ukl[k,l,m]);
########### constraint 13-14 #############
@constraint(scnd0, sum(y[j,t] for j=1:dt.N["plant"] for t=1:dt.m) <= dt.upl);
@constraint(scnd0, sum(y[j,t] for j=dt.N["plant"]+1:dt.N["distribution"]+dt.N["plant"] for t=1:dt.m) <= dt.udc);

optimize!(scnd0)
println(objective_value(scnd0))
1
# with supplivers availability optsol obj1: 1.351098484864179e8
# with BigM obj1: 1.3414953963385555e8
# without BigM obj1: 1.3414953963385555e8
# sum(value.(y))
# sum(value.(uij))
# value.(h)
# dot(dtt.C[2,:][1+54+145+269+2549+725+1345+12745:54+145+269+2549+725+1345+12745+270],value.(h)[:])
# dot(dtt.C[1,:][1+54+145:54+145+269],value.(ujk))
# (dtt.C[2,:])
########### products can be delivered only by chosen transportation mode #############
# BigM = sum(sum(dt.N["demand"]))
# for i=1:dt.N["supplier"]
#     for j=1:dt.N["plant"]
#         for m=1:dt.Mij[i,j]
#             @constraint(scnd0, sum(xij[i,j,m,:] ) <= BigM*uij[i,j,m]);
#         end
#     end
# end
# for j=1:dt.N["plant"]
#     for k=1:dt.N["distribution"]
#         for m=1:dt.Mjk[j,k]
#             @constraint(scnd0, sum(xjk[j,k,m,:] ) <= BigM*ujk[j,k,m]);
#         end
#     end
# end
# for k=1:dt.N["distribution"]
#     for l=1:dt.N["customer"]
#         for m=1:dt.Mkl[k,l]
#             @constraint(scnd0, sum(xkl[k,l,m,:] ) <= BigM*ukl[k,l,m]);
#         end
#     end
# end
########### Suppliers Availibility  ################
# for i=1:dt.N["supplier"]
#     for p=1:5
#         if dt.N["SuppliersAvailibility"][i][p]==0
#             @constraint(scnd0, sum(xij[i,:,:,p])==0)
#         end
#     end
# end
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
