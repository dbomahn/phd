using DataFrames,DelimitedFiles,JuMP,LinearAlgebra,MathProgBase#,CPLEX
const MPB = MathProgBase

function loadlp(filename,solver=CplexSolver(CPX_PARAM_SCRIND=0))
    model=buildlp([-1,0],[2 1],'<',1.5, solver) # create dummy model with correct solver
    MPB.loadproblem!(model,filename) # load what we actually want
    return model
end

mutable struct Data
    filepath::String; N::Dict{}; d::Array{}; m::Int; c::Array{}; a::Array{}; e::Array{}; cap::Array{};
    b::Array{}; q::Array{}; Mij::Array{}; Mjk::Array{}; Mkl::Array{};
    # gij::Array{}; gjk::Array{}; gkl::Array{}; vij::Array{}; vjk::Array{}; vkl::Array{};
    Vij::Array{}; Vjk::Array{}; Vkl::Array{};
    # rij::Array{}; rjk::Array{}; rkl::Array{};
    function Data(filepath)
        dt = readdlm(filepath);
        notafile = readdlm("E:/scnd/Notations.txt", '=');
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
        d = N["demand"]; m = N["transportation"]
        c = append!(N["fcp"],N["fcd"]);
        a = N["vcs"]; e = append!(N["vcp"],N["vcd"]);
        cap = append!(N["cas"],N["cap"],N["cad"]);
        Mij = transpose(reshape(N["ModeIJ"], (N["plant"],N["supplier"])));
        Mjk = transpose(reshape(N["ModeJK"], (N["distribution"],N["plant"])));
        Mkl = transpose(reshape(N["ModeKL"], (N["customer"],N["distribution"])));

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

        new(filepath,N,d,m,c,a,e,cap,b,q,Mij,Mjk,Mkl,Vij,Vjk,Vkl);
    end
end
# gij,gjk,gkl,vij,vjk,vkl,,rij,rjk,rkl
# fpath = "/home/ak121396/Desktop/instances/SCND/"
fpath = "E:/scnd/"
dt = Data(fpath*"Test1S1")
# arc = 1:N["supplier"]*N["plant"]+N["plant"]*N["distribution"]+N["distribution"]*N["customer"]

##########################  Mathematical model  #########################
scnd = Model()#CPLEX.Optimizer
# @variable(scnd, x[1:N["supplier"][1], 1:N["plant"][1],1:N["mode"][1]])
@variable(scnd, y[1:dt.N["plant"]+dt.N["distribution"],1:2] )
@variable(scnd, xij[1:dt.N["supplier"],1:dt.N["plant"],1:dt.m,1:5] )
@variable(scnd, xjk[1:dt.N["plant"],1:dt.N["distribution"],1:dt.m,1:5] )
@variable(scnd, xkl[1:dt.N["distribution"],1:dt.N["customer"],1:dt.m,1:5] )

@variable(scnd, h[1:dt.N["plant"]+dt.N["distribution"],1:5,1:2] )

@variable(scnd, uij[1:dt.N["supplier"],1:dt.N["plant"],1:dt.m] )
@variable(scnd, ujk[1:dt.N["plant"],1:dt.N["distribution"],1:dt.m] )
@variable(scnd, ukl[1:dt.N["distribution"],1:dt.N["customer"],1:dt.m] )
@variable(scnd, )

#g_ijm*u_ijm expression
exg = AffExpr(0)
for i=1:dt.N["supplier"]
    idx = 1;
    for j=1:dt.N["plant"]
        if dt.Mij[i,j]<dt.m
            for d=1:dt.m-dt.Mij[i,j]
                delete(scnd,uij[i,j,dt.m-d+1])
            end
        end
        for m=1:dt.Mij[i,j]
            add_to_expression!(exg, dt.N["fixedcostModesp"][i][idx]*uij[i,j,m])
            idx+=1
        end
    end
end
for j=1:dt.N["plant"]
    idx = 1;
    for k=1:dt.N["distribution"]
        if dt.Mjk[j,k]<dt.m
            for d=1:dt.m-dt.Mjk[j,k]
                delete(scnd,ujk[j,k,dt.m-d+1])
            end
        end
        for m=1:dt.Mjk[j,k]
            add_to_expression!(exg, dt.N["fixedcostModepd"][j][idx]*ujk[j,k,m])
            idx+=1
        end
    end
end
for k=1:dt.N["distribution"]
    idx = 1;
    for l=1:dt.N["customer"]
        if dt.Mkl[k,l]<dt.m
            for d=1:dt.m-dt.Mkl[k,l]
                delete(scnd,ukl[k,l,dt.m-d+1])
            end
        end
        for m=1:dt.Mkl[k,l]
            add_to_expression!(exg,dt.N["fixedcostModedc"][k][idx]*ukl[k,l,m])
            idx+=1
        end
    end
end

#v_ijmp*x_ijmp expression
exv = AffExpr(0)
for i=1:dt.N["supplier"]
    idx = 1;
    for j=1:dt.N["plant"]
        if dt.Mij[i,j]<dt.m
            for d=1:dt.m-dt.Mij[i,j]
                delete(scnd,xij[i,j,dt.m-d+1,1:5])
            end
        end
        for m=1:dt.Mij[i,j]
            add_to_expression!(exv,sum(dot.(dt.N["tcp"][i][idx:idx+4],xij[i,j,m,1:5])) )
            idx+=5
        end
    end
end
for j=1:dt.N["plant"]
    idx = 1;
    for k=1:dt.N["distribution"]
        if dt.Mjk[j,k]<dt.m
            for d=1:dt.m-dt.Mjk[j,k]
                delete(scnd,xjk[j,k,dt.m-d+1,1:5])
            end
        end
        for m=1:dt.Mjk[j,k]
            add_to_expression!(exv,sum(dot.(dt.N["tcd"][j][idx:idx+4],xjk[j,k,m,1:5])) )
            idx+=5
        end
    end
end
for k=1:dt.N["distribution"]
    idx = 1;
    for l=1:dt.N["customer"]
        if dt.Mkl[k,l]<dt.m
            for d=1:dt.m-dt.Mkl[k,l]
                delete(scnd,xkl[k,l,dt.m-d+1,1:5])
            end
        end
        for m=1:dt.Mkl[k,l]
            add_to_expression!(exv,sum(dot.(dt.N["tcc"][k][idx:idx+4],xkl[k,l,m,1:5])) )
            idx+=5
        end
    end
end
#a_ip*x_ijmp
exa = AffExpr(0)
for i=1:dt.N["supplier"]
    for j=1:dt.N["plant"]
        for m=1:dt.m
            for p=1:5
                if is_valid(scnd,xij[i,j,m,p])==true
                    add_to_expression!(exa,sum(dt.a[i][p]*xij[i,j,m,p]) )
                end
            end
        end
    end
end

#b_ip*x_ijmp
exb = AffExpr(0)
for i=1:dt.N["supplier"]
    for j=1:dt.N["plant"]
        for m=1:dt.m
            for p=1:5
                if is_valid(scnd,xij[i,j,m,p])==true
                    add_to_expression!(exb,sum(dt.b[i,p]*xij[i,j,m,p]) )
                end
            end
        end
    end
end

exr = AffExpr(0)
for i=1:dt.N["supplier"]
    idx = 1;
    for j=1:dt.N["plant"]
        for m=1:dt.Mij[i,j]
            add_to_expression!(exr,sum(dot.(dt.N["cep"][i][idx:idx+4],xij[i,j,m,1:5])))
            idx+=5
        end
    end
end
for j=1:dt.N["plant"]
    idx = 1;
    for k=1:dt.N["distribution"]
        for m=1:dt.Mjk[j,k]
            add_to_expression!(exr,sum(dot.(dt.N["ced"][j][idx:idx+4],xjk[j,k,m,1:5])))
            idx+=5
        end
    end
end
for k=1:dt.N["distribution"]
    idx = 1;
    for l=1:dt.N["customer"]
        for m=1:dt.Mkl[k,l]
            add_to_expression!(exr,sum(dot.(dt.N["cec"][k][idx:idx+4],xkl[k,l,m,1:5])))
            idx+=5
        end
    end
end

#1st obj
@constraint(scnd, obj1,
    sum(dt.c[j][t]*y[j,t] for j=1:1:dt.N["plant"]+dt.N["distribution"] for t=1:2)+
    sum(dt.e[j][(p-1)*2+t]*h[j,p,t] for j=1:dt.N["plant"]+dt.N["distribution"] for p=1:5 for t=1:2)+
    exg + exv + exa <= 0)
#2nd obj
@constraint(scnd, obj2, exb+sum(dt.q[j][(p-1)*2+t]*h[j,p,t] for j=1:dt.N["plant"]+dt.N["distribution"] for p=1:5 for t=1:2) +exr <=0)


@constraint(scnd, sum(xij[i,:,:,m]  for i=1:dt.N["supplier"] for m=1:dt.m) == sum(xjk[:,k,:,m]  for k=1:dt.N["distribution"] for m=1:dt.m) )
@constraint(scnd, sum(xjk[j,:,:,m]  for j=1:dt.N["plant"] for m=1:dt.m) == sum(xkl[:,l,:,m]  for l=1:dt.N["customer"] for m=1:dt.m)


@constraint(scnd, dot(data.B[k,:],x) <= dot(xij[1:N[]]))
@constraint(scnd, dot(data.B[k,:],x) == data.RHS[k])






all_variables(scnd)
is_valid(scnd,xjk[1,2,2,1])
1
