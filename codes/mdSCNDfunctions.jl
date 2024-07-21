using CPUTime,DataFrames,DelimitedFiles,JuMP,LinearAlgebra,CPLEX,vOptGeneric,SparseArrays,StatsBase,CSV,JLD2,LazySets,PolygonInbounds
struct Data
    filepath::String; N::Dict{}; d::Array{}; c::Array{};  Mij::Array{}; Mjk::Array{}; Mkl::Array{};
    vij::Array{}; vjk::Array{}; vkl::Array{}; rij::Array{}; rjk::Array{}; rkl::Array{}; #gij::Array{}; gjk::Array{}; gkl::Array{};
    Vij::Array{}; Vjk::Array{}; b::Array{}; upl::Int; udc::Int; bigM::Int # e::Array{};q::Array{};
    function Data(filepath)
        dt = readdlm(filepath);
        # notafile = readdlm("/home/ak121396/Desktop/instances/scnd/Notations.txt", '=');
        # notafile = readdlm("F:/scnd/Notations.txt", '=');
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
        d = N["demand"];  c = append!(N["fcp"],N["fcd"]);
        Mij = transpose(reshape(N["ModeIJ"], (N["plant"],N["supplier"])));
        Mjk = transpose(reshape(N["ModeJK"], (N["distribution"],N["plant"])));
        Mkl = transpose(reshape(N["ModeKL"], (N["customer"],N["distribution"])));
        vij = [];
        for i=1:N["supplier"]
            tmp = [];
            for j=1:N["plant"]
                tmp2 = []
                if j==1
                    for m=1:Mij[i,1]
                        push!(tmp2, N["tcp"][i][5*(m-1)+1:5*(m-1)+5])
                    end
                else
                    for m=1:Mij[i,j]
                        push!(tmp2, N["tcp"][i][5*sum(Mij[i,1:j-1])+5*(m-1)+1:5*sum(Mij[i,1:j-1])+5*(m-1)+5]);
                    end
                end
                push!(tmp,tmp2);
            end
            push!(vij,tmp);
        end
        vjk = [];
        for j=1:N["plant"]
            tmp = [];
            for k=1:N["distribution"]
                tmp2 = []
                if k==1
                    for m=1:Mjk[j,1]
                        push!(tmp2, N["tcd"][j][5*(m-1)+1:5*(m-1)+5])
                    end
                else
                    for m=1:Mjk[j,k]
                        push!(tmp2, N["tcd"][j][5*sum(Mjk[j,1:k-1])+5*(m-1)+1:5*sum(Mjk[j,1:k-1])+5*(m-1)+5]);
                    end
                end
                push!(tmp,tmp2);
            end
            push!(vjk,tmp);
        end
        vkl = [];
        for k=1:N["distribution"]
            tmp = [];
            for l=1:N["customer"]
                tmp2 = []
                if l==1
                    for m=1:Mkl[k,1]
                        push!(tmp2, N["tcc"][k][5*(m-1)+1:5*(m-1)+5])
                    end
                else
                    for m=1:Mkl[k,l]
                        push!(tmp2, N["tcc"][k][5*sum(Mkl[k,1:l-1])+5*(m-1)+1:5*sum(Mkl[k,1:l-1])+5*(m-1)+5]);
                    end
                end
                push!(tmp,tmp2);
            end
            push!(vkl,tmp);
        end

        rij = [];
        for i=1:N["supplier"]
            tmp = [];
            for j=1:N["plant"]
                tmp2 = []
                if j==1
                    for m=1:Mij[i,1]
                        push!(tmp2, unique(N["cep"][i][5*(m-1)+1:5*(m-1)+5])[1])
                    end
                else
                    for m=1:Mij[i,j]
                        push!(tmp2, unique(N["cep"][i][5*sum(Mij[i,1:j-1])+5*(m-1)+1:5*sum(Mij[i,1:j-1])+5*(m-1)+5])[1]);
                    end
                end
                push!(tmp,tmp2);
            end
            push!(rij,tmp);
        end
        rjk = [];
        for j=1:N["plant"]
            tmp = [];
            for k=1:N["distribution"]
                tmp2 = []
                if k==1
                    for m=1:Mjk[j,1]
                        push!(tmp2, unique(N["ced"][j][5*(m-1)+1:5*(m-1)+5])[1])
                    end
                else
                    for m=1:Mjk[j,k]
                        push!(tmp2, unique(N["ced"][j][5*sum(Mjk[j,1:k-1])+5*(m-1)+1:5*sum(Mjk[j,1:k-1])+5*(m-1)+5])[1]);
                    end
                end
                push!(tmp,tmp2);
            end
            push!(rjk,tmp);
        end
        rkl = [];
        for k=1:N["distribution"]
            tmp = [];
            for l=1:N["customer"]
                tmp2 = []
                if l==1
                    for m=1:Mkl[k,1]
                        push!(tmp2, unique(N["cec"][k][5*(m-1)+1:5*(m-1)+5])[1])
                    end
                else
                    for m=1:Mkl[k,l]
                        push!(tmp2, unique(N["cec"][k][5*sum(Mkl[k,1:l-1])+5*(m-1)+1:5*sum(Mkl[k,1:l-1])+5*(m-1)+5])[1]);
                    end
                end
                push!(tmp,tmp2);
            end
            push!(rkl,tmp);
        end
        Vij = [0, maximum(maximum(N["LcapacityModesp"][i] for i=1:N["supplier"]))]
        Vjk = [0, maximum(maximum(N["LcapacityModepd"][i] for i=1:N["plant"]))]
        # Vij = [];
        # for i=1:N["supplier"]
        #     idx = 1; push!(Vij,[]);
        #     for j=1:N["plant"]
        #         th = []
        #         for m=1:Mij[i,j]
        #             push!(th, N["LcapacityModesp"][i][idx]);
        #             idx+=1
        #         end
        #         push!(Vij[i],th);
        #     end
        # end
        # Vjk = [];
        # for j=1:N["plant"]
        #     idx = 1; push!(Vjk,[]);
        #     for k=1:N["distribution"]
        #         th = []
        #         for m=1:Mjk[j,k]
        #             push!(th, N["LcapacityModepd"][j][idx]);
        #             idx+=1
        #         end
        #         push!(Vjk[j],th);
        #     end
        # end
        # Vkl = [];
        # for k=1:N["distribution"]
        #     idx = 1; push!(Vkl,[]);
        #     for l=1:N["customer"]
        #         th= []
        #         for m=1:Mkl[k,l]
        #             push!(th, N["LcapacityModedc"][k][idx]);
        #             idx+=1
        #         end
        #         push!(Vkl[k],th);
        #     end
        # end
        b = reshape(N["ves"],N["supplier"],5);
        # q = append!(N["vep"],N["ved"]);
        upl = N["upperpants"]; udc = N["upperdistribution"]
        bigM = sum(sum(N["demand"]))
        new(filepath,N,d,c,Mij,Mjk,Mkl,vij,vjk,vkl,rij,rjk,rkl,Vij,Vjk,b,upl,udc,bigM); #cap,Mij,Mjk,Mkl,
    end
end
function SCND_LP()
    ##########################  Mathematical model  #########################
    scnd = vModel(CPLEX.Optimizer); set_silent(scnd);
    optimizer_with_attributes(
            CPLEX.Optimizer,
            "CPX_PARAM_EPGAP" => 1e-4);
    @variable(scnd, 0<= y[1:dt.N["plant"]+dt.N["distribution"],1:2] <= 1);
    @variable(scnd, 0<= uij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:2] <= 1)
    @variable(scnd, 0<= ujk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:2] <= 1)
    @variable(scnd, 0<= ukl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:2] <= 1)
    ############
    @variable(scnd, 0<= xij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:2,1:5] );
    @variable(scnd, 0<= xjk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:dt.Mjk[j,k],1:5] );
    @variable(scnd, 0<= xkl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:dt.Mkl[k,l],1:5] );
    @variable(scnd, 0<= h[1:dt.N["plant"]+dt.N["distribution"],1:5,1:2] );
    exg = AffExpr(0);
    for i=1:dt.N["supplier"]
        for j=1:dt.N["plant"]
            add_to_expression!(exg, 10000*uij[i,j,1]);
        end
    end
    for j=1:dt.N["plant"]
        for k=1:dt.N["distribution"]
            add_to_expression!(exg, 10000*ujk[j,k,1]);
        end
    end
    for k=1:dt.N["distribution"]
        for l=1:dt.N["customer"]
            add_to_expression!(exg,10000*ukl[k,l,1]);
        end
    end
    @addobjective(scnd, Min, sum(dt.c[j][t]*y[j,t] for j=1:1:dt.N["plant"]+dt.N["distribution"] for t=1:2) + exg +
        sum(dt.N["vcs"][i][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) +
        sum(dt.vij[i][j][m][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
        sum(dt.vjk[j][k][m][p]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
        sum(dt.vkl[k][l][m][p]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5)+
        sum(dt.N["vcp"][j][2*(p-1)+t]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
        sum(dt.N["vcd"][k][2*(p-1)+t]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2));
    @addobjective(scnd, Min, sum(dt.b[i,p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
        sum(dt.N["vep"][j][2*(p-1)+t]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
        sum(dt.N["ved"][k][2*(p-1)+t]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2)+
        sum(dt.rij[i][j][m]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
        sum(dt.rjk[j][k][m]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
        sum(dt.rkl[k][l][m]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5));
    @constraint(scnd,[j=1:dt.N["plant"],p=1:5], sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j]) == sum(xjk[j,k,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]));
    @constraint(scnd,[k=1:dt.N["distribution"],p=1:5], sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k]) == sum(xkl[k,l,m,p] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l]));
    ########### constraint 4-6 #############
    @constraint(scnd,[j=1:dt.N["plant"],p=1:5], sum(h[j,p,:]) == sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j]));
    @constraint(scnd,[k=1:dt.N["distribution"],p=1:5], sum(h[k+dt.N["plant"],p,:] ) == sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k]));
    @constraint(scnd, [l=1:dt.N["customer"],p=1:5], sum(xkl[k,l,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mkl[k,l]) >= dt.d[l][p]);
    ########### constraint 7-9 #############
    @constraint(scnd,[i=1:dt.N["supplier"]], sum(xij[i,j,m,p] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) <= dt.N["cas"][i]);
    @constraint(scnd,[j=1:dt.N["plant"]+dt.N["distribution"], t=1:2], sum(h[j,:,t]) <=[dt.N["cap"];dt.N["cad"]][j]*y[j,t]);
    @constraint(scnd,[j=1:dt.N["plant"]+dt.N["distribution"]], sum(y[j,:]) <= 1);
    ########### constraint 10 #############
    @constraint(scnd,[i=1:dt.N["supplier"],j=1:dt.N["plant"]], sum(uij[i,j,m] for m=1:dt.Mij[i,j]) <= 1);
    @constraint(scnd,[j=1:dt.N["plant"],k=1:dt.N["distribution"]], sum(ujk[j,k,m] for m=1:dt.Mjk[j,k]) <= 1);
    @constraint(scnd,[k=1:dt.N["distribution"],l=1:dt.N["customer"]], sum(ukl[k,l,m] for m=1:dt.Mkl[k,l]) <= 1);    
    ########### constraint 11 #############
    @constraint(scnd,[i=1:dt.N["supplier"],j=1:dt.N["plant"],m=1:dt.Mij[i,j]], sum(xij[i,j,m,p] for p=1:5) <= dt.bigM*uij[i,j,m]);
    @constraint(scnd,[j=1:dt.N["plant"],k=1:dt.N["distribution"],m=1:dt.Mjk[j,k]] ,sum(xjk[j,k,m,p] for p=1:5) <= dt.bigM*ujk[j,k,m]);
    @constraint(scnd,[k=1:dt.N["distribution"],l=1:dt.N["customer"],m=1:dt.Mkl[k,l]], sum(xkl[k,l,m,p] for p=1:5) <= dt.bigM*ukl[k,l,m]);
    # ########### constraint 12 #############
    @constraint(scnd,[i=1:dt.N["supplier"], j=1:dt.N["plant"], m=1:dt.Mij[i,j]], sum(xij[i,j,m,p] for p=1:5) >= dt.Vij[m]*uij[i,j,m] );
    @constraint(scnd,[j=1:dt.N["plant"], k=1:dt.N["distribution"], m=1:dt.Mjk[j,k]], sum(xjk[j,k,m,p] for p=1:5) >= dt.Vjk[m]*ujk[j,k,m]);
    # ########### constraint 13-14 #############
    @constraint(scnd,sum(y[j,t] for j=1:dt.N["plant"] for t=1:2) <= dt.upl);
    @constraint(scnd,sum(y[j,t] for j=dt.N["plant"]+1:dt.N["distribution"]+dt.N["plant"] for t=1:2) <= dt.udc);
    return scnd
end
function SCND_MIP()
    ##########################  Mathematical model  #########################
    scnd = vModel(CPLEX.Optimizer); 
    # optimizer_with_attributes(
    #         CPLEX.Optimizer,
    #         "CPX_PARAM_EPGAP" => 1e-8);
    set_silent(scnd);
    MOI.set(scnd, MOI.NumberOfThreads(), 1);
    @variable(scnd, y[1:dt.N["plant"]+dt.N["distribution"],1:2], Bin);
    @variable(scnd, uij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:2], Bin);
    @variable(scnd, ujk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:2], Bin);
    @variable(scnd, ukl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:2], Bin);
    ############
    @variable(scnd, 0<= xij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:dt.Mij[i,j],1:5] );
    @variable(scnd, 0<= xjk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:dt.Mjk[j,k],1:5] );
    @variable(scnd, 0<= xkl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:dt.Mkl[k,l],1:5] );
    @variable(scnd, 0<= h[1:dt.N["plant"]+dt.N["distribution"],1:5,1:2] );
    exg = AffExpr(0);
    for i=1:dt.N["supplier"]
        for j=1:dt.N["plant"]
            add_to_expression!(exg, 10000*uij[i,j,1]);
        end
    end
    for j=1:dt.N["plant"]
        for k=1:dt.N["distribution"]
            add_to_expression!(exg, 10000*ujk[j,k,1]);
        end
    end
    for k=1:dt.N["distribution"]
        for l=1:dt.N["customer"]
            add_to_expression!(exg,10000*ukl[k,l,1]);
        end
    end
    @addobjective(scnd, Min, sum(dt.c[j][t]*y[j,t] for j=1:1:dt.N["plant"]+dt.N["distribution"] for t=1:2) + exg +
        sum(dt.N["vcs"][i][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) +
        sum(dt.vij[i][j][m][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
        sum(dt.vjk[j][k][m][p]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
        sum(dt.vkl[k][l][m][p]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5)+
        sum(dt.N["vcp"][j][2*(p-1)+t]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
        sum(dt.N["vcd"][k][2*(p-1)+t]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2));
    @addobjective(scnd, Min, sum(dt.b[i,p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
        sum(dt.N["vep"][j][2*(p-1)+t]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
        sum(dt.N["ved"][k][2*(p-1)+t]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2)+
        sum(dt.rij[i][j][m]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
        sum(dt.rjk[j][k][m]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
        sum(dt.rkl[k][l][m]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5));
    ######### constraint 3 #############
    @constraint(scnd,[j=1:dt.N["plant"],p=1:5], sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j]) == sum(xjk[j,k,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]));
    @constraint(scnd,[k=1:dt.N["distribution"],p=1:5], sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k]) == sum(xkl[k,l,m,p] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l]));
    ########### constraint 4-6 #############
    @constraint(scnd,[j=1:dt.N["plant"],p=1:5], sum(h[j,p,:]) == sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j]));
    @constraint(scnd,[k=1:dt.N["distribution"],p=1:5], sum(h[k+dt.N["plant"],p,:] ) == sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k]));
    @constraint(scnd, [l=1:dt.N["customer"],p=1:5], sum(xkl[k,l,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mkl[k,l]) >= dt.d[l][p]);
    ########### constraint 7-9 #############
    @constraint(scnd,[i=1:dt.N["supplier"]], sum(xij[i,j,m,p] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) <= dt.N["cas"][i]);
    @constraint(scnd,[j=1:dt.N["plant"]+dt.N["distribution"], t=1:2], sum(h[j,:,t]) <=[dt.N["cap"];dt.N["cad"]][j]*y[j,t]);
    @constraint(scnd,[j=1:dt.N["plant"]+dt.N["distribution"]], sum(y[j,:]) <= 1);
    ############ constraint 10 #############
    @constraint(scnd,[i=1:dt.N["supplier"],j=1:dt.N["plant"]], sum(uij[i,j,m] for m=1:dt.Mij[i,j]) <= 1);
    @constraint(scnd,[j=1:dt.N["plant"],k=1:dt.N["distribution"]], sum(ujk[j,k,m] for m=1:dt.Mjk[j,k]) <= 1);
    @constraint(scnd,[k=1:dt.N["distribution"],l=1:dt.N["customer"]], sum(ukl[k,l,m] for m=1:dt.Mkl[k,l]) <= 1);
    ########### constraint 11 #############
    @constraint(scnd,[i=1:dt.N["supplier"],j=1:dt.N["plant"],m=1:dt.Mij[i,j]], sum(xij[i,j,m,p] for p=1:5) <= dt.bigM*uij[i,j,m]);
    @constraint(scnd,[j=1:dt.N["plant"],k=1:dt.N["distribution"],m=1:dt.Mjk[j,k]] ,sum(xjk[j,k,m,p] for p=1:5) <= dt.bigM*ujk[j,k,m]);
    @constraint(scnd,[k=1:dt.N["distribution"],l=1:dt.N["customer"],m=1:dt.Mkl[k,l]], sum(xkl[k,l,m,p] for p=1:5) <= dt.bigM*ukl[k,l,m]);
    # ########### constraint 12 #############
    @constraint(scnd,[i=1:dt.N["supplier"], j=1:dt.N["plant"], m=1:dt.Mij[i,j]], sum(xij[i,j,m,p] for p=1:5) >= dt.Vij[m]*uij[i,j,m] );
    @constraint(scnd,[j=1:dt.N["plant"], k=1:dt.N["distribution"], m=1:dt.Mjk[j,k]], sum(xjk[j,k,m,p] for p=1:5) >= dt.Vjk[m]*ujk[j,k,m]);
    # ########### constraint 13-14 #############
    @constraint(scnd,sum(y[j,t] for j=1:dt.N["plant"] for t=1:2) <= dt.upl);
    @constraint(scnd,sum(y[j,t] for j=dt.N["plant"]+1:dt.N["distribution"]+dt.N["plant"] for t=1:2) <= dt.udc);

    return scnd
end
function lexobj1()
    ##########################  Mathematical model  #########################
    lex = Model(CPLEX.Optimizer); set_silent(lex);
    MOI.set(lex, MOI.NumberOfThreads(), 1);
    # set_optimizer_attribute(lex, "CPXPARAM_TimeLimit", TL)
    @variable(lex, y[1:dt.N["plant"]+dt.N["distribution"],1:2], Bin);
    @variable(lex, 0<= uij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:2] <= 1);
    @variable(lex, 0<= ujk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:2] <= 1);
    @variable(lex, 0<= ukl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:2] <= 1);
    ############
    @variable(lex, 0<= xij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:dt.Mij[i,j],1:5] );
    @variable(lex, 0<= xjk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:dt.Mjk[j,k],1:5] );
    @variable(lex, 0<= xkl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:dt.Mkl[k,l],1:5] );
    @variable(lex, 0<= h[1:dt.N["plant"]+dt.N["distribution"],1:5,1:2] );
    exg = AffExpr(0);
    for i=1:dt.N["supplier"]
        for j=1:dt.N["plant"]
            add_to_expression!(exg, 10000*uij[i,j,1]);
        end
    end
    for j=1:dt.N["plant"]
        for k=1:dt.N["distribution"]
            add_to_expression!(exg, 10000*ujk[j,k,1]);
        end
    end
    for k=1:dt.N["distribution"]
        for l=1:dt.N["customer"]
            add_to_expression!(exg,10000*ukl[k,l,1]);
        end
    end

    @objective(lex, Min, sum(dt.c[j][t]*y[j,t] for j=1:1:dt.N["plant"]+dt.N["distribution"] for t=1:2) + exg +
            sum(dt.N["vcs"][i][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) +
            sum(dt.vij[i][j][m][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
            sum(dt.vjk[j][k][m][p]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
            sum(dt.vkl[k][l][m][p]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5)+
            sum(dt.N["vcp"][j][2*(p-1)+t]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
            sum(dt.N["vcd"][k][2*(p-1)+t]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2));
    ########### constraint 3 #############
    @constraint(lex,[j=1:dt.N["plant"],p=1:5], sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j]) == sum(xjk[j,k,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]));
    @constraint(lex,[k=1:dt.N["distribution"],p=1:5], sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k]) == sum(xkl[k,l,m,p] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l]));
    ########### constraint 4-6 #############
    @constraint(lex,[j=1:dt.N["plant"],p=1:5], sum(h[j,p,:]) == sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j]));
    @constraint(lex,[k=1:dt.N["distribution"],p=1:5], sum(h[k+dt.N["plant"],p,:] ) == sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k]));
    @constraint(lex, [l=1:dt.N["customer"],p=1:5], sum(xkl[k,l,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mkl[k,l]) >= dt.d[l][p]);
    ########### constraint 7-9 #############
    @constraint(lex,[i=1:dt.N["supplier"]], sum(xij[i,j,m,p] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) <= dt.N["cas"][i]);
    @constraint(lex,[j=1:dt.N["plant"]+dt.N["distribution"], t=1:2], sum(h[j,:,t]) <=[dt.N["cap"];dt.N["cad"]][j]*y[j,t]);
    @constraint(lex,[j=1:dt.N["plant"]+dt.N["distribution"]], sum(y[j,:]) <= 1);
    ########### constraint 10 #############
    @constraint(lex,[i=1:dt.N["supplier"],j=1:dt.N["plant"]], sum(uij[i,j,m] for m=1:dt.Mij[i,j]) <= 1);
    @constraint(lex,[j=1:dt.N["plant"],k=1:dt.N["distribution"]], sum(ujk[j,k,m] for m=1:dt.Mjk[j,k]) <= 1);
    @constraint(lex,[k=1:dt.N["distribution"],l=1:dt.N["customer"]], sum(ukl[k,l,m] for m=1:dt.Mkl[k,l]) <= 1);  
    ########### constraint 11 #############
    @constraint(lex,[i=1:dt.N["supplier"],j=1:dt.N["plant"],m=1:dt.Mij[i,j]], sum(xij[i,j,m,p] for p=1:5) <= dt.bigM*uij[i,j,m]);
    @constraint(lex,[j=1:dt.N["plant"],k=1:dt.N["distribution"],m=1:dt.Mjk[j,k]] ,sum(xjk[j,k,m,p] for p=1:5) <= dt.bigM*ujk[j,k,m]);
    @constraint(lex,[k=1:dt.N["distribution"],l=1:dt.N["customer"],m=1:dt.Mkl[k,l]], sum(xkl[k,l,m,p] for p=1:5) <= dt.bigM*ukl[k,l,m]);
    # ########### constraint 12 #############
    @constraint(lex,[i=1:dt.N["supplier"], j=1:dt.N["plant"], m=1:dt.Mij[i,j]], sum(xij[i,j,m,p] for p=1:5) >= dt.Vij[m]*uij[i,j,m] );
    @constraint(lex,[j=1:dt.N["plant"], k=1:dt.N["distribution"], m=1:dt.Mjk[j,k]], sum(xjk[j,k,m,p] for p=1:5) >= dt.Vjk[m]*ujk[j,k,m]);
    # ########### constraint 13-14 #############
    @constraint(lex,sum(y[j,t] for j=1:dt.N["plant"] for t=1:2) <= dt.upl);
    @constraint(lex,sum(y[j,t] for j=dt.N["plant"]+1:dt.N["distribution"]+dt.N["plant"] for t=1:2) <= dt.udc);
    return lex
end
function lexobj2()
    ##########################  Mathematical model  #########################
    lex = Model(CPLEX.Optimizer); set_silent(lex);
    MOI.set(lex, MOI.NumberOfThreads(), 1);
    # set_optimizer_attribute(lex, "CPXPARAM_TimeLimit", TL)
    @variable(lex, y[1:dt.N["plant"]+dt.N["distribution"],1:2], Bin);
    @variable(lex, 0<= uij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:2] <= 1);
    @variable(lex, 0<= ujk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:2] <= 1);
    @variable(lex, 0<= ukl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:2] <= 1);
    @variable(lex, 0<= xij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:dt.Mij[i,j],1:5] );
    @variable(lex, 0<= xjk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:dt.Mjk[j,k],1:5] );
    @variable(lex, 0<= xkl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:dt.Mkl[k,l],1:5] );
    @variable(lex, 0<= h[1:dt.N["plant"]+dt.N["distribution"],1:5,1:2] );
    @objective(lex, Min, sum(dt.b[i,p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
            sum(dt.N["vep"][j][2*(p-1)+t]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
            sum(dt.N["ved"][k][2*(p-1)+t]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2)+
            sum(dt.rij[i][j][m]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
            sum(dt.rjk[j][k][m]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
            sum(dt.rkl[k][l][m]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5));
    
    @constraint(lex,[j=1:dt.N["plant"],p=1:5], sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j]) == sum(xjk[j,k,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]));
    @constraint(lex,[k=1:dt.N["distribution"],p=1:5], sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k]) == sum(xkl[k,l,m,p] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l]));
    ########### constraint 4-6 #############
    @constraint(lex,[j=1:dt.N["plant"],p=1:5], sum(h[j,p,:]) == sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j]));
    @constraint(lex,[k=1:dt.N["distribution"],p=1:5], sum(h[k+dt.N["plant"],p,:] ) == sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k]));
    @constraint(lex, [l=1:dt.N["customer"],p=1:5], sum(xkl[k,l,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mkl[k,l]) >= dt.d[l][p]);
    ########### constraint 7-9 #############
    @constraint(lex,[i=1:dt.N["supplier"]], sum(xij[i,j,m,p] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) <= dt.N["cas"][i]);
    @constraint(lex,[j=1:dt.N["plant"]+dt.N["distribution"], t=1:2], sum(h[j,:,t]) <=[dt.N["cap"];dt.N["cad"]][j]*y[j,t]);
    @constraint(lex,[j=1:dt.N["plant"]+dt.N["distribution"]], sum(y[j,:]) <= 1);
    ########### constraint 10 #############
    @constraint(lex,[i=1:dt.N["supplier"],j=1:dt.N["plant"]], sum(uij[i,j,m] for m=1:dt.Mij[i,j]) <= 1);
    @constraint(lex,[j=1:dt.N["plant"],k=1:dt.N["distribution"]], sum(ujk[j,k,m] for m=1:dt.Mjk[j,k]) <= 1);
    @constraint(lex,[k=1:dt.N["distribution"],l=1:dt.N["customer"]], sum(ukl[k,l,m] for m=1:dt.Mkl[k,l]) <= 1);  
    ########### constraint 11 #############
    @constraint(lex,[i=1:dt.N["supplier"],j=1:dt.N["plant"],m=1:dt.Mij[i,j]], sum(xij[i,j,m,p] for p=1:5) <= dt.bigM*uij[i,j,m]);
    @constraint(lex,[j=1:dt.N["plant"],k=1:dt.N["distribution"],m=1:dt.Mjk[j,k]] ,sum(xjk[j,k,m,p] for p=1:5) <= dt.bigM*ujk[j,k,m]);
    @constraint(lex,[k=1:dt.N["distribution"],l=1:dt.N["customer"],m=1:dt.Mkl[k,l]], sum(xkl[k,l,m,p] for p=1:5) <= dt.bigM*ukl[k,l,m]);
    # ########### constraint 12 #############
    @constraint(lex,[i=1:dt.N["supplier"], j=1:dt.N["plant"], m=1:dt.Mij[i,j]], sum(xij[i,j,m,p] for p=1:5) >= dt.Vij[m]*uij[i,j,m] );
    @constraint(lex,[j=1:dt.N["plant"], k=1:dt.N["distribution"], m=1:dt.Mjk[j,k]], sum(xjk[j,k,m,p] for p=1:5) >= dt.Vjk[m]*ujk[j,k,m]);
    # ########### constraint 13-14 #############
    @constraint(lex,sum(y[j,t] for j=1:dt.N["plant"] for t=1:2) <= dt.upl);
    @constraint(lex,sum(y[j,t] for j=dt.N["plant"]+1:dt.N["distribution"]+dt.N["plant"] for t=1:2) <= dt.udc);
    return lex
end
function getobjval(model)
    y = value.(model[:y]);
    uij = value.(model[:uij]);
    ujk = value.(model[:ujk]);
    ukl = value.(model[:ukl]);
    xij = value.(model[:xij]);
    xjk = value.(model[:xjk]);
    xkl = value.(model[:xkl]);
    h = value.(model[:h]);
    exg = 0
    for i=1:dt.N["supplier"]
        for j=1:dt.N["plant"]
            exg = exg+10000*uij[i,j,1];
        end
    end
    for j=1:dt.N["plant"]
        for k=1:dt.N["distribution"]
            exg = exg + 10000*ujk[j,k,1] #add_to_expression!(exg, );
        end
    end
    for k=1:dt.N["distribution"]
        for l=1:dt.N["customer"]
            exg = exg + 10000*ukl[k,l,1] #add_to_expression!(exg,);
        end
    end

    obj1 = sum(dt.c[j][t]*y[j,t] for j=1:1:dt.N["plant"]+dt.N["distribution"] for t=1:2) + exg +
        sum(dt.N["vcs"][i][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) +
        sum(dt.vij[i][j][m][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
        sum(dt.vjk[j][k][m][p]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
        sum(dt.vkl[k][l][m][p]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5)+
        sum(dt.N["vcp"][j][2*(p-1)+t]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
        sum(dt.N["vcd"][k][2*(p-1)+t]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2)
    
    obj2 = sum(dt.b[i,p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
        sum(dt.N["vep"][j][2*(p-1)+t]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
        sum(dt.N["ved"][k][2*(p-1)+t]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2)+
        sum(dt.rij[i][j][m]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
        sum(dt.rjk[j][k][m]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
        sum(dt.rkl[k][l][m]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5)
    return [obj1,obj2]#,dval
end
function dominated(x,P)
    st = false
    for k=1:length(P)
        if all( x .>= P[k])
            st=true; break
        else
            continue
        end
    end
    return st
end
function NDfilter(P,Pobj)
    copysol = Dict(); copyobj = Dict();
    for i=1:length(Pobj)
        copysol[i] = P[i]
        copyobj[i] = Pobj[i]
    end
    for i=1:length(Pobj)-1
        for j=i+1:length(Pobj)
            if all(Pobj[i] .>= Pobj[j]) == true #dominated by PF[j]
                copyobj[i]=nothing; copysol[i]=nothing; break
            elseif all(Pobj[j] .>= Pobj[i]) == true
                copyobj[j]=nothing; copysol[j]=nothing;
            end
        end
    end

    finalsol = filter!(a->a!=nothing, collect(values(copysol)))
    finalobj = filter!(a->a!=nothing, collect(values(copyobj)))
    return finalsol,finalobj
end
function NDsetfilter(ndset)
    removal = []
    for i=1:length(ndset)-1
        for j=i+1:length(ndset)
            if all(ndset[i].val .> ndset[j].val) == true #dominated by PF[j]
                push!(removal,i); @goto nextnd  
            elseif all(ndset[j].val .> ndset[i].val) == true
                push!(removal,j)
            end
        end
        @label nextnd
    end
    sort!(unique!(removal))
    filteredND = deleteat!(ndset, removal)
    return filteredND
end
function SortingSol(P,Pobj)
    copysol = Dict(); copyobj = Dict();
    for i=1:length(Pobj)
        copysol[i] = P[i]
        copyobj[i] = Pobj[i]
    end
    for i=1:length(Pobj)-1
        for j=i+1:length(Pobj)
            if all(Pobj[i] .>= Pobj[j]) == true #dominated by PF[j]
                copyobj[i]=0; copysol[i]=0; break
            elseif all(Pobj[j] .>= Pobj[i]) == true
                copyobj[j]=0; copysol[j]=0;
            end
        end
    end
    sortedsol = filter!(a->a!=0, collect(values(copysol)))
    sortedobj = filter!(a->a!=0, collect(values(copyobj)))
    df = DataFrame(X=sortedsol,Y=sortedobj);
    sort!(df,[:Y])
    return df
end
function FP_Model()
    model = Model(CPLEX.Optimizer);
    # optimizer_with_attributes(
    #     CPLEX.Optimizer,
    #     "CPX_PARAM_EPGAP" => 1e-4);
    set_silent(model)
    MOI.set(model, MOI.NumberOfThreads(), 1)
    @variable(model, y[1:dt.N["plant"]+dt.N["distribution"],1:2], Bin);
    @variable(model, uij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:2], Bin);
    @variable(model, ujk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:2], Bin);
    @variable(model, ukl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:2], Bin);
    ############
    @variable(model, 0<= xij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:dt.Mij[i,j],1:5] );
    @variable(model, 0<= xjk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:dt.Mjk[j,k],1:5] );
    @variable(model, 0<= xkl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:dt.Mkl[k,l],1:5] );
    @variable(model, 0<= h[1:dt.N["plant"]+dt.N["distribution"],1:5,1:2] );
    
    ######### constraint 3 #############
    @constraint(model,[j=1:dt.N["plant"],p=1:5], sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j]) == sum(xjk[j,k,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]));
    @constraint(model,[k=1:dt.N["distribution"],p=1:5], sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k]) == sum(xkl[k,l,m,p] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l]));
    ########### constraint 4-6 #############
    @constraint(model,[j=1:dt.N["plant"],p=1:5], sum(h[j,p,:]) == sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j]));
    @constraint(model,[k=1:dt.N["distribution"],p=1:5], sum(h[k+dt.N["plant"],p,:] ) == sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k]));
    @constraint(model, [l=1:dt.N["customer"],p=1:5], sum(xkl[k,l,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mkl[k,l]) >= dt.d[l][p]);
    ########### constraint 7-9 #############
    @constraint(model,[i=1:dt.N["supplier"]], sum(xij[i,j,m,p] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) <= dt.N["cas"][i]);
    @constraint(model,[j=1:dt.N["plant"]+dt.N["distribution"], t=1:2], sum(h[j,:,t]) <=[dt.N["cap"];dt.N["cad"]][j]*y[j,t]);
    @constraint(model,[j=1:dt.N["plant"]+dt.N["distribution"]], sum(y[j,:]) <= 1);
    ############ constraint 10 #############
    @constraint(model,[i=1:dt.N["supplier"],j=1:dt.N["plant"]], sum(uij[i,j,m] for m=1:dt.Mij[i,j]) <= 1);
    @constraint(model,[j=1:dt.N["plant"],k=1:dt.N["distribution"]], sum(ujk[j,k,m] for m=1:dt.Mjk[j,k]) <= 1);
    @constraint(model,[k=1:dt.N["distribution"],l=1:dt.N["customer"]], sum(ukl[k,l,m] for m=1:dt.Mkl[k,l]) <= 1);
    ########### constraint 11 #############
    @constraint(model,[i=1:dt.N["supplier"],j=1:dt.N["plant"],m=1:dt.Mij[i,j]], sum(xij[i,j,m,p] for p=1:5) <= dt.bigM*uij[i,j,m]);
    @constraint(model,[j=1:dt.N["plant"],k=1:dt.N["distribution"],m=1:dt.Mjk[j,k]] ,sum(xjk[j,k,m,p] for p=1:5) <= dt.bigM*ujk[j,k,m]);
    @constraint(model,[k=1:dt.N["distribution"],l=1:dt.N["customer"],m=1:dt.Mkl[k,l]], sum(xkl[k,l,m,p] for p=1:5) <= dt.bigM*ukl[k,l,m]);
    # ########### constraint 12 #############
    @constraint(model,[i=1:dt.N["supplier"], j=1:dt.N["plant"], m=1:dt.Mij[i,j]], sum(xij[i,j,m,p] for p=1:5) >= dt.Vij[m]*uij[i,j,m] );
    @constraint(model,[j=1:dt.N["plant"], k=1:dt.N["distribution"], m=1:dt.Mjk[j,k]], sum(xjk[j,k,m,p] for p=1:5) >= dt.Vjk[m]*ujk[j,k,m]);
    # ########### constraint 13-14 #############
    @constraint(model,sum(y[j,t] for j=1:dt.N["plant"] for t=1:2) <= dt.upl);
    @constraint(model,sum(y[j,t] for j=dt.N["plant"]+1:dt.N["distribution"]+dt.N["plant"] for t=1:2) <= dt.udc);
    return model
end
function LP_Model()
    ##########################  Mathematical model  #########################
    lp = Model(CPLEX.Optimizer); set_silent(lp);
    # optimizer_with_attributes(
    #         CPLEX.Optimizer,
    #         "CPX_PARAM_EPGAP" => 1e-4);
    @variable(lp, 0<= y[1:dt.N["plant"]+dt.N["distribution"],1:2] <= 1);
    @variable(lp, 0<= uij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:2] <= 1)
    @variable(lp, 0<= ujk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:2] <= 1)
    @variable(lp, 0<= ukl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:2] <= 1)
    ############
    @variable(lp, 0<= xij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:dt.Mij[i,j],1:5] );
    @variable(lp, 0<= xjk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:dt.Mjk[j,k],1:5] );
    @variable(lp, 0<= xkl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:dt.Mkl[k,l],1:5] );
    @variable(lp, 0<= h[1:dt.N["plant"]+dt.N["distribution"],1:5,1:2] );
    exg = AffExpr(0);
    for i=1:dt.N["supplier"]
        for j=1:dt.N["plant"]
            add_to_expression!(exg, 10000*uij[i,j,1]);
        end
    end
    for j=1:dt.N["plant"]
        for k=1:dt.N["distribution"]
            add_to_expression!(exg, 10000*ujk[j,k,1]);
        end
    end
    for k=1:dt.N["distribution"]
        for l=1:dt.N["customer"]
            add_to_expression!(exg,10000*ukl[k,l,1]);
        end
    end
    @constraint(lp,[j=1:dt.N["plant"],p=1:5], sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j]) == sum(xjk[j,k,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]));
    @constraint(lp,[k=1:dt.N["distribution"],p=1:5], sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k]) == sum(xkl[k,l,m,p] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l]));
    ########### constraint 4-6 #############
    @constraint(lp,[j=1:dt.N["plant"],p=1:5], sum(h[j,p,:]) == sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j]));
    @constraint(lp,[k=1:dt.N["distribution"],p=1:5], sum(h[k+dt.N["plant"],p,:] ) == sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k]));
    @constraint(lp, [l=1:dt.N["customer"],p=1:5], sum(xkl[k,l,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mkl[k,l]) >= dt.d[l][p]);
    ########### constraint 7-9 #############
    @constraint(lp,[i=1:dt.N["supplier"]], sum(xij[i,j,m,p] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) <= dt.N["cas"][i]);
    @constraint(lp,[j=1:dt.N["plant"]+dt.N["distribution"], t=1:2], sum(h[j,:,t]) <=[dt.N["cap"];dt.N["cad"]][j]*y[j,t]);
    @constraint(lp,[j=1:dt.N["plant"]+dt.N["distribution"]], sum(y[j,:]) <= 1);
    ########### constraint 10 #############
    @constraint(lp,[i=1:dt.N["supplier"],j=1:dt.N["plant"]], sum(uij[i,j,m] for m=1:dt.Mij[i,j]) <= 1);
    @constraint(lp,[j=1:dt.N["plant"],k=1:dt.N["distribution"]], sum(ujk[j,k,m] for m=1:dt.Mjk[j,k]) <= 1);
    @constraint(lp,[k=1:dt.N["distribution"],l=1:dt.N["customer"]], sum(ukl[k,l,m] for m=1:dt.Mkl[k,l]) <= 1);    
    ########### constraint 11 #############
    @constraint(lp,[i=1:dt.N["supplier"],j=1:dt.N["plant"],m=1:dt.Mij[i,j]], sum(xij[i,j,m,p] for p=1:5) <= dt.bigM*uij[i,j,m]);
    @constraint(lp,[j=1:dt.N["plant"],k=1:dt.N["distribution"],m=1:dt.Mjk[j,k]] ,sum(xjk[j,k,m,p] for p=1:5) <= dt.bigM*ujk[j,k,m]);
    @constraint(lp,[k=1:dt.N["distribution"],l=1:dt.N["customer"],m=1:dt.Mkl[k,l]], sum(xkl[k,l,m,p] for p=1:5) <= dt.bigM*ukl[k,l,m]);
    # ########### constraint 12 #############
    @constraint(lp,[i=1:dt.N["supplier"], j=1:dt.N["plant"], m=1:dt.Mij[i,j]], sum(xij[i,j,m,p] for p=1:5) >= dt.Vij[m]*uij[i,j,m] );
    @constraint(lp,[j=1:dt.N["plant"], k=1:dt.N["distribution"], m=1:dt.Mjk[j,k]], sum(xjk[j,k,m,p] for p=1:5) >= dt.Vjk[m]*ujk[j,k,m]);
    # ########### constraint 13-14 #############
    @constraint(lp,sum(y[j,t] for j=1:dt.N["plant"] for t=1:2) <= dt.upl);
    @constraint(lp,sum(y[j,t] for j=dt.N["plant"]+1:dt.N["distribution"]+dt.N["plant"] for t=1:2) <= dt.udc);
    return lp
end
function flip(x_h,j,e)
    if x_h[e[j]]==1
        x_h[e[j]] = 0
    else #make a feasible solution
        if isodd(e[j])
            x_h[e[j]+1] = 0
        else
            x_h[e[j]-1] = 0
        end
        x_h[e[j]] = 1

    end
    return x_h
end
function flipoper(Tabu,x_t,x_r)
    e = sortperm(abs.(x_t-x_r),rev=true)
    xi = []
    x_h = copy(x_r)
    j = 1
    M=length(x_t) #
    while j<=M && xi==[]
        x_h = flip(x_h,j,e)
        if x_h ∉ Tabu
            xi=x_h
        else
            j+=1
        end
    end
    if xi==[]
        while j<=M
            x_h=copy(x_r)
            Num = Int64(rand(ceil(length(x_r)/2):length(x_r)-1))
            R = StatsBase.sample(1:M,Num, replace=false)
            for r in R
                x_h = flip(x_h,r,e)
                if x_h ∉ Tabu
                    xi = x_h
                end
            end
            j+=1
        end
    end
    return xi
end
function FP_FBcheck(model,yr,iter)        
    yr = reshape(yr,(Int(len[1]/2),2))
    if isodd(iter)==true
        exg = AffExpr(0)
        for i=1:dt.N["supplier"]
            for j=1:dt.N["plant"]
                exg+10000*model[:uij][i,j,1];
            end
        end
        for j=1:dt.N["plant"]
            for k=1:dt.N["distribution"]
                exg + 10000*model[:ujk][j,k,1] #add_to_expression!(exg, );
            end
        end
        for k=1:dt.N["distribution"]
            for l=1:dt.N["customer"]
                exg + 10000*model[:ukl][k,l,1] #add_to_expression!(exg,);
            end
        end
        @objective(model, Min, 
            sum(dt.c[j][t]*model[:y][j,t] for j=1:1:dt.N["plant"]+dt.N["distribution"] for t=1:2) + exg +
            sum(dt.N["vcs"][i][p]*model[:xij][i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) +
            sum(dt.vij[i][j][m][p]*model[:xij][i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
            sum(dt.vjk[j][k][m][p]*model[:xjk][j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
            sum(dt.vkl[k][l][m][p]*model[:xkl][k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5)+
            sum(dt.N["vcp"][j][2*(p-1)+t]*model[:h][j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
            sum(dt.N["vcd"][k][2*(p-1)+t]*model[:h][k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2)) 
    else
        @objective(model, Min, 
            sum(dt.b[i,p]*model[:xij][i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
            sum(dt.N["vep"][j][2*(p-1)+t]*model[:h][j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
            sum(dt.N["ved"][k][2*(p-1)+t]*model[:h][k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2)+
            sum(dt.rij[i][j][m]*model[:xij][i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
            sum(dt.rjk[j][k][m]*model[:xjk][j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
            sum(dt.rkl[k][l][m]*model[:xkl][k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5))
    end

    JuMP.fix.(model[:y],yr; force=true)
    optimize!(model)
    st = termination_status(model)
    if st == MOI.OPTIMAL
        return true
    else
        return false
    end
end
function fbsearch(yr)#,u1r,u2r,u3r) #solveLP
    yr = reshape(yr, Int(len[1]/2),2)
    idy_0 = findall(k->k==0, yr)
    idy_1 = findall(k->k==1, yr)
    @objective( dist, Min, sum(dist[:y][i] for i in idy_0) + sum(1-(dist[:y][j]) for j in idy_1))
    optimize!(dist)
    if termination_status(dist) == MOI.OPTIMAL
        return reshape(value.(dist[:y]),len[1])  #JuMP.value.(dist[:y])
    else
        return 0
    end
end
function fbsearch2(yr,u1r,u2r,u3r) #solveLP
    yr = reshape(yr, Int(len[1]/2),2)
    u1r = reshape(u1r, dt.N["supplier"],dt.N["plant"],2)
    u2r = reshape(u2r, dt.N["plant"],dt.N["distribution"],2)
    u3r = reshape(u3r, dt.N["distribution"],dt.N["customer"],2)

    idy_0 = findall(k->k==0, yr)
    idy_1 = findall(k->k==1, yr)
    idu1_0 = findall(k->k==0, u1r)
    idu1_1 = findall(k->k==1, u1r)
    idu2_0 = findall(k->k==0, u2r)
    idu2_1 = findall(k->k==1, u2r)
    idu3_0 = findall(k->k==0, u3r)
    idu3_1 = findall(k->k==1, u3r)
    @objective( dist, Min, sum(dist[:y][i] for i in idy_0) + sum(1-(dist[:y][j]) for j in idy_1) +
        sum(dist[:uij][i] for i in idu1_0) + sum(1-(dist[:uij][j]) for j in idu1_1)+
        sum(dist[:ujk][i] for i in idu2_0) + sum(1-(dist[:ujk][j]) for j in idu2_1)+
        sum(dist[:ukl][i] for i in idu3_0) + sum(1-(dist[:ukl][j]) for j in idu3_1))
    optimize!(dist)
    if termination_status(dist) == MOI.OPTIMAL
        nx = value.(all_variables(dist))
        return nx[1:len[1]],nx[1+len[1]:len[1]+len[2]],nx[1+sum(len[i] for i=1:2):sum(len[i] for i=1:3)],nx[1+sum(len[i] for i=1:3):sum(len[i] for i=1:4)]
        # JuMP.value.(dist[:y]),JuMP.value.(dist[:uij]),JuMP.value.(dist[:ujk]),JuMP.value.(dist[:ukl])
    else
        return 0,0,0,0
    end
end
function FP(candX,len,TL,Max_iter)
    X = []; PF =[]; Tabu = []; newsol = 0; Y = [];
    candlist = StatsBase.sample(1:length(candX), length(candX), replace=false)
    t0 = time(); 
    for k in candlist  #&& time() - t0 < TL 
        x_t = candX[k]
        yt = x_t[1:len[1]]; yr = copy(yt)
        SearchDone = false; iter=0;
        while iter < Max_iter && SearchDone == false && time() - t0 < TL
            yid = findall(p->p>0.2,yt);
            for j=1:len[1]
                if j in yid
                    yr[j]=1
                else
                    yr[j]=0
                end
            end
            if FP_FBcheck(fpmodel,yr,iter) == true #,u1r,u2r,u3r)
                sol = value.(all_variables(fpmodel)); ndp = getobjval(fpmodel)
                if sol ∉ X  && dominated(ndp,collect(values(PF)))==false
                    push!(X,sol); push!(Y,yr); push!(PF,ndp) #PF[k] = ndp
                    newsol+=1; SearchDone = true; 
                    # deleteat!(candlist, k);
                end
            else
                if yr ∈ Tabu
                    yr = flipoper(Y,yt,yr); # u1r = flipoper(U1,u1t,u1r); u2r = flipoper(U2,u2t,u2r); u3r = flipoper(U3,u3t,u3r)
                    if yr == [] # if any(i->i==[], [yr,u1r,u2r,u3r])
                        SearchDone = true; 
                        # deleteat!(candlist, k);
                       
                    else
                        if FP_FBcheck(fpmodel,yr,iter) == true #,u1r,u2r,u3r)
                            sol = value.(all_variables(fpmodel)); ndp = getobjval(fpmodel)
                            if sol ∉ X && dominated(ndp,collect(values(PF)))==false
                                push!(X,sol); push!(Y,yr); push!(PF,ndp) #PF[k] = ndp
                                newsol+=1; SearchDone = true;
                            end
                        end
                    end
                end
                if time()-t0 >= TL
                    break
                end
            
                if SearchDone == false
                    push!(Tabu,yr)
                    yt = fbsearch(yr);
                    # println("new lp sol")
                    if yt==0  #when there's no new feasible lp sol
                        SearchDone = true
                    end
                end
            end
			iter+=1
        end
    end
    return X,PF
end
function FPplus(dvar,Y_N,len,TL,Max_iter) 
    X = copy(dvar); PF = copy(Y_N); Y = []; IGPair=[]; Tabu = []; t0=time();
    # U1 = []; U2= []; U3 = []; newsol=0; 
    while time()-t0 < TL && length(IGPair)<(length(PF)*(length(PF)-1))
        I,G = StatsBase.sample(1:length(X), 2, replace=false)
        x1 = X[I][1:sum(len[i] for i=1:4)]; x2 = X[G][1:sum(len[i] for i=1:4)];
        λ = round(rand(Float64, 1)[1]; digits=1)
        x_t = x1*λ + x2*(1-λ);
        # x_t = x1*.5 + x2*.5;
        yt = x_t[1:len[1]]; 
        u1t = x_t[1+len[1]:len[1]+len[2]]; u2t = x_t[1+sum(len[i] for i=1:2):sum(len[i] for i=1:3)]; u3t = x_t[1+sum(len[i] for i=1:3):sum(len[i] for i=1:4)]
        SearchDone = false; iter=0;
        while [I,G]∉IGPair && iter<Max_iter && SearchDone == false
            # x_r = round.(Int,x_t);
            yr = round.(Int, yt); 
            u1r = round.(Int, u1t);u2r = round.(Int, u2t); u3r = round.(Int, u3t);

            if FP_FBcheck(fpmodel,yr,iter) == true
                sol = value.(all_variables(fpmodel)); ndp = getobjval(fpmodel)
                if sol ∉ X  && dominated(ndp,PF)==false
                    push!(X,sol); push!(PF,ndp)
                    # push!(Y,yr); push!(U1,u1r); push!(U2,u2r); push!(U3,u3r);
                    # newsol+=1; 
                    SearchDone = true
                    # println("rounding worked")
                end
            else
                if yr ∈ Tabu #[yr;u1r;u2r;u3r] ∈ Tabu
                    yr = flipoper(Y,yt,yr);
                    # u1r = flipoper(U1,u1t,u1r); u2r = flipoper(U2,u2t,u2r); u3r = flipoper(U3,u3t,u3r)
                    if yr == [] # if any(i->i==[], [yr,u1r,u2r,u3r])
                        SearchDone = true;
                        # println("flip failed")
                    else
                        if FP_FBcheck(fpmodel,yr,iter) == true
                            sol = value.(all_variables(fpmodel)); ndp = getobjval(fpmodel)
                            if sol ∉ X && dominated(ndp,PF)==false
                                push!(X,sol); push!(PF,ndp)
                                # newsol+=1; push!(Y,yr); push!(U1,u1r); push!(U2,u2r); push!(U3,u3r);
                                SearchDone = true;
                                # println("flip worked")
                            end
                        end
                    end
                end
                if time()-t0 >= TL
                    break
                end
                if SearchDone == false
                    push!(Tabu,yr) #[yr;u1r;u2r;u3r]
                    yt,u1t,u2t,u3t = fbsearch2(yr,u1r,u2r,u3r) # yt = fbsearch(yr) 
                    # if yt==0
                    if any(i->i==0, [yt,u1t,u2t,u3t])  #when there's no new feasible lp sol
                        # println("no solution")
                        SearchDone = true
                    end
                end
            end
			iter+=1
        end
        push!(IGPair,[I,G])

    end
    return X,PF#,IGPair,newsol
end
function PR_Model()
    model = Model(CPLEX.Optimizer);
    set_silent(model)
    MOI.set(model, MOI.NumberOfThreads(), 1)
    @variable(model, y[1:dt.N["plant"]+dt.N["distribution"],1:2], Bin);
    @variable(model, uij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:2], Bin);
    @variable(model, ujk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:2], Bin);
    @variable(model, ukl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:2], Bin);
    ############
    @variable(model, 0<= xij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:dt.Mij[i,j],1:5] );
    @variable(model, 0<= xjk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:dt.Mjk[j,k],1:5] );
    @variable(model, 0<= xkl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:dt.Mkl[k,l],1:5] );
    @variable(model, 0<= h[1:dt.N["plant"]+dt.N["distribution"],1:5,1:2] );
    
    ######### constraint 3 #############
    @constraint(model,[j=1:dt.N["plant"],p=1:5], sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j]) == sum(xjk[j,k,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]));
    @constraint(model,[k=1:dt.N["distribution"],p=1:5], sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k]) == sum(xkl[k,l,m,p] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l]));
    ########### constraint 4-6 #############
    @constraint(model,[j=1:dt.N["plant"],p=1:5], sum(h[j,p,:]) == sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j]));
    @constraint(model,[k=1:dt.N["distribution"],p=1:5], sum(h[k+dt.N["plant"],p,:] ) == sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k]));
    @constraint(model, [l=1:dt.N["customer"],p=1:5], sum(xkl[k,l,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mkl[k,l]) >= dt.d[l][p]);
    ########### constraint 7-9 #############
    @constraint(model,[i=1:dt.N["supplier"]], sum(xij[i,j,m,p] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) <= dt.N["cas"][i]);
    @constraint(model,[j=1:dt.N["plant"]+dt.N["distribution"], t=1:2], sum(h[j,:,t]) <=[dt.N["cap"];dt.N["cad"]][j]*y[j,t]);
    @constraint(model,[j=1:dt.N["plant"]+dt.N["distribution"]], sum(y[j,:]) <= 1);
    ############ constraint 10 #############
    @constraint(model,[i=1:dt.N["supplier"],j=1:dt.N["plant"]], sum(uij[i,j,m] for m=1:dt.Mij[i,j]) <= 1);
    @constraint(model,[j=1:dt.N["plant"],k=1:dt.N["distribution"]], sum(ujk[j,k,m] for m=1:dt.Mjk[j,k]) <= 1);
    @constraint(model,[k=1:dt.N["distribution"],l=1:dt.N["customer"]], sum(ukl[k,l,m] for m=1:dt.Mkl[k,l]) <= 1);
    ########### constraint 11 #############
    @constraint(model,[i=1:dt.N["supplier"],j=1:dt.N["plant"],m=1:dt.Mij[i,j]], sum(xij[i,j,m,p] for p=1:5) <= dt.bigM*uij[i,j,m]);
    @constraint(model,[j=1:dt.N["plant"],k=1:dt.N["distribution"],m=1:dt.Mjk[j,k]] ,sum(xjk[j,k,m,p] for p=1:5) <= dt.bigM*ujk[j,k,m]);
    @constraint(model,[k=1:dt.N["distribution"],l=1:dt.N["customer"],m=1:dt.Mkl[k,l]], sum(xkl[k,l,m,p] for p=1:5) <= dt.bigM*ukl[k,l,m]);
    # ########### constraint 12 #############
    @constraint(model,[i=1:dt.N["supplier"], j=1:dt.N["plant"], m=1:dt.Mij[i,j]], sum(xij[i,j,m,p] for p=1:5) >= dt.Vij[m]*uij[i,j,m] );
    @constraint(model,[j=1:dt.N["plant"], k=1:dt.N["distribution"], m=1:dt.Mjk[j,k]], sum(xjk[j,k,m,p] for p=1:5) >= dt.Vjk[m]*ujk[j,k,m]);
    # ########### constraint 13-14 #############
    @constraint(model,sum(y[j,t] for j=1:dt.N["plant"] for t=1:2) <= dt.upl);
    @constraint(model,sum(y[j,t] for j=dt.N["plant"]+1:dt.N["distribution"]+dt.N["plant"] for t=1:2) <= dt.udc);
    return model
end
function SelectIG_by_Distance(ndp,IGPair) #,len
    dist = [];  connect = []
    for i=1:length(ndp.Y)-1 
        if i==1 
            d = sqrt(sum((ndp.Y[i].-ndp.Y[i+1]).^2))
            push!(dist,d); push!(connect,2)
        else
            d1 = sqrt(sum((ndp.Y[i].-ndp.Y[i-1]).^2))
            d2 = sqrt(sum((ndp.Y[i].-ndp.Y[i+1]).^2))
            push!(dist, max(d1,d2))
            if d1 >= d2
                push!(connect,1)
            else
                push!(connect,2)
            end
        end
    end
    while count(==(0), dist)< round(Int,(length(ndp.Y)-1)*1) #sum(dist) != 0
        id = findmax(dist)[2]
        if connect[id] == 1
            I,G = ndp.Y[id],ndp.Y[id-1] 
        else
            I,G = ndp.Y[id],ndp.Y[id+1]
        end
        if [I,G] in IGPair #|| round.(ndp.X[id][1:len[1]]) == round.(ndp.X[id+1][1:len[1]])
            dist[id] = 0
        else
            return I,G
        end
    end
    return 0,0
end
function createAllNB2(SI,dif)
    SI = round.(SI)
    neibour = []; 
    for i in dif
        cpSI = copy(SI)
        if cpSI[i] == 1
            cpSI[i] = 0
        else # make a feasibiel neighbour
            if isodd(i)
                cpSI[i+1] = 0
            else # i is even
                cpSI[i-1] = 0
            end
            cpSI[i] = 1
        end
        push!(neibour, cpSI);# push!(neiobj, binobj(cpSI,bvar))
    end
    return neibour
end
function nextSI(neibour,neiobj)#SI
    if length(neibour) == 1  #if there is one candiate sol
        return neibour[1][1:len[1]] #,neibour[1]
    else length(neibour) > 1 # if there are multiple candiates, check the improved ratio
        # neiobj = [getobjval(prmodel)]
        ratiotb = zeros(length(neiobj),length(neiobj[1]))
        for i=1:length(neiobj)
            ratiotb[i,:] = neiobj[i]#./SIobj
        end
        ranktb = zeros(length(neiobj),length(neiobj[1]))
        for i=1:length(neiobj[1])
            ranktb[:,i] = tiedrank(ratiotb[:,i])
        end
        ranksum = [sum(ranktb[i,:]) for i=1:length(neiobj)]
        mostimp = findall(x-> x == minimum(ranksum), ranksum)
        k = rand(mostimp)[1]
        return neibour[k][1:len[1]]#,neiobj[k]
    end
end
function PR_FBcheck_w(model,weight,yr)
    yr = reshape(yr,(Int(len[1]/2),2))
    exg = AffExpr(0)
    for i=1:dt.N["supplier"]
        for j=1:dt.N["plant"]
            exg+10000*model[:uij][i,j,1];
        end
    end
    for j=1:dt.N["plant"]
        for k=1:dt.N["distribution"]
            exg + 10000*model[:ujk][j,k,1] #add_to_expression!(exg, );
        end
    end
    for k=1:dt.N["distribution"]
        for l=1:dt.N["customer"]
            exg + 10000*model[:ukl][k,l,1] #add_to_expression!(exg,);
        end
    end
    @objective(model, Min, weight[1]*(sum(dt.c[j][t]*model[:y][j,t] for j=1:1:dt.N["plant"]+dt.N["distribution"] for t=1:2) + exg +
        sum(dt.N["vcs"][i][p]*model[:xij][i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) +
        sum(dt.vij[i][j][m][p]*model[:xij][i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
        sum(dt.vjk[j][k][m][p]*model[:xjk][j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
        sum(dt.vkl[k][l][m][p]*model[:xkl][k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5)+
        sum(dt.N["vcp"][j][2*(p-1)+t]*model[:h][j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
        sum(dt.N["vcd"][k][2*(p-1)+t]*model[:h][k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2)) 
        +
        weight[2]*(sum(dt.b[i,p]*model[:xij][i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
        sum(dt.N["vep"][j][2*(p-1)+t]*model[:h][j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
        sum(dt.N["ved"][k][2*(p-1)+t]*model[:h][k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2)+
        sum(dt.rij[i][j][m]*model[:xij][i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
        sum(dt.rjk[j][k][m]*model[:xjk][j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
        sum(dt.rkl[k][l][m]*model[:xkl][k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5))
    );
    JuMP.fix.(model[:y], yr; force=true)
    optimize!(model)
    if has_values(model) == true #
    # if termination_status(model) == MOI.OPTIMAL
        return true
    else
        return false
    end
end
############## Initialising nondominated line segments
struct node
    val::Array; arm::String; #dom::String x::Float64; y::Float64
end
function Newnodes(lsg)
    nodes = []
    if length(lsg) == 1
        push!(nodes, node(lsg[1], "null"))
    else
        for i=1:length(lsg)
            if i == 1
                connect = "R" 
            elseif i == length(lsg)
                connect = "L" 
            else
                connect = "LR" 
            end
            push!(nodes, node(lsg[i], connect)) 
        end
    end
    return nodes
end
function SectionQueue(ndset,newsg)
    # When newsg has multiple elements
    # que = Matrix(undef,length(newsg),length(ndset)) 
    # for i=1:length(newsg)
    #     for j=1:length(ndset)
    #         if newsg[i].val[1] >= ndset[j].val[1] && newsg[i].val[2] >= ndset[j].val[2] 
    #             que[i,j] = "u"
    #         elseif newsg[i].val[1] <= ndset[j].val[1] && newsg[i].val[2] <= ndset[j].val[2] 
    #             que[i,j] = "d"
    #         elseif newsg[i].val[1] >= ndset[j].val[1] && newsg[i].val[2] <= ndset[j].val[2] 
    #             que[i,j] = "r"
    #         else
    #             que[i,j] = "l"
    #         end
    #     end
    # end

    # single new point (newsg)
    que = []
    for j=1:length(ndset)
        if newsg[1].val[1] >= ndset[j].val[1] && newsg[1].val[2] >= ndset[j].val[2] 
            push!(que, "u")
        elseif newsg[1].val[1] <= ndset[j].val[1] && newsg[1].val[2] <= ndset[j].val[2] 
            push!(que, "d")
        elseif newsg[1].val[1] >= ndset[j].val[1] && newsg[1].val[2] <= ndset[j].val[2] 
            push!(que, "r")
        else
            push!(que, "l")
        end
    end
    return que
end
function dominance_count(x,P)
    st = 0
    for k=1:length(P)
        if all( x .<= P[k])
            st+=1; 
        else
            continue
        end
    end
    return st
end
function ClassifyND(ndset,dsol)
    Larms = findall(x-> ndset[x].arm == "L", 1:length(ndset))
    start = 0;  ndom = []; dom = [];

    if isempty(Larms) == false 
        for l in Larms 
            refset = ndset[1+start:l]
            polygon = [refset[i].val for i=1:length(refset)]
            insert!(polygon, 1, [refset[1].val[1], 99^30])
            push!(polygon, [99^30, refset[end].val[2]],[99^30, refset[1].val[2]])
            edges = zeros(Int,length(polygon),2)
            for i=1:length(polygon)-1
                edges[i,1] = i; edges[i,2] = i+1
            end
            edges[end,1] = length(polygon); edges[end,2] = 1
            stat = inpoly2(dsol, polygon, edges, atol = 1e-1)
            tmpnd = findall(x-> stat[x,:] == [0,0], 1:length(dsol))
            push!(ndom, tmpnd)
            union!(dom, filter( x-> x ∉ tmpnd, 1:length(dsol)))
            start = l
        end
    else
        l = length(ndset)
        refset = ndset[1+start:l]
        polygon = [refset[i].val for i=1:length(refset)]
        insert!(polygon, 1, [refset[1].val[1], 99^30])
        push!(polygon, [99^30, refset[end].val[2]],[99^30, refset[1].val[2]])
        edges = zeros(Int,length(polygon),2)
        for i=1:length(polygon)-1
            edges[i,1] = i; edges[i,2] = i+1
        end
        edges[end,1] = length(polygon); edges[end,2] = 1
        stat = inpoly2(dsol, polygon, edges, atol = 1e-1)
        tmpnd = findall(x-> stat[x,:] == [0,0], 1:length(dsol))
        push!(ndom, tmpnd)
        union!(dom, filter( x-> x ∉ tmpnd, 1:length(dsol)))
    end

    ndom2 = []
    if length(ndom)>1
        ndom2 = ndom[1]
        for i=2:length(ndom)
            ndom2 = intersect(ndom2, ndom[i])
        end
    elseif length(ndom)==1 
        ndom2 = ndom[1]
    end

    return ndom2,dom
end
function SolveLPdicho(model,sol,ndp0,ndlist)
    ndp = round.(ndp0)
    y′ = reshape(sol[1:len[1]], Int(len[1]/2),2)
    uij′ = reshape(sol[1+len[1]:sum(len[i] for i=1:2)], dt.N["supplier"], dt.N["plant"],2)
    ujk′ = reshape(sol[1+sum(len[i] for i=1:2):sum(len[i] for i=1:3)], dt.N["plant"], dt.N["distribution"],2)
    ukl′ = reshape(sol[1+sum(len[i] for i=1:3):sum(len[i] for i=1:4)], dt.N["distribution"],dt.N["customer"],2)
    
    JuMP.fix.(model[:y], y′; force = true)
    JuMP.fix.(model[:uij], uij′; force = true)
    JuMP.fix.(model[:ujk], ujk′; force = true)
    JuMP.fix.(model[:ukl], ukl′; force = true)
    vSolve(model, 5, method=:dicho, verbose=false)
    res = getvOptData(model);
    if res.Y_N == []  # no trade-off solusions
        dichoy = [ndp0]
        ddicho = Dict(ndp0 => sol)
    else
        for i=1:length(res.Y_N)
            res.Y_N[i] = round.(res.Y_N[i])
        end

        # if ndp ∈ res.Y_N
        #     id = findall(x -> x == ndp, res.Y_N)
        #     deleteat!(res.X_E,id); deleteat!(res.Y_N,id)
        # end
        dichodf = DataFrame(X = res.X_E, Y = res.Y_N) 
        sort!(dichodf, [:Y])
        # dichoy = dichodf.Y
        # for i=1:length(dichodf.Y)
        #     ddicho[dichodf.Y[i]] = dichodf.X[i]
        # end        
        
        dichoy = []; ddicho = Dict()
        for i=1:length(dichodf.Y)
            if dichodf.Y[i] ∉ ndlist
                push!(dichoy,dichodf.Y[i]); push!(ddicho, dichodf.Y[i] => dichodf.X[i])
            end
        end

    end
    return dichoy,ddicho
end
function InitialNDset(df,model)  
    pair = Dict()
    dsol1,dict1 = SolveLPdicho(model,df.X[1],df.Y[1],[])
    ndset = Newnodes(dsol1)
    merge!(pair,dict1)

    for l=2:length(df.Y)
        ndlist = [ndset[i].val for i=1:length(ndset)]
        if df.Y[l] ∉ ndlist
            dsol,dict= SolveLPdicho(model,df.X[l],df.Y[l],ndlist)
            if isempty(dsol) == true
                # println("dsol is empty")
                # break 
            else
                merge!(pair,dict)
                ndicho,ddicho = ClassifyND(ndset,dsol)
                newsg = Newnodes(dsol)
                ndomset,domset = ClassifyND(newsg,ndlist)
                ndset = UpdateNDset(ddicho,ndicho,domset,ndomset,newsg,ndset) 
            end
        end
    end

    ndf = DataFrame(X = [], Y = [])
    for y in [ndset[i].val for i=1:length(ndset)] 
        if haskey(pair,y)
            push!(ndf.X,pair[y])
            push!(ndf.Y,y) 
        end
    end
    return ndf,ndset
end
function UpdateNDset(ddicho,ndicho,domset,ndomset,newsg_og,ndset_og)
    ndset = copy(ndset_og)
    if length(newsg_og) == 1    
        # println(newsg_og ,"dsol is a point")                               # new solution is a point
        que = SectionQueue(ndset_og,newsg_og)
        nw = newsg_og[1] 
        if "u" ∈ que                                     # new solution is dominated
            return ndset
        elseif all(x->x == "d", que)                     # new solution dominate the whole ndset
            ndset = nw; return ndset
        elseif all(x->x == "l", que)      
            insert!(ndset, 1, nw); return ndset               
        elseif all(x->x == "r", que)      
            insert!(ndset, length(ndset)+1, nw); return ndset
        end
    
        d_id = findall(x->x =="d", que[1,:])
        if length(d_id) == 1 # # new solution is dominating one ndpoint 
            case = 0
            if ndset[d_id[1]].arm == "null" # one dominated point without an arm
                replace!(ndset, ndset[d_id[1]] => nw)
            else  # one dominated point with an arm
                if d_id[1] == 1                 # 1st ndset point is dominated
                    lsgx = LineSegment(nw.val,[nw.val[1],ndset[2].val[2]])
                    pjx = LazySets.isdisjoint(LineSegment(ndset[1].val,ndset[2].val), lsgx, true)
                    case = 1
                elseif d_id[1] == length(ndset) # last ndset point is dominated
                    # pjy = LazySets.isdisjoint(LineSegment(ndset[end-1].val,ndset[end].val), Line(nw.val,[1,0.]),true)
                    pjy = LazySets.isdisjoint(LineSegment(ndset[end-1].val,ndset[end].val), LineSegment(nw.val,[nw.val[1],ndset[end-1].val[2]]),true)
                    case = 2
                else                            # midle ndset point is dominated
                    lsgx = LineSegment(nw.val,[ndset[d_id[1]-1].val[1],nw.val[2]])
                    pjx = LazySets.is_intersection_empty(LineSegment(ndset[d_id[1]-1].val,ndset[d_id[1]].val), lsgx, true)
                    pjy = LazySets.isdisjoint(LineSegment(ndset[end-1].val,ndset[end].val), LineSegment(nw.val,[nw.val[1],ndset[end-1].val[2]]),true)
                    case = 3
                end
                if case == 1
                    replace!(ndset, ndset[d_id[1]] => nw); insert!(ndset, d_id[1]+1, node(pjx[2], "R"))
                elseif case == 2
                    replace!(ndset, ndset[end] => node(pjy[2], "L")); insert!(ndset, length(ndset), nw)
                else 
                    replace!(ndset, ndset[d_id[1]] => node(pjy[2], "L")); 
                    insert!(ndset, d_id[1]+1, nw)
                    insert!(ndset, d_id[1]+2, node(pjx[2], "R"))
                end
            end
        elseif length(d_id) >1 # new solution is dominating some ndpoints/linesegment
            case = 0  
            if d_id[1] == 1
                if d_id[end].arm == "R" || d_id[end].arm == "LR"
                    case = 1 # case1: only pjx calculation needed 
                else
                    case =0 # case0: replace dominated lsg and return ndset
                end
            else
                if ndset[d_id[1]].arm == "L" || ndset[d_id[1]].arm == "LR"
                    if d_id[end].arm == "R" || d_id[end].arm == "LR"
                        case = 3
                    else
                        case = 2
                    end
                else
                    if d_id[end].arm == "R" || d_id[end].arm == "LR"
                        case = 3
                    else
                        case =0
                    end
                end
            end
            deleteat!(ndset, d_id);
            if case == 0
                insert!(ndset, d_id[1], nw)
                return ndset
            elseif case ==1 
                if d_id[1] == 1
                    lsgx = LineSegment(nw.val,[nw.val[1],ndset[2].val[2]])
                    pjx = LazySets.is_intersection_empty(LineSegment(ndset[1].val,ndset[2].val), lsgx, true)
                else
                    lsgx = LineSegment(nw.val,[ndset[d_id[1]-1].val[1],nw.val[2]])
                    pjx = LazySets.is_intersection_empty(LineSegment(ndset[d_id[1]-1].val,ndset[d_id[1]].val), lsgx, true)
                end
                insert!(ndset, d_id[1], nw); insert!(ndset, d_id[1]+1, node(pjx[2], "R"))
            elseif case == 2
                pjy = LazySets.isdisjoint(LineSegment(ndset[d_id[1]-1].val,ndset[d_id[1]].val), LineSegment(nw.val,[nw.val[1],ndset[d_id[1]-1].val[2]]),true)
                insert!(ndset, d_id[1], node(pjy[2], "L")); insert!(ndset, d_id[1]+1, nw)
            else 
                insert!(ndset, d_id[1], node(pjy[2], "L")); 
                insert!(ndset, d_id[1]+1, nw)
                insert!(ndset, d_id[1]+2, node(pjx[2], "R"))
            end
        end

        rl = findall(i-> que[i] == "r" && que[i+1] == "l", 1:length(que))
        if isempty(rl) == false 
            if ndset[rl[1]].arm == "R" # connected linesegment
                if abs((ndset[rl[1]].val[2]-ndset[rl[1]+1].val[2])/(ndset[rl[1]].val[1]-ndset[rl[1]+1].val[1])) > abs((ndset[rl[1]].val[2]-nw.val[2])/(ndset[rl[1]].val[1]-nw.val[1])) #new sol is nondominated
                    if ndset[rl[1]].arm != "null"  # divide the lsg into two parts and insert the new sol
                        seg = [ndset[rl[1]].val,ndset[rl[1]+1].val]
                        p1 = LineSegment(seg[1],seg[2])
                        p2 = LineSegment(nw.val,[nw.val[1],seg[1][2]])
                        interpt1 = LazySets.isdisjoint(p1,p2,true)
                        p3 = LineSegment(nw.val,[seg[2][1],nw.val[2]])
                        interpt2 = LazySets.isdisjoint(p1,p3,true)
                        
                        #insert interpt1,nw,interpt2
                        insert!( ndset, rl[1]+1, node(interpt1[2], "L") )
                        insert!( ndset, rl[1]+2, node(nw.val, "null") )
                        insert!( ndset, rl[1]+3, node(interpt2[2], "R") )
                    else
                        insert!( ndset, rl[1]+1, node(nw.val, "null") ) #insert the new sol btw two points
                    end
                end
            else                # new solution between two points
                if dominated(nw.val,[ndset[rl[1]].val,ndset[rl[1]+1].val]) == false
                    insert!( ndset, rl[1]+1, node(nw.val, "null") )  # add new sol only if it's nondominated
                end
            end
        end
    else
        newsg = copy(newsg_og)    
        removal = []; addnw = []; ndstart = 1;
        for k =1:length(newsg)-1
            if k ∈ ddicho && k+1 ∈ ddicho
            else
                for i=ndstart:length(ndset)-1
                    # @show (k,newsg[k])
                    # @show (i,ndset[i])
                    if i ∈ ndomset && i+1 ∈ ndomset && k+1 ∉ ndicho
                        ndstart = i+1
                        # println("all nd nondominated. next: ", ndstart)
                    elseif i ∈ domset && i+1 ∈ domset
                        # println("all nd domset")
                        push!(removal,ndset[i],ndset[i+1])
                        ndstart = i+1 #MODI
                    else
                        if ndset[i].arm == "R" || ndset[i].arm =="LR"
                            nwline = LazySets.LineSegment(newsg[k].val,newsg[k+1].val)
                            ndline = LazySets.LineSegment(ndset[i].val,ndset[i+1].val)
                            inter = LazySets.isdisjoint(nwline,ndline,true)                    
                            fourpt0 = sort!([newsg[k],newsg[k+1],ndset[i],ndset[i+1]], by = x-> x.val[1]) 
                            fourpt = copy(fourpt0)
                            if isempty(inter[2]) # two line segments are disjoint
        
                                # 2nd point projected to yaxis
                                if fourpt0[2].val ∈ nwline
                                    l3 = ndline
                                else
                                    l3 = nwline
                                end
                                lsg2y = LazySets.LineSegment([fourpt0[2].val[1],0],[fourpt0[2].val[1],99^30])
                                pj2y = is_intersection_empty(lsg2y,l3,true) 
                                # pj2y = LazySets.isdisjoint(Line(fourpt0[2].val,[0.,1]), l3, true)
                                if isempty(pj2y[2]) # no intersection
                                    dom3 = dominated(fourpt0[3].val,[fourpt0[1].val,fourpt0[2].val])
                                    dom4 = dominated(fourpt0[4].val,[fourpt0[1].val,fourpt0[2].val])
        
                                    if (dom3,dom4) == (1,1) #nd points all dominated
                                        # println("dom = 1,1")
                                        push!(removal, fourpt0[3],fourpt0[4])
                                    elseif (dom3,dom4) == (0,0)  
                                        # println("dom = 0,0")
                                        if fourpt0[3].val ∈ nwline
                                            if i == length(ndset)-1 
                                                push!(addnw, newsg[k], newsg[k+1])
                                                @goto Nextnw
                                            else
                                                ndstart = i+1
                                            end
                                        else
                                            # println("add all newsg: ", newsg[k], newsg[k+1])
                                            push!(addnw, newsg[k], newsg[k+1])
                                            @goto Nextnw
                                        end
                                    elseif (dom3,dom4) == (0,1) 
                                        # println("dom = 0,1")
                                    else # (dom3,dom4) == (1,0)  #pt3 dominated
                                        # println("dom = 1,0")
                                        lsg2x = LazySets.LineSegment(fourpt0[2].val,[fourpt0[4].val[1],fourpt0[2].val[2]])
                                        pj2x = is_intersection_empty(lsg2x,l3,true)
                                        if isempty(pj2x[2])
                                            # println("pj2x is empty")
                                            push!(removal, fourpt0[3],fourpt0[4])
                                            @goto Nextnw
                                        else
                                            if fourpt0[3].val ∈ nwline
                                                push!(removal, newsg[k])
                                                push!(addnw, newsg[k+1]) 
                                                # println("pt3 on the newline, pj2x added: ", pj2x[2])
                                                if i+1 == length(ndset)
                                                    push!(addnw, node(pj2x[2], "R"))
                                                    @goto Nextnw
                                                end
                                                newsg[k] = node(pj2x[2], "R")
                                                ndstart = i+1
                                                
                                            else
                                                # println("add two new points")
                                                push!(removal, ndset[i])
                                                push!(addnw, newsg[k], node(newsg[k+1].val, "L") , node(pj2x[2], "R"))
                                                insert!(ndset, i+1, node(pj2x[2], "R"))
                                                @goto Nextnw
                                            end
                                        end
                                    end
                                elseif pj2y[2][2] > fourpt0[2].val[2]
                                    if fourpt0[2].val ∈ nwline
                                        # println("add pj 2y to ndset: ",pj2y[2])
                                        # println("add nwsg to ndset: ",newsg[k].val)
                                        # replace!(ndset, ndset[i+1] => node(ndset[i+1].val, "R"))
                                        push!(ndset, node(pj2y[2], "L"), node(newsg[k].val, "R"))
                                    else
                                        # println("only add pj 2y to ndset: ",pj2y[2])
                                        push!(addnw, newsg[k], node(pj2y[2], "L"))
                                    end
                                elseif pj2y[2][2] < fourpt0[2].val[2]
                                    # println(" no add from pj 2y")
                                    if fourpt0[2].val ∈ nwline
                                        push!(removal, newsg[k])
                                    else
                                        push!(addnw, newsg[k])
                                        push!(removal, ndset[i])
                                    end
                                end
                        
                                # 3rd point projected to xaxis
                                if fourpt0[3].val[2] > fourpt0[4].val[2] 
                                    if fourpt0[3].val ∈ nwline
                                        l3 = ndline
                                    else
                                        l3 = nwline
                                    end
                                    lsg3x = LazySets.LineSegment([fourpt0[1].val[1],fourpt0[3].val[2]],[fourpt0[4].val[1],fourpt0[3].val[2]])
                                    pj3x = is_intersection_empty(lsg3x,l3,true)
                                    
                                    if isempty(pj3x[2])  # no intersection 
                                        ndstart = i+1 
                                    elseif pj3x[2][1] < fourpt0[3].val[1]
        
                                        if fourpt0[3].val ∈ nwline
                                            if k+1 != length(newsg)
                                                # println("pj3x no need. remove newsg: ", fourpt0[3].val)
                                                # push!(removal, fourpt0[3].val)  #MODI
                                                push!(removal, newsg[k])
                                            end
                                        else
                                            # println("remove ndset ",ndset[i+1]) 
                                            # println("add newsg ", newsg[k+1])
                                            push!(addnw, newsg[k+1]); push!(removal, ndset[i+1])
                                        end
                                        ndstart = i+1
                                    elseif pj3x[2][1] > fourpt0[3].val[1]
                                        
                                        if fourpt0[3].val ∈ nwline
                                            # println("on the nwline added from proj 3: ", pj3x[2])
                                            push!(addnw, node(newsg[k+1].val, "L"))
                                            insert!(ndset, i+1, node(pj3x[2], "R"))
                                            @goto Nextnw
                                        else
                                            # println("not on the line. add newsg[k+1]: ", newsg[k+1])
                                            newsg[k] = node(pj3x[2], "R")
                                            push!(addnw, newsg[k+1])
                                            ndstart = i+1
                                        end
                                    end
                                else
                                    if fourpt0[3].val ∈ nwline
                                        # println("add from proj 3", newsg[k+1])
                                        # println("rm from proj 3", ndset[i+1])
                                        push!(addnw,newsg[k+1]); push!(removal, ndset[i+1])
                                    else
                                        # println("rm from proj 3")
                                        push!(removal, newsg[k+1])
                                    end
                                    ndstart = i+1; @goto Nextnw
                                end
                            else #two line segments intersect
                                # line segments share the point => calculate with three points
                                if inter[2] ∈ [newsg[k].val,newsg[k+1].val,ndset[i].val,ndset[i+1].val]
                                    threept0 = copy(fourpt0)
                                    interid = findall(x->x.val==inter[2],threept0)[1]
                                    deleteat!(threept0, interid)
                                    threept = copy(threept0)
                                    if interid == 1
                                        if threept0[2].val[2] < threept0[3].val[2]
                                            if threept0[3] ∈ newsg
                                                @goto Nextnw
                                            else
                                                push!(addnw, newsg[k+1]); push!(removal, ndset[i+1])
                                                ndstart = i+1
                                            end
                                        else
                                            lsg2x = LazySets.LineSegment([threept0[1].val[1],threept0[2].val[2]],[threept0[3].val[1],threept0[2].val[2]])
                                            pj2x = LazySets.isdisjoint(lsg2x, LineSegment(threept0[1].val,threept0[3].val), true)
                                            if pj2x[2][1] < threept0[2].val[1] # no intersection
                                                if threept0[2] ∈ newsg
                                                    @goto Nextnw
                                                else
                                                    push!(removal, ndset[i+1])
                                                    push!(addnw, newsg[k+1])
                                                    ndstart = i+1
                                                end
                                            else
                                                if isempty(pj2x[2])
                                                    @show( newsg[k:k+1], ndset[i:i+1])
                                                    println("FIX!! interid =1, pj2x empty")
                                                end
                                                if threept0[2] ∈ newsg
                                                    push!(addnw, node(newsg[k+1].val, "L"))
                                                    insert!(ndset, i+1, node(pj2x[2], "R"))
                                                    ndstart = i+1
                                                    @goto Nextnw    
                                                else
                                                    replace!(newsg, newsg[k] => node(pj2x[2], "R"))
                                                    push!(addnw, newsg[k+1], node(pj2x[2], "R"))
                                                    ndstart = i+1
                                                end
                                            end
                                        end
                                    elseif interid == 3
                                        if threept0[1].val[2] < threept0[2].val[2]
                                            if threept0[1] ∈ newsg
                                                push!(addnw, newsg[k]); push!(removal, ndset[i])
                                            end
                                        else
                                            lsg2y = LazySets.LineSegment([threept0[2].val[1],0],[fourpt0[2].val[1],99^30])
                                            pj2y = is_intersection_empty(lsg2y,LineSegment(threept0[1].val,threept0[3].val),true)
                                            # pj2y = LazySets.isdisjoint(Line(threept0[2].val, [0,1.]), Line(threept0[1].val,threept0[3].val), true)
                                            push!(addnw, newsg[k])
                                            if isempty(pj2y[2])
                                                # @show( newsg[k:k+1], ndset[i:i+1])
                                                # println("interid =3, pj2y empty")
                                            else
                                                push!(addnw, node(pj2y[2], "L"))
                                            end
                                        end
                                        ndstart = i+1 
                                        @goto Nextnw
                                    elseif interid == 2
                                        if threept0[2] ∈ newsg
                                            push!(addnw, newsg[k])
                                        else
                                            push!(addnw, newsg[k+1])
                                        end
                                        
                                        ndstart = i+1 
                                        @goto Nextnw
                                    end
        
                                # calculate with four points
                                else 
                                    if fourpt0[2].val ∈ nwline
                                        l3 = ndline
                                    else
                                        l3 = nwline
                                    end
                                    # 2nd point projected to yaxis
                                    lsg2y = LazySets.LineSegment([fourpt0[2].val[1],0],[fourpt0[2].val[1],99^30])
                                    pj2y = is_intersection_empty(lsg2y,l3,true) 
                                    if isempty(pj2y[2])
                                        # @show( newsg[k:k+1], ndset[i:i+1])
                                        println("FIX !! fourpt, pj2y empty")
                                    # else
                                        # println("added pj2y: ", pj2y[2])
                                    end  
                                    if pj2y[2][2] > fourpt0[2].val[2]
                                        if fourpt0[2].val ∈ nwline
                                            # fourpt[2].arm changed from "LR" to "R"
                                            # println("add pj2y: ",pj2y[2])
                                            # println("add pt2: ",fourpt[2])
                                            # println("intersection: ",inter[2])
                                            push!(addnw, node(pj2y[2], "L"), node(fourpt[2].val, "R"), node(inter[2], "LR"))
                                        else
                                            # println("add pt2: ",fourpt[2])
                                            push!(addnw, fourpt[2]) 
                                            if k ∉ ddicho 
                                                push!(addnw, node(pj2y[2], "L"))
                                            end
                                            push!(addnw, node(inter[2], "LR"))
                                        end
                                    else
                                        # println(" add inter: ", inter[2])
                                        if fourpt0[2].val ∈ nwline
                                            # println("rm newsg[k]: ", newsg[k])
                                            push!(removal, newsg[k])
                                            push!(addnw,  node(inter[2], "LR"),)
                                        else
                                            # println("rm ndset[i]: ", ndset[i])
                                            push!(removal, ndset[i])
                                            push!(addnw, newsg[k], node(inter[2], "LR"))
                                        end
                                    end
                                    
                                    #3rd point projected to xaxis
                                    #slope calculation
                                    nwdif = abs.(newsg[k].val-newsg[k+1].val)
                                    slop_newseg = abs(nwdif[2]/nwdif[1])
                                    nddif = abs.(ndset[i].val-ndset[i+1].val)
                                    slop_ndseg = abs(nddif[2]/nddif[1])
                                    if fourpt0[3].val ∈ nwline
                                        l3 = ndline
                                        pjslop = slop_newseg; theother_slop = slop_ndseg
                                    else
                                        l3 = nwline
                                        pjslop = slop_ndseg; theother_slop = slop_newseg
                                    end
                                    
                                    if fourpt0[3].val[2] > fourpt0[4].val[2] && abs(pjslop) > abs(theother_slop)  
                                        lsg3x = LazySets.LineSegment([fourpt0[1].val[1],fourpt0[3].val[2]],[fourpt0[4].val[1],fourpt0[3].val[2]])
                                        pj3x = is_intersection_empty(lsg3x,l3,true)
                                        
                                        if isempty(pj3x[2])
                                            # @show( newsg[k:k+1], ndset[i:i+1])
                                            println("FIX !! fourpt, pj3x empty")
                                            # @show ndset; 
                                            # @show newsg
                                        # else
                                        #     println("pj3x added: ", pj3x[2])    
                                        end  
                                        
                                        if fourpt0[3].val ∈ nwline
                                            # println("add newsg[k+1]: ", newsg[k+1])
                                            # println("replace ndset[i] with pj3x: ", pj3x[2])
                                            push!(addnw, newsg[k+1], node(pj3x[2], "R"))
                                            # insert!(ndset, i+1, node(pj3x[2], "R")); # ndstart = i+1 
                                            replace!(ndset, ndset[i] => node(pj3x[2], "R")) #MODI
                                            @goto Nextnw
                                        else
                                            push!(addnw, node(pj3x[2], "R"), newsg[k+1])
                                            newsg[k] = node(pj3x[2], "R")
                                            ndstart = i+1
                                        end
                                    elseif fourpt0[3].val[2] < fourpt0[4].val[2]
                                        if fourpt0[3].val ∈ nwline
                                            push!(removal, ndset[i+1])
                                            push!(addnw, newsg[k+1])
                                        end
                                        ndstart = i+1
                                        @goto Nextnw
                                    end
                                end
                            end
                        elseif ndset[i].arm =="L" #only Left arm => nextnode
                            if i != length(ndset) 
                                ndstart = i+1
                            else #last node
                                return ndset
                            end
                        else # ndset[i] is a "point"
                            # println("ndset is a point")
                            if ndset[i].val[1] < newsg[k].val[1]
                                dc = dominance_count(ndset[i].val,[newsg[j].val for j=k:k+1])
                                if dc == 0
                                    insert!(ndset, i+1, newsg[k])
                                    insert!(ndset, i+2, newsg[k+1])
                                elseif dc == 1
                                    lsgx = LineSegment(ndset[i].val,[newsg[k+1].val[1],ndset[i].val[2]])
                                    pjx = LazySets.is_intersection_empty( lsgx, LineSegment(newsg[k].val,newsg[k+1].val), true)
                                    if isempty(pjx[2])
                                        @show( newsg[k:k+1], ndset[i])
                                        println("FIX !! ndset point, pjx empty")
                                    end  
                                    insert!(ndset, i+1 , node(pjx[2], "LR")); insert!(ndset, i+2 , newsg[k+1])
                                end
                            elseif ndset[i].val[1] > newsg[k+1].val[1] && ndset[i].val[2] < newsg[k+1].val[2]
                                insert!(ndset, i, newsg[k]); insert!(ndset, i+1, newsg[k+1])
                            elseif newsg[k].val[1] < ndset[i].val[1] < newsg[k+1].val[1] 
                                ndslop = (ndset[i].val[2] - newsg[k+1].val[2])/(ndset[i].val[1] - newsg[k+1].val[1])
                                segslop = (newsg[k+1].val[2] - newsg[k].val[2])/(newsg[k].val[1] - newsg[k+1].val[1])
                                if ndslop < 0 && ndslop > segslop
                                    
                                    pjy = LazySets.isdisjoint( Line(ndset[i].val,[0,1.]), LineSegment(newsg[k].val,newsg[k+1].val), true)
                                    if isempty(pjy[2])
                                        @show( newsg[k:k+1], ndset[i])
                                        println("FIX !! ndset point-2, pjy empty")
                                    end  
                                    lsgx = LineSegment(ndset[i].val,[newsg[k+1].val[1],ndset[i].val[2]])
                                    pjx = LazySets.is_intersection_empty( lsgx, LineSegment(newsg[k].val,newsg[k+1].val), true)
                                    if isempty(pjx[2])
                                        @show( newsg[k:k+1], ndset[i])
                                        println("FIX !! ndset point-2, pjx empty")
                                    end  
                                    insert!(ndset, i, newsg[k]); insert!(ndset, i+1, node(pjy[2], "L")); 
                                    insert!(ndset, i+3, node(pjx[2], "R")); insert!(ndset, i+4, newsg[k])
                                elseif ndslop > 0
                                    if isempty(pjy[2])
                                        @show( newsg[k:k+1], ndset[i])
                                        println("FIX !! ndset point-3, pjy empty")
                                    end  
                                    pjy = LazySets.isdisjoint( Line(ndset[i].val,[0,1.]), LineSegment(newsg[k].val,newsg[k+1].val), true)
                                    insert!(ndset, i, newsg[k]); insert!(ndset, i+1, node(pjy[2], "L")); 
                                end
                            end
                        end
                    end
                end
            end
            @label Nextnw

        end
         #Post processing 
        unique!(addnw); unique!(removal)
        append!(ndset,addnw)
        append!(ndset,newsg_og[ndicho]); append!(ndset,ndset_og[ndomset])
        setdiff!(ndset,removal)
        setdiff!(ndset,ndset_og[domset])

        #Remove dominated points
        ndset = NDsetfilter(ndset)
        ndsort = DataFrame(v1 = [ndset[i].val[1] for i=1:length(ndset)], v2 = [ndset[i].val[2] for i=1:length(ndset)], p = [ndset[i].val for i=1:length(ndset)], arm = [ndset[i].arm for i=1:length(ndset)])
        sort!( ndsort, [:v1,:v2],rev = [false,true] )
        obj = [[ndsort.v1[i],ndsort.v2[i]] for i=1:length(ndsort.v1)]
        ndset = node.(obj, ndsort.arm)
     
        dup = findall(i-> i != 1, countmap(ndsort.p))
        for d in dup
            setdiff!(ndset,ndset[findall(i-> i == d, ndsort.p)[1:end-1]])
            delete!(ndsort, findall(i-> i == d, ndsort.p)[1:end-1])
        end

        Larm = findall(x -> ndset[x].arm == "LR" && ndset[x+1].arm == "R", 1:length(ndset)-1)
        Rarm = findall(x -> ndset[x-1].arm == "L" && ndset[x].arm == "LR", 2:length(ndset))

        for l in Larm
            replace!( ndset, ndset[l] => node(ndset[l].val,"L") )
        end
        for r in Rarm
            replace!( ndset, ndset[r] => node(ndset[r].val,"R") )
        end    
    end
    return ndset
end
function PRdicho(Xp,Yp,fsi0,PF0,len,TL,pc1)
    X = copy(Xp); Y = copy(Yp); PF = copy(PF0); fsi = copy(fsi0)
    IGPair=[]; infsi = []; #fsi = DataFrame(key=[],x=[],y=[]); #IGPair2=[]; Tabu = []; exploredSI = [];
    # println("#######  1st STAGE ##########")
    iter1 = 0; newsol1 = 0; iter2 = 0; newsol2 = 0; iter3 = 0; newsol3 = 0; t0 = time();
    while time()-t0 < TL*pc1
        nd = SortingSol(X,Y)
        I,G = SelectIG_by_Distance(nd,IGPair)
        # println("I G = ", I,"  ",G )
        if [I,G] == [0,0]
            allcomb = []
            for i=1:length(nd.Y)
                for j=1:length(nd.Y)
                    if i!=j
                        push!(allcomb, [nd.Y[i],nd.Y[j]])
                    end
                end
            end
            unused = setdiff(allcomb,IGPair)
            if isempty(unused) == false
                I,G = StatsBase.sample(unused, 1)[1] 
                Iid = findall(i->i==I, nd.Y)[1]; Gid = findall(i->i==G, nd.Y)[1]
                SI = round.(nd.X[Iid]); SG = round.(nd.X[Gid])
                iter2+=1
                # println("random IG selection")
            else
                # println("All IGpairs are used")
                @goto lsgCal
            end                      
        else
            # println("Diatance based I-G")
            Iid = findall(i->i==I, nd.Y)[1]; Gid = findall(i->i==G, nd.Y)[1]
            SI = round.(nd.X[Iid]); SG = round.(nd.X[Gid])
            iter1+=1
        end
        
        ct = 0; λ1 = abs(round.(I[2]-G[2])); λ2 = abs(round.(G[1]-I[1])); new_weight = [λ1,λ2]# nosol = 0; 
        while round.(SI[1:len[1]]) != round.(SG[1:len[1]]) && time()-t0 < TL*pc1 #&& nosol < 6 
            dif = findall(i-> SI[1:len[1]][i]!=SG[1:len[1]][i], 1:len[1])
            neibour = createAllNB2(SI[1:len[1]],dif)
            candSI =[]; candobj = []; 
            for l=1:length(neibour)
                if (time()-t0 >= TL)
                    @goto lsgCal
                end
                if neibour[l] ∈ fsi.key
                    id0 = findall(k -> k == neibour[l], fsi.key); id = rand(id0)[1]
                    sol = fsi.x[id]; ndp = fsi.y[id]
                    push!(candSI,sol); push!(candobj,ndp); ct+=1;
                    if ct >2
                        # println(" rep in fsikey")
                        @goto NextPair
                    end
                elseif neibour[l] ∈ infsi
                    # st = false
                else
                    st = PR_FBcheck_w(prmodel, new_weight, neibour[l])
                    if st == true
                        sol = value.(all_variables(prmodel)); ndp = getobjval(prmodel)
                        push!(fsi.key, neibour[l]); push!(fsi.x,sol); push!(fsi.y,ndp)
                        push!(candSI,sol); push!(candobj,ndp)
                        if dominated(ndp,PF.Y)==false #&& ndp ∉ Tabu #
                            # push!(candSI,sol)
                            push!(PF.X, sol); push!(PF.Y, ndp);                             # push!(new.X, sol);  push!(new.Y, ndp); push!(Tabu, ndp)
                            newsol1+=1;
                            if sol[1:len[1]] ∉ nd.X
                                # println( iter1,"  new sol ",ndp)
                                push!(X, length(X)+1 => sol[1:len[1]]); push!(Y, length(Y)+1 => ndp); 
                            end
                        end
                    else
                        push!(infsi, neibour[l])
                    end
                end
                # iter2+=1
            end
            if candSI == []
                # println("no candidates")
                break
            else
                SI = nextSI(candSI,candobj)
            end
        end
        @label NextPair
        push!(IGPair,[I,G]);
    end
    @label lsgCal
    df = SortingSol(PF.X,PF.Y)
    ndf,ndset = InitialNDset(df,LPdicho)
    return ndf,ndset,iter1,newsol1,iter2,newsol2,iter3,newsol3
end


# include("/home/ak121396/Dropbox/mdSCNDfunctions.jl")
# file = "/home/ak121396/Desktop/instances/scnd/test01S2"
# fname = file[end-7:end]


include("/home/k2g00/k2g3475/scnd/vopt/mdSCNDfunctions.jl")
#=
file = "/home/k2g00/k2g3475/scnd/instances/test01S2"
dt = Data(file);
fname = file[end-7:end]
tnum = parse(Int,fname[end-3:end-2])
=#
@show file = ARGS[1]
fname = ARGS[1][end-7:end]
dt = Data(file);
fname = ARGS[1][end-7:end]
tnum = parse(Int,fname[end-3:end-2])


if tnum ==  1
    TL = 500;     lTL = 60; l2TL = 10;
elseif tnum == 2
    TL = 600;    lTL = 10; l2TL = 5; 
elseif tnum == 3
    TL = 800;    lTL = 70; l2TL = 10; 
elseif tnum == 4
    TL = 1000;    lTL = 14; l2TL = 10; 
elseif tnum ==5
    TL = 1500;    lTL = 16; l2TL = 5; #20, 20
elseif tnum ==6
    TL = 2000;    lTL = 70; l2TL = 20; #80-120,30
elseif tnum ==7
    TL = 3500;    lTL = 35; l2TL = 30; #37 
elseif tnum ==8
    TL = 4200;    lTL = 200; l2TL = 80; #80-160,90
elseif tnum ==9
    TL = 5000;    lTL = 90; l2TL = 100; #130,80-150
elseif tnum ==10
    TL = 8000;    lTL = 160; l2TL = 30; #200-220,75
elseif tnum ==11
    TL = 10000;    lTL = 320; l2TL = 50; #60
elseif tnum ==12
    TL = 9000;    lTL = 550; l2TL = 300; #800-400
elseif tnum ==13
    TL = 11000;    lTL = 600; l2TL = 200; #800,200-500
elseif tnum ==14
    TL = 14000;    lTL = 500; l2TL = 150; #700-600,400
else
    TL = 20000;    lTL = 300; l2TL = 150; #600
end
scndlp = SCND_LP()
LPtime = @CPUelapsed vSolve(scndlp, 60, method=:dicho, verbose=false)
lp = getvOptData(scndlp);#getY_N(scndlp)
len = [length(scndlp[:y]),length(scndlp[:uij]),length(scndlp[:ujk]),length(scndlp[:ukl]),length(scndlp[:xij]),length(scndlp[:xjk]),length(scndlp[:xkl]),length(scndlp[:h])];

lex1 = lexobj1();
set_optimizer_attribute(lex1, "CPXPARAM_TimeLimit", 3); optimize!(lex1);
set_optimizer_attribute(lex1, "CPXPARAM_TimeLimit",  lTL);#3600
l1time = @CPUelapsed optimize!(lex1); lex1X = [value.(all_variables(lex1))]; lex1Y = [getobjval(lex1)];
println("l1time $fname: ", l1time)

lex2 = lexobj2();
set_optimizer_attribute(lex2, "CPXPARAM_TimeLimit", 3); optimize!(lex2);
set_optimizer_attribute(lex2, "CPXPARAM_TimeLimit", l2TL);#3600
l2time = @CPUelapsed optimize!(lex2); lex2X = [value.(all_variables(lex2))]; lex2Y = [getobjval(lex2)]
println("l2time $fname: ", l2time)

l1mip = SCND_MIP()
JuMP.fix.(l1mip[:y], value.(lex1[:y]); force = true);
l3time = @CPUelapsed vSolve(l1mip, 5, method=:dicho, verbose=false)
stg1 = getvOptData(l1mip);
stg1.Y_N

l2mip = SCND_MIP()
JuMP.fix.(l2mip[:y], value.(lex2[:y]); force = true);
l4time = @CPUelapsed vSolve(l2mip, 5, method=:dicho, verbose=false)
stg2 = getvOptData(l2mip);
stg2.Y_N

rdX = []; temp = [];
for i=1:length(lp.Y_N)
    xt = lp.X_E[i][1:len[1]]    
    yid = findall(p->p>0.2,xt);
    for j=1:len[1]
        if j in yid
            xt[j]=1
        else
            xt[j]=0
        end
    end
    if xt ∉ temp 
        push!(temp, xt); push!(rdX, lp.X_E[i])
    end
end
#FPmodel,PRmodel
fpmodel = FP_Model(); optimize!(fpmodel);
set_optimizer_attribute(fpmodel, "CPXPARAM_TimeLimit", 10)
dist = LP_Model(); optimize!(dist);
set_optimizer_attribute(dist, "CPXPARAM_TimeLimit", 10)
prmodel = PR_Model(); optimize!(prmodel);
set_optimizer_attribute(prmodel, "CPXPARAM_TimeLimit", 10);

#FP start
FP(rdX,len,3,5)
FPtime = @CPUelapsed f1x,f1y = FP(rdX,len,Inf,5);
println("FPtime: ", FPtime, " FPsol: ", length(f1y))
fpX,fpY = NDfilter([f1x;stg1.X_E;stg2.X_E],[f1y;stg1.Y_N;stg2.Y_N])
dfp = DataFrame(X=fpX,Y=fpY); sort!(dfp, [:Y])

#FP+ start
FPplus(dfp.X,dfp.Y,len,3,5)
FPPtime = @CPUelapsed fpp1x,fpp1y = FPplus(dfp.X,dfp.Y,len,round(Int,TL/5),5)
fpp2x,fpp2y = NDfilter(fpp1x,fpp1y);
println("FPPtime: ", FPPtime, " FPPsol: ", length(fpp2y))
dfpp = DataFrame(X=fpp2x,Y=fpp2y);
sort!(dfpp,[:Y])
################################ 
LPdicho = SCND_LP(); InitialNDset(DataFrame(X = [dfpp.X[1]], Y = [dfpp.Y[1]]),LPdicho) #SolveMIPdicho(MIPdicho,dfpp.X[1],dfpp.Y[1],[]);
initime = @CPUelapsed ndf0,ndy0 = InitialNDset(dfpp,LPdicho);
PF = DataFrame(X=ndf0.X, Y=ndf0.Y);
dX = Dict(); dY = Dict(); fsi0 = DataFrame(key=[],x=[],y=[]); 
for i=1:length(ndf0.Y)
    if ndf0.X[i][1:len[1]] ∉ collect(values(dX))
        push!(dX, length(dX)+1 => ndf0.X[i][1:len[1]]); push!(dY, length(dY)+1 => round.(ndf0.Y[i]))
        push!(fsi0.key, round.(ndf0.X[i][1:len[1]]));push!(fsi0.x, ndf0.X[i]); push!(fsi0.y, round.(ndf0.Y[i]))
    end
end
PRdicho(dX,dY,fsi0,PF,len,10,.7)
PRtime = @CPUelapsed  ndf,nodeset,iter1,newsol1,iter2,newsol2,iter3,newsol3 = PRdicho(dX,dY,fsi0,PF,len,TL,1) 
println("PRtime  $fname: ", PRtime)
########################## Saving the output file ###########################
otable = zeros(length(ndf.Y),2)
for i=1:length(ndf.Y)
    for j=1:2
        otable[i,j] = ndf.Y[i][j]
    end
end
CSV.write("/home/k2g00/k2g3475/scnd/vopt/Y/md/1/"*fname*"Y.log", DataFrame(otable, :auto), append=false, header=false,delim=' ')
# nodetb = DataFrame(v1 = [nodeset[i].val[1] for i=1:length(nodeset)], v2 = [nodeset[i].val[2] for i=1:length(nodeset)], arm = [nodeset[i].arm for i=1:length(nodeset)])
# CSV.write("/home/k2g00/k2g3475/scnd/vopt/nodes/1/"*fname*".csv", nodetb, header=true)
# cpuT = LPtime+FPtime+FPPtime+PRtime+l1time+l2time+l3time+l4time+initime
println("mdichotime $fname: ", cpuT," #sol: ", length(ndf.Y))

