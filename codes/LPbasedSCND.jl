using CPUTime,DataFrames,DelimitedFiles,JuMP,LinearAlgebra,CPLEX,vOptGeneric,SparseArrays,StatsBase
struct Data
    filepath::String; N::Dict{}; d::Array{}; cj::Array{}; ck::Array{};  Mij::Array{}; Mjk::Array{}; Mkl::Array{};
    vij::Array{}; vjk::Array{}; vkl::Array{}; rij::Array{}; rjk::Array{}; rkl::Array{}; #gij::Array{}; gjk::Array{}; gkl::Array{};
    Vij::Array{}; Vjk::Array{}; b::Array{}; upl::Int; udc::Int; bigM::Int # e::Array{};q::Array{};
    function Data(filepath)
        dt = readdlm(filepath);
        # notafile = readdlm("/home/ak121396/Desktop/instances/SCND/Notations.txt", '=');
        notafile = readdlm("F:/scnd/Notations.txt", '=');
        # notafile = readdlm("/home/k2g00/k2g3475/scnd/Notations.txt", '=');
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
        d = N["demand"];  cj = N["fcp"]; ck= N["fcd"];
        # e = append!(N["vcp"],N["vcd"]);

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
        b = reshape(N["ves"],N["supplier"],5);
        # q = append!(N["vep"],N["ved"]);
        upl = N["upperpants"]; udc = N["upperdistribution"]
        bigM = sum(sum(N["demand"]))
        new(filepath,N,d,cj,ck,Mij,Mjk,Mkl,vij,vjk,vkl,rij,rjk,rkl,Vij,Vjk,b,upl,udc,bigM); #cap,Mij,Mjk,Mkl,gij,gjk,gkl,
    end
end
# file = "/home/k2g00/k2g3475/scnd/instances/test01S2"
file = "F:scnd/Test1S2"

dt = Data(file);
function voptLPmodel()
    ##########################  Mathematical model  #########################
    scnd = vModel(optimizer_with_attributes(
            CPLEX.Optimizer,
            "CPX_PARAM_EPGAP" => 1e-8
          ));
    set_silent(scnd);
    # MOI.set(scnd, MOI.NumberOfThreads(), 1);
    @variable(scnd, 0<=yj[1:dt.N["plant"],1:2]<=1)
    @variable(scnd, 0<=yk[1:dt.N["distribution"],1:2]<=1)
    @variable(scnd, 0<=uij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:dt.Mij[i,j]]<=1);
    @variable(scnd, 0<=ujk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:dt.Mjk[j,k]]<=1);
    @variable(scnd, 0<=ukl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:dt.Mkl[k,l]]<=1);
    # @variable(scnd, yj[1:dt.N["plant"],1:2], Bin)
    # @variable(scnd, yk[1:dt.N["distribution"],1:2], Bin)
    # @variable(scnd, uij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:dt.Mij[i,j]], Bin);
    # @variable(scnd, ujk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:dt.Mjk[j,k]], Bin);
    # @variable(scnd, ukl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:dt.Mkl[k,l]], Bin);
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

    @addobjective(scnd, Min, sum(dt.cj[j][t]*yj[j,t] for j=1:dt.N["plant"] for t=1:2)+
        sum(dt.ck[k][t]*yk[k,t] for k=1:dt.N["distribution"] for t=1:2) + exg +
        sum(dt.N["vcs"][i][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) +
        sum(dt.vij[i][j][m][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
        sum(dt.vjk[j][k][m][p]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
        sum(dt.vkl[k][l][m][p]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5)+
        sum(dt.N["vcp"][j][2*(p-1)+t]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
        sum(dt.N["vcd"][k][2*(p-1)+t]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2));
    #2nd obj
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
    @constraint(scnd,[j=1:dt.N["plant"], t=1:2], sum(h[j,:,t]) <= dt.N["cap"][j]*yj[j,t]);
    @constraint(scnd,[k=1:dt.N["distribution"],t=1:2], sum(h[k+dt.N["plant"],:,t]) <= dt.N["cad"][k]*yk[k,t]);
    @constraint(scnd,[j=1:dt.N["plant"]], sum(yj[j,:]) <= 1);
    @constraint(scnd,[k=1:dt.N["distribution"]], sum(yk[k,:]) <= 1);

    ########### constraint 10 #############
    @constraint(scnd,[i=1:dt.N["supplier"],j=1:dt.N["plant"]], sum(uij[i,j,m] for m=1:dt.Mij[i,j]) <= 1);
    @constraint(scnd,[j=1:dt.N["plant"],k=1:dt.N["distribution"]], sum(ujk[j,k,m] for m=1:dt.Mjk[j,k]) <= 1);
    @constraint(scnd,[k=1:dt.N["distribution"],l=1:dt.N["customer"]], sum(ukl[k,l,m] for m=1:dt.Mkl[k,l]) <= 1);
    ########### constraint 11 #############
    @constraint(scnd,[i=1:dt.N["supplier"],j=1:dt.N["plant"],m=1:dt.Mij[i,j]], sum(xij[i,j,m,p] for p=1:5) <= dt.bigM*uij[i,j,m]);
    @constraint(scnd,[j=1:dt.N["plant"],k=1:dt.N["distribution"],m=1:dt.Mjk[j,k]] ,sum(xjk[j,k,m,p] for p=1:5) <= dt.bigM*ujk[j,k,m]);
    @constraint(scnd,[k=1:dt.N["distribution"],l=1:dt.N["customer"],m=1:dt.Mkl[k,l]], sum(xkl[k,l,m,p] for p=1:5) <= dt.bigM*ukl[k,l,m]);
    ########### constraint 12 #############
    @constraint(scnd,[i=1:dt.N["supplier"], j=1:dt.N["plant"], m=1:dt.Mij[i,j]], sum(xij[i,j,m,p] for p=1:5) >= dt.Vij[m]*uij[i,j,m] );
    @constraint(scnd,[j=1:dt.N["plant"], k=1:dt.N["distribution"], m=1:dt.Mjk[j,k]], sum(xjk[j,k,m,p] for p=1:5) >= dt.Vjk[m]*ujk[j,k,m]);
    ########### constraint 13-14 #############
    @constraint(scnd,sum(yj) <= dt.upl);
    @constraint(scnd,sum(yk) <= dt.udc);
    return scnd
end
voptlp = voptLPmodel()
# @CPUtime vSolve(scndm, TL, method=:dicho, verbose=false)
@CPUtime vSolve(voptlp, 500, method=:dicho, verbose=false)
LB = getvOptData(voptlp).X_E
len = [length(voptlp[:yj]),length(voptlp[:yk]),length(voptlp[:uij]),length(voptlp[:ujk]),length(voptlp[:ukl]),length(voptlp[:xij]),length(voptlp[:xjk]),length(voptlp[:xkl]),length(voptlp[:h])]

weight = round(Int,mean([getvOptData(voptlp).Y_N[i][1]/getvOptData(voptlp).Y_N[i][2] for i=1:length(getvOptData(voptlp).Y_N)]))

function FPmodel(weight)
    ##########################  Mathematical model  #########################
    scnd = Model(CPLEX.Optimizer); set_silent(scnd);
    # MOI.set(scnd, MOI.NumberOfThreads(), 1);
    @variable(scnd, yj[1:dt.N["plant"],1:2], Bin);
    @variable(scnd, yk[1:dt.N["distribution"],1:2], Bin);
    @variable(scnd, uij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:dt.Mij[i,j]], Bin);
    @variable(scnd, ujk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:dt.Mjk[j,k]], Bin);
    @variable(scnd, ukl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:dt.Mkl[k,l]], Bin);
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

    @objective(scnd, Min, sum(dt.cj[j][t]*yj[j,t] for j=1:dt.N["plant"] for t=1:2)+
        sum(dt.ck[k][t]*yk[k,t] for k=1:dt.N["distribution"] for t=1:2) + exg +
        sum(dt.N["vcs"][i][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) +
        sum(dt.vij[i][j][m][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
        sum(dt.vjk[j][k][m][p]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
        sum(dt.vkl[k][l][m][p]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5)+
        sum(dt.N["vcp"][j][2*(p-1)+t]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
        sum(dt.N["vcd"][k][2*(p-1)+t]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2)+
        weight*(sum(dt.b[i,p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
        sum(dt.N["vep"][j][2*(p-1)+t]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
        sum(dt.N["ved"][k][2*(p-1)+t]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2)+
        sum(dt.rij[i][j][m]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
        sum(dt.rjk[j][k][m]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
        sum(dt.rkl[k][l][m]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5)));
    ######### constraint 3 #############
    @constraint(scnd,[j=1:dt.N["plant"],p=1:5], sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j]) == sum(xjk[j,k,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]));
    @constraint(scnd,[k=1:dt.N["distribution"],p=1:5], sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k]) == sum(xkl[k,l,m,p] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l]));
    ########### constraint 4-6 #############
    @constraint(scnd,[j=1:dt.N["plant"],p=1:5], sum(h[j,p,:]) == sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j]));
    @constraint(scnd,[k=1:dt.N["distribution"],p=1:5], sum(h[k+dt.N["plant"],p,:] ) == sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k]));
    @constraint(scnd, [l=1:dt.N["customer"],p=1:5], sum(xkl[k,l,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mkl[k,l]) >= dt.d[l][p]);
    ########### constraint 7-9 #############
    @constraint(scnd,[i=1:dt.N["supplier"]], sum(xij[i,j,m,p] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) <= dt.N["cas"][i]);
    @constraint(scnd,[j=1:dt.N["plant"], t=1:2], sum(h[j,:,t]) <= dt.N["cap"][j]*yj[j,t]);
    @constraint(scnd,[k=1:dt.N["distribution"],t=1:2], sum(h[k+dt.N["plant"],:,t]) <= dt.N["cad"][k]*yk[k,t]);
    # @constraint(scnd,[j=1:dt.N["plant"]], sum(yj[j,:]) <= 1);
    # @constraint(scnd,[k=1:dt.N["distribution"]], sum(yk[k,:]) <= 1);
    ########### constraint 10 #############
    # @constraint(scnd,[i=1:dt.N["supplier"],j=1:dt.N["plant"]], sum(uij[i,j,m] for m=1:dt.Mij[i,j]) <= 1);
    # @constraint(scnd,[j=1:dt.N["plant"],k=1:dt.N["distribution"]], sum(ujk[j,k,m] for m=1:dt.Mjk[j,k]) <= 1);
    # @constraint(scnd,[k=1:dt.N["distribution"],l=1:dt.N["customer"]], sum(ukl[k,l,m] for m=1:dt.Mkl[k,l]) <= 1);
    ########### constraint 11 #############
    @constraint(scnd,[i=1:dt.N["supplier"],j=1:dt.N["plant"],m=1:dt.Mij[i,j]], sum(xij[i,j,m,p] for p=1:5) <= dt.bigM*uij[i,j,m]);
    @constraint(scnd,[j=1:dt.N["plant"],k=1:dt.N["distribution"],m=1:dt.Mjk[j,k]] ,sum(xjk[j,k,m,p] for p=1:5) <= dt.bigM*ujk[j,k,m]);
    @constraint(scnd,[k=1:dt.N["distribution"],l=1:dt.N["customer"],m=1:dt.Mkl[k,l]], sum(xkl[k,l,m,p] for p=1:5) <= dt.bigM*ukl[k,l,m]);
    ########### constraint 12 #############
    @constraint(scnd,[i=1:dt.N["supplier"], j=1:dt.N["plant"], m=1:dt.Mij[i,j]], sum(xij[i,j,m,p] for p=1:5) >= dt.Vij[m]*uij[i,j,m] );
    @constraint(scnd,[j=1:dt.N["plant"], k=1:dt.N["distribution"], m=1:dt.Mjk[j,k]], sum(xjk[j,k,m,p] for p=1:5) >= dt.Vjk[m]*ujk[j,k,m]);
    ########### constraint 13-14 #############
    # @constraint(scnd,sum(yj) <= dt.upl);
    # @constraint(scnd,sum(yk) <= dt.udc);
    return scnd
end
fplp = FPmodel(weight);

candlist = copy(LB)
for k=1:length(LB)
    k = rand(1:length(candlist))
    x_t = candlist[k]
    yjt = x_t[1:len[1]]; ykt = x_t[1+len[1]:sum(len[i] for i=1:2)]
    jid = findall(p->p>0.2,yjt); kid = findall(p->p>0.2,ykt);
    if Rounding(jid,kid) ==true
        println("YESSSS")
    end
end

1
# using Combinatorics
k = rand(1:length(candlist))
x_t = candlist[k]
yjt = x_t[1:len[1]]; ykt = x_t[1+len[1]:sum(len[i] for i=1:2)]
jid = findall(p->p>0.2,yjt); kid = findall(p->p>0.2,ykt);
u1t = x_t[1+sum(len[i] for i=1:2):sum(len[i] for i=1:3)]
u2t = x_t[1+sum(len[i] for i=1:3):sum(len[i] for i=1:4)]
u3t = x_t[1+sum(len[i] for i=1:4):sum(len[i] for i=1:5)]


function FixYandU(jid,kid)
    jset = collect(combinations(jid))
    for j in jset
        oddj = findall(p->isodd(p)==true && p+1 in j, j)
        deletej = [sample(oddj[i]:oddj[i]+1) for i=1:length(oddj)]
        deleteat!(j,deletej)
    end
    unique(jset)
    kset = collect(combinations(kid))
    for k in kset
        oddk = findall(p->isodd(p)==true && p+1 in k, k)
        deletek = [sample(oddk[i]:oddk[i]+1) for i=1:length(oddk)]
        deleteat!(k,deletek)
    end
    unique(kset)
    filter!(p->length(p)>dt.upl, kset)
    # jset,kset = JKset(yjt,ykt)
    # jplant = sample(jid,rand(1:min(length(jid),dt.upl)), replace=false)
    # kdc = sample(kid,rand(dt.upl:min(length(kid),dt.udc)), replace=false)
    jloca = sample(jset); kloca = sample(kset)

    for j=1:len[1]
        if j in jloca
            yjt[j]=1
        else
            yjt[j]=0
        end
    end
    for k=1:len[2]
        if k in kloca
            ykt[k]=1
        else
            ykt[k]=0
        end
    end
    yj = reshape(yjt,2,Int(length(yjt)/2))'
    yk = reshape(ykt,2,Int(length(ykt)/2))'
    findj = findall(j->j==1,yj)
    fixj = [findj[i][1] for i=1:length(findj)]
    findk = findall(j->j==1,yk)
    fixk = [findk[i][1] for i=1:length(findk)]
    fix.(fplp[:yj], yj; force=true); fix.(fplp[:yk], yk; force=true)

    u1id = []; u3id = []; u2id = []
    fix.(fplp[:uij], 0; force=true); fix.(fplp[:ujk], 0; force=true); fix.(fplp[:ukl], 0; force=true)
    for i=1:dt.N["supplier"]
        for j in fixj
            m = sample(1:dt.Mij[i,j])
            fix(fplp[:uij][i,j,m], 1; force=true)
            push!(u1id,[i,j,m])
        end
    end
    for j in fixj
        for k in fixk
            m = sample(1:dt.Mjk[j,k])
            fix(fplp[:ujk][j,k,m], 1; force=true)
            push!(u2id,[j,k,m])
        end
    end
    for k in fixk
        for l=1:dt.N["customer"]
            m = sample(1:dt.Mkl[k,l])
            fix(fplp[:ukl][k,l,m], 1; force=true)
            push!(u3id,[k,l,m])
        end
    end

    optimize!(fplp)
    if termination_status(fplp) == MOI.OPTIMAL
        return true
    else
        u1r = []
        for i=1:dt.N["supplier"]
            for j=1:dt.N["plant"]
                for m=1:dt.Mij[i,j]
                    if [i,j,m] in u1id
                        push!(u1r,1)
                    else
                        push!(u1r,0)
                    end
                end
            end
        end

        u2r = []
        for j=1:dt.N["plant"]
            for k=1:dt.N["distribution"]
                for m=1:dt.Mjk[j,k]
                    if [j,k,m] in u2id
                        push!(u2r,1)
                    else
                        push!(u2r,0)
                    end
                end
            end
        end

        u3r = []
        for k=1:dt.N["distribution"]
            for l=1:dt.N["customer"]
                for m=1:dt.Mkl[k,l]
                    if [k,l,m] in u3id
                        push!(u3r,1)
                    else
                        push!(u3r,0)
                    end
                end
            end
        end
        return false,yjt,ykt,u1r,u2r,u3r
    end
end

function FixY(jid,kid)
    jset = collect(combinations(jid))
    for j in jset
        oddj = findall(p->isodd(p)==true && p+1 in j, j)
        deletej = [sample(oddj[i]:oddj[i]+1) for i=1:length(oddj)]
        deleteat!(j,deletej)
    end
    unique(jset)
    kset = collect(combinations(kid))
    for k in kset
        oddk = findall(p->isodd(p)==true && p+1 in k, k)
        deletek = [sample(oddk[i]:oddk[i]+1) for i=1:length(oddk)]
        deleteat!(k,deletek)
    end
    unique(kset)
    filter!(p->length(p)>dt.upl, kset)
    # jset,kset = JKset(yjt,ykt)
    # jplant = sample(jid,rand(1:min(length(jid),dt.upl)), replace=false)
    # kdc = sample(kid,rand(dt.upl:min(length(kid),dt.udc)), replace=false)
    jloca = sample(jset); kloca = sample(kset)

    for j=1:len[1]
        if j in jloca
            yjt[j]=1
        else
            yjt[j]=0
        end
    end
    for k=1:len[2]
        if k in kloca
            ykt[k]=1
        else
            ykt[k]=0
        end
    end
    yj = reshape(yjt,2,Int(length(yjt)/2))'
    yk = reshape(ykt,2,Int(length(ykt)/2))'
    findj = findall(j->j==1,yj)
    fixj = [findj[i][1] for i=1:length(findj)]
    findk = findall(j->j==1,yk)
    fixk = [findk[i][1] for i=1:length(findk)]
    fix.(fplp[:yj], yj; force=true); fix.(fplp[:yk], yk; force=true)

    for i=1:dt.N["supplier"]
        for j in setdiff(1:dt.N["plant"],fixj)
            if dt.Mij[i,j] == 1
                fix(fplp[:uij][i,j,1], 0; force=true)
            else
                fix(fplp[:uij][i,j,1], 0; force=true)
                fix(fplp[:uij][i,j,2], 0; force=true)
            end
        end
    end
    for j in setdiff(1:dt.N["plant"],fixj)
        for k in setdiff(1:dt.N["distribution"],fixk)
            if dt.Mjk[j,k] ==1
                fix(fplp[:ujk][j,k,1], 0; force=true)
            else
                fix(fplp[:ujk][j,k,1], 0; force=true)
                fix(fplp[:ujk][j,k,2], 0; force=true)
            end
        end
    end
    for k in setdiff(1:dt.N["distribution"],fixk)
        for l=1:dt.N["customer"]
            if dt.Mkl[k,l] ==1
                fix(fplp[:ukl][k,l,1], 0; force=true)
            else
                fix(fplp[:ukl][k,l,1], 0; force=true)
                fix(fplp[:ukl][k,l,2], 0; force=true)
            end
        end
    end

    optimize!(fplp)
    if termination_status(fplp) == MOI.OPTIMAL
        return true
    else
        u1r = []
        for i=1:dt.N["supplier"]
            for j=1:dt.N["plant"]
                for m=1:dt.Mij[i,j]
                    if [i,j,m] in u1id
                        push!(u1r,1)
                    else
                        push!(u1r,0)
                    end
                end
            end
        end

        u2r = []
        for j=1:dt.N["plant"]
            for k=1:dt.N["distribution"]
                for m=1:dt.Mjk[j,k]
                    if [j,k,m] in u2id
                        push!(u2r,1)
                    else
                        push!(u2r,0)
                    end
                end
            end
        end

        u3r = []
        for k=1:dt.N["distribution"]
            for l=1:dt.N["customer"]
                for m=1:dt.Mkl[k,l]
                    if [k,l,m] in u3id
                        push!(u3r,1)
                    else
                        push!(u3r,0)
                    end
                end
            end
        end
        return false,yjt,ykt,u1r,u2r,u3r
    end
end

# fplp = FPmodel(weight);
fbcheck = FixYandU(jid,kid)
if fbcheck[1] == true
    sol = value.(all_variables(fblp)); ndp = getobjval(sol)
    if dominated(ndp,collect(values(PF)))==false #&& ndp ∉ PF
        push!(X,sol); PF[k] = ndp#push!(PF,ndp)
        yjr = sol[1:len[1]]; ykr = sol[1+len[1]:sum(len[i] for i=1:2)]
        u1r = sol[1+sum(len[i] for i=1:2):sum(len[i] for i=1:3)]
        u1r = sol[1+sum(len[i] for i=1:3):sum(len[i] for i=1:4)]
        u1r = sol[1+sum(len[i] for i=1:4):sum(len[i] for i=1:5)]
        push!(Yj,yjr); push!(Yk,ykr) push!(U1,u1r); push!(U2,u2r); push!(U3,u3r);
        newsol+=1; SearchDone = true
        deleteat!(candlist, k); ## deleteat!(glist, k);
    end
else
     yjr,ykr,u1r,u2r,u3r = fbcheck[2:6]
    if  [yjr;ykr;u1r;u2r;u3r] ∈ Tabu
            yr = flipoper(Y,yt,yr);
            # u1r = flipoper(U1,u1t,u1r); u2r = flipoper(U2,u2t,u2r); u3r = flipoper(U3,u3t,u3r)


1



##################### fix zero values
# fix.(fplp[:uij], 0; force=true)
for i=1:dt.N["supplier"]
    for j in setdiff(1:dt.N["plant"],fixj)
        m = sample(1:dt.Mij[i,j])
        fix(fplp[:uij][i,j,m], 0; force=true)
    end
end
# fix.(fplp[:ujk], 0; force=true)
for j in setdiff(1:dt.N["plant"],fixj)
    for k in setdiff(1:dt.N["distribution"],fixk)
        m = sample(1:dt.Mjk[j,k])
        fix(fplp[:ujk][j,k,m], 0; force=true)
    end
end
# fix.(fplp[:ukl], 0; force=true)
for k in setdiff(1:dt.N["distribution"],fixk)
    for l=1:dt.N["customer"]
        m = sample(1:dt.Mkl[k,l])
        fix(fplp[:ukl][k,l,m], 0; force=true)
    end
end
