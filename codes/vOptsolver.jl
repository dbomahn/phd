using DataFrames,DelimitedFiles,JuMP,LinearAlgebra,CPLEX,MathOptInterface,vOptGeneric
# using CPUTime
mutable struct Data2
    filepath::String; N::Dict{}; d::Array{}; c::Array{};  Mij::Array{}; Mjk::Array{}; Mkl::Array{};
    gij::Array{}; gjk::Array{}; gkl::Array{}; vij::Array{}; vjk::Array{}; vkl::Array{}; rij::Array{}; rjk::Array{}; rkl::Array{};
    Vij::Array{}; Vjk::Array{}; Vkl::Array{}; b::Array{}; upl::Int; udc::Int; bigM::Int # e::Array{};q::Array{};
    function Data2(filepath)
        dt = readdlm(filepath);
        notafile = readdlm("/home/ak121396/Desktop/instances/SCND/Notations.txt", '=');
        # notafile = readdlm("F:/scnd/Notations.txt", '=');
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
        d = N["demand"];  c = append!(N["fcp"],N["fcd"]);
        # e = append!(N["vcp"],N["vcd"]);

        Mij = transpose(reshape(N["ModeIJ"], (N["plant"],N["supplier"])));
        Mjk = transpose(reshape(N["ModeJK"], (N["distribution"],N["plant"])));
        Mkl = transpose(reshape(N["ModeKL"], (N["customer"],N["distribution"])));

        gij = [];
        for i=1:N["supplier"]
            idx = 1; push!(gij,[]);
            for j=1:N["plant"]
                th = []
                for m=1:Mij[i,j]
                    push!(th, N["fixedcostModesp"][i][idx]);
                    idx+=1
                end
                push!(gij[i],th);
            end
        end
        gjk = [];
        for j=1:N["plant"]
            idx = 1; push!(gjk,[]);
            for k=1:N["distribution"]
                th = []
                for m=1:Mjk[j,k]
                    push!(th, N["fixedcostModepd"][j][idx]);
                    idx+=1
                end
                push!(gjk[j],th);
            end
        end
        gkl = [];
        for k=1:N["distribution"]
            idx = 1; push!(gkl,[]);
            for l=1:N["customer"]
                th= []
                for m=1:Mkl[k,l]
                    push!(th, N["fixedcostModedc"][k][idx]);
                    idx+=1
                end
                push!(gkl[k],th);
            end
        end
        vij = [];
        for i=1:N["supplier"]
            tmp = [];
            for j=1:N["plant"]
                tmp2 = []
                if j==1
                    for m=Mij[i,1]:-1:1
                        push!(tmp2, N["tcp"][i][5*(m-1)+1:5*(m-1)+5])
                    end
                else
                    for m=Mij[i,j]:-1:1
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
                    for m=Mjk[j,1]:-1:1
                        push!(tmp2, N["tcd"][j][5*(m-1)+1:5*(m-1)+5])
                    end
                else
                    for m=Mjk[j,k]:-1:1
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
                    for m=Mkl[k,1]:-1:1
                        push!(tmp2, N["tcc"][k][5*(m-1)+1:5*(m-1)+5])
                    end
                else
                    for m=Mkl[k,l]:-1:1
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
                        push!(tmp2, N["cep"][i][5*(m-1)+1:5*(m-1)+5])
                    end
                else
                    for m=1:Mij[i,j]
                        push!(tmp2, N["cep"][i][5*sum(Mij[i,1:j-1])+5*(m-1)+1:5*sum(Mij[i,1:j-1])+5*(m-1)+5]);
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
                        push!(tmp2, N["ced"][j][5*(m-1)+1:5*(m-1)+5])
                    end
                else
                    for m=1:Mjk[j,k]
                        push!(tmp2, N["ced"][j][5*sum(Mjk[j,1:k-1])+5*(m-1)+1:5*sum(Mjk[j,1:k-1])+5*(m-1)+5]);
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
                        push!(tmp2, N["cec"][k][5*(m-1)+1:5*(m-1)+5])
                    end
                else
                    for m=1:Mkl[k,l]
                        push!(tmp2, N["cec"][k][5*sum(Mkl[k,1:l-1])+5*(m-1)+1:5*sum(Mkl[k,1:l-1])+5*(m-1)+5]);
                    end
                end
                push!(tmp,tmp2);
            end
            push!(rkl,tmp);
        end

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
        b = reshape(N["ves"],5,N["supplier"])';
        # q = append!(N["vep"],N["ved"]);
        upl = N["upperpants"]; udc = N["upperdistribution"]
        bigM = sum(sum(N["demand"]))
        new(filepath,N,d,c,Mij,Mjk,Mkl,gij,gjk,gkl,vij,vjk,vkl,rij,rjk,rkl,Vij,Vjk,Vkl,b,upl,udc,bigM); #cap,Mij,Mjk,Mkl,
    end
end
file = "/home/ak121396/Desktop/instances/SCND/test01S2"
# file = "/home/k2g00/k2g3475/scnd/instances/test01S2"
# file = "F:/model/Test1S2"
dt = Data2(file);
function buildmodel()
    model = vModel(CPLEX.Optimizer);
    set_silent(model)
    @variable(model, y[1:dt.N["plant"]+dt.N["distribution"],1:2], Bin);
    @variable(model, uij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:dt.Mij[i,j]], Bin);
    @variable(model, ujk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:dt.Mjk[j,k]], Bin);
    @variable(model, ukl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:dt.Mkl[k,l]], Bin);
    ############
    @variable(model, 0<= xij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:dt.Mij[i,j],1:5] );
    @variable(model, 0<= xjk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:dt.Mjk[j,k],1:5] );
    @variable(model, 0<= xkl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:dt.Mkl[k,l],1:5] );
    @variable(model, 0<= h[1:dt.N["plant"]+dt.N["distribution"],1:5,1:2] );
    @addobjective(model, Min,  sum(dt.c[j][t]*y[j,t] for j=1:1:dt.N["plant"]+dt.N["distribution"] for t=1:2) +
        sum(dt.gij[i][j][m]*uij[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j])+
        sum(dt.gjk[j][k][m]*ujk[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]) +
        sum(dt.gkl[k][l][m]*ukl[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])+
        sum(dt.N["vcs"][i][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) +
        sum(dt.vij[i][j][m][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
        sum(dt.vjk[j][k][m][p]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
        sum(dt.vkl[k][l][m][p]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5)+
        sum(dt.N["vcp"][j][2*(p-1)+t]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
        sum(dt.N["vcd"][k][2*(p-1)+t]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2));
    @addobjective(model, Min, sum(dt.b[i,p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
        sum(dt.N["vep"][j][2*(p-1)+t]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
        sum(dt.N["ved"][k][2*(p-1)+t]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2)+
        sum(dt.rij[i][j][m][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+sum(dt.rjk[j][k][m][p]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
        sum(dt.rkl[k][l][m][p]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5));
    @constraint(model,[j=1:dt.N["plant"],p=1:5], sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j]) == sum(xjk[j,k,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]))
    @constraint(model,[k=1:dt.N["distribution"],p=1:5], sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k]) == sum(xkl[k,l,m,p] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l]))
    ########### constraint 4-6 #############
    @constraint(model,[j=1:dt.N["plant"],p=1:5], sum(h[j,p,:]) == sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j]))
    @constraint(model,[k=1:dt.N["distribution"],p=1:5], sum(h[k+dt.N["plant"],p,:] ) == sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k]))
    @constraint(model, [l=1:dt.N["customer"],p=1:5], sum(xkl[k,l,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mkl[k,l]) >= dt.d[l][p])
    ########### constraint 7-9 #############
    @constraint(model,[i=1:dt.N["supplier"]], sum(xij[i,j,m,p] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) <= dt.N["cas"][i])
    @constraint(model,[j=1:dt.N["plant"]+dt.N["distribution"], t=1:2], sum(h[j,:,t]) <=[dt.N["cap"];dt.N["cad"]][j]*y[j,t])
    @constraint(model,[j=1:dt.N["plant"]+dt.N["distribution"]], sum(y[j,:]) <= 1);
    ########### constraint 10 #############
    @constraint(model,[i=1:dt.N["supplier"],j=1:dt.N["plant"]], sum(uij[i,j,m] for m=1:dt.Mij[i,j]) <= 1);
    @constraint(model,[j=1:dt.N["plant"],k=1:dt.N["distribution"]], sum(ujk[j,k,m] for m=1:dt.Mjk[j,k]) <= 1);
    @constraint(model,[k=1:dt.N["distribution"],l=1:dt.N["customer"]], sum(ukl[k,l,m] for m=1:dt.Mkl[k,l]) <= 1);
    ########### constraint 11 #############
    BigM = sum(sum(dt.N["demand"]));
    @constraint(model,[i=1:dt.N["supplier"],j=1:dt.N["plant"],m=1:dt.Mij[i,j]], sum(xij[i,j,m,p] for p=1:5) <= BigM*uij[i,j,m]);
    @constraint(model,[j=1:dt.N["plant"],k=1:dt.N["distribution"],m=1:dt.Mjk[j,k]] ,sum(xjk[j,k,m,p] for p=1:5) <= BigM*ujk[j,k,m]);
    @constraint(model,[k=1:dt.N["distribution"],l=1:dt.N["customer"],m=1:dt.Mkl[k,l]], sum(xkl[k,l,m,p] for p=1:5) <= BigM*ukl[k,l,m]);
    ########### constraint 12 #############
    @constraint(model,[i=1:dt.N["supplier"], j=1:dt.N["plant"], m=1:dt.Mij[i,j]], sum(xij[i,j,m,p] for p=1:5) >= dt.Vij[i][j][m]*uij[i,j,m] );
    @constraint(model,[j=1:dt.N["plant"], k=1:dt.N["distribution"], m=1:dt.Mjk[j,k]], sum(xjk[j,k,m,p] for p=1:5) >= dt.Vjk[j][k][m]*ujk[j,k,m]);
    @constraint(model,[k=1:dt.N["distribution"], l=1:dt.N["customer"], m=1:dt.Mkl[k,l]], sum(xkl[k,l,m,p] for p=1:5) >= dt.Vkl[k][l][m]*ukl[k,l,m]);
    ########### constraint 13-14 #############
    @constraint(model,sum(y[j,t] for j=1:dt.N["plant"] for t=1:2) <= dt.upl);
    @constraint(model,sum(y[j,t] for j=dt.N["plant"]+1:dt.N["distribution"]+dt.N["plant"] for t=1:2) <= dt.udc);
    return model
end
model = buildmodel()
TimeLimit = 50
vSolve(model, TimeLimit, method=:dicho)
# vSolve(model, method=:epsilon, step=1000000000.0, verbose=true);
# ndpoints = getY_N(model)
vd = getvOptData(model)
vd.Y_N
vd.X_E[1]

function saveSol(fname::String, Y_N)#, elapsedTime::Float64)
    way = fname*"Y_N"
    open(way, "w") do f
        # write(f, "$elapsedTime \n")
        write(f, "$(length(Y_N)) \n")
        for i in 1:length(Y_N)
            write(f, "$(Y_N[i][1]) $(Y_N[i][2]) \n")
        end
    end
end
path = "/home/k2g00/k2g3475/model/vopt/"*file[36:end]*"Y_N"
saveSol(path, ndpoints)#, elapsedTime)
# ---- Displaying the results (X_E and Y_N)
# for n = 1:length(Y_N)
#     X = value.(x, n)
#     print(findall(elt -> elt ≈ 1, X))
#     println("| z = ",Y_N[n])
# end


########################### 1-dimensional model  ###############################
struct Data1
    file::String; N::Dict{}; d::Array{}; c::Array{}; a::Array{}; gij::SparseVector{}; gjk::SparseVector{}; gkl::SparseVector{};
    Vij::SparseVector{}; Vjk::SparseVector{};  b::Array{}; upl::Int; udc::Int; bigM::Int #Vkl::SparseVector{};
    function Data1(file)
        dt1 = readdlm(file);
        # notafile = readdlm("/home/ak121396/Desktop/instances/SCND/Notations.txt", '=');
        # notafile = readdlm("F:/scnd/Notations.txt", '=');
        notafile = readdlm("/home/k2g00/k2g3475/scnd/Notations.txt", '=');
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
        d = reshape(N["demand"],5,N["customer"])'; c = append!(N["fcp"],N["fcd"]); a = reshape(N["vcs"],5,N["supplier"])';
        gij = sparse(N["fixedcostModesp"]); gjk = sparse(N["fixedcostModepd"]); gkl = sparse(N["fixedcostModedc"]);
        Vij = sparse(N["LcapacityModesp"]); Vjk = sparse(N["LcapacityModepd"]);
        # Vkl =  sparse(N["LcapacityModedc"]);
        b = reshape(N["ves"],5,N["supplier"])'; upl = N["upperpants"]; udc = N["upperdistribution"]; bigM = sum(N["demand"])

        new(file,N,d,c,a,gij,gjk,gkl,Vij,Vjk,b,upl,udc,bigM); #,Vkl
    end
end
struct Data2
    filepath::String; N::Dict{}; Mij::Array{}; Mjk::Array{}; Mkl::Array{}; vij::Array{}; vjk::Array{}; vkl::Array{}; rij::Array{}; rjk::Array{}; rkl::Array{};
    function Data2(filepath)
        dt = readdlm(filepath);
        # notafile = readdlm("/home/ak121396/Desktop/instances/SCND/Notations.txt", '=');
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

        Mij = transpose(reshape(N["ModeIJ"], (N["plant"],N["supplier"])));
        Mjk = transpose(reshape(N["ModeJK"], (N["distribution"],N["plant"])));
        Mkl = transpose(reshape(N["ModeKL"], (N["customer"],N["distribution"])));

        vij = [];
        for i=1:N["supplier"]
            tmp = [];
            for j=1:N["plant"]
                tmp2 = []
                if j==1
                    for m=Mij[i,1]:-1:1
                        push!(tmp2, N["tcp"][i][5*(m-1)+1:5*(m-1)+5])
                    end
                else
                    for m=Mij[i,j]:-1:1
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
                    for m=Mjk[j,1]:-1:1
                        push!(tmp2, N["tcd"][j][5*(m-1)+1:5*(m-1)+5])
                    end
                else
                    for m=Mjk[j,k]:-1:1
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
                    for m=Mkl[k,1]:-1:1
                        push!(tmp2, N["tcc"][k][5*(m-1)+1:5*(m-1)+5])
                    end
                else
                    for m=Mkl[k,l]:-1:1
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
                        push!(tmp2, N["cep"][i][5*(m-1)+1:5*(m-1)+5])
                    end
                else
                    for m=1:Mij[i,j]
                        push!(tmp2, N["cep"][i][5*sum(Mij[i,1:j-1])+5*(m-1)+1:5*sum(Mij[i,1:j-1])+5*(m-1)+5]);
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
                        push!(tmp2, N["ced"][j][5*(m-1)+1:5*(m-1)+5])
                    end
                else
                    for m=1:Mjk[j,k]
                        push!(tmp2, N["ced"][j][5*sum(Mjk[j,1:k-1])+5*(m-1)+1:5*sum(Mjk[j,1:k-1])+5*(m-1)+5]);
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
                        push!(tmp2, N["cec"][k][5*(m-1)+1:5*(m-1)+5])
                    end
                else
                    for m=1:Mkl[k,l]
                        push!(tmp2, N["cec"][k][5*sum(Mkl[k,1:l-1])+5*(m-1)+1:5*sum(Mkl[k,1:l-1])+5*(m-1)+5]);
                    end
                end
                push!(tmp,tmp2);
            end
            push!(rkl,tmp);
        end
        new(filepath,N,Mij,Mjk,Mkl,vij,vjk,vkl,rij,rjk,rkl); #cap,Mij,Mjk,Mkl,
    end
end
dt = Data1(file); dt2 = Data2(file);
function buildMIP()
    scnd1 = vModel(CPLEX.Optimizer); set_silent(scnd1)
    @variable(scnd1, y[1:(dt.N["plant"]+dt.N["distribution"])*2], Bin);
    @variable(scnd1, uij[1:sum(dt2.Mij)], Bin);
    @variable(scnd1, ujk[1:sum(dt2.Mjk)], Bin);
    @variable(scnd1, ukl[1:sum(dt2.Mkl)], Bin);

    @variable( scnd1, 0<= xij[1:sum(dt2.Mij)*5] );
    @variable( scnd1, 0<= xjk[1:sum(dt2.Mjk)*5] );
    @variable( scnd1, 0<= xkl[1:sum(dt2.Mkl)*5] );
    @variable( scnd1, 0<= h[1:(dt.N["plant"]+dt.N["distribution"])*5*2] );
    @addobjective(scnd1, Min, (sum(dt.c.*y) + sum(dt.gij[i]*uij[i] for i in findnz(dt.gij)[1]) + sum(dt.gjk[i]*ujk[i] for i in findnz(dt.gjk)[1]) + sum(dt.gkl[i].*ukl[i] for i in findnz(dt.gkl)[1])+
        sum(repeat(dt.a[1,:], outer=sum(dt2.Mij[1,:])).*xij[1:sum(dt2.Mij[1,:])*5])+ sum(sum(repeat(dt.a[i,:], outer=sum(dt2.Mij[i,:])).*xij[sum(dt2.Mij[1:i-1,:])*5+1:sum(dt2.Mij[1:i,:])*5]) for i=2:dt.N["supplier"])+
        sum([dt.N["vcp"];dt.N["vcd"]].*h)))
    @addobjective(scnd1, Min,(sum(dt.N["tcp"].*xij)+sum(dt.N["tcd"].*xjk)+sum(dt.N["tcc"].*xkl)+
        sum(repeat(dt.b[1,:], outer=sum(dt2.Mij[1,:])).*xij[1:sum(dt2.Mij[1,:])*5]) +
        sum(sum(repeat(dt.b[i,:], outer=sum(dt2.Mij[i,:])).*xij[sum(dt2.Mij[1:i-1,:])*5+1:sum(dt2.Mij[1:i,:])*5]) for i=2:dt.N["supplier"]) +
        sum([dt.N["vep"];dt.N["ved"]].*h) + sum(dt.N["cep"].*xij)+sum(dt.N["ced"].*xjk)+sum(dt.N["cec"].*xkl)))

    @constraint(scnd1, [p=1:5], sum(xij[5*(m-1)+p] for m=1:dt2.Mij[1,1])+sum(xij[5*(m-1)+p+(5*sum(dt2.Mij[1:i-1,:]))] for i=2:dt.N["supplier"] for m=1:dt2.Mij[i,1]) == sum(xjk[5*(m-1)+p] for m=1:sum(dt2.Mjk[1,:])) )
    @constraint(scnd1, [j=2:dt.N["plant"],p=1:5], sum(xij[5*sum(dt2.Mij[1,1:j-1])+5*(m-1)+p] for m=1:dt2.Mij[1,j])+sum(xij[5*sum(dt2.Mij[i,1:j-1])+5*(m-1)+p+(5*sum(dt2.Mij[1:i-1,:]))] for i=2:dt.N["supplier"] for m=1:dt2.Mij[i,j]) == sum(xjk[sum(dt2.Mjk[1:j-1,:])*5 + 5*(m-1)+p] for m=1:sum(dt2.Mjk[j,:])) )
    @constraint(scnd1, [p=1:5], sum(xjk[5*(m-1)+p] for m=1:dt2.Mjk[1,1])+sum(xjk[5*(m-1)+p+(5*sum(dt2.Mjk[1:j-1,:]))] for j=2:dt.N["plant"] for m=1:dt2.Mjk[j,1]) == sum(xkl[5*(m-1)+p] for m=1:sum(dt2.Mkl[1,:])) )
    @constraint(scnd1, [k=2:dt.N["distribution"],p=1:5],sum(xjk[5*sum(dt2.Mjk[1,1:k-1])+5*(m-1)+p] for m=1:dt2.Mjk[1,k])+sum(xjk[5*sum(dt2.Mjk[j,1:k-1])+5*(m-1)+p+(5*sum(dt2.Mjk[1:j-1,:]))] for j=2:dt.N["plant"] for m=1:dt2.Mjk[j,k]) == sum(xkl[sum(dt2.Mkl[1:k-1,:])*5 + 5*(m-1)+p] for m=1:sum(dt2.Mkl[k,:])) )

    @constraint(scnd1, [p=1:5],sum(h[2*(p-1)+t] for t=1:2) == sum(xij[5*(m-1)+p] for m=1:dt2.Mij[1,1])+sum(xij[5*(m-1)+p+(5*sum(dt2.Mij[1:i-1,:]))] for i=2:dt.N["supplier"] for m=1:dt2.Mij[i,1]))
    @constraint(scnd1, [j=2:dt.N["plant"],p=1:5], sum(h[5*2*(j-1)+2*(p-1)+t] for t=1:2) == sum(xij[5*sum(dt2.Mij[1,1:j-1])+5*(m-1)+p] for m=1:dt2.Mij[1,j])+sum(xij[5*sum(dt2.Mij[i,1:j-1])+5*(m-1)+p+(5*sum(dt2.Mij[1:i-1,:]))] for i=2:dt.N["supplier"] for m=1:dt2.Mij[i,j]) )
    @constraint(scnd1, [p=1:5], sum(h[5*2*dt.N["plant"]+2*(p-1)+t] for t=1:2) == sum(xjk[5*(m-1)+p] for m=1:dt2.Mjk[1,1])+sum(xjk[5*(m-1)+p+(5*sum(dt2.Mjk[1:j-1,:]))] for j=2:dt.N["plant"] for m=1:dt2.Mjk[j,1]) )
    @constraint(scnd1, [k=2:dt.N["distribution"],p=1:5], sum(h[5*2*dt.N["plant"]+5*2*(k-1)+2*(p-1)+t] for t=1:2) == sum(xjk[5*sum(dt2.Mjk[1,1:k-1])+5*(m-1)+p] for m=1:dt2.Mjk[1,k])+sum(xjk[5*sum(dt2.Mjk[j,1:k-1])+5*(m-1)+p+(5*sum(dt2.Mjk[1:j-1,:]))] for j=2:dt.N["plant"] for m=1:dt2.Mjk[j,k]))
    @constraint(scnd1, [p=1:5], sum(xkl[5*(m-1)+p] for m=1:dt2.Mkl[1,1]) +sum(xkl[5*(m-1)+p+(5*sum(dt2.Mkl[1:k-1,:]))] for k=2:dt.N["distribution"] for m=1:dt2.Mkl[k,1]) >= dt.d[1,p])
    @constraint(scnd1, [l=2:dt.N["customer"], p=1:5], sum(xkl[sum(dt2.Mkl[1,1:l-1])*5 + 5*(m-1)+p] for m=1:dt2.Mkl[1,l])+ sum(xkl[5*sum(dt2.Mkl[1:k-1,:])+5*sum(dt2.Mkl[k,1:l-1])+5*(m-1)+p] for k=2:dt.N["distribution"] for m=1:dt2.Mkl[k,l]) >= dt.d[l,p])

    @constraint(scnd1, sum(xij[1:5*sum(dt2.Mij[1,:])]) <= dt.N["cas"][1])
    @constraint(scnd1, [i=2:dt.N["supplier"]],  sum(xij[5*sum(dt2.Mij[1:i-1,:])+1:5*sum(dt2.Mij[1:i,:])]) <= dt.N["cas"][i])
    @constraint(scnd1,[j=1:dt.N["plant"]+dt.N["distribution"], t=1:2], sum(h[5*2*(j-1)+((p-1)*2)+t] for p=1:5) <= [dt.N["cap"];dt.N["cad"]][j]*y[2*(j-1)+t])
    @constraint(scnd1,[j=1:dt.N["plant"]+dt.N["distribution"]], sum(y[2*(j-1)+1:2*(j-1)+2]) <= 1);
    @constraints(scnd1,begin
            sum(uij[1:dt2.Mij[1,1]]) <= 1
            [j=2:dt.N["plant"]], sum(uij[sum(dt2.Mij[1,1:j-1])+1:sum(dt2.Mij[1,1:j-1])+dt2.Mij[1,j]]) <= 1
            [i=2:dt.N["supplier"],j=2:dt.N["plant"]],  sum(uij[sum(dt2.Mij[1:i-1,:])+sum(dt2.Mij[i,1:j-1])+1:sum(dt2.Mij[1:i-1,:])+sum(dt2.Mij[i,1:j-1])+dt2.Mij[i,j]])<= 1
            sum(ujk[1:dt2.Mjk[1,1]]) <= 1
            [k=2:dt.N["distribution"]], sum(ujk[sum(dt2.Mjk[1,1:k-1])+1:sum(dt2.Mjk[1,1:k-1])+dt2.Mjk[1,k]]) <= 1
            [j=2:dt.N["plant"],k=2:dt.N["distribution"]],  sum(ujk[sum(dt2.Mjk[1:j-1,:])+sum(dt2.Mjk[j,1:j-1])+1:sum(dt2.Mjk[1:j-1,:])+sum(dt2.Mjk[j,1:j-1])+dt2.Mjk[j,k]]) <= 1
            sum(ukl[1:dt2.Mkl[1,1]]) <= 1
            [l=2:dt.N["customer"]], sum(ukl[sum(dt2.Mkl[1,1:l-1])+1:sum(dt2.Mkl[1,1:l-1])+dt2.Mkl[1,l]]) <= 1
            [k=2:dt.N["distribution"],l=2:dt.N["customer"]],  sum(ukl[sum(dt2.Mkl[1:k-1,:])+sum(dt2.Mkl[k,1:l-1])+1:sum(dt2.Mkl[1:k-1,:])+sum(dt2.Mkl[k,1:l-1])+dt2.Mkl[k,l]])<= 1
    end);
    @constraints(scnd1, begin
        [i=1:sum(dt2.Mij)], sum(xij[5*(i-1)+1:5*i]) <= dt.bigM*uij[i]
        [i=1:sum(dt2.Mjk)], sum(xjk[5*(i-1)+1:5*i]) <= dt.bigM*ujk[i]
        [i=1:sum(dt2.Mkl)], sum(xkl[5*(i-1)+1:5*i]) <= dt.bigM*ukl[i]
    end)
    @constraints(scnd1, begin
            [i in findnz(dt.Vij)[1]], sum(xij[5*(i-1)+1:5*i]) >= dt.Vij[i]*uij[i]
            [i in findnz(dt.Vjk)[1]], sum(xjk[5*(i-1)+1:5*i]) >= dt.Vjk[i]*ujk[i]
            # [i in findnz(dt.Vkl)[1]], sum(xkl[5*(i-1)+1:5*i]) >= dt.Vkl[i]*ukl[i]
    end);
    @constraint(scnd1, sum(y[1:dt.N["plant"]*2]) <= dt.upl);
    @constraint(scnd1, sum(y[dt.N["plant"]*2+1:end]) <= dt.udc);

    return scnd1
end
scnd = buildMIP()
# AddCuts(scnd,optcuts,feasicuts) #adding cuts generated from initial Benders with obj1
vSolve(scnd, method=:dicho, verbose=true)

#################################################################################################
using vOptGeneric,DelimitedFiles,JuMP,JLD2,LinearAlgebra,CPLEX,MathProgBase,MathOptInterface
const MPB = MathProgBase;

mutable struct CallModel
    lpfile::String; m::Int; n::Int; C::Array{}; B::Array{}; RHS::Dict{}; signs::Array{}; vub::Array{}
    function CallModel(lpfile::String)
        lpmodel=buildlp([-1,0],[2 1],'<',1.5, CplexSolver(CPX_PARAM_SCRIND=0))
        # lpmodel = CPLEX.CplexMathProgModel();
        MPB.loadproblem!(lpmodel,lpfile)
        Bmtx = MPB.getconstrmatrix(lpmodel);
        B = Bmtx[3:end,:]; C = Bmtx[1:2,:]
        m,n=size(B)
        vub = MPB.getvarUB(lpmodel)
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
        new(lpfile,m,n,C,B,RHS,signs,vub)
    end
end
struct Valu
    x::String
    y::String
    dvar::Array{}
    LB::Array{}
    LBmtx::Array{}
    function Valu(x, y)
        JLD2.@load x dv
        dv0 = Array(dv)
        # dv0 = readdlm(x)
        dv1 = round.(dv0; digits = 4)
        objs = round.(readdlm(y); digits = 4)
        ind = findall(i -> 0 in objs[i, :], 1:size(objs)[1])
        dv2 = dv1[setdiff(1:end, ind), :]
        LBmtx = objs[setdiff(1:end, ind), 2:end]
        dvar = [Vector(dv2[i, :]) for i = 1:size(LBmtx)[1]]
        LB = [Vector(LBmtx[i, :]) for i = 1:size(LBmtx)[1]]
        new(x, y, dvar, LB, LBmtx)
    end
end
file = ARGS[1];

mt = CallModel("/home/ak121396/Desktop/relise/newlp.lp");
# pr = Valu(
#     "/home/ak121396/Desktop/relise/test01S2_X.jld2",
#     "/home/ak121396/Desktop/relise/test01S2_img_p.sol",
# );
#

function compute(file)

    mt = CallModel(file)

    bvar = findall(i->i==1,mt.vub);
    model = vModel( CPLEX.Optimizer );
    MOI.set(mip, MOI.RawParameter("CPX_PARAM_SCRIND"), false);
    @variable(model, x[1:mt.n] >=0)

    @addobjective( model, Min, dot(x, mt.C[1,:]) )
    @addobjective( model, Min, dot(x, mt.C[2,:]) )

    for k=1:mt.m
        if mt.signs[k] == "l"
            @constraint(model, dot(mt.B[k,:],x) >= mt.RHS[k])
        elseif mt.signs[k] == "u"
            @constraint(model, dot(mt.B[k,:],x) <= mt.RHS[k])
        else
            @constraint(model, dot(mt.B[k,:],x) == mt.RHS[k])
        end
    end
    start = time()
    vSolve(model, method=:dicho, verbose=false)
    # vSolve(model, method=:lex, verbose=false)


    Y_N = getY_N(model)
    elapsedTime = time() - start
    return Y_N, elapsedTime
end


# Y_N, elapsedTime = compute(file)
file = "/home/ak121396/Desktop/instances/SCND/test01S2"
# mt = CallModel(file)
# pr = Valu(
#     "/home/k2g00/k2g3475/model/vopt/X/test01S2_X.jld2",
#     "/home/k2g00/k2g3475/model/vopt/Y/test01S2_img_p.sol",
# );
#
# collector = []
# for j=1:length(pr.dvar)
#     @show j
#     mip = vModel(CPLEX.Optimizer);
#     MOI.set(mip, MOI.RawParameter("CPX_PARAM_SCRIND"), false);
#     bvar = findall(i -> i == 1, mt.vub);
#     @variable(mip, x[1:mt.n] >= 0);
#     for i = 1:mt.n
#         if i in bvar
#             set_binary(x[i])
#         else
#             JuMP.fix(x[i], pr.dvar[j][i]; force = true)
#         end
#     end
#     for k = 1:mt.m
#         if mt.signs[k] == "l"
#             @constraint(mip, dot(mt.B[k, :], x) >= mt.RHS[k])
#         elseif mt.signs[k] == "u"
#             @constraint(mip, dot(mt.B[k, :], x) <= mt.RHS[k])
#         else
#             @constraint(mip, dot(mt.B[k, :], x) == mt.RHS[k])
#         end
#     end
#     @addobjective( mip, Min, dot(x, mt.C[1,:]) );
#     @addobjective( mip, Min, dot(x, mt.C[2,:]) );
#
#     cput = @CPUelapsed vSolve(mip, method=:dicho, verbose=false)
#     if termination_status(mip) == MOI.OPTIMAL
#         println("solution found!")
#         Y_N = getY_N(mip)
#         push!(collector,Y_N)
#     end
# end
# println(collector)

############FPBH
# using JuMP, FPBH, GLPKMathProgInterface, CPLEX, MathProgBase,CPUTime ,Modof
# const MPB = MathProgBase;
# bvar = find(i->i==1,mt.vub); rvar = find(i->i!=1,mt.vub);

# fpbh_model = ModoModel()
# @variable(fpbh_model, yu[i in bvar], Bin);
# @variable(fpbh_model, xh[i in rvar] >=0 );
# @variable(fpbh_model, x[1:mt.n] );
# objective!(fpbh_model, 1, :Min, dot(mt.C[1,:],x));
# objective!(fpbh_model, 2, :Min, dot(mt.C[2,:],x));
#
# for i=1:mt.n
#     if i in bvar
#         @constraint(fpbh_model, x[i]==yu[i]  )
#     else
#         @constraint(fpbh_model, x[i]==xh[i] )
#     end
# end
# for k=1:mt.m
#     if mt.signs[k] == "l"
#         @constraint(fpbh_model, dot(mt.B[k,:],x) >= mt.RHS[k])
#     elseif mt.signs[k] == "u"
#         @constraint(fpbh_model, dot(mt.B[k,:],x) <= mt.RHS[k])
#     else
#         @constraint(fpbh_model, dot(mt.B[k,:],x) == mt.RHS[k])
#     end
# end
#
# runtime = @CPUelapsed sol = fpbh(fpbh_model, solution_polishing=false, timelimit=60.0)
#

#################### Mathematical Model of SCND ###############
