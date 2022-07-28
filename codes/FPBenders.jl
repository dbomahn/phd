 struct Data1
    file::String; N::Dict{}; d::Array{}; c::Array{}; a::Array{}; gij::SparseVector{}; gjk::SparseVector{}; gkl::SparseVector{};
    Vij::SparseVector{}; Vjk::SparseVector{}; Vkl::SparseVector{}; b::Array{}; upl::Int; udc::Int; bigM::Int
    function Data1(file)
        dt1 = readdlm(file);
        # notafile = readdlm("/home/ak121396/Desktop/instances/SCND/Notations.txt", '=');
        notafile = readdlm("F:/scnd/Notations.txt", '=');
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
        d = reshape(N["demand"],5,N["customer"])'; c = append!(N["fcp"],N["fcd"]); a = reshape(N["vcs"],5,N["supplier"])';
        gij = sparse(N["fixedcostModesp"]); gjk = sparse(N["fixedcostModepd"]); gkl = sparse(N["fixedcostModedc"]);
        Vij = sparse(N["LcapacityModesp"]); Vjk = sparse(N["LcapacityModepd"]); Vkl =  sparse(N["LcapacityModedc"]);
        b = reshape(N["ves"],5,N["supplier"])'; upl = N["upperpants"]; udc = N["upperdistribution"]; bigM = sum(N["demand"])

        new(file,N,d,c,a,gij,gjk,gkl,Vij,Vjk,Vkl,b,upl,udc,bigM);
    end
end

struct Data2
    filepath::String; N::Dict{}; Mij::Array{}; Mjk::Array{}; Mkl::Array{}; vij::Array{}; vjk::Array{}; vkl::Array{}; rij::Array{}; rjk::Array{}; rkl::Array{};
    function Data2(filepath)
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
file = "F:scnd/Test1S2";
dt = Data1(file); dt2 = Data2(file);
struct MasterP
    y::Vector{VariableRef}
    uij::Vector{VariableRef}
    ujk::Vector{VariableRef}
    ukl::Vector{VariableRef}
    θ::VariableRef
    m::Model
end

function MasterP(w)
    mas = Model(CPLEX.Optimizer); set_silent(mas)
    # MOI.set(mas, MOI.RawOptimizerAttribute("CPXPARAM_MIP_Strategy_Search"), 1) # conventional branch-and-cut. without dynamic search
    MOI.set(mas, MOI.NumberOfThreads(), 1)
    MOI.NodeCount()
    @variable(mas, y[1:(dt.N["plant"]+dt.N["distribution"])*2], Bin);
    @variable(mas, uij[1:sum(dt2.Mij)], Bin);
    @variable(mas, ujk[1:sum(dt2.Mjk)], Bin);
    @variable(mas, ukl[1:sum(dt2.Mkl)], Bin);
    @variable(mas, θ>=-1000);
    @constraint(mas,[j=1:dt.N["plant"]+dt.N["distribution"]], sum(y[2*(j-1)+1:2*(j-1)+2]) <= 1);
    @constraints(mas, begin
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
    @constraint(mas, sum(y[1:dt.N["plant"]*2]) <= dt.upl);
    @constraint(mas, sum(y[dt.N["plant"]*2+1:end]) <= dt.udc);
    @objective(mas, Min, w[1]*(sum(dt.c.*y) + sum(dt.gij[i]*uij[i] for i in findnz(dt.gij)[1]) + sum(dt.gjk[i]*ujk[i] for i in findnz(dt.gjk)[1]) + sum(dt.gkl[i].*ukl[i] for i in findnz(dt.gkl)[1]))+ θ );
    return MasterP(y,uij,ujk,ukl,θ,mas);
end
function dmodel(w)
    mas = Model(CPLEX.Optimizer); set_silent(mas)
    # MOI.set(mas, MOI.RawOptimizerAttribute("CPXPARAM_MIP_Strategy_Search"), 1) # conventional branch-and-cut. without dynamic search
    MOI.set(mas, MOI.NumberOfThreads(), 1)
    MOI.NodeCount()
    @variable(mas, y[1:(dt.N["plant"]+dt.N["distribution"])*2], Bin);
    @variable(mas, uij[1:sum(dt2.Mij)], Bin);
    @variable(mas, ujk[1:sum(dt2.Mjk)], Bin);
    @variable(mas, ukl[1:sum(dt2.Mkl)], Bin);
    @variable(mas, θ>=-1000);
    @constraint(mas,[j=1:dt.N["plant"]+dt.N["distribution"]], sum(y[2*(j-1)+1:2*(j-1)+2]) <= 1);
    @constraints(mas, begin
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
    @constraint(mas, sum(y[1:dt.N["plant"]*2]) <= dt.upl);
    @constraint(mas, sum(y[dt.N["plant"]*2+1:end]) <= dt.udc);
    @objective(mas, Min, w[1]*(sum(dt.c.*y) + sum(dt.gij[i]*uij[i] for i in findnz(dt.gij)[1]) + sum(dt.gjk[i]*ujk[i] for i in findnz(dt.gjk)[1]) + sum(dt.gkl[i].*ukl[i] for i in findnz(dt.gkl)[1]))+ θ );
    return MasterP(y,uij,ujk,ukl,θ,mas);
end
struct DualSP
    # data::SubProblemData
    α1::Vector{VariableRef}
    α2::Vector{VariableRef}
    α3::Vector{VariableRef}
    α4::Vector{VariableRef}
    α5::Vector{VariableRef}
    α6::Vector{VariableRef}
    α7::Vector{VariableRef}
    α8::Vector{VariableRef}
    α9::Vector{VariableRef}
    α10::Vector{VariableRef}
    α11::Vector{VariableRef}
    α12::Vector{VariableRef}
    α13::Vector{VariableRef}
    α14::Vector{VariableRef}
    m::Model
end
function DualSP(w)
    sub = direct_model(CPLEX.Optimizer());
    set_optimizer_attribute(sub, "CPX_PARAM_REDUCE", 0); # presolve must be turned off
    MOI.set(sub, MOI.RawOptimizerAttribute("CPXPARAM_MIP_Strategy_Search"), 1) # conventional branch-and-cut
    set_silent(sub) # display is turned off
    MOI.set(sub, MOI.NumberOfThreads(), 1) # of threads

    @variables(sub, begin
        α1[1:dt.N["plant"]*5]
        α2[1:dt.N["distribution"]*5]
        α3[1:dt.N["plant"]*5]
        α4[1:dt.N["distribution"]*5]
    end)
    @variable(sub, α5[1:dt.N["customer"]*5] >= 0);
    @variable(sub, α6[1:dt.N["supplier"]] >= 0);
    @variable(sub, α7[1:dt.N["plant"]*2] >= 0);
    @variable(sub, α8[1:dt.N["distribution"]*2] >= 0);
    @variable(sub, α9[1:sum(dt2.Mij)] >= 0);
    @variable(sub, α10[1:sum(dt2.Mjk)] >= 0);
    @variable(sub, α11[1:sum(dt2.Mkl)] >= 0);
    @variable(sub, α12[1:sum(dt2.Mij)] >= 0);
    @variable(sub, α13[1:sum(dt2.Mjk)] >= 0);
    @variable(sub, α14[1:sum(dt2.Mkl)] >= 0);

    @constraint(sub, [i=1,j=1,m=1:dt2.Mij[i,j],p=1:5], α1[5*(j-1)+p]-α3[5*(j-1)+p]-α6[i]+α9[m]-α12[m] <= w[1]*(dt.a[i,p]+dt2.vij[i][j][m][p])+w[2]*(dt.b[i,p]+dt2.rij[i][j][m][p]))
    @constraint(sub, [i=2:dt.N["supplier"],j=1,m=1:dt2.Mij[i,j],p=1:5], α1[5*(j-1)+p]-α3[5*(j-1)+p]-α6[i]+α9[m+sum(dt2.Mij[1:i-1,:])]-α12[m+sum(dt2.Mij[1:i-1,:])] <= w[1]*(dt.a[i,p]+dt2.vij[i][j][m][p])+w[2]*(dt.b[i,p]+dt2.rij[i][j][m][p]))
    @constraint(sub, [i=1,j=2:dt.N["plant"],m=1:dt2.Mij[i,j],p=1:5], α1[5*(j-1)+p]-α3[5*(j-1)+p]-α6[i]+α9[sum(dt2.Mij[1,1:j-1])+ m]-α12[sum(dt2.Mij[1,1:j-1])+ m] <= w[1]*(dt.a[i,p]+dt2.vij[i][j][m][p])+w[2]*(dt.b[i,p]+dt2.rij[i][j][m][p]))
    @constraint(sub, [i=2:dt.N["supplier"],j=2:dt.N["plant"], m=1:dt2.Mij[i,j],p=1:5],  α1[5*(j-1)+p]-α3[5*(j-1)+p]-α6[i]+α9[sum(dt2.Mij[1:i-1,:])+sum(dt2.Mij[i,1:j-1])+m]-α12[sum(dt2.Mij[1:i-1,:])+sum(dt2.Mij[i,1:j-1])+m] <= w[1]*(dt.a[i,p]+dt2.vij[i][j][m][p])+w[2]*(dt.b[i,p]+dt2.rij[i][j][m][p]))

    @constraint(sub, [j=1,k=1,m=1:dt2.Mjk[j,k],p=1:5], -α1[5*(j-1)+p]+α2[5*(k-1)+p]-α4[5*(k-1)+p]+α10[m]-α13[m] <= w[1]*(dt2.vjk[j][k][m][p])+w[2]*(dt2.rjk[j][k][m][p]) )
    @constraint(sub, [j=2:dt.N["plant"],k=1,m=1:dt2.Mjk[j,k],p=1:5], -α1[5*(j-1)+p]+α2[5*(k-1)+p]-α4[5*(k-1)+p]+α10[sum(dt2.Mjk[1:j-1,:])]-α13[sum(dt2.Mjk[1:j-1,:])+m] <= w[1]*(dt2.vjk[j][k][m][p])+w[2]*(dt2.rjk[j][k][m][p]));
    @constraint(sub, [j=1,k=2:dt.N["distribution"],m=1:dt2.Mjk[j,k],p=1:5], -α1[5*(j-1)+p]+α2[5*(k-1)+p]-α4[5*(k-1)+p]+α10[sum(dt2.Mjk[1,1:k-1])+ m]-α13[sum(dt2.Mjk[1,1:k-1])+ m] <= w[1]*(dt2.vjk[j][k][m][p])+w[2]*(dt2.rjk[j][k][m][p]) )
    @constraint(sub, [j=2:dt.N["plant"],k=2:dt.N["distribution"],m=1:dt2.Mjk[j,k],p=1:5], -α1[5*(j-1)+p]+α2[5*(k-1)+p]-α4[5*(k-1)+p]+α10[sum(dt2.Mjk[1:j-1,:])+sum(dt2.Mjk[j,1:k-1])+m]-α13[sum(dt2.Mjk[1:j-1,:])+sum(dt2.Mjk[j,1:k-1])+m] <= w[1]*(dt2.vjk[j][k][m][p])+w[2]*(dt2.rjk[j][k][m][p]));

    @constraint(sub, [k=1,l=1,m=1:dt2.Mkl[k,l],p=1:5], -α2[5*(k-1)+p]+α5[5*(l-1)+p]+α11[m]-α14[m] <= w[1]*(dt2.vkl[k][l][m][p])+w[2]*(dt2.rkl[k][l][m][p]))
    @constraint(sub, [k=1,l=2:dt.N["customer"],m=1:dt2.Mkl[k,l],p=1:5], -α2[5*(k-1)+p]+α5[5*(l-1)+p]+α11[sum(dt2.Mkl[1,1:l-1])+ m]-α14[sum(dt2.Mkl[1,1:l-1])+ m] <= w[1]*(dt2.vkl[k][l][m][p])+w[2]*(dt2.rkl[k][l][m][p]))
    @constraint(sub, [k=2:dt.N["distribution"],l=1,m=1:dt2.Mkl[k,l],p=1:5], -α2[5*(k-1)+p]+α5[5*(l-1)+p]+α11[sum(dt2.Mkl[1:k-1,:])+m]-α14[sum(dt2.Mkl[1:k-1,:])+m] <= w[1]*(dt2.vkl[k][l][m][p])+w[2]*(dt2.rkl[k][l][m][p]))
    @constraint(sub, [k=2:dt.N["distribution"],l=2:dt.N["customer"],m=1:dt2.Mkl[k,l],p=1:5], -α2[5*(k-1)+p]+α5[5*(l-1)+p]+α11[sum(dt2.Mkl[1:k-1,:])+sum(dt2.Mkl[k,1:l-1])+m]-α14[sum(dt2.Mkl[1:k-1,:])+sum(dt2.Mkl[k,1:l-1])+m] <= w[1]*(dt2.vkl[k][l][m][p])+w[2]*(dt2.rkl[k][l][m][p]))

    @constraint(sub, [j=1:dt.N["plant"],t=1:2,p=1:5], α3[5*(j-1)+p]-α7[2*(j-1)+t] <= w[1]*(dt2.N["vcp"][j][2*(t-1)+p])+w[2]*(dt2.N["vep"][j][2*(t-1)+p]))
    @constraint(sub, [k=1:dt.N["distribution"],t=1:2,p=1:5], α4[5*(k-1)+p]-α8[2*(k-1)+t] <= w[1]*(dt2.N["vcd"][k][2*(t-1)+p])+w[2]*(dt2.N["ved"][k][2*(t-1)+p]))
    return DualSP(α1,α2,α3,α4,α5,α6,α7,α8,α9,α10,α11,α12,α13,α14,sub)
end
function solve_dsp(dsp::DualSP,yb,ubij,ubjk,ubkl)
    @objective(dsp.m, Max, sum(dsp.α5[5*(l-1)+p]*dt.d[l,p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*dsp.α6[i] for i=1:dt.N["supplier"])-
        sum(dt.N["cap"][j]*yb[2*(j-1)+t]*dsp.α7[2*(j-1)+t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*yb[2*(dt.N["plant"]+k-1)+t]*dsp.α8[2*(k-1)+t] for k=1:dt.N["distribution"] for t=1:2)+
        sum(dt.Vij[i]*ubij[i]*dsp.α9[i] for i in findnz(dt.Vij)[1])+
        sum(dt.Vjk[i]*ubjk[i]*dsp.α10[i] for i in findnz(dt.Vjk)[1])+
        sum(dt.Vkl[i]*ubkl[i]*dsp.α11[i] for i in findnz(dt.Vkl)[1])+
        sum(dt.bigM.*ubij.*dsp.α12) - sum(dt.bigM.*ubjk.*dsp.α13) - sum(dt.bigM.*ubkl.*dsp.α14)
        )
    optimize!(dsp.m)
    st = termination_status(dsp.m)
    if st == MOI.OPTIMAL
      return (res = :OptimalityCut, obj = objective_value(dsp.m), α5 = value.(dsp.α5), α6 = value.(dsp.α6), α7 = value.(dsp.α7), α8= value.(dsp.α8), α9= value.(dsp.α9), α10= value.(dsp.α10),α11= value.(dsp.α11), α12= value.(dsp.α12),α13= value.(dsp.α13),α14= value.(dsp.α14) )
    elseif st == MOI.DUAL_INFEASIBLE
        return ( res = :FeasibilityCut, α5 = value.(dsp.α5), α6 = value.(dsp.α6), α7 = value.(dsp.α7), α8= value.(dsp.α8), α9= value.(dsp.α9), α10= value.(dsp.α10), α11= value.(dsp.α11), α12= value.(dsp.α12),α13= value.(dsp.α13),α14= value.(dsp.α14) )
    else
      error("DualSubProblem error: status $st")
    end
end



mp = MasterP(w)

function FPBenders(mp)
    undo = relax_integrality(mp.m)
    optimize!(mp.m)
    if termination_status(mp.m) == MOI.OPTIMAL
        yt = value.(mp.y); utij = value.(mp.uij); utjk = value.(mp.ujk); utkl = value.(mp.ukl);
        FP(yt,utij,utjk,utkl,TL)
end

dsp = DualSP(w)

###################    Benders with Feasibility Pump    ########################
function flip(x_h,j,e)
    if x_h[e[j]]==1
        x_h[e[j]] = 0
    else
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
            R = sample(1:M,Num, replace=false)
            for i in R
                x_h = flip(x_h,r)
                if x_h ∉ Tabu
                    xi = x_h
                end
            end
            j+=1
        end
    end
    return xi
end
function FBcheck(yr,u1r,u2r,u3r)
    undo()
    for k=1:length(yr)
        JuMP.fix(mp.y[k],yr[k]; force=true)
    end
    for k=1:length(u1r)
        JuMP.fix(mp.uij[k],u1r[k]; force=true)
    end
    for k=1:length(u2r)
        JuMP.fix(mp.ujk[k],u2r[k]; force=true)
    end
    for k=1:length(u3r)
        JuMP.fix(mp.ukl[k],u3r[k]; force=true)
    end
    optimize!(mp.m)
    if termination_status(mp.m) == MOI.OPTIMAL
        return true
    else
        return false
    end
end
function fbsearch(yr,u1r,u2r,u3r,mp) #solveLP

    idy_0 = findall(k->k==0, yr)
    idy_1 = findall(k->k==1, yr)
    idu1_0 = findall(k->k==0, u1r)
    idu1_1 = findall(k->k==1, u1r)
    idu2_0 = findall(k->k==0, u2r)
    idu2_1 = findall(k->k==1, u2r)
    idu3_0 = findall(k->k==0, u3r)
    idu3_1 = findall(k->k==1, u3r)
    relax_integrality(mp.m)
    @objective( mp.m, Min, sum(mp.y[i] for i in idy_0) + sum(1-(mp.y[j]) for j in idy_1) +
        sum(mp.uij[i] for i in idu1_0) + sum(1-(mp.uij[j]) for j in idu1_1)+
        sum(mp.ujk[i] for i in idu2_0) + sum(1-(mp.ujk[j]) for j in idu2_1)+
        sum(mp.ukl[i] for i in idu3_0) + sum(1-(mp.ukl[j]) for j in idu3_1)
    )
    optimize!(mp.m)
    if termination_status(mp.m) == MOI.OPTIMAL
        return JuMP.value.(mp.y),JuMP.value.(uij),JuMP.value.(ujk),JuMP.value.(ukl)
    else
        return 0,0,0,0
    end
end
function FP(yt,u1t,u2t,u3t,TL)
	sol = []; SearchDone = false;	Tabu = []; newsol=0; t0=time(); iter=0; Max_iter = 50 #Y = [];
    # while candlist != [] &&  time()-t0 < TL && k < length(candX)+1

    while iter<Max_iter && SearchDone == false
        yr = round.(Int,yt);u1r = round.(Int,u1t);u2r = round.(Int,u2t);u3r = round.(Int,u3t)
        if ( (FBcheck(yr,u1r,u2r,u3r) == true) && [yr;u1r;u2r;u3r] ∉sol)
			push!(sol,[yr,u1r,u2r,u3r]); newsol+=1;
            @show SearchDone = true            #push!(Y,getobjval(x_r,C))
        else
            if [yr;u1r;u2r;u3r] ∈ Tabu
                yr = flipoper(Tabu,yt,yr); u1r = flipoper(Tabu,u1t,u1r); u2r = flipoper(Tabu,u2t,u2r); u3r = flipoper(Tabu,u3t,u3r)
                if any(i->i==[], [yr,u1r,u2r,u3r])
                    @show SearchDone = true
                else
                    if ( (FBcheck(yr,u1r,u2r,u3r) == true) && [yr;u1r;u2r;u3r] ∉sol)
						push!(sol,[yr,u1r,u2r,u3r]); newsol+=1; SearchDone = true            #push!(Y,getobjval(x_r,C))
                    end
                end
            end
            if time()-t0 >= TL
                break
            end
            if SearchDone == false
                push!(Tabu,[yr;u1r;u2r;u3r])
                yt,u1t,u2t,u3t = fbsearch(yr,u1r,u2r,u3r)
                if any(i->i==0, [yt,u1t,u2t,u3t])  #when there's no new feasible lp sol
                    @show SearchDone = true
                end
            end
        end
		iter+=1
    end
    return sol
end
function benders_decomposition(w,mas::MasterP)

    @objective(mas, Min, w[1]*(sum(dt.c.*y) + sum(dt.gij[i]*uij[i] for i in findnz(dt.gij)[1]) + sum(dt.gjk[i]*ujk[i] for i in findnz(dt.gjk)[1]) + sum(dt.gkl[i].*ukl[i] for i in findnz(dt.gkl)[1]))+ θ );
    optimize!(mas); st = termination_status(mas)
    nopt_cons, nfeasi_cons = (0, 0)
    feasicuts = [];   thetas=[] # optcuts = Vector{Float64}[]

    while (st == MOI.INFEASIBLE) || (st == MOI.OPTIMAL)
        optimize!(mas);
        st = termination_status(mas)
        θ1 = value(θ); yb = value.(y); ubij = value.(uij); ubjk = value.(ujk); ubkl = value.(ukl)
        subp = solve_dsp(dsp,yb,ubij,ubjk,ubkl)
        push!(thetas,θ1)
        if subp.res == :OptimalityCut
            @info "Optimality cut found"
            println(θ1,"  and  ",round(subp.obj; digits=4))
            if round(θ1; digits=4) ≥ round(subp.obj; digits=4)
                break
            else
                nopt_cons+=1
                cut = @constraint( mas, θ ≥ sum(subp.α5[l,p]*dt.d[l][p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
                    sum(dt.N["cap"][j]*y[j,t]*subp.α7[j,t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*y[dt.N["plant"]+k,t]*subp.α8[k,t] for k=1:dt.N["distribution"] for t=1:2)+
                    sum(dt.Vij[i][j][m]*uij[i,j,m]*subp.α9[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j]) +
                    sum(dt.Vjk[j][k][m]*ujk[j,k,m]*subp.α10[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]) +
                    sum(dt.Vkl[k][l][m]*ukl[k,l,m]*subp.α11[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l]) -
                    sum(dt.bigM*uij[i,j,m]*subp.α12[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j]) -
                    sum(dt.bigM*ujk[j,k,m]*subp.α13[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]) -
                    sum(dt.bigM*ukl[k,l,m]*subp.α14[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])
                    );
                # push!(optcuts, cut)
            end
        else
            @info "Feasibility cut found"
            nfeasi_cons += 1
            cut = @constraint( mas, 0 ≥ sum(subp.α5[l,p]*dt.d[l][p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
                sum(dt.N["cap"][j]*y[j,t]*subp.α7[j,t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*y[dt.N["plant"]+k,t]*subp.α8[k,t] for k=1:dt.N["distribution"] for t=1:2)+
                sum(dt.Vij[i][j][m]*uij[i,j,m]*subp.α9[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j]) +
                sum(dt.Vjk[j][k][m]*ujk[j,k,m]*subp.α10[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]) +
                sum(dt.Vkl[k][l][m]*ukl[k,l,m]*subp.α11[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])-
                sum(dt.bigM*uij[i,j,m]*subp.α12[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j]) -
                sum(dt.bigM*ujk[j,k,m]*subp.α13[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]) -
                sum(dt.bigM*ukl[k,l,m]*subp.α14[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])
                );
            push!(feasicuts, cut)
            # push!(feasicuts, (subp.α5,subp.α6,subp.α7,subp.α8,subp.α9,subp.α10,subp.α11));
        end
        # @info "Addin g the cut $(cut)"
    end
    return (mas, y, uij, ujk, ukl, nopt_cons, nfeasi_cons, feasicuts,thetas)
end
function masterFP(mp::build_master)
    optimize!(mp.m)    # end
    yt = value.(mp.y)
    u1t = value.(mp.uij); u2t = value.(mp.ujk); u3t = value.(mp.ukl)
    sol = FP(yt,u1t,u2t,u3t,60)
    return sol
end
optimize!(mp.m)    # end
yt = value.(mp.y)
round.(Int,yt)

sol = masterFP(mp)
