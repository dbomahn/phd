cd("C:/Users/AK121396/Desktop/ProjectBenders")
using JuMP,CPLEX,LinearAlgebra,DelimitedFiles,CPUTime,SparseArrays

struct Data1
    file::String; N::Dict{}; d::Array{}; c::Array{}; a::Array{}; gij::SparseVector{}; gjk::SparseVector{}; gkl::SparseVector{};
    Vij::SparseVector{}; Vjk::SparseVector{};  b::Array{}; upl::Int; udc::Int; bigM::Int #Vkl::SparseVector{};
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
# file = "/home/ak121396/Desktop/instances/SCND/test01S2"
dt = Data1(file); dt2 = Data2(file);

function f1Master()
    mas = Model(CPLEX.Optimizer); set_silent(mas)
    # MOI.set(mas, MOI.RawOptimizerAttribute("CPXPARAM_MIP_Strategy_Search"), 1) # conventional branch-and-cut. without dynamic search
    MOI.set(mas, MOI.NumberOfThreads(), 1)
    # MOI.NodeCount()
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
    @objective(mas, Min, sum(dt.c.*y) + sum(dt.gij[i]*uij[i] for i in findnz(dt.gij)[1]) + sum(dt.gjk[i]*ujk[i] for i in findnz(dt.gjk)[1]) + sum(dt.gkl[i].*ukl[i] for i in findnz(dt.gkl)[1])+ θ );
    return MasterP(y,uij,ujk,ukl,θ,mas);
end

function solve_dsp(dsp::DualSP,yb,ubij,ubjk,ubkl)
    @objective(dsp.m, Max, sum(dsp.α5[5*(l-1)+p]*dt.d[l,p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*dsp.α6[i] for i=1:dt.N["supplier"])-
        sum(dt.N["cap"][j]*yb[2*(j-1)+t]*dsp.α7[2*(j-1)+t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*yb[2*(dt.N["plant"]+k-1)+t]*dsp.α8[2*(k-1)+t] for k=1:dt.N["distribution"] for t=1:2)+
        sum(dt.Vij[i]*ubij[i]*dsp.α9[i] for i in findnz(dt.Vij)[1])+
        sum(dt.Vjk[i]*ubjk[i]*dsp.α10[i] for i in findnz(dt.Vjk)[1])+
        # sum(dt.Vkl[i]*ubkl[i]*dsp.α11[i] for i in findnz(dt.Vkl)[1])+
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
function f1DualSP()
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

    @constraint(sub, [i=1,j=1,m=1:dt2.Mij[i,j],p=1:5], α1[5*(j-1)+p]-α3[5*(j-1)+p]-α6[i]+α9[m]-α12[m] <= dt.a[i,p]+dt2.vij[i][j][m][p])
    @constraint(sub, [i=2:dt.N["supplier"],j=1,m=1:dt2.Mij[i,j],p=1:5], α1[5*(j-1)+p]-α3[5*(j-1)+p]-α6[i]+α9[m+sum(dt2.Mij[1:i-1,:])]-α12[m+sum(dt2.Mij[1:i-1,:])] <= dt.a[i,p]+dt2.vij[i][j][m][p])
    @constraint(sub, [i=1,j=2:dt.N["plant"],m=1:dt2.Mij[i,j],p=1:5], α1[5*(j-1)+p]-α3[5*(j-1)+p]-α6[i]+α9[sum(dt2.Mij[1,1:j-1])+ m]-α12[sum(dt2.Mij[1,1:j-1])+ m] <= dt.a[i,p]+dt2.vij[i][j][m][p])
    @constraint(sub, [i=2:dt.N["supplier"],j=2:dt.N["plant"], m=1:dt2.Mij[i,j],p=1:5],  α1[5*(j-1)+p]-α3[5*(j-1)+p]-α6[i]+α9[sum(dt2.Mij[1:i-1,:])+sum(dt2.Mij[i,1:j-1])+m]-α12[sum(dt2.Mij[1:i-1,:])+sum(dt2.Mij[i,1:j-1])+m] <= dt.a[i,p]+dt2.vij[i][j][m][p])

    @constraint(sub, [j=1,k=1,m=1:dt2.Mjk[j,k],p=1:5], -α1[5*(j-1)+p]+α2[5*(k-1)+p]-α4[5*(k-1)+p]+α10[m]-α13[m] <= dt2.vjk[j][k][m][p])
    @constraint(sub, [j=2:dt.N["plant"],k=1,m=1:dt2.Mjk[j,k],p=1:5], -α1[5*(j-1)+p]+α2[5*(k-1)+p]-α4[5*(k-1)+p]+α10[sum(dt2.Mjk[1:j-1,:])]-α13[sum(dt2.Mjk[1:j-1,:])+m] <= dt2.vjk[j][k][m][p])
    @constraint(sub, [j=1,k=2:dt.N["distribution"],m=1:dt2.Mjk[j,k],p=1:5], -α1[5*(j-1)+p]+α2[5*(k-1)+p]-α4[5*(k-1)+p]+α10[sum(dt2.Mjk[1,1:k-1])+ m]-α13[sum(dt2.Mjk[1,1:k-1])+ m] <= dt2.vjk[j][k][m][p])
    @constraint(sub, [j=2:dt.N["plant"],k=2:dt.N["distribution"],m=1:dt2.Mjk[j,k],p=1:5], -α1[5*(j-1)+p]+α2[5*(k-1)+p]-α4[5*(k-1)+p]+α10[sum(dt2.Mjk[1:j-1,:])+sum(dt2.Mjk[j,1:k-1])+m]-α13[sum(dt2.Mjk[1:j-1,:])+sum(dt2.Mjk[j,1:k-1])+m] <= dt2.vjk[j][k][m][p])

    @constraint(sub, [k=1,l=1,m=1:dt2.Mkl[k,l],p=1:5], -α2[5*(k-1)+p]+α5[5*(l-1)+p]+α11[m]-α14[m] <= dt2.vkl[k][l][m][p])
    @constraint(sub, [k=1,l=2:dt.N["customer"],m=1:dt2.Mkl[k,l],p=1:5], -α2[5*(k-1)+p]+α5[5*(l-1)+p]+α11[sum(dt2.Mkl[1,1:l-1])+ m]-α14[sum(dt2.Mkl[1,1:l-1])+ m] <= dt2.vkl[k][l][m][p])
    @constraint(sub, [k=2:dt.N["distribution"],l=1,m=1:dt2.Mkl[k,l],p=1:5], -α2[5*(k-1)+p]+α5[5*(l-1)+p]+α11[sum(dt2.Mkl[1:k-1,:])+m]-α14[sum(dt2.Mkl[1:k-1,:])+m] <= dt2.vkl[k][l][m][p])
    @constraint(sub, [k=2:dt.N["distribution"],l=2:dt.N["customer"],m=1:dt2.Mkl[k,l],p=1:5], -α2[5*(k-1)+p]+α5[5*(l-1)+p]+α11[sum(dt2.Mkl[1:k-1,:])+sum(dt2.Mkl[k,1:l-1])+m]-α14[sum(dt2.Mkl[1:k-1,:])+sum(dt2.Mkl[k,1:l-1])+m] <= dt2.vkl[k][l][m][p])

    @constraint(sub, [j=1:dt.N["plant"],t=1:2,p=1:5], α3[5*(j-1)+p]-α7[2*(j-1)+t] <= dt2.N["vcp"][j][2*(t-1)+p])
    @constraint(sub, [k=1:dt.N["distribution"],t=1:2,p=1:5], α4[5*(k-1)+p]-α8[2*(k-1)+t] <= dt2.N["vcd"][k][2*(t-1)+p])
    return DualSP(α1,α2,α3,α4,α5,α6,α7,α8,α9,α10,α11,α12,α13,α14,sub)
end

function lazy_callback(cb_data)
    yb = callback_value.(cb_data, f1mp.y);    ubij = callback_value.(cb_data, f1mp.uij);    ubjk = callback_value.(cb_data, f1mp.ujk)
    ubkl = callback_value.(cb_data, f1mp.ukl);    θb = callback_value(cb_data, f1mp.θ)
    subp = solve_dsp(f1dsp,yb,ubij,ubjk,ubkl)

    if subp.res == :OptimalityCut
        if round(θb; digits=4) ≥ round(subp.obj; digits=4)
            return
        else
            cut = @build_constraint( f1mp.θ ≥ sum(subp.α5[5*(l-1)+p]*dt.d[l,p] for l=1:dt.N["customer"] for p=1:5)-
                sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
                sum(dt.N["cap"][j]*f1mp.y[2*(j-1)+t]*subp.α7[2*(j-1)+t] for j=1:dt.N["plant"] for t=1:2)-
                sum(dt.N["cad"][k]*f1mp.y[2*(dt.N["plant"]+k-1)+t]*subp.α8[2*(k-1)+t] for k=1:dt.N["distribution"] for t=1:2)+
                sum(dt.Vij[i]*f1mp.uij[i]*subp.α9[i] for i in findnz(dt.Vij)[1])+
                sum(dt.Vjk[i]*f1mp.ujk[i]*subp.α10[i] for i in findnz(dt.Vjk)[1])-
                # sum(dt.Vkl[i]*f1mp.ukl[i]*subp.α11[i] for i in findnz(dt.Vkl)[1])+
                sum(dt.bigM.*f1mp.uij.*subp.α12) - sum(dt.bigM.*f1mp.ujk.*subp.α13) - sum(dt.bigM.*f1mp.ukl.*subp.α14)
                )
            MOI.submit(f1mp.m, MOI.LazyConstraint(cb_data), cut)
            push!(optcuts, (α5=subp.α5,α6=subp.α6,α7=subp.α7,α8=subp.α8,α9=subp.α9,α10=subp.α10,α11=subp.α11,α12=subp.α12,α13=subp.α13,α14=subp.α14))
        end
    else
        cut = @build_constraint( 0 ≥ sum(subp.α5[5*(l-1)+p]*dt.d[l,p] for l=1:dt.N["customer"] for p=1:5)-
            sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
            sum(dt.N["cap"][j]*f1mp.y[2*(j-1)+t]*subp.α7[2*(j-1)+t] for j=1:dt.N["plant"] for t=1:2)-
            sum(dt.N["cad"][k]*f1mp.y[2*(dt.N["plant"]+k-1)+t]*subp.α8[2*(k-1)+t] for k=1:dt.N["distribution"] for t=1:2)+
            sum(dt.Vij[i]*f1mp.uij[i]*subp.α9[i] for i in findnz(dt.Vij)[1])+
            sum(dt.Vjk[i]*f1mp.ujk[i]*subp.α10[i] for i in findnz(dt.Vjk)[1])-
            # sum(dt.Vkl[i]*f1mp.ukl[i]*subp.α11[i] for i in findnz(dt.Vkl)[1])+
            sum(dt.bigM.*f1mp.uij.*subp.α12) - sum(dt.bigM.*f1mp.ujk.*subp.α13) - sum(dt.bigM.*f1mp.ukl.*subp.α14)
            )
        MOI.submit(f1mp.m, MOI.LazyConstraint(cb_data), cut)
        push!(feasicuts, (α5=subp.α5,α6=subp.α6,α7=subp.α7,α8=subp.α8,α9=subp.α9,α10=subp.α10,α11=subp.α11,α12=subp.α12,α13=subp.α13,α14=subp.α14))
    end
    return
end
f1mp = f1Master(); f1dsp = f1DualSP(); optcuts = []; feasicuts = [];
MOI.set(f1mp.m, MOI.LazyConstraintCallback(), lazy_callback);
# MOI.set(f1mp.m, MOI.UserCutCallback(), rootfrac_callback)
set_time_limit_sec(f1mp.m, 60)
optimize!(f1mp.m)
# solve_time(f1mp.m)
length(optcuts)
length(feasicuts),length(unique(feasicuts))
1

function MOLP(w)
    molp = Model(CPLEX.Optimizer); set_silent(molp)
    @variable(molp, 0<= y[1:(dt.N["plant"]+dt.N["distribution"])*2]<=1);
    @variable(molp, 0<=uij[1:sum(dt2.Mij)]<=1);
    @variable(molp, 0<=ujk[1:sum(dt2.Mjk)]<=1);
    @variable(molp, 0<=ukl[1:sum(dt2.Mkl)]<=1);
    @variable(molp, θ>= -1000);
    @variable( molp, 0<= xij[1:sum(dt2.Mij)*5] )
    @variable( molp, 0<= xjk[1:sum(dt2.Mjk)*5] )
    @variable( molp, 0<= xkl[1:sum(dt2.Mkl)*5] )
    @variable( molp, 0<= h[1:(dt.N["plant"]+dt.N["distribution"])*5*2] )

    @constraint(molp, [p=1:5], sum(xij[5*(m-1)+p] for m=1:dt2.Mij[1,1])+sum(xij[5*(m-1)+p+(5*sum(dt2.Mij[1:i-1,:]))] for i=2:dt.N["supplier"] for m=1:dt2.Mij[i,1]) == sum(xjk[5*(m-1)+p] for m=1:sum(dt2.Mjk[1,:])) )
    @constraint(molp, [j=2:dt.N["plant"],p=1:5], sum(xij[5*sum(dt2.Mij[1,1:j-1])+5*(m-1)+p] for m=1:dt2.Mij[1,j])+sum(xij[5*sum(dt2.Mij[i,1:j-1])+5*(m-1)+p+(5*sum(dt2.Mij[1:i-1,:]))] for i=2:dt.N["supplier"] for m=1:dt2.Mij[i,j]) == sum(xjk[sum(dt2.Mjk[1:j-1,:])*5 + 5*(m-1)+p] for m=1:sum(dt2.Mjk[j,:])) )
    @constraint(molp, [p=1:5], sum(xjk[5*(m-1)+p] for m=1:dt2.Mjk[1,1])+sum(xjk[5*(m-1)+p+(5*sum(dt2.Mjk[1:j-1,:]))] for j=2:dt.N["plant"] for m=1:dt2.Mjk[j,1]) == sum(xkl[5*(m-1)+p] for m=1:sum(dt2.Mkl[1,:])) )
    @constraint(molp, [k=2:dt.N["distribution"],p=1:5],sum(xjk[5*sum(dt2.Mjk[1,1:k-1])+5*(m-1)+p] for m=1:dt2.Mjk[1,k])+sum(xjk[5*sum(dt2.Mjk[j,1:k-1])+5*(m-1)+p+(5*sum(dt2.Mjk[1:j-1,:]))] for j=2:dt.N["plant"] for m=1:dt2.Mjk[j,k]) == sum(xkl[sum(dt2.Mkl[1:k-1,:])*5 + 5*(m-1)+p] for m=1:sum(dt2.Mkl[k,:])) )

    @constraint(molp, [p=1:5],sum(h[2*(p-1)+t] for t=1:2) == sum(xij[5*(m-1)+p] for m=1:dt2.Mij[1,1])+sum(xij[5*(m-1)+p+(5*sum(dt2.Mij[1:i-1,:]))] for i=2:dt.N["supplier"] for m=1:dt2.Mij[i,1]))
    @constraint(molp, [j=2:dt.N["plant"],p=1:5], sum(h[5*2*(j-1)+2*(p-1)+t] for t=1:2) == sum(xij[5*sum(dt2.Mij[1,1:j-1])+5*(m-1)+p] for m=1:dt2.Mij[1,j])+sum(xij[5*sum(dt2.Mij[i,1:j-1])+5*(m-1)+p+(5*sum(dt2.Mij[1:i-1,:]))] for i=2:dt.N["supplier"] for m=1:dt2.Mij[i,j]) )
    @constraint(molp, [p=1:5], sum(h[5*2*dt.N["plant"]+2*(p-1)+t] for t=1:2) == sum(xjk[5*(m-1)+p] for m=1:dt2.Mjk[1,1])+sum(xjk[5*(m-1)+p+(5*sum(dt2.Mjk[1:j-1,:]))] for j=2:dt.N["plant"] for m=1:dt2.Mjk[j,1]) )
    @constraint(molp, [k=2:dt.N["distribution"],p=1:5], sum(h[5*2*dt.N["plant"]+5*2*(k-1)+2*(p-1)+t] for t=1:2) == sum(xjk[5*sum(dt2.Mjk[1,1:k-1])+5*(m-1)+p] for m=1:dt2.Mjk[1,k])+sum(xjk[5*sum(dt2.Mjk[j,1:k-1])+5*(m-1)+p+(5*sum(dt2.Mjk[1:j-1,:]))] for j=2:dt.N["plant"] for m=1:dt2.Mjk[j,k]))
    @constraint(molp, [p=1:5], sum(xkl[5*(m-1)+p] for m=1:dt2.Mkl[1,1]) +sum(xkl[5*(m-1)+p+(5*sum(dt2.Mkl[1:k-1,:]))] for k=2:dt.N["distribution"] for m=1:dt2.Mkl[k,1]) >= dt.d[1,p])
    @constraint(molp, [l=2:dt.N["customer"], p=1:5], sum(xkl[sum(dt2.Mkl[1,1:l-1])*5 + 5*(m-1)+p] for m=1:dt2.Mkl[1,l])+ sum(xkl[5*sum(dt2.Mkl[1:k-1,:])+5*sum(dt2.Mkl[k,1:l-1])+5*(m-1)+p] for k=2:dt.N["distribution"] for m=1:dt2.Mkl[k,l]) >= dt.d[l,p])

    @constraint(molp, sum(xij[1:5*sum(dt2.Mij[1,:])]) <= dt.N["cas"][1])
    @constraint(molp, [i=2:dt.N["supplier"]],  sum(xij[5*sum(dt2.Mij[1:i-1,:])+1:5*sum(dt2.Mij[1:i,:])]) <= dt.N["cas"][i])
    @constraint(molp,[j=1:dt.N["plant"]+dt.N["distribution"], t=1:2], sum(h[5*2*(j-1)+((p-1)*2)+t] for p=1:5) <= [dt.N["cap"];dt.N["cad"]][j]*y[2*(j-1)+t])
    @constraint(molp,[j=1:dt.N["plant"]+dt.N["distribution"]], sum(y[2*(j-1)+1:2*(j-1)+2]) <= 1);
    @constraints(molp, begin
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
    @constraints(molp, begin
        [i=1:sum(dt2.Mij)], sum(xij[5*(i-1)+1:5*i]) <= dt.bigM*uij[i]
        [i=1:sum(dt2.Mjk)], sum(xjk[5*(i-1)+1:5*i]) <= dt.bigM*ujk[i]
        [i=1:sum(dt2.Mkl)], sum(xkl[5*(i-1)+1:5*i]) <= dt.bigM*ukl[i]
    end)
    @constraints(molp, begin
            [i in findnz(dt.Vij)[1]], sum(xij[5*(i-1)+1:5*i]) >= dt.Vij[i]*uij[i]
            [i in findnz(dt.Vjk)[1]], sum(xjk[5*(i-1)+1:5*i]) >= dt.Vjk[i]*ujk[i]
            # [i in findnz(dt.Vkl)[1]], sum(xkl[5*(i-1)+1:5*i]) >= dt.Vkl[i]*ukl[i]
    end);
    @constraint(molp, sum(y[1:dt.N["plant"]*2]) <= dt.upl);
    @constraint(molp, sum(y[dt.N["plant"]*2+1:end]) <= dt.udc);
    @objective(molp, Min, w[1]*(sum(dt.c.*y) + sum(dt.gij[i]*uij[i] for i in findnz(dt.gij)[1]) + sum(dt.gjk[i]*ujk[i] for i in findnz(dt.gjk)[1]) + sum(dt.gkl[i].*ukl[i] for i in findnz(dt.gkl)[1])+
            sum(repeat(dt.a[1,:], outer=sum(dt2.Mij[1,:])).*xij[1:sum(dt2.Mij[1,:])*5])+ sum(sum(repeat(dt.a[i,:], outer=sum(dt2.Mij[i,:])).*xij[sum(dt2.Mij[1:i-1,:])*5+1:sum(dt2.Mij[1:i,:])*5]) for i=2:dt.N["supplier"])+
            sum([dt.N["vcp"];dt.N["vcd"]].*h)) +
            w[2]*(sum(dt.N["tcp"].*xij)+sum(dt.N["tcd"].*xjk)+sum(dt.N["tcc"].*xkl)+
            sum(repeat(dt.b[1,:], outer=sum(dt2.Mij[1,:])).*xij[1:sum(dt2.Mij[1,:])*5]) +
            sum(sum(repeat(dt.b[i,:], outer=sum(dt2.Mij[i,:])).*xij[sum(dt2.Mij[1:i-1,:])*5+1:sum(dt2.Mij[1:i,:])*5]) for i=2:dt.N["supplier"]) +
            sum([dt.N["vep"];dt.N["ved"]].*h) + sum(dt.N["cep"].*xij)+sum(dt.N["ced"].*xjk)+sum(dt.N["cec"].*xkl)) + θ  )

    return molp
end
# w = [0.5,0.5]
molp = MOLP(w);
AddCuts(molp,optcuts,unique(feasicuts))
optimize!(molp)
termination_status(molp)
objective_value(molp)
yt = value.(molp[:y]); utij = value.(molp[:uij]); utjk = value.(molp[:ujk]); utkl = value.(molp[:ukl]);
sol = FP(yt,utij,utjk,utkl)



function FBcheck(yr,u1r,u2r,u3r)
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
function fbsearch(yr,u1r,u2r,u3r) #solveLP
    idy_0 = findall(k->k==0, yr)
    idy_1 = findall(k->k==1, yr)
    idu1_0 = findall(k->k==0, u1r)
    idu1_1 = findall(k->k==1, u1r)
    idu2_0 = findall(k->k==0, u2r)
    idu2_1 = findall(k->k==1, u2r)
    idu3_0 = findall(k->k==0, u3r)
    idu3_1 = findall(k->k==1, u3r)
    @objective( di.m, Min, sum(di.y[i] for i in idy_0) + sum(1-(di.y[j]) for j in idy_1) +
        sum(di.uij[i] for i in idu1_0) + sum(1-(di.uij[j]) for j in idu1_1)+
        sum(di.ujk[i] for i in idu2_0) + sum(1-(di.ujk[j]) for j in idu2_1)+
        sum(di.ukl[i] for i in idu3_0) + sum(1-(di.ukl[j]) for j in idu3_1))
    optimize!(di.m)
    if termination_status(di.m) == MOI.OPTIMAL
        return JuMP.value.(di.y),JuMP.value.(di.uij),JuMP.value.(di.ujk),JuMP.value.(di.ukl)
    else
        return 0,0,0,0
    end
end
function FP(yt,u1t,u2t,u3t)
    # yt = value.(lp.y); utij = value.(lp.uij); utjk = value.(lp.ujk); utkl = value.(lp.ukl);
	sol = []; SearchDone = false;	Tabu = []; newsol=0; t0=time(); iter=0; Max_iter = 50 #Y = [];
    # while candlist != [] &&  time()-t0 < TL && k < length(candX)+1

    while SearchDone == false && iter<Max_iter
        yr = round.(Int,yt);u1r = round.(Int,u1t);u2r = round.(Int,u2t);u3r = round.(Int,u3t)
        if ( (FBcheck(yr,u1r,u2r,u3r) == true) && [yr,u1r,u2r,u3r] ∉sol)
			push!(sol,[yr,u1r,u2r,u3r]); newsol+=1;
            SearchDone = true            #push!(Y,getobjval(x_r,C))
        else
            if [yr;u1r;u2r;u3r] ∈ Tabu
                yr = flipoper(Tabu,yt,yr); u1r = flipoper(Tabu,u1t,u1r); u2r = flipoper(Tabu,u2t,u2r); u3r = flipoper(Tabu,u3t,u3r)
                if any(i->i==[], [yr,u1r,u2r,u3r])
                    SearchDone = true
                else
                    if ( (FBcheck(yr,u1r,u2r,u3r) == true) && [yr,u1r,u2r,u3r] ∉sol)
						push!(sol,[yr,u1r,u2r,u3r]); newsol+=1; SearchDone = true            #push!(Y,getobjval(x_r,C))
                    end
                end
            end
            # if time()-t0 >= TL
            #     break
            # end
            if SearchDone == false
                push!(Tabu,[yr;u1r;u2r;u3r])
                yt,u1t,u2t,u3t = fbsearch(yr,u1r,u2r,u3r)
                if any(i->i==0, [yt,u1t,u2t,u3t])  #when there's no new feasible lp sol
                    SearchDone = true
                end
            end
        end
		@show iter+=1
    end
    return sol
end

function benders_decomposition(yt,utij,utjk,utkl)
    ocuts = []; fcuts = []; Archiv = [];
    optimize!(mp.m); st = termination_status(mp.m)

    while (st == MOI.INFEASIBLE) || (st == MOI.OPTIMAL)
        optimize!(molp)
        yt = value.(molp[:y]); utij = value.(molp[:uij]); utjk = value.(molp[:ujk]); utkl = value.(molp[:ukl]);
        sol = FP(yt,utij,utjk,utkl);
        if sol != [] &&  sol ∉ Archiv # feasible sol found by FP
            push!(Archiv, sol)
            yb = sol[1][1]; ubij  = sol[1][2]; ubjk = sol[1][3]; ubkl = sol[1][4]; θb = value(mp.θ);
            subp = solve_dsp(dsp,yb,ubij,ubjk,ubkl)
            if subp.res == :OptimalityCut
                @info "FP == Optimality cut found"
                if round(θb; digits=4) ≥ round(subp.obj; digits=4)
                    break
                else
                    # nopt_cons+=1
                    cut = @constraint( mp.m, mp.θ ≥ sum(subp.α5[5*(l-1)+p]*dt.d[l,p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
                        sum(dt.N["cap"][j]*mp.y[2*(j-1)+t]*subp.α7[2*(j-1)+t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*mp.y[2*(dt.N["plant"]+k-1)+t]*subp.α8[2*(k-1)+t] for k=1:dt.N["distribution"] for t=1:2)+
                        sum(dt.Vij[i]*mp.uij[i]*subp.α9[i] for i in findnz(dt.Vij)[1])+
                        sum(dt.Vjk[i]*mp.ujk[i]*subp.α10[i] for i in findnz(dt.Vjk)[1])-
                        # sum(dt.Vkl[i]*ukl[i]*subp.α11[i] for i in findnz(dt.Vkl)[1])+
                        sum(dt.bigM.*mp.uij.*subp.α12) - sum(dt.bigM.*mp.ujk.*subp.α13) - sum(dt.bigM.*mp.ukl.*subp.α14))

                    @constraint( molp, molp[:θ] ≥ sum(subp.α5[5*(l-1)+p]*dt.d[l,p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
                        sum(dt.N["cap"][j]*molp[:y][2*(j-1)+t]*subp.α7[2*(j-1)+t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*mp.y[2*(dt.N["plant"]+k-1)+t]*subp.α8[2*(k-1)+t] for k=1:dt.N["distribution"] for t=1:2)+
                        sum(dt.Vij[i]*molp[:uij][i]*subp.α9[i] for i in findnz(dt.Vij)[1])+
                        sum(dt.Vjk[i]*molp[:ujk][i]*subp.α10[i] for i in findnz(dt.Vjk)[1])-
                        # sum(dt.Vkl[i]*ukl[i]*subp.α11[i] for i in findnz(dt.Vkl)[1])+
                        sum(dt.bigM.*molp[:uij].*subp.α12) - sum(dt.bigM.*molp[:ujk].*subp.α13) - sum(dt.bigM.*molp[:ukl].*subp.α14))
                    push!(ocuts, (α5=subp.α5,α6=subp.α6,α7=subp.α7,α8=subp.α8,α9=subp.α9,α10=subp.α10,α11=subp.α11,α12=subp.α12,α13=subp.α13,α14=subp.α14))
                end
            else
                @info "FP == Feasibility cut found"
                # nfeasi_cons += 1
                cut = @constraint( mp.m, 0 ≥ sum(subp.α5[5*(l-1)+p]*dt.d[l,p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
                    sum(dt.N["cap"][j]*mp.y[2*(j-1)+t]*subp.α7[2*(j-1)+t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*mp.y[2*(dt.N["plant"]+k-1)+t]*subp.α8[2*(k-1)+t] for k=1:dt.N["distribution"] for t=1:2)+
                    sum(dt.Vij[i]*mp.uij[i]*subp.α9[i] for i in findnz(dt.Vij)[1])+
                    sum(dt.Vjk[i]*mp.ujk[i]*subp.α10[i] for i in findnz(dt.Vjk)[1])-
                    # sum(dt.Vkl[i]*ukl[i]*subp.α11[i] for i in findnz(dt.Vkl)[1])+
                    sum(dt.bigM.*mp.uij.*subp.α12) - sum(dt.bigM.*mp.ujk.*subp.α13) - sum(dt.bigM.*mp.ukl.*subp.α14))
                @constraint( molp, 0 ≥ sum(subp.α5[5*(l-1)+p]*dt.d[l,p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
                    sum(dt.N["cap"][j]*molp[:y][2*(j-1)+t]*subp.α7[2*(j-1)+t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*mp.y[2*(dt.N["plant"]+k-1)+t]*subp.α8[2*(k-1)+t] for k=1:dt.N["distribution"] for t=1:2)+
                    sum(dt.Vij[i]*molp[:uij][i]*subp.α9[i] for i in findnz(dt.Vij)[1])+
                    sum(dt.Vjk[i]*molp[:ujk][i]*subp.α10[i] for i in findnz(dt.Vjk)[1])-
                    # sum(dt.Vkl[i]*ukl[i]*subp.α11[i] for i in findnz(dt.Vkl)[1])+
                    sum(dt.bigM.*molp[:uij].*subp.α12) - sum(dt.bigM.*molp[:ujk].*subp.α13) - sum(dt.bigM.*molp[:ukl].*subp.α14))
                push!(fcuts, (α5=subp.α5,α6=subp.α6,α7=subp.α7,α8=subp.α8,α9=subp.α9,α10=subp.α10,α11=subp.α11,α12=subp.α12,α13=subp.α13,α14=subp.α14))
            end
        else #no solution found by FP
            optimize!(mp.m)
            @show st = termination_status(mp.m)
            yb = value.(mp.y); ubij = value.(mp.uij); ubjk = value.(mp.ujk); ubkl = value.(mp.ukl); θb = value(mp.θ)
            subp = solve_dsp(dsp,yb,ubij,ubjk,ubkl)
            if subp.res == :OptimalityCut
                @info "BD == Optimality cut found"
                if round(θb; digits=4) ≥ round(subp.obj; digits=4)
                    return
                else
                    cut = @constraint(mp.m, mp.θ ≥ sum(subp.α5[5*(l-1)+p]*dt.d[l,p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
                        sum(dt.N["cap"][j]*mp.y[2*(j-1)+t]*subp.α7[2*(j-1)+t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*mp.y[2*(dt.N["plant"]+k-1)+t]*subp.α8[2*(k-1)+t] for k=1:dt.N["distribution"] for t=1:2)+
                        sum(dt.Vij[i]*mp.uij[i]*subp.α9[i] for i in findnz(dt.Vij)[1])+
                        sum(dt.Vjk[i]*mp.ujk[i]*subp.α10[i] for i in findnz(dt.Vjk)[1])-
                        # sum(dt.Vkl[i]*ukl[i]*subp.α11[i] for i in findnz(dt.Vkl)[1])+
                        sum(dt.bigM.*mp.uij.*subp.α12) - sum(dt.bigM.*mp.ujk.*subp.α13) - sum(dt.bigM.*mp.ukl.*subp.α14))
                    @constraint( molp, molp[:θ] ≥ sum(subp.α5[5*(l-1)+p]*dt.d[l,p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
                        sum(dt.N["cap"][j]*molp[:y][2*(j-1)+t]*subp.α7[2*(j-1)+t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*mp.y[2*(dt.N["plant"]+k-1)+t]*subp.α8[2*(k-1)+t] for k=1:dt.N["distribution"] for t=1:2)+
                        sum(dt.Vij[i]*molp[:uij][i]*subp.α9[i] for i in findnz(dt.Vij)[1])+
                        sum(dt.Vjk[i]*molp[:ujk][i]*subp.α10[i] for i in findnz(dt.Vjk)[1])-
                        # sum(dt.Vkl[i]*ukl[i]*subp.α11[i] for i in findnz(dt.Vkl)[1])+
                        sum(dt.bigM.*molp[:uij].*subp.α12) - sum(dt.bigM.*molp[:ujk].*subp.α13) - sum(dt.bigM.*molp[:ukl].*subp.α14))
                    push!(ocuts, (α5=subp.α5,α6=subp.α6,α7=subp.α7,α8=subp.α8,α9=subp.α9,α10=subp.α10,α11=subp.α11,α12=subp.α12,α13=subp.α13,α14=subp.α14))
                end
            else
                @info "BD == Feasibility cut found"
                cut = @constraint(mp.m, 0  ≥ sum(subp.α5[5*(l-1)+p]*dt.d[l,p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
                    sum(dt.N["cap"][j]*mp.y[2*(j-1)+t]*subp.α7[2*(j-1)+t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*mp.y[2*(dt.N["plant"]+k-1)+t]*subp.α8[2*(k-1)+t] for k=1:dt.N["distribution"] for t=1:2)+
                    sum(dt.Vij[i]*mp.uij[i]*subp.α9[i] for i in findnz(dt.Vij)[1])+
                    sum(dt.Vjk[i]*mp.ujk[i]*subp.α10[i] for i in findnz(dt.Vjk)[1])-
                    # sum(dt.Vkl[i]*ukl[i]*subp.α11[i] for i in findnz(dt.Vkl)[1])+
                    sum(dt.bigM.*mp.uij.*subp.α12) - sum(dt.bigM.*mp.ujk.*subp.α13) - sum(dt.bigM.*mp.ukl.*subp.α14))
                @constraint( molp, 0 ≥ sum(subp.α5[5*(l-1)+p]*dt.d[l,p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
                    sum(dt.N["cap"][j]*molp[:y][2*(j-1)+t]*subp.α7[2*(j-1)+t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*mp.y[2*(dt.N["plant"]+k-1)+t]*subp.α8[2*(k-1)+t] for k=1:dt.N["distribution"] for t=1:2)+
                    sum(dt.Vij[i]*molp[:uij][i]*subp.α9[i] for i in findnz(dt.Vij)[1])+
                    sum(dt.Vjk[i]*molp[:ujk][i]*subp.α10[i] for i in findnz(dt.Vjk)[1])-
                    # sum(dt.Vkl[i]*ukl[i]*subp.α11[i] for i in findnz(dt.Vkl)[1])+
                    sum(dt.bigM.*molp[:uij].*subp.α12) - sum(dt.bigM.*molp[:ujk].*subp.α13) - sum(dt.bigM.*molp[:ukl].*subp.α14))
                push!(fcuts, (α5=subp.α5,α6=subp.α6,α7=subp.α7,α8=subp.α8,α9=subp.α9,α10=subp.α10,α11=subp.α11,α12=subp.α12,α13=subp.α13,α14=subp.α14))
            end
        end
    end
    return (mp,dsp,sol,ocuts,fcuts,Archiv)
end

mp = MasterP(w)
dsp = DualSP(w)

mp2,dsp2,sol,ocuts,fcuts,archiv = benders_decomposition(yt,utij,utjk,utkl)
1
