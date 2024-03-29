cd("C:/Users/AK121396/Desktop/ProjectBenders")
# activate .
using JuMP,CPLEX,LinearAlgebra,DelimitedFiles,CPUTime,SparseArrays
# import Printf SparseArrays,MathProgBase,
"If we use weighted sum for BD, obj values must be normalised."
# "Maybe we can use the weighted sum code by Gandibluex (Vopt)? and just solve sub problems with BD->no dichotomic search available for bi-obj mip"
# https://github.com/matbesancon/SimpleBenders.jl/blob/master/test/runtests.jl
# https://matbesancon.xyz/post/2019-05-08-simple-benders/ #matrix form
# https://co-at-work.zib.de/slides/Donnerstag_24.9/Benders_decomposition-Fundamentals.pdf

##############  Benders decomposition using Mathematical models  ###############
struct Data3
    filepath::String; N::Dict{}; d::Array{}; c::Array{};  Mij::Array{}; Mjk::Array{}; Mkl::Array{};
    gij::Array{}; gjk::Array{}; gkl::Array{}; vij::Array{}; vjk::Array{}; vkl::Array{}; rij::Array{}; rjk::Array{}; rkl::Array{};
    Vij::Array{}; Vjk::Array{}; Vkl::Array{}; b::Array{}; upl::Int; udc::Int; bigM::Int # e::Array{};q::Array{};
    function Data3(filepath)
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
# file = "F:scnd/Test1S2"
dt = Data3(file); # for Benders
##########################  Mathematical model  #########################
#Master Problem
struct build_master
    y::Matrix{VariableRef}
    uij::JuMP.Containers.SparseAxisArray{VariableRef}
    ujk::JuMP.Containers.SparseAxisArray{VariableRef}
    ukl::JuMP.Containers.SparseAxisArray{VariableRef}
    θ::VariableRef
    m::Model
end

function build_master(w)
    mas = Model(CPLEX.Optimizer);
    MOI.set(mas, MOI.RawOptimizerAttribute("CPXPARAM_MIP_Strategy_Search"), 1) # conventional branch-and-cut. without dynamic search
    MOI.set(mas, MOI.NumberOfThreads(), 1)
    set_silent(mas)
    MOI.NodeCount()
    @variable(mas, y[1:dt.N["plant"]+dt.N["distribution"],1:2], Bin);
    @variable(mas, uij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:dt.Mij[i,j]], Bin);
    @variable(mas, ujk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:dt.Mjk[j,k]], Bin);
    @variable(mas, ukl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:dt.Mkl[k,l]], Bin);
    @variable(mas, θ>= -1000);

    @constraint(mas,[j=1:dt.N["plant"]+dt.N["distribution"]], sum(y[j,:]) <= 1);
    @constraint(mas,[i=1:dt.N["supplier"],j=1:dt.N["plant"]], sum(uij[i,j,m] for m=1:dt.Mij[i,j]) <= 1);
    @constraint(mas,[j=1:dt.N["plant"],k=1:dt.N["distribution"]], sum(ujk[j,k,m] for m=1:dt.Mjk[j,k]) <= 1);
    @constraint(mas,[k=1:dt.N["distribution"],l=1:dt.N["customer"]], sum(ukl[k,l,m] for m=1:dt.Mkl[k,l]) <= 1);
    @constraint(mas, sum(y[j,t] for j=1:dt.N["plant"] for t=1:2) <= dt.upl);
    @constraint(mas, sum(y[k+dt.N["plant"],t] for k=1:dt.N["distribution"] for t=1:2) <= dt.udc);
    @objective(mas, Min, w[1]*(sum(dt.c[j][t]*y[j,t] for j=1:dt.N["plant"]+dt.N["distribution"] for t=1:2)+
        sum(dt.gij[i][j][m]*uij[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j]) +
        sum(dt.gjk[j][k][m]*ujk[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]) +
        sum(dt.gkl[k][l][m]*ukl[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])) + θ );
    return build_master(y,uij,ujk,ukl,θ,mas)
end

struct DualSP
    # data::SubProblemData
    α1::Matrix{VariableRef}
    α2::Matrix{VariableRef}
    α3::Matrix{VariableRef}
    α4::Matrix{VariableRef}
    α5::Matrix{VariableRef}
    α6::Vector{VariableRef}
    α7::Matrix{VariableRef}
    α8::Matrix{VariableRef}
    α9::JuMP.Containers.SparseAxisArray{VariableRef}
    α10::JuMP.Containers.SparseAxisArray{VariableRef}
    α11::JuMP.Containers.SparseAxisArray{VariableRef}
    α12::JuMP.Containers.SparseAxisArray{VariableRef}
    α13::JuMP.Containers.SparseAxisArray{VariableRef}
    α14::JuMP.Containers.SparseAxisArray{VariableRef}
    m::Model
end

function DualSP(w)
    sub = direct_model(CPLEX.Optimizer());
    set_optimizer_attribute(sub, "CPX_PARAM_REDUCE", 0); # presolve must be turned off
    MOI.set(sub, MOI.RawOptimizerAttribute("CPXPARAM_MIP_Strategy_Search"), 1) # conventional branch-and-cut
    set_silent(sub) # display is turned off
    MOI.set(sub, MOI.NumberOfThreads(), 1) # of threads

    @variables(sub, begin
        α1[1:dt.N["plant"],1:5]
        α2[1:dt.N["distribution"],1:5]
        α3[1:dt.N["plant"],1:5]
        α4[1:dt.N["distribution"],1:5]
    end)
    @variable(sub, α5[1:dt.N["customer"],1:5] >= 0);
    @variable(sub, α6[1:dt.N["supplier"]] >= 0);
    @variable(sub, α7[1:dt.N["plant"],1:2] >= 0);
    @variable(sub, α8[1:dt.N["distribution"],1:2] >= 0);
    @variable(sub, α9[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:dt.Mij[i,j]] >= 0);
    @variable(sub, α10[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:dt.Mjk[j,k]] >= 0);
    @variable(sub, α11[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:dt.Mkl[k,l]] >= 0);
    @variable(sub, α12[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:dt.Mij[i,j]] >= 0);
    @variable(sub, α13[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:dt.Mjk[j,k]] >= 0);
    @variable(sub, α14[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:dt.Mkl[k,l]] >= 0);

    @constraint(sub, [i=1:dt.N["supplier"],j=1:dt.N["plant"],m=1:dt.Mij[i,j],p=1:5], α1[j,p]-α3[j,p]-α6[i]+α9[i,j,m]-α12[i,j,m] <= w[1]*(dt.N["vcs"][i][p]+dt.vij[i][j][m][p])+w[2]*(dt.b[i,p]+dt.rij[i][j][m][p]));
    @constraint(sub, [j=1:dt.N["plant"],k=1:dt.N["distribution"],m=1:dt.Mjk[j,k],p=1:5], -α1[j,p]+α2[k,p]-α4[k,p]+α10[j,k,m]-α13[j,k,m] <= w[1]*(dt.vjk[j][k][m][p])+w[2]*(dt.rjk[j][k][m][p]));
    @constraint(sub, [k=1:dt.N["distribution"],l=1:dt.N["customer"],m=1:dt.Mkl[k,l],p=1:5], -α2[k,p]+α5[l,p]+α11[k,l,m]-α14[k,l,m] <= w[1]*(dt.vkl[k][l][m][p])+w[2]*(dt.rkl[k][l][m][p]));
    @constraint(sub, [j=1:dt.N["plant"],t=1:2,p=1:5], α3[j,p]-α7[j,t] <= w[1]*(dt.N["vcp"][j][2*(t-1)+p])+w[2]*(dt.N["vep"][j][2*(t-1)+p]));
    @constraint(sub, [k=1:dt.N["distribution"],t=1:2,p=1:5], α4[k,p]-α8[k,t] <= w[1]*(dt.N["vcd"][k][2*(t-1)+p])+w[2]*(dt.N["ved"][k][2*(t-1)+p]));
    return DualSP(α1,α2,α3,α4,α5,α6,α7,α8,α9,α10,α11,α12,α13,α14,sub)
end
function solve_dsp(dsp::DualSP,yb,ubij,ubjk,ubkl)
    @objective(dsp.m, Max, sum(dsp.α5[l,p]*dt.d[l][p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*dsp.α6[i] for i=1:dt.N["supplier"])-
        sum(dt.N["cap"][j]*yb[j,t]*dsp.α7[j,t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*yb[dt.N["plant"]+k,t]*dsp.α8[k,t] for k=1:dt.N["distribution"] for t=1:2)+
        sum(dt.Vij[i][j][m]*ubij[i,j,m]*dsp.α9[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j]) +
        sum(dt.Vjk[j][k][m]*ubjk[j,k,m]*dsp.α10[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]) +
        # sum(dt.Vkl[k][l][m]*ubkl[k,l,m]*dsp.α11[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])-
        sum(dt.bigM*ubij[i,j,m]*dsp.α12[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j]) -
        sum(dt.bigM*ubjk[j,k,m]*dsp.α13[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]) -
        sum(dt.bigM*ubkl[k,l,m]*dsp.α14[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])
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
function lazy_callback(cb_data)
    yb = callback_value.(cb_data, mp.y);    ubij = callback_value.(cb_data, mp.uij);    ubjk = callback_value.(cb_data, mp.ujk)
    ubkl = callback_value.(cb_data, mp.ukl);    θb = callback_value(cb_data, mp.θ)
    subp = solve_dsp(dsp,yb,ubij,ubjk,ubkl)
    if subp.res == :OptimalityCut
        if round(θb; digits=4) ≥ round(subp.obj; digits=4)
            return
        else
            cut = @build_constraint( mp.θ ≥ sum(subp.α5[l,p]*dt.d[l][p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
                sum(dt.N["cap"][j]*mp.y[j,t]*subp.α7[j,t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*mp.y[dt.N["plant"]+k,t]*subp.α8[k,t] for k=1:dt.N["distribution"] for t=1:2)+
                sum(dt.Vij[i][j][m]*mp.uij[i,j,m]*subp.α9[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j]) +
                sum(dt.Vjk[j][k][m]*mp.ujk[j,k,m]*subp.α10[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k])
                # sum(dt.Vkl[k][l][m]*mp.ukl[k,l,m]*subp.α11[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l]) -
                - sum(dt.bigM*mp.uij[i,j,m]*subp.α12[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j]) -
                sum(dt.bigM*mp.ujk[j,k,m]*subp.α13[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]) -
                sum(dt.bigM*mp.ukl[k,l,m]*subp.α14[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])
                );
            MOI.submit(mp.m, MOI.LazyConstraint(cb_data), cut)
        end
    else
        cut = @build_constraint( 0 ≥ sum(subp.α5[l,p]*dt.d[l][p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
            sum(dt.N["cap"][j]*mp.y[j,t]*subp.α7[j,t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*mp.y[dt.N["plant"]+k,t]*subp.α8[k,t] for k=1:dt.N["distribution"] for t=1:2)+
            sum(dt.Vij[i][j][m]*mp.uij[i,j,m]*subp.α9[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j]) +
            sum(dt.Vjk[j][k][m]*mp.ujk[j,k,m]*subp.α10[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k])
            # sum(dt.Vkl[k][l][m]*mp.ukl[k,l,m]*subp.α11[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])-
            - sum(dt.bigM*mp.uij[i,j,m]*subp.α12[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j]) -
            sum(dt.bigM*mp.ujk[j,k,m]*subp.α13[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]) -
            sum(dt.bigM*mp.ukl[k,l,m]*subp.α14[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])
            );
        # push!(feasicuts, cut)
        MOI.submit(mp.m, MOI.LazyConstraint(cb_data), cut)
        push!(feasicuts, (α5=subp.α5,α6=subp.α6,α7=subp.α7,α8=subp.α8,α9=subp.α9,α10=subp.α10,α11=subp.α11,α12=subp.α12,α13=subp.α13,α14=subp.α14))
    end
    return
end
function everyfrac_callback(cb_data)
    yb = callback_value.(cb_data, mp.y);    ubij = callback_value.(cb_data, mp.uij);    ubjk = callback_value.(cb_data, mp.ujk)
    ubkl = callback_value.(cb_data, mp.ukl);    θb = callback_value(cb_data, mp.θ)
    subp = solve_dsp(dsp,yb,ubij,ubjk,ubkl)
    if subp.res == :OptimalityCut
        if round(θb; digits=4) ≥ round(subp.obj; digits=4)
            return
        else
            cut = @build_constraint( mp.θ ≥ sum(subp.α5[l,p]*dt.d[l][p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
                sum(dt.N["cap"][j]*mp.y[j,t]*subp.α7[j,t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*mp.y[dt.N["plant"]+k,t]*subp.α8[k,t] for k=1:dt.N["distribution"] for t=1:2)+
                sum(dt.Vij[i][j][m]*mp.uij[i,j,m]*subp.α9[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j]) +
                sum(dt.Vjk[j][k][m]*mp.ujk[j,k,m]*subp.α10[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k])
                # sum(dt.Vkl[k][l][m]*mp.ukl[k,l,m]*subp.α11[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l]) -
                -sum(dt.bigM*mp.uij[i,j,m]*subp.α12[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j]) -
                sum(dt.bigM*mp.ujk[j,k,m]*subp.α13[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]) -
                sum(dt.bigM*mp.ukl[k,l,m]*subp.α14[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])
                );
            MOI.submit(mp.m, MOI.UserCut(cb_data), cut)
        end
    else
        cut = @build_constraint( 0 ≥ sum(subp.α5[l,p]*dt.d[l][p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
            sum(dt.N["cap"][j]*mp.y[j,t]*subp.α7[j,t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*mp.y[dt.N["plant"]+k,t]*subp.α8[k,t] for k=1:dt.N["distribution"] for t=1:2)+
            sum(dt.Vij[i][j][m]*mp.uij[i,j,m]*subp.α9[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j]) +
            sum(dt.Vjk[j][k][m]*mp.ujk[j,k,m]*subp.α10[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k])
            # sum(dt.Vkl[k][l][m]*mp.ukl[k,l,m]*subp.α11[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])-
            - sum(dt.bigM*mp.uij[i,j,m]*subp.α12[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j]) -
            sum(dt.bigM*mp.ujk[j,k,m]*subp.α13[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]) -
            sum(dt.bigM*mp.ukl[k,l,m]*subp.α14[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])
            );
        MOI.submit(mas.m, MOI.UserCut(cb_data), cut)
        # push!(mas.feasicuts, cut)
        #Should we keep those Feasibility cuts as well?
        # push!(feasicuts, (α5=subp.α5,α6=subp.α6,α7=subp.α7,α8=subp.α8,α9=subp.α9,α10=subp.α10,α11=subp.α11,α12=subp.α12,α13=subp.α13,α14=subp.α14))
    end
    return
end

function rootfrac_callback(cb_data)
    ndepth = Ref{CPXLONG}()
    CPXcallbackgetinfolong(cb_data, CPXCALLBACKINFO_NODEDEPTH, ndepth)
    if ndepth[] == 0
        yb = callback_value.(cb_data, mp.y);    ubij = callback_value.(cb_data, mp.uij);    ubjk = callback_value.(cb_data, mp.ujk)
        ubkl = callback_value.(cb_data, mp.ukl);    θb = callback_value(cb_data, mp.θ)
        subp = solve_dsp(dsp,yb,ubij,ubjk,ubkl)
        if subp.res == :OptimalityCut
            if round(θb; digits=4) ≥ round(subp.obj; digits=4)
                println("finish")
                return
            else
                cut = @build_constraint( mp.θ ≥ sum(subp.α5[l,p]*dt.d[l][p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
                    sum(dt.N["cap"][j]*mp.y[j,t]*subp.α7[j,t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*mp.y[dt.N["plant"]+k,t]*subp.α8[k,t] for k=1:dt.N["distribution"] for t=1:2)+
                    sum(dt.Vij[i][j][m]*mp.uij[i,j,m]*subp.α9[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j]) +
                    sum(dt.Vjk[j][k][m]*mp.ujk[j,k,m]*subp.α10[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k])
                    # sum(dt.Vkl[k][l][m]*mp.ukl[k,l,m]*subp.α11[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l]) -
                    -sum(dt.bigM*mp.uij[i,j,m]*subp.α12[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j]) -
                    sum(dt.bigM*mp.ujk[j,k,m]*subp.α13[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]) -
                    sum(dt.bigM*mp.ukl[k,l,m]*subp.α14[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])
                    );
                MOI.submit(mp.m, MOI.UserCut(cb_data), cut)
                push!(optcuts, (α5=subp.α5,α6=subp.α6,α7=subp.α7,α8=subp.α8,α9=subp.α9,α10=subp.α10,α11=subp.α11,α12=subp.α12,α13=subp.α13,α14=subp.α14))

            end
        else
            cut = @build_constraint( 0 ≥ sum(subp.α5[l,p]*dt.d[l][p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
                sum(dt.N["cap"][j]*mp.y[j,t]*subp.α7[j,t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*mp.y[dt.N["plant"]+k,t]*subp.α8[k,t] for k=1:dt.N["distribution"] for t=1:2)+
                sum(dt.Vij[i][j][m]*mp.uij[i,j,m]*subp.α9[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j]) +
                sum(dt.Vjk[j][k][m]*mp.ujk[j,k,m]*subp.α10[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k])
                # sum(dt.Vkl[k][l][m]*mp.ukl[k,l,m]*subp.α11[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])-
                - sum(dt.bigM*mp.uij[i,j,m]*subp.α12[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j]) -
                sum(dt.bigM*mp.ujk[j,k,m]*subp.α13[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]) -
                sum(dt.bigM*mp.ukl[k,l,m]*subp.α14[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])
                );
            MOI.submit(mp.m, MOI.UserCut(cb_data), cut)
            # push!(mas.feasicuts, cut)
            push!(feasicuts, (α5=subp.α5,α6=subp.α6,α7=subp.α7,α8=subp.α8,α9=subp.α9,α10=subp.α10,α11=subp.α11,α12=subp.α12,α13=subp.α13,α14=subp.α14))
        end
        return
    end
end


w = [0.5,0.5]; feasicuts = []; optcuts = []; nopt_cons=0; nfeasi_cons=0; cint = 0; cfrac = 0;
mp = build_master(w);
dsp = DualSP(w);
MOI.set(mp.m, MOI.LazyConstraintCallback(), lazy_callback)
# MOI.set(mp.m, MOI.UserCutCallback(), everyfrac_callback)
MOI.set(mp.m, MOI.UserCutCallback(), rootfrac_callback)
optimize!(mp.m)
# termination_status(mp.m)
node_count(mp.m),solve_time(mp.m)
# nopt,nfeasi
value.(mp.y)
1
# function benders_with_callback(cb_data)#benders_decomposition(w,mas::Model, y::Matrix{VariableRef},uij::JuMP.Containers.SparseAxisArray{VariableRef},ujk::JuMP.Containers.SparseAxisArray{VariableRef},ukl::JuMP.Containers.SparseAxisArray{VariableRef})
#     # ,uij::Array{VariableRef},ujk::Array{VariableRef},ukl::Array{VariableRef})
#     yb = callback_value.(cb_data, mp.y);    ubij = callback_value.(cb_data, mp.uij);    ubjk = callback_value.(cb_data, mp.ujk)
#     ubkl = callback_value.(cb_data, mp.ukl);    θb = callback_value(cb_data, mp.θ)
#     subp = solve_dsp(dsp,yb,ubij,ubjk,ubkl)
#
#     if subp.res == :OptimalityCut
#         # @info "Optimality cut found"
#         # println(θb,"  and  ",round(subp.obj; digits=4))
#         if round(θb; digits=4) ≥ round(subp.obj; digits=4)
#             return
#         else
#             global nopt_cons+=1
#             cut = @build_constraint( mp.θ ≥ sum(subp.α5[l,p]*dt.d[l][p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
#                 sum(dt.N["cap"][j]*mp.y[j,t]*subp.α7[j,t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*mp.y[dt.N["plant"]+k,t]*subp.α8[k,t] for k=1:dt.N["distribution"] for t=1:2)+
#                 sum(dt.Vij[i][j][m]*mp.uij[i,j,m]*subp.α9[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j]) +
#                 sum(dt.Vjk[j][k][m]*mp.ujk[j,k,m]*subp.α10[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]) +
#                 sum(dt.Vkl[k][l][m]*mp.ukl[k,l,m]*subp.α11[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l]) -
#                 sum(dt.bigM*mp.uij[i,j,m]*subp.α12[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j]) -
#                 sum(dt.bigM*mp.ujk[j,k,m]*subp.α13[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]) -
#                 sum(dt.bigM*mp.ukl[k,l,m]*subp.α14[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])
#                 );
#             MOI.submit(mp.m, MOI.LazyConstraint(cb_data), cut)
#         end
#     else
#         # @info "Feasibility cut found"
#         global nfeasi_cons += 1
#         cut = @build_constraint( 0 ≥ sum(subp.α5[l,p]*dt.d[l][p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
#             sum(dt.N["cap"][j]*mp.y[j,t]*subp.α7[j,t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*mp.y[dt.N["plant"]+k,t]*subp.α8[k,t] for k=1:dt.N["distribution"] for t=1:2)+
#             sum(dt.Vij[i][j][m]*mp.uij[i,j,m]*subp.α9[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j]) +
#             sum(dt.Vjk[j][k][m]*mp.ujk[j,k,m]*subp.α10[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]) +
#             sum(dt.Vkl[k][l][m]*mp.ukl[k,l,m]*subp.α11[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])-
#             sum(dt.bigM*mp.uij[i,j,m]*subp.α12[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j]) -
#             sum(dt.bigM*mp.ujk[j,k,m]*subp.α13[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]) -
#             sum(dt.bigM*mp.ukl[k,l,m]*subp.α14[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])
#             );
#         # push!(feasicuts, cut)
#         MOI.submit(mp.m, MOI.LazyConstraint(cb_data), cut)
#         push!(feasicuts, (α5=subp.α5,α6=subp.α6,α7=subp.α7,α8=subp.α8,α9=subp.α9,α10=subp.α10,α11=subp.α11,α12=subp.α12,α13=subp.α13,α14=subp.α14))
#     end
#     # @info "Addin g the cut $(cut)"
#     # end
#     return
# end
function benders_decomposition(w,mas::Model, y::Matrix{VariableRef},uij::JuMP.Containers.SparseAxisArray{VariableRef},ujk::JuMP.Containers.SparseAxisArray{VariableRef},ukl::JuMP.Containers.SparseAxisArray{VariableRef})
    # ,uij::Array{VariableRef},ujk::Array{VariableRef},ukl::Array{VariableRef})
    sub = direct_model(CPLEX.Optimizer()); dsp = DualSubP(sub,w)

    @variable(mas, θ>= -1000);
    @objective(mas, Min, w[1]*(sum(dt.c[j][t]*y[j,t] for j=1:dt.N["plant"]+dt.N["distribution"] for t=1:2)+
        sum(dt.gij[i][j][m]*uij[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j]) +
        sum(dt.gjk[j][k][m]*ujk[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]) +
        sum(dt.gkl[k][l][m]*ukl[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])) + θ );
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
            if round(θ1; digits=4) ≥ round(subp.obj; digits=4)
                break
            else
                nopt_cons+=1
                cut = @constraint( mas, θ ≥ sum(subp.α5[l,p]*dt.d[l][p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
                    sum(dt.N["cap"][j]*y[j,t]*subp.α7[j,t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*y[dt.N["plant"]+k,t]*subp.α8[k,t] for k=1:dt.N["distribution"] for t=1:2)+
                    sum(dt.Vij[i][j][m]*uij[i,j,m]*subp.α9[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j]) +
                    sum(dt.Vjk[j][k][m]*ujk[j,k,m]*subp.α10[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]) +
                    # sum(dt.Vkl[k][l][m]*ukl[k,l,m]*subp.α11[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l]) -
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
                # sum(dt.Vkl[k][l][m]*ukl[k,l,m]*subp.α11[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])-
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
Sub problem (primal)
sub = Model(CPLEX.Optimizer);
set_optimizer_attribute(sub, "CPX_PARAM_REDUCE", 0)
MOI.set(sub, MOI.NumberOfThreads(), 1);set_silent(sub)
@variable(sub, 0<= xij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:dt.Mij[i,j],1:5] );
@variable(sub, 0<= xjk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:dt.Mjk[j,k],1:5] );
@variable(sub, 0<= xkl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:dt.Mkl[k,l],1:5] );
@variable(sub, 0<= h[1:dt.N["plant"]+dt.N["distribution"],1:5,1:2] );
@variable(sub, yb[1:dt.N["plant"]+dt.N["distribution"],1:2], Bin);
@variable(sub, ubij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:dt.Mij[i,j]], Bin);
@variable(sub, ubjk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:dt.Mjk[j,k]], Bin);
@variable(sub, ubkl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:dt.Mkl[k,l]], Bin);
@constraints(sub, begin
    [j=1:dt.N["plant"],p=1:5], sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j]) == sum(xjk[j,k,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k])
    [k=1:dt.N["distribution"],p=1:5], sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k]) == sum(xkl[k,l,m,p] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])
end)
@constraints(sub, begin
    [j=1:dt.N["plant"],p=1:5], sum(h[j,p,:]) == sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j])
    [k=1:dt.N["distribution"],p=1:5], sum(h[k+dt.N["plant"],p,:] ) == sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k])
    [l=1:dt.N["customer"],p=1:5], sum(xkl[k,l,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mkl[k,l]) >= dt.d[l][p]
end )
con6 = @constraint(sub,[i=1:dt.N["supplier"]], sum(xij[i,j,m,p] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) <= dt.N["cas"][i]);
con7 = @constraint(sub,[j=1:dt.N["plant"], t=1:2], sum(h[j,1:5,t]) <=dt.N["cap"][j]*yb[j,t]);
con8 = @constraint(sub,[k=1:J+dt.N["distribution"], t=1:2], sum(h[k,1:5,t]) >= dt.N["cad"][k]*yb[J+k,t]);
con9 = @constraint(sub,[i=1:dt.N["supplier"], j=1:dt.N["plant"], m=1:dt.Mij[i,j]], sum(xij[i,j,m,p] for p=1:5) >= dt.Vij[i][j][m]*ubij[i,j,m] );
con10 = @constraint(sub,[j=1:dt.N["plant"], k=1:dt.N["distribution"], m=1:dt.Mjk[j,k]], sum(xjk[j,k,m,p] for p=1:5) >= dt.Vjk[j][k][m]*ubjk[j,k,m]);
# con11 = @constraint(sub,[k=1:dt.N["distribution"], l=1:dt.N["customer"], m=1:dt.Mkl[k,l]], sum(xkl[k,l,m,p] for p=1:5) >= dt.Vkl[k][l][m]*ubkl[k,l,m]);
@objective( sub, Min, w[1]*( exa + sum(dt.e[j][(p-1)*2+t]*h[j,p,t] for j=1:dt.N["plant"] + dt.N["distribution"] for p=1:5 for t=1:2)) +
    exv + w[2]*( exb + sum(dt.q[j][(p-1)*2+t]*h[j,p,t] for j=1:dt.N["plant"]+dt.N["distribution"] for p=1:5 for t=1:2) +exr) );

mp = build_master(w);
dsp = DualSP(w);
MOI.set(mp.m, MOI.LazyConstraintCallback(), benders_with_callback)
node_count(mp.m),solve_time(mp.m)


mas = Model(CPLEX.Optimizer); set_silent(mas)
MOI.set(mas, MOI.NumberOfThreads(), 1);set_silent(mas)
@variable(mas, y[1:dt.N["plant"]+dt.N["distribution"],1:2], Bin);
@variable(mas, uij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:dt.Mij[i,j]], Bin);
@variable(mas, ujk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:dt.Mjk[j,k]], Bin);
@variable(mas, ukl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:dt.Mkl[k,l]], Bin);
@variable(mas, θ>= -1000);
@constraint(mas,[j=1:dt.N["plant"]+dt.N["distribution"]], sum(y[j,:]) <= 1);
@constraint(mas,[i=1:dt.N["supplier"],j=1:dt.N["plant"]], sum(uij[i,j,m] for m=1:dt.Mij[i,j]) <= 1);
@constraint(mas,[j=1:dt.N["plant"],k=1:dt.N["distribution"]], sum(ujk[j,k,m] for m=1:dt.Mjk[j,k]) <= 1);
@constraint(mas,[k=1:dt.N["distribution"],l=1:dt.N["customer"]], sum(ukl[k,l,m] for m=1:dt.Mkl[k,l]) <= 1);
@constraint(mas, sum(y[j,t] for j=1:dt.N["plant"] for t=1:2) <= dt.upl);
@constraint(mas, sum(y[k+dt.N["plant"],t] for k=1:dt.N["distribution"] for t=1:2) <= dt.udc);
@objective(mas, Min, w[1]*(sum(dt.c[j][t]*y[j,t] for j=1:dt.N["plant"]+dt.N["distribution"] for t=1:2)+
    sum(dt.gij[i][j][m]*uij[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j]) +
    sum(dt.gjk[j][k][m]*ujk[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]) +
    sum(dt.gkl[k][l][m]*ukl[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])) + θ );

#####obj calculator
y = round.(value.(scnd1[:y1]));
uij = round.(value.(scnd1[:uij1]));
ujk = round.(value.(scnd1[:ujk1]));
ukl = round.(value.(scnd1[:ukl1]));
xij = round.(value.(scnd1[:xij1]); digits=2);
xjk = round.(value.(scnd1[:xjk1]); digits=2);
xkl = round.(value.(scnd1[:xkl1]); digits=2);
h = round.(value.(scnd1[:h1]); digits=2);

y = value.(scnd1[:y1]);
uij = value.(scnd1[:uij1]);
ujk = value.(scnd1[:ujk1]);
ukl = value.(scnd1[:ukl1]);
xij = value.(scnd1[:xij1]);
xjk = value.(scnd1[:xjk1]);
xkl = value.(scnd1[:xkl1]);
h = value.(scnd1[:h1]) ;

obj1 = sum(dt1.c.*y)+sum(repeat(dt1.a[1,:], outer=sum(dt1.Mij[1,:])).*xij[1:sum(dt1.Mij[1,:])*5])+
        sum(sum(repeat(dt1.a[i,:], outer=sum(dt1.Mij[i,:])).*xij[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"])+
        sum(dt1.e.*h) + sum(dt1.gij[i]*uij[i] for i in findnz(dt1.gij)[1]) + sum(dt1.gjk[i]*ujk[i] for i in findnz(dt1.gjk)[1]) + sum(dt1.gkl[i].*ukl[i] for i in findnz(dt1.gkl)[1])+
        sum(dt1.vij.*xij)+sum(dt1.vjk.*xjk)+sum(dt1.vkl.*xkl)
obj2 = sum(repeat(dt1.b[1,:], outer=sum(dt1.Mij[1,:])).*xij[1:sum(dt1.Mij[1,:])*5]) +
        sum(sum(repeat(dt1.b[i,:], outer=sum(dt1.Mij[i,:])).*xij[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"]) +
        sum(dt1.q.*h) + sum(dt1.rij.*xij)+sum(dt1.rjk.*xjk)+sum(dt1.rkl.*xkl)




y = value.(mymodel[:y]);
uij = value.(mymodel[:uij])
ujk = value.(mymodel[:ujk]);
ukl = value.(mymodel[:ukl]);
xij = value.(mymodel[:xij]);
xjk = value.(mymodel[:xjk]);
xkl = value.(mymodel[:xkl]) ;
h = value.(mymodel[:h]) ;
sum(round.(reshape(h,270,1)))

hcat(h)

sum(value.(all_variables(scnd)))

y = round.(value.(mymodel[:y]))
findall(i->0<i<1,value.(mymodel[:ukl]))

uij = round.(value.(mymodel[:uij]))
ujk = round.(value.(mymodel[:ujk]))
ukl = round.(value.(mymodel[:ukl]))
xij = round.(value.(mymodel[:xij]); digits=4)
xjk = round.(value.(mymodel[:xjk]); digits=4)
xkl = round.(value.(mymodel[:xkl]); digits=4)
h = round.(value.(mymodel[:h]); digits=4)

exg = 0;
for i=1:dt.N["supplier"]
    for j=1:dt.N["plant"]
        exg = exg + (10000*uij[i,j,1]);
    end
end
for j=1:dt.N["plant"]
    for k=1:dt.N["distribution"]
        exg = exg + (10000*ujk[j,k,1]);
    end
end
for k=1:dt.N["distribution"]
    for l=1:dt.N["customer"]
        exg = exg + (10000*ukl[k,l,1]);
    end
end


obj1 =   sum(dt.c[j][t]*y[j,t] for j=1:1:dt.N["plant"]+dt.N["distribution"] for t=1:2) + exg +
        sum(dt.N["vcs"][i][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) +
        sum(dt.vij[i][j][m][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
        sum(dt.vjk[j][k][m][p]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
        sum(dt.vkl[k][l][m][p]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5)+
        sum(dt.N["vcp"][j][2*(p-1)+t]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
        sum(dt.N["vcd"][k][2*(p-1)+t]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2)


obj1 = sum(dt.c[j][t]*y[j,t] for j=1:1:dt.N["plant"]+dt.N["distribution"] for t=1:2) + #exg +
    sum(dt.gij[i][j][m]*uij[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j])+
    sum(dt.gjk[j][k][m]*ujk[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]) +
    sum(dt.gkl[k][l][m]*ukl[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])+
    sum(dt.N["vcs"][i][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) +
    sum(dt.vij[i][j][m][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
    sum(dt.vjk[j][k][m][p]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
    sum(dt.vkl[k][l][m][p]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5)+
    sum(dt.N["vcp"][j][2*(p-1)+t]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
    sum(dt.N["vcd"][k][2*(p-1)+t]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2)

obj2 =  sum(dt.b[i,p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
    sum(dt.N["vep"][j][2*(p-1)+t]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
    sum(dt.N["ved"][k][2*(p-1)+t]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2)+
    sum(dt.rij[i][j][m]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
    sum(dt.rjk[j][k][m]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
    sum(dt.rkl[k][l][m]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5)

################

