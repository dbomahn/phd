using JuMP,CPLEX,LinearAlgebra,MathProgBase,DelimitedFiles,SparseArrays,CPUTime
# import Printf
"If we use weighted sum for BD, obj values must be normalised"
"Maybe we can use the weighted sum code by Gandibluex (Vopt)? and just solve sub problems with BD"
# https://github.com/matbesancon/SimpleBenders.jl/blob/master/test/runtests.jl
# https://matbesancon.xyz/post/2019-05-08-simple-benders/ #matrix form
# https://co-at-work.zib.de/slides/Donnerstag_24.9/Benders_decomposition-Fundamentals.pdf

##############  Benders decomposition using Mathematical models  ###############
mutable struct Data
    filepath::String; N::Dict{}; d::Array{}; c::Array{};  Mij::Array{}; Mjk::Array{}; Mkl::Array{};
    gij::Array{}; gjk::Array{}; gkl::Array{}; vij::Array{}; vjk::Array{}; vkl::Array{}; rij::Array{}; rjk::Array{}; rkl::Array{};
    Vij::Array{}; Vjk::Array{}; Vkl::Array{}; b::Array{}; upl::Int; udc::Int; bigM::Int # e::Array{};q::Array{};
    function Data(filepath)
        dt = readdlm(filepath);
        notafile = readdlm("/home/ak121396/Desktop/instances/SCND/Notations.txt", '=');
        # notafile = readdlm("E:/scnd/Notations.txt", '=');
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
dt = Data(file);
##########################  Mathematical model  #########################
#Master Problem
struct build_master
    y::Matrix{VariableRef}
    uij::JuMP.Containers.SparseAxisArray{VariableRef}
    ujk::JuMP.Containers.SparseAxisArray{VariableRef}
    ukl::JuMP.Containers.SparseAxisArray{VariableRef}
    m::Model
end

function build_master()
    mas = Model(CPLEX.Optimizer); set_silent(mas)
    MOI.set(mas, MOI.NumberOfThreads(), 1);set_silent(mas)
    @variable(mas, y[1:dt.N["plant"]+dt.N["distribution"],1:2], Bin);
    @variable(mas, uij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:dt.Mij[i,j]], Bin);
    @variable(mas, ujk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:dt.Mjk[j,k]], Bin);
    @variable(mas, ukl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:dt.Mkl[k,l]], Bin);
    @constraint(mas,[j=1:dt.N["plant"]+dt.N["distribution"]], sum(y[j,:]) <= 1);
    @constraint(mas,[i=1:dt.N["supplier"],j=1:dt.N["plant"]], sum(uij[i,j,m] for m=1:dt.Mij[i,j]) <= 1);
    @constraint(mas,[j=1:dt.N["plant"],k=1:dt.N["distribution"]], sum(ujk[j,k,m] for m=1:dt.Mjk[j,k]) <= 1);
    @constraint(mas,[k=1:dt.N["distribution"],l=1:dt.N["customer"]], sum(ukl[k,l,m] for m=1:dt.Mkl[k,l]) <= 1);
    @constraint(mas, sum(y[j,t] for j=1:dt.N["plant"] for t=1:2) <= dt.upl);
    @constraint(mas, sum(y[k+dt.N["plant"],t] for k=1:dt.N["distribution"] for t=1:2) <= dt.udc);
    return build_master(y,uij,ujk,ukl,mas)
end

struct DualSubP
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

function DualSubP(sub::Model,w)
    set_optimizer_attribute(sub, "CPX_PARAM_REDUCE", 0); # turning presolve off
    set_silent(sub) #MOI.set(sub, MOI.RawParameter("CPX_PARAM_SCRIND"), false )
    MOI.set(sub, MOI.NumberOfThreads(), 1) # MOI.set(sub, MOI.RawParameter("CPX_PARAM_THREADS"),1  )

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
    return DualSubP(α1,α2,α3,α4,α5,α6,α7,α8,α9,α10,α11,α12,α13,α14,sub)
end
function solve_dsp(dsp::DualSubP,yb,ubij,ubjk,ubkl)
    @objective(dsp.m, Max, sum(dsp.α5[l,p]*dt.d[l][p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*dsp.α6[i] for i=1:dt.N["supplier"])-
        sum(dt.N["cap"][j]*yb[j,t]*dsp.α7[j,t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*yb[dt.N["plant"]+k,t]*dsp.α8[k,t] for k=1:dt.N["distribution"] for t=1:2)+
        sum(dt.Vij[i][j][m]*ubij[i,j,m]*dsp.α9[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j]) +
        sum(dt.Vjk[j][k][m]*ubjk[j,k,m]*dsp.α10[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]) +
        sum(dt.Vkl[k][l][m]*ubkl[k,l,m]*dsp.α11[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])-
        sum(dt.bigM*ubij[i,j,m]*dsp.α12[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j]) -
        sum(dt.bigM*ubjk[j,k,m]*dsp.α13[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]) -
        sum(dt.bigM*ubkl[k,l,m]*dsp.α14[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])
        )
    optimize!(dsp.m)
    st = termination_status(dsp.m)

    if st == MOI.OPTIMAL
      return (res = :OptimalityCut, obj = objective_value(dsp.m), α5 = value.(dsp.α5), α6 = value.(dsp.α6), α7 = value.(dsp.α7), α8= value.(dsp.α8), α9= value.(dsp.α9), α10= value.(dsp.α10),α11= value.(dsp.α11), α12= value.(dsp.α12),α13= value.(dsp.α13),α14= value.(dsp.α14) )
    elseif st == MOI.DUAL_INFEASIBLE
        # lens = [length(dsp.α5),length(dsp.α6),length(dsp.α7),length(dsp.α8),length(dsp.α9),length(dsp.α10),length(dsp.α11),length(dsp.α1),length(dsp.α2),length(dsp.α3),length(dsp.α4)]
        # rays = Vector{Float64}(undef,sum(dsp.lens))
        # CPXgetray(backend(dsp.m).env, backend(dsp.m).lp, rays) # vr = value.(rays)
        # pvci = [sum(lens[1:l]) for l=1:length(lens)]
        # insert!(pvci,1,0)
        # slices = [rays[pvci[l]+1:pvci[l+1]] for l=1:length(pvci)-1]
        # return ( res = :FeasibilityCut, α5=slices[8], α6=slices[1], α7=slices[10], α8=slices[11], α9=slices[2], α10=slices[5], α11=slices[9] )
        return ( res = :FeasibilityCut, α5 = value.(dsp.α5), α6 = value.(dsp.α6), α7 = value.(dsp.α7), α8= value.(dsp.α8), α9= value.(dsp.α9), α10= value.(dsp.α10), α11= value.(dsp.α11), α12= value.(dsp.α12),α13= value.(dsp.α13),α14= value.(dsp.α14) )
    else
      error("DualSubProblem error: status $st")
    end
end


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
w = [0.5,0.5];
mas = build_master();
@CPUtime fmodel,fy,fuij,fujk,fukl,noptcut,nfeasicut,feasicuts,thetas = benders_decomposition(w,mas.m,mas.y,mas.uij,mas.ujk,mas.ukl)
value.(fy)
sum(value.(fuij))
sum(value.(fujk))
feasicuts[2][3]
nfeasisut[1]
noptcut
202860
1

#Sub problem (primal)
# sub = Model(CPLEX.Optimizer);
# set_optimizer_attribute(sub, "CPX_PARAM_REDUCE", 0)
# MOI.set(sub, MOI.NumberOfThreads(), 1);set_silent(sub)
# @variable(sub, 0<= xij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:dt.Mij[i,j],1:5] );
# @variable(sub, 0<= xjk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:dt.Mjk[j,k],1:5] );
# @variable(sub, 0<= xkl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:dt.Mkl[k,l],1:5] );
# @variable(sub, 0<= h[1:dt.N["plant"]+dt.N["distribution"],1:5,1:2] );
# @variable(sub, yb[1:dt.N["plant"]+dt.N["distribution"],1:2], Bin);
# @variable(sub, ubij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:dt.Mij[i,j]], Bin);
# @variable(sub, ubjk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:dt.Mjk[j,k]], Bin);
# @variable(sub, ubkl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:dt.Mkl[k,l]], Bin);
# @constraints(sub, begin
#     [j=1:dt.N["plant"],p=1:5], sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j]) == sum(xjk[j,k,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k])
#     [k=1:dt.N["distribution"],p=1:5], sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k]) == sum(xkl[k,l,m,p] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])
# end)
# @constraints(sub, begin
#     [j=1:dt.N["plant"],p=1:5], sum(h[j,p,:]) == sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j])
#     [k=1:dt.N["distribution"],p=1:5], sum(h[k+dt.N["plant"],p,:] ) == sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k])
#     [l=1:dt.N["customer"],p=1:5], sum(xkl[k,l,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mkl[k,l]) >= dt.d[l][p]
# end )
# con6 = @constraint(sub,[i=1:dt.N["supplier"]], sum(xij[i,j,m,p] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) <= dt.N["cas"][i]);
# con7 = @constraint(sub,[j=1:dt.N["plant"], t=1:2], sum(h[j,1:5,t]) <=dt.N["cap"][j]*yb[j,t]);
# con8 = @constraint(sub,[k=1:J+dt.N["distribution"], t=1:2], sum(h[k,1:5,t]) >= dt.N["cad"][k]*yb[J+k,t]);
# con9 = @constraint(sub,[i=1:dt.N["supplier"], j=1:dt.N["plant"], m=1:dt.Mij[i,j]], sum(xij[i,j,m,p] for p=1:5) >= dt.Vij[i][j][m]*ubij[i,j,m] );
# con10 = @constraint(sub,[j=1:dt.N["plant"], k=1:dt.N["distribution"], m=1:dt.Mjk[j,k]], sum(xjk[j,k,m,p] for p=1:5) >= dt.Vjk[j][k][m]*ubjk[j,k,m]);
# con11 = @constraint(sub,[k=1:dt.N["distribution"], l=1:dt.N["customer"], m=1:dt.Mkl[k,l]], sum(xkl[k,l,m,p] for p=1:5) >= dt.Vkl[k][l][m]*ubkl[k,l,m]);
# @objective( sub, Min, w[1]*( exa + sum(dt.e[j][(p-1)*2+t]*h[j,p,t] for j=1:dt.N["plant"] + dt.N["distribution"] for p=1:5 for t=1:2)) +
#     exv + w[2]*( exb + sum(dt.q[j][(p-1)*2+t]*h[j,p,t] for j=1:dt.N["plant"]+dt.N["distribution"] for p=1:5 for t=1:2) +exr) );

######################### Mathematical model ###################################
################   original SCND  Benders Decomposition     ################################
w = [0.5,0.5]
mas = Model(CPLEX.Optimizer);
MOI.set(mas, MOI.RawParameter("CPX_PARAM_SCRIND"), false )
MOI.set(mas, MOI.RawParameter("CPX_PARAM_THREADS"),1  )
@variable(mas, yjt[1:dt.N["plant"],1:2], Bin  );
@variable(mas, ykt[1:dt.N["distribution"],1:2], Bin  );
@variable(mas, uij[1:dt.N["supplier"],1:dt.N["plant"], 1:dt.m] , Bin);
@variable(mas, ujk[1:dt.N["plant"],1:dt.N["distribution"],1:dt.m] , Bin);
@variable(mas, ukl[1:dt.N["distribution"],1:dt.N["customer"],1:dt.m] , Bin);
@variable(mas, θ)
exg = AffExpr(0);
for i=1:I
    idx = 1;
    for j=1:J
        for m=1:Mij[i,j]
            add_to_expression!(exg, dt.gij[i][idx]*uij[i,j,m]);
            idx+=1
        end
    end
end
@objective(mas, Min, w[1]*(sum(dt.c[j][t]*y[j,t] for j=1:1:dt.N["plant"]+dt.N["distribution"] for t=1:2) + exg) + θ);
#####################  modelling with matrices ############################
mutable struct CallModel
    lpfile::String; m::Int; n::Int; C::Array{}; B::Array{}; RHS::Array{}; signs::Array{}; vub::Array{}
    function CallModel(lpfile::String)
        lpmodel=buildlp([-1,0],[2 1],'<',1.5, CplexSolver(CPX_PARAM_SCRIND=0))
        MPB.loadproblem!(lpmodel,lpfile)
        Bmtx = MPB.getconstrmatrix(lpmodel);
        B = Bmtx[3:end,:]; C = Bmtx[1:2,:]
        # cut = find(i-> varub[i]==1 &&varub[i+1]!=1, 1:length(varub))[end]
        # vub = varub[1:cut]; B = Bmtx[3:end,1:cut]; C = Bmtx[1:2,1:cut]
        m,n=size(B)
        vub = MPB.getvarUB(lpmodel)
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
        new(lpfile,m,n,C,B,RHS,signs,vub)
    end
end
function print_iteration(k, args...)
    f(x) = Printf.@sprintf("%12.4e", x)
    println(lpad(k, 9), " ", join(f.(args), " "))
    return
end
# mt = SCNDModel("F:/scnd/Test1S2.lp")
mt = CallModel("/home/ak121396/Desktop/instances/SCND/test01S2.lp") #small_ins.lp")
bvar = findall(i->i==1,mt.vub); rvar = findall(i->i!=1,mt.vub);
c1b = [mt.C[1,i] for i in bvar];  #c2b is zeros
c1r = [mt.C[1,i] for i in rvar]; c2r = [mt.C[2,i] for i in rvar];
range1 = findall(i-> mt.vub[i]==1 && mt.vub[i+1]!=1, 1:length(mt.vub)-1)[1]
range2 = findall(i-> mt.vub[i]!=1 && mt.vub[i+1]==1, 1:length(mt.vub)-1)[1]
firstcon = findall(i->sum(mt.B[i,range1+1:range2])==0, 1:length(mt.RHS))
secondcon = findall(i->i ∉ firstcon, 1:length(mt.RHS))
mas0 = nothing;
for i in firstcon
    mas0 = [mas0; mt.B[i,:]']
end
mas1 = mas0[2:end,:];
mascon = hcat(mas1[:,1:range1],mas1[:,range2+1:end]) #master problem constraints
mRHS = [mt.RHS[i] for i in firstcon];
msigns = [mt.signs[i] for i in firstcon];
sub0 = nothing;
for i in secondcon
    sub0 = [sub0; mt.B[i,:]']
end
sub1 = sub0[2:end,:];
# subcon = [mt.B[1:firstcon[1]-1,range1+1:range2]; mt.B[firstcon[end]+1:end,range1+1:range2]]
sRHS = [mt.RHS[i] for i in secondcon]
subsigns = [mt.signs[i] for i in secondcon]
dim_x = length(c1b); dim_y = length(c1r);
M = sum(mt.RHS[i] for i in findall(i->i<0, mt.RHS));

A_1 = hcat(sub1[:,1:range1],sub1[:,range2+1:end])
A_2 = sub1[:,range1+1:range2]
# b = [mRHS; sRHS]
inequal = findall(i->subsigns[i]!="s",1:length(subsigns))
equal = findall(i->subsigns[i]=="s",1:length(subsigns))
inRHS = [sRHS[j] for j in inequal]
eqRHS = [sRHS[j] for j in equal]
A1,A1e,A2,A2e = zeros(1,dim_x),zeros(1,dim_x),zeros(1,dim_y),zeros(1,dim_y);
for j=1:size(A_1,1)
    if j in inequal
        A1 = [A1; A_1[j,:]']; A2 = [A2; A_2[j,:]']
    else
        A1e = [A1e; A_1[j,:]']; A2e = [A2e; A_2[j,:]']
    end
end
inA1 = A1[2:end,:]; eqA1 = A1e[2:end,:];
inA2 = A2[2:end,:]; eqA2 = A2e[2:end,:];

#####################   Benders Decomposition  #######################
#Master problem
model = Model(CPLEX.Optimizer)
@variable(model, x[1:dim_x], Bin);
# @variable(model, 0<=x[1:dim_x]<=1) LP relaxation
@variable(model, θ >= M)
for k=1:length(msigns)
    if msigns[k] == "l"
        @constraint(model, dot(mascon[k,:]',x) >= mRHS[k])
    elseif msigns[k] == "u"
        @constraint(model, dot(mascon[k,:]',x) <= mRHS[k])
    else
        @constraint(model, dot(mascon[k,:]',x) == mRHS[k])
    end
end
@objective(model, Min, c1b' * x + θ)

1
function solve_subproblem(x)
    model = Model(CPLEX.Optimizer)
    @variable(model, y[1:dim_y] >= 0)
    con = @constraint(model, inA2 * y .<= inRHS - inA1 * x)
    econ = @constraint(model, eqA2 * y .== eqRHS - eqA1 *x)
    @objective(model, Min, c1r' * y)# + c2r' * y*w[2] )
    optimize!(model)
    st = termination_status(model)
    if st == MOI.OPTIMAL
        return (res = :OptimalityCut, obj = objective_value(model), y = value.(y), π = dual.(con), σ = dual.(econ) )
    elseif st == MOI.INFEASIBLE
        return (res = :FeasibilityCut, obj = objective_value(model), y = value.(y), π = dual.(con), σ = dual.(econ) )
    else
        error("DualSubProblem error: status $st")
    end
end
