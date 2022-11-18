using DataFrames,DelimitedFiles,JuMP,LinearAlgebra,CPLEX,MathProgBase,MathOptInterface,SparseArrays#,CPUTime
mutable struct Data
    file::String; N::Dict{}; d::Array{}; c::Array{}; a::Array{}; e::Array{}; gij::SparseVector{}; gjk::SparseVector{}; gkl::SparseVector{};
    vij::Array{}; vjk::Array{}; vkl::Array{}; Vij::SparseVector{}; Vjk::SparseVector{}; Vkl::SparseVector{}; Mij::Array{}; Mjk::Array{}; Mkl::Array{};
    b::Array{}; q::Array{}; rij::Array{}; rjk::Array{}; rkl::Array{}; upl::Int; udc::Int
    function Data(file)
        dt = readdlm(file);
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
                    append!(W,tmp)
                    # push!(W,tmp);
                end
                # tmp = [filter(x->x!="", dt[x,:]) for x in id1+1:id1+(id2-id1-1)]
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
# @show file = ARGS[1];
# file = "E:/scnd/Test1S2"
file = "/home/ak121396/Desktop/instances/SCND/test01S2"
# file = "/home/k2g00/k2g3475/scnd/instances/test01S1"
dt = Data(file);
############################  Benders Decomposition  ##################################
struct build_master
    y::Vector{VariableRef}
    uij::Vector{VariableRef}
    ujk::Vector{VariableRef}
    ukl::Vector{VariableRef}
    m::Model
end
function build_master()

    mas = Model(CPLEX.Optimizer); set_silent(mas)
    MOI.set(mas, MOI.NumberOfThreads(), 1); set_silent(mas)
    @variable(mas, y[1:(dt.N["plant"]+dt.N["distribution"])*2], Bin);
    @variable(mas, uij[1:sum(dt.Mij)], Bin);
    @variable(mas, ujk[1:sum(dt.Mjk)], Bin);
    @variable(mas, ukl[1:sum(dt.Mkl)], Bin);
    @constraint(mas,[j=1:dt.N["plant"]+dt.N["distribution"]], sum(y[2*(j-1)+1:2*(j-1)+2]) <= 1);
    @constraints(mas,begin
            sum(uij[1:dt.Mij[1,1]]) <= 1
            [j=2:dt.N["plant"]], sum(uij[sum(dt.Mij[1,1:j-1])+1:sum(dt.Mij[1,1:j-1])+dt.Mij[1,j]]) <= 1
            [i=2:dt.N["supplier"],j=2:dt.N["plant"]],  sum(uij[sum(dt.Mij[1:i-1,:])+sum(dt.Mij[i,1:j-1])+1:sum(dt.Mij[1:i-1,:])+sum(dt.Mij[i,1:j-1])+dt.Mij[i,j]])<= 1
            sum(ujk[1:dt.Mjk[1,1]]) <= 1
            [k=2:dt.N["distribution"]], sum(ujk[sum(dt.Mjk[1,1:k-1])+1:sum(dt.Mjk[1,1:k-1])+dt.Mjk[1,k]]) <= 1
            [j=2:dt.N["plant"],k=2:dt.N["distribution"]],  sum(ujk[sum(dt.Mjk[1:j-1,:])+sum(dt.Mjk[j,1:j-1])+1:sum(dt.Mjk[1:j-1,:])+sum(dt.Mjk[j,1:j-1])+dt.Mjk[j,k]]) <= 1
            sum(ukl[1:dt.Mkl[1,1]]) <= 1
            [l=2:dt.N["customer"]], sum(ukl[sum(dt.Mkl[1,1:l-1])+1:sum(dt.Mkl[1,1:l-1])+dt.Mkl[1,l]]) <= 1
            [k=2:dt.N["distribution"],l=2:dt.N["customer"]],  sum(ukl[sum(dt.Mkl[1:k-1,:])+sum(dt.Mkl[k,1:l-1])+1:sum(dt.Mkl[1:k-1,:])+sum(dt.Mkl[k,1:l-1])+dt.Mkl[k,l]])<= 1
    end);
    @constraint(mas, sum(y[1:dt.N["plant"]]) <= dt.upl);
    @constraint(mas, sum(y[dt.N["plant"]+1:end]) <= dt.udc);
    return (y=y,uij=uij,ujk=ujk,ukl=ukl,m=mas)
end

struct DualSubP
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
    m::Model
end

function DualSubP(sub::Model)
    set_optimizer_attribute(sub, "CPX_PARAM_REDUCE", 0); # turning presolve off
    # set_optimizer_attribute(sub, "CPX_PARAM_PREIND", 0); # turning preprocessing entirely off
    set_silent(sub) #MOI.set(sub, MOI.RawParameter("CPX_PARAM_SCRIND"), false )
    MOI.set(sub, MOI.NumberOfThreads(), 1) # MOI.set(sub, MOI.RawParameter("CPX_PARAM_THREADS"),1  )
    @variables(sub, begin
        α1[1:dt.N["plant"]*5]
        α2[1:dt.N["distribution"]*5]
        α3[1:dt.N["plant"]*5]
        α4[1:dt.N["distribution"]*5]
    end);
    @variable(sub, α5[1:dt.N["customer"]*5] >= 0);
    @variable(sub, α6[1:dt.N["supplier"]] >= 0);
    @variable(sub, α7[1:dt.N["plant"]*2] >= 0);
    @variable(sub, α8[1:dt.N["distribution"]*2] >= 0);
    @variable(sub, α9[1:sum(dt.Mij)] >= 0);
    @variable(sub, α10[1:sum(dt.Mjk)] >= 0);
    @variable(sub, α11[1:sum(dt.Mkl)] >= 0);
    ###################### constraint 39  ####################################
    @constraint(sub, [i=1,j=1,m=1:dt.Mij[i,j],p=1:5], α1[p]-α3[p]-α6[i]+α9[m] <= dt.a[i,p]+dt.vij[m]+dt.b[i,p]+dt.rij[m]);
    @constraint(sub, [i=1,j=2:dt.N["plant"],m=1:dt.Mij[i,j],p=1:5],  α1[5*(j-1)+p]-α3[5*(j-1)+p]-α6[i]+α9[sum(dt.Mij[i,1:j-1])+m] <= dt.a[i,p]+dt.vij[sum(dt.Mij[i,1:j-1])+m]+dt.b[i,p]+dt.rij[sum(dt.Mij[i,1:j-1])+m]);
    @constraint(sub, [i=2:dt.N["supplier"],j=1,m=1:dt.Mij[i,j],p=1:5],  α1[p]-α3[p]-α6[i]+α9[sum(dt.Mij[1:i-1,:])+sum(dt.Mij[i,1])+m] <= dt.a[i,p]+dt.vij[sum(dt.Mij[1:i-1,:])+m]+dt.b[i,p]+dt.rij[sum(dt.Mij[1:i-1,:])+m]);
    @constraint(sub, [i=2:dt.N["supplier"],j=2:dt.N["plant"],m=1:dt.Mij[i,j],p=1:5], α1[5*(j-1)+p]-α3[5*(j-1)+p]-α6[i]+α9[sum(dt.Mij[1:i-1,:])+sum(dt.Mij[i,1:j-1])+m] <= dt.a[i,p]+dt.vij[sum(dt.Mij[1:i-1,:])+sum(dt.Mij[i,1:j-1])+m]+dt.b[i,p]+dt.rij[sum(dt.Mij[1:i-1,:])+sum(dt.Mij[i,1:j-1])+m]);
    ###################### constraint 40  ####################################
    @constraint(sub, [j=1,k=1,m=1:dt.Mjk[j,k],p=1:5], -α1[p]-α2[p]-α4[p]+α10[m] <= dt.vjk[m]+dt.rjk[m]);
    @constraint(sub, [j=1,k=2:dt.N["distribution"],m=1:dt.Mjk[j,k],p=1:5], -α1[p]-α2[5*(k-1)+p]-α4[5*(k-1)+p]+α10[m] <= dt.vjk[sum(dt.Mjk[j,1:k-1])+m]+dt.rjk[sum(dt.Mjk[j,1:k-1])+m]);
    @constraint(sub, [j=2:dt.N["plant"],k=1,m=1:dt.Mjk[j,k],p=1:5], -α1[5*(j-1)+p]-α2[p]-α4[p]+α10[sum(dt.Mjk[1:j-1,:])+sum(dt.Mjk[j,1:j-1])+m] <= dt.vjk[sum(dt.Mjk[1:j-1,:])+m]+dt.rjk[m]);
    @constraint(sub, [j=2:dt.N["plant"],k=2:dt.N["distribution"],m=1:dt.Mjk[j,k],p=1:5], -α1[5*(j-1)+p]-α2[5*(k-1)+p]-α4[5*(k-1)+p]+α10[sum(dt.Mjk[1:j-1,:])+sum(dt.Mjk[j,1:j-1])+m] <= dt.vjk[sum(dt.Mjk[1:j-1,:])+sum(dt.Mjk[j,1:k-1])+m]+dt.rjk[sum(dt.Mjk[j,1:k-1])+m]);
    ###################### constraint 41  ####################################
    @constraint(sub, [k=1,l=1,m=1:dt.Mkl[k,l],p=1:5],-α2[p]+α5[p]+α11[m] <= dt.vkl[m]+dt.rkl[m]);
    @constraint(sub, [k=1,l=2:dt.N["customer"],m=1:dt.Mkl[k,l],p=1:5], -α2[p]+α5[5*(l-1)+p]+α11[sum(dt.Mkl[k,1:l-1])+m] <= dt.vkl[sum(dt.Mkl[k,1:l-1])+m]+dt.rkl[sum(dt.Mkl[k,1:l-1])+m]);
    @constraint(sub, [k=2:dt.N["distribution"],l=1,m=1:dt.Mkl[k,l],p=1:5], -α2[5*(k-1)+p]+α5[p]+α11[sum(dt.Mkl[1:k-1,:])+m] <= dt.vkl[sum(dt.Mkl[1:k-1,:])+m]+dt.rkl[m]);
    @constraint(sub, [k=2:dt.N["distribution"],l=2:dt.N["customer"],m=1:dt.Mkl[k,l],p=1:5],-α2[5*(k-1)+p]+α5[5*(l-1)+p]+α11[sum(dt.Mkl[1:k-1,:])+sum(dt.Mkl[k,1:l-1])+m] <= dt.vkl[sum(dt.Mkl[1:k-1,:])+sum(dt.Mkl[k,1:l-1])+m]+dt.rkl[sum(dt.Mkl[k,1:l-1])+m]);
    ###################### constraint 42  ####################################
    @constraint(sub, [t=1:2,p=1:5], α3[p]-α7[t] <= dt.e[5*(t-1)+p]);
    @constraint(sub, [j=2:dt.N["plant"],t=1:2,p=1:5], α3[5*(j-1)+p]-α7[2*(j-1)+t] <= dt.e[5*2*(j-1)+((t-1)*5)+p]);
    ###################### constraint 43  ####################################
    @constraint(sub, [t=1:2,p=1:5], α4[p]-α8[t] <= dt.e[(dt.N["plant"]*2*5)+5*(t-1)+p]);
    @constraint(sub, [k=2:dt.N["distribution"],t=1:2,p=1:5], α4[5*(k-1)+p]-α8[2*(k-1)+t] <= dt.e[(dt.N["plant"]*2*5)+5*2*(k-1)+((t-1)*5)+p]);
    return DualSubP(α1,α2,α3,α4,α5,α6,α7,α8,α9,α10,α11,sub)
end

function solve_dsp(dsp::DualSubP,yb,ubij,ubjk,ubkl)
    @objective(dsp.m, Max, sum(dsp.α5[5*(l-1)+p]*dt.d[l,p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*dsp.α6[i] for i=1:dt.N["supplier"])-
        sum(dt.N["cap"][j]*yb[2*(j-1)+t]*dsp.α7[2*(j-1)+t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*yb[2*dt.N["plant"]+k]*dsp.α8[2*(k-1)+t] for k=1:dt.N["distribution"] for t=1:2)+
        sum(dt.Vij[i]*ubij[i]*dsp.α9[i] for i in findnz(dt.Vij)[1])+sum(dt.Vjk[i]*ubjk[i]*dsp.α10[i] for i in findnz(dt.Vjk)[1])+sum(dt.Vkl[i]*ubkl[i]*dsp.α11[i] for i in findnz(dt.Vkl)[1])
        )
    optimize!(dsp.m)
    @show st = termination_status(dsp.m)

    if st == MOI.OPTIMAL
      return (res = :OptimalityCut, obj = objective_value(dsp.m), α5 = value.(dsp.α5), α6 = value.(dsp.α6), α7 = value.(dsp.α7), α8= value.(dsp.α8), α9= value.(dsp.α9), α10= value.(dsp.α10), α11= value.(dsp.α11) )
    elseif st == MOI.DUAL_INFEASIBLE
        # lens = [length(dsp.α5),length(dsp.α6),length(dsp.α7),length(dsp.α8),length(dsp.α9),length(dsp.α10),length(dsp.α11),length(dsp.α1),length(dsp.α2),length(dsp.α3),length(dsp.α4)]
        # lens = [length(dsp.α6),length(dsp.α9),length(dsp.α1),length(dsp.α3),length(dsp.α10),length(dsp.α2),length(dsp.α4),length(dsp.α5),length(dsp.α11),length(dsp.α7),length(dsp.α8)]
        # rays = Vector{Float64}(undef,sum(dsp.lens))
        # CPXgetray(backend(dsp.m).env, backend(dsp.m).lp, rays) # vr = value.(rays)
        # pvci = [sum(lens[1:l]) for l=1:length(lens)]
        # insert!(pvci,1,0)
        # slices = [rays[pvci[l]+1:pvci[l+1]] for l=1:length(pvci)-1]
        return ( res = :FeasibilityCut, α5 = value.(dsp.α5), α6 = value.(dsp.α6), α7 = value.(dsp.α7), α8= value.(dsp.α8), α9= value.(dsp.α9), α10= value.(dsp.α10), α11= value.(dsp.α11) )
        # return ( res = :FeasibilityCut, α5=slices[8], α6=slices[1], α7=slices[10], α8=slices[11], α9=slices[2], α10=slices[5], α11=slices[9] )
    else
      error("DualSubProblem error: status $st")
    end
end

function benders_decomposition(w,mas::Model, y::Vector{VariableRef},uij::Vector{VariableRef},ujk::Vector{VariableRef},ukl::Vector{VariableRef})
    # uij::JuMP.Containers.SparseAxisArray{VariableRef},ujk::JuMP.Containers.SparseAxisArray{VariableRef},ukl::JuMP.Containers.SparseAxisArray{VariableRef})
    sub = direct_model(CPLEX.Optimizer()); dsp = DualSubP(sub)

    @variable(mas, θ>= -1000);
    @objective(mas, Min, w[1]*(sum(dt.c.*y)+sum(dt.gij[i]*uij[i] for i in findnz(dt.gij)[1]) + sum(dt.gjk[i]*ujk[i] for i in findnz(dt.gjk)[1]) + sum(dt.gkl[i]*ukl[i] for i in findnz(dt.gkl)[1])+ θ));
    optimize!(mas);
    st = termination_status(mas)

    while (st == MOI.INFEASIBLE) || (st == MOI.OPTIMAL)
        optimize!(mas);
        @show  st = termination_status(mas)
        θ1 = value(θ); yb = value.(y); ubij = value.(uij); ubjk = value.(ujk); ubkl = value.(ukl)
        # @show iszero(yb),iszero(ubij),iszero(ubjk),iszero(ubkl)
        subp = solve_dsp(dsp,yb,ubij,ubjk,ubkl)

        if subp.res == :OptimalityCut
            @info "Optimality cut found"
            if θ1 ≥ subp.obj
                break
            else
                cut = @constraint( mas, θ ≥ sum(subp.α5[5*(l-1)+p]*dt.d[l,p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
                sum(dt.N["cap"][j]*y[2*(j-1)+t]*subp.α7[2*(j-1)+t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*y[2*dt.N["plant"]+k]*subp.α8[2*(k-1)+t] for k=1:dt.N["distribution"] for t=1:2)+
                sum(dt.Vij[i]*uij[i]*subp.α9[i] for i in findnz(dt.Vij)[1])+sum(dt.Vjk[i]*ujk[i]*subp.α10[i] for i in findnz(dt.Vjk)[1])+sum(dt.Vkl[i]*ukl[i]*subp.α11[i] for i in findnz(dt.Vkl)[1]))
            end
        else
            @info "Feasibility cut found"
            cut = @constraint( mas, 0 ≥ sum(subp.α5[5*(l-1)+p]*dt.d[l,p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
            sum(dt.N["cap"][j]*y[2*(j-1)+t]*subp.α7[2*(j-1)+t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*y[2*dt.N["plant"]+k]*subp.α8[2*(k-1)+t] for k=1:dt.N["distribution"] for t=1:2)+
            sum(dt.Vij[i]*uij[i]*subp.α9[i] for i in findnz(dt.Vij)[1])+sum(dt.Vjk[i]*ujk[i]*subp.α10[i] for i in findnz(dt.Vjk)[1])+sum(dt.Vkl[i]*ukl[i]*subp.α11[i] for i in findnz(dt.Vkl)[1]))
        end
        @info "Adding the cut $(cut)"
        # println(mas)
    end
end
# w = [1,1];
mas = build_master();
benders_decomposition(w,mas.m,mas.y,mas.uij,mas.ujk,mas.ukl)
