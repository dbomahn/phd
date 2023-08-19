using DataFrames,DelimitedFiles,JuMP,LinearAlgebra,CPLEX,MathOptInterface,CPUTime # SparseArrays
mutable struct Data
    file::String; N::Dict{}; d::Array{}; c::Array{}; a::Array{}; e::Array{}; gij::Array{}; gjk::Array{}; gkl::Array{};
    vij::Array{}; vjk::Array{}; vkl::Array{}; Vij::Array{}; Vjk::Array{}; Mij::Array{}; Mjk::Array{}; Mkl::Array{};
    b::Array{}; q::Array{}; rij::Array{}; rjk::Array{}; rkl::Array{}; upl::Int; udc::Int; bigM::Int
    function Data(file)
        dt = readdlm(file);
        notafile = readdlm("/home/ak121396/Desktop/instances/scnd/Notations.txt", '=');
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
        gij = (N["fixedcostModesp"]); gjk = (N["fixedcostModepd"]); gkl = (N["fixedcostModedc"]);
        vij = N["tcp"]; vjk = N["tcd"]; vkl = N["tcc"];
        Vij = (N["LcapacityModesp"]); Vjk = (N["LcapacityModepd"]); #Vkl = (N["LcapacityModedc"]);
        Mij = transpose(reshape(N["ModeIJ"], (N["plant"],N["supplier"])));
        Mjk = transpose(reshape(N["ModeJK"], (N["distribution"],N["plant"])));
        Mkl = transpose(reshape(N["ModeKL"], (N["customer"],N["distribution"])));
        b = reshape(N["ves"],5,N["supplier"])';  q = append!(N["vep"],N["ved"]);
        rij = N["cep"]; rjk = N["ced"]; rkl = N["cec"];
        upl = N["upperpants"]; udc = N["upperdistribution"];bigM = sum(N["demand"])

        new(file,N,d,c,a,e,gij,gjk,gkl,vij,vjk,vkl,Vij,Vjk,Mij,Mjk,Mkl,b,q,rij,rjk,rkl,upl,udc,bigM)
    end
end
# @show file = ARGS[1];
# file = "E:/scnd/Test1S2"
file = "/home/ak121396/Desktop/instances/scnd/test01S2"
# file = "/home/k2g00/k2g3475/scnd/instances/test01S1"
dt = Data(file);
############################  Benders Decomposition  ##################################
struct Master
    y::Vector{VariableRef}
    uij::Vector{VariableRef}
    ujk::Vector{VariableRef}
    ukl::Vector{VariableRef}
    θ::VariableRef
    m::Model
end
function build_master(w)
    mas = direct_model(CPLEX.Optimizer()); set_silent(mas)
    MOI.set(mas, MOI.NumberOfThreads(), 1); 
    @variable(mas, y[1:(dt.N["plant"]+dt.N["distribution"])*2], Bin);
    @variable(mas, uij[1:sum(dt.Mij)], Bin);
    @variable(mas, ujk[1:sum(dt.Mjk)], Bin);
    @variable(mas, ukl[1:sum(dt.Mkl)], Bin);
    @variable(mas, θ>= -10000);

    @constraint(mas,[j=1:dt.N["plant"]+dt.N["distribution"]], sum(y[2*(j-1)+1:2*(j-1)+2]) <= 1);
    @constraints(mas,begin
    sum(uij[1:dt.Mij[1,1]]) <= 1
    sum(uij[sum(dt.Mij[1,:])+dt.Mij[2,1]]) <= 1
        [j=2:dt.N["plant"]], sum(uij[sum(dt.Mij[1,1:j-1])+1:sum(dt.Mij[1,1:j-1])+dt.Mij[1,j]]) <= 1
        [i=2:dt.N["supplier"],j=2:dt.N["plant"]],  sum(uij[sum(dt.Mij[1:i-1,:])+sum(dt.Mij[i,1:j-1])+1:sum(dt.Mij[1:i-1,:])+sum(dt.Mij[i,1:j-1])+dt.Mij[i,j]])<= 1
        sum(ujk[1:dt.Mjk[1,1]]) <= 1
        sum(ujk[sum(dt.Mjk[1,:])+dt.Mjk[2,1]]) <= 1
        [k=2:dt.N["distribution"]], sum(ujk[sum(dt.Mjk[1,1:k-1])+1:sum(dt.Mjk[1,1:k-1])+dt.Mjk[1,k]]) <= 1
        [j=2:dt.N["plant"],k=2:dt.N["distribution"]],  sum(ujk[sum(dt.Mjk[1:j-1,:])+sum(dt.Mjk[j,1:j-1])+1:sum(dt.Mjk[1:j-1,:])+sum(dt.Mjk[j,1:j-1])+dt.Mjk[j,k]]) <= 1
        sum(ukl[1:dt.Mkl[1,1]]) <= 1 #[1,1]
        sum(ukl[sum(dt.Mkl[1,:])+dt.Mkl[2,1]]) <= 1 #[2,1]
        [l=2:dt.N["customer"]], sum(ukl[sum(dt.Mkl[1,1:l-1])+1:sum(dt.Mkl[1,1:l-1])+dt.Mkl[1,l]]) <= 1
        [k=2:dt.N["distribution"],l=2:dt.N["customer"]],  sum(ukl[sum(dt.Mkl[1:k-1,:])+sum(dt.Mkl[k,1:l-1])+1:sum(dt.Mkl[1:k-1,:])+sum(dt.Mkl[k,1:l-1])+dt.Mkl[k,l]])<= 1
    end);
    @constraint(mas, sum(y[1:dt.N["plant"]]*2) <= dt.upl);
    @constraint(mas, sum(y[dt.N["plant"]*2+1:end]) <= dt.udc);
    @objective(mas, Min, w*(sum(dt.c.*y)+dot(dt.gij,uij) + dot(dt.gjk,ujk)+dot(dt.gkl,ukl))+ θ);
    
    return Master(y,uij,ujk,ukl,θ,mas)
end

struct DualSubP
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
    α12::Vector{VariableRef}
    α13::Vector{VariableRef}
    α14::Vector{VariableRef}
    m::Model
end

function DualSubP(w,ratio)
    sub = direct_model(CPLEX.Optimizer());
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
    @variable(sub, α12[1:sum(dt.Mij)] >= 0);
    @variable(sub, α13[1:sum(dt.Mjk)] >= 0);
    @variable(sub, α14[1:sum(dt.Mkl)] >= 0);

    ###################### constraint 39  ####################################
    @constraint(sub, [i=1,j=1,m=1:dt.Mij[i,j],p=1:5], α1[p]-α3[p]-α6[i]+α9[m]-α12[m] <= w*(dt.a[i,p]+dt.vij[m])+ratio*(1-w)*(dt.b[i,p]+dt.rij[m]));
    @constraint(sub, [i=1,j=2:dt.N["plant"],m=1:dt.Mij[i,j],p=1:5],  α1[5*(j-1)+p]-α3[5*(j-1)+p]-α6[i]+α9[sum(dt.Mij[i,1:j-1])+m]-α12[sum(dt.Mij[i,1:j-1])+m] <= w*(dt.a[i,p]+dt.vij[sum(dt.Mij[i,1:j-1])+m])+ratio*(1-w)*(dt.b[i,p]+dt.rij[sum(dt.Mij[i,1:j-1])+m]));
    @constraint(sub, [i=2:dt.N["supplier"],j=1,m=1:dt.Mij[i,j],p=1:5],  α1[p]-α3[p]-α6[i]+α9[sum(dt.Mij[1:i-1,:])+sum(dt.Mij[i,1])+m]-α12[sum(dt.Mij[1:i-1,:])+sum(dt.Mij[i,1])+m] <= w*(dt.a[i,p]+dt.vij[sum(dt.Mij[1:i-1,:])+m])+ratio*(1-w)*(dt.b[i,p]+dt.rij[sum(dt.Mij[1:i-1,:])+m]));
    @constraint(sub, [i=2:dt.N["supplier"],j=2:dt.N["plant"],m=1:dt.Mij[i,j],p=1:5], α1[5*(j-1)+p]-α3[5*(j-1)+p]-α6[i]+α9[sum(dt.Mij[1:i-1,:])+sum(dt.Mij[i,1:j-1])+m]-α12[sum(dt.Mij[1:i-1,:])+sum(dt.Mij[i,1:j-1])+m] <= w*(dt.a[i,p]+dt.vij[sum(dt.Mij[1:i-1,:])+sum(dt.Mij[i,1:j-1])+m])+ratio*(1-w)*(dt.b[i,p]+dt.rij[sum(dt.Mij[1:i-1,:])+sum(dt.Mij[i,1:j-1])+m]));
    ###################### constraint 40  ####################################
    @constraint(sub, [j=1,k=1,m=1:dt.Mjk[j,k],p=1:5], -α1[p]-α2[p]-α4[p]+α10[m]-α13[m] <= w*(dt.vjk[m])+ratio*(1-w)*(dt.rjk[m]));
    @constraint(sub, [j=1,k=2:dt.N["distribution"],m=1:dt.Mjk[j,k],p=1:5], -α1[p]-α2[5*(k-1)+p]-α4[5*(k-1)+p]+α10[m]-α13[m] <= w*dt.vjk[sum(dt.Mjk[j,1:k-1])+m]+ratio*(1-w)*(dt.rjk[sum(dt.Mjk[j,1:k-1])+m]));
    @constraint(sub, [j=2:dt.N["plant"],k=1,m=1:dt.Mjk[j,k],p=1:5], -α1[5*(j-1)+p]-α2[p]-α4[p]+α10[sum(dt.Mjk[1:j-1,:])+sum(dt.Mjk[j,1:j-1])+m]-α13[sum(dt.Mjk[1:j-1,:])+sum(dt.Mjk[j,1:j-1])+m] <= w*dt.vjk[sum(dt.Mjk[1:j-1,:])+m]+ratio*(1-w)*dt.rjk[m]);
    @constraint(sub, [j=2:dt.N["plant"],k=2:dt.N["distribution"],m=1:dt.Mjk[j,k],p=1:5], -α1[5*(j-1)+p]-α2[5*(k-1)+p]-α4[5*(k-1)+p]+α10[sum(dt.Mjk[1:j-1,:])+sum(dt.Mjk[j,1:j-1])+m]-α13[sum(dt.Mjk[1:j-1,:])+sum(dt.Mjk[j,1:j-1])+m] <= w*dt.vjk[sum(dt.Mjk[1:j-1,:])+sum(dt.Mjk[j,1:k-1])+m]+ratio*(1-w)*dt.rjk[sum(dt.Mjk[j,1:k-1])+m]);
    ###################### constraint 41  ####################################
    @constraint(sub, [k=1,l=1,m=1:dt.Mkl[k,l],p=1:5],-α2[p]+α5[p]-α14[m] <= w*dt.vkl[m]+ratio*(1-w)*dt.rkl[m]);
    @constraint(sub, [k=1,l=2:dt.N["customer"],m=1:dt.Mkl[k,l],p=1:5], -α2[p]+α5[5*(l-1)+p]-α14[sum(dt.Mkl[k,1:l-1])+m] <= w*dt.vkl[sum(dt.Mkl[k,1:l-1])+m]+ratio*(1-w)*dt.rkl[sum(dt.Mkl[k,1:l-1])+m]);
    @constraint(sub, [k=2:dt.N["distribution"],l=1,m=1:dt.Mkl[k,l],p=1:5], -α2[5*(k-1)+p]+α5[p]-α14[sum(dt.Mkl[1:k-1,:])+m] <= w*dt.vkl[sum(dt.Mkl[1:k-1,:])+m]+ratio*(1-w)*dt.rkl[m]);
    @constraint(sub, [k=2:dt.N["distribution"],l=2:dt.N["customer"],m=1:dt.Mkl[k,l],p=1:5],-α2[5*(k-1)+p]+α5[5*(l-1)+p]-α14[sum(dt.Mkl[1:k-1,:])+sum(dt.Mkl[k,1:l-1])+m] <= w*dt.vkl[sum(dt.Mkl[1:k-1,:])+sum(dt.Mkl[k,1:l-1])+m]+ratio*(1-w)*dt.rkl[sum(dt.Mkl[k,1:l-1])+m]);
    ###################### constraint 42  ####################################
    @constraint(sub, [t=1:2,p=1:5], α3[p]-α7[t] <= w*dt.e[5*(t-1)+p]+ ratio*(1-w)*dt.q[5*(t-1)+p])
    @constraint(sub, [j=2:dt.N["plant"],t=1:2,p=1:5], α3[5*(j-1)+p]-α7[2*(j-1)+t] <= w*dt.e[5*2*(j-1)+((t-1)*5)+p]+ratio*(1-w)*dt.q[5*2*(j-1)+((t-1)*5)]);
    ###################### constraint 43  ####################################
    @constraint(sub, [t=1:2,p=1:5], α4[p]-α8[t] <= w*dt.e[(dt.N["plant"]*2*5)+5*(t-1)+p] + ratio*(1-w)*dt.q[(dt.N["plant"]*2*5)+5*(t-1)+p]);
    @constraint(sub, [k=2:dt.N["distribution"],t=1:2,p=1:5], α4[5*(k-1)+p]-α8[2*(k-1)+t] <= w*dt.e[(dt.N["plant"]*2*5)+5*2*(k-1)+((t-1)*5)+p] + ratio*(1-w)*dt.q[(dt.N["plant"]*2*5)+5*2*(k-1)+((t-1)*5)+p]);
    return DualSubP(α1,α2,α3,α4,α5,α6,α7,α8,α9,α10,α12,α13,α14,sub)
end

function solve_dsp(dsp::DualSubP,yb,ubij,ubjk,ubkl)
    @objective(dsp.m, Max, sum(dsp.α5[5*(l-1)+p]*dt.d[l,p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*dsp.α6[i] for i=1:dt.N["supplier"])-
        sum(dt.N["cap"][j]*yb[2*(j-1)+t]*dsp.α7[2*(j-1)+t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*yb[2*dt.N["plant"]+k]*dsp.α8[2*(k-1)+t] for k=1:dt.N["distribution"] for t=1:2)+
        dot(dt.Vij.*ubij,dsp.α9)+dot(dt.Vjk.*ubjk,dsp.α10) -
        dot(dt.bigM.*ubij,dsp.α12)-dot(dt.bigM.*ubjk,dsp.α13) -dot(dt.bigM.*ubkl,dsp.α14)
        )
    subt = @CPUelapsed optimize!(dsp.m); push!(subtime,subt)
    st = termination_status(dsp.m)

    if st == MOI.OPTIMAL
        return (res = :OptimalityCut, obj = objective_value(dsp.m), α5 = value.(dsp.α5), α6 = value.(dsp.α6), α7 = value.(dsp.α7), α8= value.(dsp.α8), α9= value.(dsp.α9), α10= value.(dsp.α10), α12= value.(dsp.α12),α13= value.(dsp.α13),α14= value.(dsp.α14) )
    elseif st == MOI.DUAL_INFEASIBLE
        return ( res = :FeasibilityCut, α5 = value.(dsp.α5), α6 = value.(dsp.α6), α7 = value.(dsp.α7), α8= value.(dsp.α8), α9= value.(dsp.α9), α10= value.(dsp.α10), α12= value.(dsp.α12),α13= value.(dsp.α13),α14= value.(dsp.α14) )
    else
        error("DualSubProblem error: status $st")
    end
end
      
function benders_decomposition(cb_data)
    yb = callback_value.(cb_data, mas.y);    ubij = callback_value.(cb_data, mas.uij);    ubjk = callback_value.(cb_data, mas.ujk)
    ubkl = callback_value.(cb_data, mas.ukl);    θb = callback_value(cb_data, mas.θ)
    subp = solve_dsp(dsp,yb,ubij,ubjk,ubkl);    
    if subp.res == :OptimalityCut
        if round(θb) ≥ round(subp.obj)#; digits=4)
            return
        else
            c1 = copy(rhs1); c2 = copy(rhs2); 
            CPLEX.CPXftran(backend(dsp.m).env, backend(dsp.m).lp, c1)
            CPLEX.CPXftran(backend(dsp.m).env, backend(dsp.m).lp, c2)
            push!(pi1, copy(c1)); push!(pi2, copy(c2)); 

            cut = @build_constraint( mas.θ ≥ sum(subp.α5[5*(l-1)+p]*dt.d[l,p] for l=1:dt.N["customer"] for p=1:5)-
                sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
                sum(dt.N["cap"][j]*mas.y[2*(j-1)+t]*subp.α7[2*(j-1)+t] for j=1:dt.N["plant"] for t=1:2)-
                sum(dt.N["cad"][k]*mas.y[2*dt.N["plant"]+k]*subp.α8[2*(k-1)+t] for k=1:dt.N["distribution"] for t=1:2)+
                dot(dt.Vij.*mas.uij,subp.α9)+dot(dt.Vjk.*mas.ujk,subp.α10) -
                dot(dt.bigM.*mas.uij,subp.α12)-dot(dt.bigM.*mas.ujk,subp.α13) -dot(dt.bigM.*mas.ukl,subp.α14)
                )
            MOI.submit(mas.m, MOI.LazyConstraint(cb_data), cut)
        end
    else
        cut = @build_constraint( 0 ≥ sum(subp.α5[5*(l-1)+p]*dt.d[l,p] for l=1:dt.N["customer"] for p=1:5)-
            sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
            sum(dt.N["cap"][j]*mas.y[2*(j-1)+t]*subp.α7[2*(j-1)+t] for j=1:dt.N["plant"] for t=1:2)-
            sum(dt.N["cad"][k]*mas.y[2*dt.N["plant"]+k]*subp.α8[2*(k-1)+t] for k=1:dt.N["distribution"] for t=1:2)+
            dot(dt.Vij.*mas.uij,subp.α9)+dot(dt.Vjk.*mas.ujk,subp.α10) -
            dot(dt.bigM.*mas.uij,subp.α12)-dot(dt.bigM.*mas.ujk,subp.α13) -dot(dt.bigM.*mas.ukl,subp.α14)
            )
        MOI.submit(mas.m, MOI.LazyConstraint(cb_data), cut)
        push!(fcuts, (α5=subp.α5,α6=subp.α6,α7=subp.α7,α8=subp.α8,α9=subp.α9,α10=subp.α10,α12=subp.α12,α13=subp.α13,α14=subp.α14))
    end
    return 
end
 
rhs1 = [[(dt.a[i,p]+dt.vij[m]) for i=1:1 for j=1:1 for m=1:dt.Mij[i,j] for p=1:5];
[dt.a[i,p]+dt.vij[sum(dt.Mij[i,1:j-1])+m] for i=1:1 for j=2:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5];
[dt.a[i,p]+dt.vij[sum(dt.Mij[1:i-1,:])+m] for i=2:dt.N["supplier"] for j=1:1 for m=1:dt.Mij[i,j] for p=1:5];
[dt.a[i,p]+dt.vij[sum(dt.Mij[1:i-1,:])+sum(dt.Mij[i,1:j-1])+m] for i=2:dt.N["supplier"] for j=2:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5];
[dt.vjk[m] for j=1 for k=1 for m=1:dt.Mjk[j,k] for p=1:5];
[dt.vjk[sum(dt.Mjk[j,1:k-1])+m] for j=1,k=2:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5];
[dt.vjk[sum(dt.Mjk[1:j-1,:])+m] for j=2:dt.N["plant"] for k=1:1 for m=1:dt.Mjk[j,k] for p=1:5];
[dt.vjk[sum(dt.Mjk[1:j-1,:])+sum(dt.Mjk[j,1:k-1])+m] for j=2:dt.N["plant"] for k=2:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5];
[dt.vkl[m] for k=1:1 for l=1:1 for m=1:dt.Mkl[k,l] for p=1:5];
[dt.vkl[sum(dt.Mkl[k,1:l-1])+m] for k=1:1 for l=2:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5];
[dt.vkl[sum(dt.Mkl[1:k-1,:])+m] for k=2:dt.N["distribution"] for l=1:1 for m=1:dt.Mkl[k,l] for p=1:5];
[dt.vkl[sum(dt.Mkl[1:k-1,:])+sum(dt.Mkl[k,1:l-1])+m] for k=2:dt.N["distribution"] for l=2:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5];
[dt.e[5*(t-1)+p] for t=1:2 for p=1:5];
[dt.e[5*2*(j-1)+((t-1)*5)+p] for j=2:dt.N["plant"] for t=1:2 for p=1:5];
[dt.e[(dt.N["plant"]*2*5)+5*(t-1)+p] for t=1:2 for p=1:5];
[dt.e[(dt.N["plant"]*2*5)+5*2*(k-1)+((t-1)*5)+p] for k=2:dt.N["distribution"] for t=1:2 for p=1:5]]

rhs2 = [[(dt.b[i,p]+dt.rij[m]) for i=1:1 for j=1:1 for m=1:dt.Mij[i,j] for p=1:5];
[dt.b[i,p]+dt.rij[sum(dt.Mij[i,1:j-1])+m] for i=1:1 for j=2:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5];
[dt.b[i,p]+dt.rij[sum(dt.Mij[1:i-1,:])+m] for i=2:dt.N["supplier"] for j=1:1 for m=1:dt.Mij[i,j] for p=1:5];
[dt.b[i,p]+dt.rij[sum(dt.Mij[1:i-1,:])+sum(dt.Mij[i,1:j-1])+m] for i=2:dt.N["supplier"] for j=2:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5];
[dt.rjk[m] for j=1 for k=1 for m=1:dt.Mjk[j,k] for p=1:5];
[dt.rjk[sum(dt.Mjk[j,1:k-1])+m] for j=1,k=2:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5];
[dt.rjk[m] for j=2:dt.N["plant"] for k=1:1 for m=1:dt.Mjk[j,k] for p=1:5];
[dt.rjk[sum(dt.Mjk[j,1:k-1])+m] for j=2:dt.N["plant"] for k=2:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5];
[dt.rkl[m] for k=1:1 for l=1:1 for m=1:dt.Mkl[k,l] for p=1:5];
[dt.rkl[sum(dt.Mkl[k,1:l-1])+m] for k=1:1 for l=2:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5];
[dt.rkl[m] for k=2:dt.N["distribution"] for l=1:1 for m=1:dt.Mkl[k,l] for p=1:5];
[dt.rkl[sum(dt.Mkl[k,1:l-1])+m] for k=2:dt.N["distribution"] for l=2:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5];
[dt.q[5*(t-1)+p] for t=1:2 for p=1:5];
[dt.q[5*2*(j-1)+((t-1)*5)] for j=2:dt.N["plant"] for t=1:2 for p=1:5];
[dt.q[(dt.N["plant"]*2*5)+5*(t-1)+p] for t=1:2 for p=1:5];
[dt.q[(dt.N["plant"]*2*5)+5*2*(k-1)+((t-1)*5)+p] for k=2:dt.N["distribution"] for t=1:2 for p=1:5]]
objc = [dt.N["cap"];dt.N["cad"];dt.Vij.-dt.bigM;dt.Vjk.-dt.bigM; [dt.bigM for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:sum(dt.Mkl)]] 


w = .5; fcuts = []; pi1,pi2 = [],[]; subtime = []; 
mas = build_master(.5); #set_optimizer_attribute(mas.m, "CPXPARAM_TimeLimit", 1);
dsp = DualSubP(.5,150); MOI.set(mas.m, MOI.LazyConstraintCallback(), benders_decomposition);
optimize!(mas.m);
termination_status(mas.m)