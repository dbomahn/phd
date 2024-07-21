using CSV,DelimitedFiles,Distributions,CPLEX,MathOptInterface,JuMP,MathOptInterface,LinearAlgebra,CPUTime,DataFrames,JLD2,vOptGeneric
import MathOptInterface as MOI

using PlotlyJS,DataFrames,Colors,CSV,JLD2
#############################        2D plot      ###########################
layout = Layout(
    title="Test",
    xaxis_title="Cost",
    yaxis_title="CO2 emission",
    legend_title="Legend Title",
    font=attr(
        family="Courier New, monospace",
        size=18
    )
)
# % $CO_2$ emission unit has changed from kg to g ($\times$ 10^3) 

struct Data1d
    insfile::String; 
    I::Int; J::Int; K::Int; L::Int; Jmax::Int; Kmax::Int; 
    Mij::Array{}; Mjk::Array{}; Mkl::Array{};
    demand::Array{}; bigM::Int; cas::Array{}; cap::Array{}; cad::Array{};  
    fc::Array{}; tcp::Array{}; tcd::Array{}; tcc::Array{}; cep::Array{}; ced::Array{}; cec::Array{}; 
    vcs::Array{}; vcp::Array{}; vcd::Array{}; ves::Array{}; vep::Array{}; ved::Array{};
    Lcapasp::Array{}; Lcapapd::Array{}; gij::Array{}; gjk::Array{}; gkl::Array{};

    function Data1d(insfile)
        dt = readdlm(insfile)
        I,J,K,L,Jmax,Kmax = dt[1:6]
        mijloc = findall(x->x=="Mij", dt)[1][1]
        mjkloc = findall(x->x=="Mjk", dt)[1][1]
        mklloc = findall(x->x=="Mkl", dt)[1][1]
        dloc = findall(x->x=="demand", dt)[1][1]
        casloc = findall(x->x=="cas", dt)[1][1]
        caploc = findall(x->x=="cap", dt)[1][1]
        cadloc = findall(x->x=="cad", dt)[1][1]
        fcloc = findall(x->x=="fixedcost", dt)[1][1]
        tcpl = findall(x->x=="tcp", dt)[1][1]
        tcdl = findall(x->x=="tcd", dt)[1][1]
        tccl = findall(x->x=="tcc", dt)[1][1]
        cepl = findall(x->x=="cep", dt)[1][1]
        cedl = findall(x->x=="ced", dt)[1][1]
        cecl = findall(x->x=="cec", dt)[1][1]
        vcsl = findall(x->x=="vcs", dt)[1][1]
        vcpl = findall(x->x=="vcp", dt)[1][1]
        vcdl = findall(x->x=="vcd", dt)[1][1]
        vesl = findall(x->x=="ves", dt)[1][1]
        vepl = findall(x->x=="vep", dt)[1][1]
        vedl = findall(x->x=="ved", dt)[1][1]
        Lcapaspl = findall(x->x=="Lcapasp", dt)[1][1]
        Lcapapdl = findall(x->x=="Lcapapd", dt)[1][1]
        gijl = findall(x->x=="fixedcostModesp", dt)[1][1]
        gjkl = findall(x->x=="fixedcostModepd", dt)[1][1]
        gkll = findall(x->x=="fixedcostModedc", dt)[1][1]
        #######  Assigning data 
        Mij = dt[mijloc+1:mjkloc-1,:][1:I,1:J]
        Mjk = dt[mjkloc+1:mklloc-1,:][1:J,1:K]
        Mkl = dt[mklloc+1:dloc-1,:][1:K,1:L]
        demand = dt[dloc+1:dloc+L]
        bigM = sum(demand)
        cas = dt[casloc+1:caploc-1]; cap = dt[caploc+1:cadloc-1]; cad = dt[cadloc+1:fcloc-1]
        fc = reshape(transpose(dt[fcloc+1:tcpl-1,1:2]),(tcpl-fcloc-1)*2,1);
        tcp0 = dt[tcpl+1:tcdl-1,:]; tcd0 = dt[tcdl+1:tccl-1,1:2]; tcc0 = dt[tccl+1:cepl-1,1:2];
        cep0 = dt[cepl+1:cedl-1,1:2]; ced0 = dt[cedl+1:cecl-1,1:2]; cec0 = dt[cecl+1:vcsl-1,:];
        vcs = dt[vcsl+1:vcpl-1]; vcp = reshape(transpose(dt[vcpl+1:vcdl-1,1:2]),(vcdl-vcpl-1)*2,1); vcd = reshape(transpose(dt[vcdl+1:vesl-1,1:2]),(vesl-vcdl-1)*2,1);
        ves = dt[vesl+1:vepl-1]*10^3; vep = reshape(transpose(dt[vepl+1:vedl-1,1:2]*10^3),(vedl-vepl-1)*2,1); ved = reshape(transpose(dt[vedl+1:Lcapaspl-1,1:2]*10^3),(Lcapaspl-vedl-1)*2,1);
        Lcapasp0 = dt[Lcapaspl+1:Lcapapdl-1,:]; Lcapapd0 = dt[Lcapapdl+1:gijl-1,:];
        gij0 = dt[gijl+1:gjkl-1,:]; gjk0 = dt[gjkl+1:gkll-1,:]; gkl0 = dt[gkll+1:end,:]
        tcp,tcd,tcc,cep,ced,cec,Lcapasp,Lcapapd,gij,gjk,gkl = [],[],[],[],[],[],[],[],[],[],[]
        ########### Adjusting the structure
        ct = 1
        for i=1:I
            tmptc = []; tmpce = []; tmpL = []; tmpg = [];
            for j=1:J
                if Mij[i,j]==1
                    append!(tmptc, tcp0[ct]); append!(tmpce, cep0[ct]);  append!(tmpL, [Lcapasp0[ct]]); append!(tmpg, [gij0[ct]]); 
                else
                    append!(tmptc, tcp0[ct,1:2]); append!(tmpce, cep0[ct,1:2]); append!(tmpL, Lcapasp0[ct,1:2]); append!(tmpg, gij0[ct,1:2]); 
                end
                ct+=1
            end
            tmpce = [round.(p*10^3; digits=1) for p in tmpce]
            append!(tcp,tmptc); append!(cep,tmpce);  append!(Lcapasp,tmpL); append!(gij,tmpg)
        end

        ct = 1
        for j=1:J
            tmptc = []; tmpce = []; tmpL = []; tmpg = [];
            for k=1:K
                if Mjk[j,k]==1
                    append!(tmptc, tcd0[ct]); append!(tmpce, ced0[ct]); append!(tmpL, Lcapapd0[ct]); append!(tmpg, gjk0[ct]); 
                else
                    append!(tmptc, tcd0[ct,1:2]); append!(tmpce, ced0[ct,1:2]); append!(tmpL, Lcapapd0[ct,1:2]); append!(tmpg, gjk0[ct,1:2]);
                end
                ct+=1
            end
            tmpce = [round.(d*10^3; digits=1) for d in tmpce]
            append!(tcd,tmptc);   append!(ced,tmpce); append!(Lcapapd,tmpL); append!(gjk,tmpg)
        end

        ct = 1
        for k=1:K
            tmptc = []; tmpce = []; tmpg = [];
            for l=1:L
                if Mkl[k,l]==1
                    append!(tmptc, tcc0[ct]);  append!(tmpce, cec0[ct]); append!(tmpg, gkl0[ct]); 
                else
                    append!(tmptc, tcc0[ct,1:2]); append!(tmpce, cec0[ct,1:2]); append!(tmpg, gkl0[ct,1:2]);
                end
                ct+=1
            end
            tmpce = [round.(c*10^3; digits=1) for c in tmpce]
            append!(tcc,tmptc); append!(cec,tmpce); append!(gkl,tmpg)
        end
        new(insfile,I,J,K,L,Jmax,Kmax,Mij,Mjk,Mkl,demand,bigM,cas,cap,cad,fc,tcp,tcd,tcc,cep,ced,cec,vcs,vcp,vcd,ves,vep,ved,Lcapasp,Lcapapd,gij,gjk,gkl)
    end
end
function NDfilter(Pobj)
    copyobj = Dict();
    for i=1:length(Pobj)
        copyobj[i] = Pobj[i]
    end
    for i=1:length(Pobj)-1
        for j=i+1:length(Pobj)
            if all(round.(Pobj[i];digits=4) .>= round.(Pobj[j];digits=4)) == true #dominated by PF[j]
                copyobj[i]=nothing; break
            elseif all(round.(Pobj[j];digits=4) .>= round.(Pobj[i];digits=4)) == true
                copyobj[j]=nothing; 
            end
        end
    end
    finalobj = filter!(a->a!=nothing, collect(values(copyobj)))
    return sort!(finalobj)
end
function dominated(x,P)
    st = false
    for k=1:length(P)
        if all( x .>= P[k])# && any(x > P[k])
            st=true; break
        else
            continue
        end
    end
    return st
end

mutable struct Master
    y::Vector{VariableRef}
    uij::Vector{VariableRef}
    ujk::Vector{VariableRef}
    ukl::Vector{VariableRef}
    θ::VariableRef
    m::Model
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
function Master(w,dt,tol)
    mas = direct_model(CPLEX.Optimizer()); set_silent(mas)
    set_optimizer_attribute(mas, "CPXPARAM_Preprocessing_Presolve", 0) # presolve must be turned off    
    MOI.set(mas, MOI.RelativeGapTolerance(), tol) 
    MOI.set(mas, MOI.RawOptimizerAttribute("CPXPARAM_MIP_Strategy_Search"), 1) # Apply traditional branch and cut strategy; disable dynamic search
    MOI.set(mas, MOI.NumberOfThreads(), 1); 
    @variable(mas, y[1:(dt.J+dt.K)*2], Bin)
    @variable(mas, uij[1:sum(dt.Mij)], Bin)
    @variable(mas, ujk[1:sum(dt.Mjk)], Bin)
    @variable(mas, ukl[1:sum(dt.Mkl)], Bin)
    @variable(mas, θ>= -10^8);

    @constraint(mas,[j=1:dt.J+dt.K], sum(y[2*(j-1)+1:2*(j-1)+2]) <= 1)
    @constraints(mas,begin
        sum(uij[1:dt.Mij[1,1]]) <= 1 #[1,1]
        [j=2:dt.J], sum(uij[sum(dt.Mij[1,1:j-1])+1:sum(dt.Mij[1,1:j-1])+dt.Mij[1,j]]) <= 1 # [1,2:6]
        [i=2:dt.I], sum(uij[sum(dt.Mij[1:i-1,:])+1:sum(dt.Mij[1:i-1,:])+dt.Mij[i,1]])<= 1 # [2:6,1]
        [i=2:dt.I,j=2:dt.J], sum(uij[sum(dt.Mij[1:i-1,:])+sum(dt.Mij[i,1:j-1])+1:sum(dt.Mij[1:i-1,:])+sum(dt.Mij[i,1:j-1])+dt.Mij[i,j]])<= 1 #sum(y[j,:])
        
        sum(ujk[1:dt.Mjk[1,1]]) <= 1 
        [k=2:dt.K], sum(ujk[sum(dt.Mjk[1,1:k-1])+1:sum(dt.Mjk[1,1:k-1])+dt.Mjk[1,k]]) <= 1 
        [j=2:dt.J], sum(ujk[sum(dt.Mjk[1:j-1,:])+1:sum(dt.Mjk[1:j-1,:])+dt.Mjk[j,1]])<= 1 # [2:6,1]
        [j=2:dt.J,k=2:dt.K],  sum(ujk[sum(dt.Mjk[1:j-1,:])+sum(dt.Mjk[j,1:j-1])+1:sum(dt.Mjk[1:j-1,:])+sum(dt.Mjk[j,1:j-1])+dt.Mjk[j,k]]) <= 1 #(sum(y[j,:])+sum(y[dt.J+k,:]))/2

        sum(ukl[1:dt.Mkl[1,1]]) <= 1 
        [l=2:dt.L], sum(ukl[sum(dt.Mkl[1,1:l-1])+1:sum(dt.Mkl[1,1:l-1])+dt.Mkl[1,l]]) <= 1 #sum(y[dt.J+1,:])
        [k=2:dt.K], sum(ukl[sum(dt.Mkl[1:k-1,:])+1:sum(dt.Mkl[1:k-1,:])+dt.Mkl[k,1]])<= 1 
        [k=2:dt.K,l=2:dt.L], sum(ukl[sum(dt.Mkl[1:k-1,:])+sum(dt.Mkl[k,1:l-1])+1:sum(dt.Mkl[1:k-1,:])+sum(dt.Mkl[k,1:l-1])+dt.Mkl[k,l]])<= 1 #sum(y[dt.J+k,:])
    end);
    @constraint(mas, sum(y[1:dt.J*2]) <= dt.Jmax);
    @constraint(mas, sum(y[dt.J*2+1:end]) <= dt.Kmax);
    @objective(mas, Min, w*(sum(dt.fc.*y)+dot(dt.gij,uij) + dot(dt.gjk,ujk) + dot(dt.gkl,ukl)) + θ)
    return Master(y,uij,ujk,ukl,θ,mas)
end
function DualSubP(w1,w2,dt,tol)
    sub = direct_model(CPLEX.Optimizer()); set_silent(sub)
    set_optimizer_attribute(sub, "CPXPARAM_Preprocessing_Presolve", 0) # presolve must be turned off
    set_optimizer_attribute(sub, "CPXPARAM_LPMethod" ,1) #simplex method
    MOI.set(sub, MOI.NumberOfThreads(), 1) # of threads
    MOI.set(sub, MOI.RelativeGapTolerance(), tol) 

    @variables(sub, begin
        α1[1:dt.J]
        α2[1:dt.K]
        α3[1:dt.J]
        α4[1:dt.K]
    end);
    @variable(sub, α5[1:dt.L] >= 0);
    @variable(sub, α6[1:dt.I] >= 0);
    @variable(sub, α7[1:dt.J*2] >= 0);
    @variable(sub, α8[1:dt.K*2] >= 0);
    @variable(sub, α9[1:sum(dt.Mij)] >= 0);
    @variable(sub, α10[1:sum(dt.Mjk)] >= 0);
    @variable(sub, α12[1:sum(dt.Mij)] >= 0);
    @variable(sub, α13[1:sum(dt.Mjk)] >= 0);
    @variable(sub, α14[1:sum(dt.Mkl)] >= 0);


    @constraints(sub, begin
        con1[i=1,j=1,m=1:dt.Mij[i,j]], α1[j]-α3[j]-α6[i]+α9[m]-α12[m] <= w1*(dt.vcs[i]+dt.tcp[m])+w2*(dt.ves[i]+dt.cep[m])
        con2[i=1,j=2:dt.J,m=1:dt.Mij[i,j]], α1[j]-α3[j]-α6[i]+α9[sum(dt.Mij[i,1:j-1])+m]-α12[sum(dt.Mij[i,1:j-1])+m] <= w1*(dt.vcs[i]+dt.tcp[sum(dt.Mij[i,1:j-1])+m])+w2*(dt.ves[i]+dt.cep[sum(dt.Mij[i,1:j-1])+m])
        con3[i=2:dt.I,j=1,m=1:dt.Mij[i,j]], α1[j]-α3[j]-α6[i]+α9[sum(dt.Mij[1:i-1,:])+sum(dt.Mij[i,1])+m]-α12[sum(dt.Mij[1:i-1,:])+sum(dt.Mij[i,1])+m] <= w1*(dt.vcs[i]+dt.tcp[sum(dt.Mij[1:i-1,:])+m])+w2*(dt.ves[i]+dt.cep[sum(dt.Mij[1:i-1,:])+m])
        con4[i=2:dt.I,j=2:dt.J,m=1:dt.Mij[i,j]], α1[j]-α3[j]-α6[i]+α9[sum(dt.Mij[1:i-1,:])+sum(dt.Mij[i,1:j-1])+m]-α12[sum(dt.Mij[1:i-1,:])+sum(dt.Mij[i,1:j-1])+m] <= w1*(dt.vcs[i]+dt.tcp[sum(dt.Mij[1:i-1,:])+sum(dt.Mij[i,1:j-1])+m])+w2*(dt.ves[i]+dt.cep[sum(dt.Mij[1:i-1,:])+sum(dt.Mij[i,1:j-1])+m])
    end)

    @constraint(sub, con5[j=1,k=1,m=1:dt.Mjk[j,k]], -α1[j]+α2[k]-α4[k]+α10[m]-α13[m] <= w1*(dt.tcd[m])+w2*(dt.ced[m]))
    @constraint(sub, con6[j=1,k=2:dt.K,m=1:dt.Mjk[j,k]], -α1[j]+α2[k]-α4[k]+α10[m]-α13[m] <= w1*dt.tcd[sum(dt.Mjk[j,1:k-1])+m]+w2*(dt.ced[sum(dt.Mjk[j,1:k-1])+m]))
    @constraint(sub, con7[j=2:dt.J,k=1,m=1:dt.Mjk[j,k]], -α1[j]+α2[k]-α4[k]+α10[sum(dt.Mjk[1:j-1,:])+sum(dt.Mjk[j,1:j-1])+m]-α13[sum(dt.Mjk[1:j-1,:])+sum(dt.Mjk[j,1:j-1])+m] <= w1*dt.tcd[sum(dt.Mjk[1:j-1,:])+m]+w2*dt.ced[m])
    @constraint(sub, con8[j=2:dt.J,k=2:dt.K,m=1:dt.Mjk[j,k]], -α1[j]+α2[k]-α4[k]+α10[sum(dt.Mjk[1:j-1,:])+sum(dt.Mjk[j,1:j-1])+m]-α13[sum(dt.Mjk[1:j-1,:])+sum(dt.Mjk[j,1:j-1])+m] <= w1*dt.tcd[sum(dt.Mjk[1:j-1,:])+sum(dt.Mjk[j,1:k-1])+m]+w2*dt.ced[sum(dt.Mjk[j,1:k-1])+m])

    @constraint(sub, con9[k=1,l=1,m=1:dt.Mkl[k,l]],-α2[k]+α5[l]-α14[m] <= w1*dt.tcc[m]+w2*dt.cec[m])
    @constraint(sub, con10[k=1,l=2:dt.L,m=1:dt.Mkl[k,l]], -α2[k]+α5[l]-α14[sum(dt.Mkl[k,1:l-1])+m] <= w1*dt.tcc[sum(dt.Mkl[k,1:l-1])+m]+w2*dt.cec[sum(dt.Mkl[k,1:l-1])+m])
    @constraint(sub, con11[k=2:dt.K,l=1,m=1:dt.Mkl[k,l]], -α2[k]+α5[l]-α14[sum(dt.Mkl[1:k-1,:])+m] <= w1*dt.tcc[sum(dt.Mkl[1:k-1,:])+m]+w2*dt.cec[m])
    @constraint(sub, con12[k=2:dt.K,l=2:dt.L,m=1:dt.Mkl[k,l]],-α2[k]+α5[l]-α14[sum(dt.Mkl[1:k-1,:])+sum(dt.Mkl[k,1:l-1])+m] <= w1*dt.tcc[sum(dt.Mkl[1:k-1,:])+sum(dt.Mkl[k,1:l-1])+m]+w2*dt.cec[sum(dt.Mkl[k,1:l-1])+m])

    @constraint(sub, con13[j=1:dt.J,t=1:2], α3[j]-α7[2*(j-1)+t] <= w1*dt.vcp[2*(j-1)+t]+w2*dt.vep[2*(j-1)+t])
    @constraint(sub, con14[k=1:dt.K,t=1:2], α4[k]-α8[2*(k-1)+t] <= w1*dt.vcd[2*(k-1)+t] + w2*dt.ved[2*(k-1)+t])
    return DualSubP(α1,α2,α3,α4,α5,α6,α7,α8,α9,α10,α12,α13,α14,sub)
end
function solve_dsp(dsp::DualSubP,yb,ubij,ubjk,ubkl)
    @objective(dsp.m, Max, dot(dt.demand,dsp.α5)-dot(dt.cas,dsp.α6)-
        sum(dt.cap[j]*yb[2*(j-1)+t]*dsp.α7[2*(j-1)+t] for j=1:dt.J for t=1:2)-sum(dt.cad[k]*yb[2*dt.J+k]*dsp.α8[2*(k-1)+t] for k=1:dt.K for t=1:2)+
        dot(dt.Lcapasp.*ubij,dsp.α9)+dot(dt.Lcapapd.*ubjk,dsp.α10) -
        dot(dt.bigM.*ubij,dsp.α12)-dot(dt.bigM.*ubjk,dsp.α13) -dot(dt.bigM.*ubkl,dsp.α14)
        );
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
function BD_lazy_callback(cb_data)
    yb = callback_value.(cb_data, mas.y);    ubij = callback_value.(cb_data, mas.uij);    ubjk = callback_value.(cb_data, mas.ujk)
    ubkl = callback_value.(cb_data, mas.ukl);    θb = callback_value(cb_data, mas.θ)
    subp = solve_dsp(dsp,yb,ubij,ubjk,ubkl);    
    if subp.res == :OptimalityCut
        # @info "Optimality cut found"
        if dt.insfile[end-11:end-4] == "test01S1" 
            num = 4
        elseif dt.insfile[end-11:end-4] == "test01S3"
            num = 4
        elseif dt.insfile[end-11:end-4] == "test01S4"
            num = 5
        elseif dt.insfile[end-11:end-4] == "test02S3"
            num = 4
        elseif dt.insfile[end-11:end-4] == "test02S4"
            num = 3
        elseif dt.insfile[end-11:end-4] == "test03S2"
            num = 4
        elseif dt.insfile[end-11:end-4] == "test03S4"
            num = 4
        elseif dt.insfile[end-11:end-4] == "test04S1"
            num = 4
        else 
            num = 6
        end
        if θb ≥ round(subp.obj; digits=num)
            return
        else
            cut = @build_constraint( mas.θ ≥ dot(subp.α5,dt.demand)-dot(dt.cas,subp.α6)-
                sum(dt.cap[j]*mas.y[2*(j-1)+t]*subp.α7[2*(j-1)+t] for j=1:dt.J for t=1:2)-
                sum(dt.cad[k]*mas.y[2*dt.J+k]*subp.α8[2*(k-1)+t] for k=1:dt.K for t=1:2)+
                dot(dt.Lcapasp.*mas.uij,subp.α9)+dot(dt.Lcapapd.*mas.ujk,subp.α10) -
                dot(dt.bigM.*mas.uij,subp.α12)-dot(dt.bigM.*mas.ujk,subp.α13) -dot(dt.bigM.*mas.ukl,subp.α14)
                )
            MOI.submit(mas.m, MOI.LazyConstraint(cb_data), cut); push!(tmpno,1)
            # push!(ocuts, (obj=subp.obj,α5=subp.α5,α6=subp.α6,α7=subp.α7,α8=subp.α8,α9=subp.α9,α10=subp.α10,α12=subp.α12,α13=subp.α13,α14=subp.α14))
        end
    elseif subp.res == :FeasibilityCut
        cut = @build_constraint( 0 ≥ dot(subp.α5,dt.demand)-dot(dt.cas,subp.α6)-
            sum(dt.cap[j]*mas.y[2*(j-1)+t]*subp.α7[2*(j-1)+t] for j=1:dt.J for t=1:2)-
            sum(dt.cad[k]*mas.y[2*dt.J+k]*subp.α8[2*(k-1)+t] for k=1:dt.K for t=1:2)+
            dot(dt.Lcapasp.*mas.uij,subp.α9)+dot(dt.Lcapapd.*mas.ujk,subp.α10) -
            dot(dt.bigM.*mas.uij,subp.α12)-dot(dt.bigM.*mas.ujk,subp.α13) -dot(dt.bigM.*mas.ukl,subp.α14)
            )
        MOI.submit(mas.m, MOI.LazyConstraint(cb_data), cut); push!(tmpnf,1)
        # push!(fcuts, (α5=subp.α5,α6=subp.α6,α7=subp.α7,α8=subp.α8,α9=subp.α9,α10=subp.α10,α12=subp.α12,α13=subp.α13,α14=subp.α14))
    # else
    #     println("Neither case")
    end
    push!(tmpit,1)
    return 
end
function Basis_lazy_callback(cb_data)
    yb = callback_value.(cb_data, mas.y);    ubij = callback_value.(cb_data, mas.uij);    ubjk = callback_value.(cb_data, mas.ujk)
    ubkl = callback_value.(cb_data, mas.ukl);    θb = callback_value(cb_data, mas.θ)
    subp = solve_dsp(dsp,yb,ubij,ubjk,ubkl,dt);    
    if subp.res == :OptimalityCut
        # @info "Optimality cut found"
        if dt.insfile[end-11:end-4] == "test01S1" 
            num = 4
        elseif dt.insfile[end-11:end-4] == "test01S3"
            num = 4
        elseif dt.insfile[end-11:end-4] == "test01S4"
            num = 5
        elseif dt.insfile[end-11:end-4] == "test02S3"
            num = 4
        elseif dt.insfile[end-11:end-4] == "test02S4"
            num = 3
        elseif dt.insfile[end-11:end-4] == "test03S2"
            num = 4
        elseif dt.insfile[end-11:end-4] == "test03S4"
            num = 4
        elseif dt.insfile[end-11:end-4] == "test04S1"
            num = 4
        else 
            num = 6
        end
        if θb ≥ round(subp.obj; digits=num)
            return
        else
            head0 = Vector{Cint}(undef, numcon); bx = Vector{Cdouble}(undef, numcon)
            CPLEX.CPXgetbhead(backend(dsp.m).env, backend(dsp.m).lp, head0, bx);head = head0.+1;
            c1 = copy(rhs1); c2 = copy(rhs2); 
            CPLEX.CPXftran(backend(dsp.m).env, backend(dsp.m).lp, c1)
            CPLEX.CPXftran(backend(dsp.m).env, backend(dsp.m).lp, c2)
            pi1,pi2 = zeros(num_variables(dsp.m)), zeros(num_variables(dsp.m))
            for (idx,i) in enumerate(head)
                if i >0
                    pi1[i] = c1[idx]
                    pi2[i] = c2[idx]
                end
            end
            cut = @build_constraint( mas.θ ≥ dot(subp.α5,dt.demand)-dot(dt.cas,subp.α6)-
                sum(dt.cap[j]*mas.y[2*(j-1)+t]*subp.α7[2*(j-1)+t] for j=1:dt.J for t=1:2)-
                sum(dt.cad[k]*mas.y[2*dt.J+k]*subp.α8[2*(k-1)+t] for k=1:dt.K for t=1:2)+
                dot(dt.Lcapasp.*mas.uij,subp.α9)+dot(dt.Lcapapd.*mas.ujk,subp.α10) -
                dot(dt.bigM.*mas.uij,subp.α12)-dot(dt.bigM.*mas.ujk,subp.α13) -dot(dt.bigM.*mas.ukl,subp.α14)
                )
            
            MOI.submit(mas.m, MOI.LazyConstraint(cb_data), cut); push!(tmpno,1)
            push!(dfo, [w,copy(pi1),copy(pi2) ,subp.α5,subp.α6])
        end
    elseif subp.res == :FeasibilityCut
        cut = @build_constraint( 0 ≥ dot(subp.α5,dt.demand)-dot(dt.cas,subp.α6)-
            sum(dt.cap[j]*mas.y[2*(j-1)+t]*subp.α7[2*(j-1)+t] for j=1:dt.J for t=1:2)-
            sum(dt.cad[k]*mas.y[2*dt.J+k]*subp.α8[2*(k-1)+t] for k=1:dt.K for t=1:2)+
            dot(dt.Lcapasp.*mas.uij,subp.α9)+dot(dt.Lcapapd.*mas.ujk,subp.α10) -
            dot(dt.bigM.*mas.uij,subp.α12)-dot(dt.bigM.*mas.ujk,subp.α13) -dot(dt.bigM.*mas.ukl,subp.α14)
            )
        MOI.submit(mas.m, MOI.LazyConstraint(cb_data), cut)
        push!(dff, [w, (α5=subp.α5,α6=subp.α6,α7=subp.α7,α8=subp.α8,α9=subp.α9,α10=subp.α10,α12=subp.α12,α13=subp.α13,α14=subp.α14)])
        push!(tmpnf,1)
    else
        println("Neither case in BD")
    end
    push!(tmpit,1)
    return 
end

function getobjdual(mp,dp,dt)
    y = value.(mp[:y])
    uij = value.(mp[:uij])
    ujk = value.(mp[:ujk])
    ukl = value.(mp[:ukl])

    obj1 = sum(dt.fc.*y)+dot(dt.gij,uij) + dot(dt.gjk,ujk) + dot(dt.gkl,ukl) +
        sum((dt.vcs[i]+dt.tcp[m])*-dual.(dp[:con1])[i,j,m] for i=1 for j=1 for m=1:dt.Mij[i,j]) +
        sum((dt.vcs[i]+dt.tcp[sum(dt.Mij[i,1:j-1])+m])*-dual.(dp[:con2])[i,j,m] for i=1 for j=2:dt.J for m=1:dt.Mij[i,j]) +
        sum((dt.vcs[i]+dt.tcp[sum(dt.Mij[1:i-1,:])+m])*-dual.(dp[:con3])[i,j,m] for i=2:dt.I for j=1 for m=1:dt.Mij[i,j]) +
        sum((dt.vcs[i]+dt.tcp[sum(dt.Mij[1:i-1,:])+sum(dt.Mij[i,1:j-1])+m])*-dual.(dp[:con4])[i,j,m] for i=2:dt.I for j=2:dt.J for m=1:dt.Mij[i,j]) +
        sum(dt.tcd[m]*-dual.(dp[:con5])[j,k,m] for j=1 for k=1 for m=1:dt.Mjk[j,k]) +
        sum(dt.tcd[sum(dt.Mjk[j,1:k-1])+m]*-dual.(dp[:con6])[j,k,m] for j=1 for k=2:dt.K for m=1:dt.Mjk[j,k]) +
        sum(dt.tcd[sum(dt.Mjk[1:j-1,:])+m]*-dual.(dp[:con7])[j,k,m] for j=2:dt.J for k=1 for m=1:dt.Mjk[j,k]) +
        sum(dt.tcd[sum(dt.Mjk[1:j-1,:])+sum(dt.Mjk[j,1:k-1])+m]*-dual.(dp[:con8])[j,k,m] for j=2:dt.J for k=2:dt.K for m=1:dt.Mjk[j,k]) +
        sum(dt.tcc[m]*-dual.(dp[:con9])[k,l,m] for k=1 for l=1 for m=1:dt.Mkl[k,l])+
        sum(dt.tcc[sum(dt.Mkl[k,1:l-1])+m]*-dual.(dp[:con10])[k,l,m] for k=1 for l=2:dt.L for m=1:dt.Mkl[k,l])+
        sum(dt.tcc[sum(dt.Mkl[1:k-1,:])+m]*-dual.(dp[:con11])[k,l,m] for k=2:dt.K for l=1 for m=1:dt.Mkl[k,l])+
        sum(dt.tcc[sum(dt.Mkl[1:k-1,:])+sum(dt.Mkl[k,1:l-1])+m]*-dual.(dp[:con12])[k,l,m] for k=2:dt.K for l=2:dt.L for m=1:dt.Mkl[k,l])+
        dot(dt.vcp, reshape(transpose(-dual.(dp[:con13])),dt.J*2,1))+
        dot(dt.vcd, reshape(transpose(-dual.(dp[:con14])),dt.K*2,1))
                
    obj2 = sum([dt.ves[i]+dt.cep[m]*-dual.(dp[:con1])[i,j,m] for i=1 for j=1 for m=1:dt.Mij[i,j]])+
        sum([dt.ves[i]+dt.cep[sum(dt.Mij[i,1:j-1])+m]*-dual.(dp[:con2])[i,j,m] for i=1 for j=2:dt.J for m=1:dt.Mij[i,j]])+
        sum([dt.ves[i]+dt.cep[sum(dt.Mij[1:i-1,:])+m]*-dual.(dp[:con3])[i,j,m] for i=2:dt.I for j=1 for m=1:dt.Mij[i,j]])+
        sum([dt.ves[i]+dt.cep[sum(dt.Mij[1:i-1,:])+sum(dt.Mij[i,1:j-1])+m]*-dual.(dp[:con4])[i,j,m] for i=2:dt.I for j=2:dt.J for m=1:dt.Mij[i,j]])+        
        sum(dt.ced[m]*-dual.(dp[:con5])[j,k,m] for j=1 for k=1 for m=1:dt.Mjk[j,k]) +
        sum(dt.ced[sum(dt.Mjk[j,1:k-1])+m]*-dual.(dp[:con6])[j,k,m] for j=1 for k=2:dt.K for m=1:dt.Mjk[j,k]) +
        sum(dt.ced[sum(dt.Mjk[1:j-1,:])+m]*-dual.(dp[:con7])[j,k,m] for j=2:dt.J for k=1 for m=1:dt.Mjk[j,k]) +
        sum(dt.ced[sum(dt.Mjk[1:j-1,:])+sum(dt.Mjk[j,1:k-1])+m]*-dual.(dp[:con8])[j,k,m] for j=2:dt.J for k=2:dt.K for m=1:dt.Mjk[j,k]) +    
        sum(dt.cec[m]*-dual.(dp[:con9])[k,l,m] for k=1 for l=1 for m=1:dt.Mkl[k,l])+
        sum(dt.cec[sum(dt.Mkl[k,1:l-1])+m]*-dual.(dp[:con10])[k,l,m] for k=1 for l=2:dt.L for m=1:dt.Mkl[k,l])+
        sum(dt.cec[sum(dt.Mkl[1:k-1,:])+m]*-dual.(dp[:con11])[k,l,m] for k=2:dt.K for l=1 for m=1:dt.Mkl[k,l])+
        sum(dt.cec[sum(dt.Mkl[1:k-1,:])+sum(dt.Mkl[k,1:l-1])+m]*-dual.(dp[:con12])[k,l,m] for k=2:dt.K for l=2:dt.L for m=1:dt.Mkl[k,l])+
        dot(dt.vep, reshape(transpose(-dual.(dp[:con13])),dt.J*2,1))+
        dot(dt.ved, reshape(transpose(-dual.(dp[:con14])),dt.K*2,1))

    return round(obj1),round(obj2)#obj1,obj2 #,y,uij,ujk,ukl,xij,xjk,xkl,h 
end
function smSCND_MIP(dt,tol)
    smscnd = vModel(CPLEX.Optimizer);  set_silent(smscnd)
    MOI.set(smscnd, MOI.RelativeGapTolerance(), tol)
    @variable(smscnd, y[1:(dt.J+dt.K)*2], Bin)
    @variable(smscnd, uij[1:sum(dt.Mij)], Bin)
    @variable(smscnd, ujk[1:sum(dt.Mjk)], Bin)
    @variable(smscnd, ukl[1:sum(dt.Mkl)], Bin)

    @variable( smscnd, 0<= xij[1:sum(dt.Mij)]);
    @variable( smscnd, 0<= xjk[1:sum(dt.Mjk)]);
    @variable( smscnd, 0<= xkl[1:sum(dt.Mkl)]);
    @variable( smscnd, 0<= h[1:(dt.J+dt.K)*2]);

    @addobjective(smscnd, Min, sum(dt.fc.*y)+dot(dt.gij,uij) + dot(dt.gjk,ujk) + dot(dt.gkl,ukl) +
        sum(repeat(dt.vcs[1,:], outer=sum(dt.Mij[1,:])).*xij[1:sum(dt.Mij[1,:])]) +
        sum(sum(repeat(dt.vcs[i,:], outer=sum(dt.Mij[i,:])).*xij[sum(dt.Mij[1:i-1,:])+1:sum(dt.Mij[1:i,:])]) for i=2:dt.I) +
        dot(dt.tcp,xij) + dot(dt.tcd,xjk) + dot(dt.tcc,xkl)+ dot([dt.vcp;dt.vcd],h)
    )

    @addobjective(smscnd, Min, sum(repeat(dt.ves[1,:], outer=sum(dt.Mij[1,:])).*xij[1:sum(dt.Mij[1,:])]) +
        sum(sum(repeat(dt.ves[i,:], outer=sum(dt.Mij[i,:])).*xij[sum(dt.Mij[1:i-1,:])+1:sum(dt.Mij[1:i,:])]) for i=2:dt.I) +
        dot(dt.cep,xij) + dot(dt.ced,xjk) + dot(dt.cec,xkl) + dot([dt.vep;dt.ved], h) 
    )
    
    ########## constraint 3 #############
    @constraint(smscnd, sum(xij[m] for m=1:dt.Mij[1,1])+ sum(xij[m+(sum(dt.Mij[1:i-1,:]))] for i=2:dt.I for m=1:dt.Mij[i,1])  == sum(xjk[m] for m=1:sum(dt.Mjk[1,:])) )
    @constraint(smscnd, [j=2:dt.J], sum(xij[sum(dt.Mij[1,1:j-1])+m] for m=1:dt.Mij[1,j])+sum(xij[sum(dt.Mij[i,1:j-1])+m+(sum(dt.Mij[1:i-1,:]))] for i=2:dt.I for m=1:dt.Mij[i,j]) == sum(xjk[sum(dt.Mjk[1:j-1,:]) + m] for m=1:sum(dt.Mjk[j,:])) )
    @constraint(smscnd, sum(xjk[m] for m=1:dt.Mjk[1,1])+sum(xjk[m+(sum(dt.Mjk[1:j-1,:]))] for j=2:dt.J for m=1:dt.Mjk[j,1]) == sum(xkl[m] for m=1:sum(dt.Mkl[1,:])) )
    @constraint(smscnd, [k=2:dt.K],sum(xjk[sum(dt.Mjk[1,1:k-1])+m] for m=1:dt.Mjk[1,k])+sum(xjk[sum(dt.Mjk[j,1:k-1])+m+(sum(dt.Mjk[1:j-1,:]))] for j=2:dt.J for m=1:dt.Mjk[j,k]) == sum(xkl[sum(dt.Mkl[1:k-1,:]) + m] for m=1:sum(dt.Mkl[k,:])) )
    ########### constraint 4-6 #############
    @constraint(smscnd, sum(h[t] for t=1:2) ==sum(xij[m] for m=1:dt.Mij[1,1])+sum(xij[m+(sum(dt.Mij[1:i-1,:]))] for i=2:dt.I for m=1:dt.Mij[i,1]));
    @constraint(smscnd, [j=2:dt.J], sum(h[2*(j-1)+t] for t=1:2) == sum(xij[sum(dt.Mij[1,1:j-1])+m] for m=1:dt.Mij[1,j])+sum(xij[sum(dt.Mij[i,1:j-1])+m+(sum(dt.Mij[1:i-1,:]))] for i=2:dt.I for m=1:dt.Mij[i,j]) );
    @constraint(smscnd,sum(h[2*dt.J+t] for t=1:2) == sum(xjk[m] for m=1:dt.Mjk[1,1])+sum(xjk[m+(sum(dt.Mjk[1:j-1,:]))] for j=2:dt.J for m=1:dt.Mjk[j,1]) );
    @constraint(smscnd, [k=2:dt.K], sum(h[2*dt.J+2*(k-1)+t] for t=1:2) == sum(xjk[sum(dt.Mjk[1,1:k-1])+m] for m=1:dt.Mjk[1,k])+sum(xjk[sum(dt.Mjk[j,1:k-1])+m+(sum(dt.Mjk[1:j-1,:]))] for j=2:dt.J for m=1:dt.Mjk[j,k]));
    @constraint(smscnd, sum(xkl[m] for m=1:dt.Mkl[1,1]) +sum(xkl[m+(sum(dt.Mkl[1:k-1,:]))] for k=2:dt.K for m=1:dt.Mkl[k,1]) >= dt.demand[1]);
    @constraint(smscnd, [l=2:dt.L], sum(xkl[sum(dt.Mkl[1,1:l-1]) + m] for m=1:dt.Mkl[1,l])+ sum(xkl[sum(dt.Mkl[1:k-1,:])+sum(dt.Mkl[k,1:l-1])+m] for k=2:dt.K for m=1:dt.Mkl[k,l]) >= dt.demand[l]);
    ########### constraint 7 #############
    @constraint(smscnd, sum(xij[1:sum(dt.Mij[1,:])]) <= dt.cas[1]);
    @constraint(smscnd, [i=2:dt.I],  sum(xij[sum(dt.Mij[1:i-1,:])+1:sum(dt.Mij[1:i,:])]) <= dt.cas[i]);
    ########### constraint 8 #############
    @constraint(smscnd,[j=1:dt.J+dt.K, t=1:2], sum(h[2*(j-1)+t]) <= [dt.cap;dt.cad][j]*y[2*(j-1)+t]);
    ########### constraint 9 #############
    @constraint(smscnd,[j=1:dt.J+dt.K], sum(y[2*(j-1)+1:2*(j-1)+2]) <= 1);
    # ########### constraint 10 #############
    @constraints(smscnd,begin
        sum(uij[1:dt.Mij[1,1]]) <= 1 #[1,1]
        [j=2:dt.J], sum(uij[sum(dt.Mij[1,1:j-1])+1:sum(dt.Mij[1,1:j-1])+dt.Mij[1,j]]) <= 1 # [1,2:6]
        [i=2:dt.I], sum(uij[sum(dt.Mij[1:i-1,:])+1:sum(dt.Mij[1:i-1,:])+dt.Mij[i,1]])<= 1 # [2:6,1]
        [i=2:dt.I,j=2:dt.J], sum(uij[sum(dt.Mij[1:i-1,:])+sum(dt.Mij[i,1:j-1])+1:sum(dt.Mij[1:i-1,:])+sum(dt.Mij[i,1:j-1])+dt.Mij[i,j]])<= 1 #sum(y[j,:])
        
        sum(ujk[1:dt.Mjk[1,1]]) <= 1 #(y[1,1]+y[dt.J+1,1])/2
        [k=2:dt.K], sum(ujk[sum(dt.Mjk[1,1:k-1])+1:sum(dt.Mjk[1,1:k-1])+dt.Mjk[1,k]]) <= 1 # (sum(y[1,:])+sum(y[dt.J+k,:]))/2
        [j=2:dt.J], sum(ujk[sum(dt.Mjk[1:j-1,:])+1:sum(dt.Mjk[1:j-1,:])+dt.Mjk[j,1]])<= 1 # [2:6,1]
        [j=2:dt.J,k=2:dt.K],  sum(ujk[sum(dt.Mjk[1:j-1,:])+sum(dt.Mjk[j,1:j-1])+1:sum(dt.Mjk[1:j-1,:])+sum(dt.Mjk[j,1:j-1])+dt.Mjk[j,k]]) <= 1 #(sum(y[j,:])+sum(y[dt.J+k,:]))/2

        sum(ukl[1:dt.Mkl[1,1]]) <= 1 #y[dt.J+1,1] #[1,1]
        [l=2:dt.L], sum(ukl[sum(dt.Mkl[1,1:l-1])+1:sum(dt.Mkl[1,1:l-1])+dt.Mkl[1,l]]) <= 1 #sum(y[dt.J+1,:])
        [k=2:dt.K], sum(ukl[sum(dt.Mkl[1:k-1,:])+1:sum(dt.Mkl[1:k-1,:])+dt.Mkl[k,1]])<= 1 
        [k=2:dt.K,l=2:dt.L],  sum(ukl[sum(dt.Mkl[1:k-1,:])+sum(dt.Mkl[k,1:l-1])+1:sum(dt.Mkl[1:k-1,:])+sum(dt.Mkl[k,1:l-1])+dt.Mkl[k,l]])<= 1 #sum(y[dt.J+k,:])
    end);
    # ########### constraint 11 #############
    @constraints(smscnd, begin
            [i=1:dt.I], xij[i] >= dt.Lcapasp[i]*uij[i]
            [j=1:dt.J], xjk[j] >= dt.Lcapapd[j]*ujk[j]
    end);
    ########### constraint 12 #############
    @constraints(smscnd, begin
        [i=1:sum(dt.Mij)], xij[i] <= dt.bigM*uij[i]
        [j=1:sum(dt.Mjk)], xjk[j] <= dt.bigM*ujk[j]
        [k=1:sum(dt.Mkl)], xkl[k] <= dt.bigM*ukl[k]
    end);
    ########### constraint 13-14 #############
    @constraint(smscnd, sum(y[1:dt.J*2]) <= dt.Jmax);
    @constraint(smscnd, sum(y[dt.J*2+1:end]) <= dt.Kmax);
    return smscnd
end

global rhs1 = append!([dt.vcs[i]+dt.tcp[m] for i=1 for j=1 for m=1:dt.Mij[i,j]],[dt.vcs[i]+dt.tcp[sum(dt.Mij[i,1:j-1])+m] for i=1 for j=2:dt.J for m=1:dt.Mij[i,j]],[dt.vcs[i]+dt.tcp[sum(dt.Mij[1:i-1,:])+m] for i=2:dt.I for j=1 for m=1:dt.Mij[i,j]],[dt.vcs[i]+dt.tcp[sum(dt.Mij[1:i-1,:])+sum(dt.Mij[i,1:j-1])+m] for i=2:dt.I for j=2:dt.J for m=1:dt.Mij[i,j]], dt.tcd,dt.tcc,dt.vcp,dt.vcd);
global rhs2 = append!([dt.ves[i]+dt.cep[m] for i=1 for j=1 for m=1:dt.Mij[i,j]],[dt.ves[i]+dt.cep[sum(dt.Mij[i,1:j-1])+m] for i=1 for j=2:dt.J for m=1:dt.Mij[i,j]],[dt.ves[i]+dt.cep[sum(dt.Mij[1:i-1,:])+m] for i=2:dt.I for j=1 for m=1:dt.Mij[i,j]],[dt.ves[i]+dt.cep[sum(dt.Mij[1:i-1,:])+sum(dt.Mij[i,1:j-1])+m] for i=2:dt.I for j=2:dt.J for m=1:dt.Mij[i,j]], dt.ced,dt.cec,dt.vep,dt.ved);
global objc = append!(-[dt.cap[j] for j=1:dt.J for t=1:2], -[dt.cad[k] for k=1:dt.K for t=1:2],dt.Lcapasp, dt.Lcapapd, fill(-dt.bigM, sum(dt.Mij)+sum(dt.Mjk)+sum(dt.Mkl) ))

dot(dff.cuts[f].α5,dt.demand)-dot(dt.cas,dff.cuts[f].α6)-sum(dt.cap[j]*vy[2*(j-1)+t]*dff.cuts[f].α7[2*(j-1)+t] for j=1:dt.J for t=1:2)-
sum(dt.cad[k]*vy[2*dt.J+k]*dff.cuts[f].α8[2*(k-1)+t] for k=1:dt.K for t=1:2)+
dot(dt.Lcapasp.*vuij,dff.cuts[f].α9)+dot(dt.Lcapapd.*vujk,dff.cuts[f].α10) -
dot(dt.bigM.*vuij,dff.cuts[f].α12)-dot(dt.bigM.*vujk,dff.cuts[f].α13) -dot(dt.bigM.*vukl,dff.cuts[f].α14)

dot(dff.cuts[f].α5,dt.demand)-dot(dt.cas,dff.cuts[f].α6)-
sum(dt.cap[j]*mas.y[2*(j-1)+t]*dff.cuts[f].α7[2*(j-1)+t] for j=1:dt.J for t=1:2)-
sum(dt.cad[k]*mas.y[2*dt.J+k]*dff.cuts[f].α8[2*(k-1)+t] for k=1:dt.K for t=1:2)+
dot(dt.Lcapasp.*mas.uij,dff.cuts[f].α9)+dot(dt.Lcapapd.*mas.ujk,dff.cuts[f].α10) -
dot(dt.bigM.*mas.uij,dff.cuts[f].α12)-dot(dt.bigM.*mas.ujk,dff.cuts[f].α13) -dot(dt.bigM.*mas.ukl,dff.cuts[f].α14) <= 0
