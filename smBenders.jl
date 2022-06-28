using StatsBase,DelimitedFiles,Distributions,CPLEX,JuMP,LinearAlgebra,CPUTime
#######################  small instance Benders Decomposition  ########################
struct build_master
    y::Matrix{VariableRef}
    uij::JuMP.Containers.SparseAxisArray{VariableRef}
    ujk::JuMP.Containers.SparseAxisArray{VariableRef}
    ukl::JuMP.Containers.SparseAxisArray{VariableRef}
    m::Model
end
function build_master()
    mas = Model(CPLEX.Optimizer);
    set_silent(mas); MOI.set(mas, MOI.NumberOfThreads(), 1);
    @variable(mas, y[1:J+K,1:2], Bin);
    # @variable(mas, ykt[1:K,1:2], Bin);
    @variable(mas, uij[i=1:I,j=1:J,1:Mij[i,j]], Bin);
    @variable(mas, ujk[j=1:J,k=1:K,1:Mjk[j,k]], Bin);
    @variable(mas, ukl[k=1:K,l=1:L,1:Mkl[k,l]], Bin);
    @constraint(mas,[j=1:J+K], sum(y[j,:]) <= 1);
    # @constraint(mas,[k=1:K], -sum(ykt[k,:]) >= -1);
    @constraint(mas,[i=1:I,j=1:J], sum(uij[i,j,m] for m=1:Mij[i,j]) <= 1);
    @constraint(mas,[j=1:J,k=1:K], sum(ujk[j,k,m] for m=1:Mjk[j,k]) <= 1);
    @constraint(mas,[k=1:K,l=1:L], sum(ukl[k,l,m] for m=1:Mkl[k,l]) <=1);
    @constraint(mas, sum(y[J+k,t] for k=1:K for t=1:2) <= Kmax);
    # optimize!(mas)
    return build_master(y,uij,ujk,ukl,mas)
end

struct DualSubP
    # data::SubProblemData
    α1::Vector{VariableRef}
    α2::Vector{VariableRef}
    α3::Vector{VariableRef}
    α4::Vector{VariableRef}
    α5::Vector{VariableRef}
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
    set_optimizer_attribute(sub, "CPX_PARAM_REDUCE", 0); # presolve must be turned off
    set_silent(sub)
    MOI.set(sub, MOI.NumberOfThreads(), 1)

    @variables(sub, begin
        α1[1:J]
        α2[1:K]
        α3[1:J]
        α4[1:K]
    end)
    @variable(sub, α5[1:L] >= 0);
    @variable(sub, α6[1:I] >= 0);
    @variable(sub, α7[1:J,1:2] >= 0);
    @variable(sub, α8[1:K,1:2] >= 0);
    @variable(sub, α9[i=1:I,j=1:J,1:Mij[i,j]] >= 0);
    @variable(sub, α10[j=1:J,k=1:K,1:Mjk[j,k]] >= 0);
    @variable(sub, α11[k=1:K,l=1:L,1:Mkl[k,l]] >= 0);
    @variable(sub, α12[i=1:I,j=1:J,1:Mij[i,j]] >= 0);
    @variable(sub, α13[j=1:J,k=1:K,1:Mjk[j,k]] >= 0);
    @variable(sub, α14[k=1:K,l=1:L,1:Mkl[k,l]] >= 0);

    @constraint(sub, [i=1:I,j=1:J,m=1:Mij[i,j]], α1[j]-α3[j]-α6[i]+α9[i,j,m]-α12[i,j,m] <= w[1]*(vcs[i]+tcp[i][j][m])+w[2]*(ves[i]+cep[i][j][m]));
    @constraint(sub, [j=1:J,k=1:K,m=1:Mjk[j,k]], -α1[j]+α2[k]-α4[k]+α10[j,k,m]-α13[j,k,m] <= w[1]*tcd[j][k][m]+w[2]*ced[j][k][m]);
    @constraint(sub, [k=1:K,l=1:L,m=1:Mkl[k,l]], -α2[k]+α5[l]+α11[k,l,m]-α14[k,l,m] <= w[1]*tcc[k][l][m]+w[2]*cec[k][l][m]);
    @constraint(sub, [j=1:J,t=1:2], α3[j]-α7[j,t] <= w[1]*vcp[j,t]+w[2]*vep[j,t]);
    @constraint(sub, [k=1:K,t=1:2], α4[k]-α8[k,t] <= w[1]*vcd[k,t]+w[2]*ved[k,t]);
    return DualSubP(α1,α2,α3,α4,α5,α6,α7,α8,α9,α10,α11,α12,α13,α14,sub)
end

function solve_dsp(dsp::DualSubP,yb,ubij,ubjk,ubkl)
    @objective(dsp.m, Max, dot(demand,dsp.α5)-dot(cas,dsp.α6)-
        sum(cap[j]*dsp.α7[j,t]*yb[j,t] for j=1:J for t=1:2) -
        sum(cad[k]*dsp.α8[k,t]*yb[J+k,t] for k=1:K for t=1:2) +
        sum(Lcapasp[i][j][m]*dsp.α9[i,j,m]*ubij[i,j,m] for i=1:I for j=1:J for m=1:Mij[i,j]) +
        sum(Lcapapd[j][k][m]*dsp.α10[j,k,m]*ubjk[j,k,m] for j=1:J for k=1:K for m=1:Mjk[j,k]) +
        sum(Lcapadc[k][l][m]*dsp.α11[k,l,m]*ubkl[k,l,m] for k=1:K for l=1:L for m=1:Mkl[k,l])-
        sum(bigM*ubij[i,j,m]*dsp.α12[i,j,m] for i=1:I for j=1:J for m=1:dt.Mij[i,j]) -
        sum(bigM*ubjk[j,k,m]*dsp.α13[j,k,m] for j=1:J for k=1:K for m=1:dt.Mjk[j,k]) -
        sum(bigM*ubkl[k,l,m]*dsp.α14[k,l,m] for k=1:K for l=1:L for m=1:dt.Mkl[k,l])
        );
    optimize!(dsp.m)
    st = termination_status(dsp.m)

    if st == MOI.OPTIMAL
      return (res = :OptimalityCut, obj = objective_value(dsp.m), α5 = value.(dsp.α5), α6 = value.(dsp.α6), α7 = value.(dsp.α7), α8= value.(dsp.α8), α9= value.(dsp.α9), α10= value.(dsp.α10), α11= value.(dsp.α11), α12= value.(dsp.α12),α13= value.(dsp.α13),α14= value.(dsp.α14) )
    elseif st == MOI.DUAL_INFEASIBLE
      # println( value.(dsp.α5),  value.(dsp.α6), value.(dsp.α7), value.(dsp.α8),  value.(dsp.α9),  value.(dsp.α10), value.(dsp.α11), value.(dsp.α12),value.(dsp.α13),value.(dsp.α14) )
      return ( res = :FeasibilityCut, α5 = value.(dsp.α5), α6 = value.(dsp.α6), α7 = value.(dsp.α7), α8= value.(dsp.α8), α9= value.(dsp.α9), α10= value.(dsp.α10), α11= value.(dsp.α11), α12= value.(dsp.α12),α13= value.(dsp.α13),α14= value.(dsp.α14) )

    else
      error("DualSubProblem error: status $st")
    end
end
function benders_decomposition(w,mas::Model, y::Matrix{VariableRef},uij::JuMP.Containers.SparseAxisArray{VariableRef},ujk::JuMP.Containers.SparseAxisArray{VariableRef},ukl::JuMP.Containers.SparseAxisArray{VariableRef})
    sub = direct_model(CPLEX.Optimizer()); dsp = DualSubP(sub,w)

    @variable(mas, θ>= -1000);
    @objective(mas, Min, w[1]*(sum(fc[j][t]*y[j,t] for j=1:J+K for t=1:2)+
        sum(gij[i][j][m]*uij[i,j,m] for i=1:I for j=1:J for m=1:Mij[i,j])+
        sum(gjk[j][k][m]*ujk[j,k,m] for j=1:J for k=1:K for m=1:Mjk[j,k])+
        sum(gkl[k][l][m]*ukl[k,l,m] for k=1:K for l=1:L for m=1:Mkl[k,l])) + θ);
    optimize!(mas); st = termination_status(mas)
    nopt_cons, nfeasi_cons = (0, 0)
    feasicuts = [];

    while (st == MOI.INFEASIBLE) || (st == MOI.OPTIMAL)
        optimize!(mas);
        st = termination_status(mas)
        θ1 = value(θ); yb = value.(y); ubij = value.(uij); ubjk = value.(ujk); ubkl = value.(ukl)
        subp = solve_dsp(dsp,yb,ubij,ubjk,ubkl)
        if subp.res == :OptimalityCut
            # @info "Optimality cut found"
            if round(θ1; digits=4) ≥ round(subp.obj; digits=4)
                break
            else
                nopt_cons+=1
                cut = @constraint( mas, θ ≥ dot(demand,subp.α5)-dot(cas,subp.α6)-
                    sum(cap[j]*subp.α7[j,t]*y[j,t] for j=1:J for t=1:2) -
                    sum(cad[k]*subp.α8[k,t]*y[k+J,t] for k=1:K for t=1:2) +
                    sum(Lcapasp[i][j][m]*subp.α9[i,j,m]*uij[i,j,m] for i=1:I for j=1:J for m=1:Mij[i,j]) +
                    sum(Lcapapd[j][k][m]*subp.α10[j,k,m]*ujk[j,k,m] for j=1:J for k=1:K for m=1:Mjk[j,k]) +
                    sum(Lcapadc[k][l][m]*subp.α11[k,l,m]*ukl[k,l,m] for k=1:K for l=1:L for m=1:Mkl[k,l])-
                    sum(bigM*uij[i,j,m]*subp.α12[i,j,m] for i=1:I for j=1:J for m=1:dt.Mij[i,j]) -
                    sum(bigM*ujk[j,k,m]*subp.α13[j,k,m] for j=1:J for k=1:K for m=1:dt.Mjk[j,k]) -
                    sum(bigM*ukl[k,l,m]*subp.α14[k,l,m] for k=1:K for l=1:L for m=1:dt.Mkl[k,l])
                    )
                    # build_constraint
                    # MOI.submit(mas, MOI.LazyConstraint(), cut)
            end
        else
            # @info "Feasibility cut found"
            nfeasi_cons += 1
            cut = @constraint( mas, 0 ≥ dot(demand,subp.α5)-dot(cas,subp.α6)-
                sum(cap[j]*subp.α7[j,t]*y[j,t] for j=1:J for t=1:2) -
                sum(cad[k]*subp.α8[k,t]*y[k+J,t] for k=1:K for t=1:2) +
                sum(Lcapasp[i][j][m]*subp.α9[i,j,m]*uij[i,j,m] for i=1:I for j=1:J for m=1:Mij[i,j]) +
                sum(Lcapapd[j][k][m]*subp.α10[j,k,m]*ujk[j,k,m] for j=1:J for k=1:K for m=1:Mjk[j,k]) +
                sum(Lcapadc[k][l][m]*subp.α11[k,l,m]*ukl[k,l,m] for k=1:K for l=1:L for m=1:Mkl[k,l])-
                sum(bigM*uij[i,j,m]*subp.α12[i,j,m] for i=1:I for j=1:J for m=1:dt.Mij[i,j]) -
                sum(bigM*ujk[j,k,m]*subp.α13[j,k,m] for j=1:J for k=1:K for m=1:dt.Mjk[j,k]) -
                sum(bigM*ukl[k,l,m]*subp.α14[k,l,m] for k=1:K for l=1:L for m=1:dt.Mkl[k,l]))
            # push!(feasicuts, cut)
            push!(feasicuts, (α5=subp.α5,α6=subp.α6,α7=subp.α7,α8=subp.α8,α9=subp.α9,α10=subp.α10,α11=subp.α11,α12=subp.α12,α13=subp.α13,α14=subp.α14))
        end
        # @info "Adding the cut $(cut)"
    end
    return (mas, y, uij, ujk, ukl, nopt_cons, nfeasi_cons, feasicuts)
end


w = [0.5,0.5];
w1buildTime = @CPUelapsed mas = build_master();
w1algoTime = @CPUelapsed fmodel,fy,fuij,fujk,fukl,fnoptcut,fnfeasicut,firstcuts = benders_decomposition(w,mas.m,mas.y,mas.uij,mas.ujk,mas.ukl);
(fnoptcut+fnfeasicut)
value.(fujk)


w=[0,1];

w2buildTime = @CPUelapsed mas = build_master();
w2algoTime = @CPUelapsed  smodel,sy,suij,sujk,sukl,snoptcut,snfeasicut,secondcuts = benders_decomposition(w,mas.m,mas.y,mas.uij,mas.ujk,mas.ukl);
(snoptcut+snfeasicut)


function build_master2(m::Model,firstcuts,secondcuts)
    mp = build_master()
    for i=1:length(firstcuts)
        @constraint(mp.m,  0 ≥ dot(demand,firstcuts[i].α5)-dot(cas,firstcuts[i].α6)-
            sum(cap[j]*firstcuts[i].α7[j,t]*mp.y[j,t] for j=1:J for t=1:2) -
            sum(cad[k]*firstcuts[i].α8[k,t]*mp.y[k+J,t] for k=1:K for t=1:2) +
            sum(Lcapasp[i][j][m]*firstcuts[i].α9[i,j,m]*mp.uij[i,j,m] for i=1:I for j=1:J for m=1:Mij[i,j]) +
            sum(Lcapapd[j][k][m]*firstcuts[i].α10[j,k,m]*mp.ujk[j,k,m] for j=1:J for k=1:K for m=1:Mjk[j,k]) +
            sum(Lcapadc[k][l][m]*firstcuts[i].α11[k,l,m]*mp.ukl[k,l,m] for k=1:K for l=1:L for m=1:Mkl[k,l])-
            sum(bigM*mp.uij[i,j,m]*firstcuts[i].α12[i,j,m] for i=1:I for j=1:J for m=1:dt.Mij[i,j]) -
            sum(bigM*mp.ujk[j,k,m]*firstcuts[i].α13[j,k,m] for j=1:J for k=1:K for m=1:dt.Mjk[j,k]) -
            sum(bigM*mp.ukl[k,l,m]*firstcuts[i].α14[k,l,m] for k=1:K for l=1:L for m=1:dt.Mkl[k,l]))
    end
    for i=1:length(secondcuts)
        @constraint(mp.m,  0 ≥ dot(demand,secondcuts[i].α5)-dot(cas,secondcuts[i].α6)-
            sum(cap[j]*secondcuts[i].α7[j,t]*mp.y[j,t] for j=1:J for t=1:2) -
            sum(cad[k]*secondcuts[i].α8[k,t]*mp.y[k+J,t] for k=1:K for t=1:2) +
            sum(Lcapasp[i][j][m]*secondcuts[i].α9[i,j,m]*mp.uij[i,j,m] for i=1:I for j=1:J for m=1:Mij[i,j]) +
            sum(Lcapapd[j][k][m]*secondcuts[i].α10[j,k,m]*mp.ujk[j,k,m] for j=1:J for k=1:K for m=1:Mjk[j,k]) +
            sum(Lcapadc[k][l][m]*secondcuts[i].α11[k,l,m]*mp.ukl[k,l,m] for k=1:K for l=1:L for m=1:Mkl[k,l])-
            sum(bigM*mp.uij[i,j,m]*secondcuts[i].α12[i,j,m] for i=1:I for j=1:J for m=1:dt.Mij[i,j]) -
            sum(bigM*mp.ujk[j,k,m]*secondcuts[i].α13[j,k,m] for j=1:J for k=1:K for m=1:dt.Mjk[j,k]) -
            sum(bigM*mp.ukl[k,l,m]*secondcuts[i].α14[k,l,m] for k=1:K for l=1:L for m=1:dt.Mkl[k,l]))
    end

    return mp
end
build2mp = @CPUelapsed mp2 = build_master2(mp.m,firstcuts,secondcuts);
algoTime = @CPUelapsed nmodel,ny,nuij,nujk,nukl,noptcut,nfeasicut,ncuts =benders_decomposition([0.5,0.5],mp2.m,mp2.y,mp2.uij,mp2.ujk,mp2.ukl);


w1buildTime+w1algoTime+w2buildTime+w2algoTime+build2mp+algoTime

basicbuildTime = @CPUelapsed bm = build_master();
basicalgoTime = @CPUelapsed bmodel,by,buij,bujk,bukl,bnoptcut,bnfeasicut,bcuts = benders_decomposition([0.5,0.5],bm.m,bm.y,bm.uij,bm.ujk,bm.ukl);
basicbuildTime+basicalgoTime


function AcceleratedBD()
    w = [1,0];
    mas = build_master();
    fmodel,fy,fuij,fujk,fukl,fnoptcut,fnfeasicut,firstcuts = benders_decomposition(w,mas.m,mas.y,mas.uij,mas.ujk,mas.ukl);
    w=[0,1];
    mas = build_master();
    smodel,sy,suij,sujk,sukl,snoptcut,snfeasicut,secondcuts = benders_decomposition(w,mas.m,mas.y,mas.uij,mas.ujk,mas.ukl);
    mp2 = build_master2(mp.m,firstcuts,secondcuts);
    w = [[0.2,0.8],[0.4,0.6],[0.5,0.5],[0.6,0.4],[0.8,0.2]]
    yvals = []
    for i=1:length(w)
        nmodel,ny,nuij,nujk,nukl,noptcut,nfeasicut,ncuts =benders_decomposition(w[i],mp2.m,mp2.y,mp2.uij,mp2.ujk,mp2.ukl)
        push!(yvals,ny)
    end
    return yvals
end


function BasicBD()
    w = [[0.2,0.8],[0.4,0.6],[0.5,0.5],[0.6,0.4],[0.8,0.2]]
    yvals = []; uijvals = []; ujkvals = []; uklvals = []
    for i=1:length(w)
        mas = build_master();
        model,y,uij,ujk,ukl,noptcut,nfeasicut,feasicuts = benders_decomposition(w[i],mas.m,mas.y,mas.uij,mas.ujk,mas.ukl);
        push!(yvals,y); push!(uijvals,uij); push!(ujkvals,ujk); push!(uklvals,ukl)
    end
    return yvals,uijvals,ujkvals,uklvals
end

@CPUelapsed byvals,buijvals,bujkvals,buklvals = BasicBD()
value.(byvals[4])
value.(buklvals[4])==value.(buklvals[2])
unique(buijvals)
1
# unionset = union(firstcuts,secondcuts);
# for i=1:length(unionset)
#     @constraint( mp.m, unionset[i] )
# end

value.(hy)
value.(huij)
value.(hukl)
hcuts
1
# sub = Model(CPLEX.Optimizer);
# set_optimizer_attribute(sub, "CPX_PARAM_REDUCE", 0)
# MOI.set(sub, MOI.NumberOfThreads(), 1);set_silent(sub)
# @variable(sub, 0<= xij[i=1:I,j=1:J,1:Mij[i,j]]);
# @variable(sub, 0<= xjk[j=1:J,k=1:K,1:Mjk[j,k]]);
# @variable(sub, 0<= xkl[k=1:K,l=1:L,1:Mkl[k,l]]);
# @variable(sub, 0<= h[1:J+K,1:2]);
# # @variable(sub, yb[1:J+K,1:2], Bin);
# # @variable(sub, ubij[1:I,1:J,1:2], Bin);
# # @variable(sub, ubjk[1:J,1:K,1:2], Bin);
# # @variable(sub, ubkl[1:K,1:L,1:2], Bin);
#
# ########### constraint 3 #############
# @constraints(sub, begin
#     [j=1:J], sum(xij[i,j,m] for i=1:I for m=1:Mij[i,j]) == sum(xjk[j,k,m] for k=1:K for m=1:Mjk[j,k])
#     [k=1:K], sum(xjk[j,k,m] for j=1:J for m=1:Mjk[j,k]) == sum(xkl[k,l,m] for l=1:L for m=1:Mkl[k,l])
# end)
# con3 = @constraint(sub, [j=1:J], sum(h[j,:]) == sum(xij[i,j,m] for i=1:I for m=1:Mij[i,j]))
# con4 = @constraint(sub, [k=1:K], sum(h[k+J,:] ) == sum(xjk[j,k,m] for j=1:J for m=1:Mjk[j,k]))
# con5 = @constraint(sub, [l=1:L], sum(xkl[k,l,m] for k=1:K for m=1:Mkl[k,l]) >= demand[l])
# con6 = @constraint(sub,[i=1:I], -sum(xij[i,j,m] for j=1:J for m=1:Mij[i,j] ) >= -cas[i]);
# # con7 = @constraint(sub,[j=1:J, t=1:2], -sum(h[j,t]) >= -cap[j]*yb[j,t]);
# # con8 = @constraint(sub,[k=1:K, t=1:2], -sum(h[k+J,t]) >= -cad[k]*yb[k+J,t]);
# # con9 = @constraint(sub,[i=1:I, j=1:J, m=1:Mij[i,j]], sum(xij[i,j,m]) >= Lcapasp[i][j][m]*ubij[i,j,m] )
# # con10 = @constraint(sub,[j=1:J, k=1:K, m=1:Mjk[j,k]], sum(xjk[j,k,m]) >= Lcapapd[j][k][m]*ubjk[j,k,m]);
# # con11 = @constraint(sub,[k=1:K, l=1:L, m=1:Mkl[k,l]], sum(xkl[k,l,m]) >= Lcapadc[k][l][m]*ubkl[k,l,m]);
#
# @objective( sub, Min, w[2]*( sum(sum(vcs[i]*xij[i,:,:] for i=1:I)[k] for k=1:I) + sum([vcp;vcd][j,t]*h[j,t] for j=1:J+K for t=1:2) +
#     sum(tcp[i][j][m]*xij[i,j,m] for i=1:I for j=1:J for m=1:Mij[i,j])+sum(tcd[j][k][m]*xjk[j,k,m] for j=1:J for k=1:K for m=1:Mjk[j,k])+
#     sum(tcc[k][l][m]*xkl[k,l,m] for k=1:K for l=1:L for m=1:Mkl[k,l])+
#     sum(ves[i]*xij[i,j,m] for i=1:I for j=1:J for m=1:Mij[i,j])+
#     sum([vep;ved][j,t]*h[j,t] for j=1:J+K for t=1:2) +
#     sum(cep[i][j][m]*xij[i,j,m] for i=1:I for j=1:J for m=1:Mij[i,j]) +
#     sum(ced[j][k][m]*xjk[j,k,m] for j=1:J for k=1:K for m=1:Mjk[j,k]) +
#     sum(cec[k][l][m]*xkl[k,l,m] for k=1:K for l=1:L for m=1:Mkl[k,l]))
# );
# function Defexp(xij,xjk,xkl)
#     exa = AffExpr(0);
#     for i=1:I
#         add_to_expression!(exa,sum(vcs[i]*xij[i,:,:]));
#     end
#     #v_ijmp*x_ijmp expression
#     exv = AffExpr(0);
#     for i=1:I
#         idx = 1;
#         for j=1:J
#             for m=1:Mij[i,j]
#                 add_to_expression!(exv,sum(dot.(tcp[i][idx],xij[i,j,m])))
#                 idx+=1
#             end
#         end
#     end
#     for j=1:J
#         idx = 1;
#         for k=1:K
#             for m=1:Mjk[j,k]
#                 add_to_expression!(exv,sum(dot.(tcd[j][idx],xjk[j,k,m])))
#                 idx+=1
#             end
#         end
#     end
#     for k=1:K
#         idx = 1;
#         for l=1:L
#             for m=1:Mkl[k,l]
#                 add_to_expression!(exv,sum(dot.(tcc[k][idx],xkl[k,l,m])))
#                 idx+=1
#             end
#         end
#     end
#     #b_ip*x_ijmp
#     exb = AffExpr(0);
#     for i=1:I
#         for j=1:J
#             for m=1:Mij[i,j]
#                 add_to_expression!(exb,sum(ves[i]*xij[i,j,m]) );
#             end
#         end
#     end
#     exr = AffExpr(0);
#     for i=1:I
#         idx = 1;
#         for j=1:J
#             for m=1:Mij[i,j]
#                 add_to_expression!(exr,sum(cep[i][j][m]*xij[i,j,m]))
#                 idx+=1
#             end
#         end
#     end
#     for j=1:J
#         idx = 1;
#         for k=1:K
#             for m=1:Mjk[j,k]
#                 add_to_expression!(exr,sum(dot.(ced[j][idx],xjk[j,k,m])))
#                 idx+=1
#             end
#         end
#     end
#     for k=1:K
#         idx = 1;
#         for l=1:L
#             for m=1:Mkl[k,l]
#                 add_to_expression!(exr,sum(dot.(cec[k][idx],xkl[k,l,m])))
#                 idx+=1
#             end
#         end
#     end
#     return exa,exv,exb,exr
# end

################## solve_dsp Function
# rays = Vector{Float64}(undef,length(dsp.α5)+length(dsp.α6)+length(dsp.α7)+length(dsp.α8)+length(dsp.α9)+length(dsp.α10)+length(dsp.α11))
# CPXgetray(backend(dsp.m).env, backend(dsp.m).lp, rays)
# lens = [length(dsp.α5),length(dsp.α6),length(dsp.α7),length(dsp.α8),length(dsp.α9),length(dsp.α10),length(dsp.α11)]
# vr = value.(rays)
# lens = [length(dsp.α5),length(dsp.α6),length(dsp.α7),length(dsp.α8),length(dsp.α9),length(dsp.α10),length(dsp.α11)]
# pvci = [sum(lens[1:l]) for l=1:length(lens)]
# insert!(pvci,1,0)
# a5,a6,a7,a8,a9,a10,a11 = [vr[pvci[l]+1:pvci[l+1]] for l=1:length(pvci)-1]
# return ( res = :FeasibilityCut, α5=reshape(a5,size(dsp.α5)), α6=reshape(a6,size(dsp.α6)), α7=reshape(a7,size(dsp.α7)), α8=reshape(a8,size(dsp.α8)),
#   α9=reshape(a9,size(dsp.α9)), α10=reshape(a10,size(dsp.α10)), α11=reshape(a11,size(dsp.α11)) )
