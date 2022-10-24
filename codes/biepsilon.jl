using DataStructures,DataFrames,DelimitedFiles,JuMP,CPLEX,LinearAlgebra,StatsBase
#,CSV,CPUTime,JLD2,,MathProgBase,MathOptInterface
# const MPB = MathProgBase;
# mutable struct importLP
#     lpfile::String; m::Int; n::Int; C::Array{}; B::Array{}; RHS::Dict{}; signs::Array{}; vub::Array{}
#     function importLP(lpfile::String)
#         lpmodel=buildlp([-1,0],[2 1],'<',1.5, CplexSolver(CPX_PARAM_SCRIND=0))
#         # lpmodel = CPLEX.CplexMathProgModel();
#         MPB.loadproblem!(lpmodel,lpfile)
#         Bmtx = MPB.getconstrmatrix(lpmodel);B = Bmtx[3:end,:]
#         C = Bmtx[1:2,:]
#         m,n=size(B)
#         vub = MPB.getvarUB(lpmodel)
#         lb = MPB.getconstrLB(lpmodel)[3:end]; ub = MPB.getconstrUB(lpmodel)[3:end]
#         RHS = Dict()
#         for i=1:m
#             if ub[i]==Inf
#                 RHS[i] = lb[i]
#             else
#                 RHS[i] = ub[i]
#             end
#         end
#         signs = []
#         for i=1:m
#             if ub[i] == Inf
#                 push!(signs,"l")
#             elseif lb[i] == -Inf
#                 push!(signs,"u")
#             else
#                 push!(signs, "s")
#             end
#         end
#         new(lpfile,m,n,C,B,RHS,signs,vub)
#     end
# end
# file = "/home/ak121396/Desktop/instances/SCND/test4s4.lp"
# file = "./lp/test4s4.lp"
# dtt = importLP(file)
# dtt = importLP("E:/scnd/Test4S3.lp")
# dtt = importLP("/home/k2g00/k2g3475/scnd/lp/"*file[36:end]*".lp")
bvar = findall(i -> i == 1, dtt.vub)
m1 = Model(CPLEX.Optimizer);
MOI.set(m1, MOI.RawParameter("CPX_PARAM_SCRIND"), false )
# MOI.set(m1, MOI.RawParameter("CPX_PARAM_THREADS"),1  )
@variable(m1, ex[1:dtt.n] >= 0);
for i = 1:dtt.n
    if i in bvar
        set_binary(ex[i])
    end
end
@variable(m1, 0 <= ep);
@constraint(m1, epcon1, dot(ex,dtt.C[2,:]) <= ep);
for k=1:dtt.m
    if dtt.signs[k] == "l"
        @constraint(m1, dot(dtt.B[k,:],ex) >= dtt.RHS[k])
    elseif dtt.signs[k] == "u"
        @constraint(m1, dot(dtt.B[k,:],ex) <= dtt.RHS[k])
    else
        @constraint(m1, dot(dtt.B[k,:],ex) == dtt.RHS[k])
    end
end
@objective(m1, Min, dot(ex,dtt.C[1,:]) );

function opt(ϵ)
    JuMP.fix(ep, ϵ; force = true);
    optimize!(m1)
    sol = JuMP.value.(ex)
    for i in bvar
        sol[i] = ceil(sol[i])
    end
    if termination_status(m1) == MOI.OPTIMAL
        return [dot(sol,dtt.C[1,:]),dot(sol,dtt.C[2,:])] #[JuMP.value.(ex),objective_value(m1)]
    else
        return nothing#,nothing
    end
end
function epsilon()
    # P = [];
    Y = []; ϵ = 19*10^(5); δ =10^(5); lb = 10*10^(5); fval = [0,ϵ]
    # Y = []; ϵ = 370*10^(6); δ =10^(7); lb = 235*10^(5); fval = [0,ϵ]

    # Test1S
    # ϵ = 5.5*10^(4); δ =10^(5); lb = 10^(6)+50; fval = [0,ϵ]
    while fval[2] >= lb
        # s,fval = opt(ϵ,C)
        fval = opt(ϵ)
        println(fval)
        if fval == nothing
            break
        # end
        # if dominated(fval,Y)==false
        else
            # push!(P,s);
            push!(Y,fval);
        end
        ϵ = ϵ-δ
    end
    return Y #P,
end
ey = epsilon() #ex

################# with built model
using DataStructures,DataFrames,DelimitedFiles,JuMP,CPLEX,SparseArrays,LinearAlgebra,StatsBase
file = "/home/k2g00/k2g3475/scnd/instances/test01S2"
# file = "F:/scnd/Test1S2"
struct Data1
    file::String; N::Dict{}; d::Array{}; c::Array{}; a::Array{}; e::Array{}; gij::SparseVector{}; gjk::SparseVector{}; gkl::SparseVector{};
    vij::Array{}; vjk::Array{}; vkl::Array{}; Vij::SparseVector{}; Vjk::SparseVector{}; Mij::Array{}; Mjk::Array{}; Mkl::Array{};
    b::Array{}; q::Array{}; rij::Array{}; rjk::Array{}; rkl::Array{}; upl::Int; udc::Int; bigM::Int
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
        d = reshape(N["demand"],5,N["customer"])'; c = append!(N["fcp"],N["fcd"]);
        a = reshape(N["vcs"],5,N["supplier"])';    e = append!(N["vcp"],N["vcd"]);
        gij = sparse(N["fixedcostModesp"]); gjk = sparse(N["fixedcostModepd"]); gkl = sparse(N["fixedcostModedc"]);
        vij = N["tcp"]; vjk = N["tcd"]; vkl = N["tcc"];
        Vij = sparse(N["LcapacityModesp"]); Vjk = sparse(N["LcapacityModepd"]); #Vkl =  sparse(N["LcapacityModedc"]);
        Mij = transpose(reshape(N["ModeIJ"], (N["plant"],N["supplier"])));
        Mjk = transpose(reshape(N["ModeJK"], (N["distribution"],N["plant"])));
        Mkl = transpose(reshape(N["ModeKL"], (N["customer"],N["distribution"])));
        b = reshape(N["ves"],N["supplier"],5);  q = append!(N["vep"],N["ved"]);
        rij = N["cep"]; rjk = N["ced"]; rkl = N["cec"];
        upl = N["upperpants"]; udc = N["upperdistribution"]; bigM = sum(N["demand"])

        new(file,N,d,c,a,e,gij,gjk,gkl,vij,vjk,vkl,Vij,Vjk,Mij,Mjk,Mkl,b,q,rij,rjk,rkl,upl,udc,bigM);
    end
end
dt1 = Data1(file)
#######################################
function dominated(x,P)
    st = false
    for k=1:length(P)
        if all( x .>= P[k])
            st=true; break
        else
            continue
        end
    end
    return st
end
function epmodel1dim()
    scnd1 = Model(optimizer_with_attributes(
            CPLEX.Optimizer,
            "CPX_PARAM_EPGAP" => 1e-8
          )); #CPLEX.Optimizer;
    set_silent(scnd1)
    # MOI.set(scnd1, MOI.RelativeGapTolerance(), 1e-8)
    # MOI.set(scnd1, MOI.NumberOfThreads(), 1);
    #########################  IP  ########################################
    @variable(scnd1, y1[1:(dt1.N["plant"]+dt1.N["distribution"])*2], Bin)
    @variable(scnd1, uij1[1:sum(dt1.Mij)], Bin);
    @variable(scnd1, ujk1[1:sum(dt1.Mjk)], Bin);
    @variable(scnd1, ukl1[1:sum(dt1.Mkl)], Bin);

    @variable( scnd1, 0<= xij1[1:sum(dt1.Mij)*5] );
    @variable( scnd1, 0<= xjk1[1:sum(dt1.Mjk)*5] );
    @variable( scnd1, 0<= xkl1[1:sum(dt1.Mkl)*5] );
    @variable( scnd1, 0<= h1[1:(dt1.N["plant"]+dt1.N["distribution"])*5*2] );
    obj1 = @expression(scnd1, sum(dt1.c.*y1) +sum(repeat(dt1.a[1,:], outer=sum(dt1.Mij[1,:])).*xij1[1:sum(dt1.Mij[1,:])*5])+
            sum(sum(repeat(dt1.a[i,:], outer=sum(dt1.Mij[i,:])).*xij1[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"])+
            sum(dt1.e.*h1) + sum(dt1.gij[i]*uij1[i] for i in findnz(dt1.gij)[1]) + sum(dt1.gjk[i]*ujk1[i] for i in findnz(dt1.gjk)[1]) + sum(dt1.gkl[i].*ukl1[i] for i in findnz(dt1.gkl)[1])+
            sum(dt1.vij.*xij1)+sum(dt1.vjk.*xjk1)+sum(dt1.vkl.*xkl1))
    obj2 = @expression(scnd1,sum(repeat(dt1.b[1,:], outer=sum(dt1.Mij[1,:])).*xij1[1:sum(dt1.Mij[1,:])*5]) +
            sum(sum(repeat(dt1.b[i,:], outer=sum(dt1.Mij[i,:])).*xij1[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"]) +
            sum(dt1.q.*h1) + sum(dt1.rij.*xij1)+sum(dt1.rjk.*xjk1)+sum(dt1.rkl.*xkl1))
    @variable(scnd1, 0 <= ep);
    @constraint(scnd1, epcon1, obj2 <= ep);
    @objective(scnd1, Min, obj1)
    # @constraint(scnd1, epcon1, obj1 <= ep);
    # @objective(scnd1, Min, obj2)
    @constraint(scnd1, [p=1:5], sum(xij1[5*(m-1)+p] for m=1:dt1.Mij[1,1])+sum(xij1[5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,1]) == sum(xjk1[5*(m-1)+p] for m=1:sum(dt1.Mjk[1,:])) );
    @constraint(scnd1, [j=2:dt1.N["plant"],p=1:5], sum(xij1[5*sum(dt1.Mij[1,1:j-1])+5*(m-1)+p] for m=1:dt1.Mij[1,j])+sum(xij1[5*sum(dt1.Mij[i,1:j-1])+5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,j]) == sum(xjk1[sum(dt1.Mjk[1:j-1,:])*5 + 5*(m-1)+p] for m=1:sum(dt1.Mjk[j,:])) );
    @constraint(scnd1, [p=1:5], sum(xjk1[5*(m-1)+p] for m=1:dt1.Mjk[1,1])+sum(xjk1[5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,1]) == sum(xkl1[5*(m-1)+p] for m=1:sum(dt1.Mkl[1,:])) );
    @constraint(scnd1, [k=2:dt1.N["distribution"],p=1:5],sum(xjk1[5*sum(dt1.Mjk[1,1:k-1])+5*(m-1)+p] for m=1:dt1.Mjk[1,k])+sum(xjk1[5*sum(dt1.Mjk[j,1:k-1])+5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,k]) == sum(xkl1[sum(dt1.Mkl[1:k-1,:])*5 + 5*(m-1)+p] for m=1:sum(dt1.Mkl[k,:])) );
    ########### constraint 4-6 #############
    @constraint(scnd1, [p=1:5],sum(h1[2*(p-1)+t] for t=1:2) == sum(xij1[5*(m-1)+p] for m=1:dt1.Mij[1,1])+sum(xij1[5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,1]));
    @constraint(scnd1, [j=2:dt1.N["plant"],p=1:5], sum(h1[5*2*(j-1)+2*(p-1)+t] for t=1:2) == sum(xij1[5*sum(dt1.Mij[1,1:j-1])+5*(m-1)+p] for m=1:dt1.Mij[1,j])+sum(xij1[5*sum(dt1.Mij[i,1:j-1])+5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,j]) );
    @constraint(scnd1, [p=1:5], sum(h1[5*2*dt1.N["plant"]+2*(p-1)+t] for t=1:2) == sum(xjk1[5*(m-1)+p] for m=1:dt1.Mjk[1,1])+sum(xjk1[5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,1]) );
    @constraint(scnd1, [k=2:dt1.N["distribution"],p=1:5], sum(h1[5*2*dt1.N["plant"]+5*2*(k-1)+2*(p-1)+t] for t=1:2) == sum(xjk1[5*sum(dt1.Mjk[1,1:k-1])+5*(m-1)+p] for m=1:dt1.Mjk[1,k])+sum(xjk1[5*sum(dt1.Mjk[j,1:k-1])+5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,k]));
    @constraint(scnd1, [p=1:5], sum(xkl1[5*(m-1)+p] for m=1:dt1.Mkl[1,1]) +sum(xkl1[5*(m-1)+p+(5*sum(dt1.Mkl[1:k-1,:]))] for k=2:dt1.N["distribution"] for m=1:dt1.Mkl[k,1]) >= dt1.d[1,p]);
    @constraint(scnd1, [l=2:dt1.N["customer"], p=1:5], sum(xkl1[sum(dt1.Mkl[1,1:l-1])*5 + 5*(m-1)+p] for m=1:dt1.Mkl[1,l])+ sum(xkl1[5*sum(dt1.Mkl[1:k-1,:])+5*sum(dt1.Mkl[k,1:l-1])+5*(m-1)+p] for k=2:dt1.N["distribution"] for m=1:dt1.Mkl[k,l]) >= dt1.d[l,p]);
    ########### constraint 7 #############
    @constraint(scnd1, sum(xij1[1:5*sum(dt1.Mij[1,:])]) <= dt1.N["cas"][1]);
    @constraint(scnd1, [i=2:dt1.N["supplier"]],  sum(xij1[5*sum(dt1.Mij[1:i-1,:])+1:5*sum(dt1.Mij[1:i,:])]) <= dt1.N["cas"][i]);
    ########### constraint 8 #############
    # @constraint(scnd1,[j=1:dt1.N["plant"]+dt1.N["distribution"], t=1:2], sum(h1[5*2*(j-1)+((t-1)*5)+1:5*2*(j-1)+((t-1)*5)+5]) <= [dt1.N["cap"];dt1.N["cad"]][j]*y1[(dt1.N["plant"]+dt1.N["distribution"])*(t-1)+j])
    @constraint(scnd1,[j=1:dt1.N["plant"]+dt1.N["distribution"], t=1:2], sum(h1[5*2*(j-1)+((p-1)*2)+t] for p=1:5) <= [dt1.N["cap"];dt1.N["cad"]][j]*y1[2*(j-1)+t]);
    ########### constraint 9 #############
    @constraint(scnd1,[j=1:dt1.N["plant"]+dt1.N["distribution"]], sum(y1[2*(j-1)+1:2*(j-1)+2]) <= 1);
    ########### constraint 10 #############
    @constraints(scnd1,begin sum(uij1[1:dt1.Mij[1,1]]) <= 1
        [j=2:dt1.N["plant"]], sum(uij1[sum(dt1.Mij[1,1:j-1])+1:sum(dt1.Mij[1,1:j-1])+dt1.Mij[1,j]]) <= 1
        [i=2:dt1.N["supplier"],j=2:dt1.N["plant"]],  sum(uij1[sum(dt1.Mij[1:i-1,:])+sum(dt1.Mij[i,1:j-1])+1:sum(dt1.Mij[1:i-1,:])+sum(dt1.Mij[i,1:j-1])+dt1.Mij[i,j]])<= 1
        sum(ujk1[1:dt1.Mjk[1,1]]) <= 1
        [k=2:dt1.N["distribution"]], sum(ujk1[sum(dt1.Mjk[1,1:k-1])+1:sum(dt1.Mjk[1,1:k-1])+dt1.Mjk[1,k]]) <= 1
        [j=2:dt1.N["plant"],k=2:dt1.N["distribution"]],  sum(ujk1[sum(dt1.Mjk[1:j-1,:])+sum(dt1.Mjk[j,1:j-1])+1:sum(dt1.Mjk[1:j-1,:])+sum(dt1.Mjk[j,1:j-1])+dt1.Mjk[j,k]]) <= 1
        sum(ukl1[1:dt1.Mkl[1,1]]) <= 1
        [l=2:dt1.N["customer"]], sum(ukl1[sum(dt1.Mkl[1,1:l-1])+1:sum(dt1.Mkl[1,1:l-1])+dt1.Mkl[1,l]]) <= 1
        [k=2:dt1.N["distribution"],l=2:dt1.N["customer"]],  sum(ukl1[sum(dt1.Mkl[1:k-1,:])+sum(dt1.Mkl[k,1:l-1])+1:sum(dt1.Mkl[1:k-1,:])+sum(dt1.Mkl[k,1:l-1])+dt1.Mkl[k,l]])<= 1
    end);
    ########### constraint 11 #############
    @constraints(scnd1, begin
        [i=1:sum(dt1.Mij)], sum(xij1[5*(i-1)+1:5*i]) <= dt1.bigM*uij1[i]
        [i=1:sum(dt1.Mjk)], sum(xjk1[5*(i-1)+1:5*i]) <= dt1.bigM*ujk1[i]
        [i=1:sum(dt1.Mkl)], sum(xkl1[5*(i-1)+1:5*i]) <= dt1.bigM*ukl1[i]
    end);
    ########### constraint 12 #############
    @constraints(scnd1, begin
        [i in findnz(dt1.Vij)[1]], sum(xij1[5*(i-1)+1:5*i]) >= dt1.Vij[i]*uij1[i]
        [i in findnz(dt1.Vjk)[1]], sum(xjk1[5*(i-1)+1:5*i]) >= dt1.Vjk[i]*ujk1[i]
        # [i in findnz(dt1.Vkl)[1]], sum(xkl1[5*(i-1)+1:5*i]) >= dt1.Vkl[i]*ukl1[i]
    end);
    ########### constraint 13-14 #############
    @constraint(scnd1, sum(y1[1:dt1.N["plant"]*2]) <= dt1.upl);
    @constraint(scnd1, sum(y1[dt1.N["plant"]*2+1:end]) <= dt1.udc);
    return scnd1
end
m1 = epmodel1dim()
function opt1dim(ϵ,TL)
    JuMP.fix(m1[:ep], ϵ; force = true);
    optimize!(m1)
    set_time_limit_sec(m1,TL)
    if termination_status(m1) == MOI.OPTIMAL
        # y1 = round.(JuMP.value.(m1[:y1]))
        # h1 = round.(JuMP.value.(m1[:h1]); digits=4)
        # uij1 =round.( JuMP.value.(m1[:uij1]))
        # ujk1 = round.(JuMP.value.(m1[:ujk1]))
        # ukl1 = round.(JuMP.value.(m1[:ukl1]))
        # xij1 = round.(JuMP.value.(m1[:xij1]); digits=4)
        # xjk1 = round.(JuMP.value.(m1[:xjk1]); digits=4)
        # xkl1 = round.(JuMP.value.(m1[:xkl1]); digits=4)
        y1 = value.(m1[:y1]);
        uij1 = value.(m1[:uij1]);
        ujk1 = value.(m1[:ujk1]);
        ukl1 = value.(m1[:ukl1]);
        xij1 = value.(m1[:xij1]);
        xjk1 = value.(m1[:xjk1]);
        xkl1 = value.(m1[:xkl1]);
        h1 = value.(m1[:h1]) ;
        obj1 = sum(dt1.c.*y1) +sum(repeat(dt1.a[1,:], outer=sum(dt1.Mij[1,:])).*xij1[1:sum(dt1.Mij[1,:])*5])+
                sum(sum(repeat(dt1.a[i,:], outer=sum(dt1.Mij[i,:])).*xij1[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"])+
                sum(dt1.e.*h1) + sum(dt1.gij[i]*uij1[i] for i in findnz(dt1.gij)[1]) + sum(dt1.gjk[i]*ujk1[i] for i in findnz(dt1.gjk)[1]) + sum(dt1.gkl[i].*ukl1[i] for i in findnz(dt1.gkl)[1])+
                sum(dt1.vij.*xij1)+sum(dt1.vjk.*xjk1)+sum(dt1.vkl.*xkl1)
        obj2 = sum(repeat(dt1.b[1,:], outer=sum(dt1.Mij[1,:])).*xij1[1:sum(dt1.Mij[1,:])*5]) +
                sum(sum(repeat(dt1.b[i,:], outer=sum(dt1.Mij[i,:])).*xij1[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"]) +
                sum(dt1.q.*h1) + sum(dt1.rij.*xij1)+sum(dt1.rjk.*xjk1)+sum(dt1.rkl.*xkl1)
        return [objective_value(m1),obj2]
        # return [objective_value(m1),obj1]
    elseif termination_status(m1) == MOI.INFEASIBLE
        return 0
    else
        return nothing
    end
end
function epsilon(TL)
    Y = []; ϵ = 1.2665377400539995e6; δ =40000.00000000074; lb = 866537.740053992; fval = [0,ϵ]
    iter = 0;
    while fval[2] >= lb && iter<11
        fval = opt1dim(ϵ,TL)
        # println(fval)
        # println([fval[2],fval[1]])
        if fval == 0
            break
        elseif fval == nothing
                fval = [0,ϵ]
        else dominated(fval,Y)==false
            push!(Y,fval);
        end
        ϵ = ϵ-δ; iter+=1;
    end
    return Y
end
ey = epsilon(300)
1
##############################

function epmodel1dim()
    scnd1 = Model(CPLEX.Optimizer); set_silent(scnd1)
    # MOI.set(scnd1, MOI.NumberOfThreads(), 1);
    #########################  IP  ########################################
    @variable(scnd1, y1[1:(dt1.N["plant"]+dt1.N["distribution"])*2], Bin)
    @variable(scnd1, uij1[1:sum(dt1.Mij)], Bin);
    @variable(scnd1, ujk1[1:sum(dt1.Mjk)], Bin);
    @variable(scnd1, ukl1[1:sum(dt1.Mkl)], Bin);

    @variable( scnd1, 0<= xij1[1:sum(dt1.Mij)*5] );
    @variable( scnd1, 0<= xjk1[1:sum(dt1.Mjk)*5] );
    @variable( scnd1, 0<= xkl1[1:sum(dt1.Mkl)*5] );
    @variable( scnd1, 0<= h1[1:(dt1.N["plant"]+dt1.N["distribution"])*5*2] );
    obj1 = @expression(scnd1, sum(dt1.c.*y1) +sum(repeat(dt1.a[1,:], outer=sum(dt1.Mij[1,:])).*xij1[1:sum(dt1.Mij[1,:])*5])+
            sum(sum(repeat(dt1.a[i,:], outer=sum(dt1.Mij[i,:])).*xij1[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"])+
            sum(dt1.e.*h1) + sum(dt1.gij[i]*uij1[i] for i in findnz(dt1.gij)[1]) + sum(dt1.gjk[i]*ujk1[i] for i in findnz(dt1.gjk)[1]) + sum(dt1.gkl[i].*ukl1[i] for i in findnz(dt1.gkl)[1])+
            sum(dt1.vij.*xij1)+sum(dt1.vjk.*xjk1)+sum(dt1.vkl.*xkl1))
    obj2 = @expression(scnd1,sum(repeat(dt1.b[1,:], outer=sum(dt1.Mij[1,:])).*xij1[1:sum(dt1.Mij[1,:])*5]) +
            sum(sum(repeat(dt1.b[i,:], outer=sum(dt1.Mij[i,:])).*xij1[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"]) +
            sum(dt1.q.*h1) + sum(dt1.rij.*xij1)+sum(dt1.rjk.*xjk1)+sum(dt1.rkl.*xkl1))
    @variable(scnd1, 0 <= ep);
    @constraint(scnd1, epcon1, obj1 <= ep);
    @objective(scnd1, Min, obj2)
    @constraint(scnd1, [p=1:5], sum(xij1[5*(m-1)+p] for m=1:dt1.Mij[1,1])+sum(xij1[5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,1]) == sum(xjk1[5*(m-1)+p] for m=1:sum(dt1.Mjk[1,:])) );
    @constraint(scnd1, [j=2:dt1.N["plant"],p=1:5], sum(xij1[5*sum(dt1.Mij[1,1:j-1])+5*(m-1)+p] for m=1:dt1.Mij[1,j])+sum(xij1[5*sum(dt1.Mij[i,1:j-1])+5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,j]) == sum(xjk1[sum(dt1.Mjk[1:j-1,:])*5 + 5*(m-1)+p] for m=1:sum(dt1.Mjk[j,:])) );
    @constraint(scnd1, [p=1:5], sum(xjk1[5*(m-1)+p] for m=1:dt1.Mjk[1,1])+sum(xjk1[5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,1]) == sum(xkl1[5*(m-1)+p] for m=1:sum(dt1.Mkl[1,:])) );
    @constraint(scnd1, [k=2:dt1.N["distribution"],p=1:5],sum(xjk1[5*sum(dt1.Mjk[1,1:k-1])+5*(m-1)+p] for m=1:dt1.Mjk[1,k])+sum(xjk1[5*sum(dt1.Mjk[j,1:k-1])+5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,k]) == sum(xkl1[sum(dt1.Mkl[1:k-1,:])*5 + 5*(m-1)+p] for m=1:sum(dt1.Mkl[k,:])) );
    ########### constraint 4-6 #############
    @constraint(scnd1, [p=1:5],sum(h1[2*(p-1)+t] for t=1:2) == sum(xij1[5*(m-1)+p] for m=1:dt1.Mij[1,1])+sum(xij1[5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,1]));
    @constraint(scnd1, [j=2:dt1.N["plant"],p=1:5], sum(h1[5*2*(j-1)+2*(p-1)+t] for t=1:2) == sum(xij1[5*sum(dt1.Mij[1,1:j-1])+5*(m-1)+p] for m=1:dt1.Mij[1,j])+sum(xij1[5*sum(dt1.Mij[i,1:j-1])+5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,j]) );
    @constraint(scnd1, [p=1:5], sum(h1[5*2*dt1.N["plant"]+2*(p-1)+t] for t=1:2) == sum(xjk1[5*(m-1)+p] for m=1:dt1.Mjk[1,1])+sum(xjk1[5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,1]) );
    @constraint(scnd1, [k=2:dt1.N["distribution"],p=1:5], sum(h1[5*2*dt1.N["plant"]+5*2*(k-1)+2*(p-1)+t] for t=1:2) == sum(xjk1[5*sum(dt1.Mjk[1,1:k-1])+5*(m-1)+p] for m=1:dt1.Mjk[1,k])+sum(xjk1[5*sum(dt1.Mjk[j,1:k-1])+5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,k]));
    @constraint(scnd1, [p=1:5], sum(xkl1[5*(m-1)+p] for m=1:dt1.Mkl[1,1]) +sum(xkl1[5*(m-1)+p+(5*sum(dt1.Mkl[1:k-1,:]))] for k=2:dt1.N["distribution"] for m=1:dt1.Mkl[k,1]) >= dt1.d[1,p]);
    @constraint(scnd1, [l=2:dt1.N["customer"], p=1:5], sum(xkl1[sum(dt1.Mkl[1,1:l-1])*5 + 5*(m-1)+p] for m=1:dt1.Mkl[1,l])+ sum(xkl1[5*sum(dt1.Mkl[1:k-1,:])+5*sum(dt1.Mkl[k,1:l-1])+5*(m-1)+p] for k=2:dt1.N["distribution"] for m=1:dt1.Mkl[k,l]) >= dt1.d[l,p]);
    ########### constraint 7 #############
    @constraint(scnd1, sum(xij1[1:5*sum(dt1.Mij[1,:])]) <= dt1.N["cas"][1]);
    @constraint(scnd1, [i=2:dt1.N["supplier"]],  sum(xij1[5*sum(dt1.Mij[1:i-1,:])+1:5*sum(dt1.Mij[1:i,:])]) <= dt1.N["cas"][i]);
    ########### constraint 8 #############
    # @constraint(scnd1,[j=1:dt1.N["plant"]+dt1.N["distribution"], t=1:2], sum(h1[5*2*(j-1)+((t-1)*5)+1:5*2*(j-1)+((t-1)*5)+5]) <= [dt1.N["cap"];dt1.N["cad"]][j]*y1[(dt1.N["plant"]+dt1.N["distribution"])*(t-1)+j])
    @constraint(scnd1,[j=1:dt1.N["plant"]+dt1.N["distribution"], t=1:2], sum(h1[5*2*(j-1)+((p-1)*2)+t] for p=1:5) <= [dt1.N["cap"];dt1.N["cad"]][j]*y1[2*(j-1)+t]);
    ########### constraint 9 #############
    @constraint(scnd1,[j=1:dt1.N["plant"]+dt1.N["distribution"]], sum(y1[2*(j-1)+1:2*(j-1)+2]) <= 1);
    ########### constraint 10 #############
    @constraints(scnd1,begin sum(uij1[1:dt1.Mij[1,1]]) <= 1
        [j=2:dt1.N["plant"]], sum(uij1[sum(dt1.Mij[1,1:j-1])+1:sum(dt1.Mij[1,1:j-1])+dt1.Mij[1,j]]) <= 1
        [i=2:dt1.N["supplier"],j=2:dt1.N["plant"]],  sum(uij1[sum(dt1.Mij[1:i-1,:])+sum(dt1.Mij[i,1:j-1])+1:sum(dt1.Mij[1:i-1,:])+sum(dt1.Mij[i,1:j-1])+dt1.Mij[i,j]])<= 1
        sum(ujk1[1:dt1.Mjk[1,1]]) <= 1
        [k=2:dt1.N["distribution"]], sum(ujk1[sum(dt1.Mjk[1,1:k-1])+1:sum(dt1.Mjk[1,1:k-1])+dt1.Mjk[1,k]]) <= 1
        [j=2:dt1.N["plant"],k=2:dt1.N["distribution"]],  sum(ujk1[sum(dt1.Mjk[1:j-1,:])+sum(dt1.Mjk[j,1:j-1])+1:sum(dt1.Mjk[1:j-1,:])+sum(dt1.Mjk[j,1:j-1])+dt1.Mjk[j,k]]) <= 1
        sum(ukl1[1:dt1.Mkl[1,1]]) <= 1
        [l=2:dt1.N["customer"]], sum(ukl1[sum(dt1.Mkl[1,1:l-1])+1:sum(dt1.Mkl[1,1:l-1])+dt1.Mkl[1,l]]) <= 1
        [k=2:dt1.N["distribution"],l=2:dt1.N["customer"]],  sum(ukl1[sum(dt1.Mkl[1:k-1,:])+sum(dt1.Mkl[k,1:l-1])+1:sum(dt1.Mkl[1:k-1,:])+sum(dt1.Mkl[k,1:l-1])+dt1.Mkl[k,l]])<= 1
    end);
    ########### constraint 11 #############
    @constraints(scnd1, begin
        [i=1:sum(dt1.Mij)], sum(xij1[5*(i-1)+1:5*i]) <= dt1.bigM*uij1[i]
        [i=1:sum(dt1.Mjk)], sum(xjk1[5*(i-1)+1:5*i]) <= dt1.bigM*ujk1[i]
        [i=1:sum(dt1.Mkl)], sum(xkl1[5*(i-1)+1:5*i]) <= dt1.bigM*ukl1[i]
    end);
    ########### constraint 12 #############
    @constraints(scnd1, begin
        [i in findnz(dt1.Vij)[1]], sum(xij1[5*(i-1)+1:5*i]) >= dt1.Vij[i]*uij1[i]
        [i in findnz(dt1.Vjk)[1]], sum(xjk1[5*(i-1)+1:5*i]) >= dt1.Vjk[i]*ujk1[i]
        # [i in findnz(dt1.Vkl)[1]], sum(xkl1[5*(i-1)+1:5*i]) >= dt1.Vkl[i]*ukl1[i]
    end);
    ########### constraint 13-14 #############
    @constraint(scnd1, sum(y1[1:dt1.N["plant"]*2]) <= dt1.upl);
    @constraint(scnd1, sum(y1[dt1.N["plant"]*2+1:end]) <= dt1.udc);
    return scnd1
end
m1 = epmodel1dim()

function opt1dim(ϵ)
    JuMP.fix(m1[:ep], ϵ; force = true);
    optimize!(m1)
    if termination_status(m1) == MOI.OPTIMAL
        y1 = value.(m1[:y1]);
        uij1 = value.(m1[:uij1]);
        ujk1 = value.(m1[:ujk1]);
        ukl1 = value.(m1[:ukl1]);
        xij1 = value.(m1[:xij1]);
        xjk1 = value.(m1[:xjk1]);
        xkl1 = value.(m1[:xkl1]);
        h1 = value.(m1[:h1]) ;
        obj1 = sum(dt1.c.*y1) +sum(repeat(dt1.a[1,:], outer=sum(dt1.Mij[1,:])).*xij1[1:sum(dt1.Mij[1,:])*5])+
                sum(sum(repeat(dt1.a[i,:], outer=sum(dt1.Mij[i,:])).*xij1[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"])+
                sum(dt1.e.*h1) + sum(dt1.gij[i]*uij1[i] for i in findnz(dt1.gij)[1]) + sum(dt1.gjk[i]*ujk1[i] for i in findnz(dt1.gjk)[1]) + sum(dt1.gkl[i].*ukl1[i] for i in findnz(dt1.gkl)[1])+
                sum(dt1.vij.*xij1)+sum(dt1.vjk.*xjk1)+sum(dt1.vkl.*xkl1)
        obj2 = sum(repeat(dt1.b[1,:], outer=sum(dt1.Mij[1,:])).*xij1[1:sum(dt1.Mij[1,:])*5]) +
                sum(sum(repeat(dt1.b[i,:], outer=sum(dt1.Mij[i,:])).*xij1[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"]) +
                sum(dt1.q.*h1) + sum(dt1.rij.*xij1)+sum(dt1.rjk.*xjk1)+sum(dt1.rkl.*xkl1)
        return [objective_value(m1),obj1]
    else
        return nothing
    end
end
function epsilon()
    # Test4S4
    Y = []; ϵ = 1.2347893923481262e8 ; δ =10; lb = 99*10^(6); fval = [0,ϵ]
    while fval[2] >= lb
        fval = opt1dim(ϵ)
        println([fval[2],fval[1]])
        if fval == nothing
            break
        # end
        # if dominated(fval,Y)==false
        else
            push!(Y,fval);
        end
        ϵ = ϵ-δ
    end
    return Y
end
ey = epsilon()


# optimize!(m1)
# objective_value(m1)
# function dominated(y,P)
#     st = false
#     for k=1:length(P)
#         if all( y .>= P[k])# && any(x > P[k])
#             st = true; break
#         else
#             continue
#         end
#     end
#     return st
# end

domFilter(ex,ey)
#
# dot(ex[1:27],dtt.C[1][1:27])
# ENV["CPLEX_STUDIO_BINARIES"] = "F:/cplex12.9/cplex/bin/x64_win64/"
# import Pkg
# Pkg.add(Pkg.PackageSpec(name = "CPLEX", version = v"0.6"))

struct MyData
    filepath::String; N::Dict{}; d::Array{}; c::Array{};  Mij::Array{}; Mjk::Array{}; Mkl::Array{};
    vij::Array{}; vjk::Array{}; vkl::Array{}; rij::Array{}; rjk::Array{}; rkl::Array{}; #gij::Array{}; gjk::Array{}; gkl::Array{};
    Vij::Array{}; Vjk::Array{}; b::Array{}; upl::Int; udc::Int; bigM::Int # e::Array{};q::Array{};
    function MyData(filepath)
        dt = readdlm(filepath);
        # notafile = readdlm("/home/ak121396/Desktop/instances/SCND/Notations.txt", '=');
        # notafile = readdlm("F:/scnd/Notations.txt", '=');
        notafile = readdlm("/home/k2g00/k2g3475/scnd/Notations.txt", '=');
        nota = notafile[1:end,1];
        N= Dict();
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
        upl = N["upperpants"]; udc = N["upperdistribution"]
        bigM = sum(sum(N["demand"]))
        new(filepath,N,d,c,Mij,Mjk,Mkl,vij,vjk,vkl,rij,rjk,rkl,Vij,Vjk,b,upl,udc,bigM); #cap,Mij,Mjk,Mkl,gij,gjk,gkl,
    end
end
# file = "/home/ak121396/Desktop/instances/SCND/test04S4"
dt = MyData(file);
function epmodel()
    ##############################  MIP   #####################################
    scnd = Model(CPLEX.Optimizer); set_silent(scnd);
    #MOI.set(scnd, MOI.NumberOfThreads(), 1);
    @variable(scnd, y[1:dt.N["plant"]+dt.N["distribution"],1:2], Bin);
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

    obj1 = @expression( scnd, sum(dt.c[j][t]*y[j,t] for j=1:1:dt.N["plant"]+dt.N["distribution"] for t=1:2) + exg +
        sum(dt.N["vcs"][i][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) +
        sum(dt.vij[i][j][m][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
        sum(dt.vjk[j][k][m][p]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
        sum(dt.vkl[k][l][m][p]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5)+
        sum(dt.N["vcp"][j][2*(p-1)+t]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
        sum(dt.N["vcd"][k][2*(p-1)+t]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2) )
    obj2 = @expression(scnd, sum(dt.b[i,p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
        sum(dt.N["vep"][j][2*(p-1)+t]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
        sum(dt.N["ved"][k][2*(p-1)+t]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2)+
        sum(dt.rij[i][j][m]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
        sum(dt.rjk[j][k][m]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
        sum(dt.rkl[k][l][m]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5))
    @variable(scnd, 0 <= ep);
    @constraint(scnd, epcon1, obj2 <= ep);
    @objective(scnd, Min, obj1)
    # 1st obj
    # @objective(scnd, Min, sum(dt.c[j][t]*y[j,t] for j=1:1:dt.N["plant"]+dt.N["distribution"] for t=1:2) + exg +
    #     # sum(dt.gij[i][j][m]*uij[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j])+
    #     # sum(dt.gjk[j][k][m]*ujk[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]) +
    #     # sum(dt.gkl[k][l][m]*ukl[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])+
    #     sum(dt.N["vcs"][i][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) +
    #     sum(dt.vij[i][j][m][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
    #     sum(dt.vjk[j][k][m][p]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
    #     sum(dt.vkl[k][l][m][p]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5)+
    #     sum(dt.N["vcp"][j][2*(p-1)+t]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
    #     sum(dt.N["vcd"][k][2*(p-1)+t]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2));
    #2nd obj
    # @objective(scnd, Min, sum(dt.b[i,p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
    #     sum(dt.N["vep"][j][2*(p-1)+t]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
    #     sum(dt.N["ved"][k][2*(p-1)+t]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2)+
    #     sum(dt.rij[i][j][m]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
    #     sum(dt.rjk[j][k][m]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
    #     sum(dt.rkl[k][l][m]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5));
    ######### constraint 3 #############
    @constraint(scnd,[j=1:dt.N["plant"],p=1:5], sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j]) == sum(xjk[j,k,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]));
    @constraint(scnd,[k=1:dt.N["distribution"],p=1:5], sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k]) == sum(xkl[k,l,m,p] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l]));
    ########### constraint 4-6 #############
    @constraints(scnd, begin
        [j=1:dt.N["plant"],p=1:5], sum(h[j,p,:]) == sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j])
        [k=1:dt.N["distribution"],p=1:5], sum(h[k+dt.N["plant"],p,:] ) == sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k])
        [l=1:dt.N["customer"],p=1:5], sum(xkl[k,l,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mkl[k,l]) >= dt.d[l][p]
    end );
    # @constraint(scnd,[j=1:dt.N["plant"],p=1:5], sum(h[j,p,:]) == sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j]));
    # @constraint(scnd,[k=1:dt.N["distribution"],p=1:5], sum(h[k+dt.N["plant"],p,:] ) == sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k]));
    # @constraint(scnd, [l=1:dt.N["customer"],p=1:5], sum(xkl[k,l,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mkl[k,l]) >= dt.d[l][p]);
    ########### constraint 7-9 #############
    @constraint(scnd,[i=1:dt.N["supplier"]], sum(xij[i,j,m,p] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) <= dt.N["cas"][i]);
    @constraint(scnd,[j=1:dt.N["plant"]+dt.N["distribution"], t=1:2], sum(h[j,:,t]) <=[dt.N["cap"];dt.N["cad"]][j]*y[j,t]);
    @constraint(scnd,[j=1:dt.N["plant"]+dt.N["distribution"]], sum(y[j,:]) <= 1);
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
    # @constraint(scnd,[k=1:dt.N["distribution"], l=1:dt.N["customer"], m=1:dt.Mkl[k,l]], sum(xkl[k,l,m,p] for p=1:5) >= dt.Vkl[k][l][m]*ukl[k,l,m]);
    ########### constraint 13-14 #############
    @constraint(scnd,sum(y[j,t] for j=1:dt.N["plant"] for t=1:2) <= dt.upl);
    @constraint(scnd,sum(y[j,t] for j=dt.N["plant"]+1:dt.N["distribution"]+dt.N["plant"] for t=1:2) <= dt.udc);
    return scnd
end
m1 = epmodel()
function opt(ϵ)
    JuMP.fix(m1[:ep], ϵ; force = true);
    optimize!(m1)
    if termination_status(m1) == MOI.OPTIMAL
        y = JuMP.value.(m1[:y]); h = JuMP.value.(m1[:h])
        uij = JuMP.value.(m1[:uij])
        ujk = JuMP.value.(m1[:ujk])
        ukl = JuMP.value.(m1[:ukl])
        xij = JuMP.value.(m1[:xij])
        xjk = JuMP.value.(m1[:xjk])
        xkl = JuMP.value.(m1[:xkl])
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
        obj1 = sum(dt.c[j][t]*y[j,t] for j=1:1:dt.N["plant"]+dt.N["distribution"] for t=1:2) + exg +
            sum(dt.N["vcs"][i][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) +
            sum(dt.vij[i][j][m][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
            sum(dt.vjk[j][k][m][p]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
            sum(dt.vkl[k][l][m][p]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5)+
            sum(dt.N["vcp"][j][2*(p-1)+t]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
            sum(dt.N["vcd"][k][2*(p-1)+t]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2)
        obj2 = sum(dt.b[i,p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
            sum(dt.N["vep"][j][2*(p-1)+t]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
            sum(dt.N["ved"][k][2*(p-1)+t]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2)+
            sum(dt.rij[i][j][m]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
            sum(dt.rjk[j][k][m]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
            sum(dt.rkl[k][l][m]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5)
        return [objective_value(m1),obj2]
    else
        return nothing
    end
end
