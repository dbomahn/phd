using CPUTime,DataFrames,DelimitedFiles,JuMP,LinearAlgebra,CPLEX,vOptGeneric,SparseArrays,StatsBase
using CSV,JLD2
#########################  1dim model  ############################
struct Data1dim
    file::String; N::Dict{}; d::Array{}; c::Array{}; a::Array{}; e::Array{}; gij::SparseVector{}; gjk::SparseVector{}; gkl::SparseVector{};
    vij::Array{}; vjk::Array{}; vkl::Array{}; Vij::SparseVector{}; Vjk::SparseVector{}; Mij::Array{}; Mjk::Array{}; Mkl::Array{};
    b::Array{}; q::Array{}; rij::Array{}; rjk::Array{}; rkl::Array{}; upl::Int; udc::Int; bigM::Int
    function Data1dim(file)
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
# file = "/home/k2g00/k2g3475/scnd/instances/test01S2"
@show file = ARGS[1]
# file = "F:scnd/Test1S2"
# file = "/home/ak121396/Desktop/instances/SCND/test01S2"

dt1 = Data1dim(file);
if dt1.N["supplier"] >=6 && dt1.N["supplier"] <= 7
    TL = dt1.N["supplier"]*60*1
elseif dt1.N["supplier"]>7 && dt1.N["supplier"] <= 9
    TL = dt1.N["supplier"]*60*2
elseif dt1.N["supplier"]>9 && dt1.N["supplier"] <= 12
    TL = dt1.N["supplier"]*60*3
elseif dt1.N["supplier"]>12 && dt1.N["supplier"] <= 16
    TL = dt1.N["supplier"]*60*4
elseif dt1.N["supplier"]>16 && dt1.N["supplier"] <= 20
    TL = dt1.N["supplier"]*60*5
elseif dt1.N["supplier"]>20 && dt1.N["supplier"] <= 26
    TL = dt1.N["supplier"]*60*6
# elseif dt1.N["supplier"] == 26
    # TL = dt1.N["supplier"]*60*7
# elseif dt1.N["supplier"]
#     TL = dt1.N["supplier"]*60*8
else
    TL = dt1.N["supplier"]*60*8
end

function SCND_LP()
    scnd1 = vModel(CPLEX.Optimizer)
    # optimizer_with_attributes(
    #         CPLEX.Optimizer,
    #         "CPX_PARAM_EPGAP" => 1e-8
    #       ));
    set_silent(scnd1)
    MOI.set(scnd1, MOI.NumberOfThreads(), 1)
    @variable(scnd1, 0<=y[1:(dt1.N["plant"]+dt1.N["distribution"])*2]<=1)
    @variable(scnd1, 0<=uij[1:sum(dt1.Mij)]<=1);
    @variable(scnd1, 0<=ujk[1:sum(dt1.Mjk)]<=1);
    @variable(scnd1, 0<=ukl[1:sum(dt1.Mkl)]<=1);
    @variable( scnd1, 0<= xij[1:sum(dt1.Mij)*5] );
    @variable( scnd1, 0<= xjk[1:sum(dt1.Mjk)*5] );
    @variable( scnd1, 0<= xkl[1:sum(dt1.Mkl)*5] );
    @variable( scnd1, 0<= h[1:(dt1.N["plant"]+dt1.N["distribution"])*5*2] );

    @addobjective(scnd1, Min, sum(dt1.c.*y) +sum(repeat(dt1.a[1,:], outer=sum(dt1.Mij[1,:])).*xij[1:sum(dt1.Mij[1,:])*5])+
            sum(sum(repeat(dt1.a[i,:], outer=sum(dt1.Mij[i,:])).*xij[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"])+
            sum(dt1.e.*h) + sum(dt1.gij[i]*uij[i] for i in findnz(dt1.gij)[1]) + sum(dt1.gjk[i]*ujk[i] for i in findnz(dt1.gjk)[1]) + sum(dt1.gkl[i].*ukl[i] for i in findnz(dt1.gkl)[1])+
            sum(dt1.vij.*xij)+sum(dt1.vjk.*xjk)+sum(dt1.vkl.*xkl));
    @addobjective(scnd1, Min, sum(repeat(dt1.b[1,:], outer=sum(dt1.Mij[1,:])).*xij[1:sum(dt1.Mij[1,:])*5]) +
            sum(sum(repeat(dt1.b[i,:], outer=sum(dt1.Mij[i,:])).*xij[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"]) +
            sum(dt1.q.*h) + sum(dt1.rij.*xij)+sum(dt1.rjk.*xjk)+sum(dt1.rkl.*xkl));

    ########## constraint 3 #############
    @constraint(scnd1, [p=1:5], sum(xij[5*(m-1)+p] for m=1:dt1.Mij[1,1])+sum(xij[5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,1]) == sum(xjk[5*(m-1)+p] for m=1:sum(dt1.Mjk[1,:])) );
    @constraint(scnd1, [j=2:dt1.N["plant"],p=1:5], sum(xij[5*sum(dt1.Mij[1,1:j-1])+5*(m-1)+p] for m=1:dt1.Mij[1,j])+sum(xij[5*sum(dt1.Mij[i,1:j-1])+5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,j]) == sum(xjk[sum(dt1.Mjk[1:j-1,:])*5 + 5*(m-1)+p] for m=1:sum(dt1.Mjk[j,:])) );
    @constraint(scnd1, [p=1:5], sum(xjk[5*(m-1)+p] for m=1:dt1.Mjk[1,1])+sum(xjk[5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,1]) == sum(xkl[5*(m-1)+p] for m=1:sum(dt1.Mkl[1,:])) );
    @constraint(scnd1, [k=2:dt1.N["distribution"],p=1:5],sum(xjk[5*sum(dt1.Mjk[1,1:k-1])+5*(m-1)+p] for m=1:dt1.Mjk[1,k])+sum(xjk[5*sum(dt1.Mjk[j,1:k-1])+5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,k]) == sum(xkl[sum(dt1.Mkl[1:k-1,:])*5 + 5*(m-1)+p] for m=1:sum(dt1.Mkl[k,:])) );
    ########### constraint 4-6 #############
    @constraint(scnd1, [p=1:5],sum(h[2*(p-1)+t] for t=1:2) == sum(xij[5*(m-1)+p] for m=1:dt1.Mij[1,1])+sum(xij[5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,1]));
    @constraint(scnd1, [j=2:dt1.N["plant"],p=1:5], sum(h[5*2*(j-1)+2*(p-1)+t] for t=1:2) == sum(xij[5*sum(dt1.Mij[1,1:j-1])+5*(m-1)+p] for m=1:dt1.Mij[1,j])+sum(xij[5*sum(dt1.Mij[i,1:j-1])+5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,j]) );
    @constraint(scnd1, [p=1:5], sum(h[5*2*dt1.N["plant"]+2*(p-1)+t] for t=1:2) == sum(xjk[5*(m-1)+p] for m=1:dt1.Mjk[1,1])+sum(xjk[5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,1]) );
    @constraint(scnd1, [k=2:dt1.N["distribution"],p=1:5], sum(h[5*2*dt1.N["plant"]+5*2*(k-1)+2*(p-1)+t] for t=1:2) == sum(xjk[5*sum(dt1.Mjk[1,1:k-1])+5*(m-1)+p] for m=1:dt1.Mjk[1,k])+sum(xjk[5*sum(dt1.Mjk[j,1:k-1])+5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,k]));
    @constraint(scnd1, [p=1:5], sum(xkl[5*(m-1)+p] for m=1:dt1.Mkl[1,1]) +sum(xkl[5*(m-1)+p+(5*sum(dt1.Mkl[1:k-1,:]))] for k=2:dt1.N["distribution"] for m=1:dt1.Mkl[k,1]) >= dt1.d[1,p]);
    @constraint(scnd1, [l=2:dt1.N["customer"], p=1:5], sum(xkl[sum(dt1.Mkl[1,1:l-1])*5 + 5*(m-1)+p] for m=1:dt1.Mkl[1,l])+ sum(xkl[5*sum(dt1.Mkl[1:k-1,:])+5*sum(dt1.Mkl[k,1:l-1])+5*(m-1)+p] for k=2:dt1.N["distribution"] for m=1:dt1.Mkl[k,l]) >= dt1.d[l,p]);
    ########### constraint 7 #############
    @constraint(scnd1, sum(xij[1:5*sum(dt1.Mij[1,:])]) <= dt1.N["cas"][1]);
    @constraint(scnd1, [i=2:dt1.N["supplier"]],  sum(xij[5*sum(dt1.Mij[1:i-1,:])+1:5*sum(dt1.Mij[1:i,:])]) <= dt1.N["cas"][i]);
    ########### constraint 8 #############
    @constraint(scnd1,[j=1:dt1.N["plant"]+dt1.N["distribution"], t=1:2], sum(h[5*2*(j-1)+((p-1)*2)+t] for p=1:5) <= [dt1.N["cap"];dt1.N["cad"]][j]*y[2*(j-1)+t]);
    ########### constraint 9 #############
    @constraint(scnd1,[j=1:dt1.N["plant"]+dt1.N["distribution"]], sum(y[2*(j-1)+1:2*(j-1)+2]) <= 1);
    ########### constraint 10 #############
    @constraints(scnd1,begin
        sum(uij[1:dt1.Mij[1,1]]) <= 1
        sum(uij[sum(dt1.Mij[1,:])+dt1.Mij[2,1]]) <= 1
        [j=2:dt1.N["plant"]], sum(uij[sum(dt1.Mij[1,1:j-1])+1:sum(dt1.Mij[1,1:j-1])+dt1.Mij[1,j]]) <= 1
        [i=2:dt1.N["supplier"],j=2:dt1.N["plant"]],  sum(uij[sum(dt1.Mij[1:i-1,:])+sum(dt1.Mij[i,1:j-1])+1:sum(dt1.Mij[1:i-1,:])+sum(dt1.Mij[i,1:j-1])+dt1.Mij[i,j]])<= 1
        sum(ujk[1:dt1.Mjk[1,1]]) <= 1
        sum(ujk[sum(dt1.Mjk[1,:])+dt1.Mjk[2,1]]) <= 1
        [k=2:dt1.N["distribution"]], sum(ujk[sum(dt1.Mjk[1,1:k-1])+1:sum(dt1.Mjk[1,1:k-1])+dt1.Mjk[1,k]]) <= 1
        [j=2:dt1.N["plant"],k=2:dt1.N["distribution"]],  sum(ujk[sum(dt1.Mjk[1:j-1,:])+sum(dt1.Mjk[j,1:j-1])+1:sum(dt1.Mjk[1:j-1,:])+sum(dt1.Mjk[j,1:j-1])+dt1.Mjk[j,k]]) <= 1
        sum(ukl[1:dt1.Mkl[1,1]]) <= 1 #[1,1]
        sum(ukl[sum(dt1.Mkl[1,:])+dt1.Mkl[2,1]]) <= 1 #[2,1]
        [l=2:dt1.N["customer"]], sum(ukl[sum(dt1.Mkl[1,1:l-1])+1:sum(dt1.Mkl[1,1:l-1])+dt1.Mkl[1,l]]) <= 1
        [k=2:dt1.N["distribution"],l=2:dt1.N["customer"]],  sum(ukl[sum(dt1.Mkl[1:k-1,:])+sum(dt1.Mkl[k,1:l-1])+1:sum(dt1.Mkl[1:k-1,:])+sum(dt1.Mkl[k,1:l-1])+dt1.Mkl[k,l]])<= 1
    end);
    ########### constraint 11 #############
    @constraints(scnd1, begin
        [i=1:sum(dt1.Mij)], sum(xij[5*(i-1)+1:5*i]) <= dt1.bigM*uij[i]
        [i=1:sum(dt1.Mjk)], sum(xjk[5*(i-1)+1:5*i]) <= dt1.bigM*ujk[i]
        [i=1:sum(dt1.Mkl)], sum(xkl[5*(i-1)+1:5*i]) <= dt1.bigM*ukl[i]
    end);
    ########### constraint 12 #############
    @constraints(scnd1, begin
            [i in findnz(dt1.Vij)[1]], sum(xij[5*(i-1)+1:5*i]) >= dt1.Vij[i]*uij[i]
            [i in findnz(dt1.Vjk)[1]], sum(xjk[5*(i-1)+1:5*i]) >= dt1.Vjk[i]*ujk[i]
            # [i in findnz(dt1.Vkl)[1]], sum(xkl[5*(i-1)+1:5*i]) >= dt1.Vkl[i]*ukl[i]
    end);
    ########### constraint 13-14 #############
    @constraint(scnd1, sum(y[1:dt1.N["plant"]*2]) <= dt1.upl);
    @constraint(scnd1, sum(y[dt1.N["plant"]*2+1:end]) <= dt1.udc);
    return scnd1
end
scndlp = SCND_LP()
LPtime = @CPUelapsed vSolve(scndlp, TL, method=:dicho, verbose=false)
lp = getvOptData(scndlp);
weight = round(Int,mean([lp.Y_N[i][1]/lp.Y_N[i][2] for i=1:length(lp.Y_N)]))
#########################  Feasibility Pump  ###########################
function FP_Model(weight)
    model = Model(CPLEX.Optimizer); set_silent(model)
    MOI.set(model, MOI.NumberOfThreads(), 1)
    # MOI.set(model, MOI.RawParameter("CPX_PARAM_SCRIND"), false);
    @variable(model, y[1:(dt1.N["plant"]+dt1.N["distribution"])*2], Bin)
    @variable(model, uij[1:sum(dt1.Mij)], Bin);
    @variable(model, ujk[1:sum(dt1.Mjk)], Bin);
    @variable(model, ukl[1:sum(dt1.Mkl)], Bin);
    @variable( model, 0<= xij[1:sum(dt1.Mij)*5] );
    @variable( model, 0<= xjk[1:sum(dt1.Mjk)*5] );
    @variable( model, 0<= xkl[1:sum(dt1.Mkl)*5] );
    @variable( model, 0<= h[1:(dt1.N["plant"]+dt1.N["distribution"])*5*2] );
    @objective(model, Min, sum(dt1.c.*y) +sum(repeat(dt1.a[1,:], outer=sum(dt1.Mij[1,:])).*xij[1:sum(dt1.Mij[1,:])*5])+
            sum(sum(repeat(dt1.a[i,:], outer=sum(dt1.Mij[i,:])).*xij[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"])+
            sum(dt1.e.*h) + sum(dt1.gij[i]*uij[i] for i in findnz(dt1.gij)[1]) + sum(dt1.gjk[i]*ujk[i] for i in findnz(dt1.gjk)[1]) + sum(dt1.gkl[i].*ukl[i] for i in findnz(dt1.gkl)[1])+
            sum(dt1.vij.*xij)+sum(dt1.vjk.*xjk)+sum(dt1.vkl.*xkl) +
            weight*(sum(repeat(dt1.b[1,:], outer=sum(dt1.Mij[1,:])).*xij[1:sum(dt1.Mij[1,:])*5]) +
            sum(sum(repeat(dt1.b[i,:], outer=sum(dt1.Mij[i,:])).*xij[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"]) +
            sum(dt1.q.*h) + sum(dt1.rij.*xij)+sum(dt1.rjk.*xjk)+sum(dt1.rkl.*xkl))
    );
    ########## constraint 3 #############
    @constraint(model, [p=1:5], sum(xij[5*(m-1)+p] for m=1:dt1.Mij[1,1])+sum(xij[5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,1]) == sum(xjk[5*(m-1)+p] for m=1:sum(dt1.Mjk[1,:])) );
    @constraint(model, [j=2:dt1.N["plant"],p=1:5], sum(xij[5*sum(dt1.Mij[1,1:j-1])+5*(m-1)+p] for m=1:dt1.Mij[1,j])+sum(xij[5*sum(dt1.Mij[i,1:j-1])+5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,j]) == sum(xjk[sum(dt1.Mjk[1:j-1,:])*5 + 5*(m-1)+p] for m=1:sum(dt1.Mjk[j,:])) );
    @constraint(model, [p=1:5], sum(xjk[5*(m-1)+p] for m=1:dt1.Mjk[1,1])+sum(xjk[5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,1]) == sum(xkl[5*(m-1)+p] for m=1:sum(dt1.Mkl[1,:])) );
    @constraint(model, [k=2:dt1.N["distribution"],p=1:5],sum(xjk[5*sum(dt1.Mjk[1,1:k-1])+5*(m-1)+p] for m=1:dt1.Mjk[1,k])+sum(xjk[5*sum(dt1.Mjk[j,1:k-1])+5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,k]) == sum(xkl[sum(dt1.Mkl[1:k-1,:])*5 + 5*(m-1)+p] for m=1:sum(dt1.Mkl[k,:])) );
    ########### constraint 4-6 #############
    @constraint(model, [p=1:5],sum(h[2*(p-1)+t] for t=1:2) == sum(xij[5*(m-1)+p] for m=1:dt1.Mij[1,1])+sum(xij[5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,1]));
    @constraint(model, [j=2:dt1.N["plant"],p=1:5], sum(h[5*2*(j-1)+2*(p-1)+t] for t=1:2) == sum(xij[5*sum(dt1.Mij[1,1:j-1])+5*(m-1)+p] for m=1:dt1.Mij[1,j])+sum(xij[5*sum(dt1.Mij[i,1:j-1])+5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,j]) );
    @constraint(model, [p=1:5], sum(h[5*2*dt1.N["plant"]+2*(p-1)+t] for t=1:2) == sum(xjk[5*(m-1)+p] for m=1:dt1.Mjk[1,1])+sum(xjk[5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,1]) );
    @constraint(model, [k=2:dt1.N["distribution"],p=1:5], sum(h[5*2*dt1.N["plant"]+5*2*(k-1)+2*(p-1)+t] for t=1:2) == sum(xjk[5*sum(dt1.Mjk[1,1:k-1])+5*(m-1)+p] for m=1:dt1.Mjk[1,k])+sum(xjk[5*sum(dt1.Mjk[j,1:k-1])+5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,k]));
    @constraint(model, [p=1:5], sum(xkl[5*(m-1)+p] for m=1:dt1.Mkl[1,1]) +sum(xkl[5*(m-1)+p+(5*sum(dt1.Mkl[1:k-1,:]))] for k=2:dt1.N["distribution"] for m=1:dt1.Mkl[k,1]) >= dt1.d[1,p]);
    @constraint(model, [l=2:dt1.N["customer"], p=1:5], sum(xkl[sum(dt1.Mkl[1,1:l-1])*5 + 5*(m-1)+p] for m=1:dt1.Mkl[1,l])+ sum(xkl[5*sum(dt1.Mkl[1:k-1,:])+5*sum(dt1.Mkl[k,1:l-1])+5*(m-1)+p] for k=2:dt1.N["distribution"] for m=1:dt1.Mkl[k,l]) >= dt1.d[l,p]);
    ########### constraint 7 #############
    @constraint(model, sum(xij[1:5*sum(dt1.Mij[1,:])]) <= dt1.N["cas"][1]);
    @constraint(model, [i=2:dt1.N["supplier"]],  sum(xij[5*sum(dt1.Mij[1:i-1,:])+1:5*sum(dt1.Mij[1:i,:])]) <= dt1.N["cas"][i]);
    ########### constraint 8 #############
    @constraint(model,[j=1:dt1.N["plant"]+dt1.N["distribution"], t=1:2], sum(h[5*2*(j-1)+((p-1)*2)+t] for p=1:5) <= [dt1.N["cap"];dt1.N["cad"]][j]*y[2*(j-1)+t]);
    ########### constraint 9 #############
    @constraint(model,[j=1:dt1.N["plant"]+dt1.N["distribution"]], sum(y[2*(j-1)+1:2*(j-1)+2]) <= 1);
    ########### constraint 10 #############
    @constraints(model,begin
        sum(uij[1:dt1.Mij[1,1]]) <= 1
        sum(uij[sum(dt1.Mij[1,:])+dt1.Mij[2,1]]) <= 1
        [j=2:dt1.N["plant"]], sum(uij[sum(dt1.Mij[1,1:j-1])+1:sum(dt1.Mij[1,1:j-1])+dt1.Mij[1,j]]) <= 1
        [i=2:dt1.N["supplier"],j=2:dt1.N["plant"]],  sum(uij[sum(dt1.Mij[1:i-1,:])+sum(dt1.Mij[i,1:j-1])+1:sum(dt1.Mij[1:i-1,:])+sum(dt1.Mij[i,1:j-1])+dt1.Mij[i,j]])<= 1
        sum(ujk[1:dt1.Mjk[1,1]]) <= 1
        sum(ujk[sum(dt1.Mjk[1,:])+dt1.Mjk[2,1]]) <= 1
        [k=2:dt1.N["distribution"]], sum(ujk[sum(dt1.Mjk[1,1:k-1])+1:sum(dt1.Mjk[1,1:k-1])+dt1.Mjk[1,k]]) <= 1
        [j=2:dt1.N["plant"],k=2:dt1.N["distribution"]],  sum(ujk[sum(dt1.Mjk[1:j-1,:])+sum(dt1.Mjk[j,1:j-1])+1:sum(dt1.Mjk[1:j-1,:])+sum(dt1.Mjk[j,1:j-1])+dt1.Mjk[j,k]]) <= 1
        sum(ukl[1:dt1.Mkl[1,1]]) <= 1 #[1,1]
        sum(ukl[sum(dt1.Mkl[1,:])+dt1.Mkl[2,1]]) <= 1 #[2,1]
        [l=2:dt1.N["customer"]], sum(ukl[sum(dt1.Mkl[1,1:l-1])+1:sum(dt1.Mkl[1,1:l-1])+dt1.Mkl[1,l]]) <= 1
        [k=2:dt1.N["distribution"],l=2:dt1.N["customer"]],  sum(ukl[sum(dt1.Mkl[1:k-1,:])+sum(dt1.Mkl[k,1:l-1])+1:sum(dt1.Mkl[1:k-1,:])+sum(dt1.Mkl[k,1:l-1])+dt1.Mkl[k,l]])<= 1
    end);
    ########### constraint 11 #############
    @constraints(model, begin
        [i=1:sum(dt1.Mij)], sum(xij[5*(i-1)+1:5*i]) <= dt1.bigM*uij[i]
        [i=1:sum(dt1.Mjk)], sum(xjk[5*(i-1)+1:5*i]) <= dt1.bigM*ujk[i]
        [i=1:sum(dt1.Mkl)], sum(xkl[5*(i-1)+1:5*i]) <= dt1.bigM*ukl[i]
    end);
    ########### constraint 12 #############
    @constraints(model, begin
            [i in findnz(dt1.Vij)[1]], sum(xij[5*(i-1)+1:5*i]) >= dt1.Vij[i]*uij[i]
            [i in findnz(dt1.Vjk)[1]], sum(xjk[5*(i-1)+1:5*i]) >= dt1.Vjk[i]*ujk[i]
            # [i in findnz(dt1.Vkl)[1]], sum(xkl[5*(i-1)+1:5*i]) >= dt1.Vkl[i]*ukl[i]
    end);
    ########### constraint 13-14 #############
    @constraint(model, sum(y[1:dt1.N["plant"]*2]) <= dt1.upl);
    @constraint(model, sum(y[dt1.N["plant"]*2+1:end]) <= dt1.udc);
    return model
end
function LP_Model(weight)
    lp = Model(CPLEX.Optimizer); set_silent(lp)
    MOI.set(lp, MOI.NumberOfThreads(), 1)
    # MOI.set(lp, MOI.RawParameter("CPX_PARAM_SCRIND"), false);
    @variable(lp, 0 <= y[1:(dt1.N["plant"]+dt1.N["distribution"])*2] <= 1)
    @variable(lp, 0 <= uij[1:sum(dt1.Mij)] <= 1);
    @variable(lp, 0 <= ujk[1:sum(dt1.Mjk)] <= 1);
    @variable(lp, 0 <= ukl[1:sum(dt1.Mkl)] <= 1);
    @variable( lp, 0 <= xij[1:sum(dt1.Mij)*5] );
    @variable( lp, 0 <= xjk[1:sum(dt1.Mjk)*5] );
    @variable( lp, 0 <= xkl[1:sum(dt1.Mkl)*5] );
    @variable( lp, 0 <= h[1:(dt1.N["plant"]+dt1.N["distribution"])*5*2] );
    @objective(lp, Min, sum(dt1.c.*y) +sum(repeat(dt1.a[1,:], outer=sum(dt1.Mij[1,:])).*xij[1:sum(dt1.Mij[1,:])*5])+
            sum(sum(repeat(dt1.a[i,:], outer=sum(dt1.Mij[i,:])).*xij[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"])+
            sum(dt1.e.*h) + sum(dt1.gij[i]*uij[i] for i in findnz(dt1.gij)[1]) + sum(dt1.gjk[i]*ujk[i] for i in findnz(dt1.gjk)[1]) + sum(dt1.gkl[i].*ukl[i] for i in findnz(dt1.gkl)[1])+
            sum(dt1.vij.*xij)+sum(dt1.vjk.*xjk)+sum(dt1.vkl.*xkl) +
            weight*(sum(repeat(dt1.b[1,:], outer=sum(dt1.Mij[1,:])).*xij[1:sum(dt1.Mij[1,:])*5]) +
            sum(sum(repeat(dt1.b[i,:], outer=sum(dt1.Mij[i,:])).*xij[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"]) +
            sum(dt1.q.*h) + sum(dt1.rij.*xij)+sum(dt1.rjk.*xjk)+sum(dt1.rkl.*xkl))
    );
    ########## constraint 3 #############
    @constraint(lp, [p=1:5], sum(xij[5*(m-1)+p] for m=1:dt1.Mij[1,1])+sum(xij[5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,1]) == sum(xjk[5*(m-1)+p] for m=1:sum(dt1.Mjk[1,:])) );
    @constraint(lp, [j=2:dt1.N["plant"],p=1:5], sum(xij[5*sum(dt1.Mij[1,1:j-1])+5*(m-1)+p] for m=1:dt1.Mij[1,j])+sum(xij[5*sum(dt1.Mij[i,1:j-1])+5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,j]) == sum(xjk[sum(dt1.Mjk[1:j-1,:])*5 + 5*(m-1)+p] for m=1:sum(dt1.Mjk[j,:])) );
    @constraint(lp, [p=1:5], sum(xjk[5*(m-1)+p] for m=1:dt1.Mjk[1,1])+sum(xjk[5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,1]) == sum(xkl[5*(m-1)+p] for m=1:sum(dt1.Mkl[1,:])) );
    @constraint(lp, [k=2:dt1.N["distribution"],p=1:5],sum(xjk[5*sum(dt1.Mjk[1,1:k-1])+5*(m-1)+p] for m=1:dt1.Mjk[1,k])+sum(xjk[5*sum(dt1.Mjk[j,1:k-1])+5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,k]) == sum(xkl[sum(dt1.Mkl[1:k-1,:])*5 + 5*(m-1)+p] for m=1:sum(dt1.Mkl[k,:])) );
    ########### constraint 4-6 #############
    @constraint(lp, [p=1:5],sum(h[2*(p-1)+t] for t=1:2) == sum(xij[5*(m-1)+p] for m=1:dt1.Mij[1,1])+sum(xij[5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,1]));
    @constraint(lp, [j=2:dt1.N["plant"],p=1:5], sum(h[5*2*(j-1)+2*(p-1)+t] for t=1:2) == sum(xij[5*sum(dt1.Mij[1,1:j-1])+5*(m-1)+p] for m=1:dt1.Mij[1,j])+sum(xij[5*sum(dt1.Mij[i,1:j-1])+5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,j]) );
    @constraint(lp, [p=1:5], sum(h[5*2*dt1.N["plant"]+2*(p-1)+t] for t=1:2) == sum(xjk[5*(m-1)+p] for m=1:dt1.Mjk[1,1])+sum(xjk[5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,1]) );
    @constraint(lp, [k=2:dt1.N["distribution"],p=1:5], sum(h[5*2*dt1.N["plant"]+5*2*(k-1)+2*(p-1)+t] for t=1:2) == sum(xjk[5*sum(dt1.Mjk[1,1:k-1])+5*(m-1)+p] for m=1:dt1.Mjk[1,k])+sum(xjk[5*sum(dt1.Mjk[j,1:k-1])+5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,k]));
    @constraint(lp, [p=1:5], sum(xkl[5*(m-1)+p] for m=1:dt1.Mkl[1,1]) +sum(xkl[5*(m-1)+p+(5*sum(dt1.Mkl[1:k-1,:]))] for k=2:dt1.N["distribution"] for m=1:dt1.Mkl[k,1]) >= dt1.d[1,p]);
    @constraint(lp, [l=2:dt1.N["customer"], p=1:5], sum(xkl[sum(dt1.Mkl[1,1:l-1])*5 + 5*(m-1)+p] for m=1:dt1.Mkl[1,l])+ sum(xkl[5*sum(dt1.Mkl[1:k-1,:])+5*sum(dt1.Mkl[k,1:l-1])+5*(m-1)+p] for k=2:dt1.N["distribution"] for m=1:dt1.Mkl[k,l]) >= dt1.d[l,p]);
    ########### constraint 7 #############
    @constraint(lp, sum(xij[1:5*sum(dt1.Mij[1,:])]) <= dt1.N["cas"][1]);
    @constraint(lp, [i=2:dt1.N["supplier"]],  sum(xij[5*sum(dt1.Mij[1:i-1,:])+1:5*sum(dt1.Mij[1:i,:])]) <= dt1.N["cas"][i]);
    ########### constraint 8 #############
    @constraint(lp,[j=1:dt1.N["plant"]+dt1.N["distribution"], t=1:2], sum(h[5*2*(j-1)+((p-1)*2)+t] for p=1:5) <= [dt1.N["cap"];dt1.N["cad"]][j]*y[2*(j-1)+t]);
    ########### constraint 9 #############
    @constraint(lp,[j=1:dt1.N["plant"]+dt1.N["distribution"]], sum(y[2*(j-1)+1:2*(j-1)+2]) <= 1);
    ########### constraint 10 #############
    @constraints(lp,begin
            sum(uij[1:dt1.Mij[1,1]]) <= 1
            sum(uij[sum(dt1.Mij[1,:])+dt1.Mij[2,1]]) <= 1
            [j=2:dt1.N["plant"]], sum(uij[sum(dt1.Mij[1,1:j-1])+1:sum(dt1.Mij[1,1:j-1])+dt1.Mij[1,j]]) <= 1
            [i=2:dt1.N["supplier"],j=2:dt1.N["plant"]],  sum(uij[sum(dt1.Mij[1:i-1,:])+sum(dt1.Mij[i,1:j-1])+1:sum(dt1.Mij[1:i-1,:])+sum(dt1.Mij[i,1:j-1])+dt1.Mij[i,j]])<= 1
            sum(ujk[1:dt1.Mjk[1,1]]) <= 1
            sum(ujk[sum(dt1.Mjk[1,:])+dt1.Mjk[2,1]]) <= 1
            [k=2:dt1.N["distribution"]], sum(ujk[sum(dt1.Mjk[1,1:k-1])+1:sum(dt1.Mjk[1,1:k-1])+dt1.Mjk[1,k]]) <= 1
            [j=2:dt1.N["plant"],k=2:dt1.N["distribution"]],  sum(ujk[sum(dt1.Mjk[1:j-1,:])+sum(dt1.Mjk[j,1:j-1])+1:sum(dt1.Mjk[1:j-1,:])+sum(dt1.Mjk[j,1:j-1])+dt1.Mjk[j,k]]) <= 1
            sum(ukl[1:dt1.Mkl[1,1]]) <= 1 #[1,1]
            sum(ukl[sum(dt1.Mkl[1,:])+dt1.Mkl[2,1]]) <= 1 #[2,1]
            [l=2:dt1.N["customer"]], sum(ukl[sum(dt1.Mkl[1,1:l-1])+1:sum(dt1.Mkl[1,1:l-1])+dt1.Mkl[1,l]]) <= 1
            [k=2:dt1.N["distribution"],l=2:dt1.N["customer"]],  sum(ukl[sum(dt1.Mkl[1:k-1,:])+sum(dt1.Mkl[k,1:l-1])+1:sum(dt1.Mkl[1:k-1,:])+sum(dt1.Mkl[k,1:l-1])+dt1.Mkl[k,l]])<= 1
    end);
    ########### constraint 11 #############
    @constraints(lp, begin
        [i=1:sum(dt1.Mij)], sum(xij[5*(i-1)+1:5*i]) <= dt1.bigM*uij[i]
        [i=1:sum(dt1.Mjk)], sum(xjk[5*(i-1)+1:5*i]) <= dt1.bigM*ujk[i]
        [i=1:sum(dt1.Mkl)], sum(xkl[5*(i-1)+1:5*i]) <= dt1.bigM*ukl[i]
    end);
    ########### constraint 12 #############
    @constraints(lp, begin
            [i in findnz(dt1.Vij)[1]], sum(xij[5*(i-1)+1:5*i]) >= dt1.Vij[i]*uij[i]
            [i in findnz(dt1.Vjk)[1]], sum(xjk[5*(i-1)+1:5*i]) >= dt1.Vjk[i]*ujk[i]
            # [i in findnz(dt1.Vkl)[1]], sum(xkl[5*(i-1)+1:5*i]) >= dt1.Vkl[i]*ukl[i]
    end);
    ########### constraint 13-14 #############
    @constraint(lp, sum(y[1:dt1.N["plant"]*2]) <= dt1.upl);
    @constraint(lp, sum(y[dt1.N["plant"]*2+1:end]) <= dt1.udc);
    return lp
end
fbmodel = FP_Model(weight); dist = LP_Model(weight)

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
function FP_FBcheck(model,yr)#,u1r,u2r,u3r)
    JuMP.fix.(model[:y],yr; force=true)
    # JuMP.fix.(model[:uij],u1r; force=true)
    # JuMP.fix.(model[:ujk],u2r; force=true)
    # JuMP.fix.(model[:ukl],u3r; force=true)
    optimize!(model)
    # @show
    st = termination_status(model)
    if st == MOI.OPTIMAL
        return true
    else
        return false
    end
end
function fbsearch(yr)#,u1r,u2r,u3r) #solveLP
    idy_0 = findall(k->k==0, yr)
    idy_1 = findall(k->k==1, yr)
    # idu1_0 = findall(k->k==0, u1r)
    # idu1_1 = findall(k->k==1, u1r)
    # idu2_0 = findall(k->k==0, u2r)
    # idu2_1 = findall(k->k==1, u2r)
    # idu3_0 = findall(k->k==0, u3r)
    # idu3_1 = findall(k->k==1, u3r)
    @objective( dist, Min, sum(dist[:y][i] for i in idy_0) + sum(1-(dist[:y][j]) for j in idy_1))
    # +
    #     sum(dist[:uij][i] for i in idu1_0) + sum(1-(dist[:uij][j]) for j in idu1_1)+
    #     sum(dist[:ujk][i] for i in idu2_0) + sum(1-(dist[:ujk][j]) for j in idu2_1)+
    #     sum(dist[:ukl][i] for i in idu3_0) + sum(1-(dist[:ukl][j]) for j in idu3_1))
    optimize!(dist)
    if termination_status(dist) == MOI.OPTIMAL
        return JuMP.value.(dist[:y])
        # ,JuMP.value.(dist[:uij]),JuMP.value.(dist[:ujk]),JuMP.value.(dist[:ukl])
    else
        return 0#,0,0,0
    end
end
len = [length(scndlp[:y]),length(scndlp[:uij]),length(scndlp[:ujk]),length(scndlp[:ukl]),length(scndlp[:xij]),length(scndlp[:xjk]),length(scndlp[:xkl]),length(scndlp[:h])]
# 𝚯 = [0,ℯ/(ℯ+ℯ^2),ℯ^2/(ℯ+ℯ^2)];
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
function getobjval(x)
    y = x[1:len[1]];
    uij = x[1+len[1]:len[1]+len[2]];
    ujk = x[1+sum(len[i] for i=1:2):sum(len[i] for i=1:3)]
    ukl = x[1+sum(len[i] for i=1:3):sum(len[i] for i=1:4)]
    xij = x[1+sum(len[i] for i=1:4):sum(len[i] for i=1:5)]
    xjk = x[1+sum(len[i] for i=1:5):sum(len[i] for i=1:6)]
    xkl = x[1+sum(len[i] for i=1:6):sum(len[i] for i=1:7)]
    h = x[1+sum(len[i] for i=1:7):sum(len[i] for i=1:8)]

    obj1 =  sum(dt1.c.*y)+sum(repeat(dt1.a[1,:], outer=sum(dt1.Mij[1,:])).*xij[1:sum(dt1.Mij[1,:])*5])+
            sum(sum(repeat(dt1.a[i,:], outer=sum(dt1.Mij[i,:])).*xij[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"])+
            sum(dt1.e.*h) + sum(dt1.gij[i]*uij[i] for i in findnz(dt1.gij)[1]) + sum(dt1.gjk[i]*ujk[i] for i in findnz(dt1.gjk)[1]) + sum(dt1.gkl[i].*ukl[i] for i in findnz(dt1.gkl)[1])+
            sum(dt1.vij.*xij)+sum(dt1.vjk.*xjk)+sum(dt1.vkl.*xkl)
    obj2 = sum(repeat(dt1.b[1,:], outer=sum(dt1.Mij[1,:])).*xij[1:sum(dt1.Mij[1,:])*5]) +
            sum(sum(repeat(dt1.b[i,:], outer=sum(dt1.Mij[i,:])).*xij[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"]) +
            sum(dt1.q.*h) + sum(dt1.rij.*xij)+sum(dt1.rjk.*xjk)+sum(dt1.rkl.*xkl)
    return [obj1,obj2]
end

function FP(candX,len,TL)
    X = []; PF =[]; Tabu = []; newsol = 0; candlist = copy(candX)
    Y = []; t0 = time();
    # U1 = []; U2= []; U3 = [];
    while candlist != [] && time() - t0 < TL
        k = rand(1:length(candlist))
        x_t = candlist[k]
        yt = x_t[1:len[1]]; yr = copy(yt)
        # u1t = x_t[1+len[1]:len[1]+len[2]];
        # u2t = x_t[1+sum(len[i] for i=1:2):sum(len[i] for i=1:3)]
        # u3t = x_t[1+sum(len[i] for i=1:3):sum(len[i] for i=1:4)]
        SearchDone = false; iter=0;
        Max_iter =  5#length(findall(i-> 0<i<1,yt))
        while iter < Max_iter && SearchDone == false && time() - t0 < TL
            # x_r = round.(Int,x_t);
            yid = findall(p->p>0.2,yt);
            for j=1:len[1]
                if j in yid
                    yr[j]=1
                else
                    yr[j]=0
                end
            end
            if FP_FBcheck(fbmodel,yr) == true #,u1r,u2r,u3r)
                sol = value.(all_variables(fbmodel)); ndp = getobjval(sol)
                if sol ∉ X  && dominated(ndp,collect(values(PF)))==false
                    push!(X,sol); push!(PF,ndp) #PF[k] = ndp
                    # u1r = sol[1+len[1]:len[1]+len[2]];
                    # u2r = sol[1+sum(len[i] for i=1:2):sum(len[i] for i=1:3)]
                    # u3r = sol[1+sum(len[i] for i=1:3):sum(len[i] for i=1:4)]
                    push!(Y,yr);
                    # push!(U1,u1r); push!(U2,u2r); push!(U3,u3r);
                    newsol+=1; SearchDone = true
                    deleteat!(candlist, k);
                    # println("rounding worked")
                end
            else
                if yr ∈ Tabu
                    yr = flipoper(Y,yt,yr); # u1r = flipoper(U1,u1t,u1r); u2r = flipoper(U2,u2t,u2r); u3r = flipoper(U3,u3t,u3r)
                    if yr == [] # if any(i->i==[], [yr,u1r,u2r,u3r])
                        SearchDone = true;
                        # println("flip failed")
                    else
                        if FP_FBcheck(fbmodel,yr) == true #,u1r,u2r,u3r)
                            sol = value.(all_variables(fbmodel)); ndp = getobjval(sol)
                            if sol ∉ X && dominated(ndp,collect(values(PF)))==false
                                push!(X,sol); push!(PF,ndp) #PF[k] = ndp
                                # u1r = sol[1+len[1]:len[1]+len[2]];
                                # u2r = sol[1+sum(len[i] for i=1:2):sum(len[i] for i=1:3)]
                                # u3r = sol[1+sum(len[i] for i=1:3):sum(len[i] for i=1:4)]
                                push!(Y,yr);
                                # push!(U1,u1r); push!(U2,u2r); push!(U3,u3r);
                                newsol+=1; SearchDone = true;
                                deleteat!(candlist, k);
                                # println("flip worked")
                            end
                        end
                    end
                end
                if time()-t0 >= TL
                    break
                end
                if SearchDone == false
                    push!(Tabu,yr)
                    yt = fbsearch(yr)
                    if yt==0  #when there's no new feasible lp sol
                        deleteat!(candlist, k);
                        # println("no solution")
                        SearchDone = true
                    end
                end
            end
			iter+=1
        end
    end
    return X,PF,newsol,candlist
end
optimize!(fbmodel); optimize!(dist);
FPtime = @CPUelapsed lx,ly,ln,candlist = FP(lp.X_E,len,round(Int,LPtime*2))
function Postpro(P,Pobj)
    copysol = Dict(); copyobj = Dict();
    for i=1:length(Pobj)
        copysol[i] = P[i]
        copyobj[i] = Pobj[i]
    end
    for i=1:length(Pobj)-1
        for j=i+1:length(Pobj)
            if all(Pobj[i] .>= Pobj[j]) == true #dominated by PF[j]
                copyobj[i]=nothing; copysol[i]=nothing; break
            elseif all(Pobj[j] .>= Pobj[i]) == true
                copyobj[j]=nothing; copysol[j]=nothing;
            end
        end
    end

    finalsol = filter!(a->a!=nothing, collect(values(copysol)))
    finalobj = filter!(a->a!=nothing, collect(values(copyobj)))

    return finalsol,finalobj
end
fx,fy = Postpro(lx,ly)

function PR_Model(weight)
    model = Model(CPLEX.Optimizer); set_silent(model)
    MOI.set(model, MOI.NumberOfThreads(), 1)
    # MOI.set(model, MOI.RawParameter("CPX_PARAM_SCRIND"), false);
    @variable(model, y[1:(dt1.N["plant"]+dt1.N["distribution"])*2], Bin)
    @variable(model, uij[1:sum(dt1.Mij)], Bin);
    @variable(model, ujk[1:sum(dt1.Mjk)], Bin);
    @variable(model, ukl[1:sum(dt1.Mkl)], Bin);
    @variable( model, 0<= xij[1:sum(dt1.Mij)*5] );
    @variable( model, 0<= xjk[1:sum(dt1.Mjk)*5] );
    @variable( model, 0<= xkl[1:sum(dt1.Mkl)*5] );
    @variable( model, 0<= h[1:(dt1.N["plant"]+dt1.N["distribution"])*5*2] );
    @objective(model, Min, sum(dt1.c.*y) +sum(repeat(dt1.a[1,:], outer=sum(dt1.Mij[1,:])).*xij[1:sum(dt1.Mij[1,:])*5])+
            sum(sum(repeat(dt1.a[i,:], outer=sum(dt1.Mij[i,:])).*xij[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"])+
            sum(dt1.e.*h) + sum(dt1.gij[i]*uij[i] for i in findnz(dt1.gij)[1]) + sum(dt1.gjk[i]*ujk[i] for i in findnz(dt1.gjk)[1]) + sum(dt1.gkl[i].*ukl[i] for i in findnz(dt1.gkl)[1])+
            sum(dt1.vij.*xij)+sum(dt1.vjk.*xjk)+sum(dt1.vkl.*xkl) +
            weight*(sum(repeat(dt1.b[1,:], outer=sum(dt1.Mij[1,:])).*xij[1:sum(dt1.Mij[1,:])*5]) +
            sum(sum(repeat(dt1.b[i,:], outer=sum(dt1.Mij[i,:])).*xij[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"]) +
            sum(dt1.q.*h) + sum(dt1.rij.*xij)+sum(dt1.rjk.*xjk)+sum(dt1.rkl.*xkl))
    );
    ########## constraint 3 #############
    @constraint(model, [p=1:5], sum(xij[5*(m-1)+p] for m=1:dt1.Mij[1,1])+sum(xij[5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,1]) == sum(xjk[5*(m-1)+p] for m=1:sum(dt1.Mjk[1,:])) );
    @constraint(model, [j=2:dt1.N["plant"],p=1:5], sum(xij[5*sum(dt1.Mij[1,1:j-1])+5*(m-1)+p] for m=1:dt1.Mij[1,j])+sum(xij[5*sum(dt1.Mij[i,1:j-1])+5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,j]) == sum(xjk[sum(dt1.Mjk[1:j-1,:])*5 + 5*(m-1)+p] for m=1:sum(dt1.Mjk[j,:])) );
    @constraint(model, [p=1:5], sum(xjk[5*(m-1)+p] for m=1:dt1.Mjk[1,1])+sum(xjk[5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,1]) == sum(xkl[5*(m-1)+p] for m=1:sum(dt1.Mkl[1,:])) );
    @constraint(model, [k=2:dt1.N["distribution"],p=1:5],sum(xjk[5*sum(dt1.Mjk[1,1:k-1])+5*(m-1)+p] for m=1:dt1.Mjk[1,k])+sum(xjk[5*sum(dt1.Mjk[j,1:k-1])+5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,k]) == sum(xkl[sum(dt1.Mkl[1:k-1,:])*5 + 5*(m-1)+p] for m=1:sum(dt1.Mkl[k,:])) );
    ########### constraint 4-6 #############
    @constraint(model, [p=1:5],sum(h[2*(p-1)+t] for t=1:2) == sum(xij[5*(m-1)+p] for m=1:dt1.Mij[1,1])+sum(xij[5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,1]));
    @constraint(model, [j=2:dt1.N["plant"],p=1:5], sum(h[5*2*(j-1)+2*(p-1)+t] for t=1:2) == sum(xij[5*sum(dt1.Mij[1,1:j-1])+5*(m-1)+p] for m=1:dt1.Mij[1,j])+sum(xij[5*sum(dt1.Mij[i,1:j-1])+5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,j]) );
    @constraint(model, [p=1:5], sum(h[5*2*dt1.N["plant"]+2*(p-1)+t] for t=1:2) == sum(xjk[5*(m-1)+p] for m=1:dt1.Mjk[1,1])+sum(xjk[5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,1]) );
    @constraint(model, [k=2:dt1.N["distribution"],p=1:5], sum(h[5*2*dt1.N["plant"]+5*2*(k-1)+2*(p-1)+t] for t=1:2) == sum(xjk[5*sum(dt1.Mjk[1,1:k-1])+5*(m-1)+p] for m=1:dt1.Mjk[1,k])+sum(xjk[5*sum(dt1.Mjk[j,1:k-1])+5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,k]));
    @constraint(model, [p=1:5], sum(xkl[5*(m-1)+p] for m=1:dt1.Mkl[1,1]) +sum(xkl[5*(m-1)+p+(5*sum(dt1.Mkl[1:k-1,:]))] for k=2:dt1.N["distribution"] for m=1:dt1.Mkl[k,1]) >= dt1.d[1,p]);
    @constraint(model, [l=2:dt1.N["customer"], p=1:5], sum(xkl[sum(dt1.Mkl[1,1:l-1])*5 + 5*(m-1)+p] for m=1:dt1.Mkl[1,l])+ sum(xkl[5*sum(dt1.Mkl[1:k-1,:])+5*sum(dt1.Mkl[k,1:l-1])+5*(m-1)+p] for k=2:dt1.N["distribution"] for m=1:dt1.Mkl[k,l]) >= dt1.d[l,p]);
    ########### constraint 7 #############
    @constraint(model, sum(xij[1:5*sum(dt1.Mij[1,:])]) <= dt1.N["cas"][1]);
    @constraint(model, [i=2:dt1.N["supplier"]],  sum(xij[5*sum(dt1.Mij[1:i-1,:])+1:5*sum(dt1.Mij[1:i,:])]) <= dt1.N["cas"][i]);
    ########### constraint 8 #############
    @constraint(model,[j=1:dt1.N["plant"]+dt1.N["distribution"], t=1:2], sum(h[5*2*(j-1)+((p-1)*2)+t] for p=1:5) <= [dt1.N["cap"];dt1.N["cad"]][j]*y[2*(j-1)+t]);
    ########### constraint 9 #############
    @constraint(model,[j=1:dt1.N["plant"]+dt1.N["distribution"]], sum(y[2*(j-1)+1:2*(j-1)+2]) <= 1);
    ########### constraint 10 #############
    @constraints(model,begin
        sum(uij[1:dt1.Mij[1,1]]) <= 1
        sum(uij[sum(dt1.Mij[1,:])+dt1.Mij[2,1]]) <= 1
        [j=2:dt1.N["plant"]], sum(uij[sum(dt1.Mij[1,1:j-1])+1:sum(dt1.Mij[1,1:j-1])+dt1.Mij[1,j]]) <= 1
        [i=2:dt1.N["supplier"],j=2:dt1.N["plant"]],  sum(uij[sum(dt1.Mij[1:i-1,:])+sum(dt1.Mij[i,1:j-1])+1:sum(dt1.Mij[1:i-1,:])+sum(dt1.Mij[i,1:j-1])+dt1.Mij[i,j]])<= 1
        sum(ujk[1:dt1.Mjk[1,1]]) <= 1
        sum(ujk[sum(dt1.Mjk[1,:])+dt1.Mjk[2,1]]) <= 1
        [k=2:dt1.N["distribution"]], sum(ujk[sum(dt1.Mjk[1,1:k-1])+1:sum(dt1.Mjk[1,1:k-1])+dt1.Mjk[1,k]]) <= 1
        [j=2:dt1.N["plant"],k=2:dt1.N["distribution"]],  sum(ujk[sum(dt1.Mjk[1:j-1,:])+sum(dt1.Mjk[j,1:j-1])+1:sum(dt1.Mjk[1:j-1,:])+sum(dt1.Mjk[j,1:j-1])+dt1.Mjk[j,k]]) <= 1
        sum(ukl[1:dt1.Mkl[1,1]]) <= 1 #[1,1]
        sum(ukl[sum(dt1.Mkl[1,:])+dt1.Mkl[2,1]]) <= 1 #[2,1]
        [l=2:dt1.N["customer"]], sum(ukl[sum(dt1.Mkl[1,1:l-1])+1:sum(dt1.Mkl[1,1:l-1])+dt1.Mkl[1,l]]) <= 1
        [k=2:dt1.N["distribution"],l=2:dt1.N["customer"]],  sum(ukl[sum(dt1.Mkl[1:k-1,:])+sum(dt1.Mkl[k,1:l-1])+1:sum(dt1.Mkl[1:k-1,:])+sum(dt1.Mkl[k,1:l-1])+dt1.Mkl[k,l]])<= 1
    end);
    ########### constraint 11 #############
    @constraints(model, begin
        [i=1:sum(dt1.Mij)], sum(xij[5*(i-1)+1:5*i]) <= dt1.bigM*uij[i]
        [i=1:sum(dt1.Mjk)], sum(xjk[5*(i-1)+1:5*i]) <= dt1.bigM*ujk[i]
        [i=1:sum(dt1.Mkl)], sum(xkl[5*(i-1)+1:5*i]) <= dt1.bigM*ukl[i]
    end);
    ########### constraint 12 #############
    @constraints(model, begin
            [i in findnz(dt1.Vij)[1]], sum(xij[5*(i-1)+1:5*i]) >= dt1.Vij[i]*uij[i]
            [i in findnz(dt1.Vjk)[1]], sum(xjk[5*(i-1)+1:5*i]) >= dt1.Vjk[i]*ujk[i]
            # [i in findnz(dt1.Vkl)[1]], sum(xkl[5*(i-1)+1:5*i]) >= dt1.Vkl[i]*ukl[i]
    end);
    ########### constraint 13-14 #############
    @constraint(model, sum(y[1:dt1.N["plant"]*2]) <= dt1.upl);
    @constraint(model, sum(y[dt1.N["plant"]*2+1:end]) <= dt1.udc);
    return model
end
function createNB(SI,dif,exploredSI)
    neibour = []; #neiobj = [];
    for i in dif
        cpSI = copy(SI)
        # if cpSI[i] == round(cpSI[i]) #if vari is int value
        if cpSI[i] == 1
            cpSI[i] = 0
        else
            cpSI[i] = 1
        end
        push!(neibour, cpSI);# push!(neiobj, binobj(cpSI,bvar))
        # else #if vari is fractional
        #     cpSI[i] = 1; push!(neibour, cpSI); push!( neiobj, binobj(cpSI,bvar) )
        #     cpSI = copy(SI)
        #     cpSI[i] = 0; push!(neibour, cpSI); push!( neiobj, binobj(cpSI,bvar) )
        # end
    end
    idx = findall(i-> neibour[i] in exploredSI, 1:length(neibour))
    neibour1 = setdiff(neibour,exploredSI)
    neibour2 = StatsBase.sample(neibour1, round(Int,length(neibour1)*0.1) , replace=false)
    # neiobj = neiobj[setdiff(1:end, idx),:]
    return neibour2#,neiobj
end
function nextSI(neibour,SI)
    SIobj = getobjval(SI)
    neiobj = [getobjval(neibour[i]) for i=1:length(neibour)]
    for i=1:length(neiobj)
        if length(neibour) == 1  #if there is one candiate sol
            return neibour[1]#,neibour[1]
        else length(neiobj) > 1 # if there are multiple candiates, check the improved ratio
            ratiotb = zeros(length(neiobj),length(neiobj[1]))
            for i=1:length(neiobj)
                ratiotb[i,:] = neiobj[i]./SIobj
            end
            ranktb = zeros(length(neiobj),length(neiobj[1]))
            for i=1:length(neiobj[1])
                ranktb[:,i] = tiedrank(ratiotb[:,i])
            end
            ranksum = [sum(ranktb[i,:]) for i=1:length(neiobj)]
            mostimp = findall(x-> x == maximum(ranksum), ranksum)
            k = rand(mostimp)[1]
            return neibour[k]#,neiobj[k]
        end
    end
end
# function First_FBcheck(neibour)
#     #constraint 13 & 14
#     neibour1 = filter!(y->(sum(y[1:dt1.N["plant"]*2])<=dt1.upl && sum(y[1:dt1.N["plant"]*2])>0) || (sum(y[1+dt1.N["plant"]*2:len[1]])<=dt1.udc && sum(y[1+dt1.N["plant"]*2:len[1]])>0), neibour)
#     #constraint 9
#     for j=1:dt1.N["plant"]+dt1.N["distribution"]
#         filter!(y-> sum(y[2*(j-1)+1:2*(j-1)+2]) <= 1, neibour1)
#     end
#     return neibour1
# end
function PR_FBcheck(model,yr)#,u1r,u2r,u3r)
    JuMP.fix.(model[:y],yr; force=true)
    # JuMP.fix.(model[:uij],u1r; force=true)
    # JuMP.fix.(model[:ujk],u2r; force=true)
    # JuMP.fix.(model[:ukl],u3r; force=true)
    optimize!(model)
    if termination_status(model) == MOI.OPTIMAL
        return true
    else
        return false
    end
end
prmodel = PR_Model(weight); optimize!(prmodel);
function PR(X,Y,len,TL)
    candX = copy(X); candY = copy(Y); bvar = sum(len[i] for i=1:4);
    IGPair=[]; exploredSI = []; t0=time(); iter=0; newsol=0;
    while time()-t0 < TL && length(IGPair)<(length(candY)*(length(candY)-1))
        @label NewIter
	    I,G = StatsBase.sample(1:length(candX), 2, replace=false);
        SI = candX[I]; SG = candX[G];
        # SI_r = round.(SI); SG_r = round.(SG)
        dif = findall(i-> SI[1:len[1]][i]!=SG[1:len[1]][i], 1:len[1])
         #Max_iter = 20;
        while length(dif)>0 && [I,G]∉IGPair && (time()-t0<TL) #&& iter<Max_iter
            neibour = createNB(SI[1:len[1]],dif,exploredSI)
            # println("# of neighbours: ", length(neibour2))
            if (length(neibour)==0) #(time()-t0 >= TL)
                @goto NewIter
            else
                candSI =[]
                # l=1;
                # newsol=0;
                # while (time()-t0<TL) && l<=length(neibour) && newsol <=1 #floor(dt1.N["supplier"]/3) &&
                # neibour2 = First_FBcheck(neibour)
                for l=1:length(neibour)
                    st = PR_FBcheck( prmodel, neibour[l])
                    # , neibour2[l][1+len[1]:sum(len[i] for i=1:2)],
                    #     neibour2[l][1+sum(len[i] for i=1:2):sum(len[i] for i=1:3)],neibour2[l][1+sum(len[i] for i=1:3):sum(len[i] for i=1:4)] )
                    # @show st
                    if st==true
                        sol = value.(all_variables(prmodel)); ndp = getobjval(sol)
                        push!(candSI,sol)
                        if sol ∉ candX && dominated(ndp,candY)==false
                            push!(candX, sol); push!(candY, ndp);
                            newsol+=1;
                            # println("new sol");
                        end
                    end
                #     l+=1;
                end
            end
            if candSI == []
                push!(IGPair,[I,G]);
                # @show iter+=1
                @goto NewIter
            else
                SI = nextSI(candSI,SI)
                if SI∉candX
                    push!(exploredSI,SI);
                end
                # @show iter+=1
            end
        end
        push!(IGPair,[I,G]);
        # iter+=1
    end
    return candX,candY,newsol,IGPair
end
PRtime = @CPUelapsed px,py,pn,pairs = PR(fx,fy,len,round(Int,FPtime))
prx,pry = Postpro(px,py)
########################## Saving the output file ###########################
otable = zeros( length(pry),2)
for i=1:length(pry)
    for j=1:2
        otable[i,j] = pry[i][j]
    end
end
name = ARGS[1][end-7:end];
CSV.write("/home/k2g00/k2g3475/scnd/vopt/lpY/"*name*"lpy.log", DataFrame(otable, :auto), append=false, header=false,delim=' ')
println("time $name: ", LPtime+FPtime+PRtime,": #sol: ", length(pry))


###############################
function Grouping(LB)
    fmin = minimum([LB[i][1] for i=1:length(LB)]),minimum([LB[i][2] for i=1:length(LB)])
    fmax = maximum([LB[i][1] for i=1:length(LB)]),maximum([LB[i][2] for i=1:length(LB)])
    steps = [round.(Int,abs(fmax[k]-fmin[k])/9) for k=1:2]
    cubes = Dict();
    for iter=1:length(LB)
        loca = [round.(Int,((LB[iter][k]-fmin[k])/steps[k])+1) for k=1:2]
        if !haskey(cubes,loca)
            cubes[loca] = [iter]
        else
            push!(cubes[loca], iter)
        end
    end
    groups = collect(values(cubes)); groupkeys = collect(keys(cubes))
    return cubes,groupkeys
end

function GFP(candX,candY,len,TL)
    X = []; PF =[]; Tabu = [];  newsol = 0;
    Y = []; U1 = []; U2= []; U3 = []; t0 = time();
    cubes,gkeys = Grouping(candY)
    glist = copy(gkeys)
    while time() - t0 < TL && glist != []
        g = sample(1:length(glist))
        k = sample(cubes[glist[g]])
        x_t = candX[k]
        yt = x_t[1:len[1]]; yr = copy(yt)
        # u1t = x_t[1+len[1]:len[1]+len[2]];
        # u2t = x_t[1+sum(len[i] for i=1:2):sum(len[i] for i=1:3)]
        # u3t = x_t[1+sum(len[i] for i=1:3):sum(len[i] for i=1:4)]
        SearchDone = false; iter=0; Max_iter = 5#length(findall(i-> 0<i<1,yt))
        while iter < Max_iter && SearchDone == false && time() - t0 < TL
            yid = findall(p->p>0.5,yt);
            for j=1:len[1]
                if j in yid
                    yr[j]=1
                else
                    yr[j]=0
                end
            end
            if FP_FBcheck(fbmodel,yr) == true #,u1r,u2r,u3r)
                sol = value.(all_variables(fbmodel)); ndp = getobjval(sol)
                if sol ∉ X  && dominated(ndp,PF)==false
                    push!(X,sol); push!(PF,ndp)
                    push!(Y,yr);
                    newsol+=1; SearchDone = true
                    deleteat!(glist, g);
                    println("rounding worked")

                end
            else
                if yr ∈ Tabu
                    yr = flipoper(Y,yt,yr);
                    if yr == []
                        SearchDone = true;
                        deleteat!(cubes[glist[g]], cubes[glist[g]] .== k);
                        println("flip failed")
                    else
                        if FP_FBcheck(fbmodel,yr) == true #,u1r,u2r,u3r)
                            sol = value.(all_variables(fbmodel)); ndp = getobjval(sol)
                            if sol ∉ X && dominated(ndp,PF)==false
                                push!(X,sol); push!(PF,ndp)
                                push!(Y,yr);
                                newsol+=1; SearchDone = true;
                                deleteat!(glist, g);
                                println("flip worked")

                            end
                        end
                    end
                end
                if time()-t0 >= TL
                    break
                end
                if SearchDone == false
                    push!(Tabu,yr)
                    yt = fbsearch(yr)
                    if yt==0  #when there's no new feasible lp sol
                        deleteat!(cubes[glist[g]], cubes[glist[g]] .== k);
                        println("no solution")
                        SearchDone = true
                    end
                end
            end
			iter+=1
        end
    end
    return X,PF,newsol,glist
end
FPtime = @CPUelapsed lx,ly,ln,glist = GFP(lp.X_E,lp.Y_N,len,round(Int,LPtime*3))
gx,gy = Postpro(lx,ly)
PRtime = @CPUelapsed px,py,pn,pairs = PR(lx,ly,len,round(Int,FPtime))
prx,pry = Postpro(px,py)


# len = [length(fbm[:yj]),length(fbm[:yk]),length(fbm[:uij]),length(fbm[:ujk]),length(fbm[:ukl]),length(fbm[:xij]),length(fbm[:xjk]),length(fbm[:xkl]),length(fbm[:h])]
# function FP2(candX,len,TL)
#     X = []; PF =Dict(); Tabu = []; t0 = time(); newsol = 0; candlist = copy(candX)
#     Y = []; U1 = []; U2= []; U3 = [];
#     while candlist != [] && time() - t0 < TL
#         @label NewStart
#         k = rand(1:length(candlist))
#         x_t = candlist[k]
#         yjt = x_t[1:len[1]]; ykt = x_t[1+len[1]:sum(len[i] for i=1:2)];
#         jset,kset = JKset(yjt,ykt)
#         if jset == [] ||kset == []
#             deleteat!(candlist,k)
#             println("empty set")
#             @goto NewStart
#         else
#             jid = sample(jset); kid = sample(kset);
#
#         # yt = x_t[1:len[1]];
#             u1t = x_t[1+sum(len[i] for i=1:2):sum(len[i] for i=1:3)]
#             u2t = x_t[1+sum(len[i] for i=1:3):sum(len[i] for i=1:4)]
#             u3t = x_t[1+sum(len[i] for i=1:4):sum(len[i] for i=1:5)]
#             SearchDone = false; iter=0;
#             Max_iter =  5#length(findall(i-> 0<i<1,yt))
#             while iter < Max_iter && SearchDone == false #&& time() - t0 < TL
#                 # yid = findall(p->p>0.3,yt);
#                 for j=1:len[1]
#                     if j in jid
#                         yjt[j]=1
#                     else
#                         yjt[j]=0
#                     end
#                 end
#                 for k=1:len[2]
#                     if k in kid
#                         ykt[k]=1
#                     else
#                         ykt[k]=0
#                     end
#                 end
#                 yr = [yjt;ykt]
#                 u1r = round.(Int, u1t); u2r = round.(Int, u2t); u3r = round.(Int, u3t);
#                 if FP_FBcheck(fbmodel,yr,u1r,u2r,u3r) == true #,u1r,u2r,u3r)
#                     sol = value.(all_variables(fbmodel)); ndp = getobjval(sol)
#                     if sol ∉ X  && dominated(ndp,collect(values(PF)))==false
#                         push!(X,sol);PF[k]=ndp
#                         push!(Y,yr); push!(U1,u1r); push!(U2,u2r); push!(U3,u3r);
#                         newsol+=1; SearchDone = true
#                         deleteat!(candlist, k);
#                         println("rounding worked")
#                     end
#                 else
#                     if yr ∈ Tabu
#                         yr = flipoper(Y,yt,yr);
#                          # u1r = flipoper(U1,u1t,u1r); u2r = flipoper(U2,u2t,u2r); u3r = flipoper(U3,u3t,u3r)
#                         # if any(i->i==[], [yr,u1r,u2r,u3r])
#                         if yr == []
#                             SearchDone = true;
#                         else
#                             if FP_FBcheck(fbmodel,yr,u1r,u2r,u3r) == true #,u1r,u2r,u3r)
#                                 sol = value.(all_variables(fbmodel)); ndp = getobjval(sol)
#                                 if sol ∉ X && dominated(ndp,collect(values(PF)))==false
#                                     push!(X,sol); PF[k]=ndp #push!(PF,ndp)
#                                     push!(Y,yr); push!(U1,u1r); push!(U2,u2r); push!(U3,u3r);
#                                     newsol+=1; SearchDone = true;
#                                     deleteat!(candlist, k);
#                                     println("flip worked")
#                                 end
#                             end
#                         end
#                     end
#                     # if time()-t0 >= TL
#                     #     break
#                     # end
#                     if SearchDone == false
#                         push!(Tabu,[yr;u1r;u2r;u3r])
#                         yt,u1t,u2t,u3t = fbsearch(yr,u1r,u2r,u3r)
#                         if any(i->i==0, [yt,u1t,u2t,u3t])  #when there's no new feasible lp sol
#                             deleteat!(candlist, k);
#                             # println("no solution")
#                             SearchDone = true
#                         end
#                     end
#                 end
#     			iter+=1
#             end
#         end
#     end
#     return X,PF,newsol,candlist
# end
# FPtime = @CPUelapsed lx,ly,ln,candlist = FP2(lp.X_E,len,50)
