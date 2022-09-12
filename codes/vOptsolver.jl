using CPUTime,DataFrames,DelimitedFiles,JuMP,LinearAlgebra,CPLEX,vOptGeneric,SparseArrays,StatsBase
# using MathOptInterface,LazySets
#########################  1dim model  ############################
struct Data1dim
    file::String; N::Dict{}; d::Array{}; c::Array{}; a::Array{}; e::Array{}; gij::SparseVector{}; gjk::SparseVector{}; gkl::SparseVector{};
    vij::Array{}; vjk::Array{}; vkl::Array{}; Vij::SparseVector{}; Vjk::SparseVector{}; Mij::Array{}; Mjk::Array{}; Mkl::Array{};
    b::Array{}; q::Array{}; rij::Array{}; rjk::Array{}; rkl::Array{}; upl::Int; udc::Int; bigM::Int
    function Data1dim(file)
        dt1 = readdlm(file);
        notafile = readdlm("/home/ak121396/Desktop/instances/SCND/Notations.txt", '=');
        # notafile = readdlm("F:/scnd/Notations.txt", '=');
        # notafile = readdlm("/home/k2g00/k2g3475/scnd/Notations.txt", '=');
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
file = "/home/k2g00/k2g3475/scnd/instances/test04S4"
# @show file = ARGS[1]
# file = "F:scnd/Test4S4"
file = "/home/ak121396/Desktop/instances/SCND/test04S4"
dt1 = Data1dim(file);
function SCND1dim()
    scnd1 = vModel(CPLEX.Optimizer); set_silent(scnd1)
    MOI.set(scnd1, MOI.NumberOfThreads(), 1)
    @variable(scnd1, y[1:(dt1.N["plant"]+dt1.N["distribution"])*2], Bin)
    @variable(scnd1, uij[1:sum(dt1.Mij)], Bin);
    @variable(scnd1, ujk[1:sum(dt1.Mjk)], Bin);
    @variable(scnd1, ukl[1:sum(dt1.Mkl)], Bin);
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
            [j=2:dt1.N["plant"]], sum(uij[sum(dt1.Mij[1,1:j-1])+1:sum(dt1.Mij[1,1:j-1])+dt1.Mij[1,j]]) <= 1
            [i=2:dt1.N["supplier"],j=2:dt1.N["plant"]],  sum(uij[sum(dt1.Mij[1:i-1,:])+sum(dt1.Mij[i,1:j-1])+1:sum(dt1.Mij[1:i-1,:])+sum(dt1.Mij[i,1:j-1])+dt1.Mij[i,j]])<= 1
            sum(ujk[1:dt1.Mjk[1,1]]) <= 1
            [k=2:dt1.N["distribution"]], sum(ujk[sum(dt1.Mjk[1,1:k-1])+1:sum(dt1.Mjk[1,1:k-1])+dt1.Mjk[1,k]]) <= 1
            [j=2:dt1.N["plant"],k=2:dt1.N["distribution"]],  sum(ujk[sum(dt1.Mjk[1:j-1,:])+sum(dt1.Mjk[j,1:j-1])+1:sum(dt1.Mjk[1:j-1,:])+sum(dt1.Mjk[j,1:j-1])+dt1.Mjk[j,k]]) <= 1
            sum(ukl[1:dt1.Mkl[1,1]]) <= 1
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
scnd = SCND1dim()
TL = 1200
@CPUtime vSolve(scnd,TL, method=:dicho, verbose=false)
vd = getvOptData(scnd);
vd.Y_N
vd.X_E
@CPUtime vSolve(scnd, method=:epsilon, step=10^5.0, verbose=true);
@CPUtime vSolve(scnd, method=:epsilon, step=10^7.0, verbose=true);
@CPUtime vSolve(scnd, method=:lex, verbose=true);

getY_N(scnd)
ndpoints = getY_N(scnd)

#########################  Feasibility Pump+   ###########################
function FeasibilityModel()
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
            sum(repeat(dt1.b[1,:], outer=sum(dt1.Mij[1,:])).*xij[1:sum(dt1.Mij[1,:])*5]) +
            sum(sum(repeat(dt1.b[i,:], outer=sum(dt1.Mij[i,:])).*xij[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"]) +
            sum(dt1.q.*h) + sum(dt1.rij.*xij)+sum(dt1.rjk.*xjk)+sum(dt1.rkl.*xkl)
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
            [j=2:dt1.N["plant"]], sum(uij[sum(dt1.Mij[1,1:j-1])+1:sum(dt1.Mij[1,1:j-1])+dt1.Mij[1,j]]) <= 1
            [i=2:dt1.N["supplier"],j=2:dt1.N["plant"]],  sum(uij[sum(dt1.Mij[1:i-1,:])+sum(dt1.Mij[i,1:j-1])+1:sum(dt1.Mij[1:i-1,:])+sum(dt1.Mij[i,1:j-1])+dt1.Mij[i,j]])<= 1
            sum(ujk[1:dt1.Mjk[1,1]]) <= 1
            [k=2:dt1.N["distribution"]], sum(ujk[sum(dt1.Mjk[1,1:k-1])+1:sum(dt1.Mjk[1,1:k-1])+dt1.Mjk[1,k]]) <= 1
            [j=2:dt1.N["plant"],k=2:dt1.N["distribution"]],  sum(ujk[sum(dt1.Mjk[1:j-1,:])+sum(dt1.Mjk[j,1:j-1])+1:sum(dt1.Mjk[1:j-1,:])+sum(dt1.Mjk[j,1:j-1])+dt1.Mjk[j,k]]) <= 1
            sum(ukl[1:dt1.Mkl[1,1]]) <= 1
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
function LPrelaxation()
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
    # @objective(lp, Min, sum(dt1.c.*y) +sum(repeat(dt1.a[1,:], outer=sum(dt1.Mij[1,:])).*xij[1:sum(dt1.Mij[1,:])*5])+
    #         sum(sum(repeat(dt1.a[i,:], outer=sum(dt1.Mij[i,:])).*xij[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"])+
    #         sum(dt1.e.*h) + sum(dt1.gij[i]*uij[i] for i in findnz(dt1.gij)[1]) + sum(dt1.gjk[i]*ujk[i] for i in findnz(dt1.gjk)[1]) + sum(dt1.gkl[i].*ukl[i] for i in findnz(dt1.gkl)[1])+
    #         sum(dt1.vij.*xij)+sum(dt1.vjk.*xjk)+sum(dt1.vkl.*xkl) +
            # sum(repeat(dt1.b[1,:], outer=sum(dt1.Mij[1,:])).*xij[1:sum(dt1.Mij[1,:])*5]) +
            # sum(sum(repeat(dt1.b[i,:], outer=sum(dt1.Mij[i,:])).*xij[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"]) +
            # sum(dt1.q.*h) + sum(dt1.rij.*xij)+sum(dt1.rjk.*xjk)+sum(dt1.rkl.*xkl)
    # );
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
            [j=2:dt1.N["plant"]], sum(uij[sum(dt1.Mij[1,1:j-1])+1:sum(dt1.Mij[1,1:j-1])+dt1.Mij[1,j]]) <= 1
            [i=2:dt1.N["supplier"],j=2:dt1.N["plant"]],  sum(uij[sum(dt1.Mij[1:i-1,:])+sum(dt1.Mij[i,1:j-1])+1:sum(dt1.Mij[1:i-1,:])+sum(dt1.Mij[i,1:j-1])+dt1.Mij[i,j]])<= 1
            sum(ujk[1:dt1.Mjk[1,1]]) <= 1
            [k=2:dt1.N["distribution"]], sum(ujk[sum(dt1.Mjk[1,1:k-1])+1:sum(dt1.Mjk[1,1:k-1])+dt1.Mjk[1,k]]) <= 1
            [j=2:dt1.N["plant"],k=2:dt1.N["distribution"]],  sum(ujk[sum(dt1.Mjk[1:j-1,:])+sum(dt1.Mjk[j,1:j-1])+1:sum(dt1.Mjk[1:j-1,:])+sum(dt1.Mjk[j,1:j-1])+dt1.Mjk[j,k]]) <= 1
            sum(ukl[1:dt1.Mkl[1,1]]) <= 1
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
fbmodel = FeasibilityModel() #obj not attached yet
dist = LPrelaxation()
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
function FBcheck(model,yr,u1r,u2r,u3r)
    JuMP.fix.(model[:y],yr; force=true)
    JuMP.fix.(model[:uij],u1r; force=true)
    JuMP.fix.(model[:ujk],u2r; force=true)
    JuMP.fix.(model[:ukl],u3r; force=true)
    optimize!(model)
    if termination_status(model) == MOI.OPTIMAL
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
    @objective( dist, Min, sum(dist[:y][i] for i in idy_0) + sum(1-(dist[:y][j]) for j in idy_1) +
        sum(dist[:uij][i] for i in idu1_0) + sum(1-(dist[:uij][j]) for j in idu1_1)+
        sum(dist[:ujk][i] for i in idu2_0) + sum(1-(dist[:ujk][j]) for j in idu2_1)+
        sum(dist[:ukl][i] for i in idu3_0) + sum(1-(dist[:ukl][j]) for j in idu3_1))
    optimize!(dist)
    if termination_status(dist) == MOI.OPTIMAL
        return JuMP.value.(dist[:y]),JuMP.value.(dist[:uij]),JuMP.value.(dist[:ujk]),JuMP.value.(dist[:ukl])
    else
        return 0,0,0,0
    end
end
len = [length(scnd[:y]),length(scnd[:uij]),length(scnd[:ujk]),length(scnd[:ukl]),length(scnd[:xij]),length(scnd[:xjk]),length(scnd[:xkl]),length(scnd[:h])]
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

    obj1 = sum(dt1.c.*y)+sum(repeat(dt1.a[1,:], outer=sum(dt1.Mij[1,:])).*xij[1:sum(dt1.Mij[1,:])*5])+
            sum(sum(repeat(dt1.a[i,:], outer=sum(dt1.Mij[i,:])).*xij[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"])+
            sum(dt1.e.*h) + sum(dt1.gij[i]*uij[i] for i in findnz(dt1.gij)[1]) + sum(dt1.gjk[i]*ujk[i] for i in findnz(dt1.gjk)[1]) + sum(dt1.gkl[i].*ukl[i] for i in findnz(dt1.gkl)[1])+
            sum(dt1.vij.*xij)+sum(dt1.vjk.*xjk)+sum(dt1.vkl.*xkl)
    obj2 = sum(repeat(dt1.b[1,:], outer=sum(dt1.Mij[1,:])).*xij[1:sum(dt1.Mij[1,:])*5]) +
            sum(sum(repeat(dt1.b[i,:], outer=sum(dt1.Mij[i,:])).*xij[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"]) +
            sum(dt1.q.*h) + sum(dt1.rij.*xij)+sum(dt1.rjk.*xjk)+sum(dt1.rkl.*xkl)
    return [obj1,obj2]
end
function FPplus(dvar,len,TL,Y_N) #𝚯
    X = copy(dvar); IGPair=[]; Tabu = []; newsol=0; t0=time(); PF = copy(Y_N);
    Y = []; U1 = []; U2= []; U3 = [];
    while time()-t0 < TL && length(IGPair)<(length(PF)*(length(PF)-1))
        I,G = StatsBase.sample(1:length(X), 2, replace=false)
        x1 = X[I][1:sum(len[i] for i=1:4)]; x2 = X[G][1:sum(len[i] for i=1:4)];
        λ = round(rand(Float64, 1)[1]; digits=1)
        x_t = x1*λ + x2*(1-λ);
        # x_t = x1*.5 + x2*.5;
        yt = x_t[1:len[1]]; u1t = x_t[1+len[1]:len[1]+len[2]];
        u2t = x_t[1+sum(len[i] for i=1:2):sum(len[i] for i=1:3)]
        u3t = x_t[1+sum(len[i] for i=1:3):sum(len[i] for i=1:4)]
        SearchDone = false; iter=0;
        Max_iter = 5 # length(findall(i-> 0<i<1,x_t))
        while [I,G]∉IGPair && iter<Max_iter && SearchDone == false
            # x_r = round.(Int,x_t);
            yr = round.(Int, yt); u1r = round.(Int, u1t);
            u2r = round.(Int, u2t); u3r = round.(Int, u3t);

            if FBcheck(fbmodel,yr,u1r,u2r,u3r) == true
                sol = value.(all_variables(fbmodel)); ndp = getobjval(sol)
                if sol ∉ X  && dominated(ndp,PF)==false
                    push!(X,sol); push!(PF,ndp)
                    push!(Y,yr); push!(U1,u1r); push!(U2,u2r); push!(U3,u3r);
                    newsol+=1; SearchDone = true
                    # println("rounding worked")
                end
            else
                if [yr;u1r;u2r;u3r] ∈ Tabu
                    yr = flipoper(Y,yt,yr); u1r = flipoper(U1,u1t,u1r); u2r = flipoper(U2,u2t,u2r); u3r = flipoper(U3,u3t,u3r)
                    if any(i->i==[], [yr,u1r,u2r,u3r])
                        SearchDone = true;
                        # println("flip failed")
                    else
                        if FBcheck(fbmodel,yr,u1r,u2r,u3r) == true
                            sol = value.(all_variables(fbmodel)); ndp = getobjval(sol)
                            if sol ∉ X && dominated(ndp,PF)==false
                                push!(X,sol); push!(PF,ndp)
                                push!(Y,yr); push!(U1,u1r); push!(U2,u2r); push!(U3,u3r);
                                newsol+=1; SearchDone = true;
                                # println("flip worked")
                            end
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
                        # println("no solution")
                        SearchDone = true
                    end
                end
            end
			iter+=1
        end
        push!(IGPair,[I,G])

    end
    return X,PF,newsol,IGPair
end
# optimize!(fbmodel); optimize!(dist);
FPtime = @CPUelapsed fx,fy,fn,usedPairs = FPplus(vd.X_E,len,300,vd.Y_N)
vd.Y_N
usedPairs

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
    neibour = setdiff(neibour,exploredSI)
    # neiobj = neiobj[setdiff(1:end, idx),:]

    return neibour#,neiobj
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
function PR(X,Y,len,TL)
    candX = copy(X); candY = copy(Y); bvar = sum(len[i] for i=1:4)
    IGPair=[]; exploredSI = []; t0=time();
    while time()-t0 < TL && length(IGPair)<(length(candY)*(length(candY)-1))
	    I,G = StatsBase.sample(1:length(candX), 2, replace=false)
        SI = candX[I]; SG = candX[G];
        # SI_r = round.(SI); SG_r = round.(SG)
        dif = findall(i-> SI[1:bvar]!=SG[1:bvar], 1:bvar)
        # println("dif is: ", length(dif))
        Max_iter = 10; iter=0;
        while length(dif)>0 && [I,G]∉IGPair && iter<Max_iter && (time()-t0<TL)
            neibour = createNB(SI[1:bvar],dif,exploredSI)

            if (length(neibour)==0) #(time()-t0 >= TL)
                break
            else
                candSI =[]
                l=1; newsol=0;
                while (time()-t0<TL) && l<=length(neibour) && newsol <=1 #floor(dt1.N["supplier"]/3) &&
                # for l=1:length(neibour)
                    st = FBcheck( fbmodel, neibour[l][1:len[1]],neibour[l][1+len[1]:sum(len[i] for i=1:2)],
                        neibour[l][1+sum(len[i] for i=1:2):sum(len[i] for i=1:3)],neibour[l][1+sum(len[i] for i=1:3):sum(len[i] for i=1:4)] )
                    if st==true
                        sol = value.(all_variables(fbmodel)); ndp = getobjval(sol)
                        push!(candSI,sol)
                        if sol ∉ candX && dominated(ndp,candY)==false
                            push!(candX, sol); push!(candY, ndp);
                            newsol+=1;
                            # println("new sol");
                        end
                    end
                    l+=1;
                end
            end
            SI = nextSI(candSI,SI)
            if SI∉candX
                push!(exploredSI,SI);
            end
            @show iter+=1
        end
        push!(IGPair,[I,G])
    end
    return candX,candY,IGPair
end
PRtime = @CPUelapsed prx,pry,prn = PR(fpx,fpy,len,1300)

px,py,pairs = PR(vd.X_E,vd.Y_N,len,1400)


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
fpx,fpy = Postpro(fx,fy)
px = prx; py = pry
prx,pry = Postpro(px,py)
setdiff(vd.Y_N,fpy)
setdiff(fpy,vd.Y_N)
setdiff(py,pry)

# function Postpro(P)
#     #Filter fractional solutions from LB
#     # initdv = candX[1:end-newsol]
#     # frac = findall(j-> trunc.(initdv[j])!= initdv[j], 1:length(initdv))
#     # dv2 = initdv[setdiff(1:end,frac)];
#     # P = union(dv2,candX[end-newsol+1:end])
#     Pobj =[getobjval(P[i]) for i=1:length(P)]
#     #Filter dominated solutions
#
#
#     copysol = Dict(); copyobj = Dict();
#     for i=1:length(Pobj)
#         copysol[i] = P[i]
#         copyobj[i] = Pobj[i]
#     end
#     for i=1:length(Pobj)-1
#         for j=i+1:length(Pobj)
#             if all(Pobj[i] .>= Pobj[j]) == true #dominated by PF[j]
#                 copyobj[i]=nothing; copysol[i]=nothing; break
#             elseif all(Pobj[j] .>= Pobj[i]) == true
#                 copyobj[j]=nothing; copysol[j]=nothing;
#             end
#         end
#     end
#
#     finalsol = filter!(a->a!=nothing, collect(values(copysol)))
#     finalobj = filter!(a->a!=nothing, collect(values(copyobj)))
#
#     return finalsol,finalobj
# end
# finalx,finaly = Postpro(fx)
1
######### check if the point on the convex combination (line) #########
ontheLine = []
for i=1:length(vd.Y_N)-2
    j=i+1; k=i+2
    Line = LineSegment(vd.Y_N[i],vd.Y_N[k])
    if vd.Y_N[j] ∈ Line
        push!(ontheLine, j)
    end
end
# ontheLine



########################## Multi-dimensional model
using CPUTime,DataFrames,DelimitedFiles,JuMP,LinearAlgebra,CPLEX,vOptGeneric,SparseArrays,StatsBase
struct Data
    filepath::String; N::Dict{}; d::Array{}; c::Array{};  Mij::Array{}; Mjk::Array{}; Mkl::Array{};
    vij::Array{}; vjk::Array{}; vkl::Array{}; rij::Array{}; rjk::Array{}; rkl::Array{}; #gij::Array{}; gjk::Array{}; gkl::Array{};
    Vij::Array{}; Vjk::Array{}; b::Array{}; upl::Int; udc::Int; bigM::Int # e::Array{};q::Array{};
    function Data(filepath)
        dt = readdlm(filepath);
        notafile = readdlm("/home/ak121396/Desktop/instances/SCND/Notations.txt", '=');
        # notafile = readdlm("F:/scnd/Notations.txt", '=');
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
        # q = append!(N["vep"],N["ved"]);
        upl = N["upperpants"]; udc = N["upperdistribution"]
        bigM = sum(sum(N["demand"]))
        new(filepath,N,d,c,Mij,Mjk,Mkl,vij,vjk,vkl,rij,rjk,rkl,Vij,Vjk,b,upl,udc,bigM); #cap,Mij,Mjk,Mkl,gij,gjk,gkl,
    end
end

file = "/home/k2g00/k2g3475/scnd/instances/test04S4"
dt = Data(file);

function voptmodel()
    ##########################  Mathematical model  #########################
    scnd = vModel(CPLEX.Optimizer); set_silent(scnd);
    # MOI.set(scnd, MOI.NumberOfThreads(), 1);
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

    @addobjective(scnd, Min, sum(dt.c[j][t]*y[j,t] for j=1:1:dt.N["plant"]+dt.N["distribution"] for t=1:2) + exg +
        sum(dt.N["vcs"][i][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) +
        sum(dt.vij[i][j][m][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
        sum(dt.vjk[j][k][m][p]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
        sum(dt.vkl[k][l][m][p]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5)+
        sum(dt.N["vcp"][j][2*(p-1)+t]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
        sum(dt.N["vcd"][k][2*(p-1)+t]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2));
    #2nd obj
    @addobjective(scnd, Min, sum(dt.b[i,p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
        sum(dt.N["vep"][j][2*(p-1)+t]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
        sum(dt.N["ved"][k][2*(p-1)+t]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2)+
        sum(dt.rij[i][j][m]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
        sum(dt.rjk[j][k][m]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
        sum(dt.rkl[k][l][m]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5));
    ######### constraint 3 #############
    @constraint(scnd,[j=1:dt.N["plant"],p=1:5], sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j]) == sum(xjk[j,k,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]));
    @constraint(scnd,[k=1:dt.N["distribution"],p=1:5], sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k]) == sum(xkl[k,l,m,p] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l]));
    ########### constraint 4-6 #############
    @constraint(scnd,[j=1:dt.N["plant"],p=1:5], sum(h[j,p,:]) == sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j]));
    @constraint(scnd,[k=1:dt.N["distribution"],p=1:5], sum(h[k+dt.N["plant"],p,:] ) == sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k]));
    @constraint(scnd, [l=1:dt.N["customer"],p=1:5], sum(xkl[k,l,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mkl[k,l]) >= dt.d[l][p]);
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
    ########### constraint 13-14 #############
    @constraint(scnd,sum(y[j,t] for j=1:dt.N["plant"] for t=1:2) <= dt.upl);
    @constraint(scnd,sum(y[j,t] for j=dt.N["plant"]+1:dt.N["distribution"]+dt.N["plant"] for t=1:2) <= dt.udc);
    return scnd
end
scndm = voptmodel()
# using vOptGeneric
@CPUtime vSolve(scndm, TL, method=:dicho, verbose=false)
res = getvOptData(scndm);
res.Y_N

vd = getob
############################   wrong model    ###############################

function saveSol(fname::String, Y_N)#, elapsedTime::Float64)
    way = fname*"Y_N"
    open(way, "w") do f
        # write(f, "$elapsedTime \n")
        write(f, "$(length(Y_N)) \n")
        for i in 1:length(Y_N)
            write(f, "$(Y_N[i][1]) $(Y_N[i][2]) \n")
        end
    end
end
path = "/home/k2g00/k2g3475/model/vopt/"*file[36:end]*"Y_N"
saveSol(path, ndpoints)#, elapsedTime)
# ---- Displaying the results (X_E and Y_N)
# for n = 1:length(Y_N)
#     X = value.(x, n)
#     print(findall(elt -> elt ≈ 1, X))
#     println("| z = ",Y_N[n])
# end

# vSolve(model, method=:epsilon, step=1000000000.0, verbose=true);
# ndpoints = getY_N(model)

using vOptGeneric,DelimitedFiles,JuMP,JLD2,LinearAlgebra,CPLEX,MathProgBase,MathOptInterface
const MPB = MathProgBase;

mutable struct CallModel
    lpfile::String; m::Int; n::Int; C::Array{}; B::Array{}; RHS::Dict{}; signs::Array{}; vub::Array{}
    function CallModel(lpfile::String)
        lpmodel=buildlp([-1,0],[2 1],'<',1.5, CplexSolver(CPX_PARAM_SCRIND=0))
        # lpmodel = CPLEX.CplexMathProgModel();
        MPB.loadproblem!(lpmodel,lpfile)
        Bmtx = MPB.getconstrmatrix(lpmodel);
        B = Bmtx[3:end,:]; C = Bmtx[1:2,:]
        m,n=size(B)
        vub = MPB.getvarUB(lpmodel)
        lb = MPB.getconstrLB(lpmodel)[3:end]
        ub = MPB.getconstrUB(lpmodel)[3:end]
        RHS = Dict()
        for i=1:m
            if ub[i]==Inf
                RHS[i] = lb[i]
            else
                RHS[i] = ub[i]
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
struct Valu
    x::String
    y::String
    dvar::Array{}
    LB::Array{}
    LBmtx::Array{}
    function Valu(x, y)
        JLD2.@load x dv
        dv0 = Array(dv)
        # dv0 = readdlm(x)
        dv1 = round.(dv0; digits = 4)
        objs = round.(readdlm(y); digits = 4)
        ind = findall(i -> 0 in objs[i, :], 1:size(objs)[1])
        dv2 = dv1[setdiff(1:end, ind), :]
        LBmtx = objs[setdiff(1:end, ind), 2:end]
        dvar = [Vector(dv2[i, :]) for i = 1:size(LBmtx)[1]]
        LB = [Vector(LBmtx[i, :]) for i = 1:size(LBmtx)[1]]
        new(x, y, dvar, LB, LBmtx)
    end
end
file = ARGS[1];

mt = CallModel("/home/ak121396/Desktop/relise/newlp.lp");
# pr = Valu(
#     "/home/ak121396/Desktop/relise/test01S2_X.jld2",
#     "/home/ak121396/Desktop/relise/test01S2_img_p.sol",
# );
#

function compute(file)

    mt = CallModel(file)

    bvar = findall(i->i==1,mt.vub);
    model = vModel( CPLEX.Optimizer );
    MOI.set(mip, MOI.RawParameter("CPX_PARAM_SCRIND"), false);
    @variable(model, x[1:mt.n] >=0)

    @addobjective( model, Min, dot(x, mt.C[1,:]) )
    @addobjective( model, Min, dot(x, mt.C[2,:]) )

    for k=1:mt.m
        if mt.signs[k] == "l"
            @constraint(model, dot(mt.B[k,:],x) >= mt.RHS[k])
        elseif mt.signs[k] == "u"
            @constraint(model, dot(mt.B[k,:],x) <= mt.RHS[k])
        else
            @constraint(model, dot(mt.B[k,:],x) == mt.RHS[k])
        end
    end
    start = time()
    vSolve(model, method=:dicho, verbose=false)
    # vSolve(model, method=:lex, verbose=false)


    Y_N = getY_N(model)
    elapsedTime = time() - start
    return Y_N, elapsedTime
end


# Y_N, elapsedTime = compute(file)
file = "/home/ak121396/Desktop/instances/SCND/test01S2"
# mt = CallModel(file)
# pr = Valu(
#     "/home/k2g00/k2g3475/model/vopt/X/test01S2_X.jld2",
#     "/home/k2g00/k2g3475/model/vopt/Y/test01S2_img_p.sol",
# );
#
# collector = []
# for j=1:length(pr.dvar)
#     @show j
#     mip = vModel(CPLEX.Optimizer);
#     MOI.set(mip, MOI.RawParameter("CPX_PARAM_SCRIND"), false);
#     bvar = findall(i -> i == 1, mt.vub);
#     @variable(mip, x[1:mt.n] >= 0);
#     for i = 1:mt.n
#         if i in bvar
#             set_binary(x[i])
#         else
#             JuMP.fix(x[i], pr.dvar[j][i]; force = true)
#         end
#     end
#     for k = 1:mt.m
#         if mt.signs[k] == "l"
#             @constraint(mip, dot(mt.B[k, :], x) >= mt.RHS[k])
#         elseif mt.signs[k] == "u"
#             @constraint(mip, dot(mt.B[k, :], x) <= mt.RHS[k])
#         else
#             @constraint(mip, dot(mt.B[k, :], x) == mt.RHS[k])
#         end
#     end
#     @addobjective( mip, Min, dot(x, mt.C[1,:]) );
#     @addobjective( mip, Min, dot(x, mt.C[2,:]) );
#
#     cput = @CPUelapsed vSolve(mip, method=:dicho, verbose=false)
#     if termination_status(mip) == MOI.OPTIMAL
#         println("solution found!")
#         Y_N = getY_N(mip)
#         push!(collector,Y_N)
#     end
# end
# println(collector)

############FPBH
# using JuMP, FPBH, GLPKMathProgInterface, CPLEX, MathProgBase,CPUTime ,Modof
# const MPB = MathProgBase;
# bvar = find(i->i==1,mt.vub); rvar = find(i->i!=1,mt.vub);

# fpbh_model = ModoModel()
# @variable(fpbh_model, yu[i in bvar], Bin);
# @variable(fpbh_model, xh[i in rvar] >=0 );
# @variable(fpbh_model, x[1:mt.n] );
# objective!(fpbh_model, 1, :Min, dot(mt.C[1,:],x));
# objective!(fpbh_model, 2, :Min, dot(mt.C[2,:],x));
#
# for i=1:mt.n
#     if i in bvar
#         @constraint(fpbh_model, x[i]==yu[i]  )
#     else
#         @constraint(fpbh_model, x[i]==xh[i] )
#     end
# end
# for k=1:mt.m
#     if mt.signs[k] == "l"
#         @constraint(fpbh_model, dot(mt.B[k,:],x) >= mt.RHS[k])
#     elseif mt.signs[k] == "u"
#         @constraint(fpbh_model, dot(mt.B[k,:],x) <= mt.RHS[k])
#     else
#         @constraint(fpbh_model, dot(mt.B[k,:],x) == mt.RHS[k])
#     end
# end
#
# runtime = @CPUelapsed sol = fpbh(fpbh_model, solution_polishing=false, timelimit=60.0)
#
