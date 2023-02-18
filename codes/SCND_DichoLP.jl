using CPUTime,DataFrames,DelimitedFiles,JuMP,LinearAlgebra,CPLEX,vOptGeneric,SparseArrays,StatsBase,CSV,JLD2,LazySets,PolygonInbounds
#########################  1dim model  ############################
struct Data1dim
    file::String; N::Dict{}; d::Array{}; c::Array{}; a::Array{}; e::Array{}; gij::SparseVector{}; gjk::SparseVector{}; gkl::SparseVector{};
    vij::Array{}; vjk::Array{}; vkl::Array{}; Vij::SparseVector{}; Vjk::SparseVector{}; Mij::Array{}; Mjk::Array{}; Mkl::Array{};
    b::Array{}; q::Array{}; rij::Array{}; rjk::Array{}; rkl::Array{}; upl::Int; udc::Int; bigM::Int
    function Data1dim(file)
        dt1 = readdlm(file);
        notafile = readdlm("/home/ak121396/Desktop/instances/scnd/Notations.txt", '=');
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
# @show file = ARGS[1]
file = "/home/ak121396/Desktop/instances/scnd/test01S2"
dt1 = Data1dim(file);

# fname = ARGS[1][end-7:end]
# testnum = parse(Int,name[end-3:end-2])
# TL = dt1.N["supplier"]*10*testnum

tnum = 1
if tnum ==  1
    TL = 500
elseif tnum == 2
    TL = 600
elseif tnum == 3
    TL = 800
elseif tnum == 4
    TL = 1200
elseif tnum ==5
    TL = 1400
elseif tnum ==6
    TL = 2000
elseif tnum ==7
    TL = 3200
elseif tnum ==8
    TL = 3600
elseif tnum ==9
    TL = 4200
elseif tnum ==10
    TL = 5500
elseif tnum ==11
    TL = 7500
elseif tnum ==12
    TL = 9000
elseif tnum ==13
    TL = 11000
elseif tnum ==14
    TL = 12000
else
    TL = 15000
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
scndlp = SCND_LP()
#COMPILE
vSolve(scndlp, 2, method=:dicho, verbose=false) 
LPtime = @CPUelapsed vSolve(scndlp, 30, method=:dicho, verbose=false)
println("LPtime: ", LPtime)
lp = getvOptData(scndlp);#getY_N(scndlp)
w1 = round(Int,mean([lp.Y_N[i][1]/lp.Y_N[i][2] for i=1:length(lp.Y_N)]))
len = [length(scndlp[:y]),length(scndlp[:uij]),length(scndlp[:ujk]),length(scndlp[:ukl]),length(scndlp[:xij]),length(scndlp[:xjk]),length(scndlp[:xkl]),length(scndlp[:h])]

function lexobj1()
    lex = Model(CPLEX.Optimizer)
    # optimizer_with_attributes(
            # ,"CPX_PARAM_EPGAP" => 1e-8);
    set_silent(lex)
    MOI.set(lex, MOI.NumberOfThreads(), 1);
    #########################  IP  ########################################
    @variable(lex, y1[1:(dt1.N["plant"]+dt1.N["distribution"])*2], Bin)
    @variable(lex, uij1[1:sum(dt1.Mij)], Bin);
    @variable(lex, ujk1[1:sum(dt1.Mjk)], Bin);
    @variable(lex, ukl1[1:sum(dt1.Mkl)], Bin);

    @variable( lex, 0<= xij1[1:sum(dt1.Mij)*5] );
    @variable( lex, 0<= xjk1[1:sum(dt1.Mjk)*5] );
    @variable( lex, 0<= xkl1[1:sum(dt1.Mkl)*5] );
    @variable( lex, 0<= h1[1:(dt1.N["plant"]+dt1.N["distribution"])*5*2] );
    obj1 = @expression(lex, sum(dt1.c.*y1) +sum(repeat(dt1.a[1,:], outer=sum(dt1.Mij[1,:])).*xij1[1:sum(dt1.Mij[1,:])*5])+
            sum(sum(repeat(dt1.a[i,:], outer=sum(dt1.Mij[i,:])).*xij1[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"])+
            sum(dt1.e.*h1) + sum(dt1.gij[i]*uij1[i] for i in findnz(dt1.gij)[1]) + sum(dt1.gjk[i]*ujk1[i] for i in findnz(dt1.gjk)[1]) + sum(dt1.gkl[i].*ukl1[i] for i in findnz(dt1.gkl)[1])+
            sum(dt1.vij.*xij1)+sum(dt1.vjk.*xjk1)+sum(dt1.vkl.*xkl1))
    obj2 = @expression(lex, sum(repeat(dt1.b[1,:], outer=sum(dt1.Mij[1,:])).*xij1[1:sum(dt1.Mij[1,:])*5]) +
            sum(sum(repeat(dt1.b[i,:], outer=sum(dt1.Mij[i,:])).*xij1[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"]) +
            sum(dt1.q.*h1) + sum(dt1.rij.*xij1)+sum(dt1.rjk.*xjk1)+sum(dt1.rkl.*xkl1))
    @objective(lex, Min, obj1+obj2)
    @constraint(lex, [p=1:5], sum(xij1[5*(m-1)+p] for m=1:dt1.Mij[1,1])+sum(xij1[5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,1]) == sum(xjk1[5*(m-1)+p] for m=1:sum(dt1.Mjk[1,:])) );
    @constraint(lex, [j=2:dt1.N["plant"],p=1:5], sum(xij1[5*sum(dt1.Mij[1,1:j-1])+5*(m-1)+p] for m=1:dt1.Mij[1,j])+sum(xij1[5*sum(dt1.Mij[i,1:j-1])+5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,j]) == sum(xjk1[sum(dt1.Mjk[1:j-1,:])*5 + 5*(m-1)+p] for m=1:sum(dt1.Mjk[j,:])) );
    @constraint(lex, [p=1:5], sum(xjk1[5*(m-1)+p] for m=1:dt1.Mjk[1,1])+sum(xjk1[5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,1]) == sum(xkl1[5*(m-1)+p] for m=1:sum(dt1.Mkl[1,:])) );
    @constraint(lex, [k=2:dt1.N["distribution"],p=1:5],sum(xjk1[5*sum(dt1.Mjk[1,1:k-1])+5*(m-1)+p] for m=1:dt1.Mjk[1,k])+sum(xjk1[5*sum(dt1.Mjk[j,1:k-1])+5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,k]) == sum(xkl1[sum(dt1.Mkl[1:k-1,:])*5 + 5*(m-1)+p] for m=1:sum(dt1.Mkl[k,:])) );
    ########### constraint 4-6 #############
    @constraint(lex, [p=1:5],sum(h1[2*(p-1)+t] for t=1:2) == sum(xij1[5*(m-1)+p] for m=1:dt1.Mij[1,1])+sum(xij1[5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,1]));
    @constraint(lex, [j=2:dt1.N["plant"],p=1:5], sum(h1[5*2*(j-1)+2*(p-1)+t] for t=1:2) == sum(xij1[5*sum(dt1.Mij[1,1:j-1])+5*(m-1)+p] for m=1:dt1.Mij[1,j])+sum(xij1[5*sum(dt1.Mij[i,1:j-1])+5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,j]) );
    @constraint(lex, [p=1:5], sum(h1[5*2*dt1.N["plant"]+2*(p-1)+t] for t=1:2) == sum(xjk1[5*(m-1)+p] for m=1:dt1.Mjk[1,1])+sum(xjk1[5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,1]) );
    @constraint(lex, [k=2:dt1.N["distribution"],p=1:5], sum(h1[5*2*dt1.N["plant"]+5*2*(k-1)+2*(p-1)+t] for t=1:2) == sum(xjk1[5*sum(dt1.Mjk[1,1:k-1])+5*(m-1)+p] for m=1:dt1.Mjk[1,k])+sum(xjk1[5*sum(dt1.Mjk[j,1:k-1])+5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,k]));
    @constraint(lex, [p=1:5], sum(xkl1[5*(m-1)+p] for m=1:dt1.Mkl[1,1]) +sum(xkl1[5*(m-1)+p+(5*sum(dt1.Mkl[1:k-1,:]))] for k=2:dt1.N["distribution"] for m=1:dt1.Mkl[k,1]) >= dt1.d[1,p]);
    @constraint(lex, [l=2:dt1.N["customer"], p=1:5], sum(xkl1[sum(dt1.Mkl[1,1:l-1])*5 + 5*(m-1)+p] for m=1:dt1.Mkl[1,l])+ sum(xkl1[5*sum(dt1.Mkl[1:k-1,:])+5*sum(dt1.Mkl[k,1:l-1])+5*(m-1)+p] for k=2:dt1.N["distribution"] for m=1:dt1.Mkl[k,l]) >= dt1.d[l,p]);
    ########### constraint 7 #############
    @constraint(lex, sum(xij1[1:5*sum(dt1.Mij[1,:])]) <= dt1.N["cas"][1]);
    @constraint(lex, [i=2:dt1.N["supplier"]],  sum(xij1[5*sum(dt1.Mij[1:i-1,:])+1:5*sum(dt1.Mij[1:i,:])]) <= dt1.N["cas"][i]);
    ########### constraint 8 #############
    # @constraint(lex,[j=1:dt1.N["plant"]+dt1.N["distribution"], t=1:2], sum(h1[5*2*(j-1)+((t-1)*5)+1:5*2*(j-1)+((t-1)*5)+5]) <= [dt1.N["cap"];dt1.N["cad"]][j]*y1[(dt1.N["plant"]+dt1.N["distribution"])*(t-1)+j])
    @constraint(lex,[j=1:dt1.N["plant"]+dt1.N["distribution"], t=1:2], sum(h1[5*2*(j-1)+((p-1)*2)+t] for p=1:5) <= [dt1.N["cap"];dt1.N["cad"]][j]*y1[2*(j-1)+t]);
    ########### constraint 9 #############
    @constraint(lex,[j=1:dt1.N["plant"]+dt1.N["distribution"]], sum(y1[2*(j-1)+1:2*(j-1)+2]) <= 1);
    ########### constraint 10 #############
    @constraints(lex,begin sum(uij1[1:dt1.Mij[1,1]]) <= 1
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
    @constraints(lex, begin
        [i=1:sum(dt1.Mij)], sum(xij1[5*(i-1)+1:5*i]) <= dt1.bigM*uij1[i]
        [i=1:sum(dt1.Mjk)], sum(xjk1[5*(i-1)+1:5*i]) <= dt1.bigM*ujk1[i]
        [i=1:sum(dt1.Mkl)], sum(xkl1[5*(i-1)+1:5*i]) <= dt1.bigM*ukl1[i]
    end);
    ########### constraint 12 #############
    @constraints(lex, begin
        [i in findnz(dt1.Vij)[1]], sum(xij1[5*(i-1)+1:5*i]) >= dt1.Vij[i]*uij1[i]
        [i in findnz(dt1.Vjk)[1]], sum(xjk1[5*(i-1)+1:5*i]) >= dt1.Vjk[i]*ujk1[i]
        # [i in findnz(dt1.Vkl)[1]], sum(xkl1[5*(i-1)+1:5*i]) >= dt1.Vkl[i]*ukl1[i]
    end);
    ########### constraint 13-14 #############
    @constraint(lex, sum(y1[1:dt1.N["plant"]*2]) <= dt1.upl);
    @constraint(lex, sum(y1[dt1.N["plant"]*2+1:end]) <= dt1.udc);
    return lex
end
nwline = lexobj1()
set_optimizer_attribute(nwline, "CPXPARAM_TimeLimit", 3); optimize!(nwline);
set_optimizer_attribute(nwline, "CPXPARAM_TimeLimit", LPtime*max(tnum*8,50));
l1time = @CPUelapsed optimize!(nwline); lex1X = [value.(all_variables(nwline))]; lex1Y = [getobjval(value.(all_variables(nwline)))]

function lexobj2()
    lex = Model(CPLEX.Optimizer)
    # optimizer_with_attributes(
            # ,"CPX_PARAM_EPGAP" => 1e-8);
    set_silent(lex)
    MOI.set(lex, MOI.NumberOfThreads(), 1);
    
    #########################  IP  ########################################
    @variable(lex, y1[1:(dt1.N["plant"]+dt1.N["distribution"])*2], Bin)
    @variable(lex, uij1[1:sum(dt1.Mij)], Bin);
    @variable(lex, ujk1[1:sum(dt1.Mjk)], Bin);
    @variable(lex, ukl1[1:sum(dt1.Mkl)], Bin);

    @variable( lex, 0<= xij1[1:sum(dt1.Mij)*5] );
    @variable( lex, 0<= xjk1[1:sum(dt1.Mjk)*5] );
    @variable( lex, 0<= xkl1[1:sum(dt1.Mkl)*5] );
    @variable( lex, 0<= h1[1:(dt1.N["plant"]+dt1.N["distribution"])*5*2] );
    obj2 = @expression(lex, sum(repeat(dt1.b[1,:], outer=sum(dt1.Mij[1,:])).*xij1[1:sum(dt1.Mij[1,:])*5]) +
            sum(sum(repeat(dt1.b[i,:], outer=sum(dt1.Mij[i,:])).*xij1[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"]) +
            sum(dt1.q.*h1) + sum(dt1.rij.*xij1)+sum(dt1.rjk.*xjk1)+sum(dt1.rkl.*xkl1))
        obj1 = @expression(lex, sum(dt1.c.*y1) +sum(repeat(dt1.a[1,:], outer=sum(dt1.Mij[1,:])).*xij1[1:sum(dt1.Mij[1,:])*5])+
    sum(sum(repeat(dt1.a[i,:], outer=sum(dt1.Mij[i,:])).*xij1[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"])+
    sum(dt1.e.*h1) + sum(dt1.gij[i]*uij1[i] for i in findnz(dt1.gij)[1]) + sum(dt1.gjk[i]*ujk1[i] for i in findnz(dt1.gjk)[1]) + sum(dt1.gkl[i].*ukl1[i] for i in findnz(dt1.gkl)[1])+
    sum(dt1.vij.*xij1)+sum(dt1.vjk.*xjk1)+sum(dt1.vkl.*xkl1))
    # @objective(lex, Min, obj2)
    @objective(lex, Min, obj2+(obj1/(w1*w1)))
    @constraint(lex, [p=1:5], sum(xij1[5*(m-1)+p] for m=1:dt1.Mij[1,1])+sum(xij1[5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,1]) == sum(xjk1[5*(m-1)+p] for m=1:sum(dt1.Mjk[1,:])) );
    @constraint(lex, [j=2:dt1.N["plant"],p=1:5], sum(xij1[5*sum(dt1.Mij[1,1:j-1])+5*(m-1)+p] for m=1:dt1.Mij[1,j])+sum(xij1[5*sum(dt1.Mij[i,1:j-1])+5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,j]) == sum(xjk1[sum(dt1.Mjk[1:j-1,:])*5 + 5*(m-1)+p] for m=1:sum(dt1.Mjk[j,:])) );
    @constraint(lex, [p=1:5], sum(xjk1[5*(m-1)+p] for m=1:dt1.Mjk[1,1])+sum(xjk1[5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,1]) == sum(xkl1[5*(m-1)+p] for m=1:sum(dt1.Mkl[1,:])) );
    @constraint(lex, [k=2:dt1.N["distribution"],p=1:5],sum(xjk1[5*sum(dt1.Mjk[1,1:k-1])+5*(m-1)+p] for m=1:dt1.Mjk[1,k])+sum(xjk1[5*sum(dt1.Mjk[j,1:k-1])+5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,k]) == sum(xkl1[sum(dt1.Mkl[1:k-1,:])*5 + 5*(m-1)+p] for m=1:sum(dt1.Mkl[k,:])) );
    ########### constraint 4-6 #############
    @constraint(lex, [p=1:5],sum(h1[2*(p-1)+t] for t=1:2) == sum(xij1[5*(m-1)+p] for m=1:dt1.Mij[1,1])+sum(xij1[5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,1]));
    @constraint(lex, [j=2:dt1.N["plant"],p=1:5], sum(h1[5*2*(j-1)+2*(p-1)+t] for t=1:2) == sum(xij1[5*sum(dt1.Mij[1,1:j-1])+5*(m-1)+p] for m=1:dt1.Mij[1,j])+sum(xij1[5*sum(dt1.Mij[i,1:j-1])+5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,j]) );
    @constraint(lex, [p=1:5], sum(h1[5*2*dt1.N["plant"]+2*(p-1)+t] for t=1:2) == sum(xjk1[5*(m-1)+p] for m=1:dt1.Mjk[1,1])+sum(xjk1[5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,1]) );
    @constraint(lex, [k=2:dt1.N["distribution"],p=1:5], sum(h1[5*2*dt1.N["plant"]+5*2*(k-1)+2*(p-1)+t] for t=1:2) == sum(xjk1[5*sum(dt1.Mjk[1,1:k-1])+5*(m-1)+p] for m=1:dt1.Mjk[1,k])+sum(xjk1[5*sum(dt1.Mjk[j,1:k-1])+5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,k]));
    @constraint(lex, [p=1:5], sum(xkl1[5*(m-1)+p] for m=1:dt1.Mkl[1,1]) +sum(xkl1[5*(m-1)+p+(5*sum(dt1.Mkl[1:k-1,:]))] for k=2:dt1.N["distribution"] for m=1:dt1.Mkl[k,1]) >= dt1.d[1,p]);
    @constraint(lex, [l=2:dt1.N["customer"], p=1:5], sum(xkl1[sum(dt1.Mkl[1,1:l-1])*5 + 5*(m-1)+p] for m=1:dt1.Mkl[1,l])+ sum(xkl1[5*sum(dt1.Mkl[1:k-1,:])+5*sum(dt1.Mkl[k,1:l-1])+5*(m-1)+p] for k=2:dt1.N["distribution"] for m=1:dt1.Mkl[k,l]) >= dt1.d[l,p]);
    ########### constraint 7 #############
    @constraint(lex, sum(xij1[1:5*sum(dt1.Mij[1,:])]) <= dt1.N["cas"][1]);
    @constraint(lex, [i=2:dt1.N["supplier"]],  sum(xij1[5*sum(dt1.Mij[1:i-1,:])+1:5*sum(dt1.Mij[1:i,:])]) <= dt1.N["cas"][i]);
    ########### constraint 8 #############
    # @constraint(lex,[j=1:dt1.N["plant"]+dt1.N["distribution"], t=1:2], sum(h1[5*2*(j-1)+((t-1)*5)+1:5*2*(j-1)+((t-1)*5)+5]) <= [dt1.N["cap"];dt1.N["cad"]][j]*y1[(dt1.N["plant"]+dt1.N["distribution"])*(t-1)+j])
    @constraint(lex,[j=1:dt1.N["plant"]+dt1.N["distribution"], t=1:2], sum(h1[5*2*(j-1)+((p-1)*2)+t] for p=1:5) <= [dt1.N["cap"];dt1.N["cad"]][j]*y1[2*(j-1)+t]);
    ########### constraint 9 #############
    @constraint(lex,[j=1:dt1.N["plant"]+dt1.N["distribution"]], sum(y1[2*(j-1)+1:2*(j-1)+2]) <= 1);
    ########### constraint 10 #############
    @constraints(lex,begin sum(uij1[1:dt1.Mij[1,1]]) <= 1
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
    @constraints(lex, begin
        [i=1:sum(dt1.Mij)], sum(xij1[5*(i-1)+1:5*i]) <= dt1.bigM*uij1[i]
        [i=1:sum(dt1.Mjk)], sum(xjk1[5*(i-1)+1:5*i]) <= dt1.bigM*ujk1[i]
        [i=1:sum(dt1.Mkl)], sum(xkl1[5*(i-1)+1:5*i]) <= dt1.bigM*ukl1[i]
    end);
    ########### constraint 12 #############
    @constraints(lex, begin
        [i in findnz(dt1.Vij)[1]], sum(xij1[5*(i-1)+1:5*i]) >= dt1.Vij[i]*uij1[i]
        [i in findnz(dt1.Vjk)[1]], sum(xjk1[5*(i-1)+1:5*i]) >= dt1.Vjk[i]*ujk1[i]
        # [i in findnz(dt1.Vkl)[1]], sum(xkl1[5*(i-1)+1:5*i]) >= dt1.Vkl[i]*ukl1[i]
    end);
    ########### constraint 13-14 #############
    @constraint(lex, sum(y1[1:dt1.N["plant"]*2]) <= dt1.upl);
    @constraint(lex, sum(y1[dt1.N["plant"]*2+1:end]) <= dt1.udc);
    return lex
end
ndline = lexobj2()
set_optimizer_attribute(ndline, "CPXPARAM_TimeLimit", 3); optimize!(ndline)
set_optimizer_attribute(ndline, "CPXPARAM_TimeLimit", LPtime*max(tnum*8,50));
@show l2time = @CPUelapsed optimize!(ndline);
lex2X = [value.(all_variables(ndline))]; lex2Y = [getobjval(value.(all_variables(ndline)))]

# jldsave("/home/k2g00/k2g3475/scnd/vopt/lpY/Lex/"*fname*"leX.jld2"; leX1=lex1X, leY1=lex1Y,leX2=lex2X, leY2=lex2Y)
# load("/home/k2g00/k2g3475/scnd/vopt/lpY/Lex/"*fname*"leX.jld2")
len = [length(scndlp[:y]),length(scndlp[:uij]),length(scndlp[:ujk]),length(scndlp[:ukl]),length(scndlp[:xij]),length(scndlp[:xjk]),length(scndlp[:xkl]),length(scndlp[:h])]

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
            R = StatsBase.sample(1:M,Num, replace=false)
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

function FP_FBcheck(model,yr,iter)

    if sum(yr[1:dt1.N["plant"]*2])>dt1.upl || sum(yr[1+dt1.N["plant"]*2:len[1]])>dt1.udc 
        return false
    else
        for j=1:dt1.N["plant"]+dt1.N["distribution"]
            if sum(yr[2*(j-1)+1:2*(j-1)+2]) <= 1
                continue
            else
                return false
            end
        end
    end
        
    if isodd(iter)==true
        @objective(model, Min, sum(dt1.c.*model[:y]) +sum(repeat(dt1.a[1,:], outer=sum(dt1.Mij[1,:])).*model[:xij][1:sum(dt1.Mij[1,:])*5])+
        sum(sum(repeat(dt1.a[i,:], outer=sum(dt1.Mij[i,:])).*model[:xij][sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"])+
        sum(dt1.e.*model[:h]) + sum(dt1.gij[i]*model[:uij][i] for i in findnz(dt1.gij)[1]) + sum(dt1.gjk[i]*model[:ujk][i] for i in findnz(dt1.gjk)[1]) + sum(dt1.gkl[i].*model[:ukl][i] for i in findnz(dt1.gkl)[1])+
        sum(dt1.vij.*model[:xij])+sum(dt1.vjk.*model[:xjk])+sum(dt1.vkl.*model[:xkl])) 
    else
        @objective(model, Min,sum(repeat(dt1.b[1,:], outer=sum(dt1.Mij[1,:])).*model[:xij][1:sum(dt1.Mij[1,:])*5]) +
        sum(sum(repeat(dt1.b[i,:], outer=sum(dt1.Mij[i,:])).*model[:xij][sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"]) +
        sum(dt1.q.*model[:h]) + sum(dt1.rij.*model[:xij])+sum(dt1.rjk.*model[:xjk])+sum(dt1.rkl.*model[:xkl]))
    end

    JuMP.fix.(model[:y],yr; force=true)
    optimize!(model)
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
function fbsearch2(yr,u1r,u2r,u3r) #solveLP
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
function NDfilter(P,Pobj)
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
function SortingSol(P,Pobj)
    copysol = Dict(); copyobj = Dict();
    for i=1:length(Pobj)
        copysol[i] = P[i]
        copyobj[i] = Pobj[i]
    end
    for i=1:length(Pobj)-1
        for j=i+1:length(Pobj)
            if all(Pobj[i] .>= Pobj[j]) == true #dominated by PF[j]
                copyobj[i]=0; copysol[i]=0; break
            elseif all(Pobj[j] .>= Pobj[i]) == true
                copyobj[j]=0; copysol[j]=0;
            end
        end
    end
    sortedsol = filter!(a->a!=0, collect(values(copysol)))
    sortedobj = filter!(a->a!=0, collect(values(copyobj)))
    df = DataFrame(X=sortedsol,Y=sortedobj);
    sort!(df,[:Y])
    return df
end
function FP_Model()
    model = Model(CPLEX.Optimizer); set_silent(model)
    MOI.set(model, MOI.NumberOfThreads(), 1)
    @variable(model, y[1:(dt1.N["plant"]+dt1.N["distribution"])*2], Bin)
    @variable(model, uij[1:sum(dt1.Mij)], Bin);
    @variable(model, ujk[1:sum(dt1.Mjk)], Bin);
    @variable(model, ukl[1:sum(dt1.Mkl)], Bin);
    @variable( model, 0<= xij[1:sum(dt1.Mij)*5] );
    @variable( model, 0<= xjk[1:sum(dt1.Mjk)*5] );
    @variable( model, 0<= xkl[1:sum(dt1.Mkl)*5] );
    @variable( model, 0<= h[1:(dt1.N["plant"]+dt1.N["distribution"])*5*2] );

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
    # @constraint(model,[j=1:dt1.N["plant"]+dt1.N["distribution"]], sum(y[2*(j-1)+1:2*(j-1)+2]) <= 1);
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
    # @constraint(model, sum(y[1:dt1.N["plant"]*2]) <= dt1.upl);
    # @constraint(model, sum(y[dt1.N["plant"]*2+1:end]) <= dt1.udc);
    return model
end
function LP_Model()
    lp = Model(CPLEX.Optimizer); set_silent(lp)
    MOI.set(lp, MOI.NumberOfThreads(), 1)
    @variable(lp, 0 <= y[1:(dt1.N["plant"]+dt1.N["distribution"])*2] <= 1)
    @variable(lp, 0 <= uij[1:sum(dt1.Mij)] <= 1);
    @variable(lp, 0 <= ujk[1:sum(dt1.Mjk)] <= 1);
    @variable(lp, 0 <= ukl[1:sum(dt1.Mkl)] <= 1);
    @variable( lp, 0 <= xij[1:sum(dt1.Mij)*5] );
    @variable( lp, 0 <= xjk[1:sum(dt1.Mjk)*5] );
    @variable( lp, 0 <= xkl[1:sum(dt1.Mkl)*5] );
    @variable( lp, 0 <= h[1:(dt1.N["plant"]+dt1.N["distribution"])*5*2] );
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
function PR_Model()
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
    # @constraint(model,[j=1:dt1.N["plant"]+dt1.N["distribution"]], sum(y[2*(j-1)+1:2*(j-1)+2]) <= 1);
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
    # @constraint(model, sum(y[1:dt1.N["plant"]*2]) <= dt1.upl);
    # @constraint(model, sum(y[dt1.N["plant"]*2+1:end]) <= dt1.udc);
    return model
end
function createAllNB(SI,dif,exploredSI)
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
    end
    # idx = findall(i-> neibour[i] in exploredSI, 1:length(neibour))
    return neibour
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
    # idx = findall(i-> neibour[i] in exploredSI, 1:length(neibour))
    # neibour1 = setdiff(neibour,exploredSI)
    neibour2 = StatsBase.StatsBase.sample(neibour, round(Int,length(neibour)*0.1) , replace=false)
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

function PR_preFBcheck(neibour)
    #constraint 13 & 14
    neibour2 = filter!(y->(sum(y[1:dt1.N["plant"]*2])<=dt1.upl && sum(y[1:dt1.N["plant"]*2])>0) || (sum(y[1+dt1.N["plant"]*2:len[1]])<=dt1.udc && sum(y[1+dt1.N["plant"]*2:len[1]])>0), neibour)
    #constraint 9
    for j=1:dt1.N["plant"]+dt1.N["distribution"]
        filter!(y-> sum(y[2*(j-1)+1:2*(j-1)+2]) <= 1, neibour2)
    end
    return neibour2
end

function PR_FBcheck_w(model,weight,yr)
    @objective(model, Min, sum(dt1.c.*model[:y]) +sum(repeat(dt1.a[1,:], outer=sum(dt1.Mij[1,:])).*model[:xij][1:sum(dt1.Mij[1,:])*5])+
    sum(sum(repeat(dt1.a[i,:], outer=sum(dt1.Mij[i,:])).*model[:xij][sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"])+
    sum(dt1.e.*model[:h]) + sum(dt1.gij[i]*model[:uij][i] for i in findnz(dt1.gij)[1]) + sum(dt1.gjk[i]*model[:ujk][i] for i in findnz(dt1.gjk)[1]) + sum(dt1.gkl[i].*model[:ukl][i] for i in findnz(dt1.gkl)[1])+
    sum(dt1.vij.*model[:xij])+sum(dt1.vjk.*model[:xjk])+sum(dt1.vkl.*model[:xkl]) +
    weight*(sum(repeat(dt1.b[1,:], outer=sum(dt1.Mij[1,:])).*model[:xij][1:sum(dt1.Mij[1,:])*5]) +
    sum(sum(repeat(dt1.b[i,:], outer=sum(dt1.Mij[i,:])).*model[:xij][sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"]) +
    sum(dt1.q.*model[:h]) + sum(dt1.rij.*model[:xij])+sum(dt1.rjk.*model[:xjk])+sum(dt1.rkl.*model[:xkl])));

    JuMP.fix.(model[:y],yr; force=true)
    optimize!(model)
    if has_values(model) == true #termination_status(model) == MOI.OPTIMAL
        return true
    else
        return false
    end
end
function PR_FBcheck(model,iter,yr)
    if isodd(iter)==true
        @objective(model, Min, sum(dt1.c.*model[:y]) +sum(repeat(dt1.a[1,:], outer=sum(dt1.Mij[1,:])).*model[:xij][1:sum(dt1.Mij[1,:])*5])+
        sum(sum(repeat(dt1.a[i,:], outer=sum(dt1.Mij[i,:])).*model[:xij][sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"])+
        sum(dt1.e.*model[:h]) + sum(dt1.gij[i]*model[:uij][i] for i in findnz(dt1.gij)[1]) + sum(dt1.gjk[i]*model[:ujk][i] for i in findnz(dt1.gjk)[1]) + sum(dt1.gkl[i].*model[:ukl][i] for i in findnz(dt1.gkl)[1])+
        sum(dt1.vij.*model[:xij])+sum(dt1.vjk.*model[:xjk])+sum(dt1.vkl.*model[:xkl])) 
    else
        @objective(model, Min,sum(repeat(dt1.b[1,:], outer=sum(dt1.Mij[1,:])).*model[:xij][1:sum(dt1.Mij[1,:])*5]) +
        sum(sum(repeat(dt1.b[i,:], outer=sum(dt1.Mij[i,:])).*model[:xij][sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"]) +
        sum(dt1.q.*model[:h]) + sum(dt1.rij.*model[:xij])+sum(dt1.rjk.*model[:xjk])+sum(dt1.rkl.*model[:xkl]))
    end
    JuMP.fix.(model[:y],yr; force=true)
    optimize!(model)
    if has_values(model) == true #termination_status(model) == MOI.OPTIMAL
        return true
    else
        return false
    end
end

fpmodel = FP_Model(); optimize!(fpmodel);
set_optimizer_attribute(fpmodel, "CPXPARAM_TimeLimit", 10)
dist = LP_Model(); optimize!(dist);
set_optimizer_attribute(dist, "CPXPARAM_TimeLimit", 10)
prmodel = PR_Model(); optimize!(prmodel);
set_optimizer_attribute(prmodel, "CPXPARAM_TimeLimit", 10);

function FP(candX,len,TL)
    X = []; PF =[]; Tabu = []; newsol = 0; Y = [];
    candlist = StatsBase.sample(1:length(candX), length(candX), replace=false)
    t0 = time(); 
    for k in candlist  #&& time() - t0 < TL 
        x_t = candX[k]
        yt = x_t[1:len[1]]; yr = copy(yt)
        SearchDone = false; iter=0;
        Max_iter = 5#length(findall(i-> 0<i<1,yt))
        while iter < Max_iter && SearchDone == false && time() - t0 < TL
            yid = findall(p->p>0.2,yt);
            for j=1:len[1]
                if j in yid
                    yr[j]=1
                else
                    yr[j]=0
                end
            end
            if FP_FBcheck(fpmodel,yr,iter) == true #,u1r,u2r,u3r)
                sol = value.(all_variables(fpmodel)); ndp = getobjval(sol)
                if sol ∉ X  && dominated(ndp,collect(values(PF)))==false
                    push!(X,sol); push!(PF,ndp) #PF[k] = ndp
                    push!(Y,yr);
                    newsol+=1; SearchDone = true; 
                    # deleteat!(candlist, k);
                end
            else
                if yr ∈ Tabu
                    yr = flipoper(Y,yt,yr); # u1r = flipoper(U1,u1t,u1r); u2r = flipoper(U2,u2t,u2r); u3r = flipoper(U3,u3t,u3r)
                    if yr == [] # if any(i->i==[], [yr,u1r,u2r,u3r])
                        SearchDone = true; 
                        # deleteat!(candlist, k);
                       
                    else
                        if FP_FBcheck(fpmodel,yr,iter) == true #,u1r,u2r,u3r)
                            sol = value.(all_variables(fpmodel)); ndp = getobjval(sol)
                            if sol ∉ X && dominated(ndp,collect(values(PF)))==false
                                push!(X,sol); push!(Y,yr); push!(PF,ndp) #PF[k] = ndp
                                newsol+=1; SearchDone = true;
                            end
                        end
                    end
                end
                if time()-t0 >= TL
                    break
                end
            
                if SearchDone == false
                    push!(Tabu,yr)
                    yt = fbsearch(yr);
                    # println("new lp sol")
                    if yt==0  #when there's no new feasible lp sol
                        SearchDone = true
                    end
                end
            end
			iter+=1
        end
    end
    return X,PF
end
FP(lp.X_E,len,3)
FPtime = @CPUelapsed f1x,f1y = FP(lp.X_E,len,Inf)
f2x,f2y = NDfilter(f1x,f1y);
println("FPtime: ", FPtime, " FPsol: ", length(f2y))
fpX,fpY = NDfilter([f2x;lex1X;lex2X],[f2y;lex1Y;lex2Y])
dfp = DataFrame(X=fpX,Y=fpY);
sort!(dfp, [:Y])


function FPplus(dvar,Y_N,len,TL) 
    X = copy(dvar); PF = copy(Y_N); Y = []; IGPair=[]; Tabu = []; t0=time();
    # U1 = []; U2= []; U3 = []; newsol=0; 
    while time()-t0 < TL && length(IGPair)<(length(PF)*(length(PF)-1))
        I,G = StatsBase.StatsBase.sample(1:length(X), 2, replace=false)
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

            if FP_FBcheck(fpmodel,yr,iter) == true
                sol = value.(all_variables(fpmodel)); ndp = getobjval(sol)
                if sol ∉ X  && dominated(ndp,PF)==false
                    push!(X,sol); push!(PF,ndp)
                    # push!(Y,yr); push!(U1,u1r); push!(U2,u2r); push!(U3,u3r);
                    # newsol+=1; 
                    SearchDone = true
                    # println("rounding worked")
                end
            else
                if yr ∈ Tabu #[yr;u1r;u2r;u3r] ∈ Tabu
                    yr = flipoper(Y,yt,yr);
                    # u1r = flipoper(U1,u1t,u1r); u2r = flipoper(U2,u2t,u2r); u3r = flipoper(U3,u3t,u3r)
                    if yr == [] # if any(i->i==[], [yr,u1r,u2r,u3r])
                        SearchDone = true;
                        # println("flip failed")
                    else
                        if FP_FBcheck(fpmodel,yr,iter) == true
                            sol = value.(all_variables(fpmodel)); ndp = getobjval(sol)
                            if sol ∉ X && dominated(ndp,PF)==false
                                push!(X,sol); push!(PF,ndp)
                                # newsol+=1; push!(Y,yr); push!(U1,u1r); push!(U2,u2r); push!(U3,u3r);
                                SearchDone = true;
                                # println("flip worked")
                            end
                        end
                    end
                end
                if time()-t0 >= TL
                    break
                end
                if SearchDone == false
                    push!(Tabu,yr) #[yr;u1r;u2r;u3r]
                    yt,u1t,u2t,u3t = fbsearch2(yr,u1r,u2r,u3r)
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
    return X,PF#,IGPair,newsol
end
FPplus(dfp.X,dfp.Y,len,3)
FPPtime = @CPUelapsed fpp1x,fpp1y = FPplus(dfp.X,dfp.Y,len,round(Int,TL/3))
fpp2x,fpp2y = NDfilter(fpp1x,fpp1y);
println("FPPtime: ", FPPtime, " FPPsol: ", length(fpp2y))
dfpp = DataFrame(X=fpp2x,Y=fpp2y);
sort!(dfpp,[:Y])
dX = Dict(i=>dfpp.X[i] for i=1:length(dfpp.Y)); dY = Dict(i=>dfpp.Y[i] for i=1:length(dfpp.Y));

############## Initialising nondominated line segments
struct node
    val::Array; arm::String; #dom::String x::Float64; y::Float64
end
function Newnodes(lsg)
    nodes = []
    for i=1:length(lsg)
        if i == 1
            connect = "R" 
        elseif i == length(lsg)
            connect = "L" 
        else
            connect = "LR" 
        end
        push!(nodes, node(lsg[i], connect)) #, "null"))
    end
    return nodes
end

function SolveLPdicho(sol,ndp)
    JuMP.fix.(LPdicho[:y], sol[1:len[1]]; force = true)
    JuMP.fix.(LPdicho[:uij], sol[1+len[1]:sum(len[i] for i=1:2)]; force = true)
    JuMP.fix.(LPdicho[:ujk], sol[1+sum(len[i] for i=1:2):sum(len[i] for i=1:3)]; force = true)
    JuMP.fix.(LPdicho[:ukl], sol[1+sum(len[i] for i=1:3):sum(len[i] for i=1:4)]; force = true)
    # dtime = @CPUelapsed 
    vSolve(LPdicho, 5, method=:dicho, verbose=false)
    res = getvOptData(LPdicho);
    if res == []
        dichosol = [ndp]
    else
        dichosol = sort!(res.Y_N)
    end
    return dichosol #dtime
end

function SectionQueue(ndset,newsg)
    que = Matrix(undef,length(newsg),length(ndset)) 
    for i=1:length(newsg)
        for j=1:length(ndset)
            if newsg[i].val[1] >= ndset[j].val[1] && newsg[i].val[2] >= ndset[j].val[2] 
                que[i,j] = "u"
            elseif newsg[i].val[1] <= ndset[j].val[1] && newsg[i].val[2] <= ndset[j].val[2] 
                que[i,j] = "d"
            elseif newsg[i].val[1] >= ndset[j].val[1] && newsg[i].val[2] <= ndset[j].val[2] 
                que[i,j] = "r"
            else
                que[i,j] = "l"
            end
        end
    end
    return que
end
function dominance_count(x,P)
    st = 0
    for k=1:length(P)
        if all( x .<= P[k])
            st+=1; 
        else
            continue
        end
    end
    return st
end
function ClassifyND(ndset,dsol)
    Larms = findall(x-> ndset[x].arm == "L", 1:length(ndset))
    start = 0;  ndom = []; dom = [];

    for l in Larms 
        refset = ndset[1+start:l]
        polygon = [refset[i].val for i=1:length(refset)]
        push!(polygon, [99^30, refset[end].val[2]],[99^30, refset[1].val[2]])
        edges = zeros(Int,length(polygon),2)
        for i=1:length(polygon)-1
            edges[i,1] = i; edges[i,2] = i+1
        end
        edges[end,1] = length(polygon); edges[end,2] = 1
        stat = inpoly2(dsol, polygon, edges, atol = 1e-1)
        tmpnd = findall(x-> stat[x,:] == [0,0], 1:length(dsol))
        push!(ndom, tmpnd)
        union!(dom, filter( x-> x ∉ tmpnd, 1:length(dsol)))
        start = l
    end
    
    ndom2 = nothing
    if length(ndom)>1
        ndom2 = ndom[1]
        for i=2:length(ndom)
            ndom2 = intersect(ndom2, ndom[i])
        end
    elseif length(ndom)==1 
        ndom2 = ndom[1]
    end

    return ndom2,dom
end
function UpdateNDset(ddicho,ndicho,domset,ndomset,newsg_og,ndset_og)
    ndset = copy(ndset_og)
    if length(newsg_og) == 1                                   # new solution is a point
        que = SectionQueue(ndset_og,newsg_og)
        nw = newsg_og[1] 
        if "u" ∈ que                                     # new solution is dominated
            return ndset
        elseif all(x->x == "d", que)                     # new solution dominate the whole ndset
            ndset = nw; return ndset
        elseif all(x->x == "l", que)      
            insert!(ndset, 1, nw); return ndset               
        elseif all(x->x == "r", que)      
            insert!(ndset, length(ndset)+1, nw); return ndset
        end
    
        d_id = findall(x->x =="d", que)
        if length(d_id) == 1 # # new solution is dominating one ndpoint 
            case = 0
            if ndset[d_id[1]].arm == "null" # one dominated point without an arm
                replace!(ndset, ndset[d_id[1]] => nw)
            else  # one dominated point with an arm
                if d_id[1] == 1                 # 1st ndset point is dominated
                    lsgx = LineSegment(nw.val,[nw.val[1],ndset[2].val[2]])
                    pjx = LazySets.isdisjoint(LineSegment(ndset[1].val,ndset[2].val), lsgx, true)
                    case = 1
                elseif d_id[1] == length(ndset) # last ndset point is dominated
                    pjy = LazySets.isdisjoint(LineSegment(ndset[end-1].val,ndset[end].val), Line(nw.val,[1,0.]))
                    case = 2
                else                            # midle ndset point is dominated
                    lsgx = LineSegment(nw.val,[ndset[d_id[1]-1].val[1],nw.val[2]])
                    pjx = LazySets.is_intersection_empty(LineSegment(ndset[d_id[1]-1].val,ndset[d_id[1]].val), lsgx, true)
                    pjy = LazySets.isdisjoint(LineSegment(ndset[end-1].val,ndset[end].val), Line(nw.val,[1,0.]))
                    case = 3
                end
                if case == 1
                    replace!(ndset, ndset[d_id[1]] => nw); insert!(ndset, d_id[1]+1, node(pjx[2], "R"))
                elseif case == 2
                    replace!(ndset, ndset[end] => node(pjy[2], "L")); insert!(ndset, length(ndset), nw)
                else 
                    replace!(ndset, ndset[d_id[1]] => node(pjy[2], "L")); 
                    insert!(ndset, d_id[1]+1, nw)
                    insert!(ndset, d_id[1]+2, node(pjx[2], "R"))
                end
            end
        elseif length(d_id) >1 # new solution is dominating some ndpoints/linesegment
            case = 0  
            if d_id[1] == 1
                if d_id[end].arm == "R" || d_id[end].arm == "LR"
                    case = 1 # case1: only pjx calculation needed 
                else
                    case =0 # case0: replace dominated lsg and return ndset
                end
            else
                if ndset[d_id[1]].arm == "L" || ndset[d_id[1]].arm == "LR"
                    if d_id[end].arm == "R" || d_id[end].arm == "LR"
                        case = 3
                    else
                        case = 2
                    end
                else
                    if d_id[end].arm == "R" || d_id[end].arm == "LR"
                        case = 3
                    else
                        case =0
                    end
                end
            end
            deleteat!(ndset, d_id);
            if case == 0
                insert(ndset, d_id[1], nw)
                return ndset
            elseif case ==1 
                if d_id[1] == 1
                    lsgx = LineSegment(nw.val,[nw.val[1],ndset[2].val[2]])
                    pjx = LazySets.is_intersection_empty(LineSegment(ndset[1].val,ndset[2].val), lsgx, true)
                else
                    lsgx = LineSegment(nw.val,[ndset[d_id[1]-1].val[1],nw.val[2]])
                    pjx = LazySets.is_intersection_empty(LineSegment(ndset[d_id[1]-1].val,ndset[d_id[1]].val), lsgx, true)
                end
                insert!(ndset, d_id[1], nw); insert!(ndset, d_id[1]+1, node(pjx[2], "R"))
            elseif case == 2
                pjy = LazySets.isdisjoint(LineSegment(ndset[d_id[1]-1].val,ndset[d_id[1]].val), Line(nw.val,[1.,0.]))
                insert!(ndset, d_id[1], node(pjy[2], "L")); insert!(ndset, d_id[1]+1, nw)
            else 
                insert!(ndset, d_id[1], node(pjy[2], "L")); 
                insert!(ndset, d_id[1]+1, nw)
                insert!(ndset, d_id[1]+2, node(pjx[2], "R"))
            end
        end
    
        if "r,l"∈ que 
            rl = findall( x->que[x]=="r" && que[x+1] == "l"  ,1:length(que))[1]
            if ndset[rl].arm == "R" # connected linesegment
                
                if abs((ndset[rl].val[2]-ndset[rl+1].val[2])/(ndset[rl].val[1]-ndset[rl+1].val[1])) > abs((ndset[rl].val[2]-nw.val[2])/(ndset[rl].val[1]-nw.val[1])) #new sol is nondominated
                    if ndset[rl].arm != "null"  # divide the lsg into two parts and insert the new sol
                        seg = [ndset[rl].val,ndset[rl+1].val]
                        p1 = LineSegment(id[end](seg[1]),id[end](seg[2]))
                        p2 = LineSegment(id[end](nw.val),id[end](nw.val[1],seg[1][2]))
                        interpt1 = intersect(p1,p2)[2]
                        p3 = LineSegment(id[end](nw.val),id[end](seg[2][1],nw.val[2]))
                        interpt2 = intersect(p1,p3)[2]
    
                        p1 = LineSegment(id[end](seg[1]),id[end](seg[2]))
                        p2 = LineSegment(id[end](nw[u_id-1].val),id[end](nw[u_id].val))
                    
                        #insert interpt1,nw,interpt2
                        insert!( ndset, rl+1, node(interpt1[1], interpt1[2], "L") )
                        insert!( ndset, rl+2, node(nw.val, "null") )
                        insert!( ndset, rl+3, node(interpt2[1], interpt2[2], "R") )
                    else
                        insert!( ndset, rl+1, node(nw.val, "null") ) #insert the new sol btw two points
                    end
                end
            else                # new solution between two points
                if dominated(nw.val,[ndset[rl].val,ndset[rl+1].val]) == false
                    insert!( ndset, rl+1, node(nw.val, "null") )  # add new sol only if it's nondominated
                end
            end
        end

    else    
        newsg = copy(newsg_og)    
        removal = []; addnw = []; ndstart = 1;
        for k =1:length(newsg)-1
            if k ∈ ddicho && k+1 ∈ ddicho
            else
                for i=ndstart:length(ndset)-1
                    @show (k,i)
                    if i ∈ ndomset && i+1 ∈ ndomset  
                        println("all nd nondominated. next: ", ndstart)
                        ndstart = i+1
                    elseif i ∈ domset && i+1 ∈ domset
                        println("all nd domset")
                        push!(removal,ndset[i,i+1]) 
                    else
                        
                        @show (newsg[k], ndset[i])
                        if ndset[i].arm == "R" || ndset[i].arm =="LR"

                            nwline = LazySets.LineSegment(newsg[k].val,newsg[k+1].val)
                            ndline = LazySets.LineSegment(ndset[i].val,ndset[i+1].val)
                            inter = LazySets.isdisjoint(nwline,ndline,true)                    
                            fourpt0 = sort!([newsg[k],newsg[k+1],ndset[i],ndset[i+1]], by = x-> x.val[1]) 
                            fourpt = copy(fourpt0)
                            if isempty(inter[2]) # two line segments are disjoint

                                # 2nd point projected to yaxis
                                if fourpt0[2].val ∈ nwline
                                    l3 = ndline
                                else
                                    l3 = nwline
                                end
                                pj2y = LazySets.isdisjoint(Line(fourpt0[2].val,[0.,1]), l3, true)
                                if isempty(pj2y[2]) # no intersection
                                    dom3 = dominated(fourpt0[3].val,[fourpt0[1].val,fourpt0[2].val])
                                    dom4 = dominated(fourpt0[4].val,[fourpt0[1].val,fourpt0[2].val])

                                    if (dom3,dom4) == (1,1) #nd points all dominated
                                        println("dom 1,1::: shouldn't be...")

                                    elseif (dom3,dom4) == (0,0)  
                                        println("dom = 0,0")
                                        if fourpt0[3].val ∈ nwline
                                            if i == length(ndset)-1 
                                                push!(addnw, newsg[k], newsg[k+1])
                                                @goto NextPair
                                            else
                                                ndstart = i+1
                                            end
                                        else
                                            push!(addnw, newsg[k], newsg[k+1])
                                            @goto NextPair
                                        end
                                    else # (dom3,dom4) == (1,0)  #pt3 dominated
                                        println("dom = 1,0")
                                        lsg2x = LazySets.LineSegment(fourpt0[2].val,[fourpt0[4].val[1],fourpt0[2].val[2]])
                                        pj2x = is_intersection_empty(lsg2x,l3,true)
                                        if fourpt0[3].val ∈ nwline
                                            push!(removal, newsg[k])
                                            push!(addnw, node(pj2x[2], "R"), newsg[k+1])
                                            newsg[k] = node(pj2x[2], "R")
                                            ndstart = i+1
                                        else
                                            push!(removal, ndset[i])
                                            push!(addnw, newsg[k], node(newsg[k+1].val, "L") , node(pj2x[2], "R"))
                                            insert!(ndset, i+1, node(pj2x[2], "R"))
                                            @goto NextPair
                                        end
                                        
                                    end

                                elseif pj2y[2][2] > fourpt0[2].val[2]
                                    println((k,i),"add from pj 2y")
                                    if fourpt0[2].val ∈ nwline
                                        # replace!(ndset, ndset[i+1] => node(ndset[i+1].val, "R"))
                                        push!(ndset, node(pj2y[2], "L"), node(newsg[i].val, "R"))
                                    else
                                        
                                        push!(addnw, newsg[k], node(pj2y[2], "L"))
                                    end
                                elseif pj2y[2][2] < fourpt0[2].val[2]
                                    println(" no add from pj 2y")

                                    if fourpt0[2].val ∈ nwline
                                        push!(removal, newsg[k])
                                    else
                                        push!(addnw, newsg[k])
                                        push!(removal, ndset[i])
                                    end
                                end
                        
                                # 3rd point projected to xaxis
                                if fourpt0[3].val[2] > fourpt0[4].val[2] 
                                    if fourpt0[3].val ∈ nwline
                                        l3 = ndline
                                    else
                                        l3 = nwline
                                    end
                                    lsg3x = LazySets.LineSegment([fourpt0[1].val[1],fourpt0[3].val[2]],[fourpt0[4].val[1],fourpt0[3].val[2]])
                                    pj3x = is_intersection_empty(lsg3x,l3,true)
                                    
                                    if isempty(pj3x[2]) 
                                        ndstart = i+1
                                    elseif pj3x[2][1] < fourpt0[3].val[1] # no intersection
                                        println("no add from proj 3")
                                        if fourpt0[3].val ∈ nwline
                                            push!(removal, newsg[k])
                                        else
                                            push!(addnw, newsg[k+1]); push!(removal, ndset[i+1])
                                        end
                                        ndstart = i+1
                                    elseif pj3x[2][1] > fourpt0[3].val[1]
                                        println("added from proj 3")
                                        if fourpt0[3].val ∈ nwline
                                            push!(addnw, node(newsg[k+1], "L"))
                                            insert!(ndset, i+1, node(pj3x[2], "R"))
                                            @goto NextPair
                                        else
                                            newsg[k] = node(pj3x[2], "R")
                                            push!(addnw, newsg[k+1])
                                            ndstart = i+1
                                        end
                                    end
                                else
                                    if fourpt0[3].val ∈ nwline
                                        push!(addnw,newsg[k+1]); push!(removal, ndset[i+1])
                                    else
                                        push!(removal, newsg[k+1])
                                    end
                                    ndstart = i+1; @goto NextPair
                                end
                            else #two line segments intersect
                                # line segments share the point => calculate with three points
                                if inter[2] ∈ [newsg[k].val,newsg[k+1].val,ndset[i].val,ndset[i+1].val]
                                    println("sharing a point")
                                    threept0 = copy(fourpt0)
                                    interid = findall(x->x.val==inter[2],threept0)[1]
                                    deleteat!(threept0, interid)
                                    threept = copy(threept0)
                                    if interid == 1
                                        if threept0[2].val[2] < threept0[3].val[2]
                                            if threept0[3] ∈ newsg
                                                @goto NextPair
                                            else
                                                push!(addnw, newsg[k+1]); push!(removal, ndset[i+1])
                                                ndstart = i+1
                                            end
                                        else
                                            lsg2x = LazySets.LineSegment([threept0[1].val[1],threept0[2].val[2]],[threept0[3].val[1],threept0[2].val[2]])
                                            pj2x = LazySets.isdisjoint(lsg2x, LineSegment(threept0[1].val,threept0[3].val), true)
                                            if pj2x[2][1] < threept0[2].val[1] # no intersection
                                                if threept0[2] ∈ newsg
                                                    @goto NextPair
                                                else
                                                    push!(removal, ndset[i+1])
                                                    push!(addnw, newsg[k+1])
                                                    ndstart = i+1
                                                end
                                            else
                                                if threept0[2] ∈ newsg
                                                    push!(addnw, node(newsg[k+1], "L"))
                                                    insert!(ndset, i+1, node(pj3x[2], "R"))
                                                    ndstart = i+1
                                                    @goto NextPair    
                                                else
                                                    replace!(newsg, newsg[k] => node(pj3x[2], "R"))
                                                    push!(addnw, newsg[k+1],node(pj2x[2], "R"))
                                                    ndstart = i+1
                                                end
                                            end
                                        end
                                    elseif interid == 3
                                        if threept0[1].val[2] < threept0[2].val[2]
                                            if threept0[1] ∈ newsg
                                                push!(addnw, newsg[k]); push!(removal, ndset[i])
                                            end
                                        else
                                            pj2y = LazySets.isdisjoint(Line(threept0[2].val, [1.,0.]), Line(threept0[1].val,threept0[3].val), true)
                                            push!(addnw, newsg[k], node(pj2y[2], "L"))
                                        end
                                        ndstart = i+1 
                                        @goto NextPair
                                    elseif interid == 2
                                        if threept0[2] ∈ newsg
                                            push!(addnw, newsg[k])
                                        else
                                            push!(addnw, newsg[k+1])
                                        end
                                        ndstart = i+1 
                                        @goto NextPair
                                    end

                                # calculate with four points
                                else 
                                    if fourpt0[2].val ∈ nwline
                                        l3 = ndline
                                    else
                                        l3 = nwline
                                    end
                                    # 2nd point projected to yaxis
                                    pj2y = LazySets.isdisjoint(Line(fourpt0[2].val,[0.,1.]), l3, true)
                                    if pj2y[2][2] > fourpt0[2].val[2]
                                        if fourpt0[2].val ∈ nwline
                                            println(" new pj2y and inter")
                                            # fourpt[2].arm changed from "LR" to "R"
                                            push!(addnw, node(pj2y[2], "L"), node(fourpt[2].val, "R"), node(inter[2], "LR"))
                                        else
                                            push!(addnw, fourpt[2], node(pj2y[2], "L"), node(inter[2], "LR"))
                                        end
                                    else
                                        println(" replace fourpt[2] with inter")
                                        if fourpt0[2].val ∈ nwline
                                            push!(removal, newsg[k])
                                            push!(addnw,  node(inter[2], "LR"))
                                        else
                                            push!(removal, ndset[i])
                                            push!(addnw, newsg[k], node(inter[2], "LR"))
                                        end
                                    end
                                    
                                    #3rd point projected to xaxis
                                    #slope calculation
                                    nwdif = abs.(newsg[k].val-newsg[k+1].val)
                                    slop_newseg = abs(nwdif[2]/nwdif[1])
                                    nddif = abs.(ndset[i].val-ndset[i+1].val)
                                    slop_ndseg = abs(nddif[2]/nddif[1])
                                    if fourpt0[3].val ∈ nwline
                                        l3 = ndline
                                        pjslop = slop_newseg; theother_slop = slop_ndseg
                                    else
                                        l3 = nwline
                                        pjslop = slop_ndseg; theother_slop = slop_newseg
                                    end
                                    
                                    if fourpt0[3].val[2] > fourpt0[4].val[2] && abs(pjslop) > abs(theother_slop)  
                                        lsg3x = LazySets.LineSegment([fourpt0[1].val[1],fourpt0[3].val[2]],[fourpt0[4].val[1],fourpt0[3].val[2]])
                                        pj3x = is_intersection_empty(lsg3x,l3,true)
                                        if fourpt0[3].val ∈ nwline
                                            println(i," new pj3x ")
                                            push!(addnw, newsg[k+1], node(pj3x[2], "R"))
                                            insert!(ndset, i+1, node(pj3x[2], "R"))
                                            ndstart = i+1
                                            @goto NextPair
                                        else
                                            push!(addnw, node(pj3x[2], "R"), newsg[k+1])
                                            newsg[k] = node(pj3x[2], "R")
                                            ndstart = i+1
                                        end
                                    elseif fourpt0[3].val[2] < fourpt0[4].val[2]
                                        if fourpt0[3].val ∈ nwline
                                            push!(removal, ndset[i+1])
                                            push!(addnw, newsg[k+1])
                                        end
                                        ndstart = i+1
                                        @goto NextPair
                                    end
                                end
                            end
                        elseif ndset[i].arm =="L"
                            if i != length(ndset) 
                                ndstart = i+1
                                println("only Left arm => nextnode", ndstart)
                            else #last node
                                return ndset
                            end
                        else # ndset[i] is a "point"
                            if ndset[i].val[1] < newsg[k].val[1]
                                dc = dominance_count(ndset[i].val,[newsg[j].val for j=k:k+1])
                                if dc == 0
                                    insert!(ndset, i+1, newsg[k])
                                    insert!(ndset, i+2, newsg[k+1])
                                elseif dc == 1
                                    lsgx = LineSegment(ndset[i].val,[newsg[k+1].val[1],ndset[i].val[2]])
                                    pjx = LazySets.is_intersection_empty( lsgx, LineSegment(newsg[k].val,newsg[k+1].val), true)
                                    insert(ndset, i+1 , node(pjx[2], "LR")); insert(ndset, i+2 , newsg[k+1])
                                end
                            elseif ndset[i].val[1] > newsg[k+1].val[1] && ndset[i].val[2] < newsg[k+1].val[2]
                                insert!(ndset, i, newsg[k]); insert!(ndset, i+1, newsg[k+1])
                            elseif newsg[k].val[1] < ndset[i].val[1] < newsg[k+1].val[1] 
                                ndslop = (ndset[i].val[2] - newsg[k+1].val[2])/(ndset[i].val[1] - newsg[k+1].val[1])
                                segslop = (newsg[k+1].val[2] - newsg[k].val[2])/(newsg[k].val[1] - newsg[k+1].val[1])
                                if ndslop < 0 && ndslop > segslop
                                    pjy = LazySets.isdisjoint( Line(ndset[i].val,[1.,0.]), LineSegment(newsg[k].val,newsg[k+1].val), true)
                                    lsgx = LineSegment(ndset[i].val,[newsg[k+1].val[1],ndset[i].val[2]])
                                    pjx = LazySets.is_intersection_empty( lsgx, LineSegment(newsg[k].val,newsg[k+1].val), true)
                                    insert!(ndset, i, newsg[k]); insert!(ndset, i+1, node(pjy[2], "L")); 
                                    insert!(ndset, i+3, node(pjx[2], "R")); insert!(ndset, i+4, newsg[k])
                                elseif ndslop > 0
                                    pjy = LazySets.isdisjoint( Line(ndset[i].val,[1.,0.]), LineSegment(newsg[k].val,newsg[k+1].val), true)
                                    insert!(ndset, i, newsg[k]); insert!(ndset, i+1, node(pjy[2], "L")); 
                                end
                            end
                        end
                    end
                end
            end
            @label NextPair
        end
        unique!(addnw); 
        unique!(removal)
        append!(ndset,addnw)
        append!(ndset,newsg_og[ndicho])
        setdiff!(ndset,removal)
        setdiff!(ndset,ndset_og[domset])
    end
    return ndset
end
########### Initialise the first nondominated set with the first FPP solution 
function InitialNDset(dfpp)
    LPdicho = SCND_LP()
    dsol1 = SolveLPdicho(dfpp.X[1],dfpp.Y[1])
    # for i=1:length(dfpp.Y)
    #     dfpp1[i] = round.(dfpp1[i])
    # end
    if length(dsol1) > 1
        ndset = Newnodes(dsol1)
    else
        ndset = node(dsol1[1], "null")
    end
    for l=2:length(dfpp.Y)
        dsol= SolveLPdicho(dfpp.X[l],dfpp.Y[l])
        for i=1:length(dsol)
            dsol[i] = round.(dsol[i])
        end
        ndicho, ddicho = ClassifyND(ndset,dsol)
        ndlist = [ndset[i].val for i=1:length(ndset)]
        newsg = Newnodes(dsol)
        ndomset, domset = ClassifyND(newsg,ndlist)
        ndset = UpdateNDset(ddicho,ndicho,domset,ndomset,newsg,ndset)
    end
    return ndset    
end

ndset0 = InitialNDset(dfpp)

dsol,dtime = SolveLPdicho(dfpp.X[3],dfpp.Y[3])
ndset = Newnodes(sort(dfpp.Y))
set2 = copy(ndset)
newset = UpdateNDset(dsol,ndset,que)

function PR(dX,dY,len,TL)
    X = copy(dX); Y = copy(dY); IGPair=[]; 
    exploredSI = Dict(i=>dX[i] for i=1:length(dX)); exploredval = Dict(i=>dY[i] for i=1:length(dY)); t0=time();    

    println("now randomly choose IG pairs")
    main_iter = 0; main_infeasi = 0;
    pc1 = 0.7
    while time()-t0 < TL*pc1
        nd = SortingSol(X,Y)
        if length(IGPair)>=(length(nd.Y)*(length(nd.Y)-1))
            return nd
        end
        @label NewIter1
        I0,G0 = StatsBase.StatsBase.sample(nd.Y, 2, replace=false);
        I = findall(i->i==I0,Y)[1]; G = findall(i->i==G0,Y)[1]
        # if [I,G] ∈ IGPair || abs(Y[I][1]-Y[G][1]) < 10^4|| abs(Y[I][2]-Y[G][2]) < 10^4 
        #     @goto NewIter1
        # end
        SI = X[I]; SG = X[G];
        nosol = 0; #exploredSI = [] #copy(collect(values(X)));
        while nosol < 5 && all.(SI != SG) && time()-t0 < TL*pc1
            dif = findall(i-> SI[1:len[1]][i]!=SG[1:len[1]][i], 1:len[1])
            # weight = round(Int,mean([nd.Y[i][1]/nd.Y[i][2] for i=1:length(nd.Y)]))
            new_weight = round(Int, (Y[I][1]+Y[G][1])/(Y[I][2]+Y[G][2]))
            neibour = createNB(SI[1:len[1]],dif,exploredSI)
            newsol = 0; 
            if (length(neibour)==0) #(time()-t0 >= TL)
                @goto NewIter1
            else
                candSI =[]
                neibour2 = PR_preFBcheck(neibour)
                # println("removed neighbour: ", length(neibour)-length(neibour2))
                for l=1:length(neibour2)
                    main_iter+=1
                    if neibour2[l] ∈ collect(values(exploredSI))
                        # keyy = findall(k-> exploredSI[k] == neibour2[l], 1:length(exploredSI))[1]; ndp = exploredval[keyy]; 
                        push!(candSI,sol) 
                    else
                        # st = PR_FBcheck( prmodel, main_iter, neibour[l])
                        st = PR_FBcheck_w( prmodel, new_weight, neibour2[l])
                        if st==true
                            sol = value.(all_variables(prmodel)); ndp = getobjval(sol)
                            push!(candSI,sol); 
                            push!(exploredSI,length(exploredSI)+1 => sol); push!(exploredval, length(exploredval)+1 => ndp)
                            if sol ∉ nd.X && dominated(ndp,nd.Y)==false
                                push!(X, length(X)+1 => sol); push!(Y, length(Y)+1 => ndp)                            
                                newsol+=1;
                                # println(main_iter, " new sol ", ndp);
                            end
                        else
                            main_infeasi+=1
                        end
                        if time()-t0 >= TL*pc1
                            println("__breaking point__!")
                            break
                        end 
                    end
                end
            end
            if candSI == []
                push!(IGPair,[I,G]);
                @goto NewIter1
            else
                SI = nextSI(candSI,SI)
                # if SI∉collect(values(exploredSI))
                #     push!(exploredSI,length(exploredSI)+1 => SI); push!(exploredval, length(exploredval)+1 => SIobj)
                # end
            end
            if newsol==0;
                nosol+=1
            end
        end
        push!(IGPair,[I,G]);
    end

    println("____1st Left____")
    left1_iter = 0; left1_infeasi = 0;
    pc2 =0.8
    while time()-t0 < TL*pc2
        nd = SortingSol(X,Y)
        if length(IGPair)>=(length(nd.Y)*(length(nd.Y)-1))
            return nd
        end
        @label NewIter2
        I0,G0 = StatsBase.StatsBase.sample(1:round(Int,length(nd.Y)*0.2), 2, replace=false);
        I = findall(i->i==nd.Y[I0],Y)[1]; G = findall(i->i==nd.Y[G0],Y)[1]
        # if [I,G] ∈ IGPair || abs(Y[I][1]-Y[G][1]) < 10^5|| abs(Y[I][2]-Y[G][2]) < 10^5
        #     @goto NewIter2
        # end
        SI = X[I]; SG = X[G];
        nosol = 0; #exploredSI = []#copy(collect(values(X)));
        while nosol < 5 && all.(SI != SG) && time()-t0 < TL*pc2
            dif = findall(i-> SI[1:len[1]][i]!=SG[1:len[1]][i], 1:len[1])
            new_weight = round(Int, (Y[I][1]+Y[G][1])/(Y[I][2]+Y[G][2]))
            neibour = createNB(SI[1:len[1]],dif,exploredSI)
            newsol = 0;
            if (length(neibour)==0) 
                @goto NewIter2
            else
                candSI =[]
                neibour2 = PR_preFBcheck(neibour)
                # println("removed neighbour: ", length(neibour)-length(neibour2))
                for l=1:length(neibour2)
                    left1_iter+=1
                    if neibour2[l] ∈ collect(values(exploredSI))
                        # keyy = findall(k-> exploredSI[k] == neibour2[l], 1:length(exploredSI))[1]; ndp = exploredval[keyy]; 
                        push!(candSI,sol) 
                    else
                        st = PR_FBcheck_w( prmodel, new_weight, neibour2[l])
                        if st==true
                            sol = value.(all_variables(prmodel)); ndp = getobjval(sol)
                            push!(candSI,sol); 
                            push!(exploredSI,length(exploredSI)+1 => sol); push!(exploredval, length(exploredval)+1 => ndp)
                            if sol ∉ nd.X && dominated(ndp,nd.Y)==false
                                push!(X, length(X)+1 => sol); push!(Y, length(Y)+1 => ndp)                            
                                newsol+=1;
                                # println(main_iter, " new sol ", ndp);
                            end
                        else
                            left1_infeasi+=1
                        end
                        if time()-t0 >= TL*pc2
                            println("__breaking point__!")
                            break
                        end 
                    end
                end
            end  
            if candSI == []
                push!(IGPair,[I,G]); 
                @goto NewIter2
            else
                SI = nextSI(candSI,SI)
                # if SI∉collect(values(X))
                #     push!(exploredSI,SI);
                # end
            end
            if newsol == 0
                nosol+=1
            end
        end
        push!(IGPair,[I,G]);
    end

    println("____2nd Left____")
    left2_iter = 0; left2_infeasi = 0;
    while time()-t0 < TL
        nd = SortingSol(X,Y)
        if length(IGPair)>=(length(nd.Y)*(length(nd.Y)-1))
            return nd
        end
        @label NewIter3
        buffer = round(Int,length(nd.Y)*0.15)
        I0,G0 = StatsBase.StatsBase.sample(buffer:round(Int,length(nd.Y)*0.3), 2, replace=false);
        I = findall(i->i==nd.Y[I0],Y)[1]; G = findall(i->i==nd.Y[G0],Y)[1]
        # if [I,G] ∈ IGPair || abs(Y[I][1]-Y[G][1]) < 10^5|| abs(Y[I][2]-Y[G][2]) < 10^5 
        #     @goto NewIter3
        # end
        SI = X[I]; SG = X[G];
        nosol = 0; #exploredSI = []#copy(collect(values(X)));
        while nosol < 5 && all.(SI != SG)&& time()-t0 < TL
            dif = findall(i-> SI[1:len[1]][i]!=SG[1:len[1]][i], 1:len[1])
            # weight = round(Int,mean([nd.Y[i][1]/nd.Y[i][2] for i=1:length(nd.Y)]))
            new_weight = round(Int, (Y[I][1]+Y[G][1])/(Y[I][2]+Y[G][2]))
            neibour = createNB(SI[1:len[1]],dif,exploredSI)
            newsol = 0;
            if (length(neibour)==0) 
                @goto NewIter3
            else
                candSI =[]
                neibour2 = PR_preFBcheck(neibour)
                # println("removed neighbour: ", length(neibour)-length(neibour2))
                for l=1:length(neibour2)
                    left2_iter+=1
                    if neibour2[l] ∈ collect(values(exploredSI))
                        # keyy = findall(k-> exploredSI[k] == neibour2[l], 1:length(exploredSI))[1]; ndp = exploredval[keyy]; 
                        push!(candSI,sol) 
                    else
                        st = PR_FBcheck_w( prmodel, new_weight, neibour2[l])
                        
                        if st==true
                            sol = value.(all_variables(prmodel)); ndp = getobjval(sol)
                            push!(candSI,sol); 
                            push!(exploredSI,length(exploredSI)+1 => sol); push!(exploredval, length(exploredval)+1 => ndp)
                            if sol ∉ nd.X && dominated(ndp,nd.Y)==false
                                push!(X, length(X)+1 => sol); push!(Y, length(Y)+1 => ndp)                            
                                newsol+=1;
                                # println(main_iter, " new sol ", ndp);
                            end 
                        else
                            left2_infeasi+=1
                        end
                        if time()-t0 >= TL
                            println("__breaking point__!")
                            break
                        end 
                    end
                end
            end  
            if candSI == []
                push!(IGPair,[I,G]); 
                @goto NewIter3
            else
                SI = nextSI(candSI,SI)
                # if SI∉collect(values(X))
                #     push!(exploredSI,SI);
                # end
            end
            if newsol == 0
                nosol+=1
            end
        end
        push!(IGPair,[I,G]);
    end
    return X,Y,main_iter,main_infeasi,left1_iter,left1_infeasi,left2_iter,left2_infeasi     #collect(values(X)),collect(values(Y))
end
PR(dX,dY,len,3);
PRtime = @CPUelapsed px,py,miter,minf,l1iter,l1inf,l2iter,l2inf = PR(dX,dY,len,TL)
minf/miter,l1inf/l1iter,l2inf/l2iter
nd = SortingSol(px,py)
# println("PRtime: ", PRtime, "PRsol: ", length(nd.Y))
py1= [nd.Y[i][1] for i=1:length(nd.Y)]; py2 = [nd.Y[i][2] for i=1:length(nd.Y)]
t5 = scatter(x=py1, y=py2, fname="LP+FP+FPP+PR", mode="markers", marker=attr(color="blue"))
plot([t1,t5],layout)

function SelectIG(ndlist)
    dist = []
    for i=1:length(ndlist)-1 
        if i==1 || i==length(ndlist)-1
            d = sqrt(sum((ndlist[i].-ndlist[i+1]).^2))
            push!(dist,d)
        else
            d1 = sqrt(sum((ndlist[i].-ndlist[i-1]).^2))
            d2 = sqrt(sum((ndlist[i].-ndlist[i+1]).^2))
            push!(dist, max(d1,d2))
        end
    end
    id = findmax(dist)[2]
    I0 = ndlist[id]; #G0 = ndlist[id+1]
    
    return I0#,G0
    
end

function PR2(dX,dY,len,TL)
    X = copy(dX); Y = copy(dY); 
    IGPair=[]; t0=time(); 
    println("randomly choose IG pairs")
    pc1 =0.3; iter1 = 0; infeasi1 = 0;
    while time()-t0 < TL*pc1
        nd = SortingSol(X,Y)
        if length(IGPair)>=(length(nd.Y)*(length(nd.Y)-1))
            return nd
        end
        @label NewIter1
        I0,G0 = StatsBase.StatsBase.sample(nd.Y, 2, replace=false);
        I = findall(i->i==I0,Y)[1]; G = findall(i->i==G0,Y)[1]
        # if [I,G] ∈ IGPair || abs(Y[I][1]-Y[G][1]) < 10^4|| abs(Y[I][2]-Y[G][2]) < 10^4 
        #     @goto NewIter1
        # end
        SI = X[I]; SG = X[G]; exploredSI = []; nosol = 0
        while nosol < 5 && all.(SI != SG) && time()-t0 < TL*pc1
            dif = findall(i-> SI[1:len[1]][i]!=SG[1:len[1]][i], 1:len[1])
            # weight = round(Int,mean([nd.Y[i][1]/nd.Y[i][2] for i=1:length(nd.Y)]))
            new_weight = round(Int, (Y[I][1]+Y[G][1])/(Y[I][2]+Y[G][2]))
            neibour = createAllNB(SI[1:len[1]],dif,exploredSI)
            newsol = 0; 
            if (length(neibour)==0) #(time()-t0 >= TL)
                @goto NewIter1
            else
                candSI =[]
                neibour2 = PR_preFBcheck(neibour)
                # println("removed neighbour: ", length(neibour)-length(neibour2))
                for l=1:length(neibour2)
                    iter1+=1
                    # st = PR_FBcheck( prmodel, main_iter, neibour[l])
                    st = PR_FBcheck_w( prmodel, new_weight, neibour2[l])
                    if st==true
                        sol = value.(all_variables(prmodel)); ndp = getobjval(sol)
                        push!(candSI,sol)
                        if sol ∉ nd.X && dominated(ndp,nd.Y)==false
                            push!(X, length(X)+1 => sol); push!(Y, length(Y)+1 => ndp)                            
                            newsol+=1;
                            # println(main_iter, " new sol ", ndp);
                        end
                    else
                        infeasi1+=1
                    end
                    if time()-t0 >= TL*pc1
                        println("__breaking point__!")
                        break
                    end
                end
            end
            if candSI == []
                push!(IGPair,[I,G]);
                @goto NewIter1
            else
                SI = nextSI(candSI,SI)
                if SI∉collect(values(X))
                    push!(exploredSI,SI);
                end
            end
            if newsol==0;
                nosol+=1
            end
        end
        push!(IGPair,[I,G]);
    end

    println("Choose I & G based on distance")
    pc = 0.8; iter2 = 0; infeasi2 = 0;

    while time()-t0 < TL*pc
        nd = SortingSol(X,Y)
        if length(IGPair)>=(length(nd.Y)*(length(nd.Y)-1))
            return X,Y
        end
        @label NewIter2
        I0 = SelectIG(nd.Y)
        I = findall(i->i==I0,Y)[1]; G = I+1 #findall(i->i==G0,Y)[1]
        # if [I,G] ∈ IGPair || abs(Y[I][1]-Y[G][1]) < 10^4|| abs(Y[I][2]-Y[G][2]) < 10^4 
        #     @goto NewIter1
        # end
        SI = X[I]; SG = X[G];
        nosol = 0; exploredSI = [] #copy(collect(values(X)));

        while nosol < 5 && all.(SI != SG) && time()-t0 < TL*pc
            dif = findall(i-> SI[1:len[1]][i]!=SG[1:len[1]][i], 1:len[1])
            # weight = round(Int,mean([nd.Y[i][1]/nd.Y[i][2] for i=1:length(nd.Y)]))
            new_weight = round(Int, (Y[I][1]+Y[G][1])/(Y[I][2]+Y[G][2]))
            neibour = createAllNB(SI[1:len[1]],dif,exploredSI)
            newsol = 0; 
            if (length(neibour)==0) #(time()-t0 >= TL)
                @goto NewIter2
            else
                candSI =[]
                neibour2 = PR_preFBcheck(neibour)
                # println("removed neighbour: ", length(neibour)-length(neibour2))
                for l=1:length(neibour2)
                    iter2+=1
                    # st = PR_FBcheck( prmodel, main_iter, neibour[l])
                    st = PR_FBcheck_w( prmodel, new_weight, neibour2[l])
                    if st==true
                        sol = value.(all_variables(prmodel)); ndp = getobjval(sol)
                        push!(candSI,sol)
                        if sol ∉ nd.X && dominated(ndp,nd.Y)==false
                            push!(X, length(X)+1 => sol); push!(Y, length(Y)+1 => ndp)                            
                            newsol+=1;
                            # println(main_iter, " new sol ", ndp);
                        end
                    else
                        infeasi2+=1
                    end
                    if time()-t0 >= TL*pc
                        println("__breaking point__!")
                        break
                    end
                end
            end
            if candSI == []
                push!(IGPair,[I,G]);
                @goto NewIter2
            else
                SI = nextSI(candSI,SI)
                if SI∉collect(values(X))
                    push!(exploredSI,SI);
                end
            end
            if newsol==0;
                nosol+=1
            end
        end
        push!(IGPair,[I,G]);
    end
    println("____1st Left____")
    left1_iter = 0; left1_infeasi = 0;
    pc2 =0.9
    while time()-t0 < TL*pc2
        nd = SortingSol(X,Y)
        if length(IGPair)>=(length(nd.Y)*(length(nd.Y)-1))
            return nd
        end
        @label Newiter3
        I0,G0 = StatsBase.StatsBase.sample(1:round(Int,length(nd.Y)*0.2), 2, replace=false);
        I = findall(i->i==nd.Y[I0],Y)[1]; G = findall(i->i==nd.Y[G0],Y)[1]
        # if [I,G] ∈ IGPair || abs(Y[I][1]-Y[G][1]) < 10^5|| abs(Y[I][2]-Y[G][2]) < 10^5
        #     @goto Newiter3
        # end
        SI = X[I]; SG = X[G];
        nosol = 0; exploredSI = []#copy(collect(values(X)));
        while nosol < 5 && all.(SI != SG) && time()-t0 < TL*pc2
            dif = findall(i-> SI[1:len[1]][i]!=SG[1:len[1]][i], 1:len[1])
            new_weight = round(Int, (Y[I][1]+Y[G][1])/(Y[I][2]+Y[G][2]))
            neibour = createNB(SI[1:len[1]],dif,exploredSI)
            newsol = 0;
            if (length(neibour)==0) 
                @goto Newiter3
            else
                candSI =[]
                neibour2 = PR_preFBcheck(neibour)
                # println("removed neighbour: ", length(neibour)-length(neibour2))
                for l=1:length(neibour2)
                    left1_iter+=1
                    # st = PR_FBcheck( prmodel, left1_iter, neibour[l])
                    st = PR_FBcheck_w( prmodel, new_weight, neibour2[l])
                    if st==true
                        sol = value.(all_variables(prmodel)); ndp = getobjval(sol)
                        push!(candSI,sol)
                        if sol ∉ nd.X && dominated(ndp,nd.Y)==false
                            push!(X, length(X)+1 => sol); push!(Y, length(Y)+1 => ndp)                            
                            newsol+=1; 
                            # println(iter, " new sol ", ndp);

                        end
                    else
                        left1_infeasi+=1
                    end
                    if time()-t0 >= TL*pc2
                        println("__breaking point__!")
                        break
                    end
                end
            end  
            if candSI == []
                push!(IGPair,[I,G]); 
                @goto Newiter3
            else
                SI = nextSI(candSI,SI)
                if SI∉collect(values(X))
                    push!(exploredSI,SI);
                end
            end
            if newsol == 0
                nosol+=1
            end
        end
        push!(IGPair,[I,G]);
    end

    println("____2nd Left____")
    left2_iter = 0; left2_infeasi = 0;
    while time()-t0 < TL
        nd = SortingSol(X,Y)
        if length(IGPair)>=(length(nd.Y)*(length(nd.Y)-1))
            return nd
        end
        @label Newiter4
        buffer = round(Int,length(nd.Y)*0.15)
        I0,G0 = StatsBase.StatsBase.sample(buffer:round(Int,length(nd.Y)*0.3), 2, replace=false);
        I = findall(i->i==nd.Y[I0],Y)[1]; G = findall(i->i==nd.Y[G0],Y)[1]
        # if [I,G] ∈ IGPair || abs(Y[I][1]-Y[G][1]) < 10^5|| abs(Y[I][2]-Y[G][2]) < 10^5 
        #     @goto Newiter4
        # end
        SI = X[I]; SG = X[G];
        nosol = 0; exploredSI = []#copy(collect(values(X)));
        while nosol < 5 && all.(SI != SG)&& time()-t0 < TL
            dif = findall(i-> SI[1:len[1]][i]!=SG[1:len[1]][i], 1:len[1])
            # weight = round(Int,mean([nd.Y[i][1]/nd.Y[i][2] for i=1:length(nd.Y)]))
            new_weight = round(Int, (Y[I][1]+Y[G][1])/(Y[I][2]+Y[G][2]))
            neibour = createNB(SI[1:len[1]],dif,exploredSI)
            newsol = 0;
            if (length(neibour)==0) 
                @goto Newiter4
            else
                candSI =[]
                neibour2 = PR_preFBcheck(neibour)
                # println("removed neighbour: ", length(neibour)-length(neibour2))
                for l=1:length(neibour2)
                    left2_iter+=1
                    # st = PR_FBcheck( prmodel, left1_iter, neibour[l])
                    st = PR_FBcheck_w( prmodel, new_weight, neibour2[l])
                    if st==true
                        sol = value.(all_variables(prmodel)); ndp = getobjval(sol)
                        push!(candSI,sol)
                        if sol ∉ nd.X && dominated(ndp,nd.Y)==false
                            push!(X, length(X)+1 => sol); push!(Y, length(Y)+1 => ndp)                            
                            newsol+=1; 
                            # println(iter, " new sol ", ndp);

                        end
                    else
                        left2_infeasi+=1
                    end
                    if time()-t0 >= TL
                        println("__breaking point__!")
                        break
                    end
                end
            end  
            if candSI == []
                push!(IGPair,[I,G]); 
                @goto Newiter4
            else
                SI = nextSI(candSI,SI)
                if SI∉collect(values(X))
                    push!(exploredSI,SI);
                end
            end
            if newsol == 0
                nosol+=1
            end
        end
        push!(IGPair,[I,G]);
    end
  
    return X,Y,iter1,infeasi1,iter2,infeasi2 #collect(values(X)),collect(values(Y))
end
############################################################################## 
PR2(dX,dY,len,3);
PPtime = @CPUelapsed tx,ty,it1,infi1,it2,infi2 = PR2(dX,dY,len,TL)
nd2 = SortingSol(tx,ty)
py21= [nd2.Y[i][1] for i=1:length(nd2.Y)]; py22 = [nd2.Y[i][2] for i=1:length(nd2.Y)]
t6 = scatter(x=py21, y=py22, fname="LP+FP+FPP+PR", mode="markers", marker=attr(color="green"))
plot([t1,t6],layout)


########################## Saving the output file ###########################
otable = zeros(length(nd.Y),2)
for i=1:length(nd.Y)
    for j=1:2
        otable[i,j] = nd.Y[i][j]
    end
end

CSV.write("/home/k2g00/k2g3475/scnd/vopt/lpY/"*fname*"lpY.log", DataFrame(otable, :auto), append=false, header=false,delim=' ')
CSV.write("/home/ak121396/Desktop/relise/lpY/ndp/"*fname*"lpY.log", DataFrame(otable, :auto), append=false, header=false,delim=' ')
println("algotime $fname: ", LPtime+FPtime+FPPtime+PRtime+l1time+l2time," #sol: ", length(pry))

dv = sparse.(nd.X)
JLD2.@save "/home/k2g00/k2g3475/scnd/vopt/lpY/X/"*fname*"X.jld2" dv
# lexsol = load("/home/k2g00/k2g3475/scnd/vopt/lpY/X/"*fname*"X.jld2")
# fname = file[end-7:end]
JLD2.@save "/home/ak121396/Desktop/relise/lpY/X/"*fname*"X.jld2" dv 
JLD2.@load "/home/ak121396/Desktop/relise/lpY/X/"*fname*"X.jld2" dv
###############################
  # println("____1st Left____")
    # pc2 = 0.8
    # while time()-t0 < TL*pc2
    #     nd = SortingSol(X,Y)
    #     if length(IGPair)>=(length(nd.Y)*(length(nd.Y)-1))
    #         return nd
    #     end
    #     @label NewIter2
    #     I0,G0 = StatsBase.StatsBase.sample(1:round(Int,length(nd.Y)*0.2), 2, replace=false);
    #     I = findall(i->i==nd.Y[I0],Y)[1]; G = findall(i->i==nd.Y[G0],Y)[1]
    #     # if [I,G] ∈ IGPair || abs(Y[I][1]-Y[G][1]) < 10^5|| abs(Y[I][2]-Y[G][2]) < 10^5
    #     #     @goto NewIter2
    #     # end
    #     SI = X[I]; SG = X[G];
    #     iter=0; nosol = 0
    #     while nosol < 5 && all.(SI != SG) && time()-t0 < TL*pc2
    #         dif = findall(i-> SI[1:len[1]][i]!=SG[1:len[1]][i], 1:len[1])
    #         new_weight = round(Int, (Y[I][1]+Y[G][1])/(Y[I][2]+Y[G][2]))
    #         neibour = createNB(SI[1:len[1]],dif,exploredSI)
    #         newsol = 0;
    #         if (length(neibour)==0) 
    #             @goto NewIter2
    #         else
    #             candSI =[]
    #             for l=1:length(neibour)
    #                 # st = PR_FBcheck( prmodel, iter, neibour[l])
    #                 st = PR_FBcheck_w( prmodel, new_weight, neibour[l])
    #                 if st==true
    #                     sol = value.(all_variables(prmodel)); ndp = getobjval(sol)
    #                     push!(candSI,sol)
    #                     if sol ∉ nd.X && dominated(ndp,nd.Y)==false
    #                         push!(X, length(X)+1 => sol); push!(Y, length(Y)+1 => ndp)                            
    #                         newsol+=1; 
    #                         # println(iter, " new sol ", ndp);

    #                     end
    #                 end
    #                 if time()-t0 >= TL*pc2
    #                     println("__breaking point__!")
    #                     break
    #                 end
    #             end
    #         end  
    #         if candSI == []
    #             @goto NewIter2
    #         else
    #             SI = nextSI(candSI,SI)
    #             if SI∉collect(values(X))
    #                 push!(exploredSI,SI);
    #             end
    #         end
    #         if newsol == 0
    #             nosol+=1
    #         end
    #         iter+=1
    #     end
    #     push!(IGPair,[I,G]);
    # end

    # println("____2nd Left____")
    # while time()-t0 < TL
    #     nd = SortingSol(X,Y)
    #     if length(IGPair)>=(length(nd.Y)*(length(nd.Y)-1))
    #         return nd
    #     end
    #     @label NewIter3
    #     buffer = round(Int,length(nd.Y)*0.15)
    #     I0,G0 = StatsBase.StatsBase.sample(buffer:round(Int,length(nd.Y)*0.3), 2, replace=false);
    #     I = findall(i->i==nd.Y[I0],Y)[1]; G = findall(i->i==nd.Y[G0],Y)[1]
    #     # if [I,G] ∈ IGPair || abs(Y[I][1]-Y[G][1]) < 10^5|| abs(Y[I][2]-Y[G][2]) < 10^5 
    #     #     @goto NewIter3
    #     # end
    #     SI = X[I]; SG = X[G];
    #     iter=0; nosol = 0
    #     while nosol < 5 && all.(SI != SG)&& time()-t0 < TL
    #         dif = findall(i-> SI[1:len[1]][i]!=SG[1:len[1]][i], 1:len[1])
    #         # weight = round(Int,mean([nd.Y[i][1]/nd.Y[i][2] for i=1:length(nd.Y)]))
    #         new_weight = round(Int, (Y[I][1]+Y[G][1])/(Y[I][2]+Y[G][2]))
    #         neibour = createNB(SI[1:len[1]],dif,exploredSI)
    #         newsol = 0;
    #         if (length(neibour)==0) 
    #             @goto NewIter3
    #         else
    #             candSI =[]
    #             for l=1:length(neibour)
    #                 # st = PR_FBcheck( prmodel, iter, neibour[l])
    #                 st = PR_FBcheck_w( prmodel, new_weight, neibour[l])
    #                 if st==true
    #                     sol = value.(all_variables(prmodel)); ndp = getobjval(sol)
    #                     push!(candSI,sol)
    #                     if sol ∉ nd.X && dominated(ndp,nd.Y)==false
    #                         push!(X, length(X)+1 => sol); push!(Y, length(Y)+1 => ndp)                            
    #                         newsol+=1; 
    #                         # println(iter, " new sol ", ndp);

    #                     end
    #                 end
    #                 if time()-t0 >= TL
    #                     println("__breaking point__!")
    #                     break
    #                 end
    #             end
    #         end  
    #         if candSI == []
    #             @goto NewIter3
    #         else
    #             SI = nextSI(candSI,SI)
    #             if SI∉collect(values(X))
    #                 push!(exploredSI,SI);
    #             end
    #         end
    #         if newsol == 0
    #             nosol+=1
    #         end
    #         iter+=1
    #     end
    #     push!(IGPair,[I,G]);

    # end
###############################
# function Grouping(LB)
#     fmin = minimum([LB[i][1] for i=1:length(LB)]),minimum([LB[i][2] for i=1:length(LB)])
#     fmax = maximum([LB[i][1] for i=1:length(LB)]),maximum([LB[i][2] for i=1:length(LB)])
#     steps = [round.(Int,abs(fmax[k]-fmin[k])/9) for k=1:2]
#     cubes = Dict();
#     for iter=1:length(LB)
#         loca = [round.(Int,((LB[iter][k]-fmin[k])/steps[k])+1) for k=1:2]
#         if !haskey(cubes,loca)
#             cubes[loca] = [iter]
#         else
#             push!(cubes[loca], iter)
#         end
#     end
#     groups = collect(values(cubes)); groupkeys = collect(keys(cubes))
#     return cubes,groupkeys
# end

# function GFP(candX,candY,len,TL)
#     X = []; PF =[]; Tabu = [];  newsol = 0;
#     Y = []; U1 = []; U2= []; U3 = []; t0 = time();
#     cubes,gkeys = Grouping(candY)
#     glist = copy(gkeys)
#     while time() - t0 < TL && glist != []
#         g = StatsBase.sample(1:length(glist))
#         k = StatsBase.sample(cubes[glist[g]])
#         x_t = candX[k]
#         yt = x_t[1:len[1]]; yr = copy(yt)
#         # u1t = x_t[1+len[1]:len[1]+len[2]];
#         # u2t = x_t[1+sum(len[i] for i=1:2):sum(len[i] for i=1:3)]
#         # u3t = x_t[1+sum(len[i] for i=1:3):sum(len[i] for i=1:4)]
#         SearchDone = false; iter=0; Max_iter = 5#length(findall(i-> 0<i<1,yt))
#         while iter < Max_iter && SearchDone == false && time() - t0 < TL
#             yid = findall(p->p>0.5,yt);
#             for j=1:len[1]
#                 if j in yid
#                     yr[j]=1
#                 else
#                     yr[j]=0
#                 end
#             end
#             if FP_FBcheck(fpmodel,yr) == true #,u1r,u2r,u3r)
#                 sol = value.(all_variables(fpmodel)); ndp = getobjval(sol)
#                 if sol ∉ X  && dominated(ndp,PF)==false
#                     push!(X,sol); push!(PF,ndp)
#                     push!(Y,yr);
#                     newsol+=1; SearchDone = true
#                     deleteat!(glist, g);
#                     println("rounding worked")

#                 end
#             else
#                 if yr ∈ Tabu
#                     yr = flipoper(Y,yt,yr);
#                     if yr == []
#                         SearchDone = true;
#                         deleteat!(cubes[glist[g]], cubes[glist[g]] .== k);
#                         println("flip failed")
#                     else
#                         if FP_FBcheck(fpmodel,yr) == true #,u1r,u2r,u3r)
#                             sol = value.(all_variables(fpmodel)); ndp = getobjval(sol)
#                             if sol ∉ X && dominated(ndp,PF)==false
#                                 push!(X,sol); push!(PF,ndp)
#                                 push!(Y,yr);
#                                 newsol+=1; SearchDone = true;
#                                 deleteat!(glist, g);
#                                 println("flip worked")

#                             end
#                         end
#                     end
#                 end
#                 if time()-t0 >= TL
#                     break
#                 end
#                 if SearchDone == false
#                     push!(Tabu,yr)
#                     yt = fbsearch(yr)
#                     if yt==0  #when there's no new feasible lp sol
#                         deleteat!(cubes[glist[g]], cubes[glist[g]] .== k);
#                         println("no solution")
#                         SearchDone = true
#                     end
#                 end
#             end
# 			iter+=1
#         end
#     end
#     return X,PF,newsol,glist
# end
# FPtime = @CPUelapsed lx,ly,ln,glist = GFP(lp.X_E,lp.Y_N,len,round(Int,LPtime*3))
# gx,gy = Postpro(lx,ly)
# PRtime = @CPUelapsed px,py,pn,pairs = PR(lx,ly,len,round(Int,FPtime))
# prx,pry = Postpro(px,py)


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
#             jid = StatsBase.sample(jset); kid = StatsBase.sample(kset);
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
#                 if FP_FBcheck(fpmodel,yr,u1r,u2r,u3r) == true #,u1r,u2r,u3r)
#                     sol = value.(all_variables(fpmodel)); ndp = getobjval(sol)
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
#                             if FP_FBcheck(fpmodel,yr,u1r,u2r,u3r) == true #,u1r,u2r,u3r)
#                                 sol = value.(all_variables(fpmodel)); ndp = getobjval(sol)
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
