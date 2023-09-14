# cd("C:/Users/AK121396/Desktop/ProjectBenders")
using DataFrames,DelimitedFiles,JuMP,LinearAlgebra,CPLEX,MathOptInterface,SparseArrays
struct Data1d
    file::String; N::Dict{}; d::Array{}; c::Array{}; a::Array{}; e::Array{}; gij::SparseVector{}; gjk::SparseVector{}; gkl::SparseVector{};
    vij::Array{}; vjk::Array{}; vkl::Array{}; Vij::SparseVector{}; Vjk::SparseVector{}; Mij::Array{}; Mjk::Array{}; Mkl::Array{};
    b::Array{}; q::Array{}; rij::Array{}; rjk::Array{}; rkl::Array{}; upl::Int; udc::Int; bigM::Int
    function Data1d(file)
        dt1 = readdlm(file);
        # notafile = readdlm("/home/ak121396/Desktop/instances/scnd/Notations.txt", '=');
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
# @show file = ARGS[1];
# file = "F:/scnd/Test4S4"
# file = "/home/ak121396/Desktop/instances/scnd/test04S4"
file = "/home/k2g00/k2g3475/scnd/instances/test01S2"
dt1 = Data1d(file);
function SCND1dim(dt1,nobj)
    ##########################  Mathematical model  #########################
    # scnd1 = Model(CPLEX.Optimizer);
    scnd1 = Model(optimizer_with_attributes( CPLEX.Optimizer, "CPX_PARAM_EPGAP" => 1e-8 ));
    set_silent(scnd1)
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
    # @constraint(scnd1, obj1, sum(dt1.c.*y1) +sum(repeat(dt1.a[1,:], outer=sum(dt1.Mij[1,:])).*xij1[1:sum(dt1.Mij[1,:])*5])+
    #         sum(sum(repeat(dt1.a[i,:], outer=sum(dt1.Mij[i,:])).*xij1[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"])+
    #         sum(dt1.e.*h1) + sum(dt1.gij[i]*uij1[i] for i in findnz(dt1.gij)[1]) + sum(dt1.gjk[i]*ujk1[i] for i in findnz(dt1.gjk)[1]) + sum(dt1.gkl[i].*ukl1[i] for i in findnz(dt1.gkl)[1])+
    #         sum(dt1.vij.*xij1)+sum(dt1.vjk.*xjk1)+sum(dt1.vkl.*xkl1) <= 0)
    if nobj == 1
        @objective(scnd1, Min, sum(dt1.c.*y1) +
            sum(repeat(dt1.a[1,:], outer=sum(dt1.Mij[1,:])).*xij1[1:sum(dt1.Mij[1,:])*5])+
            sum(sum(repeat(dt1.a[i,:], outer=sum(dt1.Mij[i,:])).*xij1[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"])+
            sum(dt1.gij[i]*uij1[i] for i in findnz(dt1.gij)[1]) + sum(dt1.gjk[i]*ujk1[i] for i in findnz(dt1.gjk)[1]) + sum(dt1.gkl[i].*ukl1[i] for i in findnz(dt1.gkl)[1])+
            sum(dt1.e.*h1) + sum(dt1.vij.*xij1)+sum(dt1.vjk.*xjk1)+sum(dt1.vkl.*xkl1)
        )
    else
    # @constraint(scnd1, obj2, sum(repeat(dt1.b[1,:], outer=sum(dt1.Mij[1,:])).*xij1[1:sum(dt1.Mij[1,:])*5]) +
    #         sum(sum(repeat(dt1.b[i,:], outer=sum(dt1.Mij[i,:])).*xij1[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"]) +
    #         sum(dt1.q.*h1) + sum(dt1.rij.*xij1)+sum(dt1.rjk.*xjk1)+sum(dt1.rkl.*xkl1) <= 1.1826921599716e6 )
        @objective(scnd1, Min, sum(repeat(dt1.b[1,:], outer=sum(dt1.Mij[1,:])).*xij1[1:sum(dt1.Mij[1,:])*5]) +
            sum(sum(repeat(dt1.b[i,:], outer=sum(dt1.Mij[i,:])).*xij1[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"]) +
            sum(dt1.q.*h1) + sum(dt1.rij.*xij1)+sum(dt1.rjk.*xjk1)+sum(dt1.rkl.*xkl1));
    end
    # @objective(scnd1, Min, w[1]*(sum(dt1.c.*y1) +sum(repeat(dt1.a[1,:], outer=sum(dt1.Mij[1,:])).*xij1[1:sum(dt1.Mij[1,:])*5])+
    #     sum(sum(repeat(dt1.a[i,:], outer=sum(dt1.Mij[i,:])).*xij1[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"])+
    #     sum(dt1.e.*h1) + sum(dt1.gij[i]*uij1[i] for i in findnz(dt1.gij)[1]) + sum(dt1.gjk[i]*ujk1[i] for i in findnz(dt1.gjk)[1]) + sum(dt1.gkl[i].*ukl1[i] for i in findnz(dt1.gkl)[1])+
    #     sum(dt1.vij.*xij1)+sum(dt1.vjk.*xjk1)+sum(dt1.vkl.*xkl1)) +
    #     w[2]*(sum(repeat(dt1.b[1,:], outer=sum(dt1.Mij[1,:])).*xij1[1:sum(dt1.Mij[1,:])*5]) +
    #     sum(sum(repeat(dt1.b[i,:], outer=sum(dt1.Mij[i,:])).*xij1[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"]) +
    #     sum(dt1.q.*h1) + sum(dt1.rij.*xij1)+sum(dt1.rjk.*xjk1)+sum(dt1.rkl.*xkl1))
    # );
    ########## constraint 3 #############
    @constraint(scnd1, [p=1:5], sum(xij1[5*(m-1)+p] for m=1:dt1.Mij[1,1])+sum(xij1[5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,1]) == sum(xjk1[5*(m-1)+p] for m=1:sum(dt1.Mjk[1,:])) );
    @constraint(scnd1, [j=2:dt1.N["plant"],p=1:5], sum(xij1[5*sum(dt1.Mij[1,1:j-1])+5*(m-1)+p] for m=1:dt1.Mij[1,j])+sum(xij1[5*sum(dt1.Mij[i,1:j-1])+5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,j]) == sum(xjk1[sum(dt1.Mjk[1:j-1,:])*5 + 5*(m-1)+p] for m=1:sum(dt1.Mjk[j,:])) );
    @constraint(scnd1, [p=1:5], sum(xjk1[5*(m-1)+p] for m=1:dt1.Mjk[1,1])+sum(xjk1[5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,1]) == sum(xkl1[5*(m-1)+p] for m=1:sum(dt1.Mkl[1,:])) );
    @constraint(scnd1, [k=2:dt1.N["distribution"],p=1:5],sum(xjk1[5*sum(dt1.Mjk[1,1:k-1])+5*(m-1)+p] for m=1:dt1.Mjk[1,k])+sum(xjk1[5*sum(dt1.Mjk[j,1:k-1])+5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,k]) == sum(xkl1[sum(dt1.Mkl[1:k-1,:])*5 + 5*(m-1)+p] for m=1:sum(dt1.Mkl[k,:])) );
    ########### constraint 4-6 #############
    # @constraints(scnd1, begin
    # @constraint(scnd1, [p=1:5],sum(h1[5*(t-1)+p] for t=1:2) == sum(xij1[5*(m-1)+p] for m=1:dt1.Mij[1,1])+sum(xij1[5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,1]))
    # @constraint(scnd1, [j=2:dt1.N["plant"],p=1:5], sum(h1[5*2*(j-1)+5*(t-1)+p] for t=1:2) == sum(xij1[5*(j-1)+5*(m-1)+p] for m=1:dt1.Mij[1,j])+sum(xij1[5*(j-1)+5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,j]))
    # @constraint(scnd1, [p=1:5], sum(h1[5*2*dt1.N["plant"]+5*(t-1)+p] for t=1:2) == sum(xjk1[5*(m-1)+p] for m=1:dt1.Mjk[1,1])+sum(xjk1[5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,1]));
    # @constraint(scnd1, [k=2:dt1.N["distribution"],p=1:5], sum(h1[5*2*dt1.N["plant"]+5*2*(k-1)+5*(t-1)+p] for t=1:2) == sum(xjk1[5*(k-1)+5*(m-1)+p] for m=1:dt1.Mjk[1,k])+sum(xjk1[5*(k-1)+5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,k]))
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
    @constraints(scnd1, begin
        sum(uij1[1:dt1.Mij[1,1]]) <= sum(y1[1:2])
        sum(uij1[sum(dt1.Mij[1,:])+dt1.Mij[2,1]]) <= sum(y1[3:4])
        [j=2:dt1.N["plant"]], sum(uij1[sum(dt1.Mij[1,1:j-1])+1:sum(dt1.Mij[1,1:j-1])+dt1.Mij[1,j]]) <= sum(y1[2*(j-1)+1:2*(j-1)+2])
        [i=2:dt1.N["supplier"],j=2:dt1.N["plant"]],  sum(uij1[sum(dt1.Mij[1:i-1,:])+sum(dt1.Mij[i,1:j-1])+1:sum(dt1.Mij[1:i-1,:])+sum(dt1.Mij[i,1:j-1])+dt1.Mij[i,j]])<= sum(y1[2*(j-1)+1:2*(j-1)+2])
        sum(ujk1[1:dt1.Mjk[1,1]]) <= (sum(y1[1:2])+sum(y1[dt1.N["plant"]+1:dt1.N["plant"]+2]))/2
        sum(ujk1[sum(dt1.Mjk[1,:])+dt1.Mjk[2,1]]) <= (sum(y1[3:4])+sum(y1[dt1.N["plant"]+1:dt1.N["plant"]+2]))/2
        [k=2:dt1.N["distribution"]], sum(ujk1[sum(dt1.Mjk[1,1:k-1])+1:sum(dt1.Mjk[1,1:k-1])+dt1.Mjk[1,k]]) <= (sum(y1[1:2])+sum(y1[dt1.N["plant"] + 2*(k-1)+1:dt1.N["plant"] + 2*(k-1)+2]))/2
        [j=2:dt1.N["plant"],k=2:dt1.N["distribution"]],  sum(ujk1[sum(dt1.Mjk[1:j-1,:])+sum(dt1.Mjk[j,1:j-1])+1:sum(dt1.Mjk[1:j-1,:])+sum(dt1.Mjk[j,1:j-1])+dt1.Mjk[j,k]]) <= (sum(y1[2*(j-1)+1:2*(j-1)+2])+sum(y1[dt1.N["plant"] + 2*(k-1)+1:dt1.N["plant"] + 2*(k-1)+2]))/2
        sum(ukl1[1:dt1.Mkl[1,1]]) <= sum(y1[dt1.N["plant"]+1:dt1.N["plant"]+2]) 
        sum(ukl1[sum(dt1.Mkl[1,:])+dt1.Mkl[2,1]]) <= sum(y1[dt1.N["plant"]+3:dt1.N["plant"]+4])
        [l=2:dt1.N["customer"]], sum(ukl1[sum(dt1.Mkl[1,1:l-1])+1:sum(dt1.Mkl[1,1:l-1])+dt1.Mkl[1,l]]) <= sum(y1[dt1.N["plant"]+1:dt1.N["plant"]+2])
        [k=2:dt1.N["distribution"],l=2:dt1.N["customer"]],  sum(ukl1[sum(dt1.Mkl[1:k-1,:])+sum(dt1.Mkl[k,1:l-1])+1:sum(dt1.Mkl[1:k-1,:])+sum(dt1.Mkl[k,1:l-1])+dt1.Mkl[k,l]])<= sum(y1[dt1.N["plant"]+ 2*(k-1)+1:dt1.N["plant"]+ 2*(k-1)+2])
    end);
    # @constraints(scnd1, begin
    #     sum(uij1[1:dt1.Mij[1,1]]) <= 1
    #     sum(uij1[sum(dt1.Mij[1,:])+dt1.Mij[2,1]]) <= 1
    #     [j=2:dt1.N["plant"]], sum(uij1[sum(dt1.Mij[1,1:j-1])+1:sum(dt1.Mij[1,1:j-1])+dt1.Mij[1,j]]) <= 1
    #     [i=2:dt1.N["supplier"],j=2:dt1.N["plant"]],  sum(uij1[sum(dt1.Mij[1:i-1,:])+sum(dt1.Mij[i,1:j-1])+1:sum(dt1.Mij[1:i-1,:])+sum(dt1.Mij[i,1:j-1])+dt1.Mij[i,j]])<= 1
    #     sum(ujk1[1:dt1.Mjk[1,1]]) <= 1
    #     sum(ujk1[sum(dt1.Mjk[1,:])+dt1.Mjk[2,1]]) <= 1
    #     [k=2:dt1.N["distribution"]], sum(ujk1[sum(dt1.Mjk[1,1:k-1])+1:sum(dt1.Mjk[1,1:k-1])+dt1.Mjk[1,k]]) <= 1
    #     [j=2:dt1.N["plant"],k=2:dt1.N["distribution"]],  sum(ujk1[sum(dt1.Mjk[1:j-1,:])+sum(dt1.Mjk[j,1:j-1])+1:sum(dt1.Mjk[1:j-1,:])+sum(dt1.Mjk[j,1:j-1])+dt1.Mjk[j,k]]) <= 1
    #     sum(ukl1[1:dt1.Mkl[1,1]]) <= 1
    #     sum(ukl1[sum(dt1.Mkl[1,:])+dt1.Mkl[2,1]]) <= 1
    #     [l=2:dt1.N["customer"]], sum(ukl1[sum(dt1.Mkl[1,1:l-1])+1:sum(dt1.Mkl[1,1:l-1])+dt1.Mkl[1,l]]) <= 1
    #     [k=2:dt1.N["distribution"],l=2:dt1.N["customer"]],  sum(ukl1[sum(dt1.Mkl[1:k-1,:])+sum(dt1.Mkl[k,1:l-1])+1:sum(dt1.Mkl[1:k-1,:])+sum(dt1.Mkl[k,1:l-1])+dt1.Mkl[k,l]])<= 1
    # end);
    ########### constraint 11 ############# This causes a different obj values
    @constraint(scnd1, [i=1:sum(dt1.Mij)], sum(xij1[5*(i-1)+1:5*i]) <= dt1.bigM*uij1[i])
    @constraint(scnd1, [j=1:sum(dt1.Mjk)], sum(xjk1[5*(j-1)+1:5*j]) <= dt1.bigM*ujk1[j])
    @constraint(scnd1, [k=1:sum(dt1.Mkl)], sum(xkl1[5*(k-1)+1:5*k]) <= dt1.bigM*ukl1[k])
    # ########### constraint 12 #############
    @constraints(scnd1, begin
        [i in findnz(dt1.Vij)[1]], sum(xij1[5*(i-1)+1:5*i]) >= dt1.Vij[i]*uij1[i]
        [i in findnz(dt1.Vjk)[1]], sum(xjk1[5*(i-1)+1:5*i]) >= dt1.Vjk[i]*ujk1[i]
        # [i in findnz(dt1.Vkl)[1]], sum(xkl1[5*(i-1)+1:5*i]) >= dt1.Vkl[i]*ukl1[i]
    end);
    # ########### constraint 13-14 #############
    @constraint(scnd1, sum(y1[1:dt1.N["plant"]*2]) <= dt1.upl);
    @constraint(scnd1, sum(y1[dt1.N["plant"]*2+1:end]) <= dt1.udc);
    return scnd1
end
scnd1 = SCND1dim(dt1,1);
optimize!(scnd1); objective_value(scnd1)
#solve_time(scnd1)
# termination_status(scnd1)
#write_to_file(scnd1, "/home/ak121396/Desktop/instances/SCND/small/test1s2_obj1.lp")



###############################   small SCND model: Data from smSCNDGenerator.jl   ####################################
include("./smSCNDGenerator.jl")
function get1d_objval(model)
    y = value.(model[:y])
    uij = value.(model[:uij])
    ujk = value.(model[:ujk])
    ukl = value.(model[:ukl])
    xij = value.(model[:xij])
    xjk = value.(model[:xjk])
    xkl = value.(model[:xkl])
    h = value.(model[:h])

    obj1 =  sum(c2.*y) +sum(repeat(vcs[1,:], outer=sum(Mij[1,:])).*xij[1:sum(Mij[1,:])]) +
        sum(sum(repeat(vcs[i,:], outer=sum(Mij[i,:])).*xij[sum(Mij[1:i-1,:])+1:sum(Mij[1:i,:])]) for i=2:I) +
        dot(e2,h) +sum(gij2[i]*uij[i] for i in findnz(gij2)[1]) + sum(gjk2[i]*ujk[i] for i in findnz(gjk2)[1]) + sum(gkl2[i].*ukl[i] for i in findnz(gkl2)[1]) +
        sum(tcp2.*xij)+sum(tcd2.*xjk)+sum(tcc2.*xkl)
    obj2 = sum(repeat(ves[1,:], outer=sum(Mij[1,:])).*xij[1:sum(Mij[1,:])]) +
        sum(sum(repeat(ves[i,:], outer=sum(Mij[i,:])).*xij[sum(Mij[1:i-1,:])+1:sum(Mij[1:i,:])]) for i=2:I) +
        dot(q2,h) + sum(rij2.*xij)+sum(rjk2.*xjk)+sum(rkl2.*xkl)
    return [obj1,obj2]
end

function sm1dModel(w)
    model = Model(CPLEX.Optimizer); set_silent(model)
    @variable(model, y[1:(J+K)*2], Bin)
    @variable(model, uij[1:sum(Mij)], Bin);
    @variable(model, ujk[1:sum(Mjk)], Bin);
    @variable(model, ukl[1:sum(Mkl)], Bin);
    @variable( model, 0<= xij[1:sum(Mij)] );
    @variable( model, 0<= xjk[1:sum(Mjk)] );
    @variable( model, 0<= xkl[1:sum(Mkl)] );
    @variable( model, 0<= h[1:(J+K)*2] );
    @objective(model, Min,  w[1]*(sum(c2.*y) +sum(repeat(vcs[1,:], outer=sum(Mij[1,:])).*xij[1:sum(Mij[1,:])]) +
        sum(sum(repeat(vcs[i,:], outer=sum(Mij[i,:])).*xij[sum(Mij[1:i-1,:])+1:sum(Mij[1:i,:])]) for i=2:I) +
        dot(e2,h) +sum(gij2[i]*uij[i] for i in findnz(gij2)[1]) + sum(gjk2[i]*ujk[i] for i in findnz(gjk2)[1]) + sum(gkl2[i].*ukl[i] for i in findnz(gkl2)[1]) +
        sum(tcp2.*xij)+sum(tcd2.*xjk)+sum(tcc2.*xkl))    +
            w[2]*(sum(repeat(ves[1,:], outer=sum(Mij[1,:])).*xij[1:sum(Mij[1,:])]) +
        sum(sum(repeat(ves[i,:], outer=sum(Mij[i,:])).*xij[sum(Mij[1:i-1,:])+1:sum(Mij[1:i,:])]) for i=2:I) +
        dot(q2,h) + sum(rij2.*xij)+sum(rjk2.*xjk)+sum(rkl2.*xkl))
    );
    ########## constraint 3 #############
    @constraint(model, sum(xij[m] for m=1:Mij[1,1])+ sum(xij[m+(sum(Mij[1:i-1,:]))] for i=2:I for m=1:Mij[i,1])  == sum(xjk[m] for m=1:sum(Mjk[1,:])) )
    @constraint(model, [j=2:J], sum(xij[sum(Mij[1,1:j-1])+m] for m=1:Mij[1,j])+sum(xij[sum(Mij[i,1:j-1])+m+(sum(Mij[1:i-1,:]))] for i=2:I for m=1:Mij[i,j]) == sum(xjk[sum(Mjk[1:j-1,:]) + m] for m=1:sum(Mjk[j,:])) )
    @constraint(model, sum(xjk[m] for m=1:Mjk[1,1])+sum(xjk[m+(sum(Mjk[1:j-1,:]))] for j=2:J for m=1:Mjk[j,1]) == sum(xkl[m] for m=1:sum(Mkl[1,:])) )
    @constraint(model, [k=2:K],sum(xjk[sum(Mjk[1,1:k-1])+m] for m=1:Mjk[1,k])+sum(xjk[sum(Mjk[j,1:k-1])+m+(sum(Mjk[1:j-1,:]))] for j=2:J for m=1:Mjk[j,k]) == sum(xkl[sum(Mkl[1:k-1,:]) + m] for m=1:sum(Mkl[k,:])) )
    ########### constraint 4-6 #############
    @constraint(model, sum(h[t] for t=1:2) ==sum(xij[m] for m=1:Mij[1,1])+sum(xij[m+(sum(Mij[1:i-1,:]))] for i=2:I for m=1:Mij[i,1]));
    @constraint(model, [j=2:J], sum(h[2*(j-1)+t] for t=1:2) == sum(xij[sum(Mij[1,1:j-1])+m] for m=1:Mij[1,j])+sum(xij[sum(Mij[i,1:j-1])+m+(sum(Mij[1:i-1,:]))] for i=2:I for m=1:Mij[i,j]) );
    @constraint(model,sum(h[2*J+t] for t=1:2) == sum(xjk[m] for m=1:Mjk[1,1])+sum(xjk[m+(sum(Mjk[1:j-1,:]))] for j=2:J for m=1:Mjk[j,1]) );
    @constraint(model, [k=2:K], sum(h[2*J+2*(k-1)+t] for t=1:2) == sum(xjk[sum(Mjk[1,1:k-1])+m] for m=1:Mjk[1,k])+sum(xjk[sum(Mjk[j,1:k-1])+m+(sum(Mjk[1:j-1,:]))] for j=2:J for m=1:Mjk[j,k]));
    @constraint(model, sum(xkl[m] for m=1:Mkl[1,1]) +sum(xkl[m+(sum(Mkl[1:k-1,:]))] for k=2:K for m=1:Mkl[k,1]) >= demand[1]);
    @constraint(model, [l=2:L], sum(xkl[sum(Mkl[1,1:l-1]) + m] for m=1:Mkl[1,l])+ sum(xkl[sum(Mkl[1:k-1,:])+sum(Mkl[k,1:l-1])+m] for k=2:K for m=1:Mkl[k,l]) >= demand[l]);
    ########### constraint 7 #############
    @constraint(model, sum(xij[1:sum(Mij[1,:])]) <= cas1[1]);
    @constraint(model, [i=2:I],  sum(xij[sum(Mij[1:i-1,:])+1:sum(Mij[1:i,:])]) <= cas1[i]);
    ########### constraint 8 #############
    @constraint(model,[j=1:J+K, t=1:2], sum(h[2*(j-1)+t]) <= capd1[j]*y[2*(j-1)+t]);
    ########### constraint 9 #############
    @constraint(model,[j=1:J+K], sum(y[2*(j-1)+1:2*(j-1)+2]) <= 1);
    ########### constraint 10 #############
    @constraints(model,begin
        # [i=1,j=1], 
        sum(uij[1:Mij[1,1]]) <= sum(y[1:2])
        # [i=2,j=1], 
        sum(uij[sum(Mij[1,:])+Mij[2,1]]) <= sum(y[3:4]) #sum(y[2*(j-1)+1:2*(j-1)+2])
        [j=2:J], sum(uij[sum(Mij[1,1:j-1])+1:sum(Mij[1,1:j-1])+Mij[1,j]]) <= sum(y[2*(j-1)+1:2*(j-1)+2])
        [i=2:I,j=2:J],  sum(uij[sum(Mij[1:i-1,:])+sum(Mij[i,1:j-1])+1:sum(Mij[1:i-1,:])+sum(Mij[i,1:j-1])+Mij[i,j]])<= sum(y[2*(j-1)+1:2*(j-1)+2])
        
        # [j=1,k=1], 
        sum(ujk[1:Mjk[1,1]]) <= (sum(y[1:2])+sum(y[J+1:J+2]))/2
        # [j=2,k=1], 
        sum(ujk[sum(Mjk[1,:])+Mjk[2,1]]) <= (sum(y[3:4])+sum(y[J+1:J+2]))/2
        [k=2:K], sum(ujk[sum(Mjk[1,1:k-1])+1:sum(Mjk[1,1:k-1])+Mjk[1,k]]) <= (sum(y[1:2])+sum(y[J + 2*(k-1)+1:J + 2*(k-1)+2]))/2
        [j=2:J,k=2:K],  sum(ujk[sum(Mjk[1:j-1,:])+sum(Mjk[j,1:j-1])+1:sum(Mjk[1:j-1,:])+sum(Mjk[j,1:j-1])+Mjk[j,k]]) <= (sum(y[2*(j-1)+1:2*(j-1)+2])+sum(y[J + 2*(k-1)+1:J + 2*(k-1)+2]))/2
        
        # [k=1,l=1], 
        sum(ukl[1:Mkl[1,1]]) <= sum(y[J+1:J+2]) #sum(y[J+ 2*(k-1)+1:J+ 2*(k-1)+2])
        # [k=2,l=1], 
        sum(ukl[sum(Mkl[1,:])+Mkl[2,1]]) <= sum(y[J+3:J+4]) #sum(y[J+ 2*(k-1)+1:J+ 2*(k-1)+2])
        [l=2:L], sum(ukl[sum(Mkl[1,1:l-1])+1:sum(Mkl[1,1:l-1])+Mkl[1,l]]) <= sum(y[J+1:J+2])
        [k=2:K,l=2:L],  sum(ukl[sum(Mkl[1:k-1,:])+sum(Mkl[k,1:l-1])+1:sum(Mkl[1:k-1,:])+sum(Mkl[k,1:l-1])+Mkl[k,l]])<= sum(y[J+ 2*(k-1)+1:J+ 2*(k-1)+2])
    end); 
    ########### constraint 11 #############
    @constraints(model, begin
        [i=1:sum(Mij)], xij[i] <= bigM*uij[i]
        [i=1:sum(Mjk)], xjk[i] <= bigM*ujk[i]
        [i=1:sum(Mkl)], xkl[i] <= bigM*ukl[i]
    end);
    ########### constraint 12 #############
    @constraints(model, begin
            [i in findnz(Vij2)[1]], xij[i] >= Vij2[i]*uij[i]
            [i in findnz(Vjk2)[1]], xjk[i] >= Vjk2[i]*ujk[i]
            # [i in findnz(Vkl)[1]], xkl[(i-1)+1:i] >= Vkl[i]*ukl[i]
    end);
    ########### constraint 13-14 #############
    @constraint(model, sum(y[1:J*2]) <= Jmax);
    @constraint(model, sum(y[J*2+1:end]) <= Kmax);
    return model
end


m2 = sm1dModel([1,150])
optimize!(m2); termination_status(m2)
get1d_objval(m2)

value.(m2[:y])
value.(m2[:h])
value.(m2[:uij])
value.(m2[:xij])
value.(m2[:ujk])
value.(m2[:xjk])
value.(m2[:ukl])
value.(m2[:xkl])

