cd("C:/Users/AK121396/Desktop/ProjectBenders")
using DataFrames,DelimitedFiles,JuMP,LinearAlgebra,CPLEX,MathOptInterface,SparseArrays#,CPUTime

mutable struct Data1
    file::String; N::Dict{}; d::Array{}; c::Array{}; a::Array{}; e::Array{}; gij::SparseVector{}; gjk::SparseVector{}; gkl::SparseVector{};
    vij::Array{}; vjk::Array{}; vkl::Array{}; Vij::SparseVector{}; Vjk::SparseVector{}; Vkl::SparseVector{}; Mij::Array{}; Mjk::Array{}; Mkl::Array{};
    b::Array{}; q::Array{}; rij::Array{}; rjk::Array{}; rkl::Array{}; upl::Int; udc::Int; bigM::Int
    function Data1(file)
        dt1 = readdlm(file);
        # notafile = readdlm("/home/ak121396/Desktop/instances/SCND/Notations.txt", '=');
        notafile = readdlm("F:/scnd/Notations.txt", '=');
        # notafile = readdlm("/home/k2g00/k2g3475/scnd1/Notations.txt", '=');
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
        Vij = sparse(N["LcapacityModesp"]); Vjk = sparse(N["LcapacityModepd"]); Vkl =  sparse(N["LcapacityModedc"]);
        Mij = transpose(reshape(N["ModeIJ"], (N["plant"],N["supplier"])));
        Mjk = transpose(reshape(N["ModeJK"], (N["distribution"],N["plant"])));
        Mkl = transpose(reshape(N["ModeKL"], (N["customer"],N["distribution"])));
        b = reshape(N["ves"],5,N["supplier"])';  q = append!(N["vep"],N["ved"]);
        rij = N["cep"]; rjk = N["ced"]; rkl = N["cec"];
        upl = N["upperpants"]; udc = N["upperdistribution"]; bigM = sum(N["demand"])

        new(file,N,d,c,a,e,gij,gjk,gkl,vij,vjk,vkl,Vij,Vjk,Vkl,Mij,Mjk,Mkl,b,q,rij,rjk,rkl,upl,udc,bigM);
    end
end

# @show file = ARGS[1];
file = "F:/scnd/Test1S2"
# file = "/home/ak121396/Desktop/instances/SCND/test01S2"
# file = "/home/k2g00/k2g3475/scnd1/instances/test01S2"
dt1 = Data1(file);
##########################  Mathematical model  #########################
scnd1 = Model(CPLEX.Optimizer); set_silent(scnd1)
# MOI.set(scnd1, MOI.NumberOfThreads(), 1);
# @variable(scnd1, 0<= y1[1:(dt1.N["plant"]+dt1.N["distribution"])*2] <= 1);
# @variable(scnd1, 0<= uij1[1:sum(dt1.Mij)] <= 1);
# @variable(scnd1, 0<= ujk1[1:sum(dt1.Mjk)] <= 1);
# @variable(scnd1, 0<= ukl1[1:sum(dt1.Mkl)] <= 1);
#########################  IP  ########################################
@variable(scnd1, y1[1:(dt1.N["plant"]+dt1.N["distribution"])*2], Bin)
@variable(scnd1, uij1[1:sum(dt1.Mij)], Bin)
@variable(scnd1, ujk1[1:sum(dt1.Mjk)], Bin)
@variable(scnd1, ukl1[1:sum(dt1.Mkl)], Bin)

@variable( scnd1, 0<= xij1[1:sum(dt1.Mij)*5] )
@variable( scnd1, 0<= xjk1[1:sum(dt1.Mjk)*5] )
@variable( scnd1, 0<= xkl1[1:sum(dt1.Mkl)*5] )
@variable( scnd1, 0<= h1[1:(dt1.N["plant"]+dt1.N["distribution"])*5*2] )


# @constraint(scnd1, obj1, sum(dt1.c.*y1) +sum(repeat(dt1.a[1,:], outer=sum(dt1.Mij[1,:])).*xij1[1:sum(dt1.Mij[1,:])*5])+
#         sum(sum(repeat(dt1.a[i,:], outer=sum(dt1.Mij[i,:])).*xij1[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"])+
#         sum(dt1.e.*h1) + sum(dt1.gij[i]*uij1[i] for i in findnz(dt1.gij)[1]) + sum(dt1.gjk[i]*ujk1[i] for i in findnz(dt1.gjk)[1]) + sum(dt1.gkl[i].*ukl1[i] for i in findnz(dt1.gkl)[1])+
#         sum(dt1.vij.*xij1)+sum(dt1.vjk.*xjk1)+sum(dt1.vkl.*xkl1) <= 0)
@objective(scnd1, Min, sum(dt1.c.*y1) +sum(repeat(dt1.a[1,:], outer=sum(dt1.Mij[1,:])).*xij1[1:sum(dt1.Mij[1,:])*5])+
        sum(sum(repeat(dt1.a[i,:], outer=sum(dt1.Mij[i,:])).*xij1[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"])+
        sum(dt1.e.*h1) + sum(dt1.gij[i]*uij1[i] for i in findnz(dt1.gij)[1]) + sum(dt1.gjk[i]*ujk1[i] for i in findnz(dt1.gjk)[1]) + sum(dt1.gkl[i].*ukl1[i] for i in findnz(dt1.gkl)[1])+
        sum(dt1.vij.*xij1)+sum(dt1.vjk.*xjk1)+sum(dt1.vkl.*xkl1));
# @constraint(scnd1, obj2, sum(repeat(dt1.b[1,:], outer=sum(dt1.Mij[1,:])).*xij1[1:sum(dt1.Mij[1,:])*5]) +
#         sum(sum(repeat(dt1.b[i,:], outer=sum(dt1.Mij[i,:])).*xij1[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"]) +
#         sum(dt1.q.*h1) + sum(dt1.rij.*xij1)+sum(dt1.rjk.*xjk1)+sum(dt1.rkl.*xkl1) <= 0)
@objective(scnd1, Min, sum(repeat(dt1.b[1,:], outer=sum(dt1.Mij[1,:])).*xij1[1:sum(dt1.Mij[1,:])*5]) +
        sum(sum(repeat(dt1.b[i,:], outer=sum(dt1.Mij[i,:])).*xij1[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"]) +
        sum(dt1.q.*h1) + sum(dt1.rij.*xij1)+sum(dt1.rjk.*xjk1)+sum(dt1.rkl.*xkl1));

########## constraint 3 #############
# @constraints(scnd1, begin
@constraint(scnd1, [p=1:5], sum(xij1[5*(m-1)+p] for m=1:dt1.Mij[1,1])+sum(xij1[5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,1]) == sum(xjk1[5*(m-1)+p] for m=1:sum(dt1.Mjk[1,:])) )
        # [j=1,p=1:5], sum(xij1[5*(j-1)+5*(m-1)+p] for m=1:dt1.Mij[1,j])+sum(xij1[5*(j-1)+5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,j]) == sum(xjk1[5*(m-1)+p] for m=1:sum(dt1.Mjk[j,:]))
@constraint(scnd1, [j=2:dt1.N["plant"],p=1:5], sum(xij1[5*(j-1)+5*(m-1)+p] for m=1:dt1.Mij[1,j])+sum(xij1[5*(j-1)+5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,j]) == sum(xjk1[sum(dt1.Mjk[1:j-1,:])*5 + 5*(m-1)+p] for m=1:sum(dt1.Mjk[j,:])) )
@constraint(scnd1, [p=1:5], sum(xjk1[5*(m-1)+p] for m=1:dt1.Mjk[1,1])+sum(xjk1[5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,1]) == sum(xkl1[5*(m-1)+p] for m=1:sum(dt1.Mkl[1,:])) )
@constraint(scnd1, [k=2:dt1.N["distribution"],p=1:5],sum(xjk1[5*(k-1)+5*(m-1)+p] for m=1:dt1.Mjk[1,k])+sum(xjk1[5*(k-1)+5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,k]) == sum(xkl1[sum(dt1.Mkl[1:k-1,:])*5 + 5*(m-1)+p] for m=1:sum(dt1.Mkl[k,:])) )
# end)

########### constraint 4-6 #############
# @constraints(scnd1, begin
@constraint(scnd1, [p=1:5],sum(h1[5*(t-1)+p] for t=1:2) == sum(xij1[5*(m-1)+p] for m=1:dt1.Mij[1,1])+sum(xij1[5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,1]))
@constraint(scnd1, [j=2:dt1.N["plant"],p=1:5], sum(h1[5*2*(j-1)+5*(t-1)+p] for t=1:2) == sum(xij1[5*(j-1)+5*(m-1)+p] for m=1:dt1.Mij[1,j])+sum(xij1[5*(j-1)+5*(m-1)+p+(5*sum(dt1.Mij[1:i-1,:]))] for i=2:dt1.N["supplier"] for m=1:dt1.Mij[i,j]))
@constraint(scnd1, [p=1:5], sum(h1[5*2*dt1.N["plant"]+5*(t-1)+p] for t=1:2) == sum(xjk1[5*(m-1)+p] for m=1:dt1.Mjk[1,1])+sum(xjk1[5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,1]));
@constraint(scnd1, [k=2:dt1.N["distribution"],p=1:5], sum(h1[5*2*dt1.N["plant"]+5*2*(k-1)+5*(t-1)+p] for t=1:2) == sum(xjk1[5*(k-1)+5*(m-1)+p] for m=1:dt1.Mjk[1,k])+sum(xjk1[5*(k-1)+5*(m-1)+p+(5*sum(dt1.Mjk[1:j-1,:]))] for j=2:dt1.N["plant"] for m=1:dt1.Mjk[j,k]))


@constraint(scnd1, [p=1:5], sum(xkl1[5*(m-1)+p+ sum(dt1.Mkl[1,1:l-1])*5] for m=1:dt1.Mkl[1,1]) +sum(xkl1[5*(m-1)+p+(5*sum(dt1.Mkl[1:k-1,:]))] for k=2:dt1.N["distribution"] for m=1:dt1.Mkl[k,1]) >= dt1.d[1,p])
@constraint(scnd1, [l=2:dt1.N["customer"], p=1:5], sum(xkl1[sum(dt1.Mkl[1,1:l-1])*5 + 5*(m-1)+p] for m=1:dt1.Mkl[1,l])+ sum(xkl1[5*sum(dt1.Mkl[1:k-1,:])+5*sum(dt1.Mkl[k,1:l-1])+5*(m-1)+p] for k=2:dt1.N["distribution"] for m=1:dt1.Mkl[k,l]) >= dt1.d[l,p])

########### constraint 7 #############
@constraints(scnd1, begin
        sum(xij1[1:5*sum(dt1.Mij[1,:])]) <= dt1.N["cas"][1]
        [i=2:dt1.N["supplier"]],  sum(xij1[5*sum(dt1.Mij[i-1,:])+1:5*sum(dt1.Mij[i,:])]) <= dt1.N["cas"][i]
end);
########### constraint 8 #############
@constraint(scnd1,[j=1:dt1.N["plant"]+dt1.N["distribution"], t=1:2], sum(h1[5*2*(j-1)+((t-1)*5)+1:5*2*(j-1)+((t-1)*5)+5]) <= [dt1.N["cad"];dt1.N["cas"]][j]*y1[(dt1.N["plant"]+dt1.N["distribution"])*(t-1)+j]);
########### constraint 9 #############
@constraint(scnd1,[j=1:dt1.N["plant"]+dt1.N["distribution"]], sum(y1[2*(j-1)+1:2*(j-1)+2]) <= 1);
########### constraint 10 #############
@constraints(scnd1,begin
        sum(uij1[1:dt1.Mij[1,1]]) <= 1
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
end)
########### constraint 12 #############
@constraints(scnd1, begin
        [i in findnz(dt1.Vij)[1]], sum(xij1[5*(i-1)+1:5*i]) >= dt1.Vij[i]*uij1[i]
        [i in findnz(dt1.Vjk)[1]], sum(xjk1[5*(i-1)+1:5*i]) >= dt1.Vjk[i]*ujk1[i]
        [i in findnz(dt1.Vkl)[1]], sum(xkl1[5*(i-1)+1:5*i]) >= dt1.Vkl[i]*ukl1[i]
end);
########### constraint 13-14 #############
@constraint(scnd1, sum(y1[1:dt1.N["plant"]]) <= dt1.upl);
@constraint(scnd1, sum(y1[dt1.N["plant"]+1:end]) <= dt1.udc);
scnd1
optimize!(scnd1)
objective_value(scnd1)
