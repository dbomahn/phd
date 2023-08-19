# cd("C:/Users/AK121396/Desktop/ProjectBenders")
using DataFrames,DelimitedFiles,JuMP,LinearAlgebra,CPLEX,MathOptInterface,SparseArrays

struct Data1d
    file::String; N::Dict{}; d::Array{}; c::Array{}; a::Array{}; e::Array{}; gij::SparseVector{}; gjk::SparseVector{}; gkl::SparseVector{};
    vij::Array{}; vjk::Array{}; vkl::Array{}; Vij::SparseVector{}; Vjk::SparseVector{}; Mij::Array{}; Mjk::Array{}; Mkl::Array{};
    b::Array{}; q::Array{}; rij::Array{}; rjk::Array{}; rkl::Array{}; upl::Int; udc::Int; bigM::Int
    function Data1d(file)
        dt1 = readdlm(file);
        notafile = readdlm("/home/ak121396/Desktop/instances/scnd/Notations.txt", '=');
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

@show file = ARGS[1];
# file = "F:/scnd/Test4S4"
file = "/home/ak121396/Desktop/instances/scnd/test01S2"
# file = "/home/k2g00/k2g3475/scnd/instances/test01S2"
dt1 = Data1d(file);
function SCND1dim()
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
    @objective(scnd1, Min, sum(dt1.c.*y1) +
        sum(repeat(dt1.a[1,:], outer=sum(dt1.Mij[1,:])).*xij1[1:sum(dt1.Mij[1,:])*5])+
        sum(sum(repeat(dt1.a[i,:], outer=sum(dt1.Mij[i,:])).*xij1[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"])+
        sum(dt1.gij[i]*uij1[i] for i in findnz(dt1.gij)[1]) + sum(dt1.gjk[i]*ujk1[i] for i in findnz(dt1.gjk)[1]) + sum(dt1.gkl[i].*ukl1[i] for i in findnz(dt1.gkl)[1])+
        sum(dt1.e.*h1) + sum(dt1.vij.*xij1)+sum(dt1.vjk.*xjk1)+sum(dt1.vkl.*xkl1)
    )
    # @constraint(scnd1, obj2, sum(repeat(dt1.b[1,:], outer=sum(dt1.Mij[1,:])).*xij1[1:sum(dt1.Mij[1,:])*5]) +
    #         sum(sum(repeat(dt1.b[i,:], outer=sum(dt1.Mij[i,:])).*xij1[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"]) +
    #         sum(dt1.q.*h1) + sum(dt1.rij.*xij1)+sum(dt1.rjk.*xjk1)+sum(dt1.rkl.*xkl1) <= 1.1826921599716e6 )
    # @objective(scnd1, Min, sum(repeat(dt1.b[1,:], outer=sum(dt1.Mij[1,:])).*xij1[1:sum(dt1.Mij[1,:])*5]) +
    #         sum(sum(repeat(dt1.b[i,:], outer=sum(dt1.Mij[i,:])).*xij1[sum(dt1.Mij[1:i-1,:])*5+1:sum(dt1.Mij[1:i,:])*5]) for i=2:dt1.N["supplier"]) +
    #         sum(dt1.q.*h1) + sum(dt1.rij.*xij1)+sum(dt1.rjk.*xjk1)+sum(dt1.rkl.*xkl1));
    # w=[0.5,0.5]
    # w=[1,1]
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
    # # ########### constraint 10 #############
    @constraint(scnd1, sum(uij1[1:dt1.Mij[1,1]]) <= 1)
    @constraint(scnd1, sum(ujk1[1:dt1.Mjk[1,1]]) <= 1)
    @constraint(scnd1, sum(ukl1[1:dt1.Mkl[1,1]]) <= 1)
    @constraints(scnd1,begin
        sum(uij[1:dt.Mij[1,1]]) <= 1
        [j=2:dt1.N["plant"]], sum(uij1[sum(dt1.Mij[1,1:j-1])+1:sum(dt1.Mij[1,1:j-1])+dt1.Mij[1,j]]) <= 1
        [i=2:dt1.N["supplier"],j=2:dt1.N["plant"]],  sum(uij1[sum(dt1.Mij[1:i-1,:])+sum(dt1.Mij[i,1:j-1])+1:sum(dt1.Mij[1:i-1,:])+sum(dt1.Mij[i,1:j-1])+dt1.Mij[i,j]])<= 1
        sum(ujk[1:dt.Mjk[1,1]]) <= 1
        [k=2:dt1.N["distribution"]], sum(ujk1[sum(dt1.Mjk[1,1:k-1])+1:sum(dt1.Mjk[1,1:k-1])+dt1.Mjk[1,k]]) <= 1
        [j=2:dt1.N["plant"],k=2:dt1.N["distribution"]],  sum(ujk1[sum(dt1.Mjk[1:j-1,:])+sum(dt1.Mjk[j,1:j-1])+1:sum(dt1.Mjk[1:j-1,:])+sum(dt1.Mjk[j,1:j-1])+dt1.Mjk[j,k]]) <= 1
        sum(ukl[1:dt.Mkl[1,1]]) <= 1
        [l=2:dt1.N["customer"]], sum(ukl1[sum(dt1.Mkl[1,1:l-1])+1:sum(dt1.Mkl[1,1:l-1])+dt1.Mkl[1,l]]) <= 1
        [k=2:dt1.N["distribution"],l=2:dt1.N["customer"]],  sum(ukl1[sum(dt1.Mkl[1:k-1,:])+sum(dt1.Mkl[k,1:l-1])+1:sum(dt1.Mkl[1:k-1,:])+sum(dt1.Mkl[k,1:l-1])+dt1.Mkl[k,l]])<= 1
    end);
    ########### constraint 11 ############# This causes a different obj values
    @constraint(scnd1, [i=1:sum(dt1.Mij)], sum(xij1[5*(i-1)+1:5*i]) <= dt1.bigM*uij1[i])
    @constraint(scnd1, [j=1:sum(dt1.Mjk)], sum(xjk1[5*(j-1)+1:5*j]) <= dt1.bigM*ujk1[j])
    # @constraint(scnd1, [k=1:sum(dt1.Mkl)], sum(xkl1[5*(k-1)+1:5*k]) <= dt1.bigM*ukl1[k])
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

@constraints(scnd1,[i in findnz(dt1.Vij)[1]], sum(xij1[5*(i-1)+1:5*i]) >= dt1.Vij[i]*uij1[i])
scnd1 = SCND1dim()
#optimize!(scnd1); objective_value(scnd1)
#solve_time(scnd1)
# termination_status(scnd1)
#write_to_file(scnd1, "/home/ak121396/Desktop/instances/SCND/small/test1s2_obj1.lp")

###########################    Make vlp file   #########################
using MathProgBase
const MPB = MathProgBase
function loadlp(filename,solver=CplexSolver(CPX_PARAM_SCRIND=0))
    model=buildlp([-1,0],[2 1],'<',1.5, solver) # create dummy model with correct solver
    MPB.loadproblem!(model,filename) # load what we actually want
    return model
end
write_to_file(scnd1, file*"1dim.lp")
lpmodel = loadlp(file*"1dim.lp") 
Bmtx = MPB.getconstrmatrix(lpmodel);
# cut = findall(i-> varub[i]==1 && varub[i+1]!=1, 1:length(varub))[end]
# B = Bmtx[3:end,1:cut];P = Bmtx[1:2,1:cut]; vub = varub[1:cut]
B = Bmtx[3:end,:];P = Bmtx[1:2,:]; vub = MPB.getvarUB(lpmodel)
m,n=size(B)
lb = MPB.getconstrLB(lpmodel)[3:end]
ub = MPB.getconstrUB(lpmodel)[3:end]
RHS = []
for i=1:m
    if ub[i]==Inf
        push!(RHS,lb[i])
    else
        push!(RHS,ub[i])
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
########################## Make vlp file for Bensolve #######################
nz = count(i->(i!=0),B)
objnz = count(i->(i!=0),P)
obj=size(P)[1]
wholearray=[];
arr=["p vlp min",m,n,nz,obj,objnz]
push!(wholearray,arr)

for i=1:m
    for j=1:n
        if (B[i,j]!=0)
            if (B[i,j]%1) == 0 #if B[i,j] is Int
                push!(wholearray,("a",i,j,Int128(B[i,j])))
            else# B[i,j] is Float
                push!(wholearray,("a",i,j,Float64(B[i,j])))
            end
        end
    end
end
for i=1:obj
    for j=1:n
        if P[i,j]!=0
            push!(wholearray,("o",i,j,P[i,j]))
        end
    end
end
for i=1:m
    push!(wholearray,("i",i,signs[i],RHS[i]))
end
# for j=1:n
#     if j in bvar
#         push!(wholearray,("j",j,"s",fpx[1][j])) #assign FP int var values
#     else
#         push!(wholearray,("j", j,'l',0))
#     end
# end
for j=1:n
    if vub[j]==1
        push!(wholearray,("j",j,"d",0,1))
    else
        push!(wholearray,("j", j,'l',0))
    end
end
push!(wholearray,"e")

ins = open("/home/k2g00/k2g3475/scnd/vlp/"*file[36:end]*"1d.vlp","w")
#ins = open(file*".vlp","w")
writedlm(ins,wholearray)
close(ins)
