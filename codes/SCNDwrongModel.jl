using DataFrames,DelimitedFiles,JuMP,LinearAlgebra,CPLEX

struct AuthorData
    filepath::String; N::Dict{}; d::Array{}; c::Array{};  Mij::Array{}; Mjk::Array{}; Mkl::Array{};
    vij::Array{}; vjk::Array{}; vkl::Array{}; rij::Array{}; rjk::Array{}; rkl::Array{};
    Vij::Int; Vjk::Int; b::Array{}; upl::Int; udc::Int; bigM::Int # e::Array{};q::Array{};
    function AuthorData(filepath)
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

        Vij = findmax(findmax(N["LcapacityModesp"])[1])[1]
        Vjk = findmax(findmax(N["LcapacityModepd"])[1])[1]

        b = reshape(N["ves"],N["supplier"],5);
        # q = append!(N["vep"],N["ved"]);
        upl = N["upperpants"]; udc = N["upperdistribution"]
        bigM = sum(sum(N["demand"]))
        new(filepath,N,d,c,Mij,Mjk,Mkl,vij,vjk,vkl,rij,rjk,rkl,Vij,Vjk,b,upl,udc,bigM); #cap,Mij,Mjk,Mkl,
    end
end
file = "/home/ak121396/Desktop/instances/SCND/test01S2"
file = "/home/ak121396/Desktop/instances/SCND/Test1S1"
dt = AuthorData(file);
# @show file = ARGS[1]


function Authormodel()
    scnd = Model(CPLEX.Optimizer); set_silent(scnd);
    @variable(scnd, y[1:dt.N["plant"]+dt.N["distribution"],1:2], Bin);
    @variable(scnd, uij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:2], Bin);
    @variable(scnd, ujk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:2], Bin);
    @variable(scnd, ukl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:2], Bin);
    ############
    @variable(scnd, 0<= xij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:2,1:5] );
    @variable(scnd, 0<= xjk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:2,1:5] );
    @variable(scnd, 0<= xkl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:2,1:5] );
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
    # 1st obj
    # @constraint(scnd, obj1,
    #     sum(dt.c[j][t]*y[j,t] for j=1:1:dt.N["plant"]+dt.N["distribution"] for t=1:2) +
    #     sum(dt.gij[i][j][m]*uij[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j])+
    #     sum(dt.gjk[j][k][m]*ujk[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]) +
    #     sum(dt.gkl[k][l][m]*ukl[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])+
    #     sum(dt.N["vcs"][i][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) +
    #     sum(dt.vij[i][j][m][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
    #     sum(dt.vjk[j][k][m][p]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
    #     sum(dt.vkl[k][l][m][p]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5)+
    #     sum(dt.N["vcp"][j][5*(t-1)+p]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
    #     sum(dt.N["vcd"][k][5*(t-1)+p]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2) <=0);
    # #2nd obj
    # @constraint(scnd, obj2, sum(dt.b[i,p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
    #     sum(dt.N["vep"][j][2*(p-1)+t]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
    #     sum(dt.N["ved"][k][2*(p-1)+t]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2)+
    #     sum(dt.rij[i][j][m][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+sum(dt.rjk[j][k][m][p]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
    #     sum(dt.rkl[k][l][m][p]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5)<=0);

    @objective(scnd, Min, sum(dt.c[j][t]*y[j,t] for j=1:1:dt.N["plant"]+dt.N["distribution"] for t=1:2) +exg+
        sum(dt.gij[i][j][m]*uij[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j])+
        sum(dt.gjk[j][k][m]*ujk[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]) +
        sum(dt.gkl[k][l][m]*ukl[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])+
        sum(dt.N["vcs"][i][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:2 for p=1:5) +
        sum(dt.vij[i][j][m][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
        sum(dt.vjk[j][k][m][p]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
        sum(dt.vkl[k][l][m][p]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5)+
        sum(dt.N["vcp"][j][2*(p-1)+t]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
        sum(dt.N["vcd"][k][2*(p-1)+t]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2));
    # @objective(scnd,Min,sum(dt.b[i,p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:2 for p=1:5)+
    #     sum(dt.N["vep"][j][2*(p-1)+t]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
    #     sum(dt.N["ved"][k][2*(p-1)+t]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2)+
    #     sum(dt.rij[i][j][m]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
    #     sum(dt.rjk[j][k][m]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
    #     sum(dt.rkl[k][l][m]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5));
    @constraint(scnd,[j=1:dt.N["plant"],p=1:5], sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:2) == sum(xjk[j,k,m,p] for k=1:dt.N["distribution"] for m=1:2));
    @constraint(scnd,[k=1:dt.N["distribution"],p=1:5], sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:2) == sum(xkl[k,l,m,p] for l=1:dt.N["customer"] for m=1:2));
    @constraint(scnd,[j=1:dt.N["plant"],p=1:5], sum(h[j,p,:]) == sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:2));
    @constraint(scnd,[k=1:dt.N["distribution"],p=1:5], sum(h[k+dt.N["plant"],p,:] ) == sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:2));
    @constraint(scnd, [l=1:dt.N["customer"],p=1:5], sum(xkl[k,l,m,p] for k=1:dt.N["distribution"] for m=1:2) >= dt.d[l][p])
    ########### constraint 7-9 #############
    @constraint(scnd,[i=1:dt.N["supplier"]], sum(xij[i,j,m,p] for j=1:dt.N["plant"] for m=1:2 for p=1:5) <= dt.N["cas"][i]);
    @constraint(scnd,[j=1:dt.N["plant"]+dt.N["distribution"], t=1:2], sum(h[j,:,t]) <=[dt.N["cap"];dt.N["cad"]][j]*y[j,t])
    @constraint(scnd,[j=1:dt.N["plant"]+dt.N["distribution"]], sum(y[j,:]) <= 1)
    ########### constraint 10 #############
    @constraint(scnd,[i=1:dt.N["supplier"],j=1:dt.N["plant"]], sum(uij[i,j,:]) <= 1)
    @constraint(scnd,[j=1:dt.N["plant"],k=1:dt.N["distribution"]], sum(ujk[j,k,:] ) <= 1)
    @constraint(scnd,[k=1:dt.N["distribution"],l=1:dt.N["customer"]], sum(ukl[k,l,:] ) <= 1)
    ########### constraint 12 #############
    @constraint(scnd,[i=1:dt.N["supplier"], j=1:dt.N["plant"]], sum(xij[i,j,2,p] for p=1:5) >= dt.Vij*uij[i,j,2] );
    @constraint(scnd,[j=1:dt.N["plant"], k=1:dt.N["distribution"]], sum(xjk[j,k,2,p] for p=1:5) >= dt.Vjk*ujk[j,k,2]);
    # ########## constraint 13-14 #############
    @constraint(scnd,sum(y[j,t] for j=1:dt.N["plant"] for t=1:2) <= dt.upl)
    @constraint(scnd,sum(y[j,t] for j=dt.N["plant"]+1:dt.N["distribution"]+dt.N["plant"] for t=1:2) <= dt.udc)
    return scnd
end

amodel = Authormodel()
optimize!(amodel); objective_value(amodel)
sum(value.(amodel[:xkl]))
1
# write_to_file(scnd, "/home/ak121396/Desktop/relise/Mytest1S1.lp")
###################
write_to_file(amodel, file*".lp")
using MathProgBase
const MPB = MathProgBase
function loadlp(filename,solver=CplexSolver(CPX_PARAM_SCRIND=0))
    model=buildlp([-1,0],[2 1],'<',1.5, solver) # create dummy model with correct solver
    MPB.loadproblem!(model,filename) # load what we actually want
    return model
end

lpmodel = loadlp(file*".lp")
Bmtx = MPB.getconstrmatrix(lpmodel);
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

for j=1:n
    if vub[j]==1
        push!(wholearray,("j",j,"d",0,1))
    else
        push!(wholearray,("j", j,'l',0))
    end
end
push!(wholearray,"e")

ins = open("/home/k2g00/k2g3475/scnd/vlp/"*file[36:end]*".vlp","w")
writedlm(ins,wholearray)
close(ins)
