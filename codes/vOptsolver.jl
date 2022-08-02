using DataFrames,DelimitedFiles,JuMP,LinearAlgebra,CPLEX,MathOptInterface,vOptGeneric
# using CPUTime
mutable struct Data2
    filepath::String; N::Dict{}; d::Array{}; c::Array{};  Mij::Array{}; Mjk::Array{}; Mkl::Array{};
    gij::Array{}; gjk::Array{}; gkl::Array{}; vij::Array{}; vjk::Array{}; vkl::Array{}; rij::Array{}; rjk::Array{}; rkl::Array{};
    Vij::Array{}; Vjk::Array{}; Vkl::Array{}; b::Array{}; upl::Int; udc::Int; bigM::Int # e::Array{};q::Array{};
    function Data2(filepath)
        dt = readdlm(filepath);
        # notafile = readdlm("/home/ak121396/Desktop/instances/SCND/Notations.txt", '=');
        # notafile = readdlm("F:/scnd/Notations.txt", '=');
        notafile = readdlm("/home/k2g00/k2g3475/scnd/Notations.txt", '=');
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

        gij = [];
        for i=1:N["supplier"]
            idx = 1; push!(gij,[]);
            for j=1:N["plant"]
                th = []
                for m=1:Mij[i,j]
                    push!(th, N["fixedcostModesp"][i][idx]);
                    idx+=1
                end
                push!(gij[i],th);
            end
        end
        gjk = [];
        for j=1:N["plant"]
            idx = 1; push!(gjk,[]);
            for k=1:N["distribution"]
                th = []
                for m=1:Mjk[j,k]
                    push!(th, N["fixedcostModepd"][j][idx]);
                    idx+=1
                end
                push!(gjk[j],th);
            end
        end
        gkl = [];
        for k=1:N["distribution"]
            idx = 1; push!(gkl,[]);
            for l=1:N["customer"]
                th= []
                for m=1:Mkl[k,l]
                    push!(th, N["fixedcostModedc"][k][idx]);
                    idx+=1
                end
                push!(gkl[k],th);
            end
        end
        vij = [];
        for i=1:N["supplier"]
            tmp = [];
            for j=1:N["plant"]
                tmp2 = []
                if j==1
                    for m=Mij[i,1]:-1:1
                        push!(tmp2, N["tcp"][i][5*(m-1)+1:5*(m-1)+5])
                    end
                else
                    for m=Mij[i,j]:-1:1
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
                    for m=Mjk[j,1]:-1:1
                        push!(tmp2, N["tcd"][j][5*(m-1)+1:5*(m-1)+5])
                    end
                else
                    for m=Mjk[j,k]:-1:1
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
                    for m=Mkl[k,1]:-1:1
                        push!(tmp2, N["tcc"][k][5*(m-1)+1:5*(m-1)+5])
                    end
                else
                    for m=Mkl[k,l]:-1:1
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
                        push!(tmp2, N["cep"][i][5*(m-1)+1:5*(m-1)+5])
                    end
                else
                    for m=1:Mij[i,j]
                        push!(tmp2, N["cep"][i][5*sum(Mij[i,1:j-1])+5*(m-1)+1:5*sum(Mij[i,1:j-1])+5*(m-1)+5]);
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
                        push!(tmp2, N["ced"][j][5*(m-1)+1:5*(m-1)+5])
                    end
                else
                    for m=1:Mjk[j,k]
                        push!(tmp2, N["ced"][j][5*sum(Mjk[j,1:k-1])+5*(m-1)+1:5*sum(Mjk[j,1:k-1])+5*(m-1)+5]);
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
                        push!(tmp2, N["cec"][k][5*(m-1)+1:5*(m-1)+5])
                    end
                else
                    for m=1:Mkl[k,l]
                        push!(tmp2, N["cec"][k][5*sum(Mkl[k,1:l-1])+5*(m-1)+1:5*sum(Mkl[k,1:l-1])+5*(m-1)+5]);
                    end
                end
                push!(tmp,tmp2);
            end
            push!(rkl,tmp);
        end

        Vij = [];
        for i=1:N["supplier"]
            idx = 1; push!(Vij,[]);
            for j=1:N["plant"]
                th = []
                for m=1:Mij[i,j]
                    push!(th, N["LcapacityModesp"][i][idx]);
                    idx+=1
                end
                push!(Vij[i],th);
            end
        end
        Vjk = [];
        for j=1:N["plant"]
            idx = 1; push!(Vjk,[]);
            for k=1:N["distribution"]
                th = []
                for m=1:Mjk[j,k]
                    push!(th, N["LcapacityModepd"][j][idx]);
                    idx+=1
                end
                push!(Vjk[j],th);
            end
        end
        Vkl = [];
        for k=1:N["distribution"]
            idx = 1; push!(Vkl,[]);
            for l=1:N["customer"]
                th= []
                for m=1:Mkl[k,l]
                    push!(th, N["LcapacityModedc"][k][idx]);
                    idx+=1
                end
                push!(Vkl[k],th);
            end
        end
        b = reshape(N["ves"],5,N["supplier"])';
        # q = append!(N["vep"],N["ved"]);
        upl = N["upperpants"]; udc = N["upperdistribution"]
        bigM = sum(sum(N["demand"]))
        new(filepath,N,d,c,Mij,Mjk,Mkl,gij,gjk,gkl,vij,vjk,vkl,rij,rjk,rkl,Vij,Vjk,Vkl,b,upl,udc,bigM); #cap,Mij,Mjk,Mkl,
    end
end
# file = "/home/ak121396/Desktop/instances/SCND/test01S2"
file = "/home/k2g00/k2g3475/scnd/instances/test01S2"
# file = "F:/model/Test1S2"
dt = Data2(file);
function buildmodel()
    model = vModel(CPLEX.Optimizer);
    set_silent(model)
    @variable(model, y[1:dt.N["plant"]+dt.N["distribution"],1:2], Bin);
    @variable(model, uij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:dt.Mij[i,j]], Bin);
    @variable(model, ujk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:dt.Mjk[j,k]], Bin);
    @variable(model, ukl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:dt.Mkl[k,l]], Bin);
    ############
    @variable(model, 0<= xij[i=1:dt.N["supplier"],j=1:dt.N["plant"],1:dt.Mij[i,j],1:5] );
    @variable(model, 0<= xjk[j=1:dt.N["plant"],k=1:dt.N["distribution"],1:dt.Mjk[j,k],1:5] );
    @variable(model, 0<= xkl[k=1:dt.N["distribution"],l=1:dt.N["customer"],1:dt.Mkl[k,l],1:5] );
    @variable(model, 0<= h[1:dt.N["plant"]+dt.N["distribution"],1:5,1:2] );
    @addobjective(model, Min,  sum(dt.c[j][t]*y[j,t] for j=1:1:dt.N["plant"]+dt.N["distribution"] for t=1:2) +
        sum(dt.gij[i][j][m]*uij[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j])+
        sum(dt.gjk[j][k][m]*ujk[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]) +
        sum(dt.gkl[k][l][m]*ukl[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])+
        sum(dt.N["vcs"][i][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) +
        sum(dt.vij[i][j][m][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
        sum(dt.vjk[j][k][m][p]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
        sum(dt.vkl[k][l][m][p]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5)+
        sum(dt.N["vcp"][j][2*(p-1)+t]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
        sum(dt.N["vcd"][k][2*(p-1)+t]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2));
    @addobjective(model, Min, sum(dt.b[i,p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+
        sum(dt.N["vep"][j][2*(p-1)+t]*h[j,p,t] for j=1:dt.N["plant"] for p=1:5 for t=1:2)+
        sum(dt.N["ved"][k][2*(p-1)+t]*h[k+dt.N["plant"],p,t] for k=1:dt.N["distribution"] for p=1:5 for t=1:2)+
        sum(dt.rij[i][j][m][p]*xij[i,j,m,p] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5)+sum(dt.rjk[j][k][m][p]*xjk[j,k,m,p] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k] for p=1:5)+
        sum(dt.rkl[k][l][m][p]*xkl[k,l,m,p] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l] for p=1:5));
    @constraint(model,[j=1:dt.N["plant"],p=1:5], sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j]) == sum(xjk[j,k,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]))
    @constraint(model,[k=1:dt.N["distribution"],p=1:5], sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k]) == sum(xkl[k,l,m,p] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l]))
    ########### constraint 4-6 #############
    @constraint(model,[j=1:dt.N["plant"],p=1:5], sum(h[j,p,:]) == sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j]))
    @constraint(model,[k=1:dt.N["distribution"],p=1:5], sum(h[k+dt.N["plant"],p,:] ) == sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k]))
    @constraint(model, [l=1:dt.N["customer"],p=1:5], sum(xkl[k,l,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mkl[k,l]) >= dt.d[l][p])
    ########### constraint 7-9 #############
    @constraint(model,[i=1:dt.N["supplier"]], sum(xij[i,j,m,p] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) <= dt.N["cas"][i])
    @constraint(model,[j=1:dt.N["plant"]+dt.N["distribution"], t=1:2], sum(h[j,:,t]) <=[dt.N["cap"];dt.N["cad"]][j]*y[j,t])
    @constraint(model,[j=1:dt.N["plant"]+dt.N["distribution"]], sum(y[j,:]) <= 1);
    ########### constraint 10 #############
    @constraint(model,[i=1:dt.N["supplier"],j=1:dt.N["plant"]], sum(uij[i,j,m] for m=1:dt.Mij[i,j]) <= 1);
    @constraint(model,[j=1:dt.N["plant"],k=1:dt.N["distribution"]], sum(ujk[j,k,m] for m=1:dt.Mjk[j,k]) <= 1);
    @constraint(model,[k=1:dt.N["distribution"],l=1:dt.N["customer"]], sum(ukl[k,l,m] for m=1:dt.Mkl[k,l]) <= 1);
    ########### constraint 11 #############
    BigM = sum(sum(dt.N["demand"]));
    @constraint(model,[i=1:dt.N["supplier"],j=1:dt.N["plant"],m=1:dt.Mij[i,j]], sum(xij[i,j,m,p] for p=1:5) <= BigM*uij[i,j,m]);
    @constraint(model,[j=1:dt.N["plant"],k=1:dt.N["distribution"],m=1:dt.Mjk[j,k]] ,sum(xjk[j,k,m,p] for p=1:5) <= BigM*ujk[j,k,m]);
    @constraint(model,[k=1:dt.N["distribution"],l=1:dt.N["customer"],m=1:dt.Mkl[k,l]], sum(xkl[k,l,m,p] for p=1:5) <= BigM*ukl[k,l,m]);
    ########### constraint 12 #############
    @constraint(model,[i=1:dt.N["supplier"], j=1:dt.N["plant"], m=1:dt.Mij[i,j]], sum(xij[i,j,m,p] for p=1:5) >= dt.Vij[i][j][m]*uij[i,j,m] );
    @constraint(model,[j=1:dt.N["plant"], k=1:dt.N["distribution"], m=1:dt.Mjk[j,k]], sum(xjk[j,k,m,p] for p=1:5) >= dt.Vjk[j][k][m]*ujk[j,k,m]);
    @constraint(model,[k=1:dt.N["distribution"], l=1:dt.N["customer"], m=1:dt.Mkl[k,l]], sum(xkl[k,l,m,p] for p=1:5) >= dt.Vkl[k][l][m]*ukl[k,l,m]);
    ########### constraint 13-14 #############
    @constraint(model,sum(y[j,t] for j=1:dt.N["plant"] for t=1:2) <= dt.upl);
    @constraint(model,sum(y[j,t] for j=dt.N["plant"]+1:dt.N["distribution"]+dt.N["plant"] for t=1:2) <= dt.udc);
    return model
end
model = buildmodel()
vSolve(model, method=:dicho)
# vSolve(model, method=:epsilon, step=1000000000.0, verbose=true);
ndpoints = getY_N(model)
function saveSol(fname::String, Y_N)#, elapsedTime::Float64)
    way = fname*"Y_N"
    open(way, "w") do f
        write(f, "$elapsedTime \n")
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

#################################################################################################
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

#################### Mathematical Model of SCND ###############
