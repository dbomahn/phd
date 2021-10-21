using DataFrames,DelimitedFiles,JuMP,LinearAlgebra,CPLEX,MathProgBase
const MPB = MathProgBase

function loadlp(filename,solver=CplexSolver(CPX_PARAM_SCRIND=0))
    model=buildlp([-1,0],[2 1],'<',1.5, solver) # create dummy model with correct solver
    MPB.loadproblem!(model,filename) # load what we actually want
    return model
end

mutable struct Data
    filepath::String; N::Dict{}; d::Array{}; c::Array{}; a::Array{}; e::Array{}; cap::Array{};
    gij::Array{}; gjk::Array{}; gkl::Array{}; vij::Array{}; vjk::Array{}; vkl::Array{};
    Vij::Array{}; Vjk::Array{}; Vkl::Array{}; b::Array{}; q::Array{};
    rij::Array{}; rjk::Array{}; rkl::Array{};
    function Data(filepath)
        dt = readdlm(filepath);
        notafile = readdlm("E:/scnd/Notations.txt", '=');
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
        d = N["demand"];
        c = append!(N["fcp"],N["fcd"]);
        a = N["vcs"];
        e = append!(N["vcp"],N["vcd"]);
        cap = append!(N["cas"],N["cap"],N["cad"]);
        Mij = transpose(reshape(N["ModeIJ"], (N["plant"],N["supplier"])));
        Mjk = transpose(reshape(N["ModeJK"], (N["distribution"],N["plant"])));
        Mkl = transpose(reshape(N["ModeKL"], (N["customer"],N["distribution"])));

        gij = [];
        for i=1:N["supplier"]
            idx = 1; push!(gij,[]);
            for j=1:N["plant"]
                fc = []
                for m=1:Mij[i,j]
                    push!(fc, N["fixedcostModesp"][i][idx]);
                    idx+=1
                end
                push!(gij[i],fc);
            end
        end
        gjk = [];
        for j=1:N["plant"]
            idx = 1; push!(gjk,[]);
            for k=1:N["distribution"]
                fc = []
                for m=1:Mjk[j,k]
                    push!(fc, N["fixedcostModepd"][j][idx]);
                    idx+=1
                end
                push!(gjk[j],fc);
            end
        end
        gkl = [];
        for k=1:N["distribution"]
            idx = 1; push!(gkl,[]);
            for l=1:N["customer"]
                fc = []
                for m=1:Mkl[k,l]
                    push!(fc, N["fixedcostModedc"][k][idx]);
                    idx+=1
                end
                push!(gkl[k],fc);
            end
        end
        vij = [];
        for i=1:N["supplier"]
            idx = 1; push!(vij,[]);
            for j=1:N["plant"]
                tc = []
                for m=1:Mij[i,j]
                    push!(tc, N["tcp"][i][idx:idx+4]);
                    idx+=5
                end
                push!(vij[i],tc);
            end
        end
        vjk = [];
        for j=1:N["plant"]
            idx = 1; push!(vjk,[]);
            for k=1:N["distribution"]
                tc = []
                for m=1:Mjk[j,k]
                    push!(tc, N["tcd"][j][idx:idx+4]);
                    idx+=5
                end
                push!(vjk[j],tc);
            end
        end
        vkl = [];
        for k=1:N["distribution"]
            idx = 1; push!(vkl,[]);
            for l=1:N["customer"]
                tc = []
                for m=1:Mkl[k,l]
                    push!(tc, N["tcc"][k][idx:idx+4]);
                    idx+=5
                end
                push!(vkl[k],tc);
            end
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
        b = reshape( N["ves"], (N["supplier"],Int(length(N["ves"])/N["supplier"])) );
        q = append!(N["vep"],N["ved"]);
        rij = [];
        for i=1:N["supplier"]
            idx = 1; push!(rij,[]);
            for j=1:N["plant"]
                em = []
                for m=1:Mij[i,j]
                    push!(em, N["cep"][i][idx:idx+4]);
                    idx+=5
                end
                push!(rij[i],em);
            end
        end
        rjk = []
        for j=1:N["plant"]
            idx = 1; push!(rjk,[]);
            for k=1:N["distribution"]
                em = []
                for m=1:Mjk[j,k]
                    push!(em, N["ced"][j][idx:idx+4]);
                    idx+=5
                end
                push!(rjk[j],em);
            end
        end
        rkl = []
        for k=1:N["distribution"]
            idx = 1; push!(rkl,[]);
            for l=1:N["customer"]
                em = []
                for m=1:Mkl[k,l]
                    push!(em, N["cec"][k][idx:idx+4]);
                    idx+=5
                end
                push!(rkl[k],em);
            end
        end

        new(filepath,N,d,c,a,e,cap,gij,gjk,gkl,vij,vjk,vkl,Vij,Vjk,Vkl,b,q,rij,rjk,rkl);
    end

end
# fpath = "/home/ak121396/Desktop/instances/SCND/"
# fpath = "E:/scnd/"
dt = Data(fpath*"Test1S1")
dt.N
arc = 1:N["supplier"]*N["plant"]+N["plant"]*N["distribution"]+N["distribution"]*N["customer"]

##########################  Mathematical model  #########################
scnd = Model()#CPLEX.Optimizer

# @variable(scnd, x[1:N["supplier"][1], 1:N["plant"][1],1:N["mode"][1]])
@variable(scnd, y[1:dt.N["plant"]+dt.N["distribution"],1:2] )
@variable(scnd, x[1:dt.N["supplier"],1:dt.N["plant"],1:2,1:5] )

@variable(scnd, h[1:N["plant"]+N["distribution"]] )
@variable(scnd, u[1:narc,] )

@constraint(scnd, sum(xij[i,m]) == sum(xjk[k,m]))
@constraint(scnd, dot(data.B[k,:],x) <= dot(xij[1:N[]]))
@constraint(scnd, dot(data.B[k,:],x) == data.RHS[k])

@objective(scnd, Min,
    sum(c[j][t]*y[j,t] for j=1:1:N["plant"]+N["distribution"] for t=1:2)+

)

optimize!(scnd);


sum(c[j][t]*y[j,t] for j=1:1:N["plant"]+N["distribution"] for t=1:2)

N["tcp"][1]

c[j][t]*y[j,t]
i=6
id1 = findall(x->x==nota[i], dt)[1][1]
id2 = findall(x->x==nota[i+1], dt)
[1][1]
if id2-id1 <3
