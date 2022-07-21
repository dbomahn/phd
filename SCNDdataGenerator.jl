using StatsBase,DelimitedFiles,Distributions,CPLEX,JuMP,LinearAlgebra,CPUTime
# "Region: divide the whole grid by 5, one region is 40*40. total 25 regions"
# "Cluster: a dense set of nodes located in the same region"
# "Choose 2(originally 4 or 5) regions (at random) and randomly assign 60% nodes there. 40% are randomly scattered"
# "Max 1 plant and 2 DCs (originally 4 plants 8 DCs) can lie in the same cluster.Candidate DCs are located in the regions where customers populated"
# "Demand ~ [100, 300]. D = ∑∑d^p_{l}"
# "Capacity of DCs~[1.1,1.5], plants~[1.1,2], supplieres~[1.1,3]"
"The fixed cost of transportation mode 1 is 5% of the total transportation costs. The variable cost of mode 2 is 20% higher than that of mode 1"
# "The whole grid (200*200) divided into 5*5 area to assign t fixed cost of facilities; very expensive areas 5%, expensive areas 50%, intermediate 40%, cheap 5%
# cost_j = 𝚽*sqrt(capacity_j). cheap~[5000,20000], intermediate[20000,35000], expensive~[35000,50000], very expensive~[50000,60000]"
# "Processing cost at suppliers~[130,150], at plants~[130,150], at DCs~[100,120], then multiply by a random noising factor [0.9,1.2]"
# "Transportation cost: Euclidean distance between the arcs was multiplied by a unitary cost~[0.8,1.2] and by
# τ at supplier-plant~[1,1.3], at plant-Dc~[1.2,1.4], at DC-customer~[1.3,1.5]"

mutable struct Data1
    filepath::String; N::Dict{}; d::Array{}; m::Int; c::Array{}; e::Array{}; gij::Array{}; gjk::Array{}; gkl::Array{};
    Mij::Array{}; Mjk::Array{}; Mkl::Array{}; Vij::Array{}; Vjk::Array{}; Vkl::Array{}; b::Array{}; q::Array{};
    upl::Int; udc::Int
    # rij::Array{}; rjk::Array{}; rkl::Array{}; vij::Array{}; vjk::Array{}; vkl::Array{};
    function Data1(filepath)
        dt = readdlm(filepath);
        # notafile = readdlm("/home/ak121396/Desktop/instances/SCND/Notations.txt", '=');
        notafile = readdlm("F:/scnd/Notations.txt", '=');
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
        d = N["demand"];  m = N["transportation"];
        c = append!(N["fcp"],N["fcd"]); e = append!(N["vcp"],N["vcd"]);
        gij = N["fixedcostModesp"]; gjk = N["fixedcostModepd"]; gkl = N["fixedcostModedc"];
        # gij = replace.(N["fixedcostModesp"], 0=>10^(-3));gjk = replace.(N["fixedcostModepd"], 0=>10^(-3)); gkl = replace.(N["fixedcostModedc"], 0=>10^(-3));
        Mij = transpose(reshape(N["ModeIJ"], (N["plant"],N["supplier"])));
        Mjk = transpose(reshape(N["ModeJK"], (N["distribution"],N["plant"])));
        Mkl = transpose(reshape(N["ModeKL"], (N["customer"],N["distribution"])));

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
        upl = N["upperpants"]; udc = N["upperdistribution"]

        new(filepath,N ,d,m,c,e,gij,gjk,gkl,Mij,Mjk,Mkl,Vij,Vjk,Vkl,b,q,upl,udc); #cap,Mij,Mjk,Mkl,
    end
end
# file = "/home/ak121396/Desktop/instances/SCND/test01S2"
file = "F:/scnd/Test6S3"
dt = Data1(file);
rgs = rand(4:5,1)[1]
xcoordi = sample(1:5, 2); ycoordi = sample(1:5, 2);
cluij = [xcoordi ycoordi];``
xcoordi = sample(1:5, 2); ycoordi = sample(1:5, 2);
clukl = [xcoordi ycoordi];
I = Int(ceil(dt.N["supplier"]/3))
if I<=3
    J=I; K=2*I; L=5*I; Jmax = 2; Kmax = 2
else
    J=I; K=2*I; L=5*I; Jmax = ceil(J/2); Kmax =I;
end
cluijkl = round.(Int,0.6*[I,J,K,L]);
Mij = dt.Mij[1:I,1:J]; Mjk = dt.Mjk[1:J,1:K]; Mkl = dt.Mkl[1:K,1:L];


demand = sample(100:300,L);
cas = round.(Int,rand(Uniform(1.1,3), I).*(sum(demand)/Kmax));
cap = round.(Int,rand(Uniform(1.1,2), J)*(sum(demand)/Kmax) );
cad = round.(Int,rand(Uniform(1.1,1.5), K)*(sum(demand)/Kmax) );
capd = [cap;cad];
function nodecoordi(cluijkl,cluij,clukl,I,J,K,L)
    interval = 40; Ipt,Jpt,Kpt,Lpt = [],[],[],[]
    for i=1:I
        if i<=cluijkl[1]
            c = rand(1:size(cluij)[1])
            x,y1 = cluij[c,:]
            nx = rand(interval*(x-1)+1:interval*x)
            ny = rand(interval*(y1-1)+1:interval*y1)
        else
            nx,ny = rand(1:200,2)
        end
        push!(Ipt,(nx,ny))
    end
    for i=1:J
        if i<=cluijkl[2]
            c = rand(1:size(cluij)[1])
            x,y1 = cluij[c,:]
            nx = rand(interval*(x-1)+1:interval*x)
            ny = rand(interval*(y1-1)+1:interval*y1)
        else
            nx,ny = rand(1:200,2)
        end
        push!(Jpt,(nx,ny))
    end
    for i=1:K
        if i<=cluijkl[3]
            c = rand(1:size(cluij)[1])
            x,y1 = cluij[c,:]
            nx = rand(interval*(x-1)+1:interval*x)
            ny = rand(interval*(y1-1)+1:interval*y1)
        else
            nx,ny = rand(1:200,2)
        end
        push!(Kpt,(nx,ny))
    end
    for i=1:L
        if i<=cluijkl[4]
            c = rand(1:size(cluij)[1])
            x,y1 = cluij[c,:]
            nx = rand(interval*(x-1)+1:interval*x)
            ny = rand(interval*(y1-1)+1:interval*y1)
        else
            nx,ny = rand(1:200,2)
        end
        push!(Lpt,(nx,ny))
    end
    return Ipt,Jpt,Kpt,Lpt
end
Ipt,Jpt,Kpt,Lpt = nodecoordi(cluijkl,cluij,clukl,I,J,K,L);

function dist(x,y1)
    return sqrt(abs(x[1]-y1[1])^2 + abs(x[2]-y1[2])^2)
end
function fac_fixedcost(cluijkl,cap,cad,J,K)
    fc = []
    for i=1:J
        if i<=cluijkl[2]
            c1 = round(sqrt(cap[i])*rand(Uniform(35000,50000),1)[1])
            c2 = round(c1*1.2)
        else
            c1 = round(sqrt(cap[i])*rand(Uniform(20000,35000),1)[1])
            c2 = round(c1*1.2)
        end
        push!(fc,[c1,c2])
    end
    for i=1:K
        if i<=cluijkl[2]
            c1 = round(sqrt(cad[i])*rand(Uniform(35000,50000),1)[1])
            c2 = round(c1*1.2)
        else
            c1 = round(sqrt(cad[i])*rand(Uniform(35000,50000),1)[1])
            c2 = round(c1*1.2)
        end
        push!(fc,[c1,c2])
    end
    return fc
end
fc = fac_fixedcost(cluijkl,cap,cad,J,K);

function trans_cost(cluijkl,Ipt,Jpt,Kpt,Lpt,Mij,Mjk,Mkl)
    tcp,tcd,tcc = [],[],[]
    for i=1:length(Ipt)
        tmp = []
        for j=1:length(Jpt)
            d = dist(Ipt[i],Jpt[j])
            if Mij[i,j]==1
                push!(tmp, [round(d*rand(Uniform(0.8,1.2),1)[1]*rand(Uniform(1,1.3),1)[1];digits=4)])
            else
                push!(tmp, [round(d*rand(Uniform(0.8,1.2),1)[1]*rand(Uniform(1,1.3),1)[1];digits=4), round(d*rand(Uniform(0.8*1.2,1.2*1.2),1)[1]*rand(Uniform(1,1.3),1)[1]; digits=4)])
            end
        end
        push!(tcp,tmp)
    end
    for j=1:length(Jpt)
        tmp = []
        for k=1:length(Kpt)
            d = dist(Jpt[j],Kpt[k])
            if Mjk[j,k]==1
                push!(tmp, [round(d*rand(Uniform(0.8,1.2),1)[1]*rand(Uniform(1.2,1.4),1)[1];digits=4)])
            else
                push!(tmp, [round(d*rand(Uniform(0.8,1.2),1)[1]*rand(Uniform(1.2,1.4),1)[1];digits=4), round(d*rand(Uniform(0.8*1.2,1.2*1.2),1)[1]*rand(Uniform(1.2,1.4),1)[1];digits=4)])
            end
        end
        push!(tcd,tmp)
    end
    for k=1:length(Kpt)
        tmp = []
        for l=1:length(Lpt)
            d = dist(Kpt[k],Lpt[l])
            if Mkl[k,l]==1
                push!(tmp, [round(d*rand(Uniform(0.8,1.2),1)[1]*rand(Uniform(1.3,1.5),1)[1];digits=4)])
            else
                push!(tmp, [round(d*rand(Uniform(0.8,1.2),1)[1]*rand(Uniform(1.3,1.5),1)[1];digits=4), round(d*rand(Uniform(0.8*1.2,1.2*1.2),1)[1]*rand(Uniform(1.3,1.5),1)[1];digits=4)])
            end
        end
        push!(tcc,tmp)
    end
    return tcp,tcd,tcc
end
tcp,tcd,tcc = trans_cost(cluijkl,Ipt,Jpt,Kpt,Lpt,Mij,Mjk,Mkl);

function trans_emission(dt,Ipt,Jpt,Kpt,Lpt)
    cep,ced,cec = [],[],[];
    for i=1:length(Ipt)
        em = []
        for j=1:length(Jpt)
            d = dist(Ipt[i],Jpt[j])
            if Mij[i,j]==1
                push!(em,[round(0.065*d; digits=4)])
            else
                push!(em,[round(0.065*d; digits=4),round(0.055*d; digits=4)]);
            end
        end
        push!(cep,em)
    end
    for j=1:length(Jpt)
        em = []
        for k=1:length(Kpt)
            d = dist(Jpt[j],Kpt[k])
            if Mjk[j,k]==1
                push!(em,[round(0.065*d; digits=4)])
            else
                push!(em,[round(0.065*d; digits=4),round(0.055*d; digits=4)]);
            end
        end
        push!(ced,em)
    end
    for k=1:length(Kpt)
        em = []
        for l=1:length(Lpt)
            d = dist(Kpt[k],Lpt[l])
            if Mkl[k,l]==1
                push!(em,[round(0.065*d; digits=4)])
            else
                push!(em,[round(0.065*d; digits=4),round(0.055*d; digits=4)]);
            end
        end
        push!(cec,em)
    end
    return cep,ced,cec
end
cep,ced,cec = trans_emission(dt,Ipt,Jpt,Kpt,Lpt);

function processing_cost(I,J,K)
    vcs = round.(rand(130:150,I).*rand(Uniform(0.9,1.2),1))
    plant1 = round.(rand(130:150,J).*rand(Uniform(0.9,1.2),1))
    plant2 = round.(rand(130*0.9:150*0.9,J).*rand(Uniform(0.9,1.2),1))
    dct1 = round.(rand(100:120,K).*rand(Uniform(0.9,1.2),1))
    dct2 = round.(rand(100*0.9:120*0.9,K).*rand(Uniform(0.9,1.2),1))
    vcp = hcat(plant1,plant2);  vcd = hcat(dct1,dct2)
    return vcs,vcp,vcd
end
vcs,vcp,vcd = processing_cost(I,J,K);
function processing_emission(I,J,K)
    ves = round.(rand(Uniform(2.5,5),I))
    p1 = round.(rand(Uniform(2.5,5),J))
    p2 = round.(rand(Uniform(2.5*0.8,5*0.8),J))
    d1 = round.(rand(Uniform(2.5,5),K))
    d2 = round.(rand(Uniform(2.5*0.8,5*0.8),K))
    vep = hcat(p1,p2);  ved = hcat(d1,d2)
    return ves,vep,ved
end
ves,vep,ved = processing_emission(I,J,K);
function Lcapa(I,J,K,L,Mij,Mjk,Mkl)
    spcost = maximum(maximum(unique(dt.N["LcapacityModesp"])))
    pdcost = maximum(maximum(unique(dt.N["LcapacityModepd"])))
    dccost = maximum(maximum(unique(dt.N["LcapacityModedc"])))
    lsp,lpd,ldc = [],[],[]
    for i=1:I
        lcp = []
        for j=1:J
            if Mij[i,j]==1
                push!(lcp,[0])
            else
                push!(lcp,[0,spcost])
            end
        end
        push!(lsp,lcp)
    end
    for j=1:J
        lcp = []
        for k=1:K
            if Mjk[j,k]==1
                push!(lcp,[0])
            else
                push!(lcp,[0,pdcost])
            end
        end
        push!(lpd,lcp)
    end
    for k=1:K
        lcp = []
        for l=1:L
            if Mkl[k,l]==1
                push!(lcp,[0])
            else
                push!(lcp,[0,dccost])
            end
        end
        push!(ldc,lcp)
    end
    return lsp,lpd,ldc
end
Lcapasp,Lcapapd,Lcapadc = Lcapa(I,J,K,L,Mij,Mjk,Mkl);

function mode_fixed_cost(I,J,K,L,Mij,Mjk,Mkl)
    fsp,fpd,fdc = [],[],[];
    for i=1:I
        fcm = []
        for j=1:J
            m = Mij[i,j]
            if m==1
                push!(fcm,[10000])
            else
                push!(fcm,[10000,0]);
            end
        end
        push!(fsp,fcm)
    end
    for j=1:J
        fcm = []
        for k=1:K
            m = Mjk[j,k]
            if m==1
                push!(fcm,[10000])
            else
                push!(fcm,[10000,0]);
            end
        end
        push!(fpd,fcm)
    end
    for k=1:K
        fcm = []
        for l=1:L
            m = Mkl[k,l]
            if m==1
                push!(fcm,[10000])
            else
                push!(fcm,[10000,0]);
            end
        end
        push!(fdc,fcm)
    end
    return fsp,fpd,fdc
end

gij,gjk,gkl = mode_fixed_cost(I,J,K,L,Mij,Mjk,Mkl);
######################### Build mathematical model #############################
w = [0.5,0.5]; bigM = sum(demand);
function build_scndmodel(w,bigM)
    scnd = direct_model(CPLEX.Optimizer());
    MOI.set(scnd, MOI.RawOptimizerAttribute("CPXPARAM_MIP_Strategy_Search"), 1) # conventional branch-and-cut
    set_silent(scnd)
    MOI.NodeCount()
    # MOI.set(scnd, MOI.RawOptimizerAttribute("CPXPARAM_MIP_Display"),2)
    MOI.set(scnd, MOI.NumberOfThreads(), 1)   #  MOI.set(scnd, MOI.RawOptimizerAttribute("CPXPARAM_Threads"), 1)
    # @variable(scnd, 0<= y1[1:J+K,1:2] <= 1 );
    # @variable(scnd, 0<= uij1[1:I,1:J,1:2] <= 1);
    # @variable(scnd, 0<= ujk1[1:J,1:K,1:2] <= 1);
    # @variable(scnd, 0<= ukl1[1:K,1:L,1:2] <= 1);
    @variable(scnd, y1[1:J+K,1:2], Bin)
    @variable(scnd, uij1[i=1:I,j=1:J,1:Mij[i,j]], Bin)
    @variable(scnd, ujk1[j=1:J,k=1:K,1:Mjk[j,k]], Bin)
    @variable(scnd, ukl1[k=1:K,l=1:L,1:Mkl[k,l]], Bin)

    @variable(scnd, 0<= xij1[i=1:I,j=1:J,1:Mij[i,j]])
    @variable(scnd, 0<= xjk1[j=1:J,k=1:K,1:Mjk[j,k]])
    @variable(scnd, 0<= xkl1[k=1:K,l=1:L,1:Mkl[k,l]])
    @variable(scnd, 0<= h1[1:J+K,1:2])


    #1st obj
    # @constraint(scnd, obj1,
    #     sum(fc[j][t]*y1[j,t] for j=1:J+K for t=1:2) + exg +
    #     exa + exv + sum(vc[j,t]*h1[j,t] for j=1:J+K for t=1:2) <=0);
    #2nd obj
    # @constraint(scnd, obj2, exb+sum(ve[j,t]*h1[j,t] for j=1:J+K for t=1:2) +exr <=0);
    # @objective(scnd,Min,exb+sum(q[j][(p-1)*2+t]*h1[j,t] for j=1:J+K  for t=1:2) +exr );

    @objective(scnd,Min, w[1]*(sum(fc[j][t]*y1[j,t] for j=1:J+K for t=1:2) +
        sum(gij[i][j][m]*uij1[i,j,m] for i=1:I for j=1:J for m=1:Mij[i,j])+
        sum(gjk[j][k][m]*ujk1[j,k,m] for j=1:J for k=1:K for m=1:Mjk[j,k]) +
        sum(gkl[k][l][m]*ukl1[k,l,m] for k=1:K for l=1:L for m=1:Mkl[k,l])) +
        w[2]*(sum([vcp;vcd][j,t]*h1[j,t] for j=1:J+K for t=1:2) +
        sum(vcs[i]*xij1[i,j,m] for i=1:I for j=1:J for m=1:Mij[i,j]) +
        sum(tcp[i][j][m]*xij1[i,j,m] for i=1:I for j=1:J for m=1:Mij[i,j] )+
        sum(tcd[j][k][m]*xjk1[j,k,m] for j=1:J for k=1:K for m=1:Mjk[j,k] )+
        sum(tcc[k][l][m]*xkl1[k,l,m] for k=1:K for l=1:L for m=1:Mkl[k,l] )+
        sum(ves[i]*xij1[i,j,m] for i=1:I for j=1:J for m=1:Mij[i,j] ) +
        sum([vep;ved][j,t]*h1[j,t] for j=1:J+K for t=1:2) +
        sum(cep[i][j][m]*xij1[i,j,m] for i=1:I for j=1:J for m=1:Mij[i,j] )+
        sum(ced[j][k][m]*xjk1[j,k,m] for j=1:J for k=1:K for m=1:Mjk[j,k] )+
        sum(cec[k][l][m]*xkl1[k,l,m] for k=1:K for l=1:L for m=1:Mkl[k,l] )) );
    ########## constraint 3 #############
    @constraints(scnd, begin
        [j=1:J], sum(xij1[i,j,m] for i=1:I for m=1:Mij[i,j]) == sum(xjk1[j,k,m] for k=1:K for m=1:Mjk[j,k]);
        [k=1:K], sum(xjk1[j,k,m] for j=1:J for m=1:Mjk[j,k]) == sum(xkl1[k,l,m] for l=1:L for m=1:Mkl[k,l]);
    end)
    ########### constraint 4-6 #############
    @constraints(scnd, begin
        [j=1:J], sum(h1[j,:]) == sum(xij1[i,j,m] for i=1:I for m=1:Mij[i,j]);
        [k=1:K], sum(h1[k+J,:] ) == sum(xjk1[j,k,m] for j=1:J for m=1:Mjk[j,k]);
        [l=1:L], sum(xkl1[k,l,m] for k=1:K for m=1:Mkl[k,l]) >= demand[l];
    end );
    ########### constraint 7-9 #############
    @constraint(scnd,[i=1:I], -sum(xij1[i,j,m] for j=1:J for m=1:Mij[i,j] ) >= -cas[i]);
    @constraint(scnd,[j=1:J+K, t=1:2], -sum(h1[j,t]) >= -capd[j]*y1[j,t]);
    @constraint(scnd,[j=1:J+K], sum(y1[j,:]) <= 1)
    ########### constraint 10 #############
    @constraint(scnd,[i=1:I,j=1:J], sum(uij1[i,j,m] for m=1:Mij[i,j]) <= 1);
    @constraint(scnd,[j=1:J,k=1:K], sum(ujk1[j,k,m] for m=1:Mjk[j,k]) <= 1);
    @constraint(scnd,[k=1:K,l=1:L], sum(ukl1[k,l,m] for m=1:Mkl[k,l]) <= 1);
    ########### constraint 11 #############
    @constraint(scnd,[i=1:I,j=1:J,m=1:Mij[i,j]], sum(xij1[i,j,m] ) <= bigM*uij1[i,j,m]);
    @constraint(scnd, [j=1:J,k=1:K,m=1:Mjk[j,k]], sum(xjk1[j,k,m] ) <= bigM*ujk1[j,k,m]);
    @constraint(scnd,[k=1:K, l=1:L, m=1:Mkl[k,l]], sum(xkl1[k,l,m] ) <= bigM*ukl1[k,l,m]);
    ########### constraint 12 #############
    @constraint(scnd,[i=1:I, j=1:J, m=1:Mij[i,j]], sum(xij1[i,j,m]) >= Lcapasp[i][j][m]*uij1[i,j,m] );
    @constraint(scnd,[j=1:J, k=1:K, m=1:Mjk[j,k]], sum(xjk1[j,k,m]) >= Lcapapd[j][k][m]*ujk1[j,k,m]);
    @constraint(scnd,[k=1:K, l=1:L, m=1:Mkl[k,l]], sum(xkl1[k,l,m]) >= Lcapadc[k][l][m]*ukl1[k,l,m]);
    ########### constraint 13-14 #############
    @constraint(scnd,sum(y1[j,t] for j=1:J for t=1:2) <= Jmax);
    @constraint(scnd,sum(y1[j,t] for j=J+1:K+J for t=1:2) <= Kmax)
    return scnd
end

scnd = build_scndmodel(w,bigM)
optimize!(scnd)
solve_time(scnd)
@show termination_status(scnd)
objective_value(scnd)
node_count(scnd)
value.(y1) == value.(mas.y)
value.(uij1) == value.(mas.uij)
value.(ujk1)# ==
value.(mas.ujk)

value.(xij1)
value.(xjk1)
value.(xkl1)
value.(h1)
1
################     Benders Decomposition     ################################
# w = [0.5,0.5]
# mas = Model(CPLEX.Optimizer); set_silent(mas)
# MOI.set(mas, MOI.RawParameter("CPX_PARAM_SCRIND"), false )
# MOI.set(mas, MOI.RawParameter("CPX_PARAM_THREADS"),1  )
# @variable(mas, y1[1:J+K,1:2], Bin);
# @variable(mas, uij1[1:I,1:J,1:2], Bin);
# @variable(mas, ujk1[1:J,1:K,1:2], Bin);
# @variable(mas, ukl1[1:K,1:L,1:2], Bin);
# @constraint(mas,[j=1:J+K], -sum(y1[j,:]) >= -1);
# @constraint(mas,[i=1:I,j=1:J], -sum(uij1[i,j,m] for m=1:Mij[i,j]) >= -1);
# @constraint(mas,[j=1:J,k=1:K], -sum(ujk1[j,k,m] for m=1:Mjk[j,k]) >= -1);
# @constraint(mas,[k=1:K,l=1:L], -sum(ukl1[k,l,m] for m=1:Mkl[k,l]) >= -1);
# @constraint(mas,-sum(y1[j,t] for j=J+1:K+J for t=1:2) >= -Kmax)
# exg = AffExpr(0);
# for i=1:I
#     idx = 1;
#     for j=1:J
#         for m=1:Mij[i,j]
#             add_to_expression!(exg, gij[i][idx]*uij1[i,j,m]);
#             idx+=1
#         end
#     end
# end
# @objective(mas, Min, w[1]*(sum(fc[j][t]*y1[j,t] for j=1:J+K for t=1:2) + exg) + );
