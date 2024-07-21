using StatsBase,DelimitedFiles,Distributions,CPLEX,JuMP,CPUTime
# "Region: divide the whole grid by 5, one region is 40*40. total 25 regions"
# "Cluster: a dense set of nodes located in the same region"
# "Choose 2(originally 4 or 5) regions (at random) and randomly assign 60% nodes there. 40% are randomly scattered"
# "Max 1 plant and 2 DCs (originally 4 plants 8 DCs) can lie in the same cluster.Candidate DCs are located in the regions where customers populated"
# "Demand ~ [100, 300]. D = ∑∑d^p_{l}"
# "Capacity of DCs~[1.1,1.5], plants~[1.1,2], supplieres~[1.1,3]"
# "The fixed cost of transportation mode 1 is 5% of the total transportation costs. The variable cost of mode 2 is 20% higher than that of mode 1"
# "The whole grid (200*200) divided into 5*5 area to assign t fixed cost of facilities; very expensive areas 5%, expensive areas 50%, intermediate 40%, cheap 5%
# cost_j = 𝚽*sqrt(capacity_j). cheap~[5000,20000], intermediate[20000,35000], expensive~[35000,50000], very expensive~[50000,60000]"
# "Processing cost at suppliers~[130,150], at plants~[130,150], at DCs~[100,120], then multiply by a random noising factor [0.9,1.2]"
# "Transportation cost: Euclidean distance between the arcs was multiplied by a unitary cost~[0.8,1.2] and by
# τ at supplier-plant~[1,1.3], at plant-Dc~[1.2,1.4], at DC-customer~[1.3,1.5]"
struct Data1
    filepath::String; N::Dict{}; d::Array{}; m::Int; c::Array{}; e::Array{}; gij::Array{}; gjk::Array{}; gkl::Array{};
    Mij::Array{}; Mjk::Array{}; Mkl::Array{}; Vij::Array{}; Vjk::Array{}; Vkl::Array{}; b::Array{}; q::Array{};
    upl::Int; udc::Int
    # rij::Array{}; rjk::Array{}; rkl::Array{}; vij::Array{}; vjk::Array{}; vkl::Array{};
    function Data1(filepath)
        data = readdlm(filepath);
        notafile = readdlm("/home/ak121396/Desktop/instances/scnd/Notations.txt", '=');
        # notafile = readdlm("F:/scnd/Notations.txt", '=');
        # notafile = readdlm("/home/k2g00/k2g3475/scnd/Notations.txt", '=');
        nota = notafile[1:end,1];  N= Dict();
        for i=1:length(nota)-1
            id1 = findall(x->x==nota[i], data)[1][1];
            id2 = findall(x->x==nota[i+1], data)[1][1];
            if id2-id1<3
                tmp = filter(x->x!="",  data[id1+(id2-id1-1),:])
                if length(tmp)<2
                    N[nota[i]] = tmp[1];
                else
                    N[nota[i]] = tmp;
                end
            else
                W = []
                for x=id1+1:id1+(id2-id1-1)
                    tmp = filter(x->x!="", data[x,:]);
                    push!(W,tmp);
                end
                # tmp = [filter(x->x!="", data[x,:]) for x in id1+1:id1+(id2-id1-1)]
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

        new(filepath,N,d,m,c,e,gij,gjk,gkl,Mij,Mjk,Mkl,Vij,Vjk,Vkl,b,q,upl,udc); #cap,Mij,Mjk,Mkl,
    end
end
file = "/home/ak121396/Desktop/instances/scnd/notused/test03S4"
# file = "F:/scnd/Test6S3"


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
function trans_cost(Ipt,Jpt,Kpt,Lpt,Mij,Mjk,Mkl)
    tcp,tcd,tcc = [],[],[]
    for i=1:length(Ipt)
        tmp = []
        for j=1:length(Jpt)
            d = dist(Ipt[i],Jpt[j])
            if Mij[i,j]==1
                push!(tmp, [round(d*rand(Uniform(0.8,1.2),1)[1]*rand(Uniform(1,1.3),1)[1]; digits=1)])
            else
                push!(tmp, [round(d*rand(Uniform(0.8,1.2),1)[1]*rand(Uniform(1,1.3),1)[1]; digits=1), round(d*rand(Uniform(0.8*1.2,1.2*1.2),1)[1]*rand(Uniform(1,1.3),1)[1]; digits=1)])
            end
        end
        push!(tcp,tmp)
    end
    for j=1:length(Jpt)
        tmp = []
        for k=1:length(Kpt)
            d = dist(Jpt[j],Kpt[k])
            if Mjk[j,k]==1
                push!(tmp, [round(d*rand(Uniform(0.8,1.2),1)[1]*rand(Uniform(1.2,1.4),1)[1]; digits=1)])
            else
                push!(tmp, [round(d*rand(Uniform(0.8,1.2),1)[1]*rand(Uniform(1.2,1.4),1)[1]; digits=1), round(d*rand(Uniform(0.8*1.2,1.2*1.2),1)[1]*rand(Uniform(1.2,1.4),1)[1]; digits=1)])
            end
        end
        push!(tcd,tmp)
    end
    for k=1:length(Kpt)
        tmp = []
        for l=1:length(Lpt)
            d = dist(Kpt[k],Lpt[l])
            if Mkl[k,l]==1
                push!(tmp, [round(d*rand(Uniform(0.8,1.2),1)[1]*rand(Uniform(1.3,1.5),1)[1]; digits=1)])
            else
                push!(tmp, [round(d*rand(Uniform(0.8,1.2),1)[1]*rand(Uniform(1.3,1.5),1)[1]; digits=1), round(d*rand(Uniform(0.8*1.2,1.2*1.2),1)[1]*rand(Uniform(1.3,1.5),1)[1]; digits=1)])
            end
        end
        push!(tcc,tmp)
    end
    return tcp,tcd,tcc
end

function trans_emission(Ipt,Jpt,Kpt,Lpt)
    cep,ced,cec = [],[],[];
    for i=1:length(Ipt)
        em = []
        for j=1:length(Jpt)
            d = dist(Ipt[i],Jpt[j])
            if Mij[i,j]==1
                push!(em,[round(0.065*d; digits=1)])
            else
                push!(em,[round(0.065*d; digits=1),round(0.055*d; digits=1)]);
            end
        end
        push!(cep,em)
    end
    for j=1:length(Jpt)
        em = []
        for k=1:length(Kpt)
            d = dist(Jpt[j],Kpt[k])
            if Mjk[j,k]==1
                push!(em,[round(0.065*d; digits=1)])
            else
                push!(em,[round(0.065*d; digits=1),round(0.055*d; digits=1)]);
            end
        end
        push!(ced,em)
    end
    for k=1:length(Kpt)
        em = []
        for l=1:length(Lpt)
            d = dist(Kpt[k],Lpt[l])
            if Mkl[k,l]==1
                push!(em,[round(0.065*d; digits=1)])
            else
                push!(em,[round(0.065*d; digits=1),round(0.055*d; digits=1)]);
            end
        end
        push!(cec,em)
    end
    return cep,ced,cec
end
function processing_cost(I,J,K)
    vcs = round.(rand(130:150,I).*rand(Uniform(0.9,1.2),1))
    plant1 = round.(rand(130:150,J).*rand(Uniform(0.9,1.2),1))
    plant2 = round.(rand(130*0.9:150*0.9,J).*rand(Uniform(0.9,1.2),1))
    dct1 = round.(rand(100:120,K).*rand(Uniform(0.9,1.2),1))
    dct2 = round.(rand(100*0.9:120*0.9,K).*rand(Uniform(0.9,1.2),1))
    vcp = hcat(plant1,plant2);  vcd = hcat(dct1,dct2)
    return vcs,vcp,vcd
end
function processing_emission(I,J,K)
    ves = round.(rand(Uniform(2.5,5),I))
    p1 = round.(rand(Uniform(2.5,5),J))
    p2 = round.(rand(Uniform(2.5*0.8,5*0.8),J))
    d1 = round.(rand(Uniform(2.5,5),K))
    d2 = round.(rand(Uniform(2.5*0.8,5*0.8),K))
    vep = hcat(p1,p2);  ved = hcat(d1,d2)
    return ves,vep,ved
end
function Lcapa(I,J,K,Mij,Mjk)
    spcost = maximum(maximum(unique(data.N["LcapacityModesp"])))
    pdcost = maximum(maximum(unique(data.N["LcapacityModepd"])))
    # dccost = maximum(maximum(unique(data.N["LcapacityModedc"])))
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
    # for k=1:K
    #     lcp = []
    #     for l=1:L
    #         if Mkl[k,l]==1
    #             push!(lcp,[0])
    #         else
    #             push!(lcp,[0,dccost])
    #         end
    #     end
    #     push!(ldc,lcp)
    # end
    return lsp,lpd#,ldc
end
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

function Generate_Data_Model(file)
    data = Data1(file);
    rgs = rand(4:5,1)[1]
    xcoordi = sample(1:5, 2); ycoordi = sample(1:5, 2);
    cluij = [xcoordi ycoordi];
    xcoordi = sample(1:5, 2); ycoordi = sample(1:5, 2);
    clukl = [xcoordi ycoordi];
    I = Int(ceil(data.N["supplier"]/3))
    if I<=3
        J=I; K=2*I; L=5*I; Jmax = 2; Kmax = 2
    else
        J=I; K=2*I; L=5*I; Jmax = ceil(J/2); Kmax =I;
    end

    cluijkl = round.(Int,0.6*[I,J,K,L]);
    Mij = data.Mij[1:I,1:J]; Mjk = data.Mjk[1:J,1:K]; Mkl = data.Mkl[1:K,1:L];


    demand = sample(100:300,L);
    cas1 = round.(Int,rand(Uniform(1.1,3), I).*(sum(demand)/Kmax));
    cap1 = round.(Int,rand(Uniform(1.1,2), J)*(sum(demand)/Kmax) );
    cad1 = round.(Int,rand(Uniform(1.1,1.5), K)*(sum(demand)/Kmax) );
    capd1 = [cap1;cad1]; 

    Ipt,Jpt,Kpt,Lpt = nodecoordi(cluijkl,cluij,clukl,I,J,K,L);
    fc = fac_fixedcost(cluijkl,cap1,cad1,J,K);
    tcp,tcd,tcc = trans_cost(Ipt,Jpt,Kpt,Lpt,Mij,Mjk,Mkl);
    cep,ced,cec = trans_emission(Ipt,Jpt,Kpt,Lpt);
    vcs,vcp,vcd = processing_cost(I,J,K);
    ves,vep,ved = processing_emission(I,J,K);
    Lcapasp,Lcapapd = Lcapa(I,J,K,Mij,Mjk);
    gij,gjk,gkl = mode_fixed_cost(I,J,K,L,Mij,Mjk,Mkl);
    bigM = sum(demand);


    w = [1,150]
    scnd = direct_model(CPLEX.Optimizer());set_silent(scnd)
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
    #     sum(fc[j][t]*y1[j,t] for j=1:J+K for t=1:2) +
    #     sum(gij[i][j][m]*uij1[i,j,m] for i=1:I for j=1:J for m=1:Mij[i,j])+
    #     sum(gjk[j][k][m]*ujk1[j,k,m] for j=1:J for k=1:K for m=1:Mjk[j,k]) +
    #     sum(gkl[k][l][m]*ukl1[k,l,m] for k=1:K for l=1:L for m=1:Mkl[k,l]) +
    #     sum(vcs[i]*xij1[i,j,m] for i=1:I for j=1:J for m=1:Mij[i,j]) +
    #     sum(tcp[i][j][m]*xij1[i,j,m] for i=1:I for j=1:J for m=1:Mij[i,j]) +
    #     sum(tcd[j][k][m]*xjk1[j,k,m] for j=1:J for k=1:K for m=1:Mjk[j,k]) +
    #     sum(tcc[k][l][m]*xkl1[k,l,m] for k=1:K for l=1:L for m=1:Mkl[k,l]) +
    #     sum([vcp;vcd][j,t]*h1[j,t] for j=1:J+K for t=1:2)  <=0 );
    # # 2nd obj
    # @constraint(scnd, obj2, sum(ves[i]*xij1[i,j,m] for i=1:I for j=1:J for m=1:Mij[i,j]) +
    #     sum([vep;ved][j,t]*h1[j,t] for j=1:J+K for t=1:2) +
    #     sum(cep[i][j][m]*xij1[i,j,m] for i=1:I for j=1:J for m=1:Mij[i,j] )+
    #     sum(ced[j][k][m]*xjk1[j,k,m] for j=1:J for k=1:K for m=1:Mjk[j,k] )+
    #     sum(cec[k][l][m]*xkl1[k,l,m] for k=1:K for l=1:L for m=1:Mkl[k,l] ) <=0);
    # @objective(scnd,Min,exb+sum(q[j][(p-1)*2+t]*h1[j,t] for j=1:J+K  for t=1:2) +exr );

    @objective(scnd,Min, w[1]*(sum(fc[j][t]*y1[j,t] for j=1:J+K for t=1:2) +
        sum(gij[i][j][m]*uij1[i,j,m] for i=1:I for j=1:J for m=1:Mij[i,j])+
        sum(gjk[j][k][m]*ujk1[j,k,m] for j=1:J for k=1:K for m=1:Mjk[j,k]) +
        sum(gkl[k][l][m]*ukl1[k,l,m] for k=1:K for l=1:L for m=1:Mkl[k,l]) +
        sum(vcs[i]*xij1[i,j,m] for i=1:I for j=1:J for m=1:Mij[i,j]) +
        sum([vcp;vcd][j,t]*h1[j,t] for j=1:J+K for t=1:2) +
        sum(tcp[i][j][m]*xij1[i,j,m] for i=1:I for j=1:J for m=1:Mij[i,j]) +
        sum(tcd[j][k][m]*xjk1[j,k,m] for j=1:J for k=1:K for m=1:Mjk[j,k]) +
        sum(tcc[k][l][m]*xkl1[k,l,m] for k=1:K for l=1:L for m=1:Mkl[k,l]) )
        +
        w[2]*(sum(ves[i]*xij1[i,j,m] for i=1:I for j=1:J for m=1:Mij[i,j]) +
        sum([vep;ved][j,t]*h1[j,t] for j=1:J+K for t=1:2) +
        sum(cep[i][j][m]*xij1[i,j,m] for i=1:I for j=1:J for m=1:Mij[i,j] )+
        sum(ced[j][k][m]*xjk1[j,k,m] for j=1:J for k=1:K for m=1:Mjk[j,k] )+
        sum(cec[k][l][m]*xkl1[k,l,m] for k=1:K for l=1:L for m=1:Mkl[k,l] )) );
    ######### constraint 3 #############
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
    @constraint(scnd,[i=1:I], sum(xij1[i,j,m] for j=1:J for m=1:Mij[i,j] ) <= cas1[i]);
    @constraint(scnd,[j=1:J+K, t=1:2], sum(h1[j,t]) <= capd1[j]*y1[j,t]);
    @constraint(scnd,[j=1:J+K], sum(y1[j,:]) <= 1)
    ########### constraint 10 #############
    @constraint(scnd,[i=1:I,j=1:J], sum(uij1[i,j,m] for m=1:Mij[i,j]) <=  1)
    @constraint(scnd,[j=1:J,k=1:K], sum(ujk1[j,k,m] for m=1:Mjk[j,k]) <= 1)
    @constraint(scnd,[k=1:K,l=1:L], sum(ukl1[k,l,m] for m=1:Mkl[k,l]) <= 1)
    # @constraint(scnd,[i=1:I,j=1:J], sum(uij1[i,j,m] for m=1:Mij[i,j]) <=  sum(y1[j,:]));
    # @constraint(scnd,[j=1:J,k=1:K], sum(ujk1[j,k,m] for m=1:Mjk[j,k]) <= (sum(y1[j,:])+sum(y1[J+k,:]))/2 );
    # @constraint(scnd,[k=1:K,l=1:L], sum(ukl1[k,l,m] for m=1:Mkl[k,l]) <= sum(y1[J+k,:]));
    ########### constraint 11 #############
    @constraint(scnd,[i=1:I,j=1:J,m=1:Mij[i,j]], xij1[i,j,m]  <= bigM*uij1[i,j,m]);
    @constraint(scnd,[j=1:J,k=1:K,m=1:Mjk[j,k]], xjk1[j,k,m]  <= bigM*ujk1[j,k,m]);
    @constraint(scnd,[k=1:K, l=1:L, m=1:Mkl[k,l]], xkl1[k,l,m]  <= bigM*ukl1[k,l,m]);
    ########### constraint 12 #############
    @constraint(scnd,[i=1:I, j=1:J, m=1:Mij[i,j]], xij1[i,j,m] >= Lcapasp[i][j][m]*uij1[i,j,m] );
    @constraint(scnd,[j=1:J, k=1:K, m=1:Mjk[j,k]], xjk1[j,k,m] >= Lcapapd[j][k][m]*ujk1[j,k,m]);
    # @constraint(scnd,[k=1:K, l=1:L, m=1:Mkl[k,l]], sum(xkl1[k,l,m]) >= Lcapadc[k][l][m]*ukl1[k,l,m]);
    # for i=1:I
    #     for j=1:J
    #         if Mij[i,j] == 2
    #             @constraint(scnd, xij1[i,j,2] >= Lcapasp[i][j][2]*uij1[i,j,2] );
    #         end
    #     end
    # end            
    
    # for j=1:J
    #     for k=1:K
    #         if Mjk[j,k] == 2
    #             @constraint(scnd, xjk1[j,k,2] >= Lcapapd[j][k][2]*ujk1[j,k,2]);
    #         end
    #     end
    # end  
    ########### constraint 13-14 #############
    @constraint(scnd,sum(y1[j,t] for j=1:J for t=1:2) <= Jmax);
    @constraint(scnd,sum(y1[j,t] for j=J+1:K+J for t=1:2) <= Kmax)


    optimize!(scnd);
    @show st = termination_status(scnd);
    if st == MOI.OPTIMAL
        arr = []
        push!(arr,I,J,K,L,Jmax,Kmax,"Mij",Mij,"Mjk",Mjk,"Mkl",Mkl,"demand",demand,"cas",cas1,"cap",cap1,"cad",cad1,
            )        
        writedlm("/home/ak121396/Desktop/ins.dat",arr)
    end

end

data = Data1(file);
rgs = rand(4:5,1)[1]
xcoordi = sample(1:5, 2); ycoordi = sample(1:5, 2);
cluij = [xcoordi ycoordi];
xcoordi = sample(1:5, 2); ycoordi = sample(1:5, 2);
clukl = [xcoordi ycoordi];
I = Int(ceil(data.N["supplier"]/3))
if I<=3
    J=I; K=2*I; L=5*I; Jmax = 2; Kmax = 2
else
    J=I; K=2*I; L=5*I; Jmax = ceil(J/2); Kmax =I;
end
# J=I; K=1*I; L=2*I; Jmax = 2; Kmax = 2

cluijkl = round.(Int,0.6*[I,J,K,L]);
Mij = data.Mij[1:I,1:J]; Mjk = data.Mjk[1:J,1:K]; Mkl = data.Mkl[1:K,1:L];

demand = sample(100:300,L);
cas1 = round.(Int,rand(Uniform(1.1,3), I).*(sum(demand)/Kmax));
cap1 = round.(Int,rand(Uniform(1.1,2), J)*(sum(demand)/Kmax) );
cad1 = round.(Int,rand(Uniform(1.1,1.5), K)*(sum(demand)/Kmax) );
capd1 = [cap1;cad1]; 

Ipt,Jpt,Kpt,Lpt = nodecoordi(cluijkl,cluij,clukl,I,J,K,L);
fc = fac_fixedcost(cluijkl,cap1,cad1,J,K);
tcp,tcd,tcc = trans_cost(Ipt,Jpt,Kpt,Lpt,Mij,Mjk,Mkl);
cep,ced,cec = trans_emission(Ipt,Jpt,Kpt,Lpt);
vcs,vcp,vcd = processing_cost(I,J,K);
ves,vep,ved = processing_emission(I,J,K);
Lcapasp,Lcapapd = Lcapa(I,J,K,Mij,Mjk);
gij,gjk,gkl = mode_fixed_cost(I,J,K,L,Mij,Mjk,Mkl);
bigM = sum(demand);


w = [1,150]
scnd = direct_model(CPLEX.Optimizer());set_silent(scnd)
@variable(scnd, y1[1:J+K,1:2], Bin)
@variable(scnd, uij1[i=1:I,j=1:J,1:Mij[i,j]], Bin)
@variable(scnd, ujk1[j=1:J,k=1:K,1:Mjk[j,k]], Bin)
@variable(scnd, ukl1[k=1:K,l=1:L,1:Mkl[k,l]], Bin)

@variable(scnd, 0<= xij1[i=1:I,j=1:J,1:Mij[i,j]])
@variable(scnd, 0<= xjk1[j=1:J,k=1:K,1:Mjk[j,k]])
@variable(scnd, 0<= xkl1[k=1:K,l=1:L,1:Mkl[k,l]])
@variable(scnd, 0<= h1[1:J+K,1:2])

@objective(scnd,Min, w[1]*(sum(fc[j][t]*y1[j,t] for j=1:J+K for t=1:2) +
    sum(gij[i][j][m]*uij1[i,j,m] for i=1:I for j=1:J for m=1:Mij[i,j])+
    sum(gjk[j][k][m]*ujk1[j,k,m] for j=1:J for k=1:K for m=1:Mjk[j,k]) +
    sum(gkl[k][l][m]*ukl1[k,l,m] for k=1:K for l=1:L for m=1:Mkl[k,l]) +
    sum(vcs[i]*xij1[i,j,m] for i=1:I for j=1:J for m=1:Mij[i,j]) +
    sum([vcp;vcd][j,t]*h1[j,t] for j=1:J+K for t=1:2) +
    sum(tcp[i][j][m]*xij1[i,j,m] for i=1:I for j=1:J for m=1:Mij[i,j]) +
    sum(tcd[j][k][m]*xjk1[j,k,m] for j=1:J for k=1:K for m=1:Mjk[j,k]) +
    sum(tcc[k][l][m]*xkl1[k,l,m] for k=1:K for l=1:L for m=1:Mkl[k,l]) )
    +
    w[2]*(sum(ves[i]*xij1[i,j,m] for i=1:I for j=1:J for m=1:Mij[i,j]) +
    sum([vep;ved][j,t]*h1[j,t] for j=1:J+K for t=1:2) +
    sum(cep[i][j][m]*xij1[i,j,m] for i=1:I for j=1:J for m=1:Mij[i,j] )+
    sum(ced[j][k][m]*xjk1[j,k,m] for j=1:J for k=1:K for m=1:Mjk[j,k] )+
    sum(cec[k][l][m]*xkl1[k,l,m] for k=1:K for l=1:L for m=1:Mkl[k,l] )) );
######### constraint 3 #############
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
@constraint(scnd,[i=1:I], sum(xij1[i,j,m] for j=1:J for m=1:Mij[i,j] ) <= cas1[i]);
@constraint(scnd,[j=1:J+K, t=1:2], sum(h1[j,t]) <= capd1[j]*y1[j,t]);
@constraint(scnd,[j=1:J+K], sum(y1[j,:]) <= 1)
########### constraint 10 #############
@constraint(scnd,[i=1:I,j=1:J], sum(uij1[i,j,m] for m=1:Mij[i,j]) <=  1)
@constraint(scnd,[j=1:J,k=1:K], sum(ujk1[j,k,m] for m=1:Mjk[j,k]) <= 1)
@constraint(scnd,[k=1:K,l=1:L], sum(ukl1[k,l,m] for m=1:Mkl[k,l]) <= 1)
########### constraint 11 #############
@constraint(scnd,[i=1:I,j=1:J,m=1:Mij[i,j]], xij1[i,j,m]  <= bigM*uij1[i,j,m]);
@constraint(scnd,[j=1:J,k=1:K,m=1:Mjk[j,k]], xjk1[j,k,m]  <= bigM*ujk1[j,k,m]);
@constraint(scnd,[k=1:K, l=1:L, m=1:Mkl[k,l]], xkl1[k,l,m]  <= bigM*ukl1[k,l,m]);
########### constraint 12 #############
@constraint(scnd,[i=1:I, j=1:J, m=1:Mij[i,j]], xij1[i,j,m] >= Lcapasp[i][j][m]*uij1[i,j,m] );
@constraint(scnd,[j=1:J, k=1:K, m=1:Mjk[j,k]], xjk1[j,k,m] >= Lcapapd[j][k][m]*ujk1[j,k,m]);
########### constraint 13-14 #############
@constraint(scnd,sum(y1[j,t] for j=1:J for t=1:2) <= Jmax);
@constraint(scnd,sum(y1[j,t] for j=J+1:K+J for t=1:2) <= Kmax)


optimize!(scnd);
@show st = termination_status(scnd);
if st == MOI.OPTIMAL

    arr = []
    push!(arr,I,J,K,L,Jmax,Kmax,"Mij")

    for i=1:size(Mij,1)
        push!(arr,Mij[i,:])
    end
    push!(arr, "Mjk") 

    for i=1:size(Mjk,1)
        push!(arr,Mjk[i,:])
    end
    push!(arr,"Mkl")
    for i=1:size(Mkl,1)
        push!(arr,Mkl[i,:])
    end
    push!(arr,"demand")
    for i=1:length(demand)
        push!(arr,demand[i])
    end

    push!(arr,"cas")
    for i=1:length(cas1)
        push!(arr,cas1[i])
    end
    push!(arr,"cap")
    for i=1:length(cap1)
        push!(arr,cap1[i])
    end
    push!(arr,"cad")
    for i=1:length(cad1)
        push!(arr,cad1[i])
    end
    push!(arr,"fixedcost")
    for i=1:length(fc)
        push!(arr,fc[i])
    end

    push!(arr,"tcp")
    for i=1:I
        for j=1:J
            push!(arr,tcp[i][j])
        end
    end

    push!(arr,"tcd")
    for i=1:J
        for j=1:K
            push!(arr,tcd[i][j])
        end
    end
    push!(arr,"tcc")
    for i=1:K
        for j=1:L
            push!(arr,tcc[i][j])
        end
    end

    push!(arr,"cep")
    for i=1:I
        for j=1:J
            push!(arr,cep[i][j])
        end
    end

    push!(arr,"ced")
    for i=1:J
        for j=1:K
            push!(arr,ced[i][j])
        end
    end
    push!(arr,"cec")
    for i=1:K
        for j=1:L
            push!(arr,cec[i][j])
        end
    end

    push!(arr,"vcs")
    for i=1:I
        push!(arr,vcs[i])
    end
    push!(arr,"vcp")
    for j=1:J
        push!(arr,vcp[j,:])
    end
    push!(arr,"vcd")
    for k=1:K
        push!(arr,vcd[k,:])
    end
    push!(arr,"ves")
    for i=1:I
        push!(arr,ves[i])
    end
    push!(arr,"vep")
    for j=1:J
        push!(arr,vep[j,:])
    end
    push!(arr,"ved")
    for k=1:K
        push!(arr,ved[k,:])
    end
    push!(arr,"Lcapasp")
    for i=1:I
        for j=1:J
            push!(arr,Lcapasp[i][j])
        end
    end
    push!(arr,"Lcapapd")
    for j=1:J
        for k=1:K
            push!(arr,Lcapapd[j][k])
        end
    end
    push!(arr,"fixedcostModesp")
    for i=1:I
        for j=1:J
            push!(arr,gij[i][j])
        end
    end
    push!(arr,"fixedcostModepd")
    for j=1:J
        for k=1:K     
            push!(arr,gjk[j][k])   
        end
    end
    push!(arr,"fixedcostModedc")
    for k=1:K        
        for l=1:L
            push!(arr,gkl[k][l])
        end
    end

    writedlm("/home/ak121396/Desktop/instances/test/"*file[end-7:end]*".dat",arr)
end

######################### Build mathematical model #############################
# function build_scndmodel(w,bigM)
   
#     return scnd
# end
# scnd = build_scndmodel([1,150],bigM)


# function getobjval(model)
#     y = value.(model[:y1])
#     uij = value.(model[:uij1])
#     ujk = value.(model[:ujk1])
#     ukl = value.(model[:ukl1])
#     xij = value.(model[:xij1])
#     xjk = value.(model[:xjk1])
#     xkl = value.(model[:xkl1])
#     h = value.(model[:h1])

    
#     obj1 = sum(fc[j][t]*y[j,t] for j=1:J+K for t=1:2) +
#     sum(gij[i][j][m]*uij[i,j,m] for i=1:I for j=1:J for m=1:Mij[i,j])+
#     sum(gjk[j][k][m]*ujk[j,k,m] for j=1:J for k=1:K for m=1:Mjk[j,k]) +
#     sum(gkl[k][l][m]*ukl[k,l,m] for k=1:K for l=1:L for m=1:Mkl[k,l]) +
#     sum(vcs[i]*xij[i,j,1] for i=1:I for j=1:J ) + sum(vcs[i]*xij[i,j,2] for i=1:I for j=1:J if Mij[i,j]==2) +
#     sum([vcp;vcd][j,t]*h[j,t] for j=1:J+K for t=1:2) +
#     sum(tcp[i][j][1]*xij[i,j,1] for i=1:I for j=1:J) + sum(tcp[i][j][2]*xij[i,j,2] for i=1:I for j=1:J if Mij[i,j]==2) +
#     sum(tcd[j][k][1]*xjk[j,k,1] for j=1:J for k=1:K) + sum(tcd[j][k][2]*xjk[j,k,2] for j=1:J for k=1:K if Mjk[j,k]==2) +
#     sum(tcc[k][l][m]*xkl[k,l,m] for k=1:K for l=1:L for m=1:Mkl[k,l]) 
    
#     obj2 = sum(ves[i]*xij[i,j,1] for i=1:I for j=1:J) + sum(ves[i]*xij[i,j,2] for i=1:I for j=1:J if Mij[i,j]==2) +
#     sum([vep;ved][j,t]*h[j,t] for j=1:J+K for t=1:2) +
#     sum(cep[i][j][1]*xij[i,j,1] for i=1:I for j=1:J)+ sum(cep[i][j][2]*xij[i,j,2] for i=1:I for j=1:J if Mij[i,j]==2 )+
#     sum(ced[j][k][1]*xjk[j,k,1] for j=1:J for k=1:K)+ sum(ced[j][k][2]*xjk[j,k,2] for j=1:J for k=1:K if Mjk[j,k]==2 )+
#     sum(cec[k][l][m]*xkl[k,l,m] for k=1:K for l=1:L for m=1:Mkl[k,l] )
#     return obj1,obj2       
# end
