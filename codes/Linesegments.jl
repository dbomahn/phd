using CPUTime,DataFrames,DelimitedFiles,JuMP,LinearAlgebra,CPLEX,vOptGeneric,SparseArrays,StatsBase,CSV,JLD2
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
file = "/home/ak121396/Desktop/instances/scnd/test01S2"
fname = file[end-7:end]
JLD2.@load "/home/ak121396/Desktop/relise/lpY/X/"*fname*"X.jld2" dv
dt1 = Data1dim(file);
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

################################## Find Line segments  ####################################
function FixedBinDicho(prx)
    linesg = Dict(); dtime = 0
    model = SCND_LP()
    len = [length(model[:y]),length(model[:uij]),length(model[:ujk]),length(model[:ukl]),length(model[:xij]),length(model[:xjk]),length(model[:xkl]),length(model[:h])]

    for k=1:length(prx)
        JuMP.fix.(model[:y], prx[k][1:len[1]]; force = true)
        JuMP.fix.(model[:uij], prx[k][1+len[1]:sum(len[i] for i=1:2)]; force = true)
        JuMP.fix.(model[:ujk], prx[k][1+sum(len[i] for i=1:2):sum(len[i] for i=1:3)]; force = true)
        JuMP.fix.(model[:ukl], prx[k][1+sum(len[i] for i=1:3):sum(len[i] for i=1:4)]; force = true)
        t = @CPUelapsed vSolve(model, 10, method=:dicho, verbose=false)
        res = getvOptData(model);
        if res!= []
            linesg[k]=res.Y_N
        end
        dtime = dtime + t;
    end
    return linesg,dtime
end

linesg,Linetime = FixedBinDicho(dv);
lsg1 = []
lsg2 = collect(values(linesg))
for i=1:length(lsg2)
    lsg1 = vcat(lsg1,lsg2[i])
end
unique(collect(values(linesg)))


1
# function LinesgPostpro(linesg)
#     lsg1 = []
#     lsg2 = collect(values(linesg))
#     for i=1:length(lsg2)
#         lsg1 = vcat(lsg1,lsg2[i])
#     end

#     copyobj = Dict();
#     for i=1:length(lsg1)
#         copyobj[i] = lsg1[i]
#     end
#     for i=1:length(lsg1)-1
#         for j=i+1:length(lsg1)
#             if all(lsg1[i] .>= lsg1[j]) == true #dominated by PF[j]
#                 copyobj[i]=0; break
#             elseif all(lsg1[j] .>= lsg1[i]) == true
#                 copyobj[j]=0;
#             end
#         end
#     end
#     ndlsg = sort!(filter!(a->a!=0, collect(values(copyobj)))) 
#     # lsgtb = hcat([ndlsg[i][1] for i=1:length(ndlsg)],[ndlsg[i][2] for i=1:length(ndlsg)])

#     finalDict = Dict()
#     for l=1:length(ndlsg)
#         for k in collect(keys(linesg))
#             if ndlsg[l] ∈ linesg[k]
#                 if haskey(finalDict,k) == true
#                     push!(finalDict[k],ndlsg[l])
#                 else
#                     push!(finalDict, k => [ndlsg[l]])
#                 end
#             end
#         end
#     end
#     return finalDict
# end
# function LinesgPostpro(lsgfile)

#     JLD2.@load lsgfile linesg;

#     # lsg1 = []
#     # lsg2 = collect(values(linesg))
#     # for i=1:length(lsg2)
#     #     lsg1 = vcat(lsg1,lsg2[i])
#     # end

#     # copyobj = Dict();
#     # for i=1:length(lsg1)
#     #     copyobj[i] = lsg1[i]
#     # end
#     # for i=1:length(lsg1)-1
#     #     for j=i+1:length(lsg1)
#     #         if all(lsg1[i] .>= lsg1[j]) == true #dominated by PF[j]
#     #             copyobj[i]=nothing; break
#     #         elseif all(lsg1[j] .>= lsg1[i]) == true
#     #             copyobj[j]=nothing;
#     #         end
#     #     end
#     # end
#     # ndlsg = sort!(filter!(a->a!=nothing, collect(values(copyobj))))

#     ndlsg = sort!(filter!(a->a!=nothing, lsg1))    
#     lsgtb = hcat([ndlsg[i][1] for i=1:length(ndlsg)],[ndlsg[i][2] for i=1:length(ndlsg)])
#     finalDict = Dict( i=>[] for i=1:length(linesg) )
#     for l=1:size(lsgtb,1)
#         for k in collect(keys(linesg))
#             if lsgtb[l,:] ∈ linesg[k]
#                 push!(finalDict[k],lsgtb[l,:])
#             end
#         end
#     end
#     return finalDict
# end
# lsg = LinesgPostpro("/home/ak121396/Desktop/relise/lpY/linesg/test01S2LS.jld2")
# JLD2.@save "/home/ak121396/Desktop/relise/lpY/linesg/test01S2LS.jld2" lsgdict=linesg;
# JLD2.@load "/home/ak121396/Desktop/relise/lpY/linesg/test01S2LS.jld2" lsgdict;

using PlotlyJS,Colors,JLD2
layout = Layout(
    title="Test",
    xaxis_title="Cost",
    yaxis_title="CO2 emission",
    legend_title="Legend Title",
    font=attr(
        family="Courier New, monospace",
        size=18
    )
)
ef = "/home/ak121396/Desktop/relise/epsilon/1/test01S2epY.log" 
nf = "/home/ak121396/Desktop/relise/lpY/ndp/test01S2lpY.log"

function MIPplot(epfile,ndfile,lsgdict)#lsgfile)
    plot_array = GenericTrace[]
    epY = readdlm(epfile)
    e1 = epY[:,1]; e2 = epY[:,2]
    lpY = readdlm(ndfile)
    l1 = lpY[:,1]; l2 = lpY[:,2]
    ep = scatter(x=e1,y=e2,name="epsilon", mode="markers", marker=attr(color="black"))
    ndp = scatter(x=l1,y=l2,name="LP+FP+FPP+PR", mode="markers", marker=attr(color="green"))
    # JLD2.@load lsgfile lsgdict

    push!(plot_array,ep,ndp) # push!(plot_array,t2)
    
    # ct = count(i->lsgdict[i]!=[],1:length(lsgdict)); cols = distinguishable_colors(ct, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)  

    cols = distinguishable_colors(length(lsgdict), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)  
    
    for i=1:length(lsgdict)
        if lsgdict[i]!=[]
            # push!(lsgdict[i],[l1[i],l2[i]]); 
            sort!(lsgdict[i])
            tradeoffs = scatter(x=[lsgdict[i][j][1] for j=1:length(lsgdict[i])],y=[lsgdict[i][j][2] for j=1:length(lsgdict[i])], mode="markers+lines", color=cols[i])
            push!(plot_array,tradeoffs)
        end
    end
    fig = plot(plot_array, layout); #savefig(fig,"/home/ak121396/Pictures/smSCNDins.png")
end
MIPplot(ef,nf,linesg)#"/home/ak121396/Desktop/relise/lpY/linesg/test01S2LS.jld2")
lpY = readdlm(nf)
l1 = lpY[:,1]; l2 = lpY[:,2]
push!(linesg[1],[l1[1],l2[1]])
cols = distinguishable_colors(length(lsgdict), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)  
plot_array = GenericTrace[]

for i in collect(keys(lsgdict))
    if lsgdict[i]!=[]
        tradeoffs = scatter(x=[lsgdict[i][j][1] for j=1:length(lsgdict[i])],y=[lsgdict[i][j][2] for j=1:length(lsgdict[i])], mode="markers+lines", color=sample(cols,1,replace=false)[1])
        push!(plot_array,tradeoffs)
    end
end
push!(plot_array,t5)
plot(plot_array,layout)

zip(collect(keys(lsgdict)),1: length(lsgdict))
range

[] in collect(values(linesg))

lsg1 = []
lsg2 = collect(values(linesg))
for i=1:length(lsg2)
    lsg1 = vcat(lsg1,lsg2[i])
end



ct = 0
for i=1:length(linesg)
    if linesg[i]!=[]
        ct = ct + length(linesg[i])
    end
end
tb = zeros(ct,2)
iter = 1
for i=1:length(linesg)
    if linesg[i]!=[]
        for j=1:length(linesg[i] )
            tb[iter,:] = linesg[i][j]#[1],linesg[i][j][2]
            iter+=1
        end
    end
end