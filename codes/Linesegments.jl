function LinesgPostpro(lsgfile)

    JLD2.@load lsgfile linesg;

    # lsg1 = []
    # lsg2 = collect(values(linesg))
    # for i=1:length(lsg2)
    #     lsg1 = vcat(lsg1,lsg2[i])
    # end

    # copyobj = Dict();
    # for i=1:length(lsg1)
    #     copyobj[i] = lsg1[i]
    # end
    # for i=1:length(lsg1)-1
    #     for j=i+1:length(lsg1)
    #         if all(lsg1[i] .>= lsg1[j]) == true #dominated by PF[j]
    #             copyobj[i]=nothing; break
    #         elseif all(lsg1[j] .>= lsg1[i]) == true
    #             copyobj[j]=nothing;
    #         end
    #     end
    # end
    # ndlsg = sort!(filter!(a->a!=nothing, collect(values(copyobj))))

    ndlsg = sort!(filter!(a->a!=nothing, lsg1))    
    lsgtb = hcat([ndlsg[i][1] for i=1:length(ndlsg)],[ndlsg[i][2] for i=1:length(ndlsg)])
    finalDict = Dict( i=>[] for i=1:length(linesg) )
    for l=1:size(lsgtb,1)
        for k in collect(keys(linesg))
            if lsgtb[l,:] ∈ linesg[k]
                push!(finalDict[k],lsgtb[l,:])
            end
        end
    end
    return finalDict
end

lsg = LinesgPostpro("/home/ak121396/Desktop/relise/lpY/linesg/test01S2LS.jld2")
JLD2.@load "/home/ak121396/Desktop/relise/lpY/linesg/test02S3LS.jld2" lsgdict;

lsgdict


function MIPplot(epfile,ndfile,lsgfile)
    plot_array = GenericTrace[]
    epY = readdlm(epfile)
    e1 = epY[:,1]; e2 = epY[:,2]
    lpY = readdlm(ndfile)
    l1 = lpY[:,1]; l2 = lpY[:,2]
    ep = scatter(x=e1,y=e2,name="epsilon", mode="markers", marker=attr(color="red"))
    ndp = scatter(x=l1,y=l2,name="LP+FP+PR", mode="markers", marker=attr(color="green"))
    JLD2.@load lsgfile lsgdict

    push!(plot_array,ep,ndp) # push!(plot_array,t2)
    
    # ct = count(i->lsgdict[i]!=[],1:length(lsgdict)); cols = distinguishable_colors(ct, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)  

    cols = distinguishable_colors(length(lsgdict), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)  
    
    for i=1:length(lsgdict)
        if lsgdict[i]!=[]
            tradeoffs = scatter(x=[lsgdict[i][j][1] for j=1:length(lsgdict[i])],y=[lsgdict[i][j][2] for j=1:length(lsgdict[i])], mode="markers+lines", color=cols[i])
            push!(plot_array,tradeoffs)
        end
    end
    fig = plot(plot_array, layout); #savefig(fig,"/home/ak121396/Pictures/smSCNDins.png")
end
MIPplot("/home/ak121396/Desktop/relise/epsilon/1/test02S3epY.log","/home/ak121396/Desktop/relise/lpY/ndp/test02S3lpY.log","/home/ak121396/Desktop/relise/lpY/linesg/test02S3LS.jld2")


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

linesg

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