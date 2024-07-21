using PlotlyJS,DataFrames,DelimitedFiles,Colors,CSV,JLD2
#############################        2D plot      ###########################
layout = Layout(
    # title="Test1",
    xaxis_title="Cost",
    yaxis_title="CO2 emission",
    legend_title="Legend Title",
    font=attr(
        family="Courier New, monospace",
        size=18
    )
)
function NDfilter2(Pobj)
    copyobj = Dict();
    for i=1:length(Pobj)
        copyobj[i] = Pobj[i]
    end
    for i=1:length(Pobj)-1
        for j=i+1:length(Pobj)
            if all(Pobj[i] .>= Pobj[j]) == true #dominated by PF[j]
                copyobj[i]=nothing; break
            elseif all(Pobj[j] .>= Pobj[i]) == true
                copyobj[j]=nothing; 
            end
        end
    end
    finalobj = filter!(a->a!=nothing, collect(values(copyobj)))
    return finalobj
end

# benobj = DataFrame(x=sol.LBmtx[:,1], y = sol.LBmtx[:,2]); sort!(benobj, :x);
function Plot_epYmatY(num,folder)
    colour = ["plum","gold","indigo","silver","royalblue","green","olive","teal","cyan", "orchid","steelblue","purple","pink"]
    # epath = "/home/ak121396/Desktop/relise/epsilon/md/2/";  eplist = readdir(epath) # epath = "/home/ak121396/Desktop/relise/epsilon/md/2/"
    lpath = "/home/ak121396/Desktop/relise/vopt/Y/dichow/$folder/"; lplist = readdir(lpath)    
    e2path = "/home/ak121396/Desktop/relise/epsilon/8/"; e2plist = readdir(e2path)
    for i=1:num
        layout = Layout(
        title="Test$i",
        xaxis_title="Cost",
        yaxis_title="CO2 emission",
        legend_title="Legend Title",
        font=attr(
            family="Courier New, monospace",
            size=18
            )
        )
        
        # ep1 = readdlm(epath*eplist[i])#[1:10,:]
        # ep1 = NDfilter2([ep1[i,:] for i=1:size(ep1,1)]);
        # t1 = scatter(x=[ep1[i][1] for i=1:length(ep1)], y=[ep1[i][2] for i=1:length(ep1)],  name="1d_eps", mode="markers", marker=attr(color = "crimson"));
        
        ep2 = readdlm(e2path*e2plist[i])
        ep2 = NDfilter2([ep2[i,:] for i=1:size(ep2,1)]);
        t2 = scatter(x=[ep2[i][1] for i=1:length(ep2)], y=[ep2[i][2] for i=1:length(ep2)],  name="1d_eps2", mode="markers", marker=attr(color = "orange"));
        
        
        # e1 = filter!(i->i!=0, ep0[:,1]); e2 = filter!(i->i!=0, ep0[:,2])
        # ep = [[e1[k],e2[k]] for k=1:length(e1)]
        # epp = NDfilter2(ep)
        # ep1 = [epp[i][1] for i=1:length(epp)]; ep2 = [epp[i][2] for i=1:length(epp)]
        # t1 = scatter(x=ep1[:,1], y=ep1[:,2],  name="epsilon", mode="markers", marker=attr(color = "crimson"))

        Y = readdlm(lpath*lplist[i])
        l1 = Y[:,1]; l2 = Y[:,2]
        t3 = scatter(x=l1,y=l2,name="LP+FP+PR", mode="markers", marker=attr(color="royalblue"))#colour[folder]
        fig = plot([t2,t3], layout)
        # savefig(fig,"/home/ak121396/Desktop/newplots/$folder/$i.png")
        # savefig(fig,"/home/ak121396/Dropbox/SCNDplots/dif/$i.png")
        savefig(fig,"/home/ak121396/Desktop/eps/$i.png")
    end
end
Plot_epYmatY(15,94) 
1
########################################################################file:///home/ak121396/Dropbox/scndfunctions.jl
pr0 =  [ ]
pr1=reshape(pr0,2,Int(length(pr0)/2))
p11,p12 =[],[]
for i=1:Int(length(pr1)/2)
    push!(p11,pr1[1,i])
    push!(p12,pr1[2,i])
end
t4 = scatter(x=p11, y=p12,  name="Matheuristic", mode="markers", marker=attr(color="royalblue")) #+lines
plot([t4,m1], layout)

tnum = 5
# mode="lines+markers+text"
fpath = "/home/ak121396/Desktop/relise/epsilon/5/"
eplist = readdir(fpath)
ep1 = readdlm(fpath*eplist[tnum]);#[1:10,:]
ep1 = NDfilter2([ep1[i,:] for i=1:size(ep1,1)]);
t1 = scatter(x=[ep1[i][1] for i=1:length(ep1)], y=[ep1[i][2] for i=1:length(ep1)],  name="eps", mode="markers", marker=attr(color = "crimson"));
plot([t1],layout)
lpath = "/home/ak121396/Desktop/relise/vopt/Y/final/1/"
llist = readdir(lpath);

# matY = readdlm(lpath*llist[tnum]);
# l1 = matY[:,1]; l2 = matY[:,2]
# t3 = scatter(x=l1,y=l2,name="LP+FP+PR", mode="markers", marker=attr(color="slateblue"));
# plot([t1,t3],layout)# savefig(fig,"/home/ak121396/Dropbox/SCNDplots/dichow/$tnum.png")

# trace0 = scatter(x=benobj[!,:x],y=benobj[!,:y],name="Bensolve",mode="line", market=attr(color="blue"))      # this sets its legend entry
# t1 = scatter(x=ep1[1:10,1], y=ep1[1:10,2],  name="epsilon", mode="markers", marker=attr(color = "crimson"))

lp1= [lp.Y_N[i][1] for i=1:length(lp.Y_N)]; lp2 = [lp.Y_N[i][2] for i=1:length(lp.Y_N)]
t1 = scatter(x=lp1, y=lp2, name="LP", mode="markers+lines", marker=attr(color="royalblue", size=10))
# fig = plot([t1],layout)
f1x,f1y = NDfilter(f1x,f1y)
fp1 = [f1y[i][1] for i=1:length(f1y)]; fp2 = [f1y[i][2] for i=1:length(f1y)]
t2 = scatter(x=fp1,y=fp2,name="FP", mode="markers", marker=attr(color="green", symbol="star-triangle-up", size=10))
# fig = plot([t1,t2],layout);savefig(fig,"/home/ak121396/Pictures/2)LP_FP.png")
t3 = scatter(x = [lex1Y[1][1];lex2Y[1][1]], y = [lex1Y[1][2]; lex2Y[1][2]], name="Lex", mode="markers", marker=attr(color="crimson", symbol="diamond", size=10))
ld11 = [stg1.Y_N[i][1] for i=1:length(stg1.Y_N)]; ld12 = [stg1.Y_N[i][2] for i=1:length(stg1.Y_N)]
ld21 = [stg2.Y_N[i][1] for i=1:length(stg2.Y_N)]; ld22 = [stg2.Y_N[i][2] for i=1:length(stg2.Y_N)]
fig = plot([t1,t2,t3],layout);savefig(fig,"/home/ak121396/Pictures/3)LP_FP_LEX.png")
t4 = scatter(x=ld11,y=ld12, name="LexD", mode="markers+lines", marker=attr(color="salmon", symbol="diamond", size=10))
t5 = scatter(x=ld21,y=ld22, name="LexD", mode="markers+lines", marker=attr(color="salmon", symbol="diamond", size=10))
fig = plot([t1,t2,t3,t4,t5],layout);savefig(fig,"/home/ak121396/Pictures/4)LP_FP_LEXD.png")
fl1 = [dfp.Y[i][1] for i=1:length(dfp.Y)]; fl2 = [dfp.Y[i][2] for i=1:length(dfp.Y)]
t6 = scatter(x=fl1, y=fl2, name="FP&LexD", mode="markers", marker=attr(color="gold", symbol="diamond", size=10))
fig = plot([t1,t2,t4,t5,t6],layout);savefig(fig,"/home/ak121396/Pictures/5)ND:FP&LEXD.png")
fpp1 = [dfpp.Y[i][1] for i=1:length(dfpp.Y)]; fpp2 = [dfpp.Y[i][2] for i=1:length(dfpp.Y)]
t7 = scatter(x=fpp1, y=fpp2, name="FP+", mode="markers", marker=attr(color="lime", symbol="cross", size=10))
plot([t6,t7], layout)
fig = plot([t1,t2,t6,t7],layout);savefig(fig,"/home/ak121396/Pictures/FP+.png")

# nd01= [ndf0.Y[i][1] for i=1:length(ndf0.Y)]; nd02 = [ndf0.Y[i][2] for i=1:length(ndf0.Y)];
# t8 = scatter(x=nd01, y=nd02, name="FP+ D", mode="markers", marker=attr(color="Turquoise", symbol = "star-square",size=10));
# fig = plot([t1,t2,t4,t5,t6,t7,t8],layout);savefig(fig,"/home/ak121396/Pictures/7)FP+ D.png")
nd11= [ndf.Y[i][1] for i=1:length(ndf.Y)]; nd12 = [ndf.Y[i][2] for i=1:length(ndf.Y)];
t9 = scatter(x=nd11, y=nd12, name="PR", mode="markers", marker=attr(color="magenta", symbol = "x",size=10));
# plot([t1,t2,t4,t5,t6,t7,t8,t9],layout)
fig = plot([t1,t2,t7,t9],layout);savefig(fig,"/home/ak121396/Pictures/7)PR.png")
######################################     Line segments       #################################
rpath = "/home/ak121396/Desktop/relise/vopt/nodes/2/"
rlsg = readdir(rpath)
ndy = CSV.read(rpath*rlsg[tnum],DataFrame)

cols = distinguishable_colors(length(ndy.arm), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
lsg_array = GenericTrace[]
start = 0
Larms = findall(x-> ndy.arm[x] == "L", 1:length(ndy.arm))
for l in Larms
    set1 = ndy.v1[1+start:l]; set2 = ndy.v2[1+start:l]
    lsg = scatter(x=[set1[j] for j=1:length(set1)],y=[set2[j] for j=1:length(set2)], mode="markers+lines", color=cols[1])    
    push!(lsg_array,lsg)
    start = l
end
plot([lsg_array; t1], layout) # ;mgplot
# mgpt = readdlm("/home/ak121396/Desktop/relise/performance/5.csv")
# mgplot = scatter(x=mgpt[:,1], y=mgpt[:,2], name="mg", mode="markers", market=attr(color="Terquios"))
################
linesg = DataFrame(nds)
cols = distinguishable_colors(size(linesg,1), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
plot_array = GenericTrace[]
start = 0
Larms = findall(x-> linesg.arm[x] == "L", 1:length(linesg.arm))
for l in Larms
    set1 = linesg.val[1+start:l]#[1]; set2 = linesg.val[1+start:l][2]
    lsg = scatter(x=[set1[j][1] for j=1:length(set1)],y=[set1[j][2] for j=1:length(set1)], mode="markers+lines", color=cols[1])    
    push!(plot_array,lsg)
    start = l
end
plot(plot_array, layout)

cols = distinguishable_colors(length(df.Y), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
plot_array2 = GenericTrace[]
for i=1:length(df.Y)
    dsol,dict= SolveLPdicho(LPdicho,df.X[i],df.Y[i],[df.Y[i]])
    if dsol == []
        # println("dsol is empty")
    else
        tradeoffs = scatter(x=[dsol[j][1] for j=1:length(dsol)],y=[dsol[j][2] for j=1:length(dsol)], mode="markers+lines", color=cols[i])
        push!(plot_array2,tradeoffs)
    end
end
plot([pt;plot_array2], layout) #fig =   savefig(fig,"/home/ak121396/Pictures/smSCNDins.png")

#############################   Benders Plot   ##############################
numcut = [20,40,60,80,100]
title = "B&CTest1"
cpuiter = [0.001	0.303
0.001	0.281
0.0005	0.165
0.001	0.240
0.002	0.527
]
title = "B&CTest4"
cpuiter = [0.021	1.174
0.024	1.177
0.017	1.143
0.023	0.959
0.022	1.078]
p = plot([
    bar(name="CPU", x=numcut, y=cpuiter[:,1],text=cpuiter[:,1], marker_color="lightsalmon"),
    bar(name="Iteration", x=numcut, y=cpuiter[:,2], marker_color="skyblue")],
    Layout(
        xaxis_title_text="Benders cuts(#)",
        yaxis_title_text="Normalised values",
        font=attr( size=22),
        yaxis_range=[0, maximum([cpuiter;1 1])]
    )
)

savefig(p,"/media/ak121396/USB DISK/plots/"*title*".png")


numcut = [20,40,60,80,100]
title = "TBD_Test1"
cpuiter = [
    0.428	0.624
    0.721	0.736
    0.432	0.505
    0.517	0.568
    0.635	0.634]
cpuiter =[0.514	0.624
0.442	0.514
0.677	0.609
0.859	0.697
0.730	0.599
]
p = plot([
    bar(name="CPU", x=numcut, y=cpuiter[:,1], marker_color="orange"),
    bar(name="Iteration", x=numcut, y=cpuiter[:,2], marker_color="darkseagreen")],
    Layout(
        # xaxis_type="category",
        # title_text=title,
        xaxis_title_text="Benders cuts(#)",
        yaxis_title_text="Normalised values",
        font=attr( size=22),
        yaxis_range=[0, maximum([cpuiter;1 1])]
    )
)

# relayout!(p, barmode="group")
savefig(p,"/media/ak121396/USB DISK/plots/"*title*".png")
# savefig(p,"/home/ak121396/Pictures/"*title*".png")


bdfig = plot(
    [
        scatter(x=[0,10,20,30], y=pd.totalCPU, name="CPUT"),
        scatter(x=[0,10,20,30], y=pd.iter, name="#iter", yaxis="y2")
    ],
    Layout(
        xaxis_type="category",
        title_text="B&C Test2_w=$w",
        # title_text="TBD Test2_w=$w",
        xaxis_title_text="%",
        yaxis_title_text="CPUtime(sec)",
        yaxis2=attr(
            title="#subp_iter",
            overlaying="y",
            side="right"
        )
    )
)
# savefig(bdfig,"/media/ak121396/USB DISK/plots/TBD_Test2_w$w.png")

######################################     3D Visualisation       #################################
# using RDatasets,GoldenSequences

# struct NDpoints
#     y::String
#     LB::Array{}
#     LBmtx::Array{}
#     function NDpoints(y)
#         objs = round.(readdlm(y); digits = 4)
#         ind = findall(i -> 0 in objs[i, :], 1:size(objs)[1])
#         LBmtx = objs[setdiff(1:end, ind), 2:end]
#         LB = [Vector(LBmtx[i, :]) for i = 1:size(LBmtx)[1]]
#         new(y, LB, LBmtx)
#     end
# end
nms = ["FPBH vs FPGPR plot"]
title = "FPBH vs FPGPR: " #this is the caption appearing in the figure
colors = [RGB(0.89, 0.1, 0.1),RGB(0.28, 0.68, 0.3),RGB(0.4,0.4,1), RGB(0.8, 0.40, 0)]
fpbh = readdlm("/media/ak121396/0526-8445/results/fpbh/KP/5/KP_p-3_n-070_ ins-03.txt")
gpr = readdlm("/media/ak121396/0526-8445/results/gpr/KP/5/_n-070_ins03kpY.log")
data = GenericTrace[]
points1 = PlotlyJS.scatter3d(;name="FPBH", mode="markers",
                   marker_size=3, marker_color=colors[1], marker_line_width=0,
                   x=fpbh[:,1], y=fpbh[:,2], z=fpbh[:,3])
points2 = PlotlyJS.scatter3d(;name="GPR", mode="markers",
                   marker_size=3, marker_color=colors[2], marker_line_width=0,
                   x=gpr[:,1], y=gpr[:,2], z=gpr[:,3])
# conv1 = mesh3d(;color=colors[1], opacity=0.3, x=fpbh[:,1], y=fpbh[:,2], z=fpbh[:,3])
# conv2 = mesh3d(;color=colors[2], opacity=0.3, x=bfp[:,1], y=bfp[:,2], z=bfp[:,3])

push!(data, points1, points2) #, conv1, conv2)

layout = Layout(width=600, height=600, autosize=true, title=title,
                scene=attr(xaxis=attr(gridcolor="rgb(255, 255, 255)",
                                      zerolinecolor="rgb(255, 255, 255)",
                                      showbackground=true,
                                      backgroundcolor="rgb(230, 230,230)"),
                           yaxis=attr(gridcolor="rgb(255, 255, 255)",
                                       zerolinecolor="rgb(255, 255, 255)",
                                       showbackground=true,
                                       backgroundcolor="rgb(230, 230,230)"),
                           zaxis=attr(gridcolor="rgb(255, 255, 255)",
                                       zerolinecolor="rgb(255, 255, 255)",
                                       showbackground=true,
                                       backgroundcolor="rgb(230, 230,230)"),
                           aspectratio=attr(x=1, y=1, z=1),
                           aspectmode = "manual"))
PlotlyJS.plot(data, layout)

###############################  Clustering  #################################
epx=points[1,:]; epy=points[2,:]; epz=points[3,:]
nms = ["Clustering LBsets"]; title = "Clustering LBsets"
colors = [RGB(0.89, 0.1, 0.1),RGB(0.28, 0.68, 0.3),RGB(0.4,0.4,1), RGB(0.8, 0.40, 0)]
data = GenericTrace[]
trace3 = PlotlyJS.scatter3d(;name=title, mode="markers",
                   marker_size=3, marker_color=colors[2], marker_line_width=0,
                   x=epx, y=epy, z=epz)
# trace2 = PlotlyJS.scatter3d(;name="FP+LS", mode="markers",
#                    marker_size=3, marker_color=colors[2], marker_line_width=0,
#                    x=PFx, y=PFy, z=PFz)
push!(data, trace3) #,trace2)
# cluster = mesh3d(;color=colors[1], opacity=0.3, x=x1, y=y1, z=z1) #for masch
# push!(data,cluster)

layout = Layout(width=800, height=600, autosize=true, title=title,
                scene=attr(xaxis=attr(gridcolor="rgb(255, 255, 255)",
                                      zerolinecolor="rgb(255, 255, 255)",
                                      showbackground=true,
                                      backgroundcolor="rgb(230, 230,230)"),
                           yaxis=attr(gridcolor="rgb(255, 255, 255)",
                                       zerolinecolor="rgb(255, 255, 255)",
                                       showbackground=true,
                                       backgroundcolor="rgb(230, 230,230)"),
                           zaxis=attr(gridcolor="rgb(255, 255, 255)",
                                       zerolinecolor="rgb(255, 255, 255)",
                                       showbackground=true,
                                       backgroundcolor="rgb(230, 230,230)"),
                           aspectratio=attr(x=1, y=1, z=1),
                           aspectmode = "manual"))
PlotlyJS.plot(data, layout)


############################## LB clustering ##############
# clr = map(x->RGB(x...), (Iterators.take(GoldenSequence(3), 20)))
function LPclustering(points)
    # load data
    nms = unique(points[:cluster])
    colors = map(x->RGB(x...), (Iterators.take(GoldenSequence(3), length(nms)+1)))
    # colors = palette("Paired",9)
    # colors = [RGB(0.89, 0.1, 0.1), RGB(0.21, 0.50, 0.72), RGB(0.28, 0.68, 0.3)]

    data = GenericTrace[]
    for (i, j) in enumerate(nms) #i is a number, and j is cluster
        df = points[points[:cluster] .== j, :]
        g1 = df[:obj1]
        g2 = df[:obj2]
        g3 = df[:obj3]

        tr = scatter3d(;name="cluster"*"$j", mode="markers", marker_color=colors[i+1],
                           marker_size=3,  marker_line_width=0,
                           x=-g1, y=-g2, z=-g3)
        push!(data, tr)

        cluster = mesh3d(;color=colors[i+1], opacity=0.3, x=-g1, y=-g2, z=-g3)
        push!(data, cluster)
    end

    # notice the nested attrs to create complex JSON objects
    layout = Layout(width=800, height=550, autosize=true, title="LB clusters",
                    scene=attr(xaxis=attr(gridcolor="rgb(255, 255, 255)",
                                          zerolinecolor="rgb(255, 255, 255)",
                                          showbackground=true,
                                          backgroundcolor="rgb(230, 230,230)"),
                               yaxis=attr(gridcolor="rgb(255, 255, 255)",
                                           zerolinecolor="rgb(255, 255, 255)",
                                           showbackground=true,
                                           backgroundcolor="rgb(230, 230,230)"),
                               zaxis=attr(gridcolor="rgb(255, 255, 255)",
                                           zerolinecolor="rgb(255, 255, 255)",
                                           showbackground=true,
                                           backgroundcolor="rgb(230, 230,230)"),
                               aspectratio=attr(x=1, y=1, z=1),
                               aspectmode = "manual"))
    plot(data, layout)
end
LPclustering(points)

##############################################################################
f= readdlm("/home/ak121396/Desktop/triflp_Y/01/05_010_01.txtFPep_2hr_Y_.csv")
x1=f[:,1]; y1=f[:,2];z1=f[:,3]

#########################  Plot solutions of 1 instance  #####################
obj = [260.1	2553.1;
    236.7	2713.6;
    238.5	2591.2;
    186.5893305391	2533.3313103322;
    206.8261297131	2510.917672409;
    220.3257868363	2499.366854594;
    178.700947814	2843.215863068;
    178.767837814	2830.384863068;
    177.6211584268	3017.524146964;
    179.330821814	2825.056863068;
    181.363630051	2820.803528352;
    177.9000429376	2969.20966591067]
trace1 = scatter(;x=obj[:,1],y=obj[:,2], mode="markers",
                    marker=attr(color="#1f77b4", size=5, symbol="circle",
                                line=attr(color="rgb(44, 160, 44)", width=0)),
                    line=attr(color="#1f77b4", width=1))
layout = Layout(autosize=false, width=500, height=500, margin=attr(l=30, r=0, b=0, t=60),
        xaxis_title = "Cost",yaxis_title = "CO2") #, xaxis_range = [220, 380], yaxis_range=[1000,3000])
plot([trace1], layout) # trace2,trace3

Plot_solution()


#######################

function CollectVal(path1,path2)
    iter = readdir(path1)
    trace1 = [];trace2 = [];
    for i=1:length(iter)
        gins = readdir(path1*iter[i]); fins = readdir(path2*iter[i]);
        for j=1:length(gins)
            gep = readdlm(path1*iter[i]*"/"*gins[j])[1]; fep = readdlm(path2*iter[i]*"/"*fins[j])[1];
            push!(trace1,gep); push!(trace2,fep)
        end
    end
    return trace1,trace2
end

ag,af = CollectVal("/media/ak121396/0526-8445/results/performance/GPR/ep/AP/","/media/ak121396/0526-8445/results/performance/FPBH/ep/AP/")
fg,ff = CollectVal("/media/ak121396/0526-8445/results/performance/GPR/ep/FLP/","/media/ak121396/0526-8445/results/performance/FPBH/ep/FLP/")
kg,kf = CollectVal("/media/ak121396/0526-8445/results/performance/GPR/ep/KP/","/media/ak121396/0526-8445/results/performance/FPBH/ep/KP/")
mg,mf = CollectVal("/media/ak121396/0526-8445/results/performance/GPR/ep/MIPLIB/","/media/ak121396/0526-8445/results/performance/FPBH/ep/MIPLIB/")

ag,af = CollectVal("F:/results/performance\\GPR/ep/AP/","F:/results/performance/FPBH/ep/AP/")
kg,kf = CollectVal("F:/results\\performance/GPR/ep/KP/","F:\\results/performance/FPBH/ep/KP/")
fg,ff = CollectVal("F:/results/performance/GPR/ep/FLP/","F:/results/performance/FPBH/ep/FLP/")
mg,mf = CollectVal("F:/results/performance\\GPR/ep/MIPLIB/","F:/results/performance/FPBH/ep/MIPLIB/")


x0 = []
for i=1:500
    push!(x0,"MOAP")
end
for i=1:500
    push!(x0,"MOKP")
end
for i=1:600
    push!(x0,"TOFLP")
end
for i=1:45
    push!(x0,"TOMIPLIB")
end

function HVbox(x0,t1,t3)
    trace1 = box(;y=t1, x=x0, name="FPBH",  marker_color="#FF4136")
    # trace2 = box(;y=t2, x=x0[1:1100], name="FFP",marker_color= "#F1A232")
    trace3 = box(;y=t3, x=x0, name="LPBM", marker_color="#3D9970")
    data = [trace1, trace3] # trace2,
    layout = Layout(;yaxis=attr(title="HV indicator value", zeroline=false),boxmode="group")
    plot(data, layout)
end
HVbox(x0,vcat(af,kf,ff,mf),vcat(ag,kg,fg,mg))

function epbox(x0,t1,t3)
    trace1 = box(;y=t1, x=x0, name="FPBH",  marker_color="#FF4136")
    # trace2 = box(;y=t2, x=x0[1:1100], name="LPBM",  marker_color= "#3D9970")#"#F1A232")
    trace3 = box(;y=t3, x=x0, name="LPBM", marker_color="#3D9970")
    data = [trace1, trace3]
    layout = Layout(;yaxis=attr(title="unary ϵ indicator value", zeroline=false),boxmode="group")
    plot(data, layout)
end
ag,af = CollectVal("/home/ak121396/Desktop/LPBMresults/performance/GPR/ep/AP/","/home/ak121396/Desktop/LPBMresults/performance/FPBH/ep/AP/")
fg,ff = CollectVal("/home/ak121396/Desktop/LPBMresults/performance/GPR/ep/FLP/","/home/ak121396/Desktop/LPBMresults/performance/FPBH/ep/FLP/")
kg,kf = CollectVal("/home/ak121396/Desktop/LPBMresults/performance/GPR/ep/KP/","/home/ak121396/Desktop/LPBMresults/performance/FPBH/ep/KP/")
mg,mf = CollectVal("/home/ak121396/Desktop/LPBMresults/performance/GPR/ep/MIPLIB/","/home/ak121396/Desktop/LPBMresults/performance/FPBH/ep/MIPLIB/")

epbox(x0,vcat(af,kf,ff,mf),vcat(ag,kg,fg,mg))
1

# function groupbox(t1,t2,t3,t4,t5,t6,t7,t8)
#     trace1 = box(;y=t1,name="AP_FPBH", marker_color="rgb(255, 76, 27)")
#     trace2 = box(;y=t2,name="AP_GPR", marker_color="#3D9970")
#     trace3 = box(;y=t3, name="KP_FPBH", marker_color="rgb(255, 76, 27)")
#     trace4 = box(;y=t4, name="KP_GPR", marker_color="#3D9970")
#     trace5 = box(;y=t5, name="FLP_FPBH", marker_color="rgb(255, 76, 27)")
#     trace6 = box(;y=t6, name="FLP_GPR", marker_color="#3D9970")
#     trace7 = box(;y=t7, name="neos151", marker_color="rgb(255, 76, 27)")
#     trace8 = box(;y=t8, name="neos151", marker_color="#3D9970")
#     data = [trace1,trace2,trace3,trace4,trace5,trace6,trace7,trace8]
#     layout = Layout(;yaxis=attr(boxmode="group", title="HV indicator value"))#unary epsilon indicator value
#     PlotlyJS.plot(data, layout)
# end
# function mipHVbox(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10)
#     trace1 = box(;y=t1,
#                   name="cvs08",
#                   marker_color="rgb(255, 76, 27)")
#     trace2 = box(;y=t2,
#                   name="cvs08",
#                   marker_color="rgb(54, 128, 56)")
#     trace3 = box(;y=t3,
#                   name="cvs16",
#                   marker_color="rgb(255, 76, 27)")
#     trace4 = box(;y=t4,
#                   name="cvs16",
#                  marker_color="rgb(54, 128, 56)")
#     trace5 = box(;y=t5,
#                 name="n2seq",
#                 marker_color="rgb(255, 76, 27)")
#     trace6 = box(;y=t6,
#                 name="n2seq",
#                 marker_color="rgb(54, 128, 56)")
#                 #               marker_color="rgb(107, 174, 214)",)
#     data = [trace1,trace2,trace3,trace4,trace5,trace6]
#     trace9 = box(;y=t9,
#                 name="neos159",
#                 marker_color="rgb(255, 76, 27)")
#     trace10 = box(;y=t10,
#                 name="neos159",
#                 marker_color="rgb(54, 128, 56)")
#     data = [trace1,trace2,trace3,trace4,trace5,trace6,trace7,trace8,trace9,trace10]
#     layout = Layout(;xaxis = attr(title = "TOMIPLIB instances"),
#         yaxis=attr(boxmode="group", title="HV indicator value", type="log"))#unary epsilon indicator value
#         # zeroline=true, zerolinecolor ="#969696", zerolinewidth= 4,
#         # paper_bgcolor="white",plot_bgcolor="white", #type="log"
#     PlotlyJS.plot(data, layout)
# end
