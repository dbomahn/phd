using PlotlyJS,DataFrames,DelimitedFiles,Colors,JLD2
#############################        2D plot      ###########################
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


function Plot_epYlpY(num)
    eplist = readdir("/home/ak121396/Desktop/relise/epsilon/")
    lplist = readdir("/home/ak121396/Desktop/relise/lpY/2/")
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
        
        ep1 = readdlm("/home/ak121396/Desktop/relise/epsilon/"*eplist[i+1])[1:10,:]

        # e1 = filter!(i->i!=0, ep0[:,1]); e2 = filter!(i->i!=0, ep0[:,2])
        # ep = [[e1[k],e2[k]] for k=1:length(e1)]
        # epp = NDfilter2(ep)
        # ep1 = [epp[i][1] for i=1:length(epp)]; ep2 = [epp[i][2] for i=1:length(epp)]

        
        lpY = readdlm("/home/ak121396/Desktop/relise/lpY/2/"*lplist[i])
        l1 = lpY[:,1]; l2 = lpY[:,2]
        t1 = scatter(x=ep1[:,1], y=ep1[:,2],  name="epsilon", mode="markers", marker=attr(color = "crimson"))
        t3 = scatter(x=l1,y=l2,name="LP+FP+PR", mode="markers", marker=attr(color="green"))
        fig = plot([t1,t3], layout)
        savefig(fig,"/home/ak121396/Dropbox/5_w/$i.png")
    end
end
Plot_epYlpY(15)
########################################################################
pr0 =  [ ]
pr1=reshape(pr0,2,Int(length(pr0)/2))
p11,p12 =[],[]
for i=1:Int(length(pr1)/2)
    push!(p11,pr1[1,i])
    push!(p12,pr1[2,i])
end
t4 = scatter(x=p11, y=p12,  name="LP+FP+PR", mode="markers", marker=attr(color="royalblue")) #
plot([t1,t4], layout)

# mip = readdlm("/home/desk/Desktop/relise/mip/test04S4mipY2.log")
m1 = mip[:,1]; m2 = mip[:,2]

# tnum = 12
fpath = "/home/ak121396/Desktop/relise/epsilon/1/"
eplist = readdir(fpath)#[2:end]
ep1 = readdlm(fpath*eplist[tnum])[1:10,:]
t1 = scatter(x=ep1[:,1], y=ep1[:,2],  name="epsilon", mode="markers", marker=attr(color = "crimson"))
# e1 = filter!(i->i!=0, ep0[:,1]); e2 = filter!(i->i!=0, ep0[:,2])
# ep = [[e1[k],e2[k]] for k=1:length(e1)]
# epp = NDfilter2(ep)
# ep1 = [epp[i][1] for i=1:length(epp)]; ep2 = [epp[i][2] for i=1:length(epp)]

lpY = readdlm("/home/ak121396/Desktop/relise/lpY/test12S4largeTL.log")
l1 = lpY[:,1]; l2 = lpY[:,2]
t3 = scatter(x=l1,y=l2,name="LP+FP+PR", mode="markers", marker=attr(color="green"))
plot([t1,t3],layout)

###########################################
# lsg = filter!(p->p!=[],collect(values(linesg)))
# lsg1 = []; lsg2 = []
# for i=1:length(lsg)
#     append!(lsg1,[lsg[i][j][1] for j=1:length(lsg[i])])
#     append!(lsg2,[lsg[i][j][2] for j=1:length(lsg[i])])
# end
# opt dominated out
# fy1= [fyy[i][1] for i=1:length(fyy)]; fy2 = [fyy[i][2] for i=1:length(fyy)]

# mode="markers+text"
# trace0 = scatter(x=benobj[!,:x],y=benobj[!,:y],name="Bensolve",mode="line", market=attr(color="blue"))      # this sets its legend entry
t1 = scatter(x=ep1[1:10,1], y=ep1[1:10,2],  name="epsilon", mode="markers", marker=attr(color = "crimson"))
t2 = scatter(x=l1,y=l2,name="LP+FP+PR", mode="markers", marker=attr(color="green"))
py1= [pry[i][1] for i=1:length(pry)]; py2 = [pry[i][2] for i=1:length(pry)]
t3 = scatter(x=[py1;], y=py2, name="LP+FP+PR", mode="markers", marker=attr(color="royalblue"))
# lp1 = [lp.Y_N[i][1] for i=1:length(lp.Y_N)]; lp2 = [lp.Y_N[i][2] for i=1:length(lp.Y_N)]


fy1= [dfpp.Y[i][1] for i=1:length(dfpp.Y)]; fy2 = [dfpp.Y[i][2] for i=1:length(dfpp.Y)]
t1 = scatter(x=fy1, y=fy2,  name="LP+FP+FPP", mode="markers", marker=attr(color="orange")) #Turquios
############### merging nd segments

ndog1 = [ndsetog[i].val[1] for i=1:length(ndsetog)]; ndog2 = [ndsetog[i].val[2] for i=1:length(ndsetog)];
t0 = scatter(x=ndog1, y=ndog2, name="ndset1", mode="markers+lines", market=attr(color="black"))
nd12x = [set2[i].val[1] for i=1:length(set2)]; nd12y = [set2[i].val[2] for i=1:length(set2)];
t12 = scatter(x=nd12x, y=nd12y, name="ndset2", mode="markers", market=attr(color="orange"))
nset3x = [set3[i].val[1] for i=1:length(set3)]; nset3y = [set3[i].val[2] for i=1:length(set3)];
t3 = scatter(x=nset3x, y=nset3y, name="ndset3", mode="markers", market=attr(color="red"))
nset4x = [set4[i].val[1] for i=1:length(set4)]; nset4y = [set4[i].val[2] for i=1:length(set4)];
t4 = scatter(x=nset4x, y=nset4y, name="ndset4", mode="markers", market=attr(color="Terquios"))
dsol = dsol2
dicho1 = [dsol[i][1] for i=1:length(dsol)]; dicho2 = [dsol[i][2] for i=1:length(dsol)]
dcho = scatter(x=dicho1, y=dicho2, name="dicho", mode="markers+lines", market=attr(color="green"))
plot([dcho,t0,t12],layout)

iid =findall( x-> x.val[2] < 892500 && x.val[2] > 892300, set3 )[1]
set3[iid]
plot([trace2,trace3,trace1,trace5], layout)

cols = distinguishable_colors(length(linesg), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
plot_array = GenericTrace[]
push!(plot_array,trace)
for i=1:length(linesg)
    if linesg[i]!=[]
        tradeoffs = scatter(x=[linesg[i][j][1] for j=1:length(linesg[i])],y=[linesg[i][j][2] for j=1:length(linesg[i])], mode="markers+lines", color=cols[i])
        push!(plot_array,tradeoffs)
    end
end
fig = plot(plot_array, layout); #savefig(fig,"/home/ak121396/Pictures/smSCNDins.png")




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
ag,af = CollectVal("/media/ak121396/0526-8445/results/performance/GPR/ep/AP/","/media/ak121396/0526-8445/results/performance/FPBH/ep/AP/")
fg,ff = CollectVal("/media/ak121396/0526-8445/results/performance/GPR/ep/FLP/","/media/ak121396/0526-8445/results/performance/FPBH/ep/FLP/")
kg,kf = CollectVal("/media/ak121396/0526-8445/results/performance/GPR/ep/KP/","/media/ak121396/0526-8445/results/performance/FPBH/ep/KP/")
mg,mf = CollectVal("/media/ak121396/0526-8445/results/performance/GPR/ep/MIPLIB/","/media/ak121396/0526-8445/results/performance/FPBH/ep/MIPLIB/")

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
