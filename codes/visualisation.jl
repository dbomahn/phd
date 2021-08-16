using PlotlyJS,DataFrames,DelimitedFiles,Colors,RDatasets
# using PlotlyJS,RDatasets,Colors,GoldenSequences,DelimitedFiles

##########################     3D Visualisation       ###########################
nms = ["FPBH vs BenFPR plot"]
title = "FPBH vs GPR: AP_n-15" #this is the caption appearing in the figure
colors = [RGB(0.89, 0.1, 0.1),RGB(0.28, 0.68, 0.3),RGB(0.4,0.4,1), RGB(0.8, 0.40, 0)]
fpbh = readdlm("/home/ak121396/Desktop/FPBH/AP/GLPK/AP_p-3_n-15_ins-10.ndf")
gpr = readdlm("/home/ak121396/Desktop/GeneralPR/goutputs/AP/GLPK/n-15_ins-10apY.log")
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
rw = [10 5; 8 13; 2 7;]
trace1 = scatter(;x=rw[:,1],y=rw[:,2], mode="markers",
                    marker=attr(color="#1f77b4", size=5, symbol="circle",
                                line=attr(color="rgb(44, 160, 44)", width=0)),
                    line=attr(color="#1f77b4", width=1))
layout = Layout(autosize=false, width=500, height=500,
                margin=attr(l=0, r=0, b=0, t=65))
plot([trace1], layout) # trace2,trace3
Plot_solution()


#######################
function boxmeasure(t1,t2,t3,t4,t5,t6,t7,t8)
    trace1 = box(;y=t1,
                  name="AP_GPR",
                  marker_color="#3D9970")
    trace2 = box(;y=t2,
                  name="AP_FPBH",
                  marker_color="#FF851B")
    trace3 = box(;y=t3,
                  name="KP_GPR",
                  marker_color="#3D9970")
    trace4 = box(;y=t4,
                  name="KP_FPBH",
                  marker_color="#FF851B")
    trace5 = box(;y=t5,
                name="FLP_GPR",
                marker_color="#3D9970")
    trace6 = box(;y=t6,
                name="FLP_FPBH",
                marker_color="#FF851B")
    # trace7 = box(;y=t7,
    #             name="MIP_FPGPR",
    #             marker_color="#3D9970")
    # trace8 = box(;y=t8,
    #             name="MIP_FPBH",
    #             marker_color="#FF851B")
                #               marker_color="rgb(107, 174, 214)" #light blue)
    data = [trace1,trace2,trace3,trace4,trace5,trace6] #,trace7,trace8]
    layout = Layout(;yaxis=attr(boxmode="group", title="CPU time (sec)", zeroline=true, #"HV indicator value"
         zerolinecolor ="#969696", zerolinewidth= 4))#paper_bgcolor="white", plot_bgcolor="white" )type="log"
    PlotlyJS.plot(data, layout)
end
1
function boxmeasure(t1,t2,t3,t4,t5,t6,t7,t8)
    trace1 = box(;y=t1,
                  name="AP_GPR",
                  marker_color="#3D9970")
    trace2 = box(;y=t2,
                  name="AP_FPBH",
                  marker_color="#FF851B")
    trace3 = box(;y=t3,
                  name="KP_GPR",
                  marker_color="#3D9970")
    trace4 = box(;y=t4,
                  name="KP_FPBH",
                  marker_color="#FF851B")
    trace5 = box(;y=t5,
                name="FLP_GPR",
                marker_color="#3D9970")
    trace6 = box(;y=t6,
                name="FLP_FPBH",
                marker_color="#FF851B")
                #               marker_color="rgb(107, 174, 214)",)
    data = [trace1,trace2,trace3,trace4,trace5,trace6]
    # layout = Layout(;yaxis=attr(type="log",boxmode="group", title="HV indicator value",zeroline=true zerolinecolor="rgb(255, 255, 255)",
    #                             zerolinewidth=2,  margin=attr(l=40, r=30, b=80, t=100) ))#paper_bgcolor="white", plot_bgcolor="white" )yaxis_type="log"
    layout = Layout(;yaxis=attr(type="log",autorange=true, showgrid=true, zeroline=true,
                                   dtick=5, gridcolor="rgb(255, 255, 255)",
                                   gridwidth=1,
                                   zerolinecolor="rgb(255, 255, 255)",
                                   zerolinewidth=2),
                        margin=attr(l=40, r=30, b=80, t=100),
                        paper_bgcolor="rgb(243, 243, 243)",
                        plot_bgcolor="rgb(243, 243, 243)",
                        showlegend=false)



    # trace7 = box(;y=t7,
    #             name="MIP_FPGPR",
    #             marker_color="#3D9970")
    # trace8 = box(;y=t8,
    #             name="MIP_FPBH",
    #             marker_color="#FF851B")
                #               marker_color="rgb(107, 174, 214)" #light blue)
    data = [trace1,trace2,trace3,trace4,trace5,trace6] #,trace7,trace8]
    layout = Layout(;yaxis=attr(boxmode="group", title="CPU time (sec)", zeroline=true, #"HV indicator value"
         zerolinecolor ="#969696", zerolinewidth= 4))#paper_bgcolor="white", plot_bgcolor="white" )type="log"
    PlotlyJS.plot(data, layout)
end
# boxmeasure(ag,af,kg,kf,fg,ff,mg,mf)
boxmeasure(gapt,fapt,gkpt,fkpt,gflpt,fflpt) #,gmipt,fmipt)

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
ag,af = CollectVal("F:/results/performance\\GPR/hv/AP/","F:/results/performance/FPBH/hv/AP/")
kg,kf = CollectVal("F:/results\\performance/GPR/hv/KP/","F:\\results/performance/FPBH/hv/KP/")
fg,ff = CollectVal("F:/results/performance/GPR/hv/FLP/","F:/results/performance/FPBH/hv/FLP/")


ag,af = CollectVal("/media/ak121396/0526-8445/results/performance/GPR/hv/AP/","/media/ak121396/0526-8445/results/performance/FPBH/hv/AP/")
kg,kf = CollectVal("/media/ak121396/0526-8445/results/performance/GPR/hv/KP/","/media/ak121396/0526-8445/results/performance/FPBH/hv/KP/")
fg,ff = CollectVal("/media/ak121396/0526-8445/results/performance/GPR/hv/FLP/","/media/ak121396/0526-8445/results/performance/GPR/hv/FLP/")


findall(x->x<0 ,kg)

boxmeasure(agtime,aftime,kgtime,kftime,fgtime,fftime)

mg,mf = CollectVal("F:/results/performance\\GPR/hv/MIPLIB/","F:/results/performance/FPBH/hv/MIPLIB/")

ag,af = CollectVal("F:/results/performance\\GPR/ep/AP/","F:/results/performance/FPBH/ep/AP/")
kg,kf = CollectVal("F:/results\\performance/GPR/ep/KP/","F:\\results/performance/FPBH/ep/KP/")
fg,ff = CollectVal("F:/results/performance/GPR/ep/FLP/","F:/results/performance/FPBH/ep/FLP/")
mg,mf = CollectVal("F:/results/performance\\GPR/ep/MIPLIB/","F:/results/performance/FPBH/ep/MIPLIB/")
boxmeasure(ag,af,kg,kf,fg,ff,mg,mf)

gflpt = vec(readdlm("F:/results/Book1.csv",',',Float64))
gkpt = vec(readdlm("F:/results/Book1.csv",',',Float64))
gapt = vec(readdlm("F:/results/Book1.csv",',',Float64))
gmipt = vec(readdlm("F:/results/Book1.csv",',',Float64))

fapt = vec(readdlm("F:/results/Book1.csv",',',Float64))
fkpt = vec(readdlm("F:/results/Book1.csv",',',Float64))
fflpt = vec(readdlm("F:/results/Book1.csv",',',Float64))
# fmipt = vec(readdlm("F:/results/Book1.csv",',',Float64))

