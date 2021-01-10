using PlotlyJS,RDatasets,Colors,GoldenSequences

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
##########################     Visualisation       ###########################
nms = ["LBsets and IP solution points"]
title = "LB sets and UB sets" #this is the caption appearing in the figure
colors = [RGB(0.89, 0.1, 0.1),RGB(0.28, 0.68, 0.3),RGB(0.4,0.4,1), RGB(0.8, 0.40, 0)]
data = GenericTrace[]
trace1 = PlotlyJS.scatter3d(;name="integer points", mode="markers",
                   marker_size=3, marker_color=colors[1], marker_line_width=0,
                   x=y[:,1], y=y[:,2], z=y[:,3])
trace2 = PlotlyJS.scatter3d(;name="LBset", mode="markers",
                   marker_size=3, marker_color=colors[2], marker_line_width=0,
                   x=-LB[:,1], y=-LB[:,2], z=-LB[:,3])
conv = mesh3d(;color=colors[2], opacity=0.3, x=-LB[:,1], y=-LB[:,2], z=-LB[:,3])
push!(data, trace2,conv) #trace1

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

##############################################################################
f= readdlm("/home/ak121396/Desktop/triflp_Y/01/05_010_01.txtFPep_2hr_Y_.csv")
x1=f[:,1]; y1=f[:,2];z1=f[:,3]

#########################  Plot solutions of 1 instance  #####################
function twoDPlot()
    rw = [10 5; 8 13; 2 7]
    trace1 = scatter(;x=rw[:,1],y=rw[:,2], z=rw[:,3], mode="markers",
                        marker=attr(color="#1f77b4", size=5, symbol="circle",
                                    line=attr(color="rgb(44, 160, 44)", width=0)),
                        line=attr(color="#1f77b4", width=1))
    layout = Layout(autosize=false, width=500, height=500,
                    margin=attr(l=0, r=0, b=0, t=65))
    plot([trace1], layout) # trace2,trace3
end
Plot_solution()
