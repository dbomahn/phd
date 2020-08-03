using Revise,DelimitedFiles,Plots,DataFrames,CSV,PlotlyJS,RDatasets,Colors
############################ Call data Ozgur (WeightDecomp) ##########################
# path = "C:\\cygwin64\\home\\AK121396\\multiobjective\\solvers\\WeightDecomp\\outputs\\"
# # path = "C:/Users/AK121396/Desktop/initial_outputs/"
# file = readdlm(path*"WD_ILP_n10_sol_var.csv", ';' ,'\n',header=true)
# path = "C:/Users/AK121396/Desktop/initial_outputs/"
# file = readdlm(path5*"Dicho_ILP_n100_sol_var.csv", ';' ,'\n',header=true)
# x = []
# y = []
# z = []
# for j=1:Int16(length(file[1])/2)
#     if file[1][j]!="sol"
#         if file[1][j]!=" "
#             a=split(file[1][j],' ')
#             b1=split(a[1],',')
#             b2 = split(b1[1],'[')
#             c1=split(a[2],',')
#             d1=split(a[3],',')
#             d2 = split(d1[1],']')
#             ps1=parse(Float64,b2[2])
#             ps2=parse(Float64,c1[1])
#             ps3=parse(Float64,d2[1])
#             push!(x,ps1)
#             push!(y,ps2)
#             push!(z,ps3)
#         else
#             break
#         end
#     end
# end

############################### Call data BENSOLVE solutions ########################
# path2 = "C:\\Users\\AK121396\\Desktop\\backupBEN\\res\\"
# fl = readdir(path2)
# p = []
# q = []
# r = []
# for u=1:10
#     file2 = readdlm(path2*fl[u], '\n' ,'\n',header=true)
#     i=9
#     while 8<i<length(file2[1])
#         # file2[1][i]
#         id = split(file2[1][i],' ')
#         idx = filter(x->x != "",id)
#         if idx[1] == "1"
#             ps1 = parse(Float64,idx[2])
#             ps2 = parse(Float64,idx[3])
#             ps3 = parse(Float64,idx[4])
#             push!(p,ps1)
#             push!(q,ps2)
#             push!(r,ps3)
#         elseif idx[1] == "0"
#         else idx[1]
#             break
#         end
#         i+=1
#     end
# end

##########################  Call data Inner Approximation   ##############################
# path3 = "C:\\Users\\AK121396\\Desktop\\backupInner\\"
# file3 = readir(path3*"AP_p-3_n-10_ins-10_t3.out", '\n', '\n',header=false)
# k = []
# n = []
# m = []
# for h=1:length(file3)
#     if split(file3[h],' ')[6]!="Elapsed:"
#         nb=split(file3[h],' ')[7:9]
#         ps1 = parse(Float64,nb[1])
#         ps2 = parse(Float64,nb[2])
#         ps3 = parse(Float64,nb[3])
#         push!(k,ps1)
#         push!(n,ps2)
#         push!(m,ps3)
#     else
#         break
#     end
# end

############################## Call data Kirlik ###########################
# path3 ="/home/ak121396/multiobjective/solvers/KirlikSayin2014/outputs/"
# result = readdir(path3)[1]
# xx = []
# yy = []
# zz = []
# global h = readdlm(path3*result)
# l = Int(length(h)/3)
#
# push!(xx,h[1:l])
# push!(yy,h[l+1: 2*l])
# push!(zz,h[(2*l)+1:3*l])
# plot( markersize = 4, gridalpha=3,xaxis="obj1",yaxis="obj2",zaxis="obj3",title = "KP_n100_ParetoFront")
# Kir= plot!(-x,-y,-z ,seriestype=:scatter, c=:bky, label=:"IP")
# WS=plot!(-p,-q,-r, seriestype=:scatter, c=:yellow, label=:"LP")
# savefig("C:\\Users\\AK121396\\Desktop\\Dropbox\\JKU\\OR2019\\KP_n100_PF.png")


#############################  Plot in 3D  ##############################
# function clustering_alpha_shapes()
#     @eval using DataFrames#, RDatasets, Colors
#
#     # load data
#     iris = dataset("datasets", "iris")
#     nms = unique(iris[:Species])
#     colors = [RGB(0.89, 0.1, 0.1), RGB(0.21, 0.50, 0.72), RGB(0.28, 0.68, 0.3)]
#
#     data = GenericTrace[]
#
#     for (i, nm) in enumerate(nms)
#         df = iris[iris[:Species] .== nm, :]
#         x=df[:SepalLength]
#         y=df[:SepalWidth]
#         z=df[:PetalLength]
#         trace = scatter3d(;name=nm, mode="markers",
#                            marker_size=3, marker_color=colors[i], marker_line_width=0,
#                            x=x, y=y, z=z)
#         push!(data, trace)
#
#         cluster = mesh3d(;color=colors[i], opacity=0.3, x=x, y=y, z=z)
#         push!(data, cluster)
#     end
#
#     # notice the nested attrs to create complex JSON objects
#     layout = Layout(width=800, height=550, autosize=false, title="Iris dataset",
#                     scene=attr(xaxis=attr(gridcolor="rgb(255, 255, 255)",
#                                           zerolinecolor="rgb(255, 255, 255)",
#                                           showbackground=true,
#                                           backgroundcolor="rgb(230, 230,230)"),
#                                yaxis=attr(gridcolor="rgb(255, 255, 255)",
#                                            zerolinecolor="rgb(255, 255, 255)",
#                                            showbackground=true,
#                                            backgroundcolor="rgb(230, 230,230)"),
#                                zaxis=attr(gridcolor="rgb(255, 255, 255)",
#                                            zerolinecolor="rgb(255, 255, 255)",
#                                            showbackground=true,
#                                            backgroundcolor="rgb(230, 230,230)"),
#                                aspectratio=attr(x=1, y=1, z=0.7),
#                                aspectmode = "manual"))
#     plot(data, layout)
# end
# clustering_alpha_shapes()
nms = ["IP","LP"]
colors = [RGB(0.89, 0.1, 0.1),RGB(0.28, 0.68, 0.3)]#[RGB(0.4,0.4,1), RGB(0.8, 0.40, 0)]

data = GenericTrace[]

# for (i, nm) in enumerate(nms)
    # df = iris[iris[:Species] .== nm, :]
    #df[:SepalLength]
    #df[:SepalWidth]
    #df[:PetalLength]
trace1 = scatter3d(;name="IP", mode="markers",
                   marker_size=3, marker_color=colors[2], marker_line_width=0,
                   x=xx, y=yy, z=zz)
trace2 = scatter3d(;name="LP", mode="markers",
                   marker_size=3, marker_color=colors[1], marker_line_width=0,
                   x=p, y=q, z=r)
push!(data, trace1,trace2)

cluster = mesh3d(;color=colors[1], opacity=0.3, x=p, y=q, z=r)
push!(data, cluster)
# end

# notice the nested attrs to create complex JSON objects
layout = Layout(width=800, height=600, autosize=true, title="AP_n10 ParetoFront",
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
                           aspectratio=attr(x=1, y=1, z=0.7),
                           aspectmode = "manual"))
plot(data, layout)
