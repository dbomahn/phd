using Polyhedra,QHull,DataFrames,DelimitedFiles
############################### Call data BENSOLVE solutions ########################
path2 = "C:\\cygwin64\\home\\AK121396\\backupBEN\\working\\"
# fl = readdir(path2)
#
# for u=1:length(fl)
    # file2 = readdlm(path2*fl[1], '\n' ,'\n',header=true)
    # x = []
    # i = 9
    # while 8<i<length(file2[1])
    #     global i
    #     id = split(file2[1][i],' ')
    #     idx = filter(x->x != "",id)
    #     if idx[1] == "1"
    #         ps1 = parse(Float64,idx[2])
    #         ps2 = parse(Float64,idx[3])
    #         ps3 = parse(Float64,idx[4])
    #         push!(x,ps1,ps2,ps3)
    #     elseif idx[1] == "0"
    #     else idx[1]
    #         break
    #     end
    #     i+=1
    # end
    # h = Int(length(x)/3)
    # x = convert(Array{Float64,1}, x)
    # p = reshape(x, h, 3)
    # n = 100 # 5 10 15 30 50 70 100
    # coeff = 80 # 20 1000 80
    # ref = [n*coeff n*coeff n*coeff]
    # p = vcat(p,ref)
    # ph = polyhedron(vrep(p), QHull.Library())
    # hv = Polyhedra.volume(ph)
    #
    # io = open(path2*"HV.txt", "a")
    # writedlm(io ,round(hv; digits=2), '\n')
    # close(io)
#
# end


for i= 1:10
    result = readdir(path2)[i]
    mtx = readdlm(path2*result)
    l = Int(length(mtx)/3)
    reshape(mtx, l,3)
    n = 5 # 5 10 15 30 50 70 100
    coeff = 20 # 20 1000 80
    ref = [n*coeff n*coeff n*coeff]
    mtx = vcat(mtx,ref)
    v2 = polyhedron(vrep(mtx), QHull.Library())
    ip_hv = Polyhedra.volume(v2)

    io = open(path2*"HV.txt", "a")
    writedlm(io ,round(ip_hv; digits=2), '\n')
    close(io)
end
