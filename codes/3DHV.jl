using DelimitedFiles,DataFrames,StatsBase,CSV
ENV["CPLEX_STUDIO_BINARIES"] = "/opt/ibm/ILOG/CPLEX_Studio221/cplex/bin/x86-64_linux/"

########################### epsilon constraint methods & matheuristic ############################
############################### Ref set of solver A + solver B   ###################################
# fp = "/home/ak121396/Desktop/FPBH/MIPLIP/GLPK/"
# fp = "F:/results/mergedMIP\\"
# fpfiles = readdir(fp)
# pr = "/home/ak121396/Desktop/GeneralPR/goutputs/MIPLIB/GLPK/2roundY/"
# prfiles = readdir(pr)

function biobjHV(ep,epfiles,mat,matfiles,num)
    tb = zeros(num,2);
    # for j=1:folder
    #     matfiles = readdir(mat)
    for i=1:num
        epobj = readdlm(ep*epfiles[i])
        matobj = readdlm(mat*matfiles[i]);
        x = [epobj[1:10,1]; matobj[:,1]] 
        y = [epobj[1:10,2]; matobj[:,2]] 
    
        ideal = [minimum(x),minimum(y)]
        nadir = [maximum(x),maximum(y)]
    
        # FPBH HV calculation
        r = size(epobj)[1]
        norm = zeros(r,2)
        for k=1:r
            norm[k,1] = (epobj[:,1][k]-ideal[1])/(nadir[1]-ideal[1])
            norm[k,2] = (epobj[:,2][k]-ideal[2])/(nadir[2]-ideal[2])
        end
        # dfE = normz,normy
        Y=DataFrame(norm, :auto);
        CSV.write(ep*"/hv/"*epfiles[i][1:end-8]*"normal_Y.csv",Y, header=false, delim=' ' )
        cd("/home/ak121396/Downloads/hv-1.3-src")
        smetric1 =readlines( pipeline(`./hv -r "2 2" $(ep*"/hv/"*epfiles[i][1:end-8]*"normal_Y.csv")`))
        tb[i,1] = parse(Float64,smetric1[1]);
    
        # mat HV calculation
        u = size(matobj)[1]
        norm = zeros(u,2)
        # normalisting
        for k=1:u
            norm[k,1] = (matobj[:,1][k]-ideal[1])/(nadir[1]-ideal[1])
            norm[k,2] = (matobj[:,2][k]-ideal[2])/(nadir[2]-ideal[2])
        end
        Y = DataFrame(norm, :auto)
        # CSV.write(mat[1:end-3]*"/hv/"*matfiles[i][1:11]*"_normal_Y.csv",Y, header=false, delim=' ' )
        # smetric2 =readlines( pipeline(`./hv -r "2 2" $(mat[1:end-3]*"/hv/"*matfiles[i][1:11]*"_normal_Y.csv")`))
        CSV.write(mat*"/hv/"*matfiles[i][1:11]*"_normal_Y.csv",Y, header=false, delim=' ' )
        smetric2 =readlines( pipeline(`./hv -r "2 2" $(mat*"/hv/"*matfiles[i][1:11]*"_normal_Y.csv")`))
        tb[i,2] = parse(Float64,smetric2[1]);
    end
    return tb
end
ep = "/home/ak121396/Desktop/relise/epsilon/"
epfiles = readdir(ep)[6:end]
mat = "/home/ak121396/Desktop/relise/lpY/"
matfiles = readdir(mat)[end-3:end]

hv2 = biobjHV(ep,epfiles,mat,matfiles,1)
ep7pr = biobjHV(ep,epfiles,mat,matfiles,1)
ep2pr = biobjHV(ep,epfiles,mat,matfiles,13)

eppr = biobjHV(ep,epfiles,mat,matfiles,10)
1
####################################LINUX###################################
ksp = "/home/ak121396/Desktop/solvers/Kirlikoutput/AP&KP/intKP_Y/"
ffp = "/home/ak121396/Desktop/clu/KP/"

ksdp = "/home/ak121396/Desktop/solvers/Kirlikoutput/FLP/ndf/"
ffp = "/home/ak121396/Desktop/clu/flpy/"

fpbh = "/home/ak121396/Desktop/FPBH/FLP/GLPK/"
pr = "/home/ak121396/Desktop/GeneralPR/goutputs/FLP/GLPK/"
ben = "/home/ak121396/Desktop/solvers/Bensolve/APoutputs/Y/"
################################Windows #####################################
# ksdir = "F://results/KS/AP/"
# pr = "F:/results/gpr/AP/"
# fpbh = "F:/results/fpbh/AP/1/"

fpfiles = readdir(fpbh)
prfiles = readdir(pr)
bfiles = readdir(ben)

######################### true PF & heuristic ##########################
function normHV(ksdir,ksfiles,dir,files,i)
  ksobj = readdlm(ksdir*ksfiles[i])
  # obj = round.(readdlm(dir*files[i]))
  # if (0;0;0 in ksobj)
  #   obj = [obj; 0 0 0]
  # end

  obj = readdlm(dir*files[i])
  x = obj[:,1]; y=obj[:,2]; z=obj[:,3];

  # Bensolve AP
  # obj2 = DataFrame(obj)
  # obj = obj2[obj2[:x1].!=0,:]
  # x = obj[:,2]; y=obj[:,3]; z=obj[:,4];

  ideal = [minimum(ksobj[:,y]) for y=1:3]
  nadir = [maximum(ksobj[:,y]) for y=1:3]
  r = length(x); norm = zeros(r,3)
  for k=1:r
      norm[k,1] = (x[k]-ideal[1])/max((nadir[1]-ideal[1]),1)
      norm[k,2] = (y[k]-ideal[2])/max((nadir[2]-ideal[2]),1)
      norm[k,3] = (z[k]-ideal[3])/max((nadir[3]-ideal[3]),1)
  end
  Y=DataFrame(norm, :auto);
  CSV.write(dir*files[i][1:end-4]*"_normal_Y.csv",Y, header=false, delim=' ' )
  cd("/home/ak121396//Downloads/hv-1.3-src")
  # cd("C:/cygwin64/home/hv-1.3-src/") #WINDOWS
  # @show
  smetric =readlines( pipeline(`./hv -r "2 2 2" $(dir*files[i][1:end-4]*"_normal_Y.csv")`))
  return parse(Float64,smetric[1])
end

tb = zeros(10,10)
for i=1:10
    for j=1:10
        k = (i-1)*10+j
        hv = normHV(ksp,kfiles,fpp,ffiles,k)
        tb[i,j] = hv
    end
    # tb = round(mean(table[:,1]),digits=3)
end
tb
round.([mean(tb[i,:]) for i=1:10],digits=3)

table = zeros(12,10)
for i=1:12
    for j=1:10
        k = (i-1)*10+j
        hv = normHV(ksp,ksfiles,fpp,ffiles,k)
        table[i,j] = hv
    end
    # tb = round(mean(table[:,1]),digits=3)
end
table
round.([mean(table[i,:]) for i=1:12],digits=3)

###############################################################################
mipdir = "F:/results/mergedMIP/";
mipdir = "/media/ak121396/0526-8445/results/mergedMIP/"
mipfiles = readdir(mipdir)

# filedir = "/media/ak121396/0526-8445/results/gpr/MIPLIB/"
filedir = "/media/ak121396/0526-8445/results/fpbh/MIPLIB/"
table = zeros(5,5)
for k=1:5
    mipfiles = readdir(mipdir*"$k")
    files = readdir(filedir*"$k")
    for i=1:5#length(mipfiles)
        hv = normHV(mipdir*"$k/",mipfiles,filedir*"$k/",files,i)
        table[i,k] = hv
    end
    # tb = round(mean(table[:,1]),digits=3)
end
round.([mean(table[i,:]) for i=1:5],digits=3)

ksp = "/home/ak121396/Desktop/solvers/Kirlikoutput/FLP/ndf/"
ffp = "/home/ak121396/Desktop/clu/KP/"
readdir(ffp)
kfiles = readdir(ksp)

function calculateHV(rep,ins,subclas,ksdir,ksfiles,path)
    allt = zeros(subclas,rep)
    for l=1:rep
        dir = readdir(path)
        files = readdir(path*dir[l])
        tb = zeros(ins,1); tt = zeros(subclas,1);
        for i=1:subclas
          for j=1:ins
            k = j+(i-1)*ins
            hv = normHV(ksdir,ksfiles,path*dir[l]*"/",files,k)
            tb[j,1] = hv
          end
          tt[i,1]=round(mean(tb[:,1]),digits=4)
        end
        allt[:,l]=tt[:,1]
    end
    return(allt)
end
table = calculateHV(5,10,10,ksp,kfiles,ffp)
round.([mean(table[i,:]) for i=1:10],digits=3)

refdir = "/media/ak121396/0526-8445/results/KS/AP/"
refiles = readdir(refdir)
ag = []
agpath = "/home/ak121396/Desktop/forHV/gpr/AP/"
af = []
afpath = "/home/ak121396/Desktop/forHV/fpbh/AP/"
for l=1:5
    refiles = readdir(refdir)
    gfiles = readdir(agpath*"$l/")
    ffiles = readdir(afpath*"$l/")

    for j=1:length(refiles)
        ghv = normHV(refdir,refiles,agpath*"$l/",gfiles,j)
        fhv = normHV(refdir,refiles,afpath*"$l/",ffiles,j)
        push!(ag,ghv); push!(af,fhv)
    end
end

refdir = "/media/ak121396/0526-8445/results/KS/KP/"
kg = []
kgpath = "/home/ak121396/Desktop/forHV/gpr/KP/"
kf = []
kfpath = "/home/ak121396/Desktop/forHV//fpbh/KP/"

for l=1:5
    refiles = readdir(refdir)
    gfiles = readdir(kgpath*"$l/")
    ffiles = readdir(kfpath*"$l/")

    for j=1:length(refiles)
        ghv = normHV(refdir,refiles,kgpath*"$l/",gfiles,j)
        fhv = normHV(refdir,refiles,kfpath*"$l/",ffiles,j)
        push!(kg,ghv); push!(kf,fhv)
    end
end

refdir = "/media/ak121396/0526-8445/results/KS/FLP/"
fg = []
fgpath = "/home/ak121396/Desktop/forHV/gpr/FLP/"
ff = []
ffpath = "/home/ak121396/Desktop/forHV/fpbh//FLP/"

for l=1:5
    refiles = readdir(refdir)
    gfiles = readdir(fgpath*"$l/")
    ffiles = readdir(ffpath*"$l/")

    for j=1:length(refiles)
        ghv = normHV(refdir,refiles,fgpath*"$l/",gfiles,j)
        fhv = normHV(refdir,refiles,ffpath*"$l/",ffiles,j)
        push!(fg,ghv); push!(ff,fhv)
    end
end


refdir = "/media/ak121396/0526-8445/results/mergedMIP/"
mg = []
mgpath = "/home/ak121396/Desktop/forHV/gpr/MIPLIB/"
mf = []
mfpath = "/home/ak121396/Desktop/forHV/fpbh/MIPLIB/"
for l=1:9
    refiles = readdir(refdir*"$l/")
    ffiles = readdir(mfpath*"$l/")
    gfiles = readdir(mgpath*"$l/")
    for j=1:length(refiles)
        ghv = normHV(refdir*"$l/",refiles,mgpath*"$l/",gfiles,j)
        fhv = normHV(refdir*"$l/",refiles,mfpath*"$l/",ffiles,j)
        push!(mg,ghv); push!(mf,fhv)
    end
end

mergedhv = reshape(miphv,(5,10))
round.([mean(mergedhv[i,1:9]) for i=1:5],digits=3)


g1 = reshape(gmip,(5,9))
round.([mean(g1[i,1:9]) for i=1:5],digits=3)
f1 = reshape(fmip,(5,9))
round.([mean(f1[i,1:9]) for i=1:5],digits=3)


tt = []
for i=1:12
    a = round(mean(table[:,i]),digits=3)
    push!(tt,a)
end
tt
function fpprHV(fp,fpfiles,pr,prfiles,num) #,pr2,prfiles2,pr3,prfiles3,ks,ksfiles,num)
    tb = zeros(num,2);
    for i=1:num
        fpobj = readdlm(fp*fpfiles[i])
        probj = readdlm(pr*prfiles[i]);
        # probj2 = readdlm(pr2*prfiles2[i]);probj3 = readdlm(pr3*prfiles3[i]);
        # ksobj = readdlm(ks*(ksfiles[i]))
        if (0;0;0 in fpobj)
          probj = [probj;0 0 0]
        end
        x = [fpobj[:,1]; probj[:,1]] #;probj2[:,1];probj3[:,1]; ksobj[:,1]]
        y = [fpobj[:,2]; probj[:,2]] #;probj2[:,1];probj3[:,1];ksobj[:,2]]
        z = [fpobj[:,3]; probj[:,3]] #;probj2[:,1];probj3[:,1];ksobj[:,3]]

        ideal = [minimum(x),minimum(y),minimum(z)]
        nadir = [maximum(x),maximum(y),maximum(z)]

        # FPBH HV calculation
        r = size(fpobj)[1]
        norm = zeros(r,3)
        for k=1:r
            norm[k,1] = (fpobj[:,1][k]-ideal[1])/(nadir[1]-ideal[1])
            norm[k,2] = (fpobj[:,2][k]-ideal[2])/(nadir[2]-ideal[2])
            norm[k,3] = (fpobj[:,3][k]-ideal[3])/(nadir[3]-ideal[3])
        end
        # dfE = normz,normy,normx;
        Y=DataFrame(norm, :auto);
        CSV.write(fp*fpfiles[i][1:end-7]*"normal_Y.csv",Y, header=false, delim=' ' )
        cd("/home/ak121396/Downloads/hv-1.3-src")
        smetric1 =readlines( pipeline(`./hv -r "2 2 2" $(fp*fpfiles[i][1:end-7]*"normal_Y.csv")`))
        tb[i,1] = parse(Float64,smetric1[1]);

        # GPR HV calculation
        u = size(probj)[1]
        norm = zeros(u,3)
        for k=1:u
            norm[k,1] = (probj[:,1][k]-ideal[1])/(nadir[1]-ideal[1])
            norm[k,2] = (probj[:,2][k]-ideal[2])/(nadir[2]-ideal[2])
            norm[k,3] = (probj[:,3][k]-ideal[3])/(nadir[3]-ideal[3])
        end
        Y = DataFrame(norm, :auto)
        CSV.write(pr*prfiles[i][1:end-6]*"_normal_Y.csv",Y, header=false, delim=' ' )
        smetric2 =readlines( pipeline(`./hv -r "2 2 2" $(pr*prfiles[i][1:end-6]*"_normal_Y.csv")`))
        tb[i,2] = parse(Float64,smetric2[1]);
    end
    return tb
end



########################  Merging GFP+Kirlik output  ######################
#NDpoint
x=[0];y=[0];z=[0]
for i=1:length(gkfiles)
  f = readdlm(gfpkirlik*"hybndf/"*gkfiles[i])
  append!(x,f[:,1]);append!(y,f[:,2]);append!(z,f[:,3])
  # for i=1:length(f)
  #   append!(obj,f)
end
###############      Adding origin[0,0,0] to ndpoints     ########################
# direc = "/home/ak121396/Desktop/epsilon2hr/Y/"
folders = readdir(direc)[3:13]
points = readdlm(direc*folders[2]*"/"*files[12])
for j in folders
  files = readdir(direc*"/"*j)

  for i in files
    ndf =[0 0 0]
    org = ndf[:,1],ndf[:,2],ndf[:,3]
    df = DataFrame(org)
    CSV.write(direc*j*"/"*i, df; append=true, writeheader=false,delim=' ')
  end
end
########################windows ##################

################## Finding nadir points: the worst values among obtained points of three algorithms ##########
# direc = "/home/ak121396/Desktop//"
# folders = readdir(direc)[2:13]
# ndpoint = Dict()
# # readdir(direc*"09")
# for j in folders[1:12]
#   output = readdir(direc*j)
#   nadirpt = []
#   for k=1:20
#     @show obj = readdlm(direc*j*"/"*output[k], ' ')
#     x=obj[:,1];y=obj[:,2];z=obj[:,3]
#     nadir = [maximum(obj[:,k]) for k=1:3]
#     push!(nadirpt,nadir)
#   end
#   uninp = [maximum(hcat(nadirpt...)[i,:]) for i=1:3]+[1,1,1]
#   np1 = uninp[1]; np2=uninp[2]; np3=uninp[3]
#   ndpoint[j] = [np1,np2,np3]
# end
# function Bensolve_normHV(ksdir,ksfiles,dir,files,i) #negative sign
#   ksobj = readdlm(ksdir*ksfiles[i])
#   obj = readdlm(dir*files[i])
#   idx = findall(x->x==1, obj[:,1])
#   #KP
#   x = -[obj[l,2] for l in idx]
#   y = -[obj[l,3] for l in idx]
#   z = -[obj[l,4] for l in idx]
#
#   ideal = [minimum(ksobj[:,y]) for y=1:3]
#   nadir = [maximum(ksobj[:,y]) for y=1:3]
#
#   r = length(x); normx = [];normy = [];normz = []
#   for k=1:r
#     push!(normx,(x[k]-ideal[1])/(nadir[1]-ideal[1]))
#     push!(normy,(y[k]-ideal[2])/(nadir[2]-ideal[2]))
#     push!(normz,(z[k]-ideal[3])/(nadir[3]-ideal[3]))
#   end
#
#   dfE = normz,normy,normx;
#   Y=DataFrame(dfE);
#   CSV.write(dir*files[i][1:end-4]*"_normal_Y.csv",Y, header=false, delim=' ' )
#   cd("/home/ak121396//Downloads/hv-1.3-src")
#   @show smetric =readlines( pipeline(`./hv -r "2 2 2" $(dir*files[i][1:end-4]*"_normal_Y.csv")`))
#   return parse(Float64,smetric[1])
# end
