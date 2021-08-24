using DelimitedFiles,DataFrames,StatsBase,CSV # Run external HV calculator

fp = "/home/ak121396/Desktop/FPBH/MIPLIP/GLPK/"
fp = "F:/results/mergedMIP\\"
fpfiles = readdir(fp)
pr = "/home/ak121396/Desktop/GeneralPR/goutputs/MIPLIB/GLPK/2roundY/"
prfiles = readdir(pr)

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
####################################LINUX###################################
# ksdir = "/home/ak121396/Desktop/solvers/Kirlikoutput/FLP/ndf/"
ksdir = "F://results/KS/AP/"
ksfiles = readdir(ksdir)
# fpbh = "/home/ak121396/Desktop/FPBH/FLP/GLPK/"
fpbh = "F:/results/fpbh/AP/1/"
fpfiles = readdir(fpbh)
# pr = "/home/ak121396/Desktop/GeneralPR/goutputs/FLP/GLPK/"
pr = "F:/results/gpr/AP/"
prfiles = readdir(pr)

###################################################
function normHV(ksdir,ksfiles,dir,files,i)
  ksobj = readdlm(ksdir*ksfiles[i])
  # obj = round.(readdlm(dir*files[i]))
  obj = readdlm(dir*files[i])
  # if (0;0;0 in ksobj)
  #   obj = [obj; 0 0 0]
  # end
  x = obj[:,1]; y=obj[:,2]; z=obj[:,3];

  # AP
  # obj2 = DataFrame(obj)
  # obj = obj2[obj2[:x1].!=0,:]
  # x = obj[:,2]; y=obj[:,3]; z=obj[:,4]; #Bensolve_AP
  ideal = [minimum(ksobj[:,y]) for y=1:3]
  nadir = [maximum(ksobj[:,y]) for y=1:3]
  r = length(x); norm = zeros(r,3)
  for k=1:r
      norm[k,1] = (x[k]-ideal[1])/max((nadir[1]-ideal[1]),1)
      norm[k,2] = (y[k]-ideal[2])/max((nadir[2]-ideal[2]),1)
      norm[k,3] = (z[k]-ideal[3])/max((nadir[3]-ideal[3]),1)
  end
  Y=DataFrame(norm) #, :auto);
  CSV.write(dir*files[i][1:end-4]*"_normal_Y.csv",Y, header=false, delim=' ' )
  cd("/home/ak121396//Downloads/hv-1.3-src")
  # cd("C:/cygwin64/home/hv-1.3-src/") #WINDOWS
  # @show
  smetric =readlines( pipeline(`./hv -r "2 2 2" $(dir*files[i][1:end-4]*"_normal_Y.csv")`))
  return parse(Float64,smetric[1])
end

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
table = calculateHV(5,10,12,ksdir,ksfiles,pr)


ksdir = "/media/ak121396/0526-8445/results/mergedMIP/"
ksfiles = readdir(ksdir)
gmip = []
path = "/media/ak121396/0526-8445/results/gpr/MIPLIB/"

fmip = []
path = "/media/ak121396/0526-8445/results/fpbh/MIPLIB/"
for l=1:5
    dir = readdir(path)
    files = readdir(path*dir[l])
    ksfiles = readdir(ksdir*dir[l])
    for j=1:length(files)-1
        hv = normHV(ksdir*dir[l]*"/",ksfiles,path*dir[l]*"/",files,j)
        # push!(gmip,hv)
        push!(fmip,hv)
    end
end

tt = []
for i=1:12
    a = round(mean(table[:,i]),digits=3)
    push!(tt,a)
end
tt


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
