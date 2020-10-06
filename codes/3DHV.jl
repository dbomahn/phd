# Run external HV calculator
using DelimitedFiles,DataFrames,Statistics,CSV
cd("/home/ak121396//Downloads/hv-1.3-src")

##################   Find nadir points: the worst value of true PF   #######################
ksdir = "/home/ak121396/Downloads/KirlikSayin2014/LPfiles//KP/ndf/"
fpbh = "/home/ak121396/Desktop/FPBH/kp/7ndf/"
# ben = "/home/ak121396/Desktop/BENoutputs/AP/Y/"
# pr50 = "/home/ak121396/Desktop/BENKP/50Y10/"

# fpbh = "/home/ak121396/Desktop/FPBH/ap_lp/ndf/"
# ksfiles = readdir(ksdir)
# benfiles = readdir(ben)
# pr50files = readdir(pr50)
# pr70files = readdir(pr70)
fpfiles = readdir(fpbh)
1
# function normHV(ksdir,ksfiles,dir,files,i)
#   ksobj = readdlm(ksdir*ksfiles[i])
#   obj = round.(readdlm(dir*files[i]))
#   #KirlikSayin,KP
#   x = obj[:,1]; y=obj[:,2]; z=obj[:,3];
#
#   # AP
#   # obj2 = DataFrame(obj)
#   # obj = obj2[obj2[:x1].!=0,:]
#   # x = obj[:,2]; y=obj[:,3]; z=obj[:,4]; #Bensolve_AP
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

table = zeros(10,10)


for i=1:10
  for j=1:10
    k = j+(i-1)*10
    ithhv = normHV(ksdir,ksfiles,ksdir,ksfiles,k)
    # ithhv = normHV(ksdir,ksfiles,ben,benfiles,k)
    # ithhv = normHV(ksdir,ksfiles,pr50,pr50files,k)
    # ithhv = normHV(ksdir,ksfiles,pr50,pr50files,k)
    table[j,i] = ithhv
  end
end

tt = []
for i=1:10
    a = round(mean(table[:,i]),digits=2)
    push!(tt,a)
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

# CPUtime
readdlm(gfpkirlik*"log/"*gklogs[10])[3,2]
gklogs = readdir(gfpkirlik*"/log/")
gkcpu = 0
for i=1:length(gklogs)
  f = readdlm(gfpkirlik*"/log/"*gklogs[i])
  global gkcpu = gkcpu+f[3,2]
end
#######################GFP+Kirlik: Normalisation and Calculate HV  ########
function GKHV(ksdir,ksfiles,dir,files,i,x,y,z)
  ksobj = readdlm(ksdir*ksfiles[i])
  # obj = readdlm(dir*files[i])
  # x=z; z=x #GFP,Laumanns

  ideal = [minimum(ksobj[:,y]) for y=1:3]
  nadir = [maximum(ksobj[:,y]) for y=1:3]

  r = length(x); normx = [];normy = [];normz = []
  for k=1:r
    push!(normx,(x[k]-ideal[1])/(nadir[1]-ideal[1]))
    push!(normy,(y[k]-ideal[2])/(nadir[2]-ideal[2]))
    push!(normz,(z[k]-ideal[3])/(nadir[3]-ideal[3]))
  end

  dfE = normz,normy,normx;
  Y=DataFrame(dfE);
  CSV.write(dir*files[i][1:9]*"_normal_Y.csv",Y, writeheader=false, delim=' ' )
  cd("/home/ak121396//Downloads/hv-1.3-src")
  @show smetric =readlines( pipeline(`./hv -r "2 2 2" $(dir*files[i][1:9]*"_normal_Y.csv")`))
  return parse(Float64,smetric[1])
end
GKHV(ksdir,ksfiles,gfpkirlik,gkfiles,36,z,y,x)

#######################     Normalisation and Calculate HV   #################


ksobj = readdlm(ksdir*ksfiles[i])
readdlm(fpbh*fpfiles[i])
obj = readdlm(gfp*gfpfiles[i])
x = obj[:,1]; y=obj[:,2]; z=obj[:,3]; #KirlikSayin,GFP

ideal = [minimum(ksobj[:,y]) for y=1:3]
nadir = [maximum(ksobj[:,y]) for y=1:3]

r = length(x); normx = [];normy = [];normz = []
for k=1:r
  push!(normx,(x[k]-ideal[1])/(nadir[1]-ideal[1]))
  push!(normy,(y[k]-ideal[2])/(nadir[2]-ideal[2]))
  push!(normz,(z[k]-ideal[3])/(nadir[3]-ideal[3]))
end

dfE = normz,normy,normx;
Y=DataFrame(dfE);
CSV.write(dir*files[i][1:end-4]*"_normal_Y.csv",Y, writeheader=false, delim=' ' )
cd("/home/ak121396//Downloads/hv-1.3-src")
@show smetric =readlines( pipeline(`./hv -r "2 2 2" $(dir*files[i][1:end-4]*"_normal_Y.csv")`))
return parse(Float64,smetric[1])





############################calculating HV #################################
direc = "/home/ak121396/Downloads/KirlikSayin2014/ndf/"

tb = zeros(10,12)
for j=1:12
  np = ndpt[j,:]
  np1=np[1];np2=np[2];np3=np[3]
  for i=1:10
    k=i+(j-1)*10
    @show smetric =readlines( pipeline(`./hv -r "$np1 $np2 $np3" $(direc*folders[k])`))
    tb[i,j] = parse(Float64,smetric[1])
  end
end


KShv = mean.(tb[:,i] for i=1:12)

##################            Calculating HV        ###########
direc = "/home/ak121396/Desktop/GFP2hr/Y/"
# direc = "/home/ak121396/Desktop/epsilon2hr/Y/"
folders = readdir(direc)
tb = zeros(10,12)
for j=1:12
  np = ndpoint[j,:]
  np1=np[1];np2=np[2];np3=np[3]
  for k=(j-1)*10+1:j*10
    k2 = k-(10*(j-1))
    # output = readdlm(direc*folders[k])
    smetric =readlines( pipeline(`./hv -r "$np1 $np2 $np3" $(direc*folders[k])`))
    tb[k2,j] = parse(Float64,smetric[1])
  end
end
ephv = mean.(tb[:,i] for i=1:12)
 # clipboard(tb)





########################No.solutions
direc = "/home/ak121396/Desktop/GFPKS/20_040/ndfile/"
# direc = "/home/ak121396/Downloads//KirlikSayin2014/ndf/"
files = readdir(direc)
KSnum = []
for i=1:length(files)
  f = readdlm(direc*files[i])
  l,m =size(f)
  push!(KSnum,l)
end
clipboard((sum(KSnum)-length(KSnum))/10)
##average solutions of subclass instances
avg = []
for i=1:12
  push!(avg,mean(launum[10*(i-1)+1:10*i]))
end

direc = "/home/ak121396/Desktop/GFP2hr/Ycubemean//"
#epsilon2hr/Y/"
files = readdir(direc)
launum = []
# f = readdlm(direc*files[1])
for i=1:120
  f = readdlm(direc*files[i])
  l,m =size(f)
  push!(launum,l)
end
clipboard(launum)

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




################## Finding nadir points: the worst values among obtained points of three algorithms ##########
# direc = "/home/ak121396/Desktop/triflp_Y/"
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
