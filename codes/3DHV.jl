using DelimitedFiles,DataFrames,StatsBase,CSV
ENV["CPLEX_STUDIO_BINARIES"] = "/opt/ibm/ILOG/CPLEX_Studio221/cplex/bin/x86-64_linux/"

########################### epsilon constraint methods & matheuristic ############################
############################### Ref set of solver A + solver B   ###################################
# fp = "/media/ak121396/0526-8445/FPBH/MIPLIP/GLPK/"
# fp = "F:/results/mergedMIP\\"
# fpfiles = readdir(fp)
# pr = "/media/ak121396/0526-8445/GeneralPR/goutputs/MIPLIB/GLPK/2roundY/"
# prfiles = readdir(pr)

function biobjHV(ep,epfiles,nd,ndfiles,lsg,lsgfiles,num)
    tb = zeros(num,2);
    # for j=1:folder
    #     ndfiles = readdir(nd)
    for i=1:num
        epobj = readdlm(ep*epfiles[i])
        ndobj = readdlm(nd*ndfiles[i])
        JLD2.@load lsg*lsgfiles[i] lsgdict

        # Line dots into matrix
        ct = 0
        for i=1:length(lsgdict)
            if lsgdict[i]!=[]
                ct = ct + length(lsgdict[i])
            end
        end
        ltb = zeros(ct,2)
        iter = 1
        for i=1:length(lsgdict)
            if lsgdict[i]!=[]
                for j=1:length(lsgdict[i] )
                    ltb[iter,:] = lsgdict[i][j]#[1],lsgdict[i][j][2]
                    iter+=1
                end
            end
        end
        
        x = [epobj[1:10,1]; ndobj[:,1]; ltb[:,1]] 
        y = [epobj[1:10,2]; ndobj[:,2]; ltb[:,2]] 
    
        ideal = [minimum(x),minimum(y)]
        nadir = [maximum(x),maximum(y)]
    
        #  HV calculation
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
    
        # nd HV calculation
        u = size(ndobj)[1]
        norm = zeros(u,2)
        # normalisting
        for k=1:u
            norm[k,1] = (ndobj[:,1][k]-ideal[1])/(nadir[1]-ideal[1])
            norm[k,2] = (ndobj[:,2][k]-ideal[2])/(nadir[2]-ideal[2])
        end
        Y = DataFrame(norm, :auto)
        CSV.write(nd[1:end-3]*"/hv/"*ndfiles[i][1:11]*"_normal_Y.csv",Y, header=false, delim=' ' )
        smetric2 =readlines( pipeline(`./hv -r "2 2" $(nd[1:end-3]*"/hv/"*ndfiles[i][1:11]*"_normal_Y.csv")`))
        # CSV.write(nd*"/hv/"*ndfiles[i][1:11]*"_normal_Y.csv",Y, header=false, delim=' ' )
        # smetric2 =readlines( pipeline(`./hv -r "2 2" $(nd*"/hv/"*ndfiles[i][1:11]*"_normal_Y.csv")`))
        tb[i,2] = parse(Float64,smetric2[1]);
    end
    return tb
end
ep = "/media/ak121396/0526-8445/relise/epsilon/"
epfiles = readdir(ep)[2:end]
nd = "/media/ak121396/0526-8445/relise/lpY/5/"
ndfiles = readdir(mat)

hv5 = biobjHV(ep,epfiles,mat,matfiles,15)

eppr = biobjHV(ep,epfiles,mat,matfiles,10)
1
####################################LINUX###################################
ksdp = "/media/ak121396/0526-8445/solvers/Kirlikoutput/FLP/"
ffp = "/media/ak121396/0526-8445/clu/flpy/"
fpbh = "/media/ak121396/0526-8445/FPBH/FLP/GLPK/"
ffp = "/media/ak121396/0526-8445/LPBMresults/ffp/FLP/"
ben = "/media/ak121396/0526-8445/LPBMresults/benFLP/"
################################Windows #####################################
# ksdir = "F://results/KS/AP/"
# pr = "F:/results/gpr/AP/"
# fpbh = "F:/results/fpbh/AP/1/"

######################### true PF & heuristic ##########################
function normHV(ksdir,ksfiles,dir,files,i)
  ksobj = readdlm(ksdir*ksfiles[i])
  # obj = round.(readdlm(dir*files[i]))
  # if (0;0;0 in ksobj)
  #   obj = [obj; 0 0 0]
  # end

  obj = readdlm(dir*files[i])
#   if (0;0;0 in ksobj)
#     obj = [obj; 0 0 0]
#   end
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
  savedir = "/home/ak121396/Desktop/forHV/"
  CSV.write(savedir*files[i][1:end-4]*"_normal_Y.csv",Y, header=false, delim=' ' )
  cd("/home/ak121396//Downloads/hv-1.3-src")
  # cd("C:/cygwin64/home/hv-1.3-src/") #WINDOWS
  # @show
  smetric =readlines( pipeline(`./hv -r "2 2 2" $(savedir*files[i][1:end-4]*"_normal_Y.csv")`))
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
        # hv = normHV(ksp,ksfiles,ben,bfiles,k)
        hv = normHV(ksp,ksfiles,pr,prfiles,k)
        table[i,j] = hv
    end
    # tb = round(mean(table[:,1]),digits=3)
end
table
round.([mean(table[i,1:5]) for i=1:12],digits=4)
########################## obj file converter: merging FPBH&GPR ################
mipdir = "/media/ak121396/0526-8445/LPBMresults/mergedMIP/"
fdir = "/media/ak121396/0526-8445/LPBMresults/fpbh/MIPLIB/"
gdir = "/media/ak121396/0526-8445/LPBMresults/gpr/MIPLIB/"
readdir(fdir)
for i=1:14
    gname = readdir(gdir*"/$i/"); fname = readdir(fdir*"/$i/")
    for j=1:5
        F = readdlm(fdir*"/$i/"*fname[j])
        G = readdlm(gdir*"/$i/"*gname[j])
        P = vcat(F,G)
        ins = fname[j][1:end-4]
        CSV.write("/media/ak121396/0526-8445/LPBMresults/mergedMIP/$i"*"/$ins"*".txt",DataFrame(P, :auto),header=false, delim=' ' )
    end
end
##############################################################
filedir = "/media/ak121396/0526-8445/LPBMresults/mergedMIP/"
gdir = "/media/ak121396/0526-8445/LPBMresults/gpr/MIPLIB/"
fdir = "/media/ak121396/0526-8445/LPBMresults/fpbh/MIPLIB/"

table = zeros(5,14)
tb = zeros(5,9)
for k=1:9
    mipfiles = readdir(mipdir*"$k")
    files = readdir(gdir*"$k")
    # files = readdir(fdir*"$k")
    for i=1:5#length(mipfiles)
        hv = normHV(mipdir*"$k/",mipfiles,gdir*"$k/",files,i)
        # hv = normHV(mipdir*"$k/",mipfiles,fdir*"$k/",files,i)
        # table[i,k] = hv
        tb[i,k] = hv
    end
    # tb = round(mean(table[:,1]),digits=3)
end
mean(table[4,:]),mean(table[5,:])

round.([mean(table[i,1:9]) for i=4:5],digits=3)
round.([mean(table[i,1:14]) for i=4:5],digits=3)

mean(tb[4,10:14]),mean(tb[5,10:14])
mean(tb[4,:]),mean(tb[5,:])
mean.([round.([mean(tb[i,10:14]) for i=4:5],digits=3),round.([mean(tb[i,1:5]) for i=4:5],digits=3)])

round.([mean(tb[i,10:14]) for i=4:5],digits=3)
round.([mean(tb[i,1:5]) for i=4:5],digits=3)

round.([mean(table[i,:]) for i=1:5],digits=3)

ksp = "/media/ak121396/0526-8445/solvers/Kirlikoutput/FLP/ndf/"
ffp = "/media/ak121396/0526-8445/clu/KP/"
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
ksp = "/media//ak121396/0526-8445/LPBMresults/KS/FLP/"
fpbh = "/media/ak121396/0526-8445//LPBMresults/fpbh/FLP/"
ben = "/media/ak121396/0526-8445//LPBMresults/bensolve/FLP/"
pr = "/media/ak121396/0526-8445//LPBMresults/ffp/FLP/"
ksfiles = readdir(ksp)
fpfiles = readdir(fpbh)
bfiles = readdir(ben)
prfiles = readdir(pr)

table = calculateHV(5,10,12,ksp,ksfiles,pr)
round.([mean(table[i,:]) for i=1:12],digits=4)

refdir = "/media/ak121396/0526-8445/LPBMresults/KS/AP/"
refiles = readdir(refdir)
ag = []
agpath = "/media/ak121396/0526-8445/LPBMresults/gpr/AP/"
af = []
afpath = "/media/ak121396/0526-8445/LPBMresults/fpbh/AP/"
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

refdir = "/media/ak121396/0526-8445/LPBMresults/KS/KP/"
kg = []
kgpath = "/media/ak121396/0526-8445/LPBMresults/gpr/KP/"
kf = []
kfpath = "/media/ak121396/0526-8445/LPBMresults//fpbh/KP/"

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

refdir = "/media/ak121396/0526-8445/LPBMresults/KS/FLP/"
fg = []
fgpath = "/media/ak121396/0526-8445/LPBMresults/gpr/FLP/"
ff = []
ffpath = "/media/ak121396/0526-8445/LPBMresults/fpbh/FLP/"
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


refdir = "/home/ak121396/Desktop/LPBMresults/mergedMIP/"
mg = []
mgpath = "/home/ak121396/Desktop/LPBMresults/gpr/MIPLIB/"
mf = []
mfpath = "/home/ak121396/Desktop/LPBMresults/fpbh/MIPLIB/"
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
# direc = "/media/ak121396/0526-8445/epsilon2hr/Y/"
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
# direc = "/media/ak121396/0526-8445//"
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
