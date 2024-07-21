using DelimitedFiles,DataFrames,StatsBase,CSV
ENV["CPLEX_STUDIO_BINARIES"] = "/opt/ibm/ILOG/CPLEX_Studio221/cplex/bin/x86-64_linux/"

########################### epsilon constraint methods & matheuristic ############################
############################### Ref set of solver A + solver B   ###################################
# fp = "/media/ak121396/0526-8445/FPBH/MIPLIP/GLPK/"
# fp = "F:/results/mergedMIP\\"
# fpfiles = readdir(fp)
# pr = "/media/ak121396/0526-8445/GeneralPR/goutputs/MIPLIB/GLPK/2roundY/"
# prfiles = readdir(pr)
########################      Linesegment's Hyper Volumn      ###########################
function LinesegHV(ndy,Y) #lsgpath,num,Y
    larms = findall(i -> i == "L", ndy.arm)
    segsets = []
    start = 0; lsghv = 0;
    for l in larms
        if 1+start != l
            push!(segsets, [1+start, l])
        end
        start = l
    end
    for pair in segsets
        for i=pair[1]:pair[2]-1
            lsghv  = lsghv + (abs(Y[i,1] - Y[i+1,1])*abs(Y[i,2] - Y[i+1,2]))/2
        end     
    end
    return lsghv
end
# lsg = readdir(lsgpath)
# ndy = CSV.read(lsgpath*lsg[num],DataFrame)
# (abs(ndy.v1[i][1] - ndy.v1[i+1][1])*abs(ndy.v1[i][2] - ndy.v1[i+1][2]))/2
function normHV2(optdir,dir,nfolder,nfile,nobj)
    tb = zeros(nfile,nfolder)
    folders = readdir(dir)  
    for j=1:nfolder
        files = readdir(dir*"$j")#folders[j]
        for i=1:nfile
            optobj = readdlm(optdir*"/$i.csv")
            ##############   nondominated points & linesegments    ###############
            # merged points 
            obj = readdlm(optdir*"/$i.csv")
            #epsilon 
            # obj = readdlm(dir*"$j"*"/"*files[i])
            #matheuristic
            nodes = CSV.read(dir*folders[j]*"/"*files[i],DataFrame); obj = hcat(nodes.v1,nodes.v2)
            if nobj == 2
                x = obj[:,1]; y=obj[:,2];
            elseif nobj == 3
                x = obj[:,1]; y=obj[:,2]; z=obj[:,3];
            end
        
            ideal = [minimum(optobj[:,y]) for y=1:nobj]
            nadir = [maximum(optobj[:,y]) for y=1:nobj]
            r = length(x); norm = zeros(r,nobj)
            for k=1:r
                norm[k,1] = (x[k]-ideal[1])/max((nadir[1]-ideal[1]),1)
                norm[k,2] = (y[k]-ideal[2])/max((nadir[2]-ideal[2]),1)
            #   norm[k,3] = (z[k]-ideal[3])/max((nadir[3]-ideal[3]),1)
            end
            Y=DataFrame(norm, :auto);
            # savedir = "/home/ak121396/Desktop/forHV/"*folders[j]*"/"
            # CSV.write(savedir*files[i][1:end-4]*"_normal_Y.csv",Y, header=false, delim=' ' )
            #merged points
            savedir = "/home/ak121396/Desktop/forHV/$i"
            CSV.write(savedir*"$i"*"_normal_Y.csv",Y, header=false, delim=' ' )
            cd("/home/ak121396//Downloads/hv-1.3-src")
            if nobj == 2
                # smetric =readlines( pipeline(`./hv -r "2 2" $(savedir*files[i][1:end-4]*"_normal_Y.csv")`))
                smetric =readlines( pipeline(`./hv -r "2 2" $(savedir*"$i"*"_normal_Y.csv")`))
            else
                smetric =readlines( pipeline(`./hv -r "2 2 2" $(savedir*files[i][1:end-4]*"_normal_Y.csv")`))
            end
            tb[i,j] = parse(Float64,smetric[1])
            lsghv = LinesegHV(nodes,Y) #lsgpath,num,Y
            tb[i,j] = tb[i,j] + lsghv
        end
    end
    return tb
end
ta = normHV2(od,dd,10,15,2)
for i=1:15
    println(mean(ta[i,:]))
end
od = "/home/ak121396/Desktop/relise/performance/"
dd = "/home/ak121396/Desktop/relise/vopt/nodes/"
ed = "/home/ak121396/Desktop/relise/epsilon/"

function biobjHV(ep,epfiles,mh,mhfiles,hvpath,num)
    tb = zeros(num,2);
    # for j=1:folder
    #     mhfiles = readdir(mh)
    for i=1:num
        epobj = readdlm(ep*epfiles[i])
        # mhobj = readdlm(mh*mhfiles[i]) # matheuristic nondominated points
        
        ##############   nondominated points & linesegments    ###############
        nodes = CSV.read(mh*mhfiles[i],DataFrame) 
        mhobj = hcat(nodes.v1,nodes.v2)
        ######################################################################

        # x = [epobj[1:10,1]; mhobj[:,1]] 
        # y = [epobj[1:10,2]; mhobj[:,2]] 
    
        # ideal = [minimum(x),minimum(y)]
        # nadir = [maximum(x),maximum(y)]
    
        #  HV calculation
        r = size(epobj)[1]
        norm = zeros(r,2)
        for k=1:r
            norm[k,1] = (epobj[:,1][k]-ideal[1])/(nadir[1]-ideal[1])
            norm[k,2] = (epobj[:,2][k]-ideal[2])/(nadir[2]-ideal[2])
        end
        # dfE = normz,normy
        Y=DataFrame(norm, :auto);
        CSV.write(ep[1:end-3]*"/hv/"*epfiles[i][1:end-8]*"normal_Y.csv",Y, header=false, delim=' ' )
        cd("/home/ak121396/Downloads/hv-1.3-src")
        smetric1 =readlines( pipeline(`./hv -r "2 2" $(ep[1:end-3]*"/hv/"*epfiles[i][1:end-8]*"normal_Y.csv")`))
        tb[i,1] = parse(Float64,smetric1[1]);
    
        # mh nd HV calculation
        u = size(mhobj)[1]
        norm = zeros(u,2)
        # normalisting
        for k=1:u
            norm[k,1] = (mhobj[:,1][k]-ideal[1])/(nadir[1]-ideal[1])
            norm[k,2] = (mhobj[:,2][k]-ideal[2])/(nadir[2]-ideal[2])
        end
        Y = DataFrame(norm, :auto)
        
        CSV.write(hvpath*mhfiles[i][1:11]*"_normal_Y.csv",Y, header=false, delim=' ' )
        smetric2 =readlines( pipeline(`./hv -r "2 2" $(hvpath*mhfiles[i][1:11]*"_normal_Y.csv")`))
        # CSV.write(mh*"/hv/"*mhfiles[i][1:11]*"_normal_Y.csv",Y, header=false, delim=' ' )
        # smetric2 =readlines( pipeline(`./hv -r "2 2" $(mh*"/hv/"*mhfiles[i][1:11]*"_normal_Y.csv")`))
        tb[i,2] = parse(Float64,smetric2[1]);
        ############################       Add Linesegment HV       #############################
        lsghv = LinesegHV(nodes,Y) #lsgpath,num,Y
        tb[i,2] = tb[i,2] + lsghv
    end
    return tb
end
function MergeNDPoints(epath,mpath,nfolder,filenum)
    obj1 = []; obj2 = []
    for i=1:nfolder
        efiles = readdir(epath*"$i"); 
        mfiles = readdir(mpath*"$i")
        eobj = readdlm(epath*"$i/"*efiles[filenum]) 
        mdt = CSV.read(mpath*"$i/"*mfiles[filenum], DataFrame; header=[:v1,:v2]); mobj = hcat(mdt.v1,mdt.v2)
        append!(obj1, eobj[:,1],mobj[:,1]); append!(obj2, eobj[:,2],mobj[:,2])
    end
    apo = hcat(obj1,obj2)        
    merged = NDfilter2([apo[i,:] for i=1:size(apo,1)])
    mgobj = zeros(length(merged),2)
    for i=1:size(mgobj,1)
        mgobj[i,:] = merged[i]
    end
    CSV.write("/home/ak121396/Desktop/relise/performance/$filenum.csv", DataFrame(mgobj, :auto), append=false, header=false, delim=' ')
end
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
for f=1:15
    MergeNDPoints(ep,mat,10,f)
end

ep = "/home/ak121396/Desktop/relise/epsilon/"
epfiles = readdir(ep)#[1:end]
mat = "/home/ak121396/Desktop/relise/vopt/Y/dichow/"
matfiles = readdir(mat)[1:10]
# hv1 = biobjHV(ep,epfiles,mat,matfiles,15)
hvpath = "/home/ak121396/Desktop/relise/vopt/Y/hv/8/"
lsgpath = "/home/ak121396/Desktop/relise/vopt/nodes/8/"
lsgfiles = readdir(lsgpath)
hv2 = biobjHV(ep,epfiles,lsgpath,lsgfiles,hvpath,15)
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
function normHV(optdir,optfiles,dir,files,i)
  optobj = readdlm(optdir*optfiles[i])
  # obj = round.(readdlm(dir*files[i]))
  # if (0;0;0 in optobj)
  #   obj = [obj; 0 0 0]
  # end

#   obj = readdlm(dir*files[i])
#   if (0;0;0 in optobj)
#     obj = [obj; 0 0 0]
#   end
    ##############   nondominated points & linesegments    ###############
    nodes = CSV.read(dir*files[i],DataFrame) 
    obj = hcat(nodes.v1,nodes.v2)


  if nobj == 2
    x = obj[:,1]; y=obj[:,2];
  elseif nobj == 3
    x = obj[:,1]; y=obj[:,2]; z=obj[:,3];
  end

  # Bensolve AP
  # obj2 = DataFrame(obj)
  # obj = obj2[obj2[:x1].!=0,:]
  # x = obj[:,2]; y=obj[:,3]; z=obj[:,4];

  ideal = [minimum(optobj[:,y]) for y=1:nobj]
  nadir = [maximum(optobj[:,y]) for y=1:nobj]
  r = length(x); norm = zeros(r,nobj)
  for k=1:r
      norm[k,1] = (x[k]-ideal[1])/max((nadir[1]-ideal[1]),1)
      norm[k,2] = (y[k]-ideal[2])/max((nadir[2]-ideal[2]),1)
    #   norm[k,3] = (z[k]-ideal[3])/max((nadir[3]-ideal[3]),1)
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

tb = zeros(15,10)
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

ksp = "/media/ak121396/0526-8445/solvers/Kirlikoutput/FLP/ndf/"
ffp = "/media/ak121396/0526-8445/clu/KP/"
readdir(ffp)
kfiles = readdir(ksp)

function calculateHV(rep,ins,subclas,optdir,optfiles,path)
    allt = zeros(subclas,rep)
    for l=1:rep
        dir = readdir(path)
        files = readdir(path*dir[l])
        tb = zeros(ins,1); tt = zeros(subclas,1);
        for i=1:subclas
          for j=1:ins
            k = j+(i-1)*ins
            hv = normHV(optdir,optfiles,path*dir[l]*"/",files,k)
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
