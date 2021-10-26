###########################  General Record mean calulation ######################
using Statistics,DataFrames,DelimitedFiles,CSV

dir1 = "/home/ak121396/Desktop/FPBH/MIPLIP/collection/Y/"
dir1 = "/home/ak121396/Desktop/solvers/Bensolve/MIPLIB//Y/"
ndf = readdir(dir1)

for i=1:length(ndf)
    yval = readdlm(dir1*ndf[i])
    numsol = size(yval)[1]
    print(numsol,"\n")
    # if numsol<=4
    #     print(ndf[i],": ""\n") #,numsol,
    # end
end


readdir("/home/ak121396/Desktop/instances/MIPLIB(official)/LP/")
numsol = size(yval)


# direc = "/home/ak121396/Desktop/GeneralPR/goutputs/FLP/GLPK//"
# direc = "/home/ak121396/Desktop/FPBH/AP/2minTL/time/"
dir1 = "/home/ak121396/Desktop/FPBH/FLP/GLPK/"
# cput = readdir(direc)
sol = readdir(dir1)
# tb = zeros(10,2); tt = zeros(10,2); #AP/KP
tb = zeros(10,2); tt = zeros(12,2); #FLP
# readdlm(direc*cput[1], ',' ,Float64)
# size(readdlm(dir1*sol[1] ,Float64))
1
for i=1:12
    for j=1:10
        k=(i-1)*10+j
        # t = readdlm(direc*cput[k], ',' ,Float64)
        s = readdlm(dir1*sol[k] ,Float64)
        # tb[j,1] = t[1];
        tb[j,2] = size(s)[1];
    end
    # tt[i,1] = round(mean(tb[:,1]),digits=1);
    tt[i,2] = round(mean(tb[:,2]),digits=1)
end
print("solutions \n"); sol=tt[:,2]
CPUtime=tt[:,1]
r = DataFrame(sol=tt[:,2],CPUtime=tt[:,1])

for i=1:7
    f = readdlm(direc*files[i], ',' ,Float64)
    init = round(mean(f[:,1]),digits=2)
    newsol = round(mean(f[:,2]),digits=2)
    totalsol = round(init+newsol, digits=2)
    feasitime = round(mean(f[:,end-2]),digits=2)
    cputime = round(mean(f[:,end]),digits=2)

    record = DataFrame(totalsol=totalsol, newsol=newsol, feasitime=feasitime,CPUtime=cputime)
    CSV.write("/home/ak121396/Documents/GPR.ods",record,append=true,header=false)
end
#FPBH records
cput = readdir(dir2)
dir1 = "/media/ak121396/0526-8445/results/performance/FPBH/ep/MIPLIB/"

function calculate(nins,rep,subclas,dir)
    tb = zeros(nins,1); tt = zeros(subclas,1);
    allt = zeros(subclas,rep)
    for l=1:rep
        pth = dir*"$l/"
        files = readdir(pth)
        for i=1:subclas
            for j=1:nins
                k=(i-1)*nins+j
                s = readdlm(pth*files[k] ,Float64)[1]
                tb[j,1] = s #t[1];
            end
            tt[i,1]=round(mean(tb[:,1]),digits=3)
        end
        allt[:,l]=tt[:,1]
    end
    avgy = [mean(allt[i,:]) for i=1:subclas]
    return round.(avgy,digits=3)
end

avgy = calculate(10,5,10,dir1)

dir1 = "/media/ak121396/0526-8445/results/performance/GPR/ep/MIPLIB/neos/"
dir2 = "/media/ak121396/0526-8445/results/performance/FPBH/ep/MIPLIB/neos/"
tb = zeros(10,2);
files = readdir(dir1)
files2 = readdir(dir2)
for i=1:10
    tb[i,1] = readdlm(dir1*files[i])[1]
    tb[i,2] = readdlm(dir2*files2[i])[1]
end

round(mean(tb[1:5,1]),digits=2)
round(mean(tb[5:10,1]),digits=2)

round(mean(tb[1:5,2]),digits=2)
round(mean(tb[5:10,2]),digits=2)

lendir = readdir(dir1)
tb = zeros(5,1); tb2 = zeros(5,9)
for i=1:9#length(lendir)
    records = readdir(dir1*"/$i/")
    for j=1:length(records)
        c = readdlm(dir1*"/$i/"*records[j])[1]
        tb2[j,i] = c
    end
    # tb2[:,i-1] = tb
end
tb2
my = [mean(tb2[i,1:9]) for i=1:5]
round.(my,digits=3)

for i=1:12
    f = readdlm(dir*f2[i], ',' ,Float64)
    tt = zeros(10,2)
    for i=1:10
        r = f[(i-1)*10+1:i*10,1:2]
        tt[i,1] = round(mean(r[:,1]),digits=2)
        tt[i,2] = round(mean(r[:,2]),digits=2)
        # tt[i,3] = round(mean(r[:,3]), digits=4)
        # tt[i,4] = round(mean(r[:,4]), digits=2)
    end

    record = [tt[:,1:2][2:end,:]; transpose(tt[:,1:2][1,:])]
    avsol[:,i] = record[:,1]
    avtime[:,i] = record[:,2]
end
tb1 = zeros(10,1); tb2 = zeros(10,1);
for i=1:10
    tb1[i,1] = round(mean(avsol[i,:]),digits=2)
    tb2[i,1] = round(mean(avtime[i,:]),digits=4)
end
CSV.write(, df)
tb1
tb2

# ep2 = readdlm("/home/ak121396/Downloads/KirlikSayin2014/ndf/10_020_02.txt.lp.ndf")
# fs = size(ep)[1]
# ep[:,1] = ep2[:,3]; ep[:,3] = ep2[:,1]
# GFP = readdlm("/home/ak121396/Desktop/GFP2hr/Ycubemean/10_020_02.txt_GFP_Y_.csv")
# gs = size(GFP)[1]
#
# epdom=0;gfpdom=0;
# for i=1:gs
#     for j=1:fs
#         if [GFP[i,:][k]> ep[j,:][k] for k=1:3]==[1,1,1]#dominated by PF[j]
#             print("GFP", i," dom by EP ",j,"\n");
#             print("\n==========\n",GFP[i,:],"dom by ", ep[j,:],"\n==========\n")
#             global gfpdom+=1
#             break
#         elseif [ep[j,:][k]> GFP[i,:][k] for k=1:3]==[1,1,1]
#             print("EP",j," dom by GFP ", i,"\n")
#             print("\n==========\n",ep[j,:],"dom by ", GFP[i,:],"\n==========\n")
#             global epdom+=1
#         end
#     end
# end
#
#
# for i=1:fs
#     for j=1:gs
#         if [ep[i,:][k]> GFP[j,:][k] for k=1:3]==[1,1,1]#dominated by PF[j]
#             print("GFP", i," dom by EP ",j,"\n");
#             print("\n==========\n",ep[i,:][1],"dom by ", GFP[j,:],"\n==========\n")
#             global epdom+=1
#             break
#         end
#     end
# end

##########################KP
path = "/home/ak121396/Desktop/FPBH/kp/time/1/"
timere = readdir(path)
avtb = zeros(100,1)
for i=1:100

    avtb[i]  = readdlm(path*timere[i])[1]

end

tbb = zeros(10,1)
for i=1:10
    tbb[i] =  round(mean(avtb[(i-1)*10+1:i*10]), digits=2)
end

path = "/home/ak121396/Desktop/FPBH/kp/10ndf/"
ndf = readdir(path)
avtb = zeros(100,1)
for i=1:100

    f = readdlm(path*ndf[i])
    a,b = size(f)
    avtb[i] = a
end
tbb = zeros(10,1)

for i=1:10
    tbb[i] =  round(mean(avtb[(i-1)*10+1:i*10]), digits=2)
end


##################AP

timere = readdir("/home/ak121396/Desktop/FPBH/ap_lp/time/")
avtb = zeros(100,1)
for i=1:100
    avtb[i]  = readdlm("/home/ak121396/Desktop/FPBH/ap_lp/time/"*timere[i])[1]
end

ndf = readdir("/home/ak121396/Desktop/FPBH/ap_lp/ndf/")
avtb = zeros(100,1)
for i=1:100

    f = readdlm("/home/ak121396/Desktop/FPBH/ap_lp/ndf/"*ndf[i])
    a,b = size(f)
    avtb[i] = a
end

a = [
    1
    2
]
tbb = zeros(10,1)
for i=1:10
    tbb[i] =  round(mean(a[(i-1)*10+1:i*10]), digits=2)
end


####################GFP
timere = readdir("/home/ak121396/Desktop/phd//APresults//")

##########Kirlik

direc = "/home/ak121396/Downloads/KirlikSayin2014/LPfiles/KP/ndf/"
sol = readdir(direc)

avtb = zeros(100,1)
for i=1:100
    tsol = readdlm(direc*sol[i])
    a,b = size(tsol)
    avtb[i] = a
end
