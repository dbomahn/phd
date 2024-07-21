###########################  General Record mean calulation ######################
using Statistics,DataFrames,DelimitedFiles,CSV

# direc = "/home/ak121396/Desktop/GeneralPR/goutputs/FLP/GLPK//"
dir1 = "/media/ak121396/0526-8445/results/gpr/AP/"
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

tb = zeros(15,10)
dir1 = "/home/ak121396/Desktop/relise/epsilon/"
folders = readdir(dir1)
readdir(dir1*folders[1])
readdlm(dir1*folders[1]*"/"*files[1])
for i=1:10
    files = readdir(dir1*folders[i])
    for j=1:15
        @show ep1 = readdlm(dir1*folders[i]*"/"*files[j])
        ep1 = NDfilter2([ep1[i,:] for i=1:size(ep1,1)]);
        tb[j,i] = length(ep1)    
    end
end
round.(mean.(tb[i,:] for i=1:15),digits=2)


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
