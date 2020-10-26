###########################  General Record mean calulation ######################
using Statistics,DataFrames,DelimitedFiles
# direc = "C:\\Users\\AK121396\\Desktop\\iteratio_KP_record.csv"
direc = "/home/ak121396/Desktop/PR_KP/diff/record/"
files = readdir(direc)
avsol = zeros(10,10);avtime = zeros(10,10);
for i=1:5
    f = readdlm(direc*files[i], ',' ,Float64)
    tt = zeros(10,4)
    for i=1:10
        r = f[(i-1)*10+1:i*10,1:4]
        # tt[i,1] = round(mean(r[:,1]),digits=2)
        tt[i,2] = round(mean(r[:,2]),digits=2)
        tt[i,3] = round(mean(r[:,3]), digits=4)
        tt[i,4] = round(mean(r[:,4]), digits=2)
    end

    record = [tt[:,2:3][2:end,:]; transpose(tt[:,2:3][1,:])]
    avsol[:,i] = record[:,1]
    avtime[:,i] = record[:,2]
end
tb1 = zeros(10,1); tb2 = zeros(10,1);
for i=1:10
    tb1[i,1] = round(mean(avsol[i,:]),digits=2)
    tb2[i,1] = round(mean(avtime[i,:]),digits=4)
end
tb1
tb2

#######################################################

##########################################################################
direc = "/home/ak121396/multiobjective/solvers/ep+FP/FPepresults/"
# "/home/ak121396/Desktop/FPep2hr/record/"
folders = readdir(direc)
tb = zeros(7,3)
for i=1:7
    f = readdlm(direc*folders[i],',',Float64)
    for k=1:3
        tb[i,k] = mean(f[:,k])
    end
end

direc = "/home/ak121396/multiobjective/solvers/ep+FP/ep_results/"  # "/home/ak121396/Desktop/epsilon2hr/record/"

folders = readdir(direc)#[4:14]
tb = zeros(12,8)
for i=1:12
    f = readdlm(direc*folders[i],',',Float64)
    for k=1:8
        tb[i,k] = mean(f[:,k])
    end
end
################################





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
