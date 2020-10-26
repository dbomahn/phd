using Statistics,DataFrames,DelimitedFiles
# direc = "C:\\Users\\AK121396\\Desktop\\iteratio_KP_record.csv"
direc = "/home/ak121396/Desktop/PR_KP/basic//record/"
#records
avsol = zeros(10,5);avtime = zeros(10,5);
for i=1:5
    files = readdir(direc)
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
    tb1[i,1] = round(mean(avsol[i,:]),digits=1)
    tb2[i,1] = round(mean(avtime[i,:]),digits=3)
end
tb1
tb2


# HV
pr = "/home/ak121396/Desktop/PR_KP/basic//"

avhv = zeros(5,10)
for u=1:5
    pr = "/home/ak121396/Desktop/PR_KP/basic/"*"$u/"

    # pr = "/home/ak121396/Desktop/PR_KP/repetition/"*"$u"*"/imp/"
    prfiles = readdir(pr)

    table = zeros(10,10);
    for i=1:10
      for j=1:10
        k = j+(i-1)*10
        # ithhv = normHV(ksdir,ksfiles,ksdir,ksfiles,k)
        # ithhv = normHV(ksdir,ksfiles,ben,benfiles,k)
        ithhv = normHV(ksdir,ksfiles,pr,prfiles,k)
        # ithhv = normHV(ksdir,ksfiles,pr50,pr50files,k)
        table[j,i] = ithhv
      end
    end


    for i=1:10
        avhv[u,i] = round(mean(table[:,i]),digits=2)
    end
end

ftb = zeros(10,1)
for i=1:10
     ftb[i,1] = round(mean(avhv[:,i]),digits=2)

end

tt= [ftb[2:10];ftb[1]]
tt #This is the final avg.HV of multiple runs
