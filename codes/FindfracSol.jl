using DataStructures, DelimitedFiles, DataFrames, CSV, Statistics,LinearAlgebra

function domFilter(sol,obj) #for Maximisation
    copysol = Dict(); copyobj = Dict();
    for i=1:size(obj)[1]
        copysol[i] = sol[i,:]
        copyobj[i] = obj[i,:]
    end

    for i=1:size(obj)[1]-1
        for j=i+1:size(obj)[1]
            if all(obj[i,:] .<= obj[j,:]) == true #dominated by PF[j]
                copyobj[i]=nothing; copysol[i]=nothing;break
            elseif all(obj[j,:] .<= obj[i,:]) == true
                copyobj[j]=nothing; copysol[j]=nothing;
            end
        end
    end

    finalsol = filter!(a->a!=nothing, collect(values(copysol)))
    finalobj = filter!(a->a!=nothing, collect(values(copyobj)))

    return finalsol,finalobj
end

d1 = "/home/ak121396/Desktop/KP/RoundDown//X/"; d2 = "/home/ak121396/Desktop/KP/data/"
file=readdir(d1); ins = readdir(d2)

global rdnum = zeros(100,1)
f = readdlm(d1*file[11])
ct = count(x->(x==1),f)

intcount = []
for k=1:length(file)
    f = readdlm(d1*file[k])
    # rdownsol = floor.(round.(f,digits=5))
    ct = count(x->(x==1),f)

    push!(intcount,ct)
end

rdownsol = floor.(round.(f,digits=5)); rdownsol = Int.(rdownsol)

###################       1-variables      ##############


global intcount = []
for k=1:length(file)
    f = readdlm(d1*file[k])
    numsol = size(f)[1]
    rdownsol = floor.(round.(f,digits=5)); rdownsol = Int.(rdownsol)

    d = readdlm(d2*ins[k])
    data = readdlm(d2*ins[k],'\t', String, '\n')
    b = data[4:length(data)-1]
    n=parse(Int,data[2])
    C= ones(length(b),n)
    C = round.(Int,C)
    tb = zeros(Int64,numsol,3)
    for x=1:length(b)
        a = split(b[x],('[',']',','))
        aa = filter!(e->!(e in [ "" ,"[","," ,"]"]) ,a)
        for y=1:length(aa)
            c=parse(Int64,aa[y])
            C[x,y] = c
            if c==0
                global ct = ct+1;
            end
        end
    end
    for i=1:numsol
        for j=1:3
            global objval = Int(dot(rdownsol[i,:],C[j,:]))
            tb[i,j] = objval
        end
    end

    ndsol,ndpoint = domFilter(rdownsol,tb)
    push!(intcount,length(ndsol))

    otable = ones(Int,length(ndsol),3)
    rdnum[k] = length(ndpoint)
    for i=1:length(ndsol)
        for j=1:3
            otable[i,j] = -ndpoint[i][j]
        end
    end
    o1,o2,o3=[hcat(-ndpoint...)[i,:] for i=1:3]

    # otable = ones(Int,length(ndsol),n)
    # for i=1:length(ndsol)
    #     for j=1:n
    #         otable[i,j] = ndsol[i][j]
    #     end
    # end
    # dv = DataFrame(otable);
    fname = file[k][1:19]
    CSV.write("/home/ak121396/Desktop/KP/RoundDown/X/"*"$fname"*"_roundown_X.log",dv,header=false, delim=' ' )
end

otable = DataFrame(rdownsol)
fname = file[11][1:19]

CSV.write("/home/ak121396/Desktop/KP/RoundDown/"*"$fname"*"_roundownY.log",otable,head=false, delim=' ' )


avsol = reshape(pc,10,10)
tt = []
for i=1:10
    a = round(mean(avsol[:,i]),digits=4)
    push!(tt,a)
end

tb = zeros(100,3)
for k=1:length(file)
    f = readdlm(d1*file[k])
    r,c = size(f)
    rcount = 0;  track = Dict();
    for i=1:r
        rf = round.(f[i,:],digits=4)
        if all(x->x ∈ [0,1], rf ) == true #if all variables are integer
            rcount+=1
        end
        ccount = 0;
        for j=1:c
            if rf[j] ∉ [0,1]
                ccount+=1
            end
        end
        track[i] = ccount
    end
    tb[k,1] = r; tb[k,2] = r-rcount; tb[k,3] = sum(collect(values(track)))/(tb[k,2]);
    # solfile = file[k][1:19]
    # avginfo = DataFrame( Bensol=r, frSol = r-rcount, avg_num_frVar = mean(collect(values(track))) )
    # insertcols!(avginfo, 4, Symbol("$solfile")=> 0)
    # CSV.write("Desktop/FPBH/Bensol_Info.csv", avginfo, append=true, header=true)
end

# info = DataFrame( no_frac = sort!(collect(values(track))) )
# insertcols!(info,2, Symbol("$solfile")=> sort!(collect(keys(track))))
#
# CSV.write("Desktop/FPBH/Bensol_Info.csv",info ,header=true, append=true)
