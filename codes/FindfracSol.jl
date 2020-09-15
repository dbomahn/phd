using DataStructures, DelimitedFiles, DataFrames, CSV, Statistics,LinearAlgebra


struct readfile{S<:String}
    d1::S;
    function readfile(d1::S)
        file = readdir(d1)
        for i=1:length(file)
            f = readdlm(d1*file[i])

        end
        return files
    end

end

# d1 = "/home/ak121396/Desktop/KP/X/"
d1 = "E:\\Bensolve\\X\\"; d2 ="E:\\Bensolve\\data\\"
file=readdir(d1); ins = readdir(d2)
for k=1:1#length(file)
k=21

f = readdlm(d1*file[k])
numsol = size(f)[1]
rdownsol =  floor.(f)
d = readdlm(d2*ins[k])
data = readdlm(d2*ins[k],'\t', String, '\n')
b = data[4:length(data)-1]
n=parse(Int,data[2])
C= ones(length(b),n)
C = round.(Int,C)
tb = zeros(numsol,3)
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
        objval = dot(rdownsol[i,:],C[j,:])
        @show tb[i,j] = objval
    end
end
fname = file[k][1:19]
CSV.write("E:\\Bensolve\\RoundDownsol\\"*"$fname"*".csv",tb,writeheader=false, delim=' ' )


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
