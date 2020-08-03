using DataStructures, DelimitedFiles, DataFrames, CSV, Statistics


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

d1 = "/home/ak121396/Desktop/KP/X/"
file=readdir("/home/ak121396/Desktop/KP/X/")
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
