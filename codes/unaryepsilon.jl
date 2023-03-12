using DelimitedFiles,DataFrames,CSV

########################## obj file converter: merging FPBH&GPR ################
fdir = "E:/LPBMresults/fpbh/MIPLIB/"
gdir = "E:/LPBMresults/gpr/MIPLIB/"
for i=1:14
    gname = readdir(gdir*"/$i/"); fname = readdir(fdir*"/$i/")
    for j=1:5
        F = readdlm(fdir*"/$i/"*fname[j])
        G = readdlm(gdir*"/$i/"*gname[j])
        P = vcat(F,G)
        ins = fname[j][1:end-4]
        CSV.write("E:/LPBMresults/mergedMIP/$i"*"/$ins"*".txt",DataFrame(P, :auto),header=false, delim=' ' )
    end
end
# edir = "/home/ak121396/Desktop/relise/epsilon/"
# mdir = "/home/ak121396/Desktop/relise/lpY/5/"
# # for i=1:9
# mname = readdir(mdir); ename = readdir(edir)[2:end]
# for j=1:5
#     E = readdlm(edir*ename[j])
#     M = readdlm(mdir*mname[j])
#     P = vcat(E,M)
#     ins = ename[j][1:end-7]
#     CSV.write("/home/ak121396/Desktop/relise/merged/$ins"*"merged.txt",DataFrame(P, :auto),header=false, delim=' ' )
# end
# end
######################################  Find Bounds  ######################################################
cd("C:\\Users\\AK121396\\Downloads\\performance_indi\\indicators_win\\")

cd("/home/ak121396/Downloads/performance_indi/indicators_linux/")
# find bounds
os = "_win" #"_linux"
function Findbounds(refpath,newpath)
    files = readdir(refpath)
        for j=1:length(files)
            @show file = files[j]
            new = newpath*file[1:end-4]*".txt"
            # ./bound [<paramfile>] <datafile> <outfile>
            open("$new","w") do io end
            run(pipeline(`./tools*"$os"*/bound ./tools*"$os"*/bound_param.txt
                $refpath$file $new`) )
        end
end

cd("../")
function Findbounds(refpath,newpath) #MIPLIB
    for i=1:14
        @show refpath2 = refpath*"$i"
        files = readdir(refpath2)
        for j=1:length(files)
            @show file = files[j]
            new = newpath*"$i/"*file[1:end-4]*".txt"
            # ./bound [<paramfile>] <datafile> <outfile>
            open("$new","w") do io end
            run(pipeline(`./tools$os/bound ./tools$os/bound_param.txt
                $refpath2/$file $new`) )
        end
    end
end
Findbounds("E:/LPBMresults/mergedMIP/","E:/LPBMresults/performance/bounds/MIPLIB/")
# cd("/home/ak121396/Downloads/performance_indi/") #indicators_linux/
# # find bounds
# os = "_linux" #_win
# function Findbounds(refpath,newpath) #MIPLIB
#     # for i=6:10
#         @show refpath2 = refpath#*"$i"
#         files = readdir(refpath2)
#         for j=1:length(files)
#             @show file = files[j]
#             new = newpath*file[1:end-10]*".txt"
#             # ./bound [<paramfile>] <datafile> <outfile>
#             open("$new","w") do io end
#             run(pipeline(`./tools$os/bound $refpath2/$file $new`) )
#         end
#     # end
# end
# Findbounds(,"/home/ak121396/Desktop/relise/performance/bounds/test01S2.txt")


################################  (Reference file) obj normalisation  ##################################
function Ref_normalise(refpath,boundpath)
    prob = readdir(refpath);
    for i=1:14 #length(prob)
        # dirn = prob[i]
        insts = readdir(refpath*"$i/")#dirn)
        for j=1:length(insts)
            file = refpath*"$i/"*insts[j]
            f = readdlm(file)
            new = refpath*"norm/$i/"*insts[j][1:end-4]*"_norm"*".txt"
            open("$new","w") do io
            end
            # bounds = readdir(boundpath*prob[i]);  boundf = bounds[j] #for AP,KP,FLP
            bounds = readdir(boundpath*"$i/") #for MIPLiB
            @show bpath = boundpath*"$i/"*bounds[j] #*num[i]

            # ./normalize [<paramfile>] <boundfile> <datafile> <outfile>
            run(pipeline(`./tools$os/normalize ./tools$os/normalize_param.txt
                $bpath $file $new`))

        end
    end
end

cd("C:\\Users\\AK121396\\Downloads\\performance_indi\\")

Ref_normalise("E:/LPBMresults/mergedMIP/","E:\\LPBMresults\\performance\\bounds\\MIPLIB\\")
Ref_normalise("E:/LPBMresults/KS/",)


# merged Ref file filter
function Filter(refn,reff)
    for j=1:14
        ref2 = refn*"$j/"
        ins = readdir(ref2)
        for i=1:length(ins)
            new = reff*"/$j/"*ins[i][1:end-13]*"_filter"*".txt"
            open("$new","w") do io
            end
            file = ins[i]
            run(pipeline(`./tools$os/filter ./tools$os/filter_param.txt
                $ref2$/$file $new`) )
        end
    end
end
Filter("E:/LPBMresults/mergedMIP/norm/","E:/LPBMresults/mergedMIP/filtered/")
readdir("E:/LPBMresults/mergedMIP/norm/1/")[1][1:end-9]


# normalise obj values
function Normalise(probpath,boundpath)
    # num = readdir(probpath)
    for i=1:14#length(num)-1
        insts = readdir(probpath*"$i/")
        for j=1:length(insts)
            file = probpath*"$i/"*insts[j]
            # new = probpath*"norm/$dirn/"*insts[j][1:end-4]*"_norm"*".txt"
            new = probpath*"norm/$i/"*insts[j][1:end-4]*"_norm"*".txt"
            open("$new","w") do io
            end
            # bounds = readdir(boundpath);   #for AP,KP,FLP
            bounds = readdir(boundpath*"$i/") #for MIPLiB
            bfile = boundpath*"$i/"*bounds[j] #*num[i] "$i/"*
            # ./normalize [<paramfile>] <boundfile> <datafile> <outfile>
            run(pipeline(`./tools$os/normalize ./tools$os/normalize_param.txt
                $bfile $file $new`) )

        end
    end
end
cd("C:\\Users\\AK121396\\Downloads\\performance_indi\\")
# cd("/home/ak121396/Downloads/performance_indi/")
Normalise("E:/LPBMresults/gpr/MIPLIB/","E:/LPBMresults/performance/bounds/MIPLIB/")
Normalise("E:/LPBMresults/fpbh/MIPLIB/","E:/LPBMresults/performance/bounds/MIPLIB/")

Normalise("E:/LPBMresults/gpr/AP/","E:/LPBMresults/performance/bounds/AP/")
Normalise("E:/LPBMresults/gpr/FLP/","E:/LPBMresults/performance/bounds/FLP/")

Normalise("E:/LPBMresults/gpr/KP/","E:/LPBMresults/performance/bounds/KP/")

# Measure Uep
function Measures(normalpath,refpath,eppath)
    for i=1:14#length(num)
        files = readdir(normalpath*"$i")
        for j=1:length(files)
            epstore = eppath*"$i/"*files[j][1:end-9]*"_ep.txt"
            open("$epstore","w") do io end;
            data = normalpath*"$i/"*files[j]
            # refs = readdir(refpath); ref = refpath*"/"*refs[j]; #for AP,KP,FLP
            refs = readdir(refpath*"$i/"); ref = refpath*"$i/"*refs[j]; #for MIPLIB
            # The order of inputs: .exe file    parameter file    approximation-set file    ref file      outputfile location
            run(pipeline(`./indicators$os/eps_ind ./indicators$os/eps_ind_param.txt $data $ref $epstore`))
        end
    end
end
############################    WINDOWS     ####################################
Measures("E:/LPBMresults/gpr/MIPLIB/norm/","E:/LPBMresults/mergedMIP/filtered/", "E:/LPBMresults/performance/GPR/ep/MIPLIB/")
Measures("E:\\LPBMresults\\fpbh/MIPLIB/norm\\","E:\\LPBMresults/mergedMIP/filtered\\", "E:/LPBMresults/performance/FPBH/ep/MIPLIB/")

Measures("E:/LPBMresults/gpr/AP/norm/","E:/LPBMresults/KS/norm/AP/", "E:/LPBMresults/performance/GPR/ep/AP/")
Measures("E:/LPBMresults/gpr/FLP/norm/","E:/LPBMresults/KS/norm/FLP/","E:/LPBMresults/performance/GPR/ep/FLP/")
Measures("E:/LPBMresults/gpr/KP/norm/","E:/LPBMresults/KS/norm/KP/", "E:/LPBMresults/performance/GPR/ep/KP/")

###############################  Linux    ######################################
cd("/home/ak121396/Downloads/performance_indi/indicators_linux/")
readdir("E:/LPBMresults/gpr//MIPLIB/norm/2/")

Measures("E:/LPBMresults/gpr/MIPLIB/neos/norm/", "E:/LPBMresults/mergedMIP/neos/norm/", "E:/LPBMresults/performance/GPR/ep/MIPLIB/neos/")
Measures("E:/LPBMresults/fpbh/MIPLIB/neos/norm/", "E:/LPBMresults/mergedMIP/neos/norm/", "E:/LPBMresults/performance/FPBH/ep/MIPLIB/neos/")

# neos w/ shorter TL
function Measures(normalpath,refpath,eppath)
        files = readdir(normalpath)
        for j=1:length(files)
            epstore = eppath*files[j][1:end-9]*"_ep.txt"
            open("$epstore","w") do io end;
            data = normalpath*files[j]
            # refs = readdir(refpath); ref = refpath*"/"*refs[j]; #for AP,KP,FLP
            refs = readdir(refpath); ref = refpath*refs[j];
            # The order of inputs: .exe file    parameter file    approximation-set file    ref file      outputfile location
            run(pipeline(`./eps_ind ./eps_ind_param.txt $data $ref $epstore`))
        end
    end
end


# The order of inputs: .exe file    parameter file    approximation-set file    ref file      outputfile location
run(pipeline(`./indicators_win/eps_ind ./indicators_win/eps_ind_param.txt
    $normalpath"*"$f
    E:/LPBMresults/KS/KP/KP_p-3_n-010_ins-01.txt
    C:/Users/AK121396/Desktop/performance/ex.txt`))
run(pipeline(`./indicators_win/hyp_ind ./indicators_win/hyp_ind_param.txt
    E:/LPBMresults/GPR/KP/1/n-010_ins-01kpY.log
    E:/LPBMresults/KS/KP/KP_p-3_n-010_ins-01.txt
    C:/Users/AK121396/Desktop/performance/ex.txt`))


######################### Linux #######################################
# cd("/home/ak121396//Downloads/performance_indi/")
# run(pipeline(`cd /home/ak121396/Downloads/performance_indi/`))
# run(`cd /home/ak121396/Downloads/performance_indi/`)

run(pipeline(`cat '>' /home/ak121396/Desktop/performance/ex.txt`))

# SAMPLE
run(pipeline(`./indicators_linux/eps_ind ./indicators_linux/eps_ind_param.txt
    /home/ak121396/Desktop/relise/epsilon/test01S2epY.log  /home/ak121396/Desktop/ref-1.txt /home/ak121396/Desktop/ex.txt`) )
run(pipeline(`./indicators_linux/hyp_ind ./indicators_linux/hyp_ind_param.txt /home/ak121396/Desktop/data.txt
        /home/ak121396/Desktop/ref-1.txt /home/ak121396/Desktop/performance/ex.txt`) )


run(pipeline(`./indicators_linux/eps_ind ./indicators_linux/eps_ind_param.txt
    /home/ak121396/Desktop/data.txt /home/ak121396/Desktop/ref-1.txt /home/ak121396/Desktop/ex.txt`) )



############################# Convert FPBH file ###############################
for j=6:10#length(fdir)-1
    files = readdir("E:/LPBMresults/fpbh\\MIPLIB/$j/")
    for i=1:length(files)
        fname = files[i]
        P = readdlm("E:/LPBMresults/fpbh/MIPLIB/$j/"*fname)
        ins = fname[1:end-7]
        CSV.write("E:/LPBMresults/fpbh/MIPLIB/$j/"*"$ins"*".txt",DataFrame(P, :auto),header=false, delim=' ' )
    end
end
