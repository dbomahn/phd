using DelimitedFiles,DataFrames

########################## obj file converter: merging FPBH&GPR ################
fdir = "/media/ak121396/0526-8445/results/fpbh/MIPLIB/"
gdir = "/media/ak121396/0526-8445/results/gpr/MIPLIB/"
for i=1:9
    gname = readdir(gdir*"/$i/"); fname = readdir(fdir*"/$i/")
    for j=1:5
        F = readdlm(fdir*"/$i/"*fname[j])
        G = readdlm(gdir*"/$i/"*gname[j])
        P = vcat(F,G)
        ins = fname[j][1:end-4]
        CSV.write("/media/ak121396/0526-8445/results/mergedMIP/$i"*"/$ins"*".txt",DataFrame(P, :auto),header=false, delim=' ' )
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
os = "_linux" #_win
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
    for i=6:10
        @show refpath2 = refpath*"$i"
        files = readdir(refpath2)
        for j=1:length(files)
            @show file = files[j]
            new = newpath*"$i/"*file[1:end-4]*".txt"
            # ./bound [<paramfile>] <datafile> <outfile>
            open("$new","w") do io end
            run(pipeline(`./tools*"$os"*/bound ./tools*"$os"*/bound_param.txt
                $refpath2/$file $new`) )
        end
    end
end
Findbounds("F:/results/mergedMIP/","F:/results/performance/bounds/MIPLIB/")
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
    for i=5:10 #length(prob)
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
            run(pipeline(`./tools_linux/normalize ./tools_linux/normalize_param.txt
                $bpath $file $new`))

        end
    end
end

cd("C:\\Users\\AK121396\\Downloads\\performance_indi\\")

Ref_normalise("F:\\results/mergedMIP/","F:\\results\\performance\\bounds\\MIPLIB\\")
Ref_normalise("/media/ak121396/0526-8445/results/KS/",)


# merged Ref file filter
function Filter(ref)
    for j=6:10
        ref2 = ref*"$j/"
        ins = readdir(ref2)
        for i=1:length(ins)
            new = ref2*ins[i][1:end-13]*"_filter"*".txt"
            open("$new","w") do io
            end
            file = ins[i]
            run(pipeline(`./tools*"$os"*/filter ./tools*"$os"*/filter_param.txt
                $ref2$/$file $new`) )
        end
    end
end
Filter("F:/results/mergedMIP/norm/")
readdir("F:/results/mergedMIP/norm/1/")[1][1:end-9]


# normalise obj values
function Normalise(probpath,boundpath)
    # num = readdir(probpath)
    for i=1:5 #length(num)-1
        insts = readdir(probpath*"$i/")
        for j=1:length(insts)
            file = probpath*"$i/"*insts[j]
            # new = probpath*"norm/$dirn/"*insts[j][1:end-4]*"_norm"*".txt"
            new = probpath*"norm/$i/"*insts[j][1:end-4]*"_norm"*".txt"
            open("$new","w") do io
            end
            bounds = readdir(boundpath);   #for AP,KP,FLP
            # bounds = readdir(boundpath*"$i/") #for MIPLiB
            bfile = boundpath*bounds[j] #*num[i] "$i/"*
            # ./normalize [<paramfile>] <boundfile> <datafile> <outfile>
            run(pipeline(`./tools_linux/normalize ./tools_linux/normalize_param.txt
                $bfile $file $new`) )

        end
    end
end
cd("/home/ak121396/Downloads/performance_indi/")
Normalise("F:/results/gpr/MIPLIB\\","F:/results/performance/bounds/MIPLIB/")
Normalise("/media/ak121396/0526-8445/results/gpr/AP/","/media/ak121396/0526-8445/results/performance/bounds/AP/")
Normalise("/media/ak121396/0526-8445/results/gpr/FLP/","/media/ak121396/0526-8445/results/performance/bounds/FLP/")

Normalise("/media/ak121396/0526-8445/results/gpr/KP/","/media/ak121396/0526-8445/results/performance/bounds/KP/")

# Measure Uep
function Measures(normalpath,refpath,eppath)
    for i=1:5#length(num)
        files = readdir(normalpath*"$i")
        for j=1:length(files)
            epstore = eppath*"$i/"*files[j][1:end-9]*"_ep.txt"
            open("$epstore","w") do io end;
            data = normalpath*"$i/"*files[j]
            refs = readdir(refpath); ref = refpath*"/"*refs[j]; #for AP,KP,FLP
            # refs = readdir(refpath*"$i/"); ref = refpath*"$i/"*refs[j];
            # The order of inputs: .exe file    parameter file    approximation-set file    ref file      outputfile location
            run(pipeline(`./indicators_linux/eps_ind ./indicators_linux/eps_ind_param.txt $data $ref $epstore`))
        end
    end
end
Measures("/media/ak121396/0526-8445/results/gpr/AP/norm/","/media/ak121396/0526-8445/results/KS/norm/AP/", "/media/ak121396/0526-8445/results/performance/GPR/ep/AP/")
Measures("/media/ak121396/0526-8445/results/gpr/FLP/norm/","/media/ak121396/0526-8445/results/KS/norm/FLP/","/media/ak121396/0526-8445/results/performance/GPR/ep/FLP/")
Measures("/media/ak121396/0526-8445/results/gpr/KP/norm/","/media/ak121396/0526-8445/results/KS/norm/KP/", "/media/ak121396/0526-8445/results/performance/GPR/ep/KP/")


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
############################    WINDOWS     ####################################
Measures("F:\\results\\gpr/MIPLIB/norm\\","F:\\results/mergedMIP/norm/", "F:/results/performance/GPR/ep/MIPLIB/")
Measures("F:\\results\\fpbh/MIPLIB/norm\\","F:\\results/mergedMIP/norm/", "F:/results/performance/FPBH/ep/MIPLIB/")
###############################  Linux    ######################################
cd("/home/ak121396/Downloads/performance_indi/indicators_linux/")
readdir("/media/ak121396/0526-8445/results/gpr//MIPLIB/norm/2/")

Measures("/media/ak121396/0526-8445/results/gpr/MIPLIB/neos/norm/", "/media/ak121396/0526-8445/results/mergedMIP/neos/norm/", "/media/ak121396/0526-8445/results/performance/GPR/ep/MIPLIB/neos/")
Measures("/media/ak121396/0526-8445/results/fpbh/MIPLIB/neos/norm/", "/media/ak121396/0526-8445/results/mergedMIP/neos/norm/", "/media/ak121396/0526-8445/results/performance/FPBH/ep/MIPLIB/neos/")


# The order of inputs: .exe file    parameter file    approximation-set file    ref file      outputfile location
run(pipeline(`./indicators_win/eps_ind ./indicators_win/eps_ind_param.txt
    $normalpath"*"$f
    F:/results/KS/KP/KP_p-3_n-010_ins-01.txt
    C:/Users/AK121396/Desktop/performance/ex.txt`))
run(pipeline(`./indicators_win/hyp_ind ./indicators_win/hyp_ind_param.txt
    F:/results/GPR/KP/1/n-010_ins-01kpY.log
    F:/results/KS/KP/KP_p-3_n-010_ins-01.txt
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
    files = readdir("F:/results/fpbh\\MIPLIB/$j/")
    for i=1:length(files)
        fname = files[i]
        P = readdlm("F:/results/fpbh/MIPLIB/$j/"*fname)
        ins = fname[1:end-7]
        CSV.write("F:/results/fpbh/MIPLIB/$j/"*"$ins"*".txt",DataFrame(P, :auto),header=false, delim=' ' )
    end
end
