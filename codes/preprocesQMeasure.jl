using DelimitedFiles,DataFrames

cd("C:\\Users\\AK121396\\Downloads\\performance_indi\\indicators_win\\")

# find bounds
function Findbounds(refpath,newpath)
    files = readdir(refpath)
        for j=1:length(files)
            @show file = files[j]
            new = newpath*file[1:end-4]*".txt"
            # ./bound [<paramfile>] <datafile> <outfile>
            open("$new","w") do io end
            run(pipeline(`./tools_win/bound ./tools_win/bound_param.txt
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
            run(pipeline(`./tools_win/bound ./tools_win/bound_param.txt
                $refpath2/$file $new`) )
        end
    end
end
Findbounds("F:/results/mergedMIP/","F:/results/performance/bounds/MIPLIB/")


#(Reference file) obj normalisation
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
            run(pipeline(`./tools_win/normalize ./tools_win/normalize_param.txt
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
            run(pipeline(`./tools_win\\filter ./tools_win\\filter_param.txt
                $ref2$/$file $new`) )
        end
    end
end
Filter("F:/results/mergedMIP/norm/")
readdir("F:/results/mergedMIP/norm/1/")[1][1:end-9]


# normalise obj values
function Normalise(probpath,boundpath)
    num = readdir(probpath)
    for i=6:10 #length(num)-1
        insts = readdir(probpath*"$i/")
        for j=1:length(insts)
            file = probpath*"$i/"*insts[j]
            f = readdlm(file)

            # new = probpath*"norm/$dirn/"*insts[j][1:end-4]*"_norm"*".txt"
            new = probpath*"norm/$i/"*insts[j][1:end-4]*"_norm"*".txt"
            open("$new","w") do io
            end
            # bounds = readdir(boundpath); boundf = bounds[j];  #for AP,KP,FLP
            bounds = readdir(boundpath*"$i/") #for MIPLiB
            bpath = boundpath*"$i/"*bounds[j] #*num[i]
            # ./normalize [<paramfile>] <boundfile> <datafile> <outfile>
            run(pipeline(`./tools_win/normalize ./tools_win/normalize_param.txt
                $bpath $file $new`) )

        end
    end
end
# cd("C:\\Users\\AK121396\\Downloads\\performance_indi\\")
Normalise("F:/results/gpr/MIPLIB\\","F:/results/performance/bounds/MIPLIB/")
Normalise("F:/results/fpbh/MIPLIB/","F:/results/performance/bounds/MIPLIB/")


# Measure Uep
function Measures(normalpath,refpath,eppath)
    for i=6:10#length(num)
        files = readdir(normalpath*"$i")
        for j=1:length(files)
            epstore = eppath*"$i/"*files[j][1:end-9]*"_ep.txt"
            open("$epstore","w") do io end; 
            data = normalpath*"$i/"*files[j]
            # refs = readdir(refpath"); ref = refpath*"/"*refs[j]; #for AP,KP,FLP
            refs = readdir(refpath*"$i/"); ref = refpath*"$i/"*refs[j];
            # The order of inputs: .exe file    parameter file    approximation-set file    ref file      outputfile location
            run(pipeline(`./eps_ind ./eps_ind_param.txt $data $ref $epstore`))
        end
    end
end

############################    WINDOWS     ####################################
cd("./indicators_win/")
Measures("F:\\results\\gpr/MIPLIB/norm\\","F:\\results/mergedMIP/norm/", "F:/results/performance/GPR/ep/MIPLIB/")
Measures("F:\\results\\fpbh/MIPLIB/norm\\","F:\\results/mergedMIP/norm/", "F:/results/performance/FPBH/ep/MIPLIB/")
###############################  Linux    ######################################
cd("/home/ak121396/Downloads/performance_indi/indicators_linux/")
readdir("/media/ak121396/0526-8445/results/oldGPR/KP/norm/1/")[1][1:end-9]
Measures("/media/ak121396/0526-8445/results/oldGPR/KP/norm/", "/media/ak121396/0526-8445/results/KS/norm/KP/", "/media/ak121396/0526-8445/results/performance/GPR/ep/KP/",)

# The order of inputs: .exe file    parameter file    approximation-set file    ref file      outputfile location
run(pipeline(`./indicators_win/eps_ind ./indicators_win/eps_ind_param.txt
    $normalpath"*"$f
    F:/results/KS/KP/KP_p-3_n-010_ins-01.txt
    C:/Users/AK121396/Desktop/performance/ex.txt`))
run(pipeline(`./indicators_win/hyp_ind ./indicators_win/hyp_ind_param.txt
    F:/results/GPR/KP/1/n-010_ins-01kpY.log
    F:/results/KS/KP/KP_p-3_n-010_ins-01.txt
    C:/Users/AK121396/Desktop/performance/ex.txt`))
