using DelimitedFiles,DataFrames

cd("C:\\Users\\AK121396\\Downloads\\performance_indi\\indicators_win\\")

# find bounds
function Findbounds(refpath,newpath)
    for i=1:length(dirs)
        dir = dirs[i]*"/"
        files = readdir(refpath*dir)
        for j=1:length(files)
            @show file = files[j]
            new = newpath*dir*file[1:end-4]*".txt"
            # ./bound [<paramfile>] <datafile> <outfile>
            open("$new","w") do io end
            run(pipeline(`./tools_win/bound ./tools_win/bound_param.txt
                $refpath$dir$file $new`) )
        end
    end
end

cd("../../performance_indi\\")
Findbounds("F:\\results\\mergedMIP\\","F:\\results\\performance\\bounds\\MIPLIB\\")

# normalise obj values
function Normalise(probpath,boundpath)
    num = readdir(probpath)
    for i=1:length(num)-1
        insts = readdir(probpath*num[i])
        for j=1:length(insts)
            file = probpath*num[i]*"/"*insts[j]
            f = readdlm(file)
            dirn = num[i];
            new = probpath*"norm/$dirn/"*insts[j][1:end-4]*"_norm"*".txt"
            open("$new","w") do io
            end
            bounds = readdir(boundpath); boundf = bounds[j];  #for AP,KP,FLP
            # bounds = readdir(boundpath*dirn) #for MIPLiB
            bpath = boundpath*num[i]*"/"*bounds[j]
            # ./normalize [<paramfile>] <boundfile> <datafile> <outfile>
            # run(pipeline(`./tools_win/normalize ./tools_win/normalize_param.txt
            #     $bpath $file $new`) )
            run(pipeline(`./tools_win/normalize ./tools_win/normalize_param.txt
                $boundpath//$boundf $file $new`) )
        end
    end
end
cd("C:\\Users\\AK121396\\Downloads\\performance_indi\\")
Normalise("F:\\results\\FPBH\\MIPLIB\\","F:\\results\\performance\\bounds\\MIPLIB\\")
Normalise("F:\\results\\FPBH\\KP\\","F:\\results\\performance\\bounds\\KP\\")

readdir("F:\\results\\FPBH\\KP\\1")[1][1:end-4]

# merged Ref file filter
function Filter(ref)
    ins = readdir(ref)
    for i=1:length(ins)
        new = ref*ins[i][1:end-9]*"_filter"*".txt"
        open("$new","w") do io
        end
        file = ins[i]
        run(pipeline(`./tools_win\\filter ./tools_win\\filter_param.txt
            $ref/$file $new`) )
    end
end
Filter("F:/results/mergedMIP/norm/1/")
readdir("F:/results/mergedMIP/norm/1/")[1][1:end-9]


#KirlikSayin (Reference file) obj normalisation
function Ref_normalise(refpath,boundpath)
    prob = readdir(refpath)
    for i=1:1#length(prob)
        insts = readdir(refpath*prob[i])
        for j=1:length(insts)
            file = refpath*prob[i]*"/"*insts[j]
            f = readdlm(file)
            dirn = prob[i]
            new = refpath*"norm/$dirn"*"/"*insts[j][1:end-4]*"_norm"*".txt"
            open("$new","w") do io
            end
            ins = insts[j]
            bounds = readdir(boundpath*dirn); boundf = bounds[j]
            # bounds = readdir(boundpath*prob[i]);  boundf = bounds[j] #for AP,KP,FLP
            # ./normalize [<paramfile>] <boundfile> <datafile> <outfile>
            run(pipeline(`./tools_win/normalize ./tools_win/normalize_param.txt
                $boundpath\\$dirn\\$boundf $file $new`))
        end
    end
end

Ref_normalise("F:\\results\\mergedMIP\\","F:\\results\\performance\\bounds\\MIPLIB\\")
readdir("F:\\results\\mergedMIP\\1\\")[1][1:end-4]
1
# Measure Uep & HV
function Measures(normalpath,refpath,eppath,hvpath)
    num = readdir(normalpath)
    for i=1:length(num)
        files = readdir(normalpath*num[i])
        for j=1:length(files)
            @show num[i],files[j]
            epstore = eppath*num[i]*"/"*files[j][1:end-9]*"_ep.txt"
            hvstore = hvpath*num[i]*"/"*files[j][1:end-9]*"_hv.txt"
            open("$epstore","w") do io end; open("$hvstore","w") do io end
            data = normalpath*num[i]*"/"*files[j]
            refs = readdir(refpath); ref = refpath*refs[j]; #for AP,KP,FLP
            # refs = readdir(refpath*num[i]); ref = refpath*num[i]*"/"*refs[j];
            # The order of inputs: .exe file    parameter file    approximation-set file    ref file      outputfile location
            run(pipeline(`./eps_ind ./eps_ind_param.txt $data $ref $epstore`))
            run(pipeline(`./hyp_ind ./hyp_ind_param.txt $data $ref $hvstore`))
        end
    end
end

############################    WINDOWS     ####################################
cd("./indicators_win/")
readdir("F:\\results\\FPBH\\FLP\\norm\\1")[1][1:end-9]
Measures("F:\\results\\FPBH\\FLP\\norm\\","F:\\results\\KS\\norm\\FLP\\",
    "F:/results/performance/ep/FLP\\", "F:/results/performance/hv/FLP\\")

Measures("F:\\results\\GPR\\AP\\norm\\","F:\\results\\KS\\norm\\AP\\")
Measures("F:\\results\\GPR\\KP\\norm\\","F:\\results\\KS\\norm\\KP\\")
Measures("F:\\results\\FPBH\\FLP\\norm\\","F:\\results\\KS\\norm\\FLP\\")
Measures("F:\\results\\FPBH\\AP\\norm\\","F:\\results\\KS\\norm\\AP\\")
Measures("F:\\results\\FPBH\\KP\\norm\\","F:\\results\\KS\\norm\\KP\\")

###############################  Linux    ######################################
cd("./indicators_linux/")
readdir("/media/ak121396/0526-8445/results/GPR/AP/norm/2/")[6][1:end-9]
Measures("/media/ak121396/0526-8445/results/GPR/AP/norm/", "/media/ak121396/0526-8445/results/KS/norm/AP/",
    "/media/ak121396/0526-8445/results/performance/ep/AP/","/media/ak121396/0526-8445/results/performance/hv/AP/")

# The order of inputs: .exe file    parameter file    approximation-set file    ref file      outputfile location
run(pipeline(`./indicators_win/eps_ind ./indicators_win/eps_ind_param.txt
    $normalpath"*"$f
    F:/results/KS/KP/KP_p-3_n-010_ins-01.txt
    C:/Users/AK121396/Desktop/performance/ex.txt`))
run(pipeline(`./indicators_win/hyp_ind ./indicators_win/hyp_ind_param.txt
    F:/results/GPR/KP/1/n-010_ins-01kpY.log
    F:/results/KS/KP/KP_p-3_n-010_ins-01.txt
    C:/Users/AK121396/Desktop/performance/ex.txt`))


readdir(normalpath*num[1])
