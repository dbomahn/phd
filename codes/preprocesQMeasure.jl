using DelimitedFiles,DataFrames

cd("C:\\Users\\AK121396\\Downloads\\performance_indi\\indicators_win\\")

# find bounds
dirs = readdir("F:\\results\\KS\\")
for i=1:length(dirs)
    dir = dirs[i]
    files = readdir("F:\\results\\KS\\"*dirs[i])
    for j=1:length(files)
        @show file = files[j]
        pa = "F:/results/performance/$dir/$dir"*"_$j"*".txt"
        # ./bound [<paramfile>] <datafile> <outfile>
        open("$pa","w") do io
        end
        run(pipeline(`./tools_win/bound ./tools_win/bound_param.txt
            F:/results/KS/$dir/$file $pa`) )
    end
end

# normalise obj values
function Normalise(probpath,boundpath)
    num = readdir(probpath)
    for i=5:5#length(num)-1
        insts = readdir(probpath*num[i])
        for j=1:length(insts)
            file = probpath*num[i]*"/"*insts[j]
            f = readdlm(file)
            dirn = num[i];
            new = probpath*"norm/$dirn/"*insts[j][1:end-4]*"_norm"*".txt"
            open("$new","w") do io
            end
            bounds = readdir(boundpath)
            boundf = bounds[j];
            # ./normalize [<paramfile>] <boundfile> <datafile> <outfile>
            run(pipeline(`./tools_linux/normalize ./tools_linux/normalize_param.txt
                $boundpath$boundf $file $new`) )
        end
    end
end
Normalise("F:\\results\\FPBH\\FLP\\","F:\\results\\performance\\bounds\\FLP\\")


#KirlikSayin obj normalisation
prob = readdir("F:\\results\\KS\\norm")
for i=1:length(prob)
    insts = readdir("F:\\results\\KS\\"*prob[i])
    for j=1:length(insts)
        file = "F:\\results\\KS\\"*prob[i]*"\\"*insts[j]
        f = readdlm(file)
        dirn = prob[i]
        new = "F:/results/KS/norm/$dirn/"*insts[j]*"_norm"*".txt"
        open("$new","w") do io
        end
        ins = insts[j]
        bounds = readdir("F:\\results\\performance\\bounds\\"*prob[i])
        boundf = bounds[j]
        # ./normalize [<paramfile>] <boundfile> <datafile> <outfile>
        run(pipeline(`./tools_win/normalize ./tools_win/normalize_param.txt
            F:\\results\\performance\\bounds\\$dirn\\$boundf $file $new`))
    end
end

# Measure Uep & HV
function Measures(normalpath,refpath,eppath,hvpath)
    num = readdir(normalpath)
    for i=1:4#length(num)
        files = readdir(normalpath*num[i])
        for j=1:length(files)
            @show num[i],files[j]
            epstore = eppath*num[i]*"/"*files[j][1:end-9]*"_ep.txt"
            hvstore = hvpath*num[i]*"/"*files[j][1:end-9]*"_hv.txt"
            open("$epstore","w") do io end; open("$hvstore","w") do io end
            data = normalpath*num[i]*"/"*files[j]
            refs = readdir(refpath); ref = refpath*refs[j];
            # The order of inputs: .exe file    parameter file    approximation-set file    ref file      outputfile location
            run(pipeline(`./eps_ind ./eps_ind_param.txt $data $ref $epstore`))
            run(pipeline(`./hyp_ind ./hyp_ind_param.txt $data $ref $hvstore`))
        end
    end
end

readdir("/media/ak121396/0526-8445/results/GPR/AP/norm/2/")[6][1:end-9]
# cd("./indicators_linux/")
# Linux
Measures("/media/ak121396/0526-8445/results/GPR/AP/norm/", "/media/ak121396/0526-8445/results/KS/norm/AP/",
    "/media/ak121396/0526-8445/results/performance/ep/AP/","/media/ak121396/0526-8445/results/performance/hv/AP/")

# WINDOWS
Measures("F:\\results\\GPR\\FLP\\norm\\","F:\\results\\KS\\norm\\FLP\\")

Measures("F:\\results\\GPR\\AP\\norm\\","F:\\results\\KS\\norm\\AP\\")
Measures("F:\\results\\GPR\\KP\\norm\\","F:\\results\\KS\\norm\\KP\\")
Measures("F:\\results\\FPBH\\FLP\\norm\\","F:\\results\\KS\\norm\\FLP\\")
Measures("F:\\results\\FPBH\\AP\\norm\\","F:\\results\\KS\\norm\\AP\\")
Measures("F:\\results\\FPBH\\KP\\norm\\","F:\\results\\KS\\norm\\KP\\")

1
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
