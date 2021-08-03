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
num = readdir("F:\\results\\GPR\\FLP\\")
for i=1:length(num)-1
    insts = readdir("F:\\results\\GPR\\FLP\\"*num[i])
    # num = readdir("F:\\results\\FPBH\\"*pth[i]*insts[j])
    for j=1:length(insts)
        # ins = readdir("F:\\results\\FPBH\\FLP\\"*num[i]*"\\"*insts[j])
        file = "F:\\results\\GPR\\FLP\\"*num[i]*"\\"*insts[j]
        f = readdlm(file)
        dirn = num[i]
        new = "F:/results/GPR/FLP/norm/$dirn/"*insts[j]*"_norm"*".txt"
        open("$new","w") do io
        end
        ins = insts[j]
        bounds = readdir("F:\\results\\performance\\bounds\\FLP\\")
        boundf = bounds[j]
        # ./normalize [<paramfile>] <boundfile> <datafile> <outfile>
        run(pipeline(`./tools_win/normalize ./tools_win/normalize_param.txt
            F:\\results\\performance\\bounds\\FLP\\$boundf $file $new`) )
    end
end

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
function Measures(normalpath,refpath)
    num = readdir(normalpath)
    for i=1:length(num)
        @show num
        files = readdir(normalpath*num[i])
        epstore = "F:\\results\\GPR\\FLP\\norm\\"*num[i]*"\\"*num[i]*"ep.txt"
        hvstore = "F:\\results\\GPR\\FLP\\norm\\"*num[i]*"\\"*num[i]*"hv.txt"
        open("$epstore","w") do io end; open("$hvstore","w") do io end

        for j=1:length(files)
            data = normalpath*num[i]*"\\"*files[j]
            refs = readdir(refpath); ref = refpath*refs[j];
            # The order of inputs: .exe file    parameter file    approximation-set file    ref file      outputfile location
            run(pipeline(`./indicators_win/eps_ind ./indicators_win/eps_ind_param.txt
                $data $ref $epstore`))
            run(pipeline(`./indicators_win/hyp_ind ./indicators_win/hyp_ind_param.txt
                $data $ref $hvstore`))
        end
    end
end

Measures("F:\\results\\GPR\\FLP\\norm\\","F:\\results\\KS\\norm\\FLP\\")

readdlm(normalpath*"\\1\\05_010_01flpY.log_norm.txt")


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
