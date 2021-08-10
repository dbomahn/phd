using DelimitedFiles,DataFrames

cd("C:\\Users\\AK121396\\Downloads\\performance_indi\\indicators_win\\")

# find bounds
function Findbounds(refpath,newpath)
    for i=1:length(refpath)
        dir = refpath[i]
        files = readdir(refpath*dir)
        for j=1:length(files)
            @show file = files[j]
            new = newpath*dir*file[1:end-4]*".txt"
            # ./bound [<paramfile>] <datafile> <outfile>
            open("$new","w") do io end
            # run(pipeline(`./tools_win/bound ./tools_win/bound_param.txt
                # $refpath$dir$file $new`) )
            run(pipeline(`./tools_linux/bound ./tools_linux/bound_param.txt
                $refpath$dir$file $new`) )
        end
    end
end

cd("../")
Findbounds("F:\\results\\mergedMIP\\","F:\\results\\performance\\bounds\\MIPLIB\\")
Findbounds("/media/ak121396/0526-8445/results/KS/KP","/media/ak121396/0526-8445/results/performance/bounds/KP")

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
            run(pipeline(`./tools_linux/normalize ./tools_linux/normalize_param.txt
                $boundpath//$boundf $file $new`) )
        end
    end
end
cd("C:\\Users\\AK121396\\Downloads\\performance_indi\\")
Normalise("F:\\results\\FPBH\\MIPLIB\\","F:\\results\\performance\\bounds\\MIPLIB\\")
Normalise("/media/ak121396/0526-8445/results/oldGPR/KP/","/media/ak121396/0526-8445/results/performance/bounds/KP/")
cd("./")
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
    for i=1:length(prob)
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
            # run(pipeline(`./tools_win/normalize ./tools_win/normalize_param.txt
            #     $boundpath\\$dirn\\$boundf $file $new`))
            run(pipeline(`../tools_linux/normalize ../tools_linux/normalize_param.txt
                $boundpath/$dirn/$boundf $file $new`))
        end
    end
end

Ref_normalise("F:\\results\\mergedMIP\\","F:\\results\\performance\\bounds\\MIPLIB\\")
readdir("F:\\results\\mergedMIP\\1\\")[1][1:end-4]
Ref_normalise("/media/ak121396/0526-8445/results/KS/",)

1
# Measure Uep & HV
function Measures(normalpath,refpath,eppath,hvpath)
    num = readdir(normalpath)
    for i=1:length(num)
        files = readdir(normalpath*num[i])
        for j=1:length(files)
            # @show num[i],files[j]
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

###############################  Linux    ######################################
cd("/home/ak121396/Downloads/performance_indi/indicators_linux/")
readdir("/media/ak121396/0526-8445/results/oldGPR/KP/norm/1/")[1][1:end-9]
Measures("/media/ak121396/0526-8445/results/oldGPR/KP/norm/", "/media/ak121396/0526-8445/results/KS/norm/KP/",
    "/media/ak121396/0526-8445/results/performance/GPR/ep/KP/","/media/ak121396/0526-8445/results/performance/GPR/hv/KP/")


Measures("/media/ak121396/0526-8445/results/oldGPR/MIPLIB/norm/", "/media/ak121396/0526-8445/results/mergedMIP/norm/1/",
    "/media/ak121396/0526-8445/results/performance/ep/MIPLIB/","/media/ak121396/0526-8445/results/performance/hv/MIPLIB/")

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

gep = [4.255319200e-02,6.470588200e-02,1.164680390e-01,5.385915700e-02]
ghv =
FPBH = []
./hv/MIPLIB/1/cvs08r139-94_hv.txt:4.655709520e-02
./hv/MIPLIB/1/n2seq36f_hv.txt:1.135397310e-02
./hv/MIPLIB/1/neos-1516309_hv.txt:3.644017939e-01
./hv/MIPLIB/1/neos-1599274_hv.txt:1.055278103e-01


function Readmeasures(mpath)
    hpath = readdir(mpath*"/hv/MIPLIB/1/"); epath = readdir(mpath*"/ep/MIPLIB/1/")
    hvarr = []; eparr = [];
    for i=1:length(hpath)
        hfig = readdlm(mpath*"/hv/MIPLIB/1/"*hpath[i]); efig = readdlm(mpath*"/ep/MIPLIB/1/"*epath[i]);
        push!(hvarr,hfig); push!(eparr,efig);
    end
    return hvarr,eparr
end
ghv, gep =  Readmeasures("/media/ak121396/0526-8445/results/performance/GPR/")
fhv, fep =  Readmeasures("/media/ak121396/0526-8445/results/performance/FPBH/")

println("HV: GPR vs FPBH \n", hcat(ghv,fhv))
println("ep: GPR vs FPBH \n", hcat(gep,fep))
