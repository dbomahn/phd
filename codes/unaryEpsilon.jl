using DelimitedFiles,DataFrames
#############################    WINDOWS    ##################################

open("C:/Users/AK121396/Desktop/performance/ex3.txt","w") do io
end

cd("C:\\Users\\AK121396\\Downloads\\performance_indi\\")


run(pipeline(`./indicators_win/hyp_ind ./indicators_win/hyp_ind_param.txt
    F:/results/GPR/KP/1/n-010_ins-01kpY.log
    F:/results/KS/KP/KP_p-3_n-010_ins-01.txt
    C:/Users/AK121396/Desktop/performance/ex.txt`) )

######################### Linux #######################################
cd("/home/ak121396//Downloads/performance_indi/")
# run(pipeline(`cd /home/ak121396/Downloads/performance_indi/`))
# run(`cd /home/ak121396/Downloads/performance_indi/`)

run(pipeline(`cat '>' /home/ak121396/Desktop/performance/ex.txt`))

# SAMPLE
run(pipeline(`./eps_ind ./eps_ind_param.txt $data $ref $epstore`))
run(pipeline(`./indicators_linux/eps_ind ./indicators_linux/eps_ind_param.txt
    /home/ak121396/Desktop/data.txt /home/ak121396/Desktop/ref-1.txt /home/ak121396/Desktop/ex.txt`) )
# run(pipeline(`./indicators_linux/hyp_ind ./indicators_linux/hyp_ind_param.txt /home/ak121396/Desktop/data.txt
#         /home/ak121396/Desktop/ref-1.txt /home/ak121396/Desktop/performance/ex.txt`) )

readdlm("/home/ak121396/Desktop/ex.txt")[1]


run(pipeline(`./indicators_linux/eps_ind ./indicators_linux/eps_ind_param.txt
    /home/ak121396/Desktop/data.txt /home/ak121396/Desktop/ref-1.txt /home/ak121396/Desktop/ex.txt`) )


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
