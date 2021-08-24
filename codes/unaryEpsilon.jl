using DelimitedFiles,DataFrames
#############################    WINDOWS    ##################################

open("C:/Users/AK121396/Desktop/performance/ex3.txt","w") do io
end

cd("C:\\Users\\AK121396\\Downloads\\performance_indi\\")


run(pipeline(`./indicators_win/hyp_ind ./indicators_win/hyp_ind_param.txt
    F:/results/GPR/KP/1/n-010_ins-01kpY.log
    F:/results/KS/KP/KP_p-3_n-010_ins-01.txt
    C:/Users/AK121396/Desktop/performance/ex.txt`) )

F:/results/GPR/KP/1/n-010_ins-01kpY.log
F:/results/KS/KP/KP_p-3_n-010_ins-01.txt
C:/Users/AK121396/Desktop/performance/ex.txt

######################### Linux #######################################
cd("/home/ak121396//Downloads/performance_indi/indicators_linux/")
run(pipeline(`cd /home/ak121396/Desktop/performance/`))
run(`cd /home/ak121396/Desktop/performance/`)

run(pipeline(`cat '>' /home/ak121396/Desktop/performance/ex.txt`))

cd("/home/ak121396/Downloads/distribution_final/")
# SAMPLE
run(pipeline(`./indicators_linux/eps_ind ./indicators_linux/eps_ind_param.txt /home/ak121396/Desktop/data.txt
        /home/ak121396/Desktop/ref-1.txt /home/ak121396/Desktop/performance/ex.txt`) )
# run(pipeline(`./indicators_linux/hyp_ind ./indicators_linux/hyp_ind_param.txt /home/ak121396/Desktop/data.txt
#         /home/ak121396/Desktop/ref-1.txt /home/ak121396/Desktop/performance/ex.txt`) )
run(pipeline(`./indicators_linux/eps_ind ./indicators_linux/eps_ind_param.txt /home/ak121396/Desktop/GeneralPR/goutputs/AP/GLPK/n-05_ins-01apY.log
        /home/ak121396/Desktop/solvers/Kirlikoutput/AP'&'KP/intAP_Y/AP_p-3_n-05_ins-01.txt /home/ak121396/Desktop/performance/ex.txt`) )
# /home/ak121396/Desktop/FPBH/AP/GLPK/AP_p-3_n-05_ins-01.txt
run(pipeline(`./indicators_linux/eps_ind ./indicators_linux/eps_ind_param.txt /home/ak121396/Desktop/FPBH/AP/GLPK/AP_p-3_n-05_ins-01.txt
        /home/ak121396/Desktop/GeneralPR/goutputs/AP/GLPK/n-05_ins-01apY.log /home/ak121396/Desktop/performance/ex.txt`) )


run(pipeline(`./indicators_linux/hyp_ind ./indicators_linux/hyp_ind_param.txt /home/ak121396/Desktop/data.txt
        /home/ak121396/Desktop/ref-1.txt /home/ak121396/Desktop/performance/ex.txt`) )



readdlm("/home/ak121396/Desktop/solvers/Kirlikoutput/AP&KP/intKP_Y/KP_p-3_n-010_ins-02.ndf")



########################## obj file converter: merging FPBH&GPR ################
fdir = "F:/results/fpbh/MIPLIB/"
gdir = "F:/results/gpr/MIPLIB/"
for i=6:10
    gname = readdir(gdir*"/$i/"); fname = readdir(fdir*"/$i/")
    for j=1:5
        F = readdlm(fdir*"/$i/"*fname[j])
        G = readdlm(gdir*"/$i/"*gname[j])
        P = vcat(F,G)
        ins = fname[j][1:end-4]
        CSV.write("F:/results/mergedMIP/$i"*"/$ins"*".txt",DataFrame(P, :auto),header=false, delim=' ' )
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
