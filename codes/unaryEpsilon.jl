# cd("/home/ak121396//Downloads/performance_indi/indicators_linux/")
cd("C:\\Users\\AK121396\\Downloads\\performance_indi\\")
pwd()
dir1 = "C:\\Users\\AK121396\\\instances\\"
# The order of inputs: .exe file    parameter file    approximation-set file    ref file      outputfile location
run(pipeline(`./indicators_win/eps_ind ./indicators_win/eps_ind_param.txt C:/Users/AK121396/Desktop/approx.txt C:/Users/AK121396/Desktop/ref.txt C:/Users/AK121396/Desktop/performance/ex.txt`) )

run(pipeline(`./hv -r "2 2 2" 1 1 1`))

# run(pipeline(`cat ">>" C:\\Users\\AK121396\\Desktop\\performance\\ex.txt `))
# run((pipeline`cat F:\\results\\gpr\\2ap\\n-05_ins-01apY.log ">>" C:\\Users\\AK121396\\Desktop\\performance\\ex.txt`))

######################### Linux #######################################
run(pipeline(`cd /home/ak121396/Desktop/performance/`))
run(`cd /home/ak121396/Desktop/performance/`)

run(pipeline(`cat '>>'/home/ak121396/Desktop/performance/ex.txt`))

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
