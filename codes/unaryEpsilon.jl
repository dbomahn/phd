# cd("/home/ak121396//Downloads/performance_indi/indicators_linux/")
cd("C:\\Users\\AK121396\\Downloads\\performance_indi\\")
pwd()

# The order of inputs: .exe file    parameter file    approximation-set file    ref file      outputfile location
run(pipeline(`./indicators_win/eps_ind ./indicators_win/eps_ind_param.txt C:/Users/AK121396/Desktop/approx.txt C:/Users/AK121396/Desktop/ref.txt C:/Users/AK121396/Desktop/ue.txt`) )

run(pipeline(`./hv -r "2 2 2" 1 1 1`))

run(mycommand)
