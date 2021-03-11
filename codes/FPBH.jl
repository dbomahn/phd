

using Modof, JuMP, FPBHCPLEX, CPUTime
# print(ARGS)
# print[end-11:end-3])
warmup_fpbhcplex(threads=1)

runtime = @CPUelapsed solutions = fpbhcplex(ARGS[1],[:Min, :Min], threads=1)

runtime = @CPUelapsed solutions = fpbhcplex("/home/k2g00/k2g3475/multiobjective/instances/triFLP/LP/05_010_01.lp",[:Min, :Min], threads=1)
@show ins = ARGS[1][end-11:end-3]
writecsv("/home/k2g00/k2g3475/multiobjective/solvers/FPBH/outputs/time/"*"$ins"*"_time.txt",runtime )#, delim=',' )
write_nondominated_frontier(solutions, "/home/k2g00/k2g3475/multiobjective/solvers/FPBH/outputs/Y/"*"$ins"*".ndf")
write_nondominated_sols(solutions, "/home/k2g00/k2g3475/multiobjective/solvers/FPBH/outputs/X/"*"$ins"*".sol")


# u=1
# while u<11
#     for i=1:100#length(lpmodel)
#         ins = lpmodel[i][1:end-7]
#         runtime = @CPUelapsed solutions = fpbh(d1*lpmodel[i],[:Min, :Min], threads=1)
#
#         writecsv(d1*"/time/"*"$u/"*"$ins"*"_time.txt",runtime )#, delim=',' )
#         write_nondominated_frontier(solutions, d1*"$u"*"ndf/"*lpmodel[i][1:end-7]*".ndf")
#         write_nondominated_sols(solutions, d1*"/sol/"*lpmodel[i][1:end-7]*".sol")
#     end
#     global u=u+1
# end



# @CPUelapsed solutions = fpbh(d1*lpmodel[100],[:Min, :Min], timelimit=1000000000000.0 )
# write_nondominated_frontier(solutions, d1*"/ndf/"*lpmodel[1][1:9]*"Inf.ndf")
# write_nondominated_sols(solutions, d1*"/sol/"*lpmodel[i][1:9]*"Inf.sol")
#
#
# solutions[1].vars
# solutions[1].obj_vals
# all(i->i==0,solutions[1].vars)
#
# wrap_sols_into_array(solutions)


# model = ModoModel()
# @variable(model, x[1:2], Bin)
# @variable(model, y[1:2] >= 0.0)
# objective!(model, 1, :Min, x[1] + x[2] + y[1] + y[2])
# objective!(model, 2, :Max, x[1] + x[2] + y[1] + y[2])
# objective!(model, 3, :Min, x[1] + 2x[2] + y[1] + 2y[2])
# @constraint(model, x[1] + x[2] <= 1)
# @constraint(model, y[1] + 2y[2] >= 1)
# @time solutions = fpbh(model, timelimit=10.0)
