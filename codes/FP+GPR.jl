include("./FFPGPRfunctions.jl")

################################  Data  ####################################
data = Data(ARGS[1],ARGS[2]); pre = Valu(ARGS[3],ARGS[4]);
######################### Mathematical Model #############################
mip = Model(GLPK.Optimizer)
@variable(mip, mx[1:data.n], Bin)
for k=1:data.m
    if data.signs[k] == "l"
        @constraint(mip, dot(data.B[k,:],mx) >= data.RHS[k])
    elseif data.signs[k] == "u"
        @constraint(mip, dot(data.B[k,:],mx) <= data.RHS[k])
    else
        @constraint(mip, dot(data.B[k,:],mx) == data.RHS[k])
    end
end
optimize!(mip);
#####################  Feasibility Search Mathematical Model  ##################
dist = Model(GLPK.Optimizer)
@variables(dist, begin
    0 <= x[1:data.n] <=1
end)
for k=1:data.m
    if data.signs[k] == "l"
        @constraint(dist, dot(data.B[k,:],x) >= data.RHS[k])
    elseif data.signs[k] == "u"
        @constraint(dist, dot(data.B[k,:],x) <= data.RHS[k])
    else
        @constraint(dist, dot(data.B[k,:],x) == data.RHS[k])
    end
end
optimize!(dist);
######################## Running Feasibility Pump ##########################
Bentime = readdlm(ARGS[5])[1]; FPTL = (TL-Bentime)/2;
weightFP(pre.dvar,data.n,data.C,5)
runtime1 = @CPUelapsed fpX,newsol = fractionalFP(pre.dvar,data.n,data.C,FPTL)

# fpX,fpY = domFilter(X2,PF2);
# clistY = []
# for i=1:length(candlist)
# 	push!(clistY,getobjval(candlist[1],data.C))
# end
# cX = [fpX;candlist]; cY = [fpY;clistY]
tx = copy(fpX); ty = copy(fpY);

######################## Running Generic Path Relinking ##########################
GPR(data.C,data.n,tx,ty,5) #compiling
runtime1 = @CPUelapsed X,nsol = fractionalFP(pre.dvar,data.n,data.C,TL)
runtime2 = @CPUelapsed candset,candobj,newsol2 = GPR(data.C,data.n,fpX,fpY,TL-FPtime)
totaltime = runtime1+Bentime+runtime2
finalX,finalY = domFilter(candset,candobj);
otable = zeros(length(finalY),3)
for i=1:length(finalY)
    for j=1:3
        otable[i,j] = finalY[i][j]
    end
end
newsol3 = length(setdiff(finalY,fpY))

ins = ARGS[2][1:end-4];
print("$ins", totaltime)
CSV.write(ins*"_1Y.log",DataFrame(otable, :auto),append=false, header=false, delim=' ' )
#record1 = DataFrame(Filename=ins[48:end],FPsol=length(fpY),PRsol=newsol2, totalsol=length(finalY)) #,Bentime=Bentime,FPtime=FPtime,PRtime=GPRtime, totaltime=totaltime )
CSV.write("/home/k2g00/k2g3475/clusterhome/multiobjective/generalPR/record.csv", record1,append=true, header=false )#, delim=','i )
# matriX = zeros(Int,length(finalX),data.n)
# for i=1:length(finalX)
#     for j=1:data.n
#         matriX[i,j] = finalX[i][j]
#     end
# end
