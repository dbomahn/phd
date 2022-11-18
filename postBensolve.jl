lpmodel = loadlp(file*".lp")
# write_to_file(scnd , "/home/k2g00/k2g3475/scnd/lp/"*file[36:end]*".lp")
# lpmodel = loadlp("/home/k2g00/k2g3475/scnd/lp/"*file[36:end]*".lp")
# varub = MPB.getvarUB(lpmodel)
Bmtx = MPB.getconstrmatrix(lpmodel);
B = Bmtx[3:end,:];P = Bmtx[1:2,:]; vub = MPB.getvarUB(lpmodel)
m,n=size(B)
lb = MPB.getconstrLB(lpmodel)[3:end]
ub = MPB.getconstrUB(lpmodel)[3:end]
RHS = Dict()
for i=1:m
    if ub[i]==Inf
        RHS[i] = lb[i]
    else
        RHS[i] = ub[i]
    end
end
signs = []
for i=1:m
    if ub[i] == Inf
        push!(signs,"l")
    elseif lb[i] == -Inf
        push!(signs,"u")
    else
        push!(signs, "s")
    end
end
nz = count(i->(i!=0),B)
objnz = count(i->(i!=0),P)
obj=size(P)[1]
wholearray=[];
arr=["p vlp min",m,n,nz,obj,objnz]
push!(wholearray,arr)

for i=1:m
   for j=1:n
       if (B[i,j]!=0)
           if (B[i,j]%1) == 0 #if B[i,j] is Int
               push!(wholearray,("a",i,j,Int128(B[i,j])))
           else# B[i,j] is Float
               push!(wholearray,("a",i,j,Float64(B[i,j])))
           end
       end
   end
end
for i=1:obj
   for j=1:n
       if P[i,j]!=0
           push!(wholearray,("o",i,j,P[i,j]))
       end
   end
end
for i=1:m
   push!(wholearray,("i",i,signs[i],RHS[i]))
end

for j=1:n
    if j in bvar
        push!(wholearray,("j",j,"s",fpx[1][j])) #assign FP int var values
    else
        push!(wholearray,("j", j,'l',0))
    end
end
push!(wholearray,"e")

# ins = open("/home/k2g00/k2g3475/scnd/vlp/"*file[36:end]*".vlp","w")
ins = open(file*"fp.vlp","w")
writedlm(ins,wholearray)
close(ins)


sum(mt.C[2,:])

sum([pr.dvar[1][k] for k in bvar])
mt.B[2,:]
dv0 = readdlm("/home/ak121396/Desktop/instances/SCND/test01S2_pre_img_p.sol")

dv1 = round.(dv0; digits=1)
(dv1[2,:])

rdvar = [round.(pr.dvar[i][k] for k in bvar) for i=1:length(pr.dvar)]
bidvar = [[pr.dvar[i][k] for k in bvar] for i=1:length(pr.dvar)]

bidvar[1]
unique(rdvar)
ap=[0.6,0.4,0.8]
dd = round.(ap)
sortperm(abs.(ap-dd),rev=true)
