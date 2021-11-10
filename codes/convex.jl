# kpp = "F:/results/KS/AP/"; gpp = "F:/results/gpr/AP/"; fpp = "F:/results/fpbh/AP/"
# kpp = "F:/results/KS/KP/"; gpp = "F:/results/gpr/KP/"; fpp = "F:/results/fpbh/KP/"
# kpp = "F:/results/KS/FLP/"; gpp = "F:/results/gpr/FLP/"; fpp = "F:/results/fpbh/FLP/"
using LazySets,Polyhedra,DelimitedFiles,StatsBase

################################################################################
ksp = "/home/ak121396/Desktop/solvers/Kirlikoutput/AP&KP/intAP_Y/"
fpp = "/home/ak121396/Desktop/clu/2minap/"

ksp = "/home/ak121396/Desktop/solvers/Kirlikoutput/AP&KP/intKP_Y/"
bep = "/home/ak121396/Desktop/solvers/Bensolve/KPoutputs/Y/"
gpp = "/home/ak121396/Desktop/GeneralPR/goutputs/KP/GLPK/"
fpp = "/home/ak121396/Desktop/FPBH/KP/GLPK/"

ksp = "/home/ak121396/Desktop/solvers/Kirlikoutput/FLP/ndf/"
gpp = "/home/ak121396/Desktop/GeneralPR/goutputs/FLP/GLPK/"
fpp = "/home/ak121396/Desktop/FPBH/FLP/GLPK/"

gpp = "/media/ak121396//0526-8445/results/gpr/FLP/2/"
fpp = "/media/ak121396//0526-8445/results/fpbh/FLP/2/"


kfiles = readdir(ksp); ffile = readdir(fpp) #gfile = readdir(gpp);
nsk = []; nsg = [];nsf = []
for k=1:length(kfiles)
    @show k
    ks = readdlm(ksp*kfiles[k]); fp = readdlm(fpp*ffile[k]) #gpr = readdlm(gpp*gfile[k]);
    kv = [ ks[i,:] for i=1:size(ks)[1] ]; fv = [ fp[i,:] for i=1:size(fp)[1] ] #gv = [ gpr[i,:] for i=1:size(gpr)[1] ];
    P = VPolytope(kv)
    constraints = constraints_list(P)
    nonk = 0;  nong = 0; nonf = 0

    for i=1:length(fv)
        for c in constraints
            if fv[i] ∈ HyperPlane(c.a,c.b)
                @goto jump3
            end
        end
        nonf+=1
        # push!(nsf,x[i,:]);
        @label jump3
    end
    push!(nsk,size(fp)[1])
    push!(nsf,nonf)

    # for i=1:length(kv)
    #     for c in constraints
    #         if kv[i] ∈ HyperPlane(c.a,c.b)
    #             @goto jump1
    #         end
    #     end
    #     nonk+=1
    #     @label jump1
    # end
    # push!(nsk,nonk)
    # for i=1:length(gv)
    #     for c in constraints
    #         if gv[i] ∈ HyperPlane(c.a,c.b)
    #             @goto jump2
    #         end
    #     end
    #     nong+=1
    #     @label jump2
    # end
    # push!(nsg,nong)

end
[mean(nsk[(i-1)*10+1:10*i,:]) for i=1:10]
[mean(nsf[(i-1)*10+1:10*i,:]) for i=1:10]



# nsg = zeros(120,5); nsf = zeros(120,5);
# gfile = readdir(gpp*"/$i/")
# ffile = readdir(fpp*"/$i/")
# ks = readdlm(kpp*kfiles[i])
# gpr = readdlm(gpp*"/$i/"*gfile[j])
# fp = readdlm(fpp*"/$i/"*ffile[j])

nsg = []; nsf = []
for k=1:length(kfiles)
    gpr = readdlm(gpp*gfile[k]);fp = readdlm(fpp*ffile[k])
    push!(nsg,size(gpr)[1]);push!(nsf,size(fp)[1]);
end

[mean(nsg[(i-1)*10+1:10*i,:]) for i=1:10]

[mean(nsk[(i-1)*10+1:10*i,:]) for i=1:12]

[mean(nsf[(i-1)*10+1:10*i,:]) for i=1:12]

[mean(nsg[(i-1)*10+1:10*i,:]) for i=1:12]













# hull = convex_hull(kv)
# kvv = setdiff(kv,hull)

# remove_redundant_constraints!(constraints)
