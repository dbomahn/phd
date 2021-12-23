# kpp = "F:/results/KS/AP/"; gpp = "F:/results/gpr/AP/"; fpp = "F:/results/fpbh/AP/"
# kpp = "F:/results/KS/KP/"; gpp = "F:/results/gpr/KP/"; fpp = "F:/results/fpbh/KP/"
# kpp = "F:/results/KS/FLP/"; gpp = "F:/results/gpr/FLP/"; fpp = "F:/results/fpbh/FLP/"
using LazySets,Polyhedra,DelimitedFiles,StatsBase

################################################################################
ksp = "/home/ak121396/Desktop/solvers/Kirlikoutput/AP&KP/intAP_Y/"
wfp = "/home/ak121396/Desktop/clu/AP/"

ksp = "/home/ak121396/Desktop/solvers/Kirlikoutput/AP&KP/intKP_Y/"
wfp = "/home/ak121396/Desktop/clu/KP/"

ksp = "/home/ak121396/Desktop/solvers/Kirlikoutput/FLP/ndf/"
wfp = "/home/ak121396/Desktop/clu/FLP/"

gpp = "/home/ak121396/Desktop/GeneralPR/goutputs/FLP/GLPK/"
fpp = "/media/ak121396/0526-8445/results/fpbh/AP/"

kfiles = readdir(ksp)#[91:100]
# wfiles = readdir(wfp)#[6:15]
ffiles = readdir(fpp)

function nspoints(ins,subclas,ksp,kfiles,wfp,wfiles)
    rep = readdir(wfp)
    allns = zeros(subclas,length(rep));  alltc = zeros(subclas,length(rep));
    for r=1:5
        @show r
        wfiles = readdir(wfp*r)
        tb = zeros(ins,2);
        for k=1:length(kfiles)
            ks = readdlm(ksp*kfiles[k]); kv = [ ks[i,:] for i=1:size(ks)[1] ];
            fp = readdlm(wfp*rep[r]*"/"*wfiles[k]); fv = [ fp[i,:] for i=1:size(fp)[1] ]
            P = VPolytope(kv)
            constraints = constraints_list(P)
            nonf = 0
            for i=1:length(fv)
                for c in constraints
                    if fv[i] ∈ HyperPlane(c.a,c.b)
                        @goto jump3
                    end
                end
                nonf+=1
                @label jump3
            end
            tb[k,1] = nonf; tb[k,2]=length(fv)
        end
        allns[:,r] = [mean(tb[(i-1)*10+1:10*i,1]) for i=1:subclas]
        alltc[:,r] = [mean(tb[(i-1)*10+1:10*i,2]) for i=1:subclas]
    end
    return allns,alltc
end

ns,tc = nspoints(100,10,ksp,kfiles,fpp,ffiles)

round.([mean(ns[i,:]) for i=1:10]; digits=1)
round.([mean(tc[i,:]) for i=1:10]; digits=1)



1
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

#gpr = readdlm(gpp*gfile[k]); #gv = [ gpr[i,:] for i=1:size(gpr)[1] ];

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

[mean(tc[(i-1)*10+1:10*i,:]) for i=1:12]

[mean(nsg[(i-1)*10+1:10*i,:]) for i=1:12]













# hull = convex_hull(kv)
# kvv = setdiff(kv,hull)

# remove_redundant_constraints!(constraints)
