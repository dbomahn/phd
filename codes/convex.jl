using LazySets,Polyhedra,DelimitedFiles,StatsBase
kpp = "F:/results/KS/AP/"; gpp = "F:/results/gpr/AP/"; fpp = "F:/results/fpbh/AP/"
kpp = "F:/results/KS/KP/"; gpp = "F:/results/gpr/KP/"; fpp = "F:/results/fpbh/KP/"
kpp = "F:/results/KS/FLP/"; gpp = "F:/results/gpr/FLP/"; fpp = "F:/results/fpbh/FLP/"


kfiles = readdir(kpp)
nsg = zeros(100,1); nsf = zeros(100,1);
# nsg = zeros(120,5); nsf = zeros(120,5);

for i=1:1
    gfile = readdir(gpp*"/$i/"); ffile = readdir(fpp*"/$i/");
    for j=1:100#120
        ks = readdlm(kpp*kfiles[j])
        gpr = readdlm(gpp*"/$i/"*gfile[j])
        fp = readdlm(fpp*"/$i/"*ffile[j])

        kv = [ ks[k,:] for k=1:size(ks)[1] ]
        gv = [ gpr[k,:] for k=1:size(gpr)[1] ]
        fv = [ fp[k,:] for k=1:size(fp)[1] ]

        hull = convex_hull(kv); P = VPolytope(hull)
        ctg = 0; ctf = 0;
        for i=1:length(gv)
            if gv[i] ∉ P
                ctg+=1
            end
        end

        for i=1:length(fv)
            if fv[i] ∉ P
                ctf+=1
            end
        end
        nsg[j,i] = ctg; nsf[j,i] = ctf;
    end
end

# gfile = readdir(gpp*"/$i/")
# ffile = readdir(fpp*"/$i/")
# ks = readdlm(kpp*kfiles[i])
# gpr = readdlm(gpp*"/$i/"*gfile[j])
# fp = readdlm(fpp*"/$i/"*ffile[j])
# nsf
# nsg
[mean(nsf[(i-1)*10+1:10*i,:]) for i=1:10]
[mean(nsg[(i-1)*10+1:10*i,:]) for i=1:10]

[mean(nsf[(i-1)*10+1:10*i,:]) for i=1:12]
[mean(nsg[(i-1)*10+1:10*i,:]) for i=1:12]

################################################################################

kpp = "F:/results/KS/AP/"; gpp = "F:/results/gpr/AP/5/"; fpp = "F:/results/fpbh/AP/5/"
kfiles = readdir(kpp); gfile = readdir(gpp); ffile = readdir(fpp);

i=j=2;
ks = readdlm(kpp*kfiles[i])
gpr = readdlm(gpp*gfile[j])
fp = readdlm(fpp*ffile[j])

kv = [ ks[i,:] for i=1:size(ks)[1] ]
gv = [ gpr[i,:] for i=1:size(gpr)[1] ]
fv = [ fp[i,:] for i=1:size(fp)[1] ]

hull = convex_hull(kv); P = VPolytope(hull)
ctg = 0; ctf = 0;
# ctg = []; ctf = [];
for i=1:length(gv)
    if gv[i] ∉ P
        # push!(ctg, gv[i])
        ctg+=1
    end
end

for i=1:length(fv)
    if fv[i] ∉ P
        ctf+=1
        # push!(ctf, fv[i])
    end
end
ctg
ctf

#################################

kf = readdir(ksp)
non = [];
for i=1:length(kf)
    @show i; ct=0;
    ks = readdlm(ksp*kf[i])
    kv = [ks[k,:] for k=1:size(ks)[1]]
    hull = convex_hull(kv); P = VPolytope(hull)
    kvv = setdiff(kv,hull)
    for l=1:length(kvv)
        if kvv[l] ∉ P
            ct+=1
        end
    end
    push!(non,ct)
end



# gp = readdlm(gpp*gf[1])
for l=1:length(kv)
    if kv[l] ∉ P
        println(kv[l])
    end
end
