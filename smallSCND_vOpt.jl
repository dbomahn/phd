using CPUTime,DataFrames,DelimitedFiles,JuMP,LinearAlgebra,CPLEX,vOptGeneric,SparseArrays,StatsBase

################# 1dim model Data & Model ##############
c2 = []
for i=1:J+K
        for j=1:2
                push!(c2,fc[i][j])
        end
end
e2 = hcat(reshape(vep',1,J*2),reshape(ved',1,K*2))
fsp,fpd,fdc = [],[],[];
for i=1:I
    for j=1:J
        m = Mij[i,j]
        if m==1
            push!(fsp,10000)
        else
            push!(fsp,10000); push!(fsp,0);
        end
    end
end
for j=1:J
    for k=1:K
        m = Mjk[j,k]
        if m==1
            push!(fpd,10000)
        else
            push!(fpd,10000);push!(fpd,0);
        end
    end
end
for k=1:K
    for l=1:L
        m = Mkl[k,l]
        if m==1
            push!(fdc,10000)
        else
            push!(fdc,10000);push!(fdc,0);
        end
    end
end
gij2 = sparse(fsp); gjk2 = sparse(fpd); gkl2 = sparse(fdc)

tcp2 = []
for i=1:length(Ipt)
    for j=1:length(Jpt)
        if Mij[i,j]==1
            push!(tcp2, tcp[i][j][1])
        else
            push!(tcp2, tcp[i][j][1],tcp[i][j][2])
        end
    end
end

tcd2 = []
for i=1:length(Jpt)
    for j=1:length(Kpt)
        if Mjk[i,j]==1
            push!(tcd2, tcd[i][j][1])
        else
            push!(tcd2, tcd[i][j][1],tcd[i][j][2])
        end
    end
end

tcc2 = []
for i=1:length(Kpt)
    for j=1:length(Lpt)
        if Mkl[i,j]==1
            push!(tcc2, tcc[i][j][1])
        else
            push!(tcc2, tcc[i][j][1],tcc[i][j][2])
        end
    end
end

q2 = hcat(reshape(vcp',1,J*2),reshape(vcd',1,K*2))

rij2 = []
for i=1:length(Ipt)
    for j=1:length(Jpt)
        if Mij[i,j]==1
            push!(rij2, cep[i][j][1])
        else
            push!(rij2, cep[i][j][1],cep[i][j][2])
        end
    end
end

rjk2 = []
for i=1:length(Jpt)
    for j=1:length(Kpt)
        if Mjk[i,j]==1
            push!(rjk2, ced[i][j][1])
        else
            push!(rjk2, ced[i][j][1],ced[i][j][2])
        end
    end
end

rkl2 = []
for i=1:length(Kpt)
    for j=1:length(Lpt)
        if Mkl[i,j]==1
            push!(rkl2, cec[i][j][1])
        else
            push!(rkl2, cec[i][j][1],cec[i][j][2])
        end
    end
end

lsp,lpd,ldc = [],[],[]

for i=1:I
    for j=1:J
        if Mij[i,j]==1
            push!(lsp,0)
        else
            push!(lsp,0); push!(lsp,maximum(maximum(maximum(Lcapasp))));
        end
    end
end
for j=1:J
    for k=1:K
        if Mjk[j,k]==1
            push!(lpd,0)
        else
            push!(lpd,0); push!(lpd,maximum(maximum(maximum(Lcapapd))));
        end
    end
end
Vij2 = sparse(lsp); Vjk2 = sparse(lpd);
# for k=1:K
#     for l=1:L
#         if Mkl[k,l]==1
#             push!(ldc,0)
#         else
#             push!(ldc,0); push!(ldc,maximum(maximum(maximum(Lcapadc))));
#         end
#     end
# end


function smscnd_lp()
    smscnd = vModel(optimizer_with_attributes(
            CPLEX.Optimizer,
            "CPX_PARAM_EPGAP" => 1e-8
          ));
    set_silent(smscnd)
    MOI.set(smscnd, MOI.NumberOfThreads(), 1)
    @variable(smscnd, 0<= y[1:(J+K)*2]<=1)
    @variable(smscnd, 0<=uij[1:sum(Mij)]<=1);
    @variable(smscnd, 0<=ujk[1:sum(Mjk)]<=1);
    @variable(smscnd, 0<=ukl[1:sum(Mkl)]<=1);
    @variable( smscnd, 0<= xij[1:sum(Mij)] );
    @variable( smscnd, 0<= xjk[1:sum(Mjk)] );
    @variable( smscnd, 0<= xkl[1:sum(Mkl)] );
    @variable( smscnd, 0<= h[1:(J+K)*2] );

    @addobjective(smscnd, Min, sum(c2.*y) +sum(repeat(vcs[1,:], outer=sum(Mij[1,:])).*xij[1:sum(Mij[1,:])]) +
        sum(sum(repeat(vcs[i,:], outer=sum(Mij[i,:])).*xij[sum(Mij[1:i-1,:])+1:sum(Mij[1:i,:])]) for i=2:I) +
        dot(e2,h) +sum(gij2[i]*uij[i] for i in findnz(gij2)[1]) + sum(gjk2[i]*ujk[i] for i in findnz(gjk2)[1]) + sum(gkl2[i].*ukl[i] for i in findnz(gkl2)[1]) +
        sum(tcp2.*xij)+sum(tcd2.*xjk)+sum(tcc2.*xkl));
    @addobjective(smscnd, Min, sum(repeat(ves[1,:], outer=sum(Mij[1,:])).*xij[1:sum(Mij[1,:])]) +
        sum(sum(repeat(ves[i,:], outer=sum(Mij[i,:])).*xij[sum(Mij[1:i-1,:])+1:sum(Mij[1:i,:])]) for i=2:I) +
        dot(q2,h) + sum(rij2.*xij)+sum(rjk2.*xjk)+sum(rkl2.*xkl));
    ########## constraint 3 #############
    @constraint(smscnd, sum(xij[m] for m=1:Mij[1,1])+ sum(xij[m+(sum(Mij[1:i-1,:]))] for i=2:I for m=1:Mij[i,1])  == sum(xjk[m] for m=1:sum(Mjk[1,:])) )
    @constraint(smscnd, [j=2:J], sum(xij[sum(Mij[1,1:j-1])+m] for m=1:Mij[1,j])+sum(xij[sum(Mij[i,1:j-1])+m+(sum(Mij[1:i-1,:]))] for i=2:I for m=1:Mij[i,j]) == sum(xjk[sum(Mjk[1:j-1,:]) + m] for m=1:sum(Mjk[j,:])) )
    @constraint(smscnd, sum(xjk[m] for m=1:Mjk[1,1])+sum(xjk[m+(sum(Mjk[1:j-1,:]))] for j=2:J for m=1:Mjk[j,1]) == sum(xkl[m] for m=1:sum(Mkl[1,:])) )
    @constraint(smscnd, [k=2:K],sum(xjk[sum(Mjk[1,1:k-1])+m] for m=1:Mjk[1,k])+sum(xjk[sum(Mjk[j,1:k-1])+m+(sum(Mjk[1:j-1,:]))] for j=2:J for m=1:Mjk[j,k]) == sum(xkl[sum(Mkl[1:k-1,:]) + m] for m=1:sum(Mkl[k,:])) )
    ########### constraint 4-6 #############
    @constraint(smscnd, sum(h[t] for t=1:2) ==sum(xij[m] for m=1:Mij[1,1])+sum(xij[m+(sum(Mij[1:i-1,:]))] for i=2:I for m=1:Mij[i,1]));
    @constraint(smscnd, [j=2:J], sum(h[2*(j-1)+t] for t=1:2) == sum(xij[sum(Mij[1,1:j-1])+m] for m=1:Mij[1,j])+sum(xij[sum(Mij[i,1:j-1])+m+(sum(Mij[1:i-1,:]))] for i=2:I for m=1:Mij[i,j]) );
    @constraint(smscnd,sum(h[2*J+t] for t=1:2) == sum(xjk[m] for m=1:Mjk[1,1])+sum(xjk[m+(sum(Mjk[1:j-1,:]))] for j=2:J for m=1:Mjk[j,1]) );
    @constraint(smscnd, [k=2:K], sum(h[2*J+2*(k-1)+t] for t=1:2) == sum(xjk[sum(Mjk[1,1:k-1])+m] for m=1:Mjk[1,k])+sum(xjk[sum(Mjk[j,1:k-1])+m+(sum(Mjk[1:j-1,:]))] for j=2:J for m=1:Mjk[j,k]));
    @constraint(smscnd, sum(xkl[m] for m=1:Mkl[1,1]) +sum(xkl[m+(sum(Mkl[1:k-1,:]))] for k=2:K for m=1:Mkl[k,1]) >= demand[1]);
    @constraint(smscnd, [l=2:L], sum(xkl[sum(Mkl[1,1:l-1]) + m] for m=1:Mkl[1,l])+ sum(xkl[sum(Mkl[1:k-1,:])+sum(Mkl[k,1:l-1])+m] for k=2:K for m=1:Mkl[k,l]) >= demand[l]);
    ########### constraint 7 #############
    @constraint(smscnd, sum(xij[1:sum(Mij[1,:])]) <= cas1[1]);
    @constraint(smscnd, [i=2:I],  sum(xij[sum(Mij[1:i-1,:])+1:sum(Mij[1:i,:])]) <= cas1[i]);
    ########### constraint 8 #############
    @constraint(smscnd,[j=1:J+K, t=1:2], sum(h[2*(j-1)+t]) <= capd1[j]*y[2*(j-1)+t]);
    ########### constraint 9 #############
    @constraint(smscnd,[j=1:J+K], sum(y[2*(j-1)+1:2*(j-1)+2]) <= 1);
    ########### constraint 10 #############
    @constraints(smscnd,begin
            sum(uij[1:Mij[1,1]]) <= 1
            [j=2:J], sum(uij[sum(Mij[1,1:j-1])+1:sum(Mij[1,1:j-1])+Mij[1,j]]) <= 1
            [i=2:I,j=2:J],  sum(uij[sum(Mij[1:i-1,:])+sum(Mij[i,1:j-1])+1:sum(Mij[1:i-1,:])+sum(Mij[i,1:j-1])+Mij[i,j]])<= 1
            sum(ujk[1:Mjk[1,1]]) <= 1
            [k=2:K], sum(ujk[sum(Mjk[1,1:k-1])+1:sum(Mjk[1,1:k-1])+Mjk[1,k]]) <= 1
            [j=2:J,k=2:K],  sum(ujk[sum(Mjk[1:j-1,:])+sum(Mjk[j,1:j-1])+1:sum(Mjk[1:j-1,:])+sum(Mjk[j,1:j-1])+Mjk[j,k]]) <= 1
            sum(ukl[1:Mkl[1,1]]) <= 1
            [l=2:L], sum(ukl[sum(Mkl[1,1:l-1])+1:sum(Mkl[1,1:l-1])+Mkl[1,l]]) <= 1
            [k=2:K,l=2:L],  sum(ukl[sum(Mkl[1:k-1,:])+sum(Mkl[k,1:l-1])+1:sum(Mkl[1:k-1,:])+sum(Mkl[k,1:l-1])+Mkl[k,l]])<= 1
    end);
    ########### constraint 11 #############
    @constraints(smscnd, begin
        [i=1:sum(Mij)], sum(xij[(i-1)+1:i]) <= bigM*uij[i]
        [i=1:sum(Mjk)], sum(xjk[(i-1)+1:i]) <= bigM*ujk[i]
        [i=1:sum(Mkl)], sum(xkl[(i-1)+1:i]) <= bigM*ukl[i]
    end);
    ########### constraint 12 #############
    @constraints(smscnd, begin
            [i in findnz(Vij2)[1]], sum(xij[(i-1)+1:i]) >= Vij2[i]*uij[i]
            [i in findnz(Vjk2)[1]], sum(xjk[(i-1)+1:i]) >= Vjk2[i]*ujk[i]
            # [i in findnz(Vkl)[1]], sum(xkl[(i-1)+1:i]) >= Vkl[i]*ukl[i]
    end);
    ########### constraint 13-14 #############
    @constraint(smscnd, sum(y[1:J*2]) <= Jmax);
    @constraint(smscnd, sum(y[J*2+1:end]) <= Kmax);
    return smscnd
end
lpm = smscnd_lp()
@CPUtime vSolve(lpm, 50,method=:dicho, verbose=false)
lyn = getvOptData(lpm).Y_N
lxe = getvOptData(lpm).X_E
lweight = round(Int,mean([lyn[i][1]/lyn[i][2] for i=1:length(lyn)]))

################################### Dicho MIP ##############################
function smscnd_1dim()
    smscnd = vModel(optimizer_with_attributes(
            CPLEX.Optimizer,
            "CPX_PARAM_EPGAP" => 1e-8
          ));
    set_silent(smscnd)
    MOI.set(smscnd, MOI.NumberOfThreads(), 1)
    @variable(smscnd, y[1:(J+K)*2], Bin)
    @variable(smscnd, uij[1:sum(Mij)], Bin);
    @variable(smscnd, ujk[1:sum(Mjk)], Bin);
    @variable(smscnd, ukl[1:sum(Mkl)], Bin);
    @variable( smscnd, 0<= xij[1:sum(Mij)] );
    @variable( smscnd, 0<= xjk[1:sum(Mjk)] );
    @variable( smscnd, 0<= xkl[1:sum(Mkl)] );
    @variable( smscnd, 0<= h[1:(J+K)*2] );

    @addobjective(smscnd, Min, sum(c2.*y) +sum(repeat(vcs[1,:], outer=sum(Mij[1,:])).*xij[1:sum(Mij[1,:])]) +
        sum(sum(repeat(vcs[i,:], outer=sum(Mij[i,:])).*xij[sum(Mij[1:i-1,:])+1:sum(Mij[1:i,:])]) for i=2:I) +
        dot(e2,h) +sum(gij2[i]*uij[i] for i in findnz(gij2)[1]) + sum(gjk2[i]*ujk[i] for i in findnz(gjk2)[1]) + sum(gkl2[i].*ukl[i] for i in findnz(gkl2)[1]) +
        sum(tcp2.*xij)+sum(tcd2.*xjk)+sum(tcc2.*xkl));
    @addobjective(smscnd, Min, sum(repeat(ves[1,:], outer=sum(Mij[1,:])).*xij[1:sum(Mij[1,:])]) +
        sum(sum(repeat(ves[i,:], outer=sum(Mij[i,:])).*xij[sum(Mij[1:i-1,:])+1:sum(Mij[1:i,:])]) for i=2:I) +
        dot(q2,h) + sum(rij2.*xij)+sum(rjk2.*xjk)+sum(rkl2.*xkl));
    ########## constraint 3 #############
    @constraint(smscnd, sum(xij[m] for m=1:Mij[1,1])+ sum(xij[m+(sum(Mij[1:i-1,:]))] for i=2:I for m=1:Mij[i,1])  == sum(xjk[m] for m=1:sum(Mjk[1,:])) )
    @constraint(smscnd, [j=2:J], sum(xij[sum(Mij[1,1:j-1])+m] for m=1:Mij[1,j])+sum(xij[sum(Mij[i,1:j-1])+m+(sum(Mij[1:i-1,:]))] for i=2:I for m=1:Mij[i,j]) == sum(xjk[sum(Mjk[1:j-1,:]) + m] for m=1:sum(Mjk[j,:])) )
    @constraint(smscnd, sum(xjk[m] for m=1:Mjk[1,1])+sum(xjk[m+(sum(Mjk[1:j-1,:]))] for j=2:J for m=1:Mjk[j,1]) == sum(xkl[m] for m=1:sum(Mkl[1,:])) )
    @constraint(smscnd, [k=2:K],sum(xjk[sum(Mjk[1,1:k-1])+m] for m=1:Mjk[1,k])+sum(xjk[sum(Mjk[j,1:k-1])+m+(sum(Mjk[1:j-1,:]))] for j=2:J for m=1:Mjk[j,k]) == sum(xkl[sum(Mkl[1:k-1,:]) + m] for m=1:sum(Mkl[k,:])) )
    ########### constraint 4-6 #############
    @constraint(smscnd, sum(h[t] for t=1:2) ==sum(xij[m] for m=1:Mij[1,1])+sum(xij[m+(sum(Mij[1:i-1,:]))] for i=2:I for m=1:Mij[i,1]));
    @constraint(smscnd, [j=2:J], sum(h[2*(j-1)+t] for t=1:2) == sum(xij[sum(Mij[1,1:j-1])+m] for m=1:Mij[1,j])+sum(xij[sum(Mij[i,1:j-1])+m+(sum(Mij[1:i-1,:]))] for i=2:I for m=1:Mij[i,j]) );
    @constraint(smscnd,sum(h[2*J+t] for t=1:2) == sum(xjk[m] for m=1:Mjk[1,1])+sum(xjk[m+(sum(Mjk[1:j-1,:]))] for j=2:J for m=1:Mjk[j,1]) );
    @constraint(smscnd, [k=2:K], sum(h[2*J+2*(k-1)+t] for t=1:2) == sum(xjk[sum(Mjk[1,1:k-1])+m] for m=1:Mjk[1,k])+sum(xjk[sum(Mjk[j,1:k-1])+m+(sum(Mjk[1:j-1,:]))] for j=2:J for m=1:Mjk[j,k]));
    @constraint(smscnd, sum(xkl[m] for m=1:Mkl[1,1]) +sum(xkl[m+(sum(Mkl[1:k-1,:]))] for k=2:K for m=1:Mkl[k,1]) >= demand[1]);
    @constraint(smscnd, [l=2:L], sum(xkl[sum(Mkl[1,1:l-1]) + m] for m=1:Mkl[1,l])+ sum(xkl[sum(Mkl[1:k-1,:])+sum(Mkl[k,1:l-1])+m] for k=2:K for m=1:Mkl[k,l]) >= demand[l]);
    ########### constraint 7 #############
    @constraint(smscnd, sum(xij[1:sum(Mij[1,:])]) <= cas1[1]);
    @constraint(smscnd, [i=2:I],  sum(xij[sum(Mij[1:i-1,:])+1:sum(Mij[1:i,:])]) <= cas1[i]);
    ########### constraint 8 #############
    @constraint(smscnd,[j=1:J+K, t=1:2], sum(h[2*(j-1)+t]) <= capd1[j]*y[2*(j-1)+t]);
    ########### constraint 9 #############
    @constraint(smscnd,[j=1:J+K], sum(y[2*(j-1)+1:2*(j-1)+2]) <= 1);
    ########### constraint 10 #############
    @constraints(smscnd,begin
            sum(uij[1:Mij[1,1]]) <= 1
            [j=2:J], sum(uij[sum(Mij[1,1:j-1])+1:sum(Mij[1,1:j-1])+Mij[1,j]]) <= 1
            [i=2:I,j=2:J],  sum(uij[sum(Mij[1:i-1,:])+sum(Mij[i,1:j-1])+1:sum(Mij[1:i-1,:])+sum(Mij[i,1:j-1])+Mij[i,j]])<= 1
            sum(ujk[1:Mjk[1,1]]) <= 1
            [k=2:K], sum(ujk[sum(Mjk[1,1:k-1])+1:sum(Mjk[1,1:k-1])+Mjk[1,k]]) <= 1
            [j=2:J,k=2:K],  sum(ujk[sum(Mjk[1:j-1,:])+sum(Mjk[j,1:j-1])+1:sum(Mjk[1:j-1,:])+sum(Mjk[j,1:j-1])+Mjk[j,k]]) <= 1
            sum(ukl[1:Mkl[1,1]]) <= 1
            [l=2:L], sum(ukl[sum(Mkl[1,1:l-1])+1:sum(Mkl[1,1:l-1])+Mkl[1,l]]) <= 1
            [k=2:K,l=2:L],  sum(ukl[sum(Mkl[1:k-1,:])+sum(Mkl[k,1:l-1])+1:sum(Mkl[1:k-1,:])+sum(Mkl[k,1:l-1])+Mkl[k,l]])<= 1
    end);
    ########### constraint 11 #############
    @constraints(smscnd, begin
        [i=1:sum(Mij)], sum(xij[(i-1)+1:i]) <= bigM*uij[i]
        [i=1:sum(Mjk)], sum(xjk[(i-1)+1:i]) <= bigM*ujk[i]
        [i=1:sum(Mkl)], sum(xkl[(i-1)+1:i]) <= bigM*ukl[i]
    end);
    ########### constraint 12 #############
    @constraints(smscnd, begin
            [i in findnz(Vij2)[1]], sum(xij[(i-1)+1:i]) >= Vij2[i]*uij[i]
            [i in findnz(Vjk2)[1]], sum(xjk[(i-1)+1:i]) >= Vjk2[i]*ujk[i]
            # [i in findnz(Vkl)[1]], sum(xkl[(i-1)+1:i]) >= Vkl[i]*ukl[i]
    end);
    ########### constraint 13-14 #############
    @constraint(smscnd, sum(y[1:J*2]) <= Jmax);
    @constraint(smscnd, sum(y[J*2+1:end]) <= Kmax);
    return smscnd
end
comp1 = smscnd_1dim()
@CPUtime vSolve(comp1,50,method=:dicho, verbose=false)
myn = getvOptData(comp1).Y_N;
mxe = getvOptData(comp1).X_E;
mweight = round(Int,mean([myn[i][1]/myn[i][2] for i=1:length(myn)]))
################
function FP_Model(weight)
    model = Model(CPLEX.Optimizer); set_silent(model)
    MOI.set(model, MOI.NumberOfThreads(), 1)
    # MOI.set(model, MOI.RawParameter("CPX_PARAM_SCRIND"), false);
    @variable(model, y[1:(J+K)*2], Bin)
    @variable(model, uij[1:sum(Mij)], Bin);
    @variable(model, ujk[1:sum(Mjk)], Bin);
    @variable(model, ukl[1:sum(Mkl)], Bin);
    @variable( model, 0<= xij[1:sum(Mij)] );
    @variable( model, 0<= xjk[1:sum(Mjk)] );
    @variable( model, 0<= xkl[1:sum(Mkl)] );
    @variable( model, 0<= h[1:(J+K)*2] );
    @objective(model, Min,  sum(c2.*y) +sum(repeat(vcs[1,:], outer=sum(Mij[1,:])).*xij[1:sum(Mij[1,:])]) +
        sum(sum(repeat(vcs[i,:], outer=sum(Mij[i,:])).*xij[sum(Mij[1:i-1,:])+1:sum(Mij[1:i,:])]) for i=2:I) +
        dot(e2,h) +sum(gij2[i]*uij[i] for i in findnz(gij2)[1]) + sum(gjk2[i]*ujk[i] for i in findnz(gjk2)[1]) + sum(gkl2[i].*ukl[i] for i in findnz(gkl2)[1]) +
        sum(tcp2.*xij)+sum(tcd2.*xjk)+sum(tcc2.*xkl)    +
            weight*(sum(repeat(ves[1,:], outer=sum(Mij[1,:])).*xij[1:sum(Mij[1,:])]) +
        sum(sum(repeat(ves[i,:], outer=sum(Mij[i,:])).*xij[sum(Mij[1:i-1,:])+1:sum(Mij[1:i,:])]) for i=2:I) +
        dot(q2,h) + sum(rij2.*xij)+sum(rjk2.*xjk)+sum(rkl2.*xkl))
    );
    ########## constraint 3 #############
    @constraint(model, sum(xij[m] for m=1:Mij[1,1])+ sum(xij[m+(sum(Mij[1:i-1,:]))] for i=2:I for m=1:Mij[i,1])  == sum(xjk[m] for m=1:sum(Mjk[1,:])) )
    @constraint(model, [j=2:J], sum(xij[sum(Mij[1,1:j-1])+m] for m=1:Mij[1,j])+sum(xij[sum(Mij[i,1:j-1])+m+(sum(Mij[1:i-1,:]))] for i=2:I for m=1:Mij[i,j]) == sum(xjk[sum(Mjk[1:j-1,:]) + m] for m=1:sum(Mjk[j,:])) )
    @constraint(model, sum(xjk[m] for m=1:Mjk[1,1])+sum(xjk[m+(sum(Mjk[1:j-1,:]))] for j=2:J for m=1:Mjk[j,1]) == sum(xkl[m] for m=1:sum(Mkl[1,:])) )
    @constraint(model, [k=2:K],sum(xjk[sum(Mjk[1,1:k-1])+m] for m=1:Mjk[1,k])+sum(xjk[sum(Mjk[j,1:k-1])+m+(sum(Mjk[1:j-1,:]))] for j=2:J for m=1:Mjk[j,k]) == sum(xkl[sum(Mkl[1:k-1,:]) + m] for m=1:sum(Mkl[k,:])) )
    ########### constraint 4-6 #############
    @constraint(model, sum(h[t] for t=1:2) ==sum(xij[m] for m=1:Mij[1,1])+sum(xij[m+(sum(Mij[1:i-1,:]))] for i=2:I for m=1:Mij[i,1]));
    @constraint(model, [j=2:J], sum(h[2*(j-1)+t] for t=1:2) == sum(xij[sum(Mij[1,1:j-1])+m] for m=1:Mij[1,j])+sum(xij[sum(Mij[i,1:j-1])+m+(sum(Mij[1:i-1,:]))] for i=2:I for m=1:Mij[i,j]) );
    @constraint(model,sum(h[2*J+t] for t=1:2) == sum(xjk[m] for m=1:Mjk[1,1])+sum(xjk[m+(sum(Mjk[1:j-1,:]))] for j=2:J for m=1:Mjk[j,1]) );
    @constraint(model, [k=2:K], sum(h[2*J+2*(k-1)+t] for t=1:2) == sum(xjk[sum(Mjk[1,1:k-1])+m] for m=1:Mjk[1,k])+sum(xjk[sum(Mjk[j,1:k-1])+m+(sum(Mjk[1:j-1,:]))] for j=2:J for m=1:Mjk[j,k]));
    @constraint(model, sum(xkl[m] for m=1:Mkl[1,1]) +sum(xkl[m+(sum(Mkl[1:k-1,:]))] for k=2:K for m=1:Mkl[k,1]) >= demand[1]);
    @constraint(model, [l=2:L], sum(xkl[sum(Mkl[1,1:l-1]) + m] for m=1:Mkl[1,l])+ sum(xkl[sum(Mkl[1:k-1,:])+sum(Mkl[k,1:l-1])+m] for k=2:K for m=1:Mkl[k,l]) >= demand[l]);
    ########### constraint 7 #############
    @constraint(model, sum(xij[1:sum(Mij[1,:])]) <= cas1[1]);
    @constraint(model, [i=2:I],  sum(xij[sum(Mij[1:i-1,:])+1:sum(Mij[1:i,:])]) <= cas1[i]);
    ########### constraint 8 #############
    @constraint(model,[j=1:J+K, t=1:2], sum(h[2*(j-1)+t]) <= capd1[j]*y[2*(j-1)+t]);
    ########### constraint 9 #############
    @constraint(model,[j=1:J+K], sum(y[2*(j-1)+1:2*(j-1)+2]) <= 1);
    ########### constraint 10 #############
    @constraints(model,begin
            sum(uij[1:Mij[1,1]]) <= 1
            [j=2:J], sum(uij[sum(Mij[1,1:j-1])+1:sum(Mij[1,1:j-1])+Mij[1,j]]) <= 1
            [i=2:I,j=2:J],  sum(uij[sum(Mij[1:i-1,:])+sum(Mij[i,1:j-1])+1:sum(Mij[1:i-1,:])+sum(Mij[i,1:j-1])+Mij[i,j]])<= 1
            sum(ujk[1:Mjk[1,1]]) <= 1
            [k=2:K], sum(ujk[sum(Mjk[1,1:k-1])+1:sum(Mjk[1,1:k-1])+Mjk[1,k]]) <= 1
            [j=2:J,k=2:K],  sum(ujk[sum(Mjk[1:j-1,:])+sum(Mjk[j,1:j-1])+1:sum(Mjk[1:j-1,:])+sum(Mjk[j,1:j-1])+Mjk[j,k]]) <= 1
            sum(ukl[1:Mkl[1,1]]) <= 1
            [l=2:L], sum(ukl[sum(Mkl[1,1:l-1])+1:sum(Mkl[1,1:l-1])+Mkl[1,l]]) <= 1
            [k=2:K,l=2:L],  sum(ukl[sum(Mkl[1:k-1,:])+sum(Mkl[k,1:l-1])+1:sum(Mkl[1:k-1,:])+sum(Mkl[k,1:l-1])+Mkl[k,l]])<= 1
    end);
    ########### constraint 11 #############
    @constraints(model, begin
        [i=1:sum(Mij)], sum(xij[(i-1)+1:i]) <= bigM*uij[i]
        [i=1:sum(Mjk)], sum(xjk[(i-1)+1:i]) <= bigM*ujk[i]
        [i=1:sum(Mkl)], sum(xkl[(i-1)+1:i]) <= bigM*ukl[i]
    end);
    ########### constraint 12 #############
    @constraints(model, begin
            [i in findnz(Vij2)[1]], sum(xij[(i-1)+1:i]) >= Vij2[i]*uij[i]
            [i in findnz(Vjk2)[1]], sum(xjk[(i-1)+1:i]) >= Vjk2[i]*ujk[i]
            # [i in findnz(Vkl)[1]], sum(xkl[(i-1)+1:i]) >= Vkl[i]*ukl[i]
    end);
    ########### constraint 13-14 #############
    @constraint(model, sum(y[1:J*2]) <= Jmax);
    @constraint(model, sum(y[J*2+1:end]) <= Kmax);
    return model
end
function LP_Model(weight)
    lp = Model(CPLEX.Optimizer); set_silent(lp)
    MOI.set(lp, MOI.NumberOfThreads(), 1)
    # MOI.set(lp, MOI.RawParameter("CPX_PARAM_SCRIND"), false);
    @variable(lp, 0<= y[1:(J+K)*2]<=1)
    @variable(lp, 0<= uij[1:sum(Mij)]<=1);
    @variable(lp, 0<= ujk[1:sum(Mjk)]<=1);
    @variable(lp, 0<= ukl[1:sum(Mkl)]<=1);
    @variable( lp, 0<= xij[1:sum(Mij)] );
    @variable( lp, 0<= xjk[1:sum(Mjk)] );
    @variable( lp, 0<= xkl[1:sum(Mkl)] );
    @variable( lp, 0<= h[1:(J+K)*2] );
    @objective(lp, Min,  sum(c2.*y) +sum(repeat(vcs[1,:], outer=sum(Mij[1,:])).*xij[1:sum(Mij[1,:])]) +
        sum(sum(repeat(vcs[i,:], outer=sum(Mij[i,:])).*xij[sum(Mij[1:i-1,:])+1:sum(Mij[1:i,:])]) for i=2:I) +
        dot(e2,h) +sum(gij2[i]*uij[i] for i in findnz(gij2)[1]) + sum(gjk2[i]*ujk[i] for i in findnz(gjk2)[1]) + sum(gkl2[i].*ukl[i] for i in findnz(gkl2)[1]) +
        sum(tcp2.*xij)+sum(tcd2.*xjk)+sum(tcc2.*xkl)    +
            weight*(sum(repeat(ves[1,:], outer=sum(Mij[1,:])).*xij[1:sum(Mij[1,:])]) +
        sum(sum(repeat(ves[i,:], outer=sum(Mij[i,:])).*xij[sum(Mij[1:i-1,:])+1:sum(Mij[1:i,:])]) for i=2:I) +
        dot(q2,h) + sum(rij2.*xij)+sum(rjk2.*xjk)+sum(rkl2.*xkl))
    );
    @constraint(lp, sum(xij[m] for m=1:Mij[1,1])+ sum(xij[m+(sum(Mij[1:i-1,:]))] for i=2:I for m=1:Mij[i,1])  == sum(xjk[m] for m=1:sum(Mjk[1,:])) )
    @constraint(lp, [j=2:J], sum(xij[sum(Mij[1,1:j-1])+m] for m=1:Mij[1,j])+sum(xij[sum(Mij[i,1:j-1])+m+(sum(Mij[1:i-1,:]))] for i=2:I for m=1:Mij[i,j]) == sum(xjk[sum(Mjk[1:j-1,:]) + m] for m=1:sum(Mjk[j,:])) )
    @constraint(lp, sum(xjk[m] for m=1:Mjk[1,1])+sum(xjk[m+(sum(Mjk[1:j-1,:]))] for j=2:J for m=1:Mjk[j,1]) == sum(xkl[m] for m=1:sum(Mkl[1,:])) )
    @constraint(lp, [k=2:K],sum(xjk[sum(Mjk[1,1:k-1])+m] for m=1:Mjk[1,k])+sum(xjk[sum(Mjk[j,1:k-1])+m+(sum(Mjk[1:j-1,:]))] for j=2:J for m=1:Mjk[j,k]) == sum(xkl[sum(Mkl[1:k-1,:]) + m] for m=1:sum(Mkl[k,:])) )
    ########### constraint 4-6 #############
    @constraint(lp, sum(h[t] for t=1:2) ==sum(xij[m] for m=1:Mij[1,1])+sum(xij[m+(sum(Mij[1:i-1,:]))] for i=2:I for m=1:Mij[i,1]));
    @constraint(lp, [j=2:J], sum(h[2*(j-1)+t] for t=1:2) == sum(xij[sum(Mij[1,1:j-1])+m] for m=1:Mij[1,j])+sum(xij[sum(Mij[i,1:j-1])+m+(sum(Mij[1:i-1,:]))] for i=2:I for m=1:Mij[i,j]) );
    @constraint(lp,sum(h[2*J+t] for t=1:2) == sum(xjk[m] for m=1:Mjk[1,1])+sum(xjk[m+(sum(Mjk[1:j-1,:]))] for j=2:J for m=1:Mjk[j,1]) );
    @constraint(lp, [k=2:K], sum(h[2*J+2*(k-1)+t] for t=1:2) == sum(xjk[sum(Mjk[1,1:k-1])+m] for m=1:Mjk[1,k])+sum(xjk[sum(Mjk[j,1:k-1])+m+(sum(Mjk[1:j-1,:]))] for j=2:J for m=1:Mjk[j,k]));
    @constraint(lp, sum(xkl[m] for m=1:Mkl[1,1]) +sum(xkl[m+(sum(Mkl[1:k-1,:]))] for k=2:K for m=1:Mkl[k,1]) >= demand[1]);
    @constraint(lp, [l=2:L], sum(xkl[sum(Mkl[1,1:l-1]) + m] for m=1:Mkl[1,l])+ sum(xkl[sum(Mkl[1:k-1,:])+sum(Mkl[k,1:l-1])+m] for k=2:K for m=1:Mkl[k,l]) >= demand[l]);
    ########### constraint 7 #############
    @constraint(lp, sum(xij[1:sum(Mij[1,:])]) <= cas1[1]);
    @constraint(lp, [i=2:I],  sum(xij[sum(Mij[1:i-1,:])+1:sum(Mij[1:i,:])]) <= cas1[i]);
    ########### constraint 8 #############
    @constraint(lp,[j=1:J+K, t=1:2], sum(h[2*(j-1)+t]) <= capd1[j]*y[2*(j-1)+t]);
    ########### constraint 9 #############
    @constraint(lp,[j=1:J+K], sum(y[2*(j-1)+1:2*(j-1)+2]) <= 1);
    ########### constraint 10 #############
    @constraints(lp,begin
            sum(uij[1:Mij[1,1]]) <= 1
            [j=2:J], sum(uij[sum(Mij[1,1:j-1])+1:sum(Mij[1,1:j-1])+Mij[1,j]]) <= 1
            [i=2:I,j=2:J],  sum(uij[sum(Mij[1:i-1,:])+sum(Mij[i,1:j-1])+1:sum(Mij[1:i-1,:])+sum(Mij[i,1:j-1])+Mij[i,j]])<= 1
            sum(ujk[1:Mjk[1,1]]) <= 1
            [k=2:K], sum(ujk[sum(Mjk[1,1:k-1])+1:sum(Mjk[1,1:k-1])+Mjk[1,k]]) <= 1
            [j=2:J,k=2:K],  sum(ujk[sum(Mjk[1:j-1,:])+sum(Mjk[j,1:j-1])+1:sum(Mjk[1:j-1,:])+sum(Mjk[j,1:j-1])+Mjk[j,k]]) <= 1
            sum(ukl[1:Mkl[1,1]]) <= 1
            [l=2:L], sum(ukl[sum(Mkl[1,1:l-1])+1:sum(Mkl[1,1:l-1])+Mkl[1,l]]) <= 1
            [k=2:K,l=2:L],  sum(ukl[sum(Mkl[1:k-1,:])+sum(Mkl[k,1:l-1])+1:sum(Mkl[1:k-1,:])+sum(Mkl[k,1:l-1])+Mkl[k,l]])<= 1
    end);
    ########### constraint 11 #############
    @constraints(lp, begin
        [i=1:sum(Mij)], sum(xij[(i-1)+1:i]) <= bigM*uij[i]
        [i=1:sum(Mjk)], sum(xjk[(i-1)+1:i]) <= bigM*ujk[i]
        [i=1:sum(Mkl)], sum(xkl[(i-1)+1:i]) <= bigM*ukl[i]
    end);
    ########### constraint 12 #############
    @constraints(lp, begin
            [i in findnz(Vij2)[1]], sum(xij[(i-1)+1:i]) >= Vij2[i]*uij[i]
            [i in findnz(Vjk2)[1]], sum(xjk[(i-1)+1:i]) >= Vjk2[i]*ujk[i]
            # [i in findnz(Vkl)[1]], sum(xkl[(i-1)+1:i]) >= Vkl[i]*ukl[i]
    end);
    ########### constraint 13-14 #############
    @constraint(lp, sum(y[1:J*2]) <= Jmax);
    @constraint(lp, sum(y[J*2+1:end]) <= Kmax);
    return lp
end
function flip(x_h,j,e)
    if x_h[e[j]]==1
        x_h[e[j]] = 0
    else
        x_h[e[j]] = 1
    end
    return x_h
end
function flipoper(Tabu,x_t,x_r)
    e = sortperm(abs.(x_t-x_r),rev=true)
    xi = []
    x_h = copy(x_r)
    j = 1
    M=length(x_t) #
    while j<=M && xi==[]
        x_h = flip(x_h,j,e)
        if x_h ∉ Tabu
            xi=x_h
        else
            j+=1
        end
    end
    if xi==[]
        while j<=M
            x_h=copy(x_r)
            Num = Int64(rand(ceil(length(x_r)/2):length(x_r)-1))
            R = sample(1:M,Num, replace=false)
            for i in R
                x_h = flip(x_h,r)
                if x_h ∉ Tabu
                    xi = x_h
                end
            end
            j+=1
        end
    end
    return xi
end
function FP_FBcheck(model,yr)#,u1r,u2r,u3r)
    JuMP.fix.(model[:y],yr; force=true)
    optimize!(model)
    if termination_status(model) == MOI.OPTIMAL
        return true
    else
        return false
    end
end

# function feasicheck(yr,u1r,u2r,u3r)
# end
function fbsearch(yr,u1r,u2r,u3r) #solveLP
    idy_0 = findall(k->k==0, yr)
    idy_1 = findall(k->k==1, yr)
    idu1_0 = findall(k->k==0, u1r)
    idu1_1 = findall(k->k==1, u1r)
    idu2_0 = findall(k->k==0, u2r)
    idu2_1 = findall(k->k==1, u2r)
    idu3_0 = findall(k->k==0, u3r)
    idu3_1 = findall(k->k==1, u3r)
    @objective( distmodel, Min, sum(distmodel[:y][i] for i in idy_0) + sum(1-(distmodel[:y][j]) for j in idy_1) +
        sum(distmodel[:uij][i] for i in idu1_0) + sum(1-(distmodel[:uij][j]) for j in idu1_1)+
        sum(distmodel[:ujk][i] for i in idu2_0) + sum(1-(distmodel[:ujk][j]) for j in idu2_1)+
        sum(distmodel[:ukl][i] for i in idu3_0) + sum(1-(distmodel[:ukl][j]) for j in idu3_1))
    optimize!(distmodel)
    if termination_status(distmodel) == MOI.OPTIMAL
        return JuMP.value.(distmodel[:y]),JuMP.value.(distmodel[:uij]),JuMP.value.(distmodel[:ujk]),JuMP.value.(distmodel[:ukl])
    else
        return 0,0,0,0
    end
end

len = [length(lpm[:y]),length(lpm[:uij]),length(lpm[:ujk]),length(lpm[:ukl]),length(lpm[:xij]),length(lpm[:xjk]),length(lpm[:xkl]),length(lpm[:h])]

function dominated(x,P)
    st = false
    for k=1:length(P)
        if all( x .>= P[k])
            st=true; break
        else
            continue
        end
    end
    return st
end
function getobjval(x)
    y = x[1:len[1]];
    uij = x[1+len[1]:len[1]+len[2]];
    ujk = x[1+sum(len[i] for i=1:2):sum(len[i] for i=1:3)]
    ukl = x[1+sum(len[i] for i=1:3):sum(len[i] for i=1:4)]
    xij = x[1+sum(len[i] for i=1:4):sum(len[i] for i=1:5)]
    xjk = x[1+sum(len[i] for i=1:5):sum(len[i] for i=1:6)]
    xkl = x[1+sum(len[i] for i=1:6):sum(len[i] for i=1:7)]
    h = x[1+sum(len[i] for i=1:7):sum(len[i] for i=1:8)]

    obj1 =  sum(c2.*y) +sum(repeat(vcs[1,:], outer=sum(Mij[1,:])).*xij[1:sum(Mij[1,:])]) +
        sum(sum(repeat(vcs[i,:], outer=sum(Mij[i,:])).*xij[sum(Mij[1:i-1,:])+1:sum(Mij[1:i,:])]) for i=2:I) +
        dot(e2,h) +sum(gij2[i]*uij[i] for i in findnz(gij2)[1]) + sum(gjk2[i]*ujk[i] for i in findnz(gjk2)[1]) + sum(gkl2[i].*ukl[i] for i in findnz(gkl2)[1]) +
        sum(tcp2.*xij)+sum(tcd2.*xjk)+sum(tcc2.*xkl)
    obj2 = sum(repeat(ves[1,:], outer=sum(Mij[1,:])).*xij[1:sum(Mij[1,:])]) +
        sum(sum(repeat(ves[i,:], outer=sum(Mij[i,:])).*xij[sum(Mij[1:i-1,:])+1:sum(Mij[1:i,:])]) for i=2:I) +
        dot(q2,h) + sum(rij2.*xij)+sum(rjk2.*xjk)+sum(rkl2.*xkl)
    return [obj1,obj2]
end

function FPplus(dvar,len,TL,Y_N) #𝚯
    X = copy(dvar); IGPair=[]; Tabu = []; newsol=0; t0=time(); PF = copy(Y_N);
    Y = []; U1 = []; U2= []; U3 = [];
    while time()-t0 < TL && length(IGPair)<(length(PF)*(length(PF)-1))
        I,G = StatsBase.sample(1:length(X), 2, replace=false)
        x1 = X[I][1:sum(len[i] for i=1:4)]; x2 = X[G][1:sum(len[i] for i=1:4)];
        λ = round(rand(Float64, 1)[1]; digits=1)
        x_t = x1*λ + x2*(1-λ);
        # x_t = x1*.5 + x2*.5;
        yt = x_t[1:len[1]]; u1t = x_t[1+len[1]:len[1]+len[2]];
        u2t = x_t[1+sum(len[i] for i=1:2):sum(len[i] for i=1:3)]
        u3t = x_t[1+sum(len[i] for i=1:3):sum(len[i] for i=1:4)]
        SearchDone = false; iter=0;
        Max_iter = 10 # length(findall(i-> 0<i<1,x_t))
        while [I,G]∉IGPair && iter<Max_iter && SearchDone == false
            # x_r = round.(Int,x_t);
            yr = round.(Int, yt); u1r = round.(Int, u1t);
            u2r = round.(Int, u2t); u3r = round.(Int, u3t);

            if FP_FBcheck(fbmodel,yr) == true
                sol = value.(all_variables(fbmodel)); ndp = getobjval(sol)
                if sol ∉ X  && dominated(ndp,PF)==false
                    push!(X,sol); push!(PF,ndp)
                    push!(Y,yr); push!(U1,u1r); push!(U2,u2r); push!(U3,u3r);
                    newsol+=1; SearchDone = true
                    println("rounding worked")
                end
            else
                if [yr;u1r;u2r;u3r] ∈ Tabu
                    yr = flipoper(Y,yt,yr); u1r = flipoper(U1,u1t,u1r); u2r = flipoper(U2,u2t,u2r); u3r = flipoper(U3,u3t,u3r)
                    if any(i->i==[], [yr,u1r,u2r,u3r])
                        SearchDone = true;
                        # println("flip failed")
                    else
                        if FP_FBcheck(fbmodel,yr) == true
                            sol = value.(all_variables(fbmodel)); ndp = getobjval(sol)
                            if sol ∉ X && dominated(ndp,PF)==false
                                push!(X,sol); push!(PF,ndp)
                                push!(Y,yr); push!(U1,u1r); push!(U2,u2r); push!(U3,u3r);
                                newsol+=1; SearchDone = true;
                                println("flip worked")
                            end
                        end
                    end
                end
                if time()-t0 >= TL
                    break
                end
                if SearchDone == false
                    push!(Tabu,[yr;u1r;u2r;u3r])
                    yt,u1t,u2t,u3t = fbsearch(yr,u1r,u2r,u3r)
                    if any(i->i==0, [yt,u1t,u2t,u3t])  #when there's no new feasible lp sol
                        # println("no solution")
                        SearchDone = true
                    end
                end
            end
			iter+=1
        end
        push!(IGPair,[I,G])

    end
    return X,PF,newsol,IGPair
end
function FP(candX,len,TL)
    X = []; PF =[]; Tabu = []; t0 = time(); newsol = 0; candlist = copy(candX)
    Y = []; U1 = []; U2= []; U3 = [];
    while candlist != [] && time() - t0 < TL
        k = rand(1:length(candlist))
        x_t = candlist[k]
        yt = x_t[1:len[1]];
        u1t = x_t[1+len[1]:len[1]+len[2]];
        u2t = x_t[1+sum(len[i] for i=1:2):sum(len[i] for i=1:3)]
        u3t = x_t[1+sum(len[i] for i=1:3):sum(len[i] for i=1:4)]
        SearchDone = false; iter=0;
        Max_iter =  5#length(findall(i-> 0<i<1,yt))
        while iter < Max_iter && SearchDone == false && time() - t0 < TL
            # x_r = round.(Int,x_t);
            yr = round.(Int, yt); u1r = round.(Int, u1t);
            u2r = round.(Int, u2t); u3r = round.(Int, u3t);

            if FP_FBcheck(fbmodel,yr) == true #,u1r,u2r,u3r)
                sol = value.(all_variables(fbmodel)); ndp = getobjval(sol)
                if sol ∉ X  && dominated(ndp,PF)==false
                    push!(X,sol); push!(PF,ndp)
                    push!(Y,yr); push!(U1,u1r); push!(U2,u2r); push!(U3,u3r);
                    newsol+=1; SearchDone = true
                    deleteat!(candlist, k);
                    println("rounding worked")
                else
                    # println("rounding failed")
                end
            else
                if [yr;u1r;u2r;u3r] ∈ Tabu
                    yr = flipoper(Y,yt,yr); u1r = flipoper(U1,u1t,u1r); u2r = flipoper(U2,u2t,u2r); u3r = flipoper(U3,u3t,u3r)
                    if any(i->i==[], [yr,u1r,u2r,u3r])
                        SearchDone = true;
                        # println("flip failed")
                    else
                        if FP_FBcheck(fbmodel,yr) == true #,u1r,u2r,u3r)
                            sol = value.(all_variables(fbmodel)); ndp = getobjval(sol)
                            if sol ∉ X && dominated(ndp,PF)==false
                                push!(X,sol); push!(PF,ndp)
                                push!(Y,yr); push!(U1,u1r); push!(U2,u2r); push!(U3,u3r);
                                newsol+=1; SearchDone = true;
                                deleteat!(candlist, k);
                                println("flip worked")
                            end
                        end
                    end
                end
                if time()-t0 >= TL
                    break
                end
                if SearchDone == false
                    push!(Tabu,[yr;u1r;u2r;u3r])
                    yt,u1t,u2t,u3t = fbsearch(yr,u1r,u2r,u3r)
                    if any(i->i==0, [yt,u1t,u2t,u3t])  #when there's no new feasible lp sol
                        deleteat!(candlist, k);
                        # println("no solution")
                        SearchDone = true
                    end
                end
            end
			iter+=1
        end
    end
    return X,PF,newsol,candlist
end

# optimize!(fbmodel); optimize!(dist);
fbmodel = FP_Model(mweight)
distmodel = LP_Model(mweight)
FPtime = @CPUelapsed fx,fy,fn,usedPairs = FPplus(mxe,len,20,myn)

fbmodel = FP_Model(lweight)
distmodel = LP_Model(lweight)
Fptime = @CPUelapsed lx,ly,ln,candlist = FP(lxe,len,20)

function Postpro(P,Pobj)
    copysol = Dict(); copyobj = Dict();
    for i=1:length(Pobj)
        copysol[i] = P[i]
        copyobj[i] = Pobj[i]
    end
    for i=1:length(Pobj)-1
        for j=i+1:length(Pobj)
            if all(Pobj[i] .>= Pobj[j]) == true #dominated by PF[j]
                copyobj[i]=nothing; copysol[i]=nothing; break
            elseif all(Pobj[j] .>= Pobj[i]) == true
                copyobj[j]=nothing; copysol[j]=nothing;
            end
        end
    end

    finalsol = filter!(a->a!=nothing, collect(values(copysol)))
    finalobj = filter!(a->a!=nothing, collect(values(copyobj)))

    return finalsol,finalobj
end
################################ Path Relinking  ###############################
function PR_Model(weight)
    model = Model(CPLEX.Optimizer); set_silent(model)
    MOI.set(model, MOI.NumberOfThreads(), 1)
    # MOI.set(model, MOI.RawParameter("CPX_PARAM_SCRIND"), false);
    @variable(model, y[1:(J+K)*2], Bin)
    @variable(model, uij[1:sum(Mij)], Bin);
    @variable(model, ujk[1:sum(Mjk)], Bin);
    @variable(model, ukl[1:sum(Mkl)], Bin);
    @variable( model, 0<= xij[1:sum(Mij)] );
    @variable( model, 0<= xjk[1:sum(Mjk)] );
    @variable( model, 0<= xkl[1:sum(Mkl)] );
    @variable( model, 0<= h[1:(J+K)*2] );
    @objective(model, Min,  sum(c2.*y) +sum(repeat(vcs[1,:], outer=sum(Mij[1,:])).*xij[1:sum(Mij[1,:])]) +
        sum(sum(repeat(vcs[i,:], outer=sum(Mij[i,:])).*xij[sum(Mij[1:i-1,:])+1:sum(Mij[1:i,:])]) for i=2:I) +
        dot(e2,h) +sum(gij2[i]*uij[i] for i in findnz(gij2)[1]) + sum(gjk2[i]*ujk[i] for i in findnz(gjk2)[1]) + sum(gkl2[i].*ukl[i] for i in findnz(gkl2)[1]) +
        sum(tcp2.*xij)+sum(tcd2.*xjk)+sum(tcc2.*xkl)    +
        weight*(sum(repeat(ves[1,:], outer=sum(Mij[1,:])).*xij[1:sum(Mij[1,:])]) +
        sum(sum(repeat(ves[i,:], outer=sum(Mij[i,:])).*xij[sum(Mij[1:i-1,:])+1:sum(Mij[1:i,:])]) for i=2:I) +
        dot(q2,h) + sum(rij2.*xij)+sum(rjk2.*xjk)+sum(rkl2.*xkl))
    );
    ########## constraint 3 #############
    @constraint(model, sum(xij[m] for m=1:Mij[1,1])+ sum(xij[m+(sum(Mij[1:i-1,:]))] for i=2:I for m=1:Mij[i,1])  == sum(xjk[m] for m=1:sum(Mjk[1,:])) )
    @constraint(model, [j=2:J], sum(xij[sum(Mij[1,1:j-1])+m] for m=1:Mij[1,j])+sum(xij[sum(Mij[i,1:j-1])+m+(sum(Mij[1:i-1,:]))] for i=2:I for m=1:Mij[i,j]) == sum(xjk[sum(Mjk[1:j-1,:]) + m] for m=1:sum(Mjk[j,:])) )
    @constraint(model, sum(xjk[m] for m=1:Mjk[1,1])+sum(xjk[m+(sum(Mjk[1:j-1,:]))] for j=2:J for m=1:Mjk[j,1]) == sum(xkl[m] for m=1:sum(Mkl[1,:])) )
    @constraint(model, [k=2:K],sum(xjk[sum(Mjk[1,1:k-1])+m] for m=1:Mjk[1,k])+sum(xjk[sum(Mjk[j,1:k-1])+m+(sum(Mjk[1:j-1,:]))] for j=2:J for m=1:Mjk[j,k]) == sum(xkl[sum(Mkl[1:k-1,:]) + m] for m=1:sum(Mkl[k,:])) )
    ########### constraint 4-6 #############
    @constraint(model, sum(h[t] for t=1:2) ==sum(xij[m] for m=1:Mij[1,1])+sum(xij[m+(sum(Mij[1:i-1,:]))] for i=2:I for m=1:Mij[i,1]));
    @constraint(model, [j=2:J], sum(h[2*(j-1)+t] for t=1:2) == sum(xij[sum(Mij[1,1:j-1])+m] for m=1:Mij[1,j])+sum(xij[sum(Mij[i,1:j-1])+m+(sum(Mij[1:i-1,:]))] for i=2:I for m=1:Mij[i,j]) );
    @constraint(model,sum(h[2*J+t] for t=1:2) == sum(xjk[m] for m=1:Mjk[1,1])+sum(xjk[m+(sum(Mjk[1:j-1,:]))] for j=2:J for m=1:Mjk[j,1]) );
    @constraint(model, [k=2:K], sum(h[2*J+2*(k-1)+t] for t=1:2) == sum(xjk[sum(Mjk[1,1:k-1])+m] for m=1:Mjk[1,k])+sum(xjk[sum(Mjk[j,1:k-1])+m+(sum(Mjk[1:j-1,:]))] for j=2:J for m=1:Mjk[j,k]));
    @constraint(model, sum(xkl[m] for m=1:Mkl[1,1]) +sum(xkl[m+(sum(Mkl[1:k-1,:]))] for k=2:K for m=1:Mkl[k,1]) >= demand[1]);
    @constraint(model, [l=2:L], sum(xkl[sum(Mkl[1,1:l-1]) + m] for m=1:Mkl[1,l])+ sum(xkl[sum(Mkl[1:k-1,:])+sum(Mkl[k,1:l-1])+m] for k=2:K for m=1:Mkl[k,l]) >= demand[l]);
    ########### constraint 7 #############
    @constraint(model, sum(xij[1:sum(Mij[1,:])]) <= cas1[1]);
    @constraint(model, [i=2:I],  sum(xij[sum(Mij[1:i-1,:])+1:sum(Mij[1:i,:])]) <= cas1[i]);
    ########### constraint 8 #############
    @constraint(model,[j=1:J+K, t=1:2], sum(h[2*(j-1)+t]) <= capd1[j]*y[2*(j-1)+t]);
    ########### constraint 9 #############
    # @constraint(model,[j=1:J+K], sum(y[2*(j-1)+1:2*(j-1)+2]) <= 1);
    ########### constraint 10 #############
    @constraints(model,begin
            sum(uij[1:Mij[1,1]]) <= 1
            [j=2:J], sum(uij[sum(Mij[1,1:j-1])+1:sum(Mij[1,1:j-1])+Mij[1,j]]) <= 1
            [i=2:I,j=2:J],  sum(uij[sum(Mij[1:i-1,:])+sum(Mij[i,1:j-1])+1:sum(Mij[1:i-1,:])+sum(Mij[i,1:j-1])+Mij[i,j]])<= 1
            sum(ujk[1:Mjk[1,1]]) <= 1
            [k=2:K], sum(ujk[sum(Mjk[1,1:k-1])+1:sum(Mjk[1,1:k-1])+Mjk[1,k]]) <= 1
            [j=2:J,k=2:K],  sum(ujk[sum(Mjk[1:j-1,:])+sum(Mjk[j,1:j-1])+1:sum(Mjk[1:j-1,:])+sum(Mjk[j,1:j-1])+Mjk[j,k]]) <= 1
            sum(ukl[1:Mkl[1,1]]) <= 1
            [l=2:L], sum(ukl[sum(Mkl[1,1:l-1])+1:sum(Mkl[1,1:l-1])+Mkl[1,l]]) <= 1
            [k=2:K,l=2:L],  sum(ukl[sum(Mkl[1:k-1,:])+sum(Mkl[k,1:l-1])+1:sum(Mkl[1:k-1,:])+sum(Mkl[k,1:l-1])+Mkl[k,l]])<= 1
    end);
    ########### constraint 11 #############
    @constraints(model, begin
        [i=1:sum(Mij)], sum(xij[(i-1)+1:i]) <= bigM*uij[i]
        [i=1:sum(Mjk)], sum(xjk[(i-1)+1:i]) <= bigM*ujk[i]
        [i=1:sum(Mkl)], sum(xkl[(i-1)+1:i]) <= bigM*ukl[i]
    end);
    ########### constraint 12 #############
    @constraints(model, begin
            [i in findnz(Vij2)[1]], sum(xij[(i-1)+1:i]) >= Vij2[i]*uij[i]
            [i in findnz(Vjk2)[1]], sum(xjk[(i-1)+1:i]) >= Vjk2[i]*ujk[i]
            # [i in findnz(Vkl)[1]], sum(xkl[(i-1)+1:i]) >= Vkl[i]*ukl[i]
    end);
    ########### constraint 13-14 #############
    # @constraint(model, sum(y[1:J*2]) <= Jmax);
    # @constraint(model, sum(y[J*2+1:end]) <= Kmax);
    return model
end
function createNB(SI,dif,exploredSI)
    neibour = []; #neiobj = [];
    for i in dif
        cpSI = copy(SI)
        # if cpSI[i] == round(cpSI[i]) #if vari is int value
        if cpSI[i] == 1
            cpSI[i] = 0
        else
            cpSI[i] = 1
        end
        push!(neibour, cpSI);# push!(neiobj, binobj(cpSI,bvar))
        # else #if vari is fractional
        #     cpSI[i] = 1; push!(neibour, cpSI); push!( neiobj, binobj(cpSI,bvar) )
        #     cpSI = copy(SI)
        #     cpSI[i] = 0; push!(neibour, cpSI); push!( neiobj, binobj(cpSI,bvar) )
        # end
    end
    idx = findall(i-> neibour[i] in exploredSI, 1:length(neibour))
    neibour1 = setdiff(neibour,exploredSI)
    neibour2 = StatsBase.sample(neibour1, round(Int,length(neibour1)*0.5) , replace=false)
    # neiobj = neiobj[setdiff(1:end, idx),:]
    return neibour2#,neiobj
end
function nextSI(neibour,SI)
    SIobj = getobjval(SI)
    neiobj = [getobjval(neibour[i]) for i=1:length(neibour)]
    for i=1:length(neiobj)
        if length(neibour) == 1  #if there is one candiate sol
            return neibour[1]#,neibour[1]
        else length(neiobj) > 1 # if there are multiple candiates, check the improved ratio
            ratiotb = zeros(length(neiobj),length(neiobj[1]))
            for i=1:length(neiobj)
                ratiotb[i,:] = neiobj[i]./SIobj
            end
            ranktb = zeros(length(neiobj),length(neiobj[1]))
            for i=1:length(neiobj[1])
                ranktb[:,i] = tiedrank(ratiotb[:,i])
            end
            ranksum = [sum(ranktb[i,:]) for i=1:length(neiobj)]
            mostimp = findall(x-> x == maximum(ranksum), ranksum)
            k = rand(mostimp)[1]
            return neibour[k]#,neiobj[k]
        end
    end
end
function First_FBcheck(neibour)
    neibour1 = filter(p->(sum(p[1:J*2]) <= Jmax && sum(p[1:J*2]) >0) || (sum(p[J*2+1:end])<= Kmax) && (sum(p[J*2+1:end]) >0), neibour)
    for j=1:J+K
        filter!(p-> sum(p[2*(j-1)+1:2*(j-1)+2])  <= 1, neibour1)
    end
    return neibour1
end
function Second_FBcheck(model,yr)#,u1r,u2r,u3r)
    JuMP.fix.(model[:y],yr; force=true)
    # JuMP.fix.(model[:uij],u1r; force=true)
    # JuMP.fix.(model[:ujk],u2r; force=true)
    # JuMP.fix.(model[:ukl],u3r; force=true)
    optimize!(model)
    if termination_status(model) == MOI.OPTIMAL
        return true
    else
        return false
    end
end
function PR(X,Y,len,TL)
    candX = copy(X); candY = copy(Y); bvar = sum(len[i] for i=1:4);
    IGPair=[]; exploredSI = []; t0=time(); iter=0;
    while time()-t0 < TL && length(IGPair)<(length(candY)*(length(candY)-1))
        @label NewIter
	    I,G = StatsBase.sample(1:length(candX), 2, replace=false);
        SI = candX[I]; SG = candX[G];
        # SI_r = round.(SI); SG_r = round.(SG)
        dif = findall(i-> SI[1:len[1]][i]!=SG[1:len[1]][i], 1:len[1])
         #Max_iter = 20;
        while length(dif)>0 && [I,G]∉IGPair && (time()-t0<TL) #&& iter<Max_iter
            neibour = createNB(SI[1:len[1]],dif,exploredSI)
            # println("# of neighbours: ", length(neibour2))
            if (length(neibour)==0) #(time()-t0 >= TL)
                @goto NewIter
            else
                candSI =[]
                # l=1;
                # newsol=0;
                # while (time()-t0<TL) && l<=length(neibour) && newsol <=1 #floor(dt1.N["supplier"]/3) &&
                neibour2 = First_FBcheck(neibour)
                for l=1:length(neibour2)
                    st = Second_FBcheck( prmodel, neibour2[l][1:len[1]])
                    # , neibour2[l][1+len[1]:sum(len[i] for i=1:2)],
                    #     neibour2[l][1+sum(len[i] for i=1:2):sum(len[i] for i=1:3)],neibour2[l][1+sum(len[i] for i=1:3):sum(len[i] for i=1:4)] )
                    # @show st
                    if st==true
                        sol = value.(all_variables(prmodel)); ndp = getobjval(sol)
                        push!(candSI,sol)
                        if sol ∉ candX && dominated(ndp,candY)==false
                            push!(candX, sol); push!(candY, ndp);
                            # newsol+=1;
                            # println("new sol");
                        end
                    end
                #     l+=1;
                end
            end
            if candSI == []
                push!(IGPair,[I,G]);
                # @show iter+=1
                @goto NewIter
            else
                SI = nextSI(candSI,SI)
                if SI∉candX
                    push!(exploredSI,SI);
                end
                # @show iter+=1
            end
        end
        push!(IGPair,[I,G]);
        iter+=1
    end
    return candX,candY,IGPair
end
prmodel = PR_Model(mweight)
fpy
PRtime = @CPUelapsed px,py,pairs = PR(fx,fy,len,15)
prx,pry = Postpro(px,py)
setdiff(py,fy)

prmodel = PR_Model(lweight)
Prtime = @CPUelapsed plx,ply,used = PR(lx,ly,len,15)
setdiff(ply,ly)
p0x,p0y = Postpro(plx,ply)
################### Fix the binary values and Solve Dicho  ####################
function smscnd_LP()
    smscnd = vModel(optimizer_with_attributes(
            CPLEX.Optimizer,
            "CPX_PARAM_EPGAP" => 1e-8
          ));
    set_silent(smscnd)
    MOI.set(smscnd, MOI.NumberOfThreads(), 1)
    @variable(smscnd, 0<= y[1:(J+K)*2]<=1)
    @variable(smscnd, 0<= uij[1:sum(Mij)]<=1);
    @variable(smscnd, 0<= ujk[1:sum(Mjk)]<=1);
    @variable(smscnd, 0<= ukl[1:sum(Mkl)]<=1);
    @variable( smscnd, 0<= xij[1:sum(Mij)] );
    @variable( smscnd, 0<= xjk[1:sum(Mjk)] );
    @variable( smscnd, 0<= xkl[1:sum(Mkl)] );
    @variable( smscnd, 0<= h[1:(J+K)*2] );


    @addobjective(smscnd, Min, sum(c2.*y) +sum(repeat(vcs[1,:], outer=sum(Mij[1,:])).*xij[1:sum(Mij[1,:])]) +
        sum(sum(repeat(vcs[i,:], outer=sum(Mij[i,:])).*xij[sum(Mij[1:i-1,:])+1:sum(Mij[1:i,:])]) for i=2:I) +
        dot(e2,h) +sum(gij2[i]*uij[i] for i in findnz(gij2)[1]) + sum(gjk2[i]*ujk[i] for i in findnz(gjk2)[1]) + sum(gkl2[i].*ukl[i] for i in findnz(gkl2)[1]) +
        sum(tcp2.*xij)+sum(tcd2.*xjk)+sum(tcc2.*xkl));
    @addobjective(smscnd, Min, sum(repeat(ves[1,:], outer=sum(Mij[1,:])).*xij[1:sum(Mij[1,:])]) +
        sum(sum(repeat(ves[i,:], outer=sum(Mij[i,:])).*xij[sum(Mij[1:i-1,:])+1:sum(Mij[1:i,:])]) for i=2:I) +
        dot(q2,h) + sum(rij2.*xij)+sum(rjk2.*xjk)+sum(rkl2.*xkl));
    ########## constraint 3 #############
    @constraint(smscnd, sum(xij[m] for m=1:Mij[1,1])+ sum(xij[m+(sum(Mij[1:i-1,:]))] for i=2:I for m=1:Mij[i,1])  == sum(xjk[m] for m=1:sum(Mjk[1,:])) )
    @constraint(smscnd, [j=2:J], sum(xij[sum(Mij[1,1:j-1])+m] for m=1:Mij[1,j])+sum(xij[sum(Mij[i,1:j-1])+m+(sum(Mij[1:i-1,:]))] for i=2:I for m=1:Mij[i,j]) == sum(xjk[sum(Mjk[1:j-1,:]) + m] for m=1:sum(Mjk[j,:])) )
    @constraint(smscnd, sum(xjk[m] for m=1:Mjk[1,1])+sum(xjk[m+(sum(Mjk[1:j-1,:]))] for j=2:J for m=1:Mjk[j,1]) == sum(xkl[m] for m=1:sum(Mkl[1,:])) )
    @constraint(smscnd, [k=2:K],sum(xjk[sum(Mjk[1,1:k-1])+m] for m=1:Mjk[1,k])+sum(xjk[sum(Mjk[j,1:k-1])+m+(sum(Mjk[1:j-1,:]))] for j=2:J for m=1:Mjk[j,k]) == sum(xkl[sum(Mkl[1:k-1,:]) + m] for m=1:sum(Mkl[k,:])) )
    ########### constraint 4-6 #############
    @constraint(smscnd, sum(h[t] for t=1:2) ==sum(xij[m] for m=1:Mij[1,1])+sum(xij[m+(sum(Mij[1:i-1,:]))] for i=2:I for m=1:Mij[i,1]));
    @constraint(smscnd, [j=2:J], sum(h[2*(j-1)+t] for t=1:2) == sum(xij[sum(Mij[1,1:j-1])+m] for m=1:Mij[1,j])+sum(xij[sum(Mij[i,1:j-1])+m+(sum(Mij[1:i-1,:]))] for i=2:I for m=1:Mij[i,j]) );
    @constraint(smscnd,sum(h[2*J+t] for t=1:2) == sum(xjk[m] for m=1:Mjk[1,1])+sum(xjk[m+(sum(Mjk[1:j-1,:]))] for j=2:J for m=1:Mjk[j,1]) );
    @constraint(smscnd, [k=2:K], sum(h[2*J+2*(k-1)+t] for t=1:2) == sum(xjk[sum(Mjk[1,1:k-1])+m] for m=1:Mjk[1,k])+sum(xjk[sum(Mjk[j,1:k-1])+m+(sum(Mjk[1:j-1,:]))] for j=2:J for m=1:Mjk[j,k]));
    @constraint(smscnd, sum(xkl[m] for m=1:Mkl[1,1]) +sum(xkl[m+(sum(Mkl[1:k-1,:]))] for k=2:K for m=1:Mkl[k,1]) >= demand[1]);
    @constraint(smscnd, [l=2:L], sum(xkl[sum(Mkl[1,1:l-1]) + m] for m=1:Mkl[1,l])+ sum(xkl[sum(Mkl[1:k-1,:])+sum(Mkl[k,1:l-1])+m] for k=2:K for m=1:Mkl[k,l]) >= demand[l]);
    ########### constraint 7 #############
    @constraint(smscnd, sum(xij[1:sum(Mij[1,:])]) <= cas1[1]);
    @constraint(smscnd, [i=2:I],  sum(xij[sum(Mij[1:i-1,:])+1:sum(Mij[1:i,:])]) <= cas1[i]);
    ########### constraint 8 #############
    @constraint(smscnd,[j=1:J+K, t=1:2], sum(h[2*(j-1)+t]) <= capd1[j]*y[2*(j-1)+t]);
    ########### constraint 9 #############
    @constraint(smscnd,[j=1:J+K], sum(y[2*(j-1)+1:2*(j-1)+2]) <= 1);
    ########### constraint 10 #############
    @constraints(smscnd,begin
            sum(uij[1:Mij[1,1]]) <= 1
            [j=2:J], sum(uij[sum(Mij[1,1:j-1])+1:sum(Mij[1,1:j-1])+Mij[1,j]]) <= 1
            [i=2:I,j=2:J],  sum(uij[sum(Mij[1:i-1,:])+sum(Mij[i,1:j-1])+1:sum(Mij[1:i-1,:])+sum(Mij[i,1:j-1])+Mij[i,j]])<= 1
            sum(ujk[1:Mjk[1,1]]) <= 1
            [k=2:K], sum(ujk[sum(Mjk[1,1:k-1])+1:sum(Mjk[1,1:k-1])+Mjk[1,k]]) <= 1
            [j=2:J,k=2:K],  sum(ujk[sum(Mjk[1:j-1,:])+sum(Mjk[j,1:j-1])+1:sum(Mjk[1:j-1,:])+sum(Mjk[j,1:j-1])+Mjk[j,k]]) <= 1
            sum(ukl[1:Mkl[1,1]]) <= 1
            [l=2:L], sum(ukl[sum(Mkl[1,1:l-1])+1:sum(Mkl[1,1:l-1])+Mkl[1,l]]) <= 1
            [k=2:K,l=2:L],  sum(ukl[sum(Mkl[1:k-1,:])+sum(Mkl[k,1:l-1])+1:sum(Mkl[1:k-1,:])+sum(Mkl[k,1:l-1])+Mkl[k,l]])<= 1
    end);
    ########### constraint 11 #############
    @constraints(smscnd, begin
        [i=1:sum(Mij)], sum(xij[(i-1)+1:i]) <= bigM*uij[i]
        [i=1:sum(Mjk)], sum(xjk[(i-1)+1:i]) <= bigM*ujk[i]
        [i=1:sum(Mkl)], sum(xkl[(i-1)+1:i]) <= bigM*ukl[i]
    end);
    ########### constraint 12 #############
    @constraints(smscnd, begin
            [i in findnz(Vij2)[1]], sum(xij[(i-1)+1:i]) >= Vij2[i]*uij[i]
            [i in findnz(Vjk2)[1]], sum(xjk[(i-1)+1:i]) >= Vjk2[i]*ujk[i]
            # [i in findnz(Vkl)[1]], sum(xkl[(i-1)+1:i]) >= Vkl[i]*ukl[i]
    end);
    ########### constraint 13-14 #############
    @constraint(smscnd, sum(y[1:J*2]) <= Jmax);
    @constraint(smscnd, sum(y[J*2+1:end]) <= Kmax);
    return smscnd
end

function FixedBinDicho(pry)
    linesg = Dict()
    for k=1:length(pry)
        model = smscnd_LP()
        JuMP.fix.(model[:y], prx[k][1:len[1]]; force = true)
        JuMP.fix.(model[:uij], prx[k][1+len[1]:sum(len[i] for i=1:2)]; force = true)
        JuMP.fix.(model[:ujk], prx[k][1+sum(len[i] for i=1:2):sum(len[i] for i=1:3)]; force = true)
        JuMP.fix.(model[:ukl], prx[k][1+sum(len[i] for i=1:3):sum(len[i] for i=1:4)]; force = true)
        @CPUtime vSolve(model, 120, method=:dicho, verbose=false)
        res = getvOptData(model);
        if res!= []
            linesg[k]=res.Y_N
        end
    end
    return linesg
end

linesg = FixedBinDicho(pry)
######################  Visualisation  #########################
using PlotlyJS,DataFrames,DelimitedFiles,JLD2,Colors

mm1 = [myn[i][1] for i=1:length(myn)]; mm2 = [myn[i][2] for i=1:length(myn)]
fy1= [fpy[i][1] for i=1:length(fpy)]; fy2 = [fpy[i][2] for i=1:length(fpy)]
py1 = [pry[i][1] for i=1:length(pry)]; py2 = [pry[i][2] for i=1:length(pry)];
p1y1 = [p0y[i][1] for i=1:length(p0y)]; p1y2 = [p0y[i][2] for i=1:length(p0y)];


layout = Layout(
    title="Plot Title",
    xaxis_title="Cost",
    yaxis_title="CO2 emission",
    legend_title="Legend Title",
    font=attr(
        family="Courier New, monospace",
        size=18
    )
)
plot_array = GenericTrace[]
trace1 = scatter(x=mm1,y=mm2,name="smSCND PF", mode="markers+lines", marker=attr(color="Crimson"))
trace2 = scatter(x=fy1, y=fy2,  name="Dicho+FFP", mode="markers", marker=attr(color="Blue"))
trace3 = scatter(x=py1, y=py2,  name="MIP+FP+PR", mode="markers", marker=attr(color="LimeGreen"))
trace4 = scatter(x=p1y1, y=p1y2,  name="LP+FP+PR", mode="markers", marker=attr(color="DarkOrange"))

push!(plot_array,trace1,trace3,trace4);
for i=1:length(linesg)
    if linesg[i]!=[]
        tradeoffs = scatter(x=[linesg[i][j][1] for j=1:length(linesg[i])],y=[linesg[i][j][2] for j=1:length(linesg[i])], mode="markers+lines", color=rand(3))
        push!(plot_array,tradeoffs)
    end
end

fig = plot(plot_array, layout)
savefig(fig,"/home/ak121396/Pictures/smSCNDins.png")
