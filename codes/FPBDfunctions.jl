struct MasterP
    y::Vector{VariableRef}
    uij::Vector{VariableRef}
    ujk::Vector{VariableRef}
    ukl::Vector{VariableRef}
    θ::VariableRef
    m::Model
end
struct DualSP
    # data::SubProblemData
    α1::Vector{VariableRef}
    α2::Vector{VariableRef}
    α3::Vector{VariableRef}
    α4::Vector{VariableRef}
    α5::Vector{VariableRef}
    α6::Vector{VariableRef}
    α7::Vector{VariableRef}
    α8::Vector{VariableRef}
    α9::Vector{VariableRef}
    α10::Vector{VariableRef}
    α11::Vector{VariableRef}
    α12::Vector{VariableRef}
    α13::Vector{VariableRef}
    α14::Vector{VariableRef}
    m::Model
end
function MasterP(w)
    mas = Model(CPLEX.Optimizer); set_silent(mas)
    # MOI.set(mas, MOI.RawOptimizerAttribute("CPXPARAM_MIP_Strategy_Search"), 1) # conventional branch-and-cut. without dynamic search
    MOI.set(mas, MOI.NumberOfThreads(), 1)
    MOI.NodeCount()
    @variable(mas, y[1:(dt.N["plant"]+dt.N["distribution"])*2], Bin);
    @variable(mas, uij[1:sum(dt2.Mij)], Bin);
    @variable(mas, ujk[1:sum(dt2.Mjk)], Bin);
    @variable(mas, ukl[1:sum(dt2.Mkl)], Bin);
    @variable(mas, θ>=-1000);
    @constraint(mas,[j=1:dt.N["plant"]+dt.N["distribution"]], sum(y[2*(j-1)+1:2*(j-1)+2]) <= 1);
    @constraints(mas, begin
            sum(uij[1:dt2.Mij[1,1]]) <= 1
            [j=2:dt.N["plant"]], sum(uij[sum(dt2.Mij[1,1:j-1])+1:sum(dt2.Mij[1,1:j-1])+dt2.Mij[1,j]]) <= 1
            [i=2:dt.N["supplier"],j=2:dt.N["plant"]],  sum(uij[sum(dt2.Mij[1:i-1,:])+sum(dt2.Mij[i,1:j-1])+1:sum(dt2.Mij[1:i-1,:])+sum(dt2.Mij[i,1:j-1])+dt2.Mij[i,j]])<= 1
            sum(ujk[1:dt2.Mjk[1,1]]) <= 1
            [k=2:dt.N["distribution"]], sum(ujk[sum(dt2.Mjk[1,1:k-1])+1:sum(dt2.Mjk[1,1:k-1])+dt2.Mjk[1,k]]) <= 1
            [j=2:dt.N["plant"],k=2:dt.N["distribution"]],  sum(ujk[sum(dt2.Mjk[1:j-1,:])+sum(dt2.Mjk[j,1:j-1])+1:sum(dt2.Mjk[1:j-1,:])+sum(dt2.Mjk[j,1:j-1])+dt2.Mjk[j,k]]) <= 1
            sum(ukl[1:dt2.Mkl[1,1]]) <= 1
            [l=2:dt.N["customer"]], sum(ukl[sum(dt2.Mkl[1,1:l-1])+1:sum(dt2.Mkl[1,1:l-1])+dt2.Mkl[1,l]]) <= 1
            [k=2:dt.N["distribution"],l=2:dt.N["customer"]],  sum(ukl[sum(dt2.Mkl[1:k-1,:])+sum(dt2.Mkl[k,1:l-1])+1:sum(dt2.Mkl[1:k-1,:])+sum(dt2.Mkl[k,1:l-1])+dt2.Mkl[k,l]])<= 1
    end);
    @constraint(mas, sum(y[1:dt.N["plant"]*2]) <= dt.upl);
    @constraint(mas, sum(y[dt.N["plant"]*2+1:end]) <= dt.udc);
    @objective(mas, Min, w[1]*(sum(dt.c.*y) + sum(dt.gij[i]*uij[i] for i in findnz(dt.gij)[1]) + sum(dt.gjk[i]*ujk[i] for i in findnz(dt.gjk)[1]) + sum(dt.gkl[i].*ukl[i] for i in findnz(dt.gkl)[1]))+ θ );
    return MasterP(y,uij,ujk,ukl,θ,mas);
end

struct distance
    y::Vector{VariableRef}
    uij::Vector{VariableRef}
    ujk::Vector{VariableRef}
    ukl::Vector{VariableRef}
    m::Model
end
function distance()
    dis = Model(CPLEX.Optimizer); set_silent(dis)
    @variable(dis, y[1:(dt.N["plant"]+dt.N["distribution"])*2], Bin);
    @variable(dis, uij[1:sum(dt2.Mij)], Bin);
    @variable(dis, ujk[1:sum(dt2.Mjk)], Bin);
    @variable(dis, ukl[1:sum(dt2.Mkl)], Bin);

    @constraint(dis,[j=1:dt.N["plant"]+dt.N["distribution"]], sum(y[2*(j-1)+1:2*(j-1)+2]) <= 1);
    @constraints(dis, begin
            sum(uij[1:dt2.Mij[1,1]]) <= 1
            [j=2:dt.N["plant"]], sum(uij[sum(dt2.Mij[1,1:j-1])+1:sum(dt2.Mij[1,1:j-1])+dt2.Mij[1,j]]) <= 1
            [i=2:dt.N["supplier"],j=2:dt.N["plant"]],  sum(uij[sum(dt2.Mij[1:i-1,:])+sum(dt2.Mij[i,1:j-1])+1:sum(dt2.Mij[1:i-1,:])+sum(dt2.Mij[i,1:j-1])+dt2.Mij[i,j]])<= 1
            sum(ujk[1:dt2.Mjk[1,1]]) <= 1
            [k=2:dt.N["distribution"]], sum(ujk[sum(dt2.Mjk[1,1:k-1])+1:sum(dt2.Mjk[1,1:k-1])+dt2.Mjk[1,k]]) <= 1
            [j=2:dt.N["plant"],k=2:dt.N["distribution"]],  sum(ujk[sum(dt2.Mjk[1:j-1,:])+sum(dt2.Mjk[j,1:j-1])+1:sum(dt2.Mjk[1:j-1,:])+sum(dt2.Mjk[j,1:j-1])+dt2.Mjk[j,k]]) <= 1
            sum(ukl[1:dt2.Mkl[1,1]]) <= 1
            [l=2:dt.N["customer"]], sum(ukl[sum(dt2.Mkl[1,1:l-1])+1:sum(dt2.Mkl[1,1:l-1])+dt2.Mkl[1,l]]) <= 1
            [k=2:dt.N["distribution"],l=2:dt.N["customer"]],  sum(ukl[sum(dt2.Mkl[1:k-1,:])+sum(dt2.Mkl[k,1:l-1])+1:sum(dt2.Mkl[1:k-1,:])+sum(dt2.Mkl[k,1:l-1])+dt2.Mkl[k,l]])<= 1
    end);
    @constraint(dis, sum(y[1:dt.N["plant"]*2]) <= dt.upl);
    @constraint(dis, sum(y[dt.N["plant"]*2+1:end]) <= dt.udc);
    return distance(y,uij,ujk,ukl,dis);
end
struct MLP
    y::Vector{VariableRef}
    uij::Vector{VariableRef}
    ujk::Vector{VariableRef}
    ukl::Vector{VariableRef}
    θ::VariableRef
    m::Model
end
function MLP(w)
    lp = Model(CPLEX.Optimizer); set_silent(lp)
    @variable(lp, 0<= y[1:(dt.N["plant"]+dt.N["distribution"])*2]<=1);
    @variable(lp, 0<=uij[1:sum(dt2.Mij)]<=1);
    @variable(lp, 0<=ujk[1:sum(dt2.Mjk)]<=1);
    @variable(lp, 0<=ukl[1:sum(dt2.Mkl)]<=1);
    @variable(lp, θ>=-1000);
    @constraint(lp,[j=1:dt.N["plant"]+dt.N["distribution"]], sum(y[2*(j-1)+1:2*(j-1)+2]) <= 1);
    @constraints(lp, begin
            sum(uij[1:dt2.Mij[1,1]]) <= 1
            [j=2:dt.N["plant"]], sum(uij[sum(dt2.Mij[1,1:j-1])+1:sum(dt2.Mij[1,1:j-1])+dt2.Mij[1,j]]) <= 1
            [i=2:dt.N["supplier"],j=2:dt.N["plant"]],  sum(uij[sum(dt2.Mij[1:i-1,:])+sum(dt2.Mij[i,1:j-1])+1:sum(dt2.Mij[1:i-1,:])+sum(dt2.Mij[i,1:j-1])+dt2.Mij[i,j]])<= 1
            sum(ujk[1:dt2.Mjk[1,1]]) <= 1
            [k=2:dt.N["distribution"]], sum(ujk[sum(dt2.Mjk[1,1:k-1])+1:sum(dt2.Mjk[1,1:k-1])+dt2.Mjk[1,k]]) <= 1
            [j=2:dt.N["plant"],k=2:dt.N["distribution"]],  sum(ujk[sum(dt2.Mjk[1:j-1,:])+sum(dt2.Mjk[j,1:j-1])+1:sum(dt2.Mjk[1:j-1,:])+sum(dt2.Mjk[j,1:j-1])+dt2.Mjk[j,k]]) <= 1
            sum(ukl[1:dt2.Mkl[1,1]]) <= 1
            [l=2:dt.N["customer"]], sum(ukl[sum(dt2.Mkl[1,1:l-1])+1:sum(dt2.Mkl[1,1:l-1])+dt2.Mkl[1,l]]) <= 1
            [k=2:dt.N["distribution"],l=2:dt.N["customer"]],  sum(ukl[sum(dt2.Mkl[1:k-1,:])+sum(dt2.Mkl[k,1:l-1])+1:sum(dt2.Mkl[1:k-1,:])+sum(dt2.Mkl[k,1:l-1])+dt2.Mkl[k,l]])<= 1
    end);
    @constraint(lp, sum(y[1:dt.N["plant"]*2]) <= dt.upl);
    @constraint(lp, sum(y[dt.N["plant"]*2+1:end]) <= dt.udc);
    @objective(lp, Min, w[1]*(sum(dt.c.*y) + sum(dt.gij[i]*uij[i] for i in findnz(dt.gij)[1]) + sum(dt.gjk[i]*ujk[i] for i in findnz(dt.gjk)[1]) + sum(dt.gkl[i].*ukl[i] for i in findnz(dt.gkl)[1]))+ θ );
    return MLP(y,uij,ujk,ukl,θ,lp);
end
function DualSP(w)
    sub = direct_model(CPLEX.Optimizer());
    set_optimizer_attribute(sub, "CPX_PARAM_REDUCE", 0); # presolve must be turned off
    MOI.set(sub, MOI.RawOptimizerAttribute("CPXPARAM_MIP_Strategy_Search"), 1) # conventional branch-and-cut
    set_silent(sub) # display is turned off
    MOI.set(sub, MOI.NumberOfThreads(), 1) # of threads

    @variables(sub, begin
        α1[1:dt.N["plant"]*5]
        α2[1:dt.N["distribution"]*5]
        α3[1:dt.N["plant"]*5]
        α4[1:dt.N["distribution"]*5]
    end)
    @variable(sub, α5[1:dt.N["customer"]*5] >= 0);
    @variable(sub, α6[1:dt.N["supplier"]] >= 0);
    @variable(sub, α7[1:dt.N["plant"]*2] >= 0);
    @variable(sub, α8[1:dt.N["distribution"]*2] >= 0);
    @variable(sub, α9[1:sum(dt2.Mij)] >= 0);
    @variable(sub, α10[1:sum(dt2.Mjk)] >= 0);
    @variable(sub, α11[1:sum(dt2.Mkl)] >= 0);
    @variable(sub, α12[1:sum(dt2.Mij)] >= 0);
    @variable(sub, α13[1:sum(dt2.Mjk)] >= 0);
    @variable(sub, α14[1:sum(dt2.Mkl)] >= 0);

    @constraint(sub, [i=1,j=1,m=1:dt2.Mij[i,j],p=1:5], α1[5*(j-1)+p]-α3[5*(j-1)+p]-α6[i]+α9[m]-α12[m] <= w[1]*(dt.a[i,p]+dt2.vij[i][j][m][p])+w[2]*(dt.b[i,p]+dt2.rij[i][j][m][p]))
    @constraint(sub, [i=2:dt.N["supplier"],j=1,m=1:dt2.Mij[i,j],p=1:5], α1[5*(j-1)+p]-α3[5*(j-1)+p]-α6[i]+α9[m+sum(dt2.Mij[1:i-1,:])]-α12[m+sum(dt2.Mij[1:i-1,:])] <= w[1]*(dt.a[i,p]+dt2.vij[i][j][m][p])+w[2]*(dt.b[i,p]+dt2.rij[i][j][m][p]))
    @constraint(sub, [i=1,j=2:dt.N["plant"],m=1:dt2.Mij[i,j],p=1:5], α1[5*(j-1)+p]-α3[5*(j-1)+p]-α6[i]+α9[sum(dt2.Mij[1,1:j-1])+ m]-α12[sum(dt2.Mij[1,1:j-1])+ m] <= w[1]*(dt.a[i,p]+dt2.vij[i][j][m][p])+w[2]*(dt.b[i,p]+dt2.rij[i][j][m][p]))
    @constraint(sub, [i=2:dt.N["supplier"],j=2:dt.N["plant"], m=1:dt2.Mij[i,j],p=1:5],  α1[5*(j-1)+p]-α3[5*(j-1)+p]-α6[i]+α9[sum(dt2.Mij[1:i-1,:])+sum(dt2.Mij[i,1:j-1])+m]-α12[sum(dt2.Mij[1:i-1,:])+sum(dt2.Mij[i,1:j-1])+m] <= w[1]*(dt.a[i,p]+dt2.vij[i][j][m][p])+w[2]*(dt.b[i,p]+dt2.rij[i][j][m][p]))

    @constraint(sub, [j=1,k=1,m=1:dt2.Mjk[j,k],p=1:5], -α1[5*(j-1)+p]+α2[5*(k-1)+p]-α4[5*(k-1)+p]+α10[m]-α13[m] <= w[1]*(dt2.vjk[j][k][m][p])+w[2]*(dt2.rjk[j][k][m][p]) )
    @constraint(sub, [j=2:dt.N["plant"],k=1,m=1:dt2.Mjk[j,k],p=1:5], -α1[5*(j-1)+p]+α2[5*(k-1)+p]-α4[5*(k-1)+p]+α10[sum(dt2.Mjk[1:j-1,:])]-α13[sum(dt2.Mjk[1:j-1,:])+m] <= w[1]*(dt2.vjk[j][k][m][p])+w[2]*(dt2.rjk[j][k][m][p]));
    @constraint(sub, [j=1,k=2:dt.N["distribution"],m=1:dt2.Mjk[j,k],p=1:5], -α1[5*(j-1)+p]+α2[5*(k-1)+p]-α4[5*(k-1)+p]+α10[sum(dt2.Mjk[1,1:k-1])+ m]-α13[sum(dt2.Mjk[1,1:k-1])+ m] <= w[1]*(dt2.vjk[j][k][m][p])+w[2]*(dt2.rjk[j][k][m][p]) )
    @constraint(sub, [j=2:dt.N["plant"],k=2:dt.N["distribution"],m=1:dt2.Mjk[j,k],p=1:5], -α1[5*(j-1)+p]+α2[5*(k-1)+p]-α4[5*(k-1)+p]+α10[sum(dt2.Mjk[1:j-1,:])+sum(dt2.Mjk[j,1:k-1])+m]-α13[sum(dt2.Mjk[1:j-1,:])+sum(dt2.Mjk[j,1:k-1])+m] <= w[1]*(dt2.vjk[j][k][m][p])+w[2]*(dt2.rjk[j][k][m][p]));

    @constraint(sub, [k=1,l=1,m=1:dt2.Mkl[k,l],p=1:5], -α2[5*(k-1)+p]+α5[5*(l-1)+p]+α11[m]-α14[m] <= w[1]*(dt2.vkl[k][l][m][p])+w[2]*(dt2.rkl[k][l][m][p]))
    @constraint(sub, [k=1,l=2:dt.N["customer"],m=1:dt2.Mkl[k,l],p=1:5], -α2[5*(k-1)+p]+α5[5*(l-1)+p]+α11[sum(dt2.Mkl[1,1:l-1])+ m]-α14[sum(dt2.Mkl[1,1:l-1])+ m] <= w[1]*(dt2.vkl[k][l][m][p])+w[2]*(dt2.rkl[k][l][m][p]))
    @constraint(sub, [k=2:dt.N["distribution"],l=1,m=1:dt2.Mkl[k,l],p=1:5], -α2[5*(k-1)+p]+α5[5*(l-1)+p]+α11[sum(dt2.Mkl[1:k-1,:])+m]-α14[sum(dt2.Mkl[1:k-1,:])+m] <= w[1]*(dt2.vkl[k][l][m][p])+w[2]*(dt2.rkl[k][l][m][p]))
    @constraint(sub, [k=2:dt.N["distribution"],l=2:dt.N["customer"],m=1:dt2.Mkl[k,l],p=1:5], -α2[5*(k-1)+p]+α5[5*(l-1)+p]+α11[sum(dt2.Mkl[1:k-1,:])+sum(dt2.Mkl[k,1:l-1])+m]-α14[sum(dt2.Mkl[1:k-1,:])+sum(dt2.Mkl[k,1:l-1])+m] <= w[1]*(dt2.vkl[k][l][m][p])+w[2]*(dt2.rkl[k][l][m][p]))

    @constraint(sub, [j=1:dt.N["plant"],t=1:2,p=1:5], α3[5*(j-1)+p]-α7[2*(j-1)+t] <= w[1]*(dt2.N["vcp"][j][2*(t-1)+p])+w[2]*(dt2.N["vep"][j][2*(t-1)+p]))
    @constraint(sub, [k=1:dt.N["distribution"],t=1:2,p=1:5], α4[5*(k-1)+p]-α8[2*(k-1)+t] <= w[1]*(dt2.N["vcd"][k][2*(t-1)+p])+w[2]*(dt2.N["ved"][k][2*(t-1)+p]))
    return DualSP(α1,α2,α3,α4,α5,α6,α7,α8,α9,α10,α11,α12,α13,α14,sub)
end
###################    Benders with Feasibility Pump    ########################
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
