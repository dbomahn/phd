mutable struct Data_1dim
    file::String; N::Dict{}; d::Array{}; c::Array{}; a::Array{}; e::Array{}; gij::SparseVector{}; gjk::SparseVector{}; gkl::SparseVector{};
    vij::Array{}; vjk::Array{}; vkl::Array{}; Vij::SparseVector{}; Vjk::SparseVector{}; Vkl::SparseVector{}; Mij::Array{}; Mjk::Array{}; Mkl::Array{};
    b::Array{}; q::Array{}; rij::Array{}; rjk::Array{}; rkl::Array{}; upl::Int; udc::Int; bigM::Int
    function Data1(file)
        dt1 = readdlm(file);
        notafile = readdlm("/home/ak121396/Desktop/instances/SCND/Notations.txt", '=');
        # notafile = readdlm("F:/scnd/Notations.txt", '=');
        # notafile = readdlm("/home/k2g00/k2g3475/scnd1/Notations.txt", '=');
        nota = notafile[1:end,1];  N= Dict();

        for i=1:length(nota)-1
            id1 = findall(x->x==nota[i], dt1)[1][1];
            id2 = findall(x->x==nota[i+1], dt1)[1][1];
            if id2-id1<3
                tmp = filter(x->x!="",  dt1[id1+(id2-id1-1),:])
                if length(tmp)<2
                    N[nota[i]] = tmp[1];
                else
                    N[nota[i]] = tmp;
                end
            else
                W = []
                for x=id1+1:id1+(id2-id1-1)
                    tmp = filter(x->x!="", dt1[x,:]);
                    append!(W,tmp)
                    # push!(W,tmp);
                end
                # tmp = [filter(x->x!="", dt1[x,:]) for x in id1+1:id1+(id2-id1-1)]
                N[nota[i]] = W;
            end
        end
        d = reshape(N["demand"],5,N["customer"])'; c = append!(N["fcp"],N["fcd"]);
        # a = reshape(N["vcs"],5,N["supplier"])';
         # e = append!(N["vcp"],N["vcd"]);
        gij = sparse(N["fixedcostModesp"]); gjk = sparse(N["fixedcostModepd"]); gkl = sparse(N["fixedcostModedc"]);
        # vij = N["tcp"]; vjk = N["tcd"]; vkl = N["tcc"];
        vij = [];
        for i=1:N["supplier"]
            tmp = [];
            for j=1:N["plant"]
                tmp2 = []
                if j==1
                    for m=1:Mij[i,1]
                        push!(tmp2, N["tcp"][i][5*(m-1)+1:5*(m-1)+5])
                    end
                else
                    for m=1:Mij[i,j]
                        push!(tmp2, N["tcp"][i][5*sum(Mij[i,1:j-1])+5*(m-1)+1:5*sum(Mij[i,1:j-1])+5*(m-1)+5]);
                    end
                end
                push!(tmp,tmp2);
            end
            push!(vij,tmp);
        end
        vjk = [];
        for j=1:N["plant"]
            tmp = [];
            for k=1:N["distribution"]
                tmp2 = []
                if k==1
                    for m=1:Mjk[j,1]
                        push!(tmp2, N["tcd"][j][5*(m-1)+1:5*(m-1)+5])
                    end
                else
                    for m=1:Mjk[j,k]
                        push!(tmp2, N["tcd"][j][5*sum(Mjk[j,1:k-1])+5*(m-1)+1:5*sum(Mjk[j,1:k-1])+5*(m-1)+5]);
                    end
                end
                push!(tmp,tmp2);
            end
            push!(vjk,tmp);
        end
        vkl = [];
        for k=1:N["distribution"]
            tmp = [];
            for l=1:N["customer"]
                tmp2 = []
                if l==1
                    for m=1:Mkl[k,1]
                        push!(tmp2, N["tcc"][k][5*(m-1)+1:5*(m-1)+5])
                    end
                else
                    for m=1:Mkl[k,l]
                        push!(tmp2, N["tcc"][k][5*sum(Mkl[k,1:l-1])+5*(m-1)+1:5*sum(Mkl[k,1:l-1])+5*(m-1)+5]);
                    end
                end
                push!(tmp,tmp2);
            end
            push!(vkl,tmp);
        end

        Vij = sparse(N["LcapacityModesp"]); Vjk = sparse(N["LcapacityModepd"]); Vkl =  sparse(N["LcapacityModedc"]);
        Mij = transpose(reshape(N["ModeIJ"], (N["plant"],N["supplier"])));
        Mjk = transpose(reshape(N["ModeJK"], (N["distribution"],N["plant"])));
        Mkl = transpose(reshape(N["ModeKL"], (N["customer"],N["distribution"])));
        b = reshape(N["ves"],5,N["supplier"])';  q = append!(N["vep"],N["ved"]);
        # rij = N["cep"]; rjk = N["ced"]; rkl = N["cec"];
        rij = [];
        for i=1:N["supplier"]
            tmp = [];
            for j=1:N["plant"]
                tmp2 = []
                if j==1
                    for m=1:Mij[i,1]
                        push!(tmp2, N["cep"][i][5*(m-1)+1:5*(m-1)+5])
                    end
                else
                    for m=1:Mij[i,j]
                        push!(tmp2, N["cep"][i][5*sum(Mij[i,1:j-1])+5*(m-1)+1:5*sum(Mij[i,1:j-1])+5*(m-1)+5]);
                    end
                end
                push!(tmp,tmp2);
            end
            push!(rij,tmp);
        end
        rjk = [];
        for j=1:N["plant"]
            tmp = [];
            for k=1:N["distribution"]
                tmp2 = []
                if k==1
                    for m=1:Mjk[j,1]
                        push!(tmp2, N["ced"][j][5*(m-1)+1:5*(m-1)+5])
                    end
                else
                    for m=1:Mjk[j,k]
                        push!(tmp2, N["ced"][j][5*sum(Mjk[j,1:k-1])+5*(m-1)+1:5*sum(Mjk[j,1:k-1])+5*(m-1)+5]);
                    end
                end
                push!(tmp,tmp2);
            end
            push!(rjk,tmp);
        end
        rkl = [];
        for k=1:N["distribution"]
            tmp = [];
            for l=1:N["customer"]
                tmp2 = []
                if l==1
                    for m=1:Mkl[k,1]
                        push!(tmp2, N["cec"][k][5*(m-1)+1:5*(m-1)+5])
                    end
                else
                    for m=1:Mkl[k,l]
                        push!(tmp2, N["cec"][k][5*sum(Mkl[k,1:l-1])+5*(m-1)+1:5*sum(Mkl[k,1:l-1])+5*(m-1)+5]);
                    end
                end
                push!(tmp,tmp2);
            end
            push!(rkl,tmp);
        end
        upl = N["upperpants"]; udc = N["upperdistribution"]; bigM = sum(N["demand"])

        new(file,N,d,c,a,e,gij,gjk,gkl,vij,vjk,vkl,Vij,Vjk,Vkl,Mij,Mjk,Mkl,b,q,rij,rjk,rkl,upl,udc,bigM);
    end
end
struct MasterP
    y::Vector{VariableRef}
    uij::Vector{VariableRef}
    ujk::Vector{VariableRef}
    ukl::Vector{VariableRef}
    θ::VariableRef
    m::Model
end

function MasterP(w)
    mas = Model(CPLEX.Optimizer); set_silent(mas)
    MOI.set(mas, MOI.RawOptimizerAttribute("CPXPARAM_MIP_Strategy_Search"), 1) # conventional branch-and-cut. without dynamic search
    MOI.set(mas, MOI.NumberOfThreads(), 1)
    MOI.NodeCount()
    @variable(mas, y[1:(dt1.N["plant"]+dt1.N["distribution"])*2], Bin)
    @variable(mas, uij[1:sum(dt1.Mij)], Bin)
    @variable(mas, ujk[1:sum(dt1.Mjk)], Bin)
    @variable(mas, ukl[1:sum(dt1.Mkl)], Bin)
    @constraint(mas,[j=1:dt1.N["plant"]+dt1.N["distribution"]], sum(y[2*(j-1)+1:2*(j-1)+2]) <= 1);
    @constraints(mas,begin
            sum(uij[1:dt1.Mij[1,1]]) <= 1
            [j=2:dt1.N["plant"]], sum(uij[sum(dt1.Mij[1,1:j-1])+1:sum(dt1.Mij[1,1:j-1])+dt1.Mij[1,j]]) <= 1
            [i=2:dt1.N["supplier"],j=2:dt1.N["plant"]],  sum(uij[sum(dt1.Mij[1:i-1,:])+sum(dt1.Mij[i,1:j-1])+1:sum(dt1.Mij[1:i-1,:])+sum(dt1.Mij[i,1:j-1])+dt1.Mij[i,j]])<= 1
            sum(ujk[1:dt1.Mjk[1,1]]) <= 1
            [k=2:dt1.N["distribution"]], sum(ujk[sum(dt1.Mjk[1,1:k-1])+1:sum(dt1.Mjk[1,1:k-1])+dt1.Mjk[1,k]]) <= 1
            [j=2:dt1.N["plant"],k=2:dt1.N["distribution"]],  sum(ujk[sum(dt1.Mjk[1:j-1,:])+sum(dt1.Mjk[j,1:j-1])+1:sum(dt1.Mjk[1:j-1,:])+sum(dt1.Mjk[j,1:j-1])+dt1.Mjk[j,k]]) <= 1
            sum(ukl[1:dt1.Mkl[1,1]]) <= 1
            [l=2:dt1.N["customer"]], sum(ukl[sum(dt1.Mkl[1,1:l-1])+1:sum(dt1.Mkl[1,1:l-1])+dt1.Mkl[1,l]]) <= 1
            [k=2:dt1.N["distribution"],l=2:dt1.N["customer"]],  sum(ukl[sum(dt1.Mkl[1:k-1,:])+sum(dt1.Mkl[k,1:l-1])+1:sum(dt1.Mkl[1:k-1,:])+sum(dt1.Mkl[k,1:l-1])+dt1.Mkl[k,l]])<= 1
    end);
    @constraint(mas, sum(y[1:dt1.N["plant"]*2]) <= dt1.upl);
    @constraint(mas, sum(y[dt1.N["plant"]*2+1:end]) <= dt1.udc);
    @objective(mas, Min, w[1]*(sum(dt1.c.*y) + sum(dt1.gij[i]*uij[i] for i in findnz(dt1.gij)[1]) + sum(dt1.gjk[i]*ujk[i] for i in findnz(dt1.gjk)[1]) + sum(dt1.gkl[i].*ukl[i] for i in findnz(dt1.gkl)[1])) );
    return MasterLP(y,uij,ujk,ukl,mas)
end
mlp = MasterLP(w)
# function getobjval(x,C)
#     obj = sum(dt.c[j][t]*y[j,t] for j=1:dt.N["plant"]+dt.N["distribution"] for t=1:2)+
#     sum(dt.gij[i][j][m]*uij[i,j,m] for i=1:dt.N["supplier"] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j]) +
#     sum(dt.gjk[j][k][m]*ujk[j,k,m] for j=1:dt.N["plant"] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k]) +
#     sum(dt.gkl[k][l][m]*ukl[k,l,m] for k=1:dt.N["distribution"] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])) + θ
#     return obj
# end
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
    @variable(sub, α7[dt.N["plant"],1:2] >= 0);
    @variable(sub, α8[1:dt.N["distribution"]*2] >= 0);
    @variable(sub, α9[1:sum(dt.Mij)] >= 0);
    @variable(sub, α10[1:sum(dt.Mjk)] >= 0);
    @variable(sub, α11[1:sum(dt.Mkl)] >= 0);
    @variable(sub, α12[1:sum(dt.Mij)] >= 0);
    @variable(sub, α13[1:sum(dt.Mjk)] >= 0);
    @variable(sub, α14[1:sum(dt.Mkl)] >= 0);

    @constraint(sub, [i=1,j=1:dt.N["plant"],m=1:dt.Mij[i,j],p=1:5], α1[5*(j-1)+p]-α3[5*(j-1)+p]-α6[i]+α9[sum(dt.Mij[1,1:j-1])+ m]-α12[sum(dt.Mij[1,1:j-1])+ m] <= w[1]*(dt.N["vcs"][i][p]+dt.vij[i][j][m][p])+w[2]*(dt.b[i,p]+dt.rij[i][j][m][p]));
    @constraint(sub, [i=2:dt.N["supplier"],j=1:dt.N["plant"], m=1:dt.Mij[i,j],p=1:5],  α1[5*(j-1)+p]-α3[5*(j-1)+p]-α6[i]+α9[sum(dt.Mij[1:i-1,:])+sum(Mij[i,1:j-1])+m]-α12[sum(dt.Mij[1:i-1,:])+sum(Mij[i,1:j-1])+m] <= w[1]*(dt.a[i,p]+dt.vij[i][j][m][p])+w[2]*(dt.b[i,p]+dt.rij[i][j][m][p]))

    @constraint(sub, [j=1,k=1:dt.N["distribution"],m=1:dt.Mjk[j,k],p=1:5], -α1[5*(j-1)+p]+α2[5*(k-1)+p]-α4[5*(k-1)+p]+α10[sum(dt.Mjk[1,1:k-1])+ m]-α13[sum(dt.Mjk[1,1:k-1])+ m] <= w[1]*(dt.vjk[j][k][m][p])+w[2]*(dt.rjk[j][k][m][p]));
    @constraint(sub, [j=2:dt.N["plant"],k=1:dt.N["distribution"],m=1:dt.Mjk[j,k],p=1:5], -α1[5*(j-1)+p]+α2[5*(k-1)+p]-α4[5*(k-1)+p]+α10[sum(dt.Mjk[1:j-1,:])+sum(Mjk[j,1:k-1])+m]-α13[sum(dt.Mjk[1:j-1,:])+sum(Mjk[j,1:k-1])+m] <= w[1]*(dt.vjk[j][k][m][p])+w[2]*(dt.rjk[j][k][m][p]));

    @constraint(sub, [k=1,l=1:dt.N["customer"],m=1:dt.Mkl[k,l],p=1:5], -α2[5*(k-1)+p]+α5[5*(l-1)+p]+α11[sum(dt.Mkl[1,1:l-1])+ m]-α14[sum(dt.Mkl[1,1:l-1])+ m] <= w[1]*(dt.vkl[k][l][m][p])+w[2]*(dt.rkl[k][l][m][p]));
    @constraint(sub, [k=2:dt.N["distribution"],l=1:dt.N["customer"],m=1:dt.Mkl[k,l],p=1:5], -α2[5*(k-1)+p]+α5[5*(l-1)+p]+α11[sum(dt.Mkl[1:k-1,:])+sum(Mkl[k,1:l-1])+m]-α14[dt.Mkl[1:k-1,:])+sum(Mkl[k,1:l-1])+m] <= w[1]*(dt.vkl[k][l][m][p])+w[2]*(dt.rkl[k][l][m][p]));

    @constraint(sub, [j=1:dt.N["plant"],t=1:2,p=1:5], α3[5*(j-1)+p]-α7[2*(j-1)+t] <= w[1]*(dt.N["vcp"][j][2*(t-1)+p])+w[2]*(dt.N["vep"][j][2*(t-1)+p]));
    @constraint(sub, [k=1:dt.N["distribution"],t=1:2,p=1:5], α4[5*(k-1)+p]-α8[2*(k-1)+t] <= w[1]*(dt.N["vcd"][k][2*(t-1)+p])+w[2]*(dt.N["ved"][k][2*(t-1)+p]));
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
function FBcheck(yr,u1r,u2r,u3r)
    # for k=1:length(yr)
    #     JuMP.fix(mp.y[k],yr[k]; force=true)
    # end
    # for k=1:length(u1r)
    #     JuMP.fix(mp.uij[k],u1r[k]; force=true)
    # end
    # for k=1:length(u2r)
    #     JuMP.fix(mp.ujk[k],u2r[k]; force=true)
    # end
    # for k=1:length(u3r)
    #     JuMP.fix(mp.ukl[k],u3r[k]; force=true)
    # end
    mp.y = yr; mp.uij=u1r; mp.ujk=u2r; mp.ukl=u3r;
    optimize!(mp.m)
    if termination_status(mp.m) == MOI.OPTIMAL
        return true
    else
        return false
    end
end
function fbsearch(yr,u1r,u2r,u3r) #solveLP
    idy_0 = findall(k->k==0, yr)
    idy_1 = findall(k->k==1, yr)
    idu1_0 = findall(k->k==0, u1r)
    idu1_1 = findall(k->k==1, u1r)
    idu2_0 = findall(k->k==0, u2r)
    idu2_1 = findall(k->k==1, u2r)
    idu3_0 = findall(k->k==0, u3r)
    idu3_1 = findall(k->k==1, u3r)
    @objective( masterlp.m, Min, sum(masterlp.y[i] for i in idy_0) + sum(1-(masterlp.y[j]) for j in idy_1) +
        sum(masterlp.uij[i] for i in idu1_0) + sum(1-(masterlp.uij[j]) for j in idu1_1)+
        sum(masterlp.ujk[i] for i in idu2_0) + sum(1-(masterlp.ujk[j]) for j in idu2_1)+
        sum(masterlp.ukl[i] for i in idu3_0) + sum(1-(masterlp.ukl[j]) for j in idu3_1)
    )
    optimize!(masterlp.m)
    if termination_status(masterlp.m) == MOI.OPTIMAL
        return JuMP.value.(masterlp.y),JuMP.value.(masterlp.uij),JuMP.value.(masterlp.ujk),JuMP.value.(masterlp.ukl)
    else
        return 0,0,0,0
    end
end
function FP(yt,u1t,u2t,u3t,TL)
	sol = []; SearchDone = false;	Tabu = []; newsol=0; k=1; t0=time(); iter=0; Max_iter = 50 #Y = [];
    # while candlist != [] &&  time()-t0 < TL && k < length(candX)+1

    while iter<Max_iter && SearchDone == false
        yr = round.(Int,yt);u1r = round.(Int,u1t);u2r = round.(Int,u2t);u3r = round.(Int,u3t)
        if ( (FBcheck(yr,u1r,u2r,u3r) == true) && [yr;u1r;u2r;u3r] ∉sol)
			push!(sol,[yr,u1r,u2r,u3r]); newsol+=1;
            @show SearchDone = true            #push!(Y,getobjval(x_r,C))
        else
            if [yr;u1r;u2r;u3r] ∈ Tabu
                yr = flipoper(Tabu,yt,yr); u1r = flipoper(Tabu,u1t,u1r); u2r = flipoper(Tabu,u2t,u2r); u3r = flipoper(Tabu,u3t,u3r)
                if any(i->i==[], [yr,u1r,u2r,u3r])
                    @show SearchDone = true
                else
                    if ( (FBcheck(yr,u1r,u2r,u3r) == true) && [yr;u1r;u2r;u3r] ∉sol)
						push!(sol,[yr,u1r,u2r,u3r]); newsol+=1; SearchDone = true            #push!(Y,getobjval(x_r,C))
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
                    @show SearchDone = true
                end
            end
        end
		iter+=1
    end
    return sol
end
function masterFP(mp::build_master)
    optimize!(mp.m)    # end
    yt = value.(mp.y)
    u1t = value.(mp.uij); u2t = value.(mp.ujk); u3t = value.(mp.ukl)
    sol = FP(yt,u1t,u2t,u3t,60)
    return sol
end
optimize!(mp.m)    # end
yt = value.(mp.y)
round.(Int,yt)

sol = masterFP(mp)
