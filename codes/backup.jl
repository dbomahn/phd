scnd = Model(CPLEX.Optimizer); set_silent(scnd)
@variable(scnd, 0<=y[1:dt.N["plant"]+dt.N["distribution"],1:2]<=1  );
@variable(scnd, 0<=uij[i=1:dt.N["supplier"],j=1:dt.N["plant"], m=1:dt.Mij[i,j]] <=1 )
@variable(scnd, 0<=ujk[j=1:dt.N["plant"],k=1:dt.N["distribution"],m=1:dt.Mjk[j,k]] <=1);
@variable(scnd, 0<=ukl[k=1:dt.N["distribution"],l=1:dt.N["customer"],m=1:dt.Mkl[k,l]] <=1);

##############################  IP   #####################################
# @variable(scnd, y[1:dt.N["plant"]+dt.N["distribution"],1:2], Bin  );
# @variable(scnd, uij[i=1:dt.N["supplier"],j=1:dt.N["plant"], m=1:dt.Mij[i,j]] , Bin);
# @variable(scnd, ujk[j=1:dt.N["plant"],k=1:dt.N["distribution"],m=1:dt.Mjk[j,k]] , Bin);
# @variable(scnd, ukl[k=1:dt.N["distribution"],l=1:dt.N["customer"],m=1:dt.Mkl[k,l]] , Bin);

@variable(scnd, 0<= xij[i=1:dt.N["supplier"],j=1:dt.N["plant"],m=1:dt.Mij[i,j],1:5] );
@variable(scnd, 0<= xjk[j=1:dt.N["plant"],k=1:dt.N["distribution"],m=1:dt.Mjk[j,k],1:5] );
@variable(scnd, 0<= xkl[k=1:dt.N["distribution"],l=1:dt.N["customer"],m=1:dt.Mkl[k,l],1:5] );
@variable(scnd, 0<= h[1:dt.N["plant"]+dt.N["distribution"],1:5,1:2] );

#a_ip*x_ijmp
exa = AffExpr(0);
for i=1:dt.N["supplier"]
    for j=1:dt.N["plant"]
        for m=1:dt.Mij[i,j]
            for p=1:5
                # if is_valid(scnd,xij[i,j,m,p])==true
                add_to_expression!(exa,sum(dt.N["vcs"][i][p]*xij[i,j,m,p]) );
                # end
            end
        end
    end
end
#g_ijm*u_ijm expression
exg = AffExpr(0);
for i=1:dt.N["supplier"]
    idx = 1;
    for j=1:dt.N["plant"]
        # if dt.Mij[i,j]<dt.m
        #     for d=1:dt.m-dt.Mij[i,j]
        #         delete(scnd,uij[i,j,dt.m-d+1])
        #     end
        # end
        for m=1:dt.Mij[i,j]
            add_to_expression!(exg, dt.gij[i][idx]*uij[i,j,m]);
            idx+=1
        end
    end
end
for j=1:dt.N["plant"]
    idx = 1;
    for k=1:dt.N["distribution"]
        # if dt.Mjk[j,k]<dt.m
        #     for d=1:dt.m-dt.Mjk[j,k]
        #         delete(scnd,ujk[j,k,dt.m-d+1])
        #     end
        # end
        for m=1:dt.Mjk[j,k]
            add_to_expression!(exg, dt.gjk[j][idx]*ujk[j,k,m]);
            idx+=1
        end
    end
end
for k=1:dt.N["distribution"]
    idx = 1;
    for l=1:dt.N["customer"]
        # if dt.Mkl[k,l]<dt.m
        #     for d=1:dt.m-dt.Mkl[k,l]
        #         delete(scnd,ukl[k,l,dt.m-d+1])
        #     end
        # end
        for m=1:dt.Mkl[k,l]
            add_to_expression!(exg,dt.gkl[k][idx]*ukl[k,l,m]);
            idx+=1
        end
    end
end

#v_ijmp*x_ijmp expression
exv = AffExpr(0);
for i=1:dt.N["supplier"]
    idx = 1;
    for j=1:dt.N["plant"]
        # if dt.Mij[i,j]<dt.m
        #     for d=1:dt.m-dt.Mij[i,j]
        #         delete(scnd,xij[i,j,dt.m-d+1,:])
        #     end
        # end
        for m=1:dt.Mij[i,j]
            add_to_expression!(exv,sum(dot.(dt.N["tcp"][i][idx:idx+4],(xij[i,j,m,p] for p=1:5))))#*sqrt((dt.N["pointsupplier"][1][i]-dt.N["pointplant"][1][j])^2+(dt.N["pointsupplier"][2][i]-dt.N["pointplant"][2][j])^2)) );
            idx+=5
        end
    end
end
for j=1:dt.N["plant"]
    idx = 1;
    for k=1:dt.N["distribution"]
        # if dt.Mjk[j,k]<dt.m
        #     for d=1:dt.m-dt.Mjk[j,k]
        #         delete(scnd,xjk[j,k,dt.m-d+1,:])
        #     end
        # end
        for m=1:dt.Mjk[j,k]
            add_to_expression!(exv,sum(dot.(dt.N["tcd"][j][idx:idx+4],xjk[j,k,m,p] for p=1:5)))#*sqrt((dt.N["pointplant"][1][j]-dt.N["pointdistribution"][1][k])^2+(dt.N["pointplant"][2][j]-dt.N["pointdistribution"][2][k])^2)) );
            idx+=5
        end
    end
end
for k=1:dt.N["distribution"]
    idx = 1;
    for l=1:dt.N["customer"]
        # if dt.Mkl[k,l]<dt.m
        #     for d=1:dt.m-dt.Mkl[k,l]
        #         delete(scnd,xkl[k,l,dt.m-d+1,1:5])
        #     end
        # end
        for m=1:dt.Mkl[k,l]
            add_to_expression!(exv,sum(dot.(dt.N["tcc"][k][idx:idx+4],xkl[k,l,m,p] for p=1:5)))#*sqrt((dt.N["pointdistribution"][1][k]-dt.N["pointCustmoer"][1][l])^2+(dt.N["pointdistribution"][2][k]-dt.N["pointCustmoer"][2][l])^2)) );
            idx+=5
        end
    end
end
#b_ip*x_ijmp
exb = AffExpr(0);
for i=1:dt.N["supplier"]
    for j=1:dt.N["plant"]
        for m=1:dt.Mij[i,j]
            for p=1:5
                # if is_valid(scnd,xij[i,j,m,p])==true
                add_to_expression!(exb,sum(dt.b[i,p]*xij[i,j,m,p]) );
                # end
            end
        end
    end
end

exr = AffExpr(0);
for i=1:dt.N["supplier"]
    idx = 1;
    for j=1:dt.N["plant"]
        for m=1:dt.Mij[i,j]
            add_to_expression!(exr,sum(dot.(dt.N["cep"][i][idx:idx+4],xij[i,j,m,p] for p=1:5)))#*sqrt((dt.N["pointsupplier"][1][i]-dt.N["pointplant"][1][j])^2+(dt.N["pointsupplier"][2][i]-dt.N["pointplant"][2][j])^2)) );
            idx+=5
        end
    end
end
for j=1:dt.N["plant"]
    idx = 1;
    for k=1:dt.N["distribution"]
        for m=1:dt.Mjk[j,k]
            add_to_expression!(exr,sum(dot.(dt.N["ced"][j][idx:idx+4],xjk[j,k,m,p] for p=1:5)))#*sqrt((dt.N["pointplant"][1][j]-dt.N["pointdistribution"][1][k])^2+(dt.N["pointplant"][2][j]-dt.N["pointdistribution"][2][k])^2)) );
            idx+=5
        end
    end
end
for k=1:dt.N["distribution"]
    idx = 1;
    for l=1:dt.N["customer"]
        for m=1:dt.Mkl[k,l]
            add_to_expression!(exr,sum(dot.(dt.N["cec"][k][idx:idx+4],xkl[k,l,m,p] for p=1:5)))#*sqrt((dt.N["pointdistribution"][1][k]-dt.N["pointCustmoer"][1][l])^2+(dt.N["pointdistribution"][2][k]-dt.N["pointCustmoer"][2][l])^2)) );
            idx+=5
        end
    end
end
#1st obj
@constraint(scnd, obj1, sum(dt.c[j][t]*y[j,t] for j=1:1:dt.N["plant"]+dt.N["distribution"] for t=1:2) + exa + sum(dt.e[j][(p-1)*2+t]*h[j,p,t] for j=1:dt.N["plant"] + dt.N["distribution"] for p=1:5 for t=1:2) + exg + exv <=0);
# @objective(scnd, Min, sum(dt.c[j][t]*y[j,t] for j=1:1:dt.N["plant"]+dt.N["distribution"] for t=1:2) + exg +exa + exv + sum(dt.e[j][(p-1)*2+t]*h[j,p,t] for j=1:(dt.N["plant"]+dt.N["distribution"]) for p=1:5 for t=1:2));
#2nd obj
@constraint(scnd, obj2, exb+sum(dt.q[j][(p-1)*2+t]*h[j,p,t] for j=1:dt.N["plant"]+dt.N["distribution"] for p=1:5 for t=1:2) +exr <=0);
# @objective(scnd,Min,exb + exr + sum(dt.q[j][(p-1)*2+t]*h[j,p,t] for j=1:dt.N["plant"]+dt.N["distribution"] for p=1:5 for t=1:2) );

########### constraint 3 ###############
@constraints(scnd, begin
    [j=1:dt.N["plant"],p=1:5], sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j]) == sum(xjk[j,k,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mjk[j,k])
    [k=1:dt.N["distribution"],p=1:5], sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k]) == sum(xkl[k,l,m,p] for l=1:dt.N["customer"] for m=1:dt.Mkl[k,l])
end);
########### constraint 4-6 #############
@constraints(scnd, begin
    [j=1:dt.N["plant"],p=1:5], sum(h[j,p,:]) == sum(xij[i,j,m,p] for i=1:dt.N["supplier"] for m=1:dt.Mij[i,j])
    [k=1:dt.N["distribution"],p=1:5], sum(h[k+dt.N["plant"],p,:] ) == sum(xjk[j,k,m,p] for j=1:dt.N["plant"] for m=1:dt.Mjk[j,k])
    [l=1:dt.N["customer"],p=1:5], sum(xkl[k,l,m,p] for k=1:dt.N["distribution"] for m=1:dt.Mkl[k,l]) >= dt.d[l][p]
end )
########### constraint 7-9 #############
@constraint(scnd,[i=1:dt.N["supplier"]], sum(xij[i,j,m,p] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] for p=1:5) <= dt.N["cas"][i]);
@constraint(scnd,[j=1:dt.N["plant"]+dt.N["distribution"], t=1:2], sum(h[j,1:5,t]) <= [dt.N["cap"];dt.N["cad"]][j]*y[j,t]);
@constraint(scnd,[j=1:dt.N["plant"]+dt.N["distribution"]], sum(y[j,:]) <= 1);
########### constraint 10 #############
@constraint(scnd,[i=1:dt.N["supplier"],j=1:dt.N["plant"]], sum(uij[i,j,m] for m=1:dt.Mij[i,j]) <= 1);
@constraint(scnd,[j=1:dt.N["plant"],k=1:dt.N["distribution"]], sum(ujk[j,k,m] for m=1:dt.Mjk[j,k]) <= 1);
@constraint(scnd,[k=1:dt.N["distribution"],l=1:dt.N["customer"]], sum(ukl[k,l,m] for m=1:dt.Mkl[k,l]) <= 1);
########### constraint 12 #############
@constraint(scnd,[i=1:dt.N["supplier"], j=1:dt.N["plant"], m=1:dt.Mij[i,j]], sum(xij[i,j,m,p] for p=1:5) >= dt.Vij[i][j][m]*uij[i,j,m] );
@constraint(scnd,[j=1:dt.N["plant"], k=1:dt.N["distribution"], m=1:dt.Mjk[j,k]], sum(xjk[j,k,m,p] for p=1:5) >= dt.Vjk[j][k][m]*ujk[j,k,m]);
@constraint(scnd,[k=1:dt.N["distribution"], l=1:dt.N["customer"], m=1:dt.Mkl[k,l]], sum(xkl[k,l,m,p] for p=1:5) >= dt.Vkl[k][l][m]*ukl[k,l,m]);
########### constraint 13-14 #############
@constraint(scnd, sum(y[j,t] for j=1:dt.N["plant"] for t=1:dt.m) <= dt.upl);
@constraint(scnd, sum(y[j,t] for j=dt.N["plant"]+1:dt.N["distribution"]+dt.N["plant"] for t=1:dt.m) <= dt.udc);
########### goods can be delivered only by chosen transportation mode #############
BigM = 10^(20)
@constraint(scnd,[i=1:dt.N["supplier"], j=1:dt.N["plant"], m=1:dt.Mij[i,j]], sum(xij[i,j,m,p] for p=1:5) <= BigM*uij[i,j,m] );
@constraint(scnd,[j=1:dt.N["plant"], k=1:dt.N["distribution"], m=1:dt.Mjk[j,k]], sum(xjk[j,k,m,p] for p=1:5) <= BigM*ujk[j,k,m]);
@constraint(scnd,[k=1:dt.N["distribution"], l=1:dt.N["customer"], m=1:dt.Mkl[k,l]], sum(xkl[k,l,m,p] for p=1:5) <= BigM*ukl[k,l,m]);
########### Suppliers Availibility  ################
for i=1:dt.N["supplier"]
    for p=1:5
        if dt.N["SuppliersAvailibility"][i][p]==0
            @constraint(scnd, sum(xij[i,j,m,p] for j=1:dt.N["plant"] for m=1:dt.Mij[i,j] )==0)
        end
    end
end
write_to_file(scnd , "/home/k2g00/k2g3475/scnd/lp/"*file[36:end]*".lp")

# optimize!(scnd)
# termination_status(scnd)
# objective_value(scnd)
# #
# value.(uij)
# sum(value.(uij))

# 7.6471829e7+140000+130000+5.307581191999919e6+
# 1.6701908309100384e7+3.1115610282100677e7+2.8889434424200263e7+
# dot(dtt.C[1,:][1+54+145+269+2549+725+1345+12745:54+145+269+2549+725+1345+12745+270],value.(h)[:])
# (dtt.C[1,:][1+54+145+269+2549+725+1345:54+145+269+2549+725+1345+12745])
# dtt.C[1,:][1+54:54+145]
