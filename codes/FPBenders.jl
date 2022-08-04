cd("C:/Users/AK121396/Desktop/ProjectBenders")


f1mp = f1Master(); f1dsp = f1DualSP(); optcuts = []; feasicuts = [];
MOI.set(f1mp.m, MOI.LazyConstraintCallback(), lazy_callback);
# MOI.set(f1mp.m, MOI.UserCutCallback(), rootfrac_callback)
set_time_limit_sec(f1mp.m, 400)
optimize!(f1mp.m)
# solve_time(f1mp.m)
length(optcuts)
length(feasicuts),length(unique(feasicuts))
w = [0.5,0.5]
# termination_status(molp)
# objective_value(molp)
# yt = value.(molp.y); utij = value.(molp.uij); utjk = value.(molp.ujk); utkl = value.(molp.ukl);
# sol = FP(yt,utij,utjk,utkl)
1

function benders_decomposition()
    # ocuts = []; fcuts = []; yArchiv = [];

    optimize!(mp.m); st = termination_status(mp.m)

    while (st == MOI.INFEASIBLE) || (st == MOI.OPTIMAL)
        optimize!(molp.m); #optimize!(mp.m)
        yt = value.(molp.y); utij = value.(molp.uij); utjk = value.(molp.ujk); utkl = value.(molp.ukl);
        sol = FP(yt,utij,utjk,utkl);
        if sol != nothing  # feasible sol found by FP
            push!(feasisol,sol)
            yb = sol[1]; ubij = sol[2]; ubjk = sol[3]; ubkl = sol[4];
            subp = solve_dsp(dsp,yb,ubij,ubjk,ubkl)
            if subp.res == :OptimalityCut
                # @info "FP == Optimality cut found"
                if round(value(mp.θ); digits=4) ≥ round(subp.obj; digits=4)
                    break
                else
                    @constraint( mp.m, mp.θ ≥ sum(subp.α5[5*(l-1)+p]*dt.d[l,p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
                        sum(dt.N["cap"][j]*mp.y[2*(j-1)+t]*subp.α7[2*(j-1)+t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*mp.y[2*(dt.N["plant"]+k-1)+t]*subp.α8[2*(k-1)+t] for k=1:dt.N["distribution"] for t=1:2)+
                        sum(dt.Vij[i]*mp.uij[i]*subp.α9[i] for i in findnz(dt.Vij)[1])+
                        sum(dt.Vjk[i]*mp.ujk[i]*subp.α10[i] for i in findnz(dt.Vjk)[1])-
                        sum(dt.bigM.*mp.uij.*subp.α12) - sum(dt.bigM.*mp.ujk.*subp.α13) - sum(dt.bigM.*mp.ukl.*subp.α14))
                    @constraint( molp.m, molp.θ ≥ sum(subp.α5[5*(l-1)+p]*dt.d[l,p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
                        sum(dt.N["cap"][j]*molp.y[2*(j-1)+t]*subp.α7[2*(j-1)+t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*mp.y[2*(dt.N["plant"]+k-1)+t]*subp.α8[2*(k-1)+t] for k=1:dt.N["distribution"] for t=1:2)+
                        sum(dt.Vij[i]*molp.uij[i]*subp.α9[i] for i in findnz(dt.Vij)[1])+
                        sum(dt.Vjk[i]*molp.ujk[i]*subp.α10[i] for i in findnz(dt.Vjk)[1])-
                        sum(dt.bigM.*molp.uij.*subp.α12) - sum(dt.bigM.*molp.ujk.*subp.α13) - sum(dt.bigM.*molp.ukl.*subp.α14))
                    push!(ocuts, (α5=subp.α5,α6=subp.α6,α7=subp.α7,α8=subp.α8,α9=subp.α9,α10=subp.α10,α11=subp.α11,α12=subp.α12,α13=subp.α13,α14=subp.α14))
                end
            else
                @info "FP == Feasibility cut found"
                @constraint( mp.m, 0 ≥ sum(subp.α5[5*(l-1)+p]*dt.d[l,p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
                    sum(dt.N["cap"][j]*mp.y[2*(j-1)+t]*subp.α7[2*(j-1)+t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*mp.y[2*(dt.N["plant"]+k-1)+t]*subp.α8[2*(k-1)+t] for k=1:dt.N["distribution"] for t=1:2)+
                    sum(dt.Vij[i]*mp.uij[i]*subp.α9[i] for i in findnz(dt.Vij)[1])+
                    sum(dt.Vjk[i]*mp.ujk[i]*subp.α10[i] for i in findnz(dt.Vjk)[1])-
                    sum(dt.bigM.*mp.uij.*subp.α12) - sum(dt.bigM.*mp.ujk.*subp.α13) - sum(dt.bigM.*mp.ukl.*subp.α14))
                @constraint( molp.m, 0 ≥ sum(subp.α5[5*(l-1)+p]*dt.d[l,p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
                    sum(dt.N["cap"][j]*molp.y[2*(j-1)+t]*subp.α7[2*(j-1)+t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*molp.y[2*(dt.N["plant"]+k-1)+t]*subp.α8[2*(k-1)+t] for k=1:dt.N["distribution"] for t=1:2)+
                    sum(dt.Vij[i]*molp.uij[i]*subp.α9[i] for i in findnz(dt.Vij)[1])+
                    sum(dt.Vjk[i]*molp.ujk[i]*subp.α10[i] for i in findnz(dt.Vjk)[1])-
                    sum(dt.bigM.*molp.uij.*subp.α12) - sum(dt.bigM.*molp.ujk.*subp.α13) - sum(dt.bigM.*molp.ukl.*subp.α14))
                push!(fcuts, (α5=subp.α5,α6=subp.α6,α7=subp.α7,α8=subp.α8,α9=subp.α9,α10=subp.α10,α11=subp.α11,α12=subp.α12,α13=subp.α13,α14=subp.α14))
            end
        else #no solution found by FP
            optimize!(mp.m)
            st = termination_status(mp.m)
            yb = value.(mp.y); ubij = value.(mp.uij); ubjk = value.(mp.ujk); ubkl = value.(mp.ukl); θb = value(mp.θ)
            subp = solve_dsp(dsp,yb,ubij,ubjk,ubkl)
            if subp.res == :OptimalityCut
                # @info "BD == Optimality cut found"
                if round(θb; digits=4) ≥ round(subp.obj; digits=4)
                    return
                else
                    @constraint(mp.m, mp.θ ≥ sum(subp.α5[5*(l-1)+p]*dt.d[l,p] for l=1:dt.N["customer"] for p=1:5)-
                        sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
                        sum(dt.N["cap"][j]*mp.y[2*(j-1)+t]*subp.α7[2*(j-1)+t] for j=1:dt.N["plant"] for t=1:2)-
                        sum(dt.N["cad"][k]*mp.y[2*(dt.N["plant"]+k-1)+t]*subp.α8[2*(k-1)+t] for k=1:dt.N["distribution"] for t=1:2)+
                        sum(dt.Vij[i]*mp.uij[i]*subp.α9[i] for i in findnz(dt.Vij)[1])+
                        sum(dt.Vjk[i]*mp.ujk[i]*subp.α10[i] for i in findnz(dt.Vjk)[1])-
                        sum(dt.bigM.*mp.uij.*subp.α12) - sum(dt.bigM.*mp.ujk.*subp.α13) - sum(dt.bigM.*mp.ukl.*subp.α14))
                    @constraint( molp.m, molp.θ ≥ sum(subp.α5[5*(l-1)+p]*dt.d[l,p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
                        sum(dt.N["cap"][j]*molp.y[2*(j-1)+t]*subp.α7[2*(j-1)+t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*mp.y[2*(dt.N["plant"]+k-1)+t]*subp.α8[2*(k-1)+t] for k=1:dt.N["distribution"] for t=1:2)+
                        sum(dt.Vij[i]*molp.uij[i]*subp.α9[i] for i in findnz(dt.Vij)[1])+
                        sum(dt.Vjk[i]*molp.ujk[i]*subp.α10[i] for i in findnz(dt.Vjk)[1])-
                        sum(dt.bigM.*molp.uij.*subp.α12) - sum(dt.bigM.*molp.ujk.*subp.α13) - sum(dt.bigM.*molp.ukl.*subp.α14))
                    push!(ocuts, (α5=subp.α5,α6=subp.α6,α7=subp.α7,α8=subp.α8,α9=subp.α9,α10=subp.α10,α11=subp.α11,α12=subp.α12,α13=subp.α13,α14=subp.α14))
                end
            else
                @info "BD == Feasibility cut found"
                @constraint(mp.m, 0  ≥ sum(subp.α5[5*(l-1)+p]*dt.d[l,p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
                    sum(dt.N["cap"][j]*mp.y[2*(j-1)+t]*subp.α7[2*(j-1)+t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*mp.y[2*(dt.N["plant"]+k-1)+t]*subp.α8[2*(k-1)+t] for k=1:dt.N["distribution"] for t=1:2)+
                    sum(dt.Vij[i]*mp.uij[i]*subp.α9[i] for i in findnz(dt.Vij)[1])+
                    sum(dt.Vjk[i]*mp.ujk[i]*subp.α10[i] for i in findnz(dt.Vjk)[1])-
                    sum(dt.bigM.*mp.uij.*subp.α12) - sum(dt.bigM.*mp.ujk.*subp.α13) - sum(dt.bigM.*mp.ukl.*subp.α14))
                @constraint( molp.m, 0 ≥ sum(subp.α5[5*(l-1)+p]*dt.d[l,p] for l=1:dt.N["customer"] for p=1:5)-sum(dt.N["cas"][i]*subp.α6[i] for i=1:dt.N["supplier"])-
                    sum(dt.N["cap"][j]*molp.y[2*(j-1)+t]*subp.α7[2*(j-1)+t] for j=1:dt.N["plant"] for t=1:2)-sum(dt.N["cad"][k]*molp.y[2*(dt.N["plant"]+k-1)+t]*subp.α8[2*(k-1)+t] for k=1:dt.N["distribution"] for t=1:2)+
                    sum(dt.Vij[i]*molp.uij[i]*subp.α9[i] for i in findnz(dt.Vij)[1])+
                    sum(dt.Vjk[i]*molp.ujk[i]*subp.α10[i] for i in findnz(dt.Vjk)[1])-
                    sum(dt.bigM.*molp.uij.*subp.α12) - sum(dt.bigM.*molp.ujk.*subp.α13) - sum(dt.bigM.*molp.ukl.*subp.α14))
                push!(fcuts, (α5=subp.α5,α6=subp.α6,α7=subp.α7,α8=subp.α8,α9=subp.α9,α10=subp.α10,α11=subp.α11,α12=subp.α12,α13=subp.α13,α14=subp.α14))
            end
        end
    end
    return (mp,dsp,sol,ocuts,fcuts,Archiv)
end
molp = MOLP(w); AddCuts(molp,optcuts,unique(feasicuts))
mp = MasterP(w); AddCuts(mp,optcuts,unique(feasicuts)); dsp = DualSP(w); dis = buildMP(); fbmodel = buildMP();
ocuts = []; fcuts = []; Archiv = []; Y = []; U1 = []; U2= []; U3 = []; feasisol = [];
mp2,dsp2,sol,ocuts,fcuts,archiv = benders_decomposition()
1
length(unique(Archiv))
length(ocuts),length(unique(fcuts))

feasisol[1]


mp = MasterP(w);
AddCuts(mp,ocuts,unique(fcuts))
optimize!(mp.m); termination_status(mp.m)

objective_value(mp.m)
sum(value.(mp.y))
sum(value.(mp.uij))
sum(value.(mp.ujk))
sum(value.(mp.ukl))
print(mp.m)
