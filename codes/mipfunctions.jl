
function mipgetlocation(x,fmin,steps)
    L = [dot(coef[1,:],x),dot(coef[2,:],x),dot(coef[3,:],x)]
    loca = [round.(Int,((L[k]-fmin[k])/steps[k])+1) for k=1:3]
    return loca
end
function fbsearch(idx,steps,fmin,txi,fmax) #solveLP
    yub = [idx[j]*steps[j]+fmin[j] for j=1:3]
    ylb = [(idx-[1,1,1])[j]*steps[j]+fmin[j] for j=1:3]
    idx0 = findall(x->x==0, txi)
    idx1 = findall(x->x==1, txi)

    ##################### Build JuMP model #################
    jump_model = Model()
    mof_model = MathOptFormat.MPS.Model()
    MOI.read_from_file(mof_model, path1*"/instances/neos-1620770.mps")
    MOI.copy_to(backend(jump_model), mof_model)

    MOI.set(jump_model, MOI.RawParameter("CPX_PARAM_SCRIND"), false)


    @constraint( lp, c1,  ylb[2]+1 <=dot(coef[1,:], x) <= yub[2]-1 ) #  fmax[2])#
    @constraint( lp, c2,  ylb[1]+1 <=dot(coef[2,:], ) <= yub[1]-1 ) #  fmax[1])#
    @constraint( lp, c3,  ylb[3]+1 <=dot(coef[3,:], ) <= yub[3]-1 ) #fmax[3]) #
    optimize!(lp)

    if termination_status(jump_model) == MOI.OPTIMAL
        return JuMP.value.(x)
    else
        return 0; #print("no new lp solution found");
    end
end



function mipFP(Xf,PFset,candX,LPcount,steps,fmin,fmax)
    Tabu = [];
    for k=1:length(candX)
        print("===============Feasi Pump",k," th candidate sol ==============\n")

        x_t = candX[k];
        SearchDone = false
        itr = 1
        Max_itr = j+i #max(count(x->0<x<1,x_t),1)  #maximum number of attempts => How to set
        while itr<Max_itr && SearchDone==false
            xi_t = round.(x_t)
            fxi = [dot(coef[1,:],xi_t),dot(coef[2,:],xi_t),dot(coef[3,:],xi_t)]
            if ( (fbcheck(xi_t) == true) && (dominated(fxi,PFset)==false) ) #checking easibility and dominance   #z[1]<=fmax[1]-1 && z[2]<=fmax[2]-1 && z[3]<=fmax[3]-1)
                idx = getlocation(xi_t,fmin,steps)
                push!(PFset,fxi) #add new solval to PFset
                if !haskey(Xf,idx)
                    # print("ROUNING => another solution added \n")
                    Xf[idx] = xi_t
                else
                    # print("ROUNING =>new intsol to a cuboid \n")
                    push!([Xf[idx]], xi_t)
                end
                SearchDone = true
            else
                if xi_t ∈ Tabu
                    xi_t = flipoper(Tabu,x_t,xi_t)
                    if xi_t==[]
                        # print("FLIP didn't work \n")
                        SearchDone = true
                    else
                        fxi = [dot(coef[1,:],xi_t),dot(coef[2,:],xi_t),dot(coef[3,:],xi_t)]
                        if ( (fbcheck(xi_t) == true) && (dominated(fxi,PFset)==false) )  #z[1]<=fmax[1]-1 && z[2]<=fmax[2]-1 && z[3]<=fmax[3]-1 )
                            idx = getlocation(xi_t,fmin,steps)
                            push!(PFset,fxi) #add new solval to PFset
                            if !haskey(Xf,idx)
                                # print("FLIP=> another solution's added \n")
                                Xf[idx] = xi_t
                            else
                                # print("FLIP=> new intsol to a cuboid \n")
                                push!([Xf[idx]], xi_t)
                            end
                            SearchDone = true
                        end
                    end
                end
                if SearchDone == false
                    push!(Tabu,xi_t) #when break
                    idx = mipgetlocation(xi_t,fmin,steps)
                    x_t = fbsearch(idx,steps,fmin,xi_t,fmax)
                    LPcount+=1
                    if x_t == 0 #when there's no new feasible lp sol
                        # print("no lp sol's found");
                        break
                    # else
                        # print("\n New lp obj: ",giveobjval(x_t), "\n")
                    end
                end
            end
            itr+=1; #@show itr
        end
    end
    # end
    return Xf,PFset,Tabu,LPcount
end
