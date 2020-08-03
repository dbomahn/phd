
for i = 1 #:length(Y)-1
    # P[i] = []
    P[i] = collect(i+1:length(S))
    P[i] = filter(x->!(x in A[i]), P[i])
    P[i] = unique(P[i])
    A[i] = []
    # PS = append!(P[i],St)
    for j=3# in P
        # if !(j in A[i]) # if !(j in keys(A[i])) #skip adjacent elements
            global ld,c=lamcost(i,j)
            # if ld[i][1]!=ld[i][2]
                # SOLVE bi-objective problem
                bi = vModel(solver = CplexSolver())
                @variable(bi, xx[1:n,1:n]>=0 )
                @addobjective(bi, Min, sum(c[1][i,j]*xx[i,j] for i=1:n, j=1:n) )
                @addobjective(bi, Min, sum(c[2][i,j]*xx[i,j] for i=1:n, j=1:n) )
                @constraint( bi, cols[i=1:n], sum(xx[i,j] for j=1:n)== 1)
                @constraint( bi, rows[j=1:n], sum(xx[i,j] for i=1:n)== 1 )
                vOptGeneric.solve(bi, method = :dicho)

                for l=1:length(getY_N(bi))
                    newx = vOptGeneric.getvalue(xx,l) #round(Integer, getvalue(x))
                    a = round(dot(newx,C[1]),RoundNearestTiesAway)
                    b = round(dot(newx,C[2]),RoundNearestTiesAway)
                    c = round(dot(newx,C[3]),RoundNearestTiesAway)
                    @show edge = [a,b,c]
                    if !(edge in S) #y1--y3 (1)(3) case
                        push!(S,edge)
                        push!(P[i], length(S))
                        # push!(P[j], length(S)+1)
                        A[length(S)+1]=[]
                        P[length(S)+1]=[]

                    else # y1--y2 or y1--y3(2) case
                        if [round(dot(ld[i][z],Y[j]),RoundNearestTiesAway) == round(dot(ld[i][z],Y[i]),RoundNearestTiesAway) for z=1:length(ld[i])] == [true for x=1:length(ld[i])]
                            push!(A[i], j)
                            push!(A[j], i)
                            P[i] = filter(x->x!=j,P[i])
                        end
                    end

                end
                #calculate <lambda,y> to figure out nonextreme point
                # if [round(dot(ld[i][z],Y[j]),RoundNearestTiesAway) == round(dot(ld[i][z],Y[i]),RoundNearestTiesAway) for z=1:length(ld[i])] == [true for x=1:length(ld[i])]
                #     push!(A[i], j)
                #     push!(A[j], i)
                #     P[i] = filter(x->x!=j,P[i])
                # end
    end
    # end
end
