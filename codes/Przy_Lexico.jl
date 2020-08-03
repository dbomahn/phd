using  vOptGeneric,JuMP, CPLEX, LinearAlgebra
 # MultiJuMP, JuMP,

n = 4
C = [[2 5 4 7; 3 3 5 7; 3 8 4 2; 6 5 2 5],[3 3 6 2; 5 3 7 3; 5 2 7 4; 4 6 3 5],[4 2 5 3; 5 3 4 3; 4 3 5 2; 6 4 7 3]]

################### LEXICOGRAPHIC OPTIMAL SOLUTIONS  ###################
# vOpt Model
lex = vModel(solver = CplexSolver())
@variable(lex, x[1:n,1:n]>=0)
@addobjective(lex, Min, sum(C[1][i,j]*x[i,j] for i=1:n, j=1:n) )
@addobjective(lex, Min, sum(C[2][i,j]*x[i,j] for i=1:n, j=1:n) )
@addobjective(lex, Min, sum(C[3][i,j]*x[i,j] for i=1:n, j=1:n) )
@constraint( lex, cols[i=1:n], sum(x[i,j] for j=1:n)== 1)
@constraint( lex, rows[j=1:n], sum(x[i,j] for i=1:n)== 1 )
# print(lex)
vOptGeneric.solve(lex, method = :lex)
# getobjectivevalue(lex)
# XE = Dict()
# Xs = []
global Y = getY_N(lex)

# for i=1:length(Y)
#     push!(Xs, vOptGeneric.getvalue(x, i))
# end
#
# for i=1:length(Y)
#     for j=1:n*n
#         Xs[i][j]=abs(round(Xs[i][j])) #,RoundNearestTiesAway
#     end
# end

for i=1:length(Y)
    for j=1:3
        Y[i][j]=round(Y[i][j],RoundNearestTiesAway) #Rounding might not be needed
    end
end


# Xs = unique!(Xs)
global Y = unique!(Y)
global S = copy(Y)
global A = Dict()
global P = Dict()
for i=1:length(Y)
    if (!(haskey(A, i)))
        A[i] = []
        P[i] = []
    end
end

function lamcost(i,j)
    bim1 = vModel(solver=CplexSolver())
    @variable(bim1, lam1[1:3]>=0)
    @expression( bim1, e1, dot(S[i],lam1) )
    @expression( bim1, e2, dot(S[j],lam1) )
    @constraint( bim1, e1 <= e2 )
    @constraint( bim1, sum{lam1[k],k=1:3} ==1 )
    @addobjective(bim1, Max, e1-e2)
    @addobjective(bim1, Min, e2-e1)
    vOptGeneric.solve(bim1, method=:dicho)
    ld1 = vOptGeneric.getvalue(lam1,1) # for k=1:length(getY_N(bim1))

    bim2 = vModel(solver=CplexSolver())
    @variable(bim2, lam2[1:3]>=0)
    @expression( bim2, ex1, dot(S[j],lam2) )
    @expression( bim2, ex2, dot(S[i],lam2) )
    @constraint( bim2, ex1 <= ex2)
    @constraint( bim2, sum{lam2[k],k=1:3} ==1 )
    @addobjective(bim2, Max, ex1-ex2)
    @addobjective(bim2, Min, ex2-ex1)
    vOptGeneric.solve(bim2, method=:dicho)
    ld2 = vOptGeneric.getvalue(lam2,1) # for k=1:length(getY_N(bim2))


    Lambda = Dict()
    c = Dict()
    Lambda[i] = [ld1,ld2]
    for y=1:length(Lambda[i])
        c[y] = sum(C[k]*Lambda[i][y][k] for k=1:3)
    end
    return Lambda,c
end

lamcost(1,2)
