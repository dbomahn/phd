#######################Lexicographic####################
using JuMP,GLPK,LinearAlgebra

###################### coefficients  ###################
n = 3
m = 4
C_1 =[2 5 4 7; 3 3 5 7; 3 8 4 2; 6 5 2 5]
C_2 = [3 3 6 2; 5 3 7 3; 5 2 7 4; 4 6 3 5]
C_3 = [4 2 5 3; 5 3 4 3; 4 3 5 2; 6 4 7 3]

###################### Lexicographic Optimisation  ###################
lex1 = Model( with_optimizer(GLPK.Optimizer))
@variable(lex1, x1[1:m,1:m], Bin)
@expression(lex1, ex1, sum(C_1[i,j]*x1[i,j] for i=1:m, j=1:m) )
@expression(lex1, ex2, sum(C_2[i,j]*x1[i,j] for i=1:m, j=1:m) )
@expression(lex1, ex3, sum(C_3[i,j]*x1[i,j] for i=1:m, j=1:m) )
@constraint( lex1, cols[i=1:n], sum(x1[i,j] for j=1:m)== 1)
@constraint( lex1, rows[j=1:m], sum(x1[i,j] for i=1:m)== 1 )
@objective(lex1, Min ,ex1)#+ex2+ex3)
# print(mm)
# termination_status(mm)
optimize!(lex1)
lex1x = value.(x1)
lex1y = [objective_value(lex1),dot(lex1x,C_2),dot(lex1x,C_3)]

lex2 = Model( with_optimizer(GLPK.Optimizer))
@variable(lex2, x2[1:m,1:m], Bin)
@expression(lex2, ex1, sum(C_1[i,j]*x2[i,j] for i=1:m, j=1:m) )
@expression(lex2, ex2, sum(C_2[i,j]*x2[i,j] for i=1:m, j=1:m) )
@expression(lex2, ex3, sum(C_3[i,j]*x2[i,j] for i=1:m, j=1:m) )
@constraint( lex2, cols[i=1:m], sum(x2[i,j] for j=1:m)== 1)
@constraint( lex2, rows[j=1:m], sum(x2[i,j] for i=1:m)== 1 )
@objective(lex2, Min ,ex2)
optimize!(lex2)
lex2x = value.(x2)
lex2y = [dot(lex2x,C_1),objective_value(lex2),dot(lex2x,C_3)]

lex3 = Model( with_optimizer(GLPK.Optimizer))
@variable(lex3, x3[1:m,1:m], Bin)
@expression(lex3, ex1, sum(C_1[i,j]*x3[i,j] for i=1:m, j=1:m) )
@expression(lex3, ex2, sum(C_2[i,j]*x3[i,j] for i=1:m, j=1:m) )
@expression(lex3, ex3, sum(C_3[i,j]*x3[i,j] for i=1:m, j=1:m) )
@constraint( lex3, cols[i=1:m], sum(x3[i,j] for j=1:m)== 1)
@constraint( lex3, rows[j=1:m], sum(x3[i,j] for i=1:m)== 1 )
@objective(lex3, Min ,ex3)
optimize!(lex3)
lex3x = value.(x3)
lex3y = [dot(lex3x,C_1),dot(lex3x,C_2),objective_value(lex3)]

# Store lexi solutions
X= []
Y= []
push!(X,lex1x)
push!(X,lex2x)
push!(X,lex3x)
push!(Y,lex1y)
push!(Y,lex2y)
push!(Y,lex3y)

# YN = objective_value()
# S = Dict()
# XE = Dict()
# Xs = []
# for i=1:length(YN)
#     for j=1:n
#         push!(Xs, value(x[i,j]))
# end
# Xs = unique!(Xs)
# YN = unique!(YN)
# for i = 1:length(YN)
#     S[i] = YN[i]
#     XE[i] = Xs[i]
# end
