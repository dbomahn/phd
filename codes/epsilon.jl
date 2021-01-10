using Revise,DelimitedFiles,JuMP,GLPK,GLPKMathProgInterface, vOptGeneric, LinearAlgebra
# competible JuMP version with vOptGeneric: 0.18.6 
# using MathProgBase
path = "C:/Users/AK121396/Desktop/Dropbox/JKU/Julia/KP/"
# path = "/home/ak121396/Dropbox/JKU/Julia/KP/"
files = readdir(path)
# f =readdlm(path*"KP_p-3_n-30_ins-"*".dat", '\t', String, '\n')
fl = readdlm(path*files[1], '\t', String, '\n')
#################################################################
obj=parse(Int,fl[1])
n=parse(Int,fl[2])
ub=parse(Int,fl[3])
########################    KP     #####################
b = fl[4:length(fl)-1]
C= ones(length(b),n)
C = round.(Int,C)
global ct=0;
for x=1:length(b)
    a = split(b[x],('[',']',','))
    aa = filter!(e->!(e in [ "" ,"[","," ,"]"]) ,a)
    for y=1:length(aa)
        c=parse(Int64,aa[y])
        C[x,y] = c
        if c==0
            global ct = ct+1;
        end
    end
end
######################## coefficient matrix (B) #############################
weight =ones(1,n)
weight = round.(Int,weight)
item = fl[length(fl)]
w1 = split(item, ('[',']',','))
w2 = filter!(e->!(e in ["" ,"[", "]"]) ,w1)
for i=1:n
    weight[i] = parse(Int64,w2[i])
end
###############################  FUNCTIONS  ##############################
function getConstraints(m,e,E)
    ϵ = Dict()
    ϵp = Dict()
    k = mod(m,length(E)+1)
    d = div(m,length(E)+1)
    ϵ[2] = e[2][k+1]
    ϵp[2] = e[2][k+2]
    ϵ[3] = e[3][d+1]
    ϵp[3] = e[3][d+2]
    return ϵ,ϵp
end
lex = vModel(solver = GLPKSolverMIP())
@variable(lex, x[1:n], Bin)
@addobjective(lex, Min, -dot(C[1,:],x[:]) )
@addobjective(lex, Min, -dot(C[2,:],x[:]) )
@addobjective(lex, Min, -dot(C[3,:],x[:]) )
@constraint(lex, -sum(weight[i]*x[i] for i=1:n) >= -ub )
@constraint(lex, ϵ[2] <= -dot(C[2,:],x[:]) <= ϵp[2])
@constraint(lex, ϵ[3] <= -dot(C[3,:],x[:]) <= ϵp[3])
status = solve(lex,method=:lex)
Y_lex = unique(getY_N(lex))

function opt(ϵ,ϵp)
    #1st phase
    lex = vModel(solver = GLPKSolverMIP())
    @variable(lex, x[1:n], Bin)
    @addobjective(lex, Min, -dot(C[1,:],x[:]) )
    @constraint(lex, -sum(weight[i]*x[i] for i=1:n) >= -ub )
    @constraint(lex, ϵp[2] <= -dot(C[2,:],x[:]) <= ϵ[2])
    @constraint(lex, ϵp[3] <= -dot(C[3,:],x[:]) <= ϵ[3])
    status = solve(lex,method=:lex)
    Y_lex = unique(getY_N(lex))
    f = [Y_lex[1],Y_lex[2],Y_lex[3]]
    #2nd phase
    m2 = Model(solver=GLPKSolverMIP(presolve=true))
    @addobjective(lex, Min, -dot(C[2,:],x[:]) )
    @addobjective(lex, Min, -dot(C[3,:],x[:]) )
    @variable( m2, 0<=x2[1:n]<=1, Int )
    @constraint( m2, -dot(C[1,:],x2[:]) >= f[1])
    @constraint(lex, ϵ[2] <= -dot(C[2,:],x[:]) <= ϵp[2])
    @constraint(lex, ϵ[3] <= -dot(C[3,:],x[:]) <= ϵp[3])
    @objective( m2, Min, -dot(C[2,:],x2[:]) - (dot(C[3,:],x2[:])/(Y_lex[3][3])*3) )
    status = solve(m2)
    f2 = -dot(C[2,:],JuMP.getvalue(x2))
    f3 = -dot(C[3,:],JuMP.getvalue(x2))
    Y = [f[1],f2,f3]
    if status == :Optimal
        return JuMP.getvalue(x2),Y
    else
        return NaN,NaN
    end
end

function updateConstraints(Y,e)
    for j=2:3
        push!(e[j],Y[j])
        sort(e[j])
    end
end

function mindominate(x,y)
    all(i -> x[i] <= y[i], eachindex(x)) && any(i -> x[i] < y[i], eachindex(x))
end
function maxdominate(x,y)
    all(i -> x[i] >= y[i], eachindex(x)) && any(i -> x[i] > y[i], eachindex(x))
end
################################################################################

E = []
F= Set()
e = Dict(2=>[-Inf,1], 3=>[-Inf,0])
copye = Dict(2=>[-Inf,1], 3=>[-Inf,0])
m = 0
# c =1
c = (length(E)+1)^2

for i=1:20
    global E,F,e,m,c
    for m = 0:((length(E)+1)^2)-1
        if m>=c
            return E
        end
        global ϵ,ϵp = getConstraints(m,e,E)
        if !(Set([[ϵ[2],ϵp[2]],[ϵ[3],ϵp[3]]]) in F) && length(E)>0
            global s,Y = opt(ϵ,ϵp)
            if Y == NaN || sum([mindominate(Y,E[i]) for i=1:length(E)]) == length(E)-1 #no feasible or dominated sol
                union( F,Set([[ϵ[2],ϵp[2]],[ϵ[3],ϵp[3]]]) )
            else
                break
            end
        end
    end
    if length(Y)>2
        push!(E,Y)
        F = union(F,Set([[Y[2],ϵp[2]],[Y[3],ϵp[3]]]))
        updateConstraints(Y,e)
    end
    c = (length(E)+1)^2
    m = 0
end






#     ϵ_ϵp = [[0,0],[0,0]] #Dict(2=>[0,0],3=>[0,0])
#     # m = ((length(E)+1)^2)-1
#     for j=2:3
#         d = mod(m,length(E)+1)
#         # m = (m-d)/(length(E)+1)
#         ϵ_ϵp[j][1] = e[j][d]
#         ϵ_ϵp[j][2] = e[j][d+1]
#     end
#     return ϵ_ϵp #values(ϵ_ϵp)
# end


##########################  Lexicographic Optimisation  #######################
# lex = vModel(solver = GLPKSolverMIP())
# @variable(lex, x[1:n], Bin)
# @addobjective(lex, Min, -dot(C[1,:],x[:]) )
# @addobjective(lex, Min, -dot(C[2,:],x[:]) )
# @addobjective(lex, Min, -dot(C[3,:],x[:]) )
# @constraint(lex, sum(weight[i]*x[i] for i=1:n) <= ub )
# status = solve(lex,method=:lex)
# Y_lex = unique(getY_N(lex))
# ini_f = [Y_lex[1][1],Y_lex[2][2],Y_lex[3][3]]
# X_lex = []
# for i=1:length(Y_lex)
#     xlex = vOptGeneric.getvalue(x, i)
#     push!(X_lex,xlex)
# end
