@enum  BCenum Dirichlet=0 Neumann=1

struct BoundaryCondition
    code::BCenum
    value::Function
end

struct BoundaryValueProblem
    p::Function # difuzní koeficient
    q::Vector{Function} # konvekce
    r::Function # reakce
    f::Function # pravá strana
    vbc::Vector{BoundaryCondition} # okrajové podmínky
end

f1(x::Float64,y::Float64)=1.0
f0(x::Float64,y::Float64)=0.0

function Poisson(f::Function,vbc::Vector{BoundaryCondition})
    BoundaryValueProblem(f1,[f0, f0],f0,f,vbc)
end

function Laplace(vbc::Vector{BoundaryCondition})
    BoundaryValueProblem(f0,[f0, f0],f0,f,vbc)
end
