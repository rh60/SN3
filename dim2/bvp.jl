struct BoundaryCondition
    code::Int # 0 Dirichlet 1 Newton
    value::Function
end

struct BoundaryValueProblem
    p::Function # difuzní koeficient
    q::Function # konvekce
    r::Function # reakce
    f::Function # pravá strana
    xLeft::BoundaryCondition # okrajová podmínka vlevo
    xRight::BoundaryCondition # vpravo
    yLeft::BoundaryCondition # okrajová podmínka dole
    yRight::BoundaryCondition # nahoře
end
