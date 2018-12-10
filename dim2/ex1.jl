include("common.jl")

function example(degree::Int, n::Int, nq::Int=4)
    msh=Mesh(Rectangle(0,1,0,1),n,n)

    f(x,y)=32*(y*(1-y)+x*(1-x))
    bc=[BoundaryCondition(Dirichlet,f0) for i=1:4]
    global bvp=Poisson(f,bc)

    U,a,A=Solve(bvp,msh,degree,nq)

    u(x,y)=16*x*(1-x)*y*(1-y)
    print("norm(u-U)="); display(norm(u.(a.msh.x,a.msh.y)-U,Inf))

    I,J,V=findnz(A)
    write_matfile("data/poisson.mat",amsh=a.msh,msh=Refine(msh,degree),U=U,b=a.b,I=I,J=J,V=V)
end

example(2,10)
