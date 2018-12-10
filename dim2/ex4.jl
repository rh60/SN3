include("common.jl")

function example(degree::Int, n::Int, nq::Int=4)
    msh=Mesh(Rectangle(-1,1,-1,1),n,n)

    f(x,y)=-2*exp(x)*(1+x+6*y)
    p(x,y)=exp(x)
    q1(x,y)=3*x
    q2(x,y)=2*y
    r(x,y)=-6

    bc=[BoundaryCondition(Dirichlet,(x,y)->x^2-2)
        BoundaryCondition(Dirichlet,(x,y)->2*y^3+1)
        BoundaryCondition(Neumann,(x,y)->6*exp(x))
        BoundaryCondition(Neumann,(x,y)->2*exp(x))]
    bvp=BoundaryValueProblem(p,[q1,q2],r,f,bc)

    U,a,A=Solve(bvp,msh,degree,nq)

    u(x,y)=x^2+2*y^3
    print("norm(u-U)="); display(norm(u.(a.msh.x,a.msh.y)-U,Inf))

    I,J,V=findnz(A)
    write_matfile("data/ex4.mat",a=a,msh=Refine(msh,degree),U=U,I=I,J=J,V=V)
end

example(4,10,6)
