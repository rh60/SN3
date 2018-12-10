include("common.jl")

function example(degree::Int, n::Int, nq::Int=4)
    msh=Mesh(Rectangle(0,2,-1,0),n,n)

    f(x,y)=2*x*exp(y)-2
    p(x,y)=exp(-x^2-y)
    q1(x,y)=exp(-x^2)
    q2(x,y)=exp(-y)
    r(x,y)=-exp(-y)

    bc=[BoundaryCondition(Neumann,(x,y)->-1)
        BoundaryCondition(Dirichlet,(x,y)->exp(y+4))
        BoundaryCondition(Dirichlet,(x,y)->exp(x^2))
        BoundaryCondition(Neumann,(x,y)->0)]
    bvp=BoundaryValueProblem(p,[q1,q2],r,f,bc)

    U,a=Solve(bvp,msh,degree,nq)

    u(x,y)=exp(x^2+y)
    print("norm(u-U)="); display(norm(u.(a.msh.x,a.msh.y)-U,Inf))

    #I,J,V=findnz(A)
    #write_matfile("data/poisson.mat",amsh=a.msh,msh=Refine(msh,degree),U=U,b=a.b,I=I,J=J,V=V)
end

example(4,10,6)
