include("common.jl")

function example(degree::Int, n::Int, nq::Int=4)
    msh=Mesh(Rectangle(0,2,0,1),n,n)

    f(x,y)=1
    p(x,y)=4
    q1(x,y)=4
    q2(x,y)=2
    r(x,y)=2
    Uboundary=[(x,y)->exp(x)+1/2,(x,y)->exp(y+2)+1/2,
        (x,y)->exp(x+1)+1/2,(x,y)->exp(y)+1/2]
    bc=[BoundaryCondition(Dirichlet,Uboundary[i]) for i=1:4]
    bvp=BoundaryValueProblem(p,[q1,q2],r,f,bc)

    U,a=Solve(bvp,msh,degree,nq)

    u(x,y)=exp(x+y)+1/2
    print("norm(u-U)="); display(norm(u.(a.msh.x,a.msh.y)-U,Inf))

    #I,J,V=findnz(A)
    #write_matfile("data/poisson.mat",amsh=a.msh,msh=Refine(msh,degree),U=U,b=a.b,I=I,J=J,V=V)
end

example(2,10)
