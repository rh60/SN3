include("common.jl")

function loadmsh(fname::String)::Mesh
    mat=read_matfile(fname)
    val=jvalue(mat["msh"])
    x=val["POS"][:,1]
    y=val["POS"][:,2]
    tri=Int.(val["TRIANGLES"][:,1:3])
    lines=Int.(val["LINES"])
    bdry=[lines[:,1:2]]
    msh=Mesh(x,y,tri,bdry)
end

function u(x::Float64,y::Float64)::Float64
    r=sqrt(x^2+y^2)
    a=atan(x/y)
    if r==0
        return 0.0
    else
        return r^(2/3)*sin(2*a/3+π/3)
    end
end

function example(degree::Int, nq::Int=4)
    msh=loadmsh("geo/L.mat")

    bc=[BoundaryCondition(Dirichlet,u)]
    bvp=Laplace(bc)

    U,a,A=Solve(bvp,msh,degree,nq)

    print("norm(u-U)="); display(norm(u.(a.msh.x,a.msh.y)-U,Inf))

    I,J,V=findnz(A)
    write_matfile("data/ex5.mat",a=a,msh=Refine(msh,degree),U=U,I=I,J=J,V=V)
end

example(4,8)
