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
    a=atan(y,x)
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
    rmsh=Refine(msh,degree)
    write_matfile("data/ex5.mat",a=a,msh=rmsh,U=U,I=I,J=J,V=V)
    write_vtkfile("data/ex5","U",rmsh.x,rmsh.y,rmsh.tri,U)
end

example(3,4)
@async run(`"C:/Program Files/ParaView 5.6.0-Windows-msvc2015-64bit/bin/paraview.exe" data/ex5.pvsm`)
