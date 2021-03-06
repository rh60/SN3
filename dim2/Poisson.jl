include("common.jl")

struct Assembly
    msh::Mesh
    bvp::BoundaryValueProblem
    L::Vector{Vector{Float64}}
    L1::Vector{Vector{Float64}}
    L2::Vector{Vector{Float64}}
    Q::QuadTri
    I::Vector{Int}
    J::Vector{Int}
    V::Vector{Float64}
    b::Vector{Float64}
    function Assembly(msh::Mesh,bvp::BoundaryValueProblem,E::Lagrange,nq::Int)
        L,L1,L2,Q=Quadrature(E,nq)
        b=zeros(length(msh.x))
        new(msh,bvp,L,L1,L2,Q,Int[],Int[],Float64[],b)
    end
end

@inline function CalculateContributionOnTriangle!(a::Assembly,k::Int)
    loc2glob=a.msh.tri[k,:]
    t=loc2glob[1:3]
    x=a.msh.x[t];y=a.msh.y[t];
    qx=affinity(x,a.Q.Points.x,a.Q.Points.y)
    qy=affinity(y,a.Q.Points.x,a.Q.Points.y)
    J,invJT=LinearTranform(x,y)
    l=length(a.L)
    p=a.bvp.p.(qx,qy)
    f=a.bvp.f.(qx,qy)
    for i=1:l
        for j=1:l
            push!(a.I,loc2glob[i]);
            push!(a.J,loc2glob[j])
            gi=invJT*[a.L1[i] a.L2[i]]'
            gj=invJT*[a.L1[j] a.L2[j]]'
            g=gi[1,:].*gj[1,:]+gi[2,:].*gj[2,:]
            #push!(a.V,intSum(a.Q,p.*g)*J)
            push!(a.V,intSum(a.Q,p.*g))
        end
        #a.b[t[i]] += intSum(a.Q,f.*a.L[i])*J
        a.b[loc2glob[i]] += intSum(a.Q,f.*a.L[i])
    end
end

function FEM(degree::Int, n::Int, nq::Int=4)
    msh=Mesh(Rectangle(0,1,0,1),n,n)

    E=Lagrange(degree)
    f(x,y)=32*(y*(1-y)+x*(1-x))
    vbc=[BoundaryCondition(Dirichlet,f0) for i=1:4]
    bvp=BoundaryValueProblem(f,vbc)

    a=Assembly(LagrangeMesh(msh,E.degree),bvp,E,nq)
    nt=size(a.msh.tri,1)
    for k=1:nt
        CalculateContributionOnTriangle!(a,k)
    end
    A=sparse(a.I,a.J,a.V)

    dirt=Int[]
    for ib=1:4
        dirt = [dirt; a.msh.boundary[ib][:,1]; a.msh.boundary[ib][:,3:end][:] ]
    end
    A[dirt,:].=0
    for i in dirt
        A[i,i]=1
    end
    a.b[dirt].=0

    U=A\a.b
    u(x,y)=16*x*(1-x)*y*(1-y)
    print("norm(u-U)="); display(norm(u.(a.msh.x,a.msh.y)-U,Inf))

    I,J,V=findnz(A)
    write_matfile("data/poisson.mat",amsh=a.msh,msh=Refine(msh,degree),U=U,b=a.b,I=I,J=J,V=V)
end
