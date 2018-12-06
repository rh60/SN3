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
    q1=a.bvp.q[1].(qx,qy)
    q2=a.bvp.q[2].(qx,qy)
    #nq=length(qx)
    r=a.bvp.r.(qx,qy)
    f=a.bvp.f.(qx,qy)
    for i=1:l
        for j=1:l
            push!(a.I,loc2glob[i]);
            push!(a.J,loc2glob[j])
            gi=invJT*[a.L1[i] a.L2[i]]'
            gj=invJT*[a.L1[j] a.L2[j]]'
            g=gj[1,:].*gi[1,:]+gj[2,:].*gi[2,:]
            #push!(a.V,intSum(a.Q,p.*g)*J)
            tmp=intSum(a.Q,p.*g+(q1.*gj[1,:]+q2.*gj[2,:]+r.*a.L[j]).*a.L[i])
            push!(a.V,tmp)
        end
        #a.b[t[i]] += intSum(a.Q,f.*a.L[i])*J
        a.b[loc2glob[i]] += intSum(a.Q,f.*a.L[i])
    end
end

function FEM(degree::Int, n::Int, nq::Int=4)
    msh=Mesh(Rectangle(0,2,0,1),n,n)
    E=Lagrange(degree)
    f(x,y)=1
    p(x,y)=4
    q1(x,y)=4
    q2(x,y)=2
    r(x,y)=2
    Uboundary=[(x,y)->exp(x)+1/2,(x,y)->exp(y+2)+1/2,
        (x,y)->exp(x+1)+1/2,(x,y)->exp(y)+1/2]
    bc=[BoundaryCondition(Dirichlet,Uboundary[i]) for i=1:4]
    bvp=BoundaryValueProblem(p,[q1,q2],r,f,bc)

    a=Assembly(LagrangeMesh(msh,E.degree),bvp,E,nq)
    nt=size(a.msh.tri,1)
    for k=1:nt
        CalculateContributionOnTriangle!(a,k)
    end
    A=sparse(a.I,a.J,a.V)

    for (ib,bdr) in enumerate(a.msh.boundary)
        inodes = bdr[:]
        A[inodes,:].=0
        for i in inodes
            A[i,i]=1
        end
        x=a.msh.x[inodes]
        y=a.msh.y[inodes]
        a.b[inodes]=bc[ib].value.(x,y)
    end

    U=A\a.b
    u(x,y)=exp(x+y)+1/2
    print("norm(u-U)="); display(norm(u.(a.msh.x,a.msh.y)-U,Inf))

    I,J,V=findnz(A)
    write_matfile("data/poisson.mat",amsh=a.msh,msh=Refine(msh,degree),U=U,b=a.b,I=I,J=J,V=V)
end

FEM(2,10)
