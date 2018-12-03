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
    t=a.msh.tri[k,1:3]
    x=a.msh.x[t];y=a.msh.y[t];
    qx=affinity(x,a.Q.Points.x)
    qy=affinity(x,a.Q.Points.y)
    J,invJT=LinearTranform(x,y)
    l=length(a.L)
    p=a.bvp.p.(qx,qy)
    for i=1:l
        for j=1:l
            push!(a.I,t[i]); push!(a.J,t[j])
            gi=invJT*[a.L1[i] a.L2[i]]'
            gj=invJT*[a.L1[j] a.L2[j]]'
            g=gi[1,:].*gj[1,:]+gi[2,:].*gj[2,:]
            push!(a.V,intSum(a.Q,p.*g)*J)
        end
        a.b[t[i]] += intSum(a.Q,a.bvp.f.(qx,qy).*a.L[i])*J
    end
end

#function FEM
msh=Mesh(Rectangle(0,1,0,1),10,10)
E=Lagrange(1)
f(x,y)=2*f1(x,y)
vbc=[BoundaryCondition(0,f0) for i=1:4]
bvp=BoundaryValueProblem(f,vbc)

a=Assembly(msh,bvp,E,4)
nt=size(a.msh.tri,1)
for k=1:nt
    CalculateContributionOnTriangle!(a,k)
end
A=sparse(a.I,a.J,a.V)
dirt=Int[]
for ib=1:4
    global dirt = [dirt; a.msh.boundary[ib][:,1]]
end
A[dirt,:].=0
for i in dirt
    A[i,i]=1
end
a.b[dirt].=0
U=A\a.b

write_matfile("data/poisson.mat",msh=a.msh,U=U)

#end
