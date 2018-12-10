struct Assembly
    msh::Mesh
    bvp::BoundaryValueProblem
    L::Vector{Vector{Float64}}
    L1::Vector{Vector{Float64}}
    L2::Vector{Vector{Float64}}
    Q::QuadTri
    Ls::Vector{Vector{Float64}}
    Qs::QuadSeg
    I::Vector{Int}
    J::Vector{Int}
    V::Vector{Float64}
    b::Vector{Float64}
    function Assembly(msh::Mesh,bvp::BoundaryValueProblem,E::Lagrange,nq::Int)
        L,L1,L2,Q,Ls,Qs=Quadrature(E,nq)
        b=zeros(length(msh.x))
        new(msh,bvp,L,L1,L2,Q,Ls,Qs,Int[],Int[],Float64[],b)
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
            tmp=intSum(a.Q,
            p.*g+(q1.*gj[1,:]+q2.*gj[2,:]+r.*a.L[j]).*a.L[i]
            )*J
            push!(a.V,tmp)
        end
        #a.b[t[i]] += intSum(a.Q,f.*a.L[i])*J
        a.b[loc2glob[i]] += intSum(a.Q,f.*a.L[i])*J
    end
end

@inline function CalculateContributionOnSegment!(a::Assembly,ib::Int,k::Int)
    l=length(a.Ls)
    loc2glob=a.msh.boundary[ib][k,:]
    t=loc2glob[1:2]
    x=a.msh.x[t];y=a.msh.y[t];
    h=norm([diff(x),diff(y)])
    qx=affinity(x,a.Qs.Points)
    qy=affinity(y,a.Qs.Points)
    g=a.bvp.vbc[ib].value.(qx,qy)
    for i=1:l
        a.b[loc2glob[i]] += intSum(a.Qs,g.*a.Ls[i])*h
    end
end

function Solve(bvp::BoundaryValueProblem,msh::Mesh,degree::Int,nq::Int)
    E=Lagrange(degree)
    a=Assembly(LagrangeMesh(msh,E.degree),bvp,E,nq)
    nt=size(a.msh.tri,1)
    for k=1:nt
        CalculateContributionOnTriangle!(a,k)
    end
    A=sparse(a.I,a.J,a.V)

    for (ib,bc) in enumerate(a.bvp.vbc)
        if bc.code==Neumann
            for k=1:size(a.msh.boundary[ib],1)
                CalculateContributionOnSegment!(a,ib,k)
            end
        end
    end

    for (ib,bc) in enumerate(a.bvp.vbc)
        if bc.code==Dirichlet
            inodes = a.msh.boundary[ib][:]
            A[inodes,:].=0
            for i in inodes
                A[i,i]=1
            end
            x=a.msh.x[inodes]
            y=a.msh.y[inodes]
            a.b[inodes]=bc.value.(x,y)
        end
    end

    U=A\a.b
    return U, a, A
end
