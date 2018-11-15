include("common.jl")

struct BoundaryCondition
    code::Int # 0 Dirichlet 1 Newton
    value::Float64
end

struct BoundaryValueProblem
    p # difuzní koeficient
    q # konvekce
    r # reakce
    f # pravá strana
    Left::BoundaryCondition # okrajová podmínka vlevo
    Right::BoundaryCondition # vpravo
end

@inline function CalculateContributionOnElement!(problem::BoundaryValueProblem,
    x::Vector{Float64}, A::Matrix{Float64}, b::Vector{Float64},
    Q::Quad, L::Vector{Vector{Float64}}, dL::Vector{Vector{Float64}})
    h=x[end]-x[1]
    t=x[1].+Q.Points*h
    p=problem.p.(t)
    q=problem.q.(t)
    r=problem.r.(t)
    f=problem.f.(t)
    nl=length(L)
    for i=1:nl
        for j=1:nl
            temp = 1/h^2 * p .* dL[i] .* dL[j] +
            1/h * q .* L[i] .* dL[j] +
            r .* L[i] .* L[j]
            A[i,j] = h * intSum(temp,Q)
        end
    end
    for i=1:nl
        temp = f .* L[i]
        b[i] = h * intSum(temp,Q)
    end    
end

@inline function CalculateContributionOnElement!(problem::BoundaryValueProblem,
    x::Vector{Float64}, globix::Int, ins::Int,
    I::Vector{Int}, J::Vector{Int}, V::Vector{Float64},
    b::Vector{Float64},
    Q::Quad, L::Vector{Vector{Float64}}, dL::Vector{Vector{Float64}})
    h=x[end]-x[1]
    t=x[1].+Q.Points*h
    p=problem.p.(t)
    q=problem.q.(t)
    r=problem.r.(t)
    f=problem.f.(t)
    nl=length(L)
    for i=1:nl
        for j=1:nl
            temp = 1/h^2 * p .* dL[i] .* dL[j] +
            1/h * q .* L[i] .* dL[j] +
            r .* L[i] .* L[j]
            I[ins] = globix+i
            J[ins] = globix+j
            V[ins] = h * intSum(temp,Q)
            ins += 1
        end
    end
    for i=1:nl
        temp = f .* L[i]
        b[i] = h * intSum(temp,Q)
    end
    ins
end

function fem(bvp::BoundaryValueProblem, E, msh::Mesh, nq::Int)
    L, dL, Q = Quadrature(E,nq)
    ndof = msh.N * msh.n + 1
    A=zeros(ndof,ndof)
    b=zeros(ndof)
    ix = 1
    n=E.degree+1
    locA=zeros(n,n)
    locb=zeros(n)
    for k=1:msh.N
        ir = ix:ix+msh.n
        x = msh.x[ir]
        CalculateContributionOnElement!(bvp,x,locA,locb,Q,L,dL)
        A[ir,ir]+=locA
        b[ir]+=locb
        ix += msh.n
    end
    if bvp.Left.code==0
        A[1,:] .= 0
        A[1,1] = 1
        b[1]=bvp.Left.value
    else
        b[1]+=bvp.Left.value
    end
    if bvp.Right.code==0
        A[end,:] .= 0
        A[end,end] = 1
        b[end] = bvp.Right.value
    else
        b[end] += bvp.Right.value
    end
    #write_matfile("data/LAb.mat", A=A, b=b)
    u=A\b
    return u
end

function femSparse(bvp::BoundaryValueProblem, E, msh::Mesh, nq::Int)
    L, dL, Q = Quadrature(E,nq)
    ndof = msh.N * msh.n + 1
    b=zeros(ndof)
    ix = 1
    n=E.degree+1
    sz=n*n*msh.N
    @show sz
    I=zeros(Int,sz)
    J=zeros(Int,sz)
    V=zeros(Float64,sz)
    locb=zeros(n)
    ins=1
    for k=1:msh.N
        ir = ix:ix+msh.n
        x = msh.x[ir]
        ins=CalculateContributionOnElement!(bvp,x,ix-1,ins,I,J,V,locb,Q,L,dL)
        b[ir]+=locb
        ix += msh.n
    end
    A=sparse(I,J,V)
    if bvp.Left.code==0
        A[1,:] .= 0
        A[1,1] = 1
        b[1]=bvp.Left.value
    else
        b[1]+=bvp.Left.value
    end
    if bvp.Right.code==0
        A[end,:] .= 0
        A[end,end] = 1
        b[end] = bvp.Right.value
    else
        b[end] += bvp.Right.value
    end
    #@show A
    u=A\b
    return u
end

function test(E,N::Int,n::Int,nq::Int=3)
    p(x) = -x
    q(x) = -(2*x+3)
    r(x) = x+2
    f(x) = x^3*exp(2*x)
    u_a = BoundaryCondition(0,exp(1))
    eta_b = BoundaryCondition(1,-5*exp(1)*exp(1))
    bvp = BoundaryValueProblem(p,q,r,f,u_a,eta_b)

    ele = E(n)

    msh = Mesh(ele,0.5,1,N)
    @time U=femSparse(bvp,ele,msh,nq)

    e=exp(1)
    E=3*e+sqrt(e)/4
    K1=3*e-32/31*E
    K2=8/31*E
    u(x)=exp(x)*(K1+K2*x^3)+exp(2*x)*(x^2-2*x+2)
    #mxcall(:plot,0,msh.x,u.(msh.x)-U)

    norm(u.(msh.x)-U)
end

N=1000; n=1
@show N n
@show test(Lagrange,N,1,3)

N=10
for n=2:8
    @show N n
    @show test(Lagrange,N,n,n+2)
    @show test(Chebyshev,N,n,n+2)
end
