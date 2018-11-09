using FastGaussQuadrature
using LinearAlgebra
using MATLAB

struct BoundaryCondition
    code::Integer # 0 Dirichlet 1 Newton
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

struct Quad
    Points::Vector{Float64}
    Weights::Vector{Float64}
end

@inline function intSum(values::Vector{Float64},Q::Quad)
    values⋅Q.Weights;
end

struct CardinalBase
    l
    dl
end

isdefined(Main,:Lagrange) ||
const Lagrange =
[
    CardinalBase(
        [x->1-x, x->x],
        [x->-1, x->1]
    )

    CardinalBase(
        [x->2 * (x - 1/2) * (x - 1), x->-4 * x * (x - 1), x->2 * x * (x - 1/2)],
        [x->4 * x - 3, x->-8 * x + 4, x->4 * x - 1]
    )

    CardinalBase(
        [x->-9/2 * (x - 1/3) * (x - 2/3) * (x - 1),
        x->27/2 * x * (x - 2/3) * (x - 1),
        x->-27/2 * x * (x - 1/3) * (x - 1),
        x->9/2 * x * (x - 1/3) * (x - 2/3)],
        [x->-9/2 * (x - 2/3) * (x - 1) - 9/2 * (x - 1/3) * (x - 1) - 9/2 * (x - 1/3) * (x - 2/3),
        x->27/2 * (x - 2/3) * (x - 1) + 27/2 * x * (x - 1) + 27/2 * x * (x - 2/3),
        x->-27/2 * (x - 1/3) * (x - 1) - 27/2 * x * (x - 1) - 27/2 * x * (x - 1/3),
        x->9/2 * (x - 1/3) * (x - 2/3) + 9/2 * x * (x - 2/3) + 9/2 * x * (x - 1/3)]
    )
]

struct Mesh
    x::Vector{Float64}
    N::Integer
    n::Integer
    function Mesh(a,b,N,n=1)
        h=(b-a)/N
        x=collect(a:h/n:b)
        return new(x,N,n)
    end
end

function prequad(nq::Integer,B::CardinalBase)
    x,w = gausslegendre(nq)
    Q=Quad( (x.+1)/2, w/2 )
    n=length(B.l)
    L=Vector{Vector{Float64}}(undef,n)
    dL=Vector{Vector{Float64}}(undef,n)
    for i = 1:n
        L[i]=B.l[i].(Q.Points)
        dL[i]=B.dl[i].(Q.Points)
    end
    return L, dL, Q
end

@inline function CalculateContributionOnElement(problem::BoundaryValueProblem,
    x::Vector{Float64},
    Q::Quad, L::Vector{Vector{Float64}}, dL::Vector{Vector{Float64}})
    h=x[end]-x[1]
    t=x[1].+Q.Points*h
    p=problem.p.(t)
    q=problem.q.(t)
    r=problem.r.(t)
    f=problem.f.(t)
    nl=length(L)
    A=zeros(nl,nl)
    for i=1:nl
        for j=1:nl
            temp = 1/h^2 * p .* dL[i] .* dL[j] +
            1/h * q .* L[i] .* dL[j] +
            r .* L[i] .* L[j]
            A[i,j] = h * intSum(temp,Q)
        end
    end
    b=zeros(nl)
    for i=1:nl
        temp = f .* L[i]
        b[i] = h * intSum(temp,Q)
    end
    #@show A b
    return A,b
end

function fem(bvp::BoundaryValueProblem, msh::Mesh, nq::Integer)
    L, dL, Q = prequad(nq,Lagrange[msh.n])
    ndof = msh.N * msh.n + 1
    A=zeros(ndof,ndof)
    b=zeros(ndof)
    ix = 1
    for k=1:msh.N
        ir = ix:ix+msh.n
        x = msh.x[ir]
        locA,locb = CalculateContributionOnElement(bvp,x,Q,L,dL)
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
    #write_matfile("data/Ab.mat", A=A, b=b)
    u=A\b
    return u
end

function test(N::Integer,n::Integer,nq::Integer=3)
    p(x) = -x
    q(x) = -(2*x+3)
    r(x) = x+2
    f(x) = x^3*exp(2*x)
    u_a = BoundaryCondition(0,exp(1))
    eta_b = BoundaryCondition(1,-5*exp(1)*exp(1))
    bvp = BoundaryValueProblem(p,q,r,f,u_a,eta_b)

    msh = Mesh(0.5,1,N,n)
    @time U=fem(bvp,msh,nq)

    e=exp(1)
    E=3*e+sqrt(e)/4
    K1=3*e-32/31*E
    K2=8/31*E
    u(x)=exp(x)*(K1+K2*x^3)+exp(2*x)*(x^2-2*x+2)
    #mxcall(:plot,0,msh.x,u.(msh.x)-U)
    norm(u.(msh.x)-U)
end

r1=test(100,1,3)
r2=test(100,2,3)
r3=test(100,3,3)
@show r1 r2 r3
