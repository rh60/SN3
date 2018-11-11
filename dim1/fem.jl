include("common.jl")

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

function fem(bvp::BoundaryValueProblem, E, msh::Mesh, nq::Integer)
    L, dL, Q = Quadrature(E,nq)
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
    write_matfile("data/LAb.mat", A=A, b=b)
    u=A\b
    return u
end

function test(E,N::Integer,n::Integer,nq::Integer=3)
    p(x) = -x
    q(x) = -(2*x+3)
    r(x) = x+2
    f(x) = x^3*exp(2*x)
    u_a = BoundaryCondition(0,exp(1))
    eta_b = BoundaryCondition(1,-5*exp(1)*exp(1))
    bvp = BoundaryValueProblem(p,q,r,f,u_a,eta_b)

    ele = E(n)

    msh = Mesh(ele,0.5,1,N)
    @time U=fem(bvp,ele,msh,nq)

    e=exp(1)
    E=3*e+sqrt(e)/4
    K1=3*e-32/31*E
    K2=8/31*E
    u(x)=exp(x)*(K1+K2*x^3)+exp(2*x)*(x^2-2*x+2)
    #mxcall(:plot,0,msh.x,u.(msh.x)-U)

    norm(u.(msh.x)-U)
end

N=10
r1=test(Lagrange,1000*N,1,3)
r2=test(Lagrange,N,6,8)
r3=test(Chebyshev,N,6,8)
@show r1 r2 r3
