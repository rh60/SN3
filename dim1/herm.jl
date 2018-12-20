using LinearAlgebra
using Polynomials
using FastGaussQuadrature
using PGFPlotsX

l=1
ϕ1=Poly([0,1])
ϕ2=Poly([1,-1])
u=16*ϕ1^2*ϕ2^2
f(x)=384
P(h) = [Poly([1,0,-3 / h ^ 2,2 / h ^ 3]),
        Poly([0,1,-2 / h,1 / h ^ 2]),
        Poly([0,0,3 / h ^ 2,-2 / h ^ 3]),
        Poly([0,0,-0.1e1 / h,0.1e1 / h ^ 2])]
K(h) = [12 / h ^ 3 6 / h ^ 2 -12 / h ^ 3 6 / h ^ 2;
        6 / h ^ 2 4 / h -6 / h ^ 2 2 / h;
        -12 / h ^ 3 -6 / h ^ 2 12 / h ^ 3 -6 / h ^ 2;
        6 / h ^ 2 2 / h -6 / h ^ 2 4 / h]
F(h) = [h / 2,h ^ 2 / 12,h / 2,-h ^ 2 / 12]

function plot(x,u,U,n)
    np=div(length(x),n)
    figure = @pgf Axis(
    {
        height = "10cm",
        width = "12cm",
        xlabel = raw"$x$",
        ylabel = raw"$u=16x^2(1-x)^2$",
    },
    Plot(
        {
            color = "black",
        },
        Coordinates(x,u),
    ),
    LegendEntry(raw"$u$"),
    Plot(
        {
            color = "red",
        },
        Coordinates(x,U),
    ),
    LegendEntry("\$u_{\\frac{1}{$n}}\$"),
    Plot(
        {
            color = "red",
            "only marks",
        },
        Coordinates(x[1:np:end],U[1:np:end]),
    ),)
end

function L2(pp::Vector{Poly{Float64}},nq=5)
    n=length(pp)
    h=l/n
    x,w=gausslegendre(nq)
    x=h*(1 .+ x)/2
    w=w*h/2
    r=0
    for i=0:n-1
        xh=i*h .+ x
        y=(u(xh)-pp[i+1](x)).^2
        r+=dot(y,w)
    end
    sqrt(r)
end

function compute(pp::Vector{Poly{Float64}},np=50)
    n=length(pp)
    h=l/n
    x=collect(0:h/np:h)
    U=pp[1].(x)
    x=x[2:end]
    for i=2:n
        U=[U;pp[i].(x)]
    end
    return collect(0:h/np:l),U
end

function femherm(n::Int)
        ndof=2*(n+1)
        A=zeros(ndof,ndof)
        b=zeros(ndof)
        h=l/n
        Ke=K(h)
        Fe=f(0)*F(h)
        for k=1:n
            r=2*k-1:2*k+2
            A[r,r]+=Ke
            b[r]+=Fe
        end
        A[1:2,:].=0;A[1,1]=1;A[2,2]=1
        A[end-1:end,:].=0;A[end-1,end-1]=1;A[end,end]=1
        b[1:2].=0;b[end-1:end].=0
        c=A\b
        B=P(h)
        pp=[c[k]*B[1]+c[k+1]*B[2]+c[k+2]*B[3]+c[k+3]*B[4] for k=1:2:2*n]
end

n=4
pp=femherm(n)
L2error=L2(pp)
x,U=compute(pp)
figure=plot(x,u.(x),U,n)
pgfsave("d:/data/obr$n.pdf", figure)
display(pp)
