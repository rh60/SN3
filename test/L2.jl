push!(LOAD_PATH,raw"D:\documents\Julia\Ming")

using Polynomials
using Ming
using GR
using LinearAlgebra

μ(i,j)=1/(i+j-1)
f(x)=1/(1+25*x^2)*sin(10*π*x)
x,fw=Ming.prequad(f,0,1,100)

n=100

M=zeros(n+1,n+1)
b=zeros(n+1)
for i=1:n+1
    b[i]=x.^(i-1)⋅fw
    for j=1:n+1
        M[i,j]=μ(i,j)
    end
end

c=M\b
P=Poly(c)

t=Ming.linspace(0,1,300)

figure(size=(1000,500))
plot(t,f.(t),t,P(t))
