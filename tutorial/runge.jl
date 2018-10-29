using Polynomials
using Plots

n=100;rn=0:n
ch2(j)=cos(j*π/n)
x=ch2.(rn)
w=ones(n+1);w[1]=w[end]=1/2;w=w.*(-1).^rn
#f=zeros(n+1); f[26]=1;
r(x)=1/(1+25*x^2); f=r.(x)
xx=-1:0.001:1
numer = zeros(size(xx))
denom = zeros(size(xx))
for j = 1:n+1
    xdiff = xx.-x[j]
    temp = w[j]./xdiff
    global numer = numer + temp*f[j]
    global denom = denom + temp
end
ff = numer./denom

plot(xx,ff)
plot!(xx,r.(xx))
