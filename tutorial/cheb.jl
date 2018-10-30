using Polynomials
using Plots

m=15; n=2*m; rn=0:n
ch2(j)=cos(j*π/n)
nodes=ch2.(rn)
w=ones(n+1);w[1]=w[end]=1/2;w=w.*(-1).^rn
dl=zeros(size(nodes))
for i=1:n+1
    if i!=m+1
        dl[i]=w[m+1]/w[i]/(nodes[i]-nodes[m+1])
    end
end
dl[m+1]=-sum(dl)

#r(x)=1/(1+25*x^2); f=r.(x)
x=-1:0.002:1
numer = zeros(size(x))
dnumer = zeros(size(x))
denom = zeros(size(x))
for j = 1:n+1
    xdiff = x.-nodes[j]
    temp = w[j]./xdiff
    if j==m+1
        global numer = numer + temp
    end
    global dnumer = dnumer + dl[j]*temp
    global denom = denom + temp
end
lx = numer./denom
dlx = dnumer./denom

plot(x,lx)
plot!(x,dlx)
