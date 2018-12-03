include("common.jl")

#function test()
E=Lagrange(5)

n=length(E.reftri.x)
A=zeros(n,n)
for i=1:n
    A[i,:]=E.l[i].(E.reftri.x,E.reftri.y)
end

display(A)
