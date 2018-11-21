include("common.jl")

#function test()
R=Rectangle(0,2,0,1)
msh=Mesh(R,10,5)
el=msh.tri[5,:]
x=msh.x[el]
y=msh.y[el]
ϕ1, ϕ2, Jϕ, invJT = LinearTranform(x,y)

E=Lagrange(4)
n=length(E.x)
A=zeros(n,n)
for i=1:n
    A[i,:]=E.l[i].(E.x,E.y)
end
display(A)
px=ϕ1.(E.x,E.y)
py=ϕ2.(E.x,E.y)
@show px py
#end
