include("common.jl")

#function test()
R=Rectangle(0,2,0,1)
msh=Mesh(R,10,5)
el=msh.tri[5,:]
elx=msh.x[el]
ely=msh.y[el]
ϕ, Jϕ = LinearTranform(elx,ely)
p=ϕ(1/2,1/2)
#end
