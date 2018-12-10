include("common.jl")

#function test()
#E=Lagrange(5)
msh=Mesh(Rectangle(0,5,0,3),5,3)
msh2=LagrangeMesh(msh,4)
write_matfile("data/msh.mat",msh=msh2)
