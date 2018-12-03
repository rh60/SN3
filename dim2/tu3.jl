include("common.jl")
msh=CircularSector(3)
LagrangeMesh!(msh,4)
write_matfile("data/msh.mat",msh=msh)
