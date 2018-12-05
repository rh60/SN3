include("common.jl")

mat=read_matfile("geo/L.mat")
val=jvalue(mat["msh"])

x=val["POS"][:,1]
y=val["POS"][:,2]
tri=Int.(val["TRIANGLES"][:,1:3])
lines=Int.(val["LINES"])
mask1=map(x->x==1,lines[:,3])
mask2=map(x->x==2,lines[:,3])
bdry=[lines[mask1,1:2], lines[mask2,1:2]]
msh=Mesh(x,y,tri,bdry)

LagrangeMesh!(msh,4)
write_matfile("data/L.mat",msh=msh)
