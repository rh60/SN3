push!(LOAD_PATH,raw"D:\documents\Julia\Ming")

using Ming

r=fem1(11)

using Plots

plot(r.mesh,r.solution,m=4,xlabel="x",label="U")
