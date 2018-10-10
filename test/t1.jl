push!(LOAD_PATH,raw"D:\documents\Julia\Ming")
using Ming

t=Ming.linspace(0,2pi,1000)
x=cos.(t)
y=sin.(t)
@msave t x y

xa = Ming.load("x")
