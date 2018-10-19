#push!(LOAD_PATH,raw"D:\documents\Julia\Ming")

using Ming

x=linspace(0.5, 1, 101)
@time U=fem2(x)

using Plots
plot(x,U,m=2,label="U")
