include("common.jl")

E=Chebyshev(30)
x = collect(-5:0.01:5)
f(x)=1/(1+x^2)
u=Interpolate(E,f,x)
mxcall(:plot,0,x,f.(x),x,u)
