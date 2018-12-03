include("common.jl")

n=4
Q=QuadTri(n)

de=", "

for i=0:n
    for j=0:n-i
        f(x,y)=x^i*y^j
            println(i,de,j,de,intSum(Q,f))
    end
end
