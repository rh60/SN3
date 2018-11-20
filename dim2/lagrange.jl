function LagrangeBase(n::Int)
        if n==1
                [(x,y) -> 1 - x - y, (x,y) -> x, (x,y) -> y],
                [(x,y) -> -1, (x,y) -> 1, (x,y) -> 0],
                [(x,y) -> -1, (x,y) -> 0, (x,y) -> 1]
        elseif n==2
                [(x,y) -> (1 - x - y) * (1 - 2 * x - 2 * y)
                (x,y) -> x * (2 * x - 1)
                (x,y) -> y * (2 * y - 1)
                (x,y) -> 2 * x * (2 - 2 * x - 2 * y)
                (x,y) -> 4 * x * y
                (x,y) -> 2 * y * (2 - 2 * x - 2 * y)],
                [(x,y) -> -3 + 4 * x + 4 * y
                (x,y) -> 4 * x - 1
                (x,y) -> 0
                (x,y) -> 4 - 8 * x - 4 * y
                (x,y) -> 4 * y
                (x,y) -> -4 * y],
                [(x,y) -> -3 + 4 * x + 4 * y
                (x,y) -> 0
                (x,y) -> 4 * y - 1
                (x,y) -> -4 * x
                (x,y) -> 4 * x
                (x,y) -> 4 - 4 * x - 8 * y]
        elseif n==3
                [(x,y) -> (1 - x - y) * (1 - 3/2 * x - 3/2 * y) * (1 - 3 * x - 3 * y)
                (x,y) -> x * (3/2 * x - 1/2) * (3 * x - 2)
                (x,y) -> y * (3/2 * y - 1/2) * (3 * y - 2)
                (x,y) -> 9/2 * x * (1 - x - y) * (2 - 3 * x - 3 * y)
                (x,y) -> 3/2 * x * (3 * x - 1) * (3 - 3 * x - 3 * y)
                (x,y) -> 9/2 * x * (3 * x - 1) * y
                (x,y) -> 9/2 * x * y * (3 * y - 1)
                (x,y) -> 3/2 * y * (3 * y - 1) * (3 - 3 * x - 3 * y)
                (x,y) -> 9/2 * y * (1 - x - y) * (2 - 3 * x - 3 * y)
                (x,y) -> 9 * x * y * (3 - 3 * x - 3 * y)],
                [(x,y) -> -(1 - 3/2 * x - 3/2 * y) * (1 - 3 * x - 3 * y) - 3/2 * (1 - x - y) * (1 - 3 * x - 3 * y) - 3 * (1 - x - y) * (1 - 3/2 * x - 3/2 * y)
                (x,y) -> (3/2 * x - 1/2) * (3 * x - 2) + 3/2 * x * (3 * x - 2) + 3 * x * (3/2 * x - 1/2)
                (x,y) -> 0
                (x,y) -> 9/2 * (1 - x - y) * (2 - 3 * x - 3 * y) - 9/2 * x * (2 - 3 * x - 3 * y) - 27/2 * x * (1 - x - y)
                (x,y) -> 3/2 * (3 * x - 1) * (3 - 3 * x - 3 * y) + 9/2 * x * (3 - 3 * x - 3 * y) - 9/2 * x * (3 * x - 1)
                (x,y) -> 9/2 * (3 * x - 1) * y + 27/2 * x * y
                (x,y) -> 9/2 * y * (3 * y - 1)
                (x,y) -> -9/2 * y * (3 * y - 1)
                (x,y) -> -9/2 * y * (2 - 3 * x - 3 * y) - 27/2 * y * (1 - x - y)
                (x,y) -> 9 * y * (3 - 3 * x - 3 * y) - 27 * x * y],
                [(x,y) -> -(1 - 3/2 * x - 3/2 * y) * (1 - 3 * x - 3 * y) - 3/2 * (1 - x - y) * (1 - 3 * x - 3 * y) - 3 * (1 - x - y) * (1 - 3/2 * x - 3/2 * y)
                (x,y) -> 0
                (x,y) -> (3/2 * y - 1/2) * (3 * y - 2) + 3/2 * y * (3 * y - 2) + 3 * y * (3/2 * y - 1/2)
                (x,y) -> -9/2 * x * (2 - 3 * x - 3 * y) - 27/2 * x * (1 - x - y)
                (x,y) -> -9/2 * x * (3 * x - 1)
                (x,y) -> 9/2 * x * (3 * x - 1)
                (x,y) -> 9/2 * x * (3 * y - 1) + 27/2 * x * y
                (x,y) -> 3/2 * (3 * y - 1) * (3 - 3 * x - 3 * y) + 9/2 * y * (3 - 3 * x - 3 * y) - 9/2 * y * (3 * y - 1)
                (x,y) -> 9/2 * (1 - x - y) * (2 - 3 * x - 3 * y) - 9/2 * y * (2 - 3 * x - 3 * y) - 27/2 * y * (1 - x - y)
                (x,y) -> 9 * x * (3 - 3 * x - 3 * y) - 27 * x * y]
        elseif n==4
                [(x,y) -> (1 - x - y) * (1 - 4/3 * x - 4/3 * y) * (1 - 2 * x - 2 * y) * (1 - 4 * x - 4 * y)
                (x,y) -> x * (4/3 * x - 1/3) * (2 * x - 1) * (4 * x - 3)
                (x,y) -> y * (4/3 * y - 1/3) * (2 * y - 1) * (4 * y - 3)
                (x,y) -> 16/3 * x * (1 - x - y) * (3/2 - 2 * x - 2 * y) * (2 - 4 * x - 4 * y)
                (x,y) -> 4 * x * (4 * x - 1) * (1 - x - y) * (3 - 4 * x - 4 * y)
                (x,y) -> 4/3 * x * (2 * x - 1/2) * (4 * x - 2) * (4 - 4 * x - 4 * y)
                (x,y) -> 16/3 * x * (2 * x - 1/2) * (4 * x - 2) * y
                (x,y) -> 4 * x * (4 * x - 1) * y * (4 * y - 1)
                (x,y) -> 16/3 * x * y * (2 * y - 1/2) * (4 * y - 2)
                (x,y) -> 4/3 * y * (2 * y - 1/2) * (4 * y - 2) * (4 - 4 * x - 4 * y)
                (x,y) -> 4 * y * (4 * y - 1) * (1 - x - y) * (3 - 4 * x - 4 * y)
                (x,y) -> 16/3 * y * (1 - x - y) * (3/2 - 2 * x - 2 * y) * (2 - 4 * x - 4 * y)
                (x,y) -> 32 * x * y * (1 - x - y) * (3 - 4 * x - 4 * y)
                (x,y) -> 8 * x * (4 * x - 1) * y * (4 - 4 * x - 4 * y)
                (x,y) -> 8 * x * y * (4 * y - 1) * (4 - 4 * x - 4 * y)],
                [(x,y) -> -(1 - 4/3 * x - 4/3 * y) * (1 - 2 * x - 2 * y) * (1 - 4 * x - 4 * y) - 4/3 * (1 - x - y) * (1 - 2 * x - 2 * y) * (1 - 4 * x - 4 * y) - 2 * (1 - x - y) * (1 - 4/3 * x - 4/3 * y) * (1 - 4 * x - 4 * y) - 4 * (1 - x - y) * (1 - 4/3 * x - 4/3 * y) * (1 - 2 * x - 2 * y)
                (x,y) -> (4/3 * x - 1/3) * (2 * x - 1) * (4 * x - 3) + 4/3 * x * (2 * x - 1) * (4 * x - 3) + 2 * x * (4/3 * x - 1/3) * (4 * x - 3) + 4 * x * (4/3 * x - 1/3) * (2 * x - 1)
                (x,y) -> 0
                (x,y) -> 16/3 * (1 - x - y) * (3/2 - 2 * x - 2 * y) * (2 - 4 * x - 4 * y) - 16/3 * x * (3/2 - 2 * x - 2 * y) * (2 - 4 * x - 4 * y) - 32/3 * x * (1 - x - y) * (2 - 4 * x - 4 * y) - 64/3 * x * (1 - x - y) * (3/2 - 2 * x - 2 * y)
                (x,y) -> 4 * (4 * x - 1) * (1 - x - y) * (3 - 4 * x - 4 * y) + 16 * x * (1 - x - y) * (3 - 4 * x - 4 * y) - 4 * x * (4 * x - 1) * (3 - 4 * x - 4 * y) - 16 * x * (4 * x - 1) * (1 - x - y)
                (x,y) -> 4/3 * (2 * x - 1/2) * (4 * x - 2) * (4 - 4 * x - 4 * y) + 8/3 * x * (4 * x - 2) * (4 - 4 * x - 4 * y) + 16/3 * x * (2 * x - 1/2) * (4 - 4 * x - 4 * y) - 16/3 * x * (2 * x - 1/2) * (4 * x - 2)
                (x,y) -> 16/3 * (2 * x - 1/2) * (4 * x - 2) * y + 32/3 * x * (4 * x - 2) * y + 64/3 * x * (2 * x - 1/2) * y
                (x,y) -> 4 * (4 * x - 1) * y * (4 * y - 1) + 16 * x * y * (4 * y - 1)
                (x,y) -> 16/3 * y * (2 * y - 1/2) * (4 * y - 2)
                (x,y) -> -16/3 * y * (2 * y - 1/2) * (4 * y - 2)
                (x,y) -> -4 * y * (4 * y - 1) * (3 - 4 * x - 4 * y) - 16 * y * (4 * y - 1) * (1 - x - y)
                (x,y) -> -16/3 * y * (3/2 - 2 * x - 2 * y) * (2 - 4 * x - 4 * y) - 32/3 * y * (1 - x - y) * (2 - 4 * x - 4 * y) - 64/3 * y * (1 - x - y) * (3/2 - 2 * x - 2 * y)
                (x,y) -> 32 * y * (1 - x - y) * (3 - 4 * x - 4 * y) - 32 * x * y * (3 - 4 * x - 4 * y) - 128 * x * y * (1 - x - y)
                (x,y) -> 8 * (4 * x - 1) * y * (4 - 4 * x - 4 * y) + 32 * x * y * (4 - 4 * x - 4 * y) - 32 * x * (4 * x - 1) * y
                (x,y) -> 8 * y * (4 * y - 1) * (4 - 4 * x - 4 * y) - 32 * x * y * (4 * y - 1)],
                [(x,y) -> -(1 - 4/3 * x - 4/3 * y) * (1 - 2 * x - 2 * y) * (1 - 4 * x - 4 * y) - 4/3 * (1 - x - y) * (1 - 2 * x - 2 * y) * (1 - 4 * x - 4 * y) - 2 * (1 - x - y) * (1 - 4/3 * x - 4/3 * y) * (1 - 4 * x - 4 * y) - 4 * (1 - x - y) * (1 - 4/3 * x - 4/3 * y) * (1 - 2 * x - 2 * y)
                (x,y) -> 0
                (x,y) -> (4/3 * y - 1/3) * (2 * y - 1) * (4 * y - 3) + 4/3 * y * (2 * y - 1) * (4 * y - 3) + 2 * y * (4/3 * y - 1/3) * (4 * y - 3) + 4 * y * (4/3 * y - 1/3) * (2 * y - 1)
                (x,y) -> -16/3 * x * (3/2 - 2 * x - 2 * y) * (2 - 4 * x - 4 * y) - 32/3 * x * (1 - x - y) * (2 - 4 * x - 4 * y) - 64/3 * x * (1 - x - y) * (3/2 - 2 * x - 2 * y)
                (x,y) -> -4 * x * (4 * x - 1) * (3 - 4 * x - 4 * y) - 16 * x * (4 * x - 1) * (1 - x - y)
                (x,y) -> -16/3 * x * (2 * x - 1/2) * (4 * x - 2)
                (x,y) -> 16/3 * x * (2 * x - 1/2) * (4 * x - 2)
                (x,y) -> 4 * x * (4 * x - 1) * (4 * y - 1) + 16 * x * (4 * x - 1) * y
                (x,y) -> 16/3 * x * (2 * y - 1/2) * (4 * y - 2) + 32/3 * x * y * (4 * y - 2) + 64/3 * x * y * (2 * y - 1/2)
                (x,y) -> 4/3 * (2 * y - 1/2) * (4 * y - 2) * (4 - 4 * x - 4 * y) + 8/3 * y * (4 * y - 2) * (4 - 4 * x - 4 * y) + 16/3 * y * (2 * y - 1/2) * (4 - 4 * x - 4 * y) - 16/3 * y * (2 * y - 1/2) * (4 * y - 2)
                (x,y) -> 4 * (4 * y - 1) * (1 - x - y) * (3 - 4 * x - 4 * y) + 16 * y * (1 - x - y) * (3 - 4 * x - 4 * y) - 4 * y * (4 * y - 1) * (3 - 4 * x - 4 * y) - 16 * y * (4 * y - 1) * (1 - x - y)
                (x,y) -> 16/3 * (1 - x - y) * (3/2 - 2 * x - 2 * y) * (2 - 4 * x - 4 * y) - 16/3 * y * (3/2 - 2 * x - 2 * y) * (2 - 4 * x - 4 * y) - 32/3 * y * (1 - x - y) * (2 - 4 * x - 4 * y) - 64/3 * y * (1 - x - y) * (3/2 - 2 * x - 2 * y)
                (x,y) -> 32 * x * (1 - x - y) * (3 - 4 * x - 4 * y) - 32 * x * y * (3 - 4 * x - 4 * y) - 128 * x * y * (1 - x - y)
                (x,y) -> 8 * x * (4 * x - 1) * (4 - 4 * x - 4 * y) - 32 * x * (4 * x - 1) * y
                (x,y) -> 8 * x * (4 * y - 1) * (4 - 4 * x - 4 * y) + 32 * x * y * (4 - 4 * x - 4 * y) - 32 * x * y * (4 * y - 1)]
        elseif n==5
                [],
                [],
                []
        elseif n==6
                [],
                [],
                []
        elseif n==7
                [],
                [],
                []
        elseif n==8
                [],
                [],
                []
        end
end

struct Lagrange
        l
        l1
        l2
        degree
        function Lagrange(degree::Int)
                l,l1,l2=LagrangeBase(degree)
                new(l,l1,l2,degree)
        end
end

function Quadrature(E::Lagrange,nq::Int)
        Q=Quad(nq)
        n=(E.degree+1)*(E.degree+2)/2
        L=Vector{Vector{Float64}}(undef,n)
        L1=Vector{Vector{Float64}}(undef,n)
        L2=Vector{Vector{Float64}}(undef,n)
        for i = 1:n
                L[i]=E.l[i].(Q.Points[:,1],Q.Points[:,2])
                L1[i]=E.l1[i].(Q.Points[:,1],Q.Points[:,2])
                L2[i]=E.l2[i].(Q.Points[:,1],Q.Points[:,2])
        end
        return L, L1, L2, Q
end

function LinearTranform(elx::Vector{Float64},ely::Vector{Float64})
        E=Lagrange(1)
        @inline ϕ(x,y)=[elx[1]*E.l[1](x,y)+elx[2]*E.l[2](x,y)+elx[3]*E.l[3](x,y)
                        ely[1]*E.l[1](x,y)+ely[2]*E.l[2](x,y)+ely[3]*E.l[3](x,y)]

        Jϕ=det([-elx[1]+elx[2] -elx[1]+elx[3]
                -ely[1]+ely[2] -ely[1]+ely[3]])

        return ϕ, Jϕ
end
