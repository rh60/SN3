struct Mesh
    x::Vector{Float64}
    N::Int
    n::Int
    function Mesh(E::Lagrange,a,b,N)
        n=E.degree
        r=r=range(a,stop=b,length=n*N+1)
        x=collect(r)
        new(x,N,n)
    end
    function Mesh(E::Chebyshev,a,b,N)
        n=length(E.nodes)-1
        if n==1
            return Mesh(a,b,N)
        end
        lx=N*n+1
        x=zeros(lx)
        r=range(a,stop=b,length=N+1)
        h=Float64(r.step)
        ix=1
        internal=E.nodes[2:end-1]
        for p in r
            x[ix]=p
            if ix==lx
                break
            end
            x[ix+1:ix+n-1]=h*internal .+ x[ix]
            ix += n
        end
        new(x,N,n)
    end
end

function testmesh()
    E=Chebyshev(4)
    m1=Mesh(E,0,10,5)
    E=Lagrange(1)
    m2=Mesh(E,0,10,5)
    m1,m2
end
