using MATLAB
using SparseArrays

struct Rectangle
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64
end

struct Mesh
    x::Vector{Float64}
    y::Vector{Float64}
    tri::Matrix{Int}
    boundary::Vector{Vector{Int}}
end

struct Pattern
    x::Vector{Float64}
    y::Vector{Float64}
end

@inline function affinity(P,ξ)
    P[1]*(1 .- ξ)+P[2]*ξ
end

@inline function affinity(P,ξ,η)
    P[1]*(1 .- ξ-η)+P[2]*ξ+P[3]*η
end

function numPointsTri(n::Int)::Int
    div((n+1)*(n+2),2)
end

function numInternalPointsTri(n::Int)::Int
    numPointsTri(n-3)
end

function numPointsSeg(n::Int)::Int
    n-1
end

function refLine(n::Int)
    h=1//n
    Float64.((1:numPointsSeg(n))*h)
end

function refTri(n::Int)
    h=1//n
    np=numInternalPointsTri(n)
    x=Float64[]; sizehint!(x,np)
    y=Float64[]; sizehint!(x,np)
    for j=1:n-1
        for i=1:n-j-1
            push!(x,Float64(i*h))
            push!(y,Float64(j*h))
        end
    end
    return Pattern(x,y)
end

@inline function appendPoints!(x,y,T,p::Vector{Float64})
    Px=x[T]; Py=y[T]
    append!(x,affinity(Px,p))
    append!(y,affinity(Py,p))
    return length(x)
end

@inline function appendPoints!(x,y,T,p::Pattern)
    Px=x[T]; Py=y[T]
    append!(x,affinity(Px,p.x,p.y))
    append!(y,affinity(Py,p.x,p.y))
    return length(x)
end

mutable struct EdgeInfo
    tri::Int
    loc::Int
    created::Bool
    function EdgeInfo(t::Int,l::Int)
        new(t,l,false)
    end
end

@inline function pushEdgeInfo!(k,l,tri,I,J,Id,efo)
    li=l;
    if l==3
        lj=1
    else
        lj=l+1
    end
    i=tri[k,li]; j=tri[k,lj]
    push!(I,i); push!(J,j); push!(Id,length(I))
    push!(efo,EdgeInfo(k,l))
end

@inline function appendEdge!(k,l,tri,E,efo,edges)
    li=l
    if l==3
        lj=1
    else
        lj=l+1
    end
    i=tri[k,li]; j=tri[k,lj];
    id1=E[i,j]; id2=E[j,i]
    if !efo[id1].created
        if  id2==0
            push!(edges,[i, j, efo[id1].tri, efo[id1].loc])
            efo[id1].created=true
        else
            push!(edges,[i, j, efo[id1].tri, efo[id1].loc, efo[id2].tri, efo[id2].loc])
            efo[id1].created=true
            efo[id2].created=true
        end
    end
end

function createEdges(tri::Matrix{Int})
    I=Int[];J=Int[];Id=Int[]
    efo=EdgeInfo[]
    nt=size(tri,1)
    for k=1:nt
        for l=1:3
            pushEdgeInfo!(k,l,tri,I,J,Id,efo)
        end
    end
    E=sparse(I,J,Id)
    edges=Vector{Int}[]
    for k=1:nt
        for l=1:3
            appendEdge!(k,l,tri,E,efo,edges)
        end
    end
    return edges
end

function Mesh(R::Rectangle,nx::Int,ny::Int,n::Int=1)
    n1=nx+1; n2=ny+1
    np=n1*n2; nt=2*nx*ny
    x=Float64[]; sizehint!(x,np)
    y=Float64[]; sizehint!(y,np)
    rx=range(R.xmin,stop=R.xmax,length=n1)
    ry=range(R.ymin,stop=R.ymax,length=n2)
    for (px,py) in Iterators.product(rx,ry)
        push!(x,px); push!(y,py)
    end
    lin(i,j)=n1*(j-1)+i
    npt=numPointsTri(n)
    tri=zeros(Int,nt,npt)
    for (k,(i,j)) in enumerate(Iterators.product(1:nx,1:ny))
        r=[lin(i,j) lin(i+1,j) lin(i+1,j+1) lin(i,j+1)]
        tri[2*k-1,1:3]=[r[1] r[2] r[4]]
        tri[2*k,1:3]=[r[2] r[3] r[4]]
    end

    if n>1
        # Create nodes of Lagrange element
        # index of the last point
        shift=length(x)
        edges=createEdges(tri)
        lp=refLine(n)
        nlp=length(lp)
        rlp=1:nlp
        ne=length(edges)
        for i=1:ne
            e=edges[i]
            start=3 + (e[4]-1)*nlp
            tri[e[3],start .+ rlp] = shift .+ rlp
            shift=appendPoints!(x,y,e[1:2],lp)
            if length(e)==6
                start=3 + (e[6]-1)*nlp
                tri[e[5],start .+ rlp] = shift .+ reverse(rlp)
            end
        end
        if n>2
            tp=refTri(n)
            ntp=length(tp.x)
            rtp=1:ntp
            start=3 + 3*nlp
            for k=1:nt
                t=tri[k,1:3]
                tri[k,start .+ rtp] = shift .+ rtp
                shift=appendPoints!(x,y,t,tp)
            end
        end
    end

    bdry=[Int[],Int[],Int[],Int[]]
    for (k,(i,j)) in enumerate(Iterators.product(1:n1,1:n2))
        if i==1 push!(bdry[4],k) end
        if i==n1 push!(bdry[2],k) end
        if j==1 push!(bdry[1],k) end
        if j==n2 push!(bdry[3],k) end
    end
    bdry[3]=reverse(bdry[3])
    bdry[4]=reverse(bdry[4])
    Mesh(x,y,tri,bdry)
end

function testMesh()
    R=Rectangle(0,5,0,3)
    msh=Mesh(R,5,3,4)
    write_matfile("data/msh.mat",msh=msh)
    #mxcall(:addpath,0,raw"d:\documents\julia\Ming\dim2")
    #mxcall(:show,0)
end

msh=testMesh()
