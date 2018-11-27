using MATLAB

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

function numPointsLine(n::Int)::Int
    n-1
end

function refTriPattern(n::Int,e1=true,e2=true,e3=true)
    np=numPointsTri(n)-3
    x=Float64[]; sizehint!(x,np)
    y=Float64[]; sizehint!(x,np)
    h=1//n; z=zeros(n-1)
    s=Float64.((1:n-1)*h); r=reverse(s)
    if e1
        append!(x,s); append!(y,z)
    end
    if e2
        append!(x,r); append!(y,s)
    end
    if e3
        append!(x,z); append!(y,r)
    end
    for j=1:n-1
        for i=1:n-j-1
            push!(x,Float64(i*h))
            push!(y,Float64(j*h))
        end
    end
    return Pattern(x,y)
end

@inline function appendPoints!(x,y,T,p::Pattern)
    Px=x[T]; Py=y[T]
    append!(x,affinity(Px,p.x,p.y))
    append!(y,affinity(Py,p.x,p.y))
    return length(x)
end

@inline function appendIndices!(indices::Vector{Int},count::Int,shift::Int)::Int
    ix=1:count
    append!(indices,shift .+ ix)
    return shift+count
end

@inline function createIndices(shift::Int,n::Int,e1=Int[],e2=Int[],e3=Int[])
    n1=numPointsLine(n)
    edges=(e1,e2,e3)
    indices=Int[];  sizehint!(indices,numPointsTri(n))
    for e in edges
        if length(e)==0
            shift=appendIndices!(indices,n1,shift)
        else
            append!(indices,e)
        end
    end
    n2=numInternalPointsTri(n)
    appendIndices!(indices,n2,shift)
    return indices, [1:n1, n1+1:2*n1, 2*n1+1:3*n1]
end

@inline function selectEdgeIndices(i::Vector{Int},k::Int)

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
    npt=numPoints(n)
    tri=zeros(Int,nt,npt)
    for (k,(i,j)) in enumerate(Iterators.product(1:nx,1:ny))
        r=[lin(i,j) lin(i+1,j) lin(i+1,j+1) lin(i,j+1)]
        tri[2*k-1,1:3]=[r[1] r[2] r[4]]
        tri[2*k,1:3]=[r[2] r[3] r[4]]
    end

    # Create nodes of Lagrange element
    # index of the last point
    shift=length(x)
    if n>1
        p1=refTriPattern(n)
        p2=refTriPattern(n,true,true,false)
        ltp=length(tp1.x)
        nl=numPointsLine(n)
        for k=1:2:nt
            temp,ie=createIndices(shift,n)
            tri[k,4:end]=temp
            temp=temp[ie[2]]
            shift=appendPoints!(x,y,tri[k,1:3],p1)
            tri[k+1,4:end]=createIndices(shift,n,[],[],temp)[1]
            shift=appendPoints!(x,y,tri[k+1,1:3],p2)
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

T=refTriPattern(4)
msh=testMesh()
