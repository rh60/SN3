using MATLAB
using SparseArrays

mutable struct Mesh
    x::Vector{Float64}
    y::Vector{Float64}
    tri::Matrix{Int}
    boundary::Vector{Matrix{Int}}
end

function Mesh(np::Int)
    piece=zeros(Int,0,2)
    bdry=fill(piece,np)
    Mesh(Float64[], Float64[], zeros(Int,0,3), bdry)
end

function addPoint!(msh::Mesh,x,y)
    push!(msh.x,x)
    push!(msh.y,y)
end

function addTriangle!(msh::Mesh,i1,i2,i3)
    msh.tri=[msh.tri; i1 i2 i3]
end

function addBoundaryEdge!(msh::Mesh,i,i1,i2)
    msh.boundary[i]=[msh.boundary[i]; i1 i2]
end

function szPoints!(msh::Mesh,np::Int)
    sizehint!(msh.x,np)
    sizehint!(msh.y,np)
end

mutable struct EdgeInfo
    dim::Int
    element::Int
    location::Int
    created::Bool
end

function EdgeInfo2(e::Int,l::Int)
    EdgeInfo(2,e,l,false)
end

function EdgeInfo1(e::Int,l::Int)
    EdgeInfo(1,e,l,false)
end

@inline function pushEdgeInfo1!(k,l,bdry,I,J,Id,efo)
    j=bdry[k][l,1]; i=bdry[k][l,2]
    push!(I,i); push!(J,j);
    push!(efo,EdgeInfo1(k,l));
    push!(Id,length(efo))
end

@inline function pushEdgeInfo2!(k,l,tri,I,J,Id,efo)
    li=l;
    if l==3
        lj=1
    else
        lj=l+1
    end
    i=tri[k,li]; j=tri[k,lj]
    push!(I,i); push!(J,j);
    push!(efo,EdgeInfo2(k,l));
    push!(Id,length(efo))
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
        push!(edges,[i, j, efo[id1].dim, efo[id1].element, efo[id1].location,
        efo[id2].dim, efo[id2].element, efo[id2].location])
        efo[id1].created=true
        efo[id2].created=true
    end
end

function createEdges(msh::Mesh)
    I=Int[];J=Int[];Id=Int[]
    efo=EdgeInfo[]
    nt=size(msh.tri,1)
    for k=1:nt
        for l=1:3
            pushEdgeInfo2!(k,l,msh.tri,I,J,Id,efo)
        end
    end
    nb=length(msh.boundary)
    for k=1:nb
        nl=size(msh.boundary[k],1)
        for l=1:nl
            pushEdgeInfo1!(k,l,msh.boundary,I,J,Id,efo)
        end
    end
    E=sparse(I,J,Id)
    edges=Vector{Int}[]
    for k=1:nt
        for l=1:3
            appendEdge!(k,l,msh.tri,E,efo,edges)
        end
    end
    return edges
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

struct Pattern
    x::Vector{Float64}
    y::Vector{Float64}
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

# Create nodes of Lagrange element
function LagrangeMesh!(msh::Mesh,n::Int)
    if n>1
        # index of the last point
        shift=length(msh.x)
        nt=size(msh.tri,1)
        edges=createEdges(msh)
        lp=refLine(n)
        np=length(lp)
        triedge=[zeros(Int,nt,np), zeros(Int,nt,np), zeros(Int,nt,np)]
        nb=length(msh.boundary)
        bdredge=fill(zeros(Int,0,0),nb)
        for i=1:nb
            ns=size(msh.boundary[i],1)
            bdredge[i]=zeros(Int,ns,np)
        end
        rp=1:np
        ne=length(edges)
        for i=1:ne
            e=edges[i]
            indices=shift .+ rp
            triedge[e[5]][e[4],:]=indices
            shift=appendPoints!(msh.x,msh.y,e[1:2],lp)
            if e[6]==2
                triedge[e[8]][e[7],:]=reverse(indices)
            else
                bdredge[e[7]][e[8],:]=indices
            end
        end
        msh.tri=[msh.tri triedge[1] triedge[2] triedge[3]]
        for i=1:nb
            msh.boundary[i]=[msh.boundary[i] bdredge[i]]
        end
        if n>2
            tp=refTri(n)
            np=length(tp.x)
            rp=1:np
            triinter=zeros(Int,nt,np)
            for k=1:nt
                triinter[k,:] = shift .+ rp
                shift=appendPoints!(msh.x,msh.y,msh.tri[k,1:3],tp)
            end
            msh.tri=[msh.tri triinter]
        end
    end
    return
end

function LagrangeMesh(msh::Mesh,n::Int)
    lmsh=deepcopy(msh)
    LagrangeMesh!(lmsh,n)
    return lmsh
end

struct Rectangle
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64
end

function Mesh(R::Rectangle,nx::Int,ny::Int)
    n1=nx+1; n2=ny+1
    np=n1*n2; nt=2*nx*ny
    msh=Mesh(4)
    szPoints!(msh,np)
    rx=range(R.xmin,stop=R.xmax,length=n1)
    ry=range(R.ymin,stop=R.ymax,length=n2)
    for (x,y) in Iterators.product(rx,ry)
        addPoint!(msh,x,y)
    end
    lin(i,j)=n1*(j-1)+i
    for (i,j) in Iterators.product(1:nx,1:ny)
        r=[lin(i,j) lin(i+1,j) lin(i+1,j+1) lin(i,j+1)]
        addTriangle!(msh,r[1],r[2],r[4])
        addTriangle!(msh,r[2],r[3],r[4])
        if j==1 addBoundaryEdge!(msh,1,r[1],r[2]) end
        if i==nx addBoundaryEdge!(msh,2,r[2],r[3]) end
        if j==ny addBoundaryEdge!(msh,3,r[3],r[4]) end
        if i==1 addBoundaryEdge!(msh,4,r[4],r[1]) end
    end
    msh.boundary[3]=reverse(msh.boundary[3],dims=(1))
    msh.boundary[4]=reverse(msh.boundary[4],dims=(1))
    return msh
end

function RefTriangle(n::Int)::Mesh
    msh=Mesh(3)
    szPoints!(msh,3)
    addPoint!(msh,0,0)
    addPoint!(msh,1,0)
    addPoint!(msh,0,1)
    addTriangle!(msh,1,2,3)
    addBoundaryEdge!(msh,1,1,2)
    addBoundaryEdge!(msh,2,2,3)
    addBoundaryEdge!(msh,3,3,1)
    LagrangeMesh!(msh,n)
    return msh
end

function tlin(i::Int,j::Int,n::Int)::Int
    if j==1
        if i==1
            return 1
        elseif i==n+1
            return 2
        else
            return i+2
        end
    elseif i==1
        if j==n+1
            return 3
        else
            return 3*n+2-j
        end
    elseif i+j==n+2
        return n+j+1
    else
        return 3*n+i-1+(j-2)*(n-1)-div((j-1)*(j-2),2)
    end
end

function Refine(msh::Mesh,n::Int)
    tl(i,j)=tlin(i,j,n)
    tripat=zeros(Int,0,3)
    for j=1:n
        for i=1:n-j+1
            tripat = [tripat; [tl(i,j) tl(i+1,j) tl(i,j+1)] ]
            if i<n-j+1
                tripat = [tripat; [tl(i+1,j) tl(i+1,j+1) tl(i,j+1)] ]
            end
        end
    end
    msh=LagrangeMesh(msh,n)
    nt=size(msh.tri,1)
    ntp=size(tripat,1)
    tri=zeros(Int,ntp*nt,3)
    ir=1
    for i=1:nt
        for j=1:ntp
            tri[ir,:] = msh.tri[i,tripat[j,:]]
            ir += 1
        end
    end
    msh.tri=tri
    for i=1:length(msh.boundary)
        bdry=zeros(Int,0,2)
        for j=1:size(msh.boundary[i],1)
            c=msh.boundary[i][j,3:end]
            a=[msh.boundary[i][j,1]; c]
            b=[c; msh.boundary[i][j,2]]
            bdry = [bdry; a b]
        end
        msh.boundary[i]=bdry
    end
    return msh
end

function CircularSector(n::Int)
    msh=Mesh(3)
    addPoint!(msh,0,0)
    addPoint!(msh,1,0)
    addPoint!(msh,1/2,sqrt(3)/2)
    addTriangle!(msh,1,2,3)
    addBoundaryEdge!(msh,1,1,2)
    addBoundaryEdge!(msh,2,2,3)
    addBoundaryEdge!(msh,3,3,1)
    for k=1:n
        msh=Refine(msh,2)
        for l=1:2:size(msh.boundary[2],1)
            i=msh.boundary[2][l,2]
            r=norm([msh.x[i],msh.y[i]])
            msh.x[i] /= r; msh.y[i] /= r
        end
    end
    return msh
end
