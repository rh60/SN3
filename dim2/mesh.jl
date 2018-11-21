using MATLAB

struct Rectangle
    xleft::Float64
    xright::Float64
    yleft::Float64
    yright::Float64
end

struct Mesh
    x::Vector{Float64}
    y::Vector{Float64}
    tri::Matrix{Int}
    xleft::Vector{Int}
    xright::Vector{Int}
    yleft::Vector{Int}
    yright::Vector{Int}
    #n::Int
    function Mesh(R::Rectangle,nx::Int,ny::Int)
        rx=range(R.xleft,stop=R.xright,length=nx+1)
        ry=range(R.xleft,stop=R.yright,length=ny+1)
        ax=zeros((nx+1)*(ny+1))
        ay=zeros((nx+1)*(ny+1))
        tri=zeros(Int,2*nx*ny,3)
        i=1
        for y in ry
            for x in rx
                    ax[i]=x
                    ay[i]=y
                    i+=1
            end
        end
        i=1
        j=1
        for iy=1:ny
            for ix=1:nx
                tri[j,:]=[i,i+1,i+nx+1]
                tri[j+1,:]=[i+1,i+nx+2,i+nx+1]
                i += 1
                j += 2
            end
            i += 1
        end
        xleft=collect(1:nx+1:(nx+1)*ny+1)
        xright=collect(nx+1:nx+1:(nx+1)*(ny+1))
        yleft=collect(1:nx+1)
        yright=collect((nx+1)*ny+1:(nx+1)*(ny+1))
        new(ax,ay,tri,xleft,xright,yleft,yright)
    end
end

function testMesh()
    R=Rectangle(0,2,0,1)
    msh=Mesh(R,10,5)
    write_matfile("data/msh.mat",msh=msh)
    mxcall(:addpath,0,raw"d:\documents\julia\Ming\dim2")
    mxcall(:show,0)
end

#testMesh()
