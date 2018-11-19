using MATLAB

struct Rectangle
    ax::Float64
    bx::Float64
    ay::Float64
    by::Float64
end

struct Mesh
    x::Vector{Float64}
    y::Vector{Float64}
    tri::Matrix{Int}
    #n::Int
    function Mesh(R::Rectangle,nx::Int,ny::Int)
        rx=range(R.ax,stop=R.bx,length=nx+1)
        ry=range(R.ay,stop=R.by,length=ny+1)
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
                tri[j,:]=[i,i+1,i+nx+2]
                tri[j+1,:]=[i,i+nx+2,i+nx+1]
                i += 1
                j += 2
            end
            i += 1
        end
        new(ax,ay,tri)
    end
end

function testMesh()
    R=Rectangle(0,2,0,1)
    msh=Mesh(R,6,3)
    write_matfile("data/msh.mat",msh=msh)
    msh
end

testMesh()
