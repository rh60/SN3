using HDF5
using FastGaussQuadrature

export @msave

linspace(a,b,n=100)=[range(a,stop=b,length=n)...]
data=raw"../data/data.h5"

function load(name::String)
    file=h5open(data,"r")
    v=read(file,name)
    close(file)
    return v
end

macro msave(vs...)
quote
    file=h5open(data,"w")
    for v in $vs
        name=string(v)
        write(file,name,Core.eval(Main, v))
    end
    close(file)
end
end

function prequad(f,a,b,n=10)
    x,w=gausslegendre(n)
    x=a.+(x.+1)*(b-a)/2
    x, f.(x).*w*(b-a)/2
end
