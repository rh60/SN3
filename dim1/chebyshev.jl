struct Chebyshev
    nodes::Vector{Float64}
    weights::Vector{Float64}
    derivatives::Matrix{Float64}
    degree::Int
    function Chebyshev(n::Int)
        r=n:-1:0
        ch2(j)=cos(j*π/n)
        nodes=(ch2.(r).+1)/2
        weights=(-1.0).^r
        weights[1]/=2;weights[end]/=2
        derivatives=zeros(n+1,n+1)
        for k=1:n+1
            for i=1:n+1
                if i!=k
                    derivatives[k,i]=weights[k]/weights[i]/(nodes[i]-nodes[k])
                end
            end
        end
        for k=1:n+1
            derivatives[k,k]=-sum(derivatives[:,k])
        end
        new(nodes,weights,derivatives,n)
    end
end

function InterpolateCardinal(E::Chebyshev,k::Int,x::Vector{Float64})
    n = length(E.nodes)
    m = length(x)
    denom = zeros(m)
    numer = zeros(m)
    dnumer = zeros(m)
    exact = zeros(Int,0)
    exactv = zeros(0)
    exactd = zeros(0)
    for j = 1:n
        xdiff = x.-E.nodes[j]
        temp = E.weights[j]./xdiff
        dnumer = dnumer + E.derivatives[k,j]*temp
        if j==k
            numer = numer + temp
        end
        denom = denom + temp
        f = findfirst(xdiff.==0)
        if f != nothing
            push!(exact,f)
            if j==k
                push!(exactv,1)
            else
                push!(exactv,0)
            end
            push!(exactd,E.derivatives[k,j])
        end
    end
    y = numer./denom
    y[exact] = exactv
    dy = dnumer./denom
    dy[exact] = exactd
    return y, dy
end

function Interpolate(E::Chebyshev,values::Vector{Float64},interpolx::Vector{Float64})
    l = interpolx[end]-interpolx[1]
    x = (interpolx.-interpolx[1])/l
    n = length(E.nodes)
    m = length(x)
    denom = zeros(m)
    numer = zeros(m)
    exact = zeros(Int,0)
    exactv = zeros(0)
    for j = 1:n
        xdiff = x.-E.nodes[j]
        temp = E.weights[j]./xdiff
        numer = numer + values[j]*temp
        denom = denom + temp
        f = findfirst(xdiff.==0)
        if f != nothing
            push!(exact,f)
            push!(exactv,values[j])
        end
    end
    temp = numer./denom
    temp[exact] = exactv
    return temp
end

function Interpolate(E::Chebyshev,f,interpolx::Vector{Float64})
    l = interpolx[end]-interpolx[1]
    t = l*E.nodes.+interpolx[1]
    Interpolate(E,f.(t),x)
end

function Quadrature(E::Chebyshev,nq::Int)
    n=E.degree
    Q=Quad(nq)
    L=Vector{Vector{Float64}}(undef,n+1)
    dL=Vector{Vector{Float64}}(undef,n+1)
    for k=1:n+1
        L[k], dL[k]  = InterpolateCardinal(E,k,Q.Points)
    end
    return L, dL, Q
end
