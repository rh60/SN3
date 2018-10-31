using FastGaussQuadrature

function prequad(n::Integer,m::Integer=2*n)
    rn=0:n
    ch2(j)=cos(j*π/n)
    nodes=sort(ch2.(rn))
    w=ones(n+1);w[1]=w[end]=1/2;w=w.*(-1).^rn
    dL=Vector{Vector{Float64}}(undef,n+1)
    for k=1:n+1
        dL[k]=zeros(n+1)
        for i=1:n+1
            if i!=k
                dL[k][i]=w[k]/w[i]/(nodes[i]-nodes[k])
            end
        end
        dL[k][k]=-sum(dL[k])
    end

    x,g=gausslegendre(m)

    Lx=Vector{Vector{Float64}}(undef,n+1)
    dLx=Vector{Vector{Float64}}(undef,n+1)
    for k=1:n+1
        numer = zeros(m)
        dnumer = zeros(m)
        denom = zeros(m)
        for j = 1:n+1
            xdiff = x.-nodes[j]
            temp = w[j]./xdiff
            if j==k
                numer = numer + temp
            end
            dnumer = dnumer + dL[k][j]*temp
            denom = denom + temp
        end
        Lx[k] = numer./denom
        dLx[k] = dnumer./denom
    end
    return (L=Lx,dL=dLx,x=x,w=g)
end

R=prequad(5)


