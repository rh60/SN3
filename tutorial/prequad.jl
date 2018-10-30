using FastGaussQuadrature

function prequad(n::Integer,m::Integer=n)
    rn=0:n
    ch2(j)=cos(j*π/n)
    nodes=sort(ch2.(rn))
    w=ones(n+1);w[1]=w[end]=1/2;w=w.*(-1).^rn
    dL=[zeros(size(nodes)) for i=1:n+1]
    for k=1:n+1
        for i=1:n+1
            if i!=k
                dL[k][i]=w[k]/w[i]/(nodes[i]-nodes[k])
            end
        end
        dL[k][k]=-sum(dL[k])
    end
    x,wg=gausslegendre(m)
    numer = zeros(m)
    dnumer = zeros(m)
    denom = zeros(m)
    for k=1:m        
        for j = 1:n+1
            xdiff = x.-nodes[j]
            temp = w[j]./xdiff
            if j==k
                numer = numer + temp
            end
            dnumer = dnumer + dL[k][j]*temp
            denom = denom + temp
        end
        lx = numer./denom
        dlx = dnumer./denom
    end
end

prequad(5)
