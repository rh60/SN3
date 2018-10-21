export fem2

@inline function CalculateContributionOnElement(x, k, p, q, r, f)
    # Calculate length of element
    # Spočti délku elementu
    h = x[k+1] - x[k]
    # Define quadrature scheme
    # Specifikuj kvadraturní metodu
    M = 3
    QuadWeights = [5/18, 4/9, 5/18]
    QuadPoints = [0.5*(1-sqrt(0.6)), 0.5, 0.5*(1+sqrt(0.6))]
    # Initialise local contributions to A and b to zero
    A_local = zeros(2,2)
    b_local = zeros(2)
    # Loop over quadrature points
    # Přes všechny kvadraturní uzly
    for m = 1:M
        X = QuadPoints[m]
        # Evaluate local basis functions and their derivatives
        phi_local = [1-X, X]
        phi_prime_local = [-1, 1]
        # Evaluate p, q, r, f at quadrature points
        t = x[k]+h*X
        p_val = p(t)
        q_val = q(t)
        r_val = r(t)
        f_val = f(t)
        # Increment entries of A_local and b_local
        # with suitably weighted function evaluations
        for i=1:2
            for j=1:2
                A_local[i,j] +=
                h*QuadWeights[m]*(
                p_val/h/h*phi_prime_local[i]*
                phi_prime_local[j] +
                q_val/h*phi_local[i]*
                phi_prime_local[j] +
                r_val*phi_local[i]*phi_local[j])
            end
            b_local[i] += h*QuadWeights[m]*
            f_val*phi_local[i]
        end
    end
    return A_local, b_local
end

function fem2(x)
    # define p(x), q(x), r(x), f(x), u_a, eta_b
    # definuj p(x), q(x), r(x), f(x), u_a, eta_b
    p(x) = -x
    q(x) = -(2*x+3)
    r(x) = x+2
    f(x) = x^3*exp(2*x)
    u_a = exp(1)
    eta_b = -5*exp(1)*exp(1)
    # Deduce number of elements
    N = length(x)-1
    # Initialise A and b to zero
    A = zeros(N+1, N+1)
    b = zeros(N+1)
    # Loop over elements calculating local contributions and
    # inserting into linear system
    for k = 1:N
        A_local, b_local = CalculateContributionOnElement(x, k, p, q, r, f)
        A[k:k+1, k:k+1] += A_local
        b[k:k+1] += b_local
    end
    # Add contributions from Neumann boundary conditions into
    # linear system
    b[N+1] += eta_b
    # Set Dirichlet boundary conditions in row 1
    # Nastav Dirichletovu okrajovou podmínku v řádku 1
    A[1,:] .= 0
    A[1,1] = 1
    b[1] = u_a
    # Solve linear system
    # Vyřeš soustavu lineárních rovnic
    return A\b
end

using Plots

function test(n::Integer)
    x=linspace(0.5, 1, n)
    @time U=fem2(x)
    gr(dpi=300)
    plot(x,U,m=1,label="U(x)")
    e=exp(1)
    E=3*e+sqrt(e)/4
    K1=3*e-32/31*E
    K2=8/31*E
    u(x)=exp(x)*(K1+K2*x^3)+exp(2*x)*(x^2-2*x+2)
    plot!(x,u.(x),label="u(x)")
    png("data/o1")
    plot(x,U-u.(x),label="chyba MKP")
    png("data/o2")
end
