export fem1

function fem1(N::Integer)
# Initialise A and b to zero
    A = zeros(N+1, N+1)
    b = zeros(N+1)
# Generate N+1 nodes, equally spaced between 0 and 1
    x = linspace(0, 1, N+1)
# Loop over elements calculating local contributions
# and incrementing to the global linear system
    for k=1:N
        # Calculate element length
        h = 1/N
        # Calculate local contributions
        A_local = [1/h -1/h; -1/h 1/h]
        b_local = [h; h]
        # Increment global system
        A[k:k+1, k:k+1] += A_local
        b[k:k+1] += b_local
    end
# Set Dirichlet boundary conditions
    A[1,:] .= 0
    A[1,1] = 1
    b[1] = 0
    A[N+1,:] .= 0
    A[N+1,N+1] = 1
    b[N+1] = 0
# Solve linear system and plot FE solution
    U = A\b;
# Plot finite element solution
# plot(x, U, '-o')
# xlabel('x')
# ylabel('U')
    return (mesh=x,solution=U)
end
