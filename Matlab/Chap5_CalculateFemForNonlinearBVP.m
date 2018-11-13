% Specify mesh x and initial guess U0
function U = Chap5_CalculateFemForNonlinearBVP(x, U0)
% Use Matlab routine to solve nonlinear algebraic system
U = fsolve(@CalculateResidual, U0, [], x);
% Plot solution
plot(x, U, '-')
xlabel('x')
ylabel('U')
% Function for calculating global residual
function Residual = CalculateResidual(U, x)
% Deduce number of elements
N = length(x) - 1;
% Initialise residual to zero
Residual = zeros(N+1, 1);
% Loop over elements, computing local contribution to
% residual and incrementing global residual
for k=1:N
    LocalResidual = CalculateResidualOnElement(U, x, k);
    Residual(k:k+1) = Residual(k:k+1) + LocalResidual;
end
% Handle Neumann boundary condition at x=2
Residual(N+1) = Residual(N+1) - 68;
% Handle Dirichlet boundary condition at x=-1
Residual(1) = U(1) - 1;
% Function for calculating local residual
function LocalResidual = CalculateResidualOnElement(U, x, k)
% Set element length h
h = x(k+1) - x(k);
% Define quadrature rule
M = 3;
QuadWeights = [5/18, 4/9, 5/18];
QuadPoints = [0.5*(1-sqrt(0.6)), 0.5, 0.5*(1+sqrt(0.6))];
% Initialise local residual to zero
LocalResidual = zeros(2,1);
% Loop over quadrature points to evaluate integration
% using Gaussian quadrature
for m=1:M
    % Set local coordinate X
    X = QuadPoints(m);
    % Calculate local values of basis functions,
    % their derivatives, x, U and the derivative of U
    phi_local = [1-X, X];
    phi_prime_local = [-1, 1];
    x_local = x(k) + h*X;
    LocalU = U(k)*phi_local(1) + U(k+1)*phi_local(2);
    LocalUPrime = U(k)*phi_prime_local(1) + ...
        U(k+1)*phi_prime_local(2);
    % Evaluate integrand and add into sum
    for i=1:2
        LocalResidual(i) = LocalResidual(i) + ...
            h*QuadWeights(m)* ...
            (1/h/h*(1+LocalU*LocalU)* ...
            LocalUPrime*phi_prime_local(i) + ...
            2*(1+x_local^4)*phi_local(i) + ...
            8*LocalU*LocalU*phi_local(i));
    end
end
