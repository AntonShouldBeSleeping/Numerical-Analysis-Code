Broyden's Method
function [x, iters] = broyden(x0, f, J, tol, maxiter)
% broyden: Broyden's method for solving systems of nonlinear equations.
%
% Inputs:
% x0 - Initial guess
% f - Function handle representing the system of equations
% J - Function handle representing the Jacobian matrix of f
% tol - Error tolerance
% maxiter - Maximum number of iterations allowed
%
% Outputs:
% x - The approximate solution of the system of equations
% iters - The number of iterations performed
%
% To run code: (for example)
% f = @(x) [x(1)^2 + x(2)^2 - 1; x(1) - x(2)];
% J = @(x) [2*x(1), 2*x(2); 1, -1];
% [x, iters] = broyden([1; 1], f, J, 1e-8, 100)
iters = 0; % Initialize iteration counter
err = 1; % Initialize error
B0 = J(x0); % Initial guess for Broyden's method
while (abs(err) > tol && iters < maxiter)
 % Compute the Newton step
 dx = B0 \ -f(x0);
 x1 = x0 + dx;
 
 % Compute the difference between function values
 y = f(x1) - f(x0);
 
12
 % Update the solution and the approximation of the Jacobian
 x0 = x1;
 B0 = B0 + ((y - B0*dx)*transpose(dx))./(transpose(dx)*dx);
 
 % Compute the error
 err = norm(dx, inf);
 
 % Increment iteration counter
 iters = iters + 1;
end
x = x0; % Final solution
