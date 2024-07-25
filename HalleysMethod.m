#Halley's Method code created for an assignment in my numerical analysis module. 

function [x, iters] = Halley(x0, f, df, dff, tol)
% Halley's method for solving f(x) = 0
%
% Inputs:
% x0 - Initial guess
% f - Function f
% df - First derivative of the function f
% dff - Second derivative of the function f
% tol - Error tolerance
%
% Outputs:
% x - The approximate root f(x) = 0
% iters - The number of iterations required
%
% To run code: (for example)
% f = @(x) (x/2)^2 - sin(x);
% df = @(x) x/2 - cos(x);
% dff = @(x) 1/2 + sin(x)
% [x, iters] = Halley(2, f, df, dff, 1e-8)

xg = x0; % Initial guess
incr = 1; % Initial increment
iters = 0; % Initialize iteration counter
10
% Define a function to avoid excessive clutter within the while loop
dividor = @(x0) df(x0) - (f(x0) .* dff(x0)) ./ (2 .* df(x0)); 
while (abs(incr) > tol)
 % Compute the increment using Halley's method
 incr = -f(xg) ./ dividor(xg);
 
 % Update the guess
 xg = xg + incr;
 
 % Increment the iteration counter
 iters = iters + 1;
end
% Final value
x = xg;
end
