%% Newthon's Method for approximating a root for a given tolerance. 
% 
function [x,iters] = newton(x0,f,df,tol)
% Newton’s method for f(x) = 0.
%
% input: x0  - initial guess
%        f   - function f
%        df  - derivative of the function f
%        tol - error tolerance
%
% output: x - the approximate root f(x) = 0
%         iters - the number of iterations required
%
% To run code: (for example)
% f = @(x) (x/2)^2-sin(x);
% df = @(x) x/2-cos(x);
% [x,iters] = newton(2,f,df,1e-8)

xg = x0; 
incr = 1;
iters = 0; 
m=1;

while (abs(incr) > tol)
    incr = -m*f(xg)/df(xg);
    xg = xg + incr;
    
    iters = iters + 1;
end

% Final value
x = xg;

end

