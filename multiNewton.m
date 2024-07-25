
% 
function [x, iters] = multiNewton(x0, f, J,tol, itersMax)

iters = 1;
incr = 1;
xg = x0;
while (abs(incr) > tol && iters < itersMax)
    
    dx = J(xg)\(-f(xg));
    xg = xg + dx;

    incr = norm(dx,inf);
    iters = iters + 1;
end

x = xg;
end
