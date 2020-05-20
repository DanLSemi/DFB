% File name: muller.m
function out = muller(f, x0, x1, x2 , epsilon, max)
% Finds zeroes using Muller¡¯s method, good for complex roots.
% Variable description:
% out - result of search
% x1,x2,x3 - previous guesses
% epsilon - tolerance
% max - max number of iterations
%
y0 = f(x0);
y1 = f(x1);
y2 = f(x2);
iter = 0;
while (iter <= max)
iter = iter + 1;
a =( (x1 - x2)*(y0 - y2) - (x0 - x2)*(y1 - y2)) / ...
( (x0 - x2)*(x1 - x2)*(x0 - x1) );
%
b = ( ( x0 - x2 )^2 *( y1 - y2 ) - ( x1 - x2 )^2 *( y0 - y2 )) / ...
( (x0 - x2)*(x1 - x2)*(x0 - x1) );
%
c = y2;
%
if (a ~= 0)
disc = b*b - 4*a*c;
q1 = b + sqrt(disc);
q2 = b - sqrt(disc);
if (abs(q1) < abs(q2))
dx = - 2*c/q2;
else
dx = - 2*c/q1;
end
elseif (b ~= 0)
dx = - c/b;
end
x3 = x2 + dx;
x0 = x1;
x1 = x2;
x2 = x3;
%
y0 = y1;
y1 = y2;
y2 = f(x2);
%
if (abs(dx) < epsilon)
out = x2; break;
end
end