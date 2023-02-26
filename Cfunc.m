% (Q1)
% Implements Equation (3.50) in Curtis. Solves C(z) circular and hyperbolic
% trig functions.
% INPUTS
%  z - usually defined as alpha*chi^2
% OUTPUTS
%  C
function[C] = Cfunc(z)
if z>0
    C = (1/z)*(1-cos(sqrt(z)));
elseif z<0
    C = (-1/z)*(cosh(sqrt(-z))-1);
else %z=0
    C = 1/2;
end
end