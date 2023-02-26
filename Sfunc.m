% (Q1)
% Implements Equation (3.49) in Curtis. Solves S(z) circular and hyperbolic
% trig functions.
% INPUTS
%  z - usually defined as alpha*chi^2
% OUTPUTS
%  S
function[S] = Sfunc(z)
if z>0
    S = (1/sqrt(z)^3)*(sqrt(z) - sin(sqrt(z)));
elseif z<0
    S = (1/sqrt(-z)^3)*(sinh(sqrt(-z)) - sqrt(-z));
else %z=0
    S = 1/6;
end
end