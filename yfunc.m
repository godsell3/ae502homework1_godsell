% (Q2)
% Implements Equation (3.38) in Curtis.
% INPUTS
%  r1mag - magnitude of vector r1
%  r2mag - magnitude of vector r2
%  A     - solution of Curtis eqn (5.35)
%  z     - usually defined as alpha*chi^2
% OUTPUTS
%  yval
%
% USES Sfunc, Cfunc
function[yval] = yfunc(r1mag, r2mag, A, z)
S = Sfunc(z);
C = Cfunc(z);
yval = r1mag + r2mag + A*(z*S - 1)/sqrt(C); %eqn 5.38
end