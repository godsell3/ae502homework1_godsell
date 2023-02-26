% (Q2)
% Implements Equation (3.40) in Curtis.
% INPUTS
%  z     - usually defined as alpha*chi^2
%  r1mag - magnitude of vector r1
%  r2mag - magnitude of vector r2
%  A     - solution of Curtis eqn (5.35)
%  dt    - delta t
%  mu    - gravitational parameter
% OUTPUTS
%  fval
%
% USES Sfunc, Cfunc, yfunc
function[fval] = Ffunc(z, r1mag, r2mag, A, dt, mu)
S = Sfunc(z);
C = Cfunc(z);
y = yfunc(r1mag, r2mag, A, z);
fval = (y/C)^(3/2)*S + A*sqrt(y) - sqrt(mu)*dt; %eqn 5.40
end