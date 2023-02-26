% (Q1)
% Implements Curtis Algorithm 3.4. Solves for final r and v given initial
% r, v, and a delta t.
% INPUTS
%  r0_ - position vector at t0
%  v0_ - velocity vector at t0
%  dt_ - delta t, the change in time
%  mu  - gravitational parameter
% OUTPUTS
%  rf - final position vector
%  vf - final velocity vector
%
% USES univAnom_solver, Sfunc, Cfunc
function[rf, vf] = findFinal_rv(r0_, v0_, dt_, mu)
r0mag = norm(r0_);
v0mag = norm(v0_);
vr0mag = dot(r0_, v0_)/r0mag;
alpha = 2/r0mag - v0mag^2/mu;

%find universal anomaly
x = univAnom_solver(r0mag, vr0mag, alpha, dt_, mu);

%define S and C
S = Sfunc(alpha*x^2);
C = Cfunc(alpha*x^2);

%obtain f and g
f = 1 - x^2/r0mag*C;
g = dt_ - x^3/sqrt(mu)*S;

%compute r
rf = f*r0_ + g*v0_;
rfmag = norm(rf);

%obtain fdot and gdot
fd = sqrt(mu)/(rfmag*r0mag)*(alpha*x^3*S - x);
gd = 1 - x^2/rfmag*C;

%compute v
vf = fd*r0_ + gd*v0_;

end