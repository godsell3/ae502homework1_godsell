% (Q1)
% Implements Curtis Algorithm 3.3 to solve for universal anomaly.
% INPUTS
%  r0_  - position vector at t0
%  vr0_ - radial velocity at t0
%  a_   - alpha, reciprocal of semimajor axis
%  dt_  - delta t, the change in time
%  mu   - gravitational parameter
% OUTPUTS
%  xresult - calculated universal anomaly
%
% USES Cfunc, Sfunc
function[xresult] = univAnom_solver(r0_, vr0_, a_, dt_, mu)
%initial estimate
x_i = sqrt(mu)*abs(a_)*dt_;

%define inputs
ratio_i = 100; %arbitrary; just needs to be >=tol
tol = 1E-8; %choose tolerance
firsttime = true;

%calculate x
while (abs(ratio_i)>=tol)
    if ~firsttime
        x_i = x_i - ratio_i;
    end

    z_i = a_*x_i^2;
    C = Cfunc(z_i);
    S = Sfunc(z_i);

    f = (r0_*vr0_/sqrt(mu))*x_i^2*C + (1-a_*r0_)*x_i^3*S + r0_*x_i - sqrt(mu)*dt_;
    fprime = (r0_*vr0_/sqrt(mu))*x_i*(1-a_*x_i^2*S) + (1-a_*r0_)*x_i^2*C + r0_;

    ratio_i = f/fprime;
    firsttime = false;
end

xresult = x_i;
end