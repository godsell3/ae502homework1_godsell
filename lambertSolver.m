% (Q2)
% Lambert's problem solver (Curtis Algorithm 5.2).
% INPUTS
%  r1   - initial position vector
%  r2   - target position vector
%  dt   - time interval delta t
%  type - either 'prograde' or 'retrograde'
%  mu   - gravitational parameter
% OUTPUTS
%  v1 - initial velocity vector
%  v2 - target position velocity vector
%
% USES Ffunc, yfunc, Sfunc, Cfunc
function[v1, v2] = lambertSolver(r1, r2, dt, type, mu)
%step 1 (eqn 5.24)
r1mag = norm(r1);
r2mag = norm(r2);

%step 2 (eqn 5.26)
dtheta = acos(dot(r1,r2)/r1mag/r2mag);
crossproduct = cross(r1, r2);
cross_z = crossproduct(3);
if strcmp(type, 'prograde')
    if cross_z<0
        dtheta = 2*pi - dtheta;
    end
elseif strcmp(type, 'retrograde')
    if cross_z>=0
        dtheta = 2*pi - dtheta;
    end
end

%step 3 (eqn 5.35)
A = sin(dtheta)*sqrt(r1mag*r2mag/(1-cos(dtheta)));

%step 4
%find good initial guess for z
z = -100;
while Ffunc(z,r1mag,r2mag,A,dt,mu)<0
    z = z + 0.1;
end
F = 1; %arbitrary; needs to be > 1
itr = 0;
while (abs(F) > 1e-10) && (itr < 5000)
    S = Sfunc(z);
    C = Cfunc(z);
    y = yfunc(r1mag, r2mag, A, z);
    F = Ffunc(z,r1mag,r2mag,A,dt,mu);
    y0 = r1mag + r2mag + A*(-1/sqrt(Cfunc(0)));
    if z==0
        Fp = sqrt(2)/40*y0^(3/2) + (A/8)*(sqrt(y0) + A*sqrt(1/2/y0));
    else
        Fp = (y/C)^(3/2)*(1/2/z*(C - 3/2*S/C) + 3/4*S^2/C) + A/8*(3*S/C*sqrt(y) + A*sqrt(C/y));
    end %eqn 5.43

    z1 = z - F/Fp;

    %update
    z = z1;
    itr = itr + 1;
end

%step 5
y = yfunc(r1mag, r2mag, A, z);

%step 6 (eqns 5.46)
f = 1 - y/r1mag;
g = A*sqrt(y/mu);
gd = 1 - y/r2mag;

%step 7
v1 = 1/g*(r2 - f*r1); %eqn 5.28
v2 = 1/g*(gd*r2 - r1); %eqn 5.29
end