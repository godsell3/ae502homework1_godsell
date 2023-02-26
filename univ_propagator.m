% (Q1)
% A Universal Variable two-body orbit propagator.
% INPUTS
%  r0    - starting position position vector
%  v0    - starting velocity vector
%  tmax  - time to propagate to, assuming start at t=0
%  steps - number of steps to take from t0 to tmax
%  mu    - gravitational parameter
% OUTPUTS
%  trajectory - array of size (steps x 6) containing state vectors
%  tvals      - array of size (steps) containing time values
%
% USES findFinal_rv
function[trajectory, tvals] = univ_propagator(r0, v0, tmax, steps, mu)

dt = tmax/steps;
i = 1;
trajectory = zeros([steps, 6]);
trajectory(i,:) = [r0 v0];
tvals = linspace(0,tmax,steps);

for i=1:steps
    [r0, v0] = findFinal_rv(r0, v0, dt, mu);

    trajectory(i,:) = [r0 v0];
end

end