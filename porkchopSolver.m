% (Q3) (Q4)
% Generates the data to be plotted to create a pork chop plot.
% INPUTS
%  E_traj    - Earth's trajectory
%  E_tvals   - the time from t0 at each step in Earth's trajectory
%  obj_traj  - the trajectory of the input object
%  obj_tvals - the time from t0 at each step of the object's trajectory
%  type      - 'prograde' or 'retrograde'
%  mu        - gravitational parameter
% OUTPUTS
%  dv    - matrix of delta v values for rendezvous
%  dv_fb - matrix of delta v values for flyby
%
% USES lambertSolver, deltaV, deltaVfb
function[dv, dv_fb] = porkchopSolver(E_traj, E_tvals, obj_traj, obj_tvals, type, mu)

dv = ones([length(E_tvals), length(obj_tvals)]);
dv = dv*-1;
dv_fb = dv;

%iterate through launch window
for i=1:length(E_tvals)
    launch_time = E_tvals(i);
    launch_state = E_traj(i,:);

    %iterate through arrival window
    for j=1:length(obj_tvals)
        [i,j];
        arr_time = obj_tvals(j);
        arr_state = obj_traj(j,:);

        %check that we launch before we arrive
        if launch_time<arr_time
            TOF = arr_time - launch_time;

            %solve lambert's problem
            [v1,v2] = lambertSolver(launch_state(1:3),arr_state(1:3),TOF,type,mu);

            %check if v1, v2 are good
            if all([~isnan(v1), ~isnan(v2)])
                %add results to data
                dv(i,j) = deltaV(v1,v2,launch_state(4:6),arr_state(4:6));
                dv_fb(i,j) = deltaVfb(v1,launch_state(4:6));
            end
        end

    end
end

end