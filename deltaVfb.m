% (Q3) (Q4)
% Calculates the delta v between two velocity vectors for a flyby.
% INPUTS
%  v1   - initial terminal velocity vector
%  vdep - departure object velocity vector
% OUTPUTS
%  deltav
function[deltav] = deltaVfb(v1, vdep)
deltav = norm(v1 - vdep);
end