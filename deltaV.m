% (Q3) (Q4)
% Calculates the delta v between two velocity vectors for a rendezvous.
% INPUTS
%  v1   - initial terminal velocity vector
%  v2   - final terminal velocity vector
%  vdep - departure object velocity vector
%  varr - arrival object velocity vector
% OUTPUTS
%  deltav
function[deltav] = deltaV(v1, v2, vdep, varr)
dv1 = norm(v1 - vdep);
dv2 = norm(varr - v2);
deltav = dv1 + dv2;
end