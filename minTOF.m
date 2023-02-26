% (Q3) (Q4)
% Calculates the minimum possible time of flight for a given r1 and r2.
% INPUTS
%  r1_    - initial position vector
%  r2_    - final position vector
%  theta_ - angle between the two position vectors (radians)
%  mu     - gravitational parameter
% OUTPUTS
%  tp - minimum TOF
function[tp] = minTOF(r1_, r2_, theta_, mu)
r1 = norm(r1_);
r2 = norm(r2_);

%calculate chord
c = sqrt(r1^2 + r2^2 - 2*r1*r2*cos(theta_));
s = (1/2)*(r1+r2+c);

tp = sqrt(2)/3/sqrt(mu)*(s^(3/2) - sign(sin(theta_))*(s-c)^(3/2));
end