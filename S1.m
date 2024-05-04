% Matlab Set 1 - Vector Analysis

% Consider these vectors as position vectors
% from the origin, no need for an origin point.

R1 = [1, 2, 3];
R2 = [3, 2, 1];

% Dot product between R1 and R2
dotP = dot(R1,R2)

% Projection of R1 onto R2 
proj = (dotP/norm(R2)^2) * R2

% Angle between R1 and R2
angle = acos(dotP / (norm(R1)*norm(R2)))



