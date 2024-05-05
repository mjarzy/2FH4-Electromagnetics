% Matlab Set 3 - Field-Based Charge Density
% Matthew Jarzynowski

clc % Clear command bar
clear % Remove prior variables

% Consider this, "variable defintion"

% Charges, Q1, Q2, in nC
Q1 = 8e-9;
Q2 = 8e-9;

pL = 4e-9; % Linear charge density
Eo = 8.8419e-12; % Permitivity of free space

P = [0 0 0]; % Observation point
A = [0 1 1]; % Q1 coordinates
B = [0 -1 1]; % Q2 coordinates

C = [3.5 3.5 0]; % Coordinates of the line charges centre, midpoint

stepL = 100000; % Step size of L

% Vector Manipulation

R1 = (P - A); % Vector from Q1 to observation point
R2 = (P - B); % Vector from Q2 to observation point

R1m = norm(R1); % R1 vector magnitude
R2m = norm(R2); % R2 vector magnitude

% Electric Field Calculation

E1 = Q1/(4*pi*Eo*R1m^3)*R1; % Field by Q1
E2 = Q2/(4*pi*Eo*R2m^3)*R2; % Field by Q2

D = norm(P - C); % Distance from observation to midpoint on line
L = sqrt(98)*D; % Length of the line, (m)

length = sqrt(98); % Relative length
dir_vec = [-7/sqrt(98) 7/sqrt(98) 0]; % Direction vector

dL = length/stepL; % Length of a segment
dL_Vector = dL*dir_vec; % A vector of a segment

% Consider changing to a vector pointing in the direction of the line

EL = [0 0 0]; % Initialize the field of the line segment

% Perpendicular centre of a segment
C_Start = C - length/2 * dir_vec;
C_Segment = C_Start;

for i=1: stepL
    
    % A vector from the observation point, P, to the centre
    % of the linearly charged line
    R = P - C_Segment;
    Rm = norm(R); % R vector magnitude

    EL = EL +dL * pL/(4*pi*Eo*(Rm)^3)*R; % Each segments contribution

    % Relative centre of the "i-th" segment
    C_Segment = C_Segment + dL_Vector; 
end

ET = E1 + E2 + EL; % Sum each electric field contribution

ET

