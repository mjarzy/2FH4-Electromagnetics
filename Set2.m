% Matlab Set 2 - Surface and Volume Integrals
% Matlab Set 2 - Volumetric Analysis

% Matthew Jarzynowski

clc;
clear;

% Initial volume, set to 0
volume = 0;

% Lower bounds of each integrand
rho_l = 0;
phi_l = (pi/4);
theta_l = (pi/4); 

% Defining the accuracy of the relative integral
% 1000, is a good approx, although 1000 should be
% used due to floating point arithmetic

rho_steps = 1000;
phi_steps = 1000;
theta_steps = 1000;

% Relative increments, for each integrand
d_rho = ((2-0)/rho_steps);
d_phi = ((pi/2 - pi/4)/phi_steps);
d_theta = ((pi/2 - pi/4)/phi_steps);

% For loop to calculate the volume of the enclosed
% surface. Order of the for loops doesnt matter
% since each calculation is independent.

for i=1:rho_steps
    for j=1:theta_steps
       for k=1:phi_steps
            
            % Adds volume contributions
            volume = volume + rho_l^2 * sin(theta_l) * d_rho * d_phi * d_theta;
        end

        % Increase theta each time phi is "traversed"
        theta_l = theta_l + d_theta;
    end

    %Reset theta to its lower bound, and increment rho
    % each time theta and phi have be "traversed"
    theta_l = (pi/4);
    rho_l = rho_l + d_rho; 
end

% Reset rho to its lower bound
rho_l = 0;

volume

% Surface Area Calculation

% The only surfaces, out of the 5 that we care about is S1
% the surface that is on the curved face, since surfaces S(2 - 5),
% are tangential to the "fluxial field array lines"

% Since we are calculating the volume of a sphere
% from its origin, there are 5 surfaces.

% Surface Area, each surface set to 0
surf1 = 0;
surf2 = 0;
surf3 = 0;
surf4 = 0;
surf5 = 0;

rho_u = 2; % Upper rho bound

% Calculating the surface area for surfaces (1 - 2), same
for i=1:rho_steps
    for j=1:theta_steps
        surf1 = surf1 + (rho_l * d_theta * d_rho);
    end

    % Increase rho each "traversal"
    rho_l = rho_l + d_rho;
end

rho_l = 0; % Reset to 0, lower bound
surf2 = surf1; % Consider them equal

rho_u = 2; % Define the upper bound of rho

% Calculating the surface area for surface 3, unique
for i=1:theta_steps
    for j=1:phi_steps
        surf3 = surf3 + (((rho_u)^2) * sin(theta_l) * d_theta * d_phi);
    end
    theta_l = theta_l + d_theta;
end

theta_l = (pi/4); % Redefine thetas lower bound
surf3;

theta_l = (pi/4);
theta_u = (pi/2);

% Calculating the surface area for surfaces (4 - 5), same
for i=1:rho_steps
    for j=1:phi_steps
        surf4 = surf4 + (sin(theta_u) * d_rho * d_phi);
        surf5 = surf5 + (sin(theta_l) * d_rho * d_phi);
    end
    rho_l = rho_l + d_rho;
end

surfaceArea = surf1 + surf2 + surf3 + surf4 +surf5; % Final sum

surfaceArea






