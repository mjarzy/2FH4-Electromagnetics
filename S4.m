% Matlab Set 4 - Disc-Based Electric Fields
% Matthew Jarzynowski

clc; % Clear the command line
clear; % Clear previous variables

Eo = 8.854e-12; % Premetivity constant, with respect to air
S = 2e-6; % Disc-surface charge density

P = [0 0 1]; % Observation point
E = [0 0 0]; % Initial electric field, 0 in all components

% Step sizes for relavent "integrals"
rho_steps = 5000;
phi_steps = 5000;

% Defining our bounds.
rho_L = 0;
rho_U = 1;

phi_L = 0;
phi_U = (2*pi);

% Relavent infinitesimally small dimension
d_rho = (rho_U - rho_L)/rho_steps;
d_phi = (phi_U - phi_L)/phi_steps;

ds = d_rho * d_phi; % Relative area of a single element
dQ = S * ds; % The charge on a single element

% Double integration, using "for loops"
for j=1:rho_steps
    for i=1:phi_steps
        
        rho = rho_L + d_rho/2+(i-1)*d_rho; % The rho component of an element
        phi = phi_L + d_phi/2+(j-1)*d_phi; % The phi component of an element
        
        % Direction vector to observation point
        R = P - [rho*cos(phi) rho*sin(phi) 0]; 
        
        Rm = norm(R); % Direction vectors magnitude
        
        % Relative electric field contribution
        E = E + (rho*dQ/(4 * Eo * pi * Rm^3))*R;
    end
end

E




