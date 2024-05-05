% Mutual Inductance, Between Two Circuits
% Matthew Jarzynowski

clc; % Clear the command bar
clear; % Remove all prior variables

% mutual_flux is afunction that calculates the mutual 
% inductance two circuits

% Circuit A, infinite size -------------------------------------

% I1, the current in the "infinite" circuit, Circuit A (A)

% a, the distance between conductors in the "infinite" circuit,
% not considered, as we can approximate the infinite circuit as a line (m)


% Circuit B, finite size ---------------------------------------

% b, distance between the "infinite" circuit, circuit A, and the
% secondary circuit, finite, Circuit B (m)

% h, the relative height of Circuit B (m)
% l, the relative length of Circuit B (m)

function M = mutual_flux(I1,a,b,h,l)
    
    % Relative permeability of free space
    u_0 = 4*pi*1e-7;
    
    % Magentic field, as a function of distance R
    B = @(r) u_0 * I1 ./ (2*pi*r); 
    
    % The integration of the magnetic field, over the height of Circuit B
    Phi = integral(@(r) B(r) * l, b, b + h); 

    % Calculate mutual inductance
    M = Phi / I1;

    % Display the calculated mutual inductance
    M;
end

% Consider these "base" parameters

I1 = 2; % Current, 2 (A)
a = 0.25; % Distance between conductors, Circuit A, 0.25 (m)
b = 0.50; % Distance between circuits, Circuit A - B, 0.50 (m)
h = 1.0; % Relative height of Circuit B, 1.0 (m)
l = 2.0; % Relative length of Circuit B, 2.0 (m)

% After calculating mutual inductance for two different currents
fprintf('------------------------------------------------------\n');
fprintf('Variable Current in the "Infinite" Circuits Conductors\n');
fprintf('------------------------------------------------------\n\n');

% Calculate mutual inductance for different currents and display results
M_current_2 = mutual_flux(2.0, a, b, h, l);
M_current_10 = mutual_flux(10.0, a, b, h, l);
fprintf('Mutual Inductance with I1 = 2 A: %0.3e Henries (H)\n', M_current_2);
fprintf('Mutual Inductance with I1 = 10 A: %0.3e Henries (H)\n\n', M_current_10);

% After calculating mutual inductance for two different sizes of Circuit B
fprintf('------------------------------------------------------\n');
fprintf('Variable Size of the Finite Circuit (B)\n');
fprintf('------------------------------------------------------\n\n');

% Calculate mutual inductance for different sizes of Circuit B and display results
M_size1 = mutual_flux(I1, a, b, 1.0, 2.0);
M_size2 = mutual_flux(I1, a, b, 1.5, 2.5);
fprintf('Mutual Inductance with Circuit B (h = 1.0 m, l = 2.0 m): %0.3e Henries (H)\n', M_size1);
fprintf('Mutual Inductance with Circuit B (h = 1.5 m, l = 2.5 m): %0.3e Henries (H)\n\n', M_size2);

% After calculating mutual inductance for two different distances
fprintf('------------------------------------------------------\n');
fprintf('Variable Distance Between Circuits (A - B)\n');
fprintf('------------------------------------------------------\n\n');

% Calculate mutual inductance for different distances and display results
M_dist_025 = mutual_flux(I1, a, 0.25, h, l);
M_dist_050 = mutual_flux(I1, a, 0.50, h, l);
fprintf('Mutual Inductance with distance b = 0.25 m: %0.3e Henries (H)\n', M_dist_025);
fprintf('Mutual Inductance with distance b = 0.50 m: %0.3e Henries (H)\n\n', M_dist_050);

% Conclusions
fprintf('------------------------------------------------------\n');
fprintf('Relavent Conclusions\n');
fprintf('------------------------------------------------------\n\n');
fprintf('The current in Circuit A, I1, does not matter, but, changing\n')
fprintf('"situtation geometry", in terms of a, b, h and l, changes\n')
fprintf('the overall "mutual inductance" relationship.\n')