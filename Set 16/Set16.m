% Matlab Set 16 - Toroid Field Plotting
% Matthew Jarzynowski

clc; % Clear the command window
clear; % Clear all previous variables

% Toroid Definition

I = 5.0; % Current inside the loops
N = 200; % Number of turns
Ri = 1.5; % Interior radii
Ro = 2.5; % Outer radii

% Plotting Points, (X,Y)
X_points = 50;
Y_points = 50;

% Regional Definition

Xmin = -4;
Xmax = 4;
Ymin = -4;
Ymax = 4;

% Step Sizes, Relative
dx = (Xmax - Xmin)/(X_points-1);
dy = (Ymax - Ymin)/(Y_points-1);

% Initial Grid Sizing
[X,Y] = meshgrid(Xmin:dx:Xmax, Ymin:dy:Ymax);

Z = zeros(size(X)); % Zeros matrix, waiting for plot values

% Magnetic Field Component Calculations
Bx = zeros(size(X)); % Relative X components
By = zeros(size(Y)); % Relative Y components

% Iterating through each point in the XY-plane using matrix dimensions
for i = 1:size(X,1)
    for j = 1:size(Y,2)
        x = X(i,j);
        y = Y(i,j);
        R = sqrt(x^2 + y^2);

        % Magnetic field, relative to the toroid's interior
        if R >= Ri && R <= Ro
            BPhi = (I*N)/(2*pi*R);
            Bx(i,j) = -BPhi * sin(atan2(y,x));
            By(i,j) = BPhi * cos(atan2(y,x));
        end
    end
end

% Plot the magnetic field, on the XY plane
quiver(X, Y, Bx, By);

% Plot Settings
xlabel('X (m)'); % Label x axis
ylabel('Y (m)', "Rotation",0); % Label y axis

title({'Toroid Magnetic Field Plot'})



