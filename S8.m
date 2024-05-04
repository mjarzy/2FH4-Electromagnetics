% Matlab Set 8 - Spherical "Fluxial" Inquiry
% Matthew Jarzynowski

clc; % Clear the command bar
clear; % Remove all prior variables

% Bound Defintions
r_upper = 3;
r_lower = 2;
phi_upper = 2*pi;
phi_lower = 0;
theta_upper = pi;
theta_lower = 0;

% Discretization Steps
r_steps = 100;
phi_steps = 100;
theta_steps = 100;

% Differential Elements
dr = (r_upper - r_lower)/r_steps;
dphi = (phi_upper - phi_lower)/phi_steps;
dtheta = (theta_upper - theta_lower)/theta_steps;

WE = 0; % Initial energy stored

% Constants
Eo = 8.85e-12;
D = 2.0e-6;


% Calculating the Relative Energy Stored (J)
for j=1:theta_steps
    for k=1:phi_steps
        for i=1:r_steps
            r = r_lower + dr/2+(i-1)*dr; % R, for current volume
            theta = theta_lower + dtheta/2+(i-1)*dtheta; % Theta, for current volume
            phi = phi_lower + dphi/2+(j-1)*dphi; % Phi, for current volume

            Emag = D/(Eo*(r*r)); % Relative magnitude
            
            dV = (r*r)*sin(theta)*dtheta*dphi*dr; % Volume of current element
            dWE = (1/2)*Eo*(Emag*Emag)*dV; % Energy stored in current element
            
            WE = WE +dWE; % Sum relavent contribution
        end
    end
end

fprintf('The energy stored is approx. equal to, %f J (Joules)\n', WE);



clc;                % Clear the command window.
clear;              % Clear all variables.

vMid = -10;         % Set the potential in the middle to -10.
vSide = 0;          % Set the potential at the sides to 0.
steps = 21;         % Define the number of steps/grid points for the matrix.
points = steps^2;   % Total points in the matrix.

A = zeros(points, points); % Initialize matrix A with zeros.
b = zeros(points,1);       % Initialize vector b with zeros.

xMid = 11;                % Define the middle point in x-direction.
yBottom = 9;              % Define the bottom limit for the inner square.
yTop = 17;                % Define the top limit for the inner square.
PointIndex = 1;           % Initialize point index to 1.

% Loop over each point in the matrix.
for x = 1:steps
    for y = 1:steps
        % Check if the point is inside the inner square.
        if((y>=yBottom && y<= yTop) && (x == xMid))
            b(PointIndex, 1) = vMid; % Set the potential at the point to vMid.
        % Check if the point is at the boundaries of the matrix.
        elseif((y==1)||(y==steps)||(x==1)||(x==steps))
            b(PointIndex, 1) = vSide; % Set the potential at the point to vSide.
        else
            % For other points, update the matrix A to represent the relationship
            % between the current point and its neighbors in the grid.
            A(PointIndex, PointIndex) = -4; 
            A(PointIndex, PointIndex - 1) = 1;
            A(PointIndex, PointIndex + 1) = 1;
            A(PointIndex, PointIndex - steps) = 1;
            A(PointIndex, PointIndex + steps) = 1;
        end
        PointIndex = PointIndex + 1; % Increment the point index.
    end
end

V = A\b;                          % Solve the system of equations A*V = b for V.
V_square = reshape(V, steps, steps); % Reshape the vector V into a square matrix.

surf(V_square);                      % Create a surface plot of the potential.
figure;                             % Open a new figure window.
contourf(V_square);                  % Create a filled contour plot of the potential.
title('show V');                     % Add a title to the contour plot.
colorbar;                            % Display a colorbar alongside the contour plot.
set(gca, 'clim',[-10, 0]);           % Set the limits for color mapping.

figure;                              % Open another new figure window.
[Ex,Ey] = gradient(-V_square);       % Compute the gradient of the potential.
contour(V_square);                   % Create a contour plot of the potential.
hold on;                             % Hold the plot for adding more data.
quiver(-Ex,-Ey);                     % Plot the electric field as arrows.
hold off;                            % Release the plot hold.

















% Constants
clc;                % Clear the command window.
clear;              % Remove all variables from the workspace.
Vin = -5;           % Set the potential at the inner conductor.
vout = 5;           % Set the potential at the outer conductor.

% Setup
NumberOfYPoints = 21; % Number of points in the Y-direction.
delta = 0.15/20;      % Cell size, computed as total length divided by the number of intervals.
NumberOfXPoints = NumberOfYPoints; % Assuming a square grid, number of X points is the same as Y.
number_of_unknowns = NumberOfXPoints * NumberOfYPoints; % Total number of unknowns in the matrix.
b = zeros(number_of_unknowns,1);   % Initialize vector b with zeros.

% Vectors for the known potential indices on the grid
y1 = ((NumberOfXPoints-1)*0.4)+1;
y2 = ((NumberOfXPoints-1)*0.6)+1;
y1and = ((NumberOfXPoints-1)*(0.4))+1;
y2and = ((NumberOfXPoints-1)*(0.6))+1;
x1 = ((NumberOfYPoints-1)*0.4)+1;
x2 = ((NumberOfYPoints-1)*0.6)+1;

EquationCounter = 1; % Initialize equation counter.

% Main loop to populate the matrix A and vector b
for j = 1:NumberOfYPoints
    for i = 1:NumberOfXPoints
        % Apply boundary conditions to the points on the edges
        if((i == 1) || (i == NumberOfXPoints) || (j == NumberOfYPoints))
            A(EquationCounter, EquationCounter) = 1; % Set matrix A's diagonal to 1 for boundary points.
            b(EquationCounter, 1) = vout;            % Set the potential at the outer boundary.
        % Check if the current point is at the inner conductor boundaries
        elseif(((i > x1)&&(i<x2))&&((j>y1)&&(j<y2)))
            if((i==x1+1)&&(j>y1and)&&(j<y2and)) || ((i==x2-1)&&(j>y1and)&&(j<y2and)) 
                || ((j==y1+1)&&(i>x1and)&&(i<x2and)) || ((j==y2-1)&&(i>x1and)&&(i<x2and))
                A(EquationCounter, EquationCounter) = 1; % Set matrix A's diagonal to 1 for boundary points.
                b(EquationCounter, 1) = Vin;             % Set the potential at the inner boundary.
            else
                b(EquationCounter, 1) = vin;             % Set the potential at the inner conductor.
            end
        else
            % Set up the Laplace equation for interior points
            A(EquationCounter, EquationCounter) = -4;  % Set the central difference coefficient.
            A(EquationCounter, EquationCounter + 1) = 1;
            A(EquationCounter, EquationCounter - 1) = 1;
            A(EquationCounter, EquationCounter + NumberOfYPoints) = 1;
            A(EquationCounter, EquationCounter - NumberOfYPoints) = 1;
        end
        EquationCounter = EquationCounter + 1; % Increment the equation counter.
    end
end

% Solve the linear system and plot the results
V = A\b;                       % Solve the matrix equation A*V = b.
V_square = reshape(V, NumberOfXPoints, NumberOfYPoints); % Reshape V into a square matrix.
surf(V_square);                % Plot the surface of V.
figure;                        % Create a new figure.
[C,h] = contour(V_square);     % Draw contour lines of V.
set(h, 'ShowText','on', 'TextStep',get(h,'LevelStep')*2);
colormap cool;                 % Set the colormap.
figure;                        % Create another new figure.
[px,py] = gradient(-V_square); % Calculate the gradient of V.
contour(V_square);             % Draw contour lines of V.
hold on;                       % Hold the current plot.
quiver(px,py);                 % Draw a quiver plot of the gradient vectors.
hold off;                      % Release the current plot hold.
























