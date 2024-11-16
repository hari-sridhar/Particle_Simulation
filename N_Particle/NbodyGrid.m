% Define the matrix size
rows = input("Enter the number of rows:");
columns = input("Enter the number of columns: ");
N = rows*columns;

% Blood Paramters
mu = 0.0012; % Viscosity of Blood Plasma
m = 27e-12; % Mass of Particles
rho = 1025; % Density of Plasma
g = 9.81; % Gravitational constant

% Time parameters
t_start = 0;
t_end = 0.25;  % Total simulation time
dt = 0.01;    % Time step
time_steps = t_start:dt:t_end;
num_time_steps = length(time_steps);

% Create a Radius matrix
 
R = zeros(N, 1); % Initialize with zeros for an Nx1 vector
min_value = 3.5e-6;
max_value = 4.1e-6;
R = min_value + (max_value - min_value) * rand(N, 1); % Populate R with random values between 3.1e-6 and 4.1e-6

% Distance separation
separation = 2.1*max_value;

% Initialize position history matrices for N particles over time
X_history = zeros(N, num_time_steps);
Y_history = zeros(N, num_time_steps);

% Initialize velocity history matrices for N particles over time
Vx_history = zeros(N, num_time_steps);
Vy_history = zeros(N, num_time_steps);

% Initialize the unit vector matrices for x and y components
ex = zeros(N, N);  % Matrix to hold x-components of unit vectors
ey = zeros(N, N);  % Matrix to hold y-components of unit vectors

D = zeros(N);

% Initial positions with non-intersection check

particle_position = zeros(N,2);
count = 1; 

for i = 1:rows
    for j = 1:columns
        X_pos = (j - 1) * separation;
        Y_pos = (i - 1) * separation;
        particle_position(count, :) = [X_pos, Y_pos];
        count = count + 1;
    end
end


% Assign initial positions to the first column of X_history and Y_history
X_history(:, 1) = particle_position(:, 1);
Y_history(:, 1) = particle_position(:, 2); 


% Initialize KL as an NxN matrix
KL = zeros(N);

% Calculate each element KL(i, j) based on the formula
for i = 1:N
    for j = 1:N
        if i ~= j % Check if i is not equal to j
            KL(i, j) = (6 * pi * mu * (R(i) + R(j))^2 / 4);
        end
    end
end

%Create a Kd matrix (Drag coefficient) 
KD = 6*pi*mu*R;

% Create a Vol. Matrix
vol = (4/3)*pi*R.^3;

% Define random colors for each particle
colors = rand(N, 3); 

video_filename = 'particle_sim.mp4';
video_writer = VideoWriter(video_filename, 'MPEG-4'); % MPEG-4 format for better compatibility
video_writer.Quality = 100; % Set quality to maximum (optional)
open(video_writer);

% Set figure properties for larger frame size
figure('Position', [100, 100, 1540, 1000]); % Change width and height as desired (e.g., 800x600)
    
%% Simulation Loop
for i = 2:num_time_steps
% Update particle positions here (if they change over time)
    % X and Y represent current positions of particles, each as an Nx1 array        
    for p = 1:N
        for q = p+1:N  % Only calculate for p < q to avoid redundant calculations
            % Calculate unit vector components for particles p and q
            
            % Distance between particles p and q
            dx = X_history(q, i-1) - X_history(p, i-1);
            dy = Y_history(q,i-1) - Y_history(p,i-1);
            D(p,q) = sqrt(dx^2 + dy^2) - (R(p)+R(q)); 
            D(q,p) = D(p,q);
            
            % Calculate unit vector components if the distance is non-zero
            if D(p,q) ~= 0
                ex(p, q) = dx / D(p,q);
                ey(p, q) = dy / D(p,q);
                ex(q, p) = -ex(p, q);  % Symmetry: opposite direction
                ey(q, p) = -ey(p, q);
            end
            

        end
    end
    % Calculate P, Q, R, S matrices
    P = zeros(N); Q = zeros(N); S = zeros(N);
    for j = 1:N
        for k = 1:N
            if j ~= k
                P(j, k) = KL(j, k) * ex(j, k)^2 ./ D(j,k);
                P(k, j) = P(j, k);
                Q(j, k) = KL(j, k) * ex(j, k) * ey(j, k) ./ D(j,k);
                Q(k, j) = Q(j, k);
                S(j, k) = KL(j, k) * ey(j, k)^2 ./ D(j,k);
                S(k, j) = S(j, k);
            end
        end
    end

    for j = 1:N
        P(j, j) = -(KD(j) + sum(P(j, :)));
        Q(j, j) = -sum(Q(j, :));
        S(j, j) = -(KD(j) + sum(S(j, :)));
    end

    % Define matrix A
    A = [P, Q; Q, S];

    % Compute result matrix for B based on mass, density, volume, and gravity
    V = (m - rho * vol) * g;

    % Define B matrix
    T = zeros(N, 1);    % Example for upper part of B (change as needed)
    B = [T; V];

    % Solve for X
    velocities = linsolve(A, B);
    
    % Extract velocities from the solution
        Vx = velocities(1:N);
        Vy = velocities(N+1:end);
    
    for p = 1:N
        for q = p+1:N  % Only calculate for p < q to avoid redundant calculations
            % Calculate unit vector components for particles p and q
        
            % Distance between particles p and q
            dx = X_history(q, i-1) - X_history(p, i-1);
            dy = Y_history(q,i-1) - Y_history(p,i-1);
            D(p,q) = sqrt(dx^2 + dy^2) - (R(p)+R(q));
            D(q,p) = D(p,q);
        
            % Otherwise, use the regular velocity update for particle movement
            X_history(p, i) = X_history(p, i-1) + Vx(p) * dt;
            Y_history(p, i) = Y_history(p, i-1) + Vy(p) * dt;
            X_history(q, i) = X_history(q, i-1) + Vx(q) * dt;
            Y_history(q, i) = Y_history(q, i-1) + Vy(q) * dt;
        end
    end
       
    
    % Update positions based on velocities
    X_history(:, i) = X_history(:, i - 1) + Vx * dt;
    Y_history(:, i) = Y_history(:, i - 1) + Vy * dt;

    % Store the velocities in history matrices
    Vx_history(:, i) = Vx;  % Store the x-component of velocities
    Vy_history(:, i) = Vy;  % Store the y-component of velocities

    % Store the velocities in history matrices
    Vx_history(:, i) = Vx;  % Store the x-component of velocities
    Vy_history(:, i) = Vy;  % Store the y-component of velocities
    
end
     %% Animation code for N-particle simulation
% Assuming X_history and Y_history contain the positions of particles

% Define zoom factor for axis limits
zoomFactor = 1.5;  % You can adjust this factor for zoom level

for i = 1:length(time_steps)
    % Clear the current figure to update for the next time step
    clf;

    % Loop through all particles and plot each one
    for j = 1:N
        % Plot each particle as a circle with radius R(j) from history
        rectangle('Position', [X_history(j, i) - R(j), Y_history(j, i) - R(j), 2*R(j), 2*R(j)], ...
                  'Curvature', [1, 1], 'EdgeColor', colors(j, :), 'LineWidth', 1.5);
        hold on;
    end

    % Calculate the centroid of all particles at current time step
    centerX = mean(X_history(:, i));  % Mean of X positions
    centerY = mean(Y_history(:, i));  % Mean of Y positions

    % Calculate the distances of all particles from the centroid
    distances = sqrt((X_history(:, i) - centerX).^2 + (Y_history(:, i) - centerY).^2);

    % Determine the maximum distance for dynamic axis limits
    maxDist = max(distances) + max(R);  % Add the maximum radius for framing

    % Set axis limits dynamically around the centroid
    xlim([centerX - zoomFactor * maxDist, centerX + zoomFactor * maxDist]);
    ylim([centerY - zoomFactor * maxDist, centerY + zoomFactor * maxDist]);
    
    % Set the aspect ratio to be equal so particles appear round
    axis equal;

    % Capture the current figure and write to video
    frame = getframe(gcf); % Capture the current figure as a frame
    writeVideo(video_writer, frame); % Write the captured frame to the video

    % Update the plot
    drawnow;
end

close(video_writer);

% Create a figure for the trajectory plot
figure;
hold on;  % Enable overlaying multiple plots on the same figure

% Loop through all particles to plot their trajectory
for j = 1:N
    plot(X_history(j, :), Y_history(j, :), 'LineWidth', 2, 'DisplayName', sprintf('Particle %d', j), 'Color', colors(j, :));  % Trajectories of particles
end
% Customize the plot
legend();  % Display legend with names of particles
xlabel('X Position');  % Label for x-axis
ylabel('Y Position');  % Label for y-axis
title('Trajectories of Particles');  % Title of the plot
grid on;  % Turn on grid for better visibility
axis equal;  % Use equal scaling on both axes for correct aspect ratio

% Loop through all particles to plot their velocity
figure;
hold on;
for j = 1:N
    plot(time_steps, Vx_history(j, :), 'LineWidth', 2, 'DisplayName', sprintf('Particle Velocity%d', j), 'Color', colors(j, :));  % X - Velocity of particles
end
% Customize the plot
legend();  % Display legend with names of particles
xlabel('Time');  % Label for x-axis
ylabel('Y Position');  % Label for y-axis
title('X-Velocity of Particles');  % Title of the plot
grid on;  % Turn on grid for better visibility
axis equal;  % Use equal scaling on both axes for correct aspect ratio

figure;
hold on;  % Enable overlaying multiple plots on the same figure

% Loop through all particles to plot their velocity
for j = 1:N
    plot(time_steps, Vy_history(j, :), 'LineWidth', 2, 'DisplayName', sprintf('Particle Velocity%d', j), 'Color', colors(j, :));  % Y - Velocity of particles
end
% Customize the plot
legend();  % Display legend with names of particles
xlabel('X Position');  % Label for x-axis
ylabel('Y Position');  % Label for y-axis
title('Y-Velocity of Particles');  % Title of the plot
grid on;  % Turn on grid for better visibility
axis equal;  % Use equal scaling on both axes for correct aspect ratio

hold off;  % Release the plot hold
