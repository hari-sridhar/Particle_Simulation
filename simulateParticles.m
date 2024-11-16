function simulateParticles(R1, R2, R3, X1, Y1, X2, Y2, X3, Y3)
    % Parameters (to be defined as per your system)
    mu = 0.0012; % Viscosity of Blood Plasma
    Kl12 = 6 * pi * mu * ((R1 + R2) / 2) ^ 2;  % Lubrication constant
    Kl13 = 6 * pi * mu * ((R1 + R3) / 2) ^ 2;  % Lubrication constant
    Kl23 = 6 * pi * mu * ((R3 + R2) / 2) ^ 2;  % Lubrication constant
    Kl21=Kl12;  % Lubrication constant
    Kl31=Kl13;  % Lubrication constant
    Kl32=Kl23;  % Lubrication constant
    Kd1 = 6 * pi * mu * R1;  % Drag coefficient for particle 1
    Kd2 = 6 * pi * mu * R2;  % Drag coefficient for particle 2
    Kd3 = 6 * pi * mu * R3;  % Drag coefficient for particle 3
    vol1 = (4 / 3) * pi * R1 ^ 3; % Volume of Particle 1
    vol2 = (4 / 3) * pi * R2 ^ 3; % Volume of Particle 2
    vol3 = (4 / 3) * pi * R3 ^ 3; % Volume of Particle 3
    m = 27e-12; % Mass of Particles
    rho = 1025; % Density of Plasma
    g = 9.81; % Gravitational constant

    % Time parameters
    t_start = 0;
    t_end = 0.5;  % Total simulation time
    dt = 0.0001;    % Time step
    time_steps = t_start:dt:t_end;

    % Initialize arrays to store particle positions over time
    X1_history = zeros(size(time_steps));
    Y1_history = zeros(size(time_steps));
    X2_history = zeros(size(time_steps));
    Y2_history = zeros(size(time_steps));
    X3_history = zeros(size(time_steps));
    Y3_history = zeros(size(time_steps));

    Vx1_history = zeros(size(time_steps));
    Vy1_history = zeros(size(time_steps));
    Vx2_history = zeros(size(time_steps));
    Vy2_history = zeros(size(time_steps));
    Vx3_history = zeros(size(time_steps));
    Vy3_history = zeros(size(time_steps));

    % Initialize arrays to store distance histories
    d12_history = zeros(size(time_steps)); % Distance between particle 1 and 2
    d13_history = zeros(size(time_steps)); % Distance between particle 1 and 3
    d23_history = zeros(size(time_steps)); % Distance between particle 2 and 3
    
    % Assign initial positions
    X1_history(1) = X1;
    Y1_history(1) = Y1;
    X2_history(1) = X2;
    Y2_history(1) = Y2;
    X3_history(1) = X3;
    Y3_history(1) = Y3;
    
video_filename = 'particle_sim.mp4';
video_writer = VideoWriter(video_filename, 'MPEG-4'); % MPEG-4 format for better compatibility
video_writer.Quality = 100; % Set quality to maximum (optional)
open(video_writer);


    
    % Simulation loop
    zoomFactor = 2; % Adjust this factor to set the zoom level
    for i = 2:length(time_steps)
        
        % Calculate distances between particles
        d12 = sqrt((X2-X1).^2 + (Y2-Y1).^2) - (R1+R2);
        d13 = sqrt((X3-X1).^2 + (Y3-Y1).^2)- (R1+R3);
        d23 = sqrt((X3-X2).^2 + (Y3-Y2).^2) - (R3+R2);
        d21 = d12;
        d31 = d13;
        d32 = d23;

    % Calculate unit vector components (assuming current positions of the particles)
    % Particle 1-2 interaction
        e12x = (X2 - X1)/d12;  
        e12y = (Y2 - Y1)/d12;  
        e21x = -e12x;  % Opposite direction
        e21y = -e12y;
    
    % Particle 1-3 interaction
        e13x = (X3 - X1)/d13;  
        e13y = (Y3 - Y1)/d13;  
        e31x = -e13x;  % Opposite direction
        e31y = -e13y;
    
    % Particle 2-3 interaction
        e23x = (X3 - X2)/d23;  
        e23y = (Y3 - Y2)/d23;  
        e32x = -e23x;  % Opposite direction
        e32y = -e23y;

     % Optional: Store distances for later analysis
        d12_history(i) = d12;
        d13_history(i) = d13;
        d23_history(i) = d23;
    
    % Define matrix A (dependent on particle positions and velocities)
        A = [(-Kd1-(Kl12*(e12x).^2/d12)-(Kl13*(e13x).^2/d13)) ,-Kl12*(e12x)*(e21x)/d12 , -Kl13*(e13x)*(e31x)/d13  , ((-Kl12*(e12x)*(e12y)/d12) - (Kl13*(e13x)*(e13y)/d13)) , -Kl12*(e21y)*(e12x)/d12, -Kl13*(e13x)*(e31y)/d13 ;...
        -Kl21*(e12x)*(e21x)/d21, -Kd2-(Kl21*(e21x).^2/d21)-(Kl23*(e23x).^2/d23) ,  -Kl23*(e23x)*(e32x)/d23 , -Kl21*(e21x)*(e12y)/d21  , -(Kl21*(e21x)*(e21y)/d21) - (Kl23*(e23x)*(e23y)/d21) , -Kl32*(e23x)*(e32y)/d32  ;...
       -Kl13*(e13x)*(e31x)/d13  ,  -Kl23*(e23x)*(e32x)/d23, -Kd3-(Kl31*(e31x).^2/d21)-(Kl32*(e32x).^2/d32)  , -Kl13*(e13x)*(e31y)/d13 , -Kl32*(e23y)*(e32x)/d32,  (-Kl31*(e31x)*(e31y)/d31) - (Kl32*(e32x)*(e32y)/d32)  ;...
        ((-Kl12*(e12x)*(e12y)/d12) - (Kl13*(e13x)*(e13y)/d13)), -Kl21*(e21x)*(e12y)/d21,  -Kl13*(e13x)*(e31y)/d13 ,  (-Kd1-(Kl12*(e12y).^2/d12)-(Kl13*(e13y).^2/d13)),-Kl21*(e21y)*(e12y)/d21 , -Kl13*(e13y)*(e31y)/d13 ;...
        -Kl12*(e21y)*(e12x)/d12 , -(Kl21*(e21x)*(e21y)/d21) - (Kl23*(e23x)*(e23y)/d21)  , -Kl32*(e23y)*(e32x)/d32, -Kl21*(e21y)*(e12y)/d21, -Kd2-(Kl21*(e21y).^2/d21)-(Kl23*(e23y).^2/d23)  ,-Kl32*(e23y)*(e32y)/d32  ;...
         -Kl13*(e13x)*(e31y)/d13, -Kl32*(e23x)*(e32y)/d32  ,  (-Kl31*(e31x)*(e31y)/d31) - (Kl32*(e32x)*(e32y)/d32)  ,-Kl13*(e13y)*(e31y)/d13 ,  -Kl32*(e23y)*(e32y)/d32 , ( -Kd3 - (Kl31*(e31y).^2/d31) - (Kl32*(e32y).^2/d32)) ];
    
    % Define matrix B (depends on mass, density, volume, and gravity)
        B = [0 ; 0; 0 ; (m - rho * vol1) * g ; (m - rho * vol2) * g  ; (m - rho * vol3) * g];
    
    % Solve for X = [Vx1 ; Vx2 ; Vx3 ; Vy1 ; Vy2 ; Vy3]
        X = linsolve(A, B);
    
    % Extract velocities from the solution
        Vx1 = X(1);
        Vx2 = X(2);
        Vx3 = X(3);
        Vy1 = X(4);
        Vy2 = X(5);
        Vy3 = X(6);
    
    % Update particle positions based on velocities
        X1 = X1 + Vx1 * dt;
        Y1 = Y1 + Vy1 * dt;
        X2 = X2 + Vx2 * dt;
        Y2 = Y2 + Vy2 * dt;
        X3 = X3 + Vx3 * dt;
        Y3 = Y3 + Vy3 * dt;
    
    % Store updated positions
        X1_history(i) = X1;
        Y1_history(i) = Y1;
        X2_history(i) = X2;
        Y2_history(i) = Y2;
        X3_history(i) = X3;
        Y3_history(i) = Y3;

     % Store updated Velocity
        Vx1_history(i) = -Vx1;
        Vy1_history(i) = -Vy1;
        Vx2_history(i) = -Vx2;
        Vy2_history(i) = -Vy2;
        Vx3_history(i) = -Vx3;
        Vy3_history(i) = -Vy3;
     

     % Combine individual histories into single matrices for X and Y positions
        X_history = [X1_history; X2_history; X3_history];
        Y_history = [Y1_history; Y2_history; Y3_history];
     
     % Combine individual histories into single matrices for X and Y positions
        D_history = [d12_history; d13_history; d23_history];

     % Combine individual histories into single matrices for X and Y positions
        Vx_history = [Vx1_history; Vx2_history; Vx3_history];
        Vy_history = [Vy1_history; Vy2_history; Vy3_history]; 

        % Animation code
        % Clear current frame
         clf;
        
        % Plot particles at their current positions with radii
        viscircles([X1_history(i), Y1_history(i)], R1, 'Color', 'r');
        hold all;
        viscircles([X2_history(i), Y2_history(i)], R2, 'Color', 'b');
        viscircles([X3_history(i), Y3_history(i)], R3, 'Color', 'g');
        
        axis equal;
        
        % Inside the loop, after updating positions
        % Calculate the centroid of all particles
        centerX = mean([X1, X2, X3]);
        centerY = mean([Y1, Y2, Y3]);

        % Determine the furthest particle distance from the centroid for consistent framing
        distances = sqrt(([X1, X2, X3] - centerX).^2 + ([Y1, Y2, Y3] - centerY).^2);
        maxDist = max(distances) + max([R1, R2, R3]);

        % Set axis limits dynamically around the centroid
        xlim([centerX - zoomFactor * maxDist, centerX + zoomFactor * maxDist]);
        ylim([centerY - zoomFactor * maxDist, centerY + zoomFactor * maxDist]);

        % Capture the current figure and write to video
    frame = getframe(gcf); % Capture the current figure as a frame
    writeVideo(video_writer, frame); % Write the captured frame to the video
        
        drawnow;
    end
close(video_writer);
    

% Plot distance histories after the simulation loop
    figure;
    hold on;  % Keep the current plot to overlay multiple lines
    plot(time_steps, d12_history, 'r', 'DisplayName', 'd12', 'LineWidth', 2);  % Distance between particle 1 and 2
    plot(time_steps, d13_history, 'g', 'DisplayName', 'd13', 'LineWidth', 2);  % Distance between particle 1 and 3
    plot(time_steps, d23_history, 'b', 'DisplayName', 'd23', 'LineWidth', 2);  % Distance between particle 2 and 3
    
    % Add the legend here
    legend('d12', 'd13', 'd23');  % Create a legend for the plot

    % Add labels and title
    xlabel('Time (s)');
    ylabel('Distance (m)');
    title('Distance Between Particles Over Time');
    grid on;  % Add a grid for better visibility

% Plot Velocity histories after the simulation loop
    figure;
    hold on;  % Keep the current plot to overlay multiple lines
    plot(time_steps, Vx1_history, 'r', 'DisplayName', 'Vx1', 'LineWidth', 2);  % Distance between particle 1 and 2
    plot(time_steps, Vx2_history, 'b', 'DisplayName', 'Vx2', 'LineWidth', 2);  % Distance between particle 1 and 3
    plot(time_steps, Vx3_history, 'g', 'DisplayName', 'Vx3', 'LineWidth', 2);  % Distance between particle 2 and 3
    
    % Add the legend here
    legend('Vx1', 'Vx2', 'Vx3');  % Create a legend for the plot

    % Add labels and title
    xlabel('Time (s)');
    ylabel('Velocity (m/s)');
    title('Particles Velocity Over Time');
    grid on;  % Add a grid for better visibility

    figure;
    hold on;  % Keep the current plot to overlay multiple lines
    plot(time_steps, Vy1_history, 'r', 'DisplayName', 'Vy1', 'LineWidth', 2);  % Distance between particle 1 and 2
    plot(time_steps, Vy2_history, 'b', 'DisplayName', 'Vy2', 'LineWidth', 2);  % Distance between particle 1 and 3
    plot(time_steps, Vy3_history, 'g', 'DisplayName', 'Vy3', 'LineWidth', 2);  % Distance between particle 2 and 3
    
    % Add the legend here
    legend('Vy1', 'Vy2', 'Vy3');  % Create a legend for the plot

    % Add labels and title
    xlabel('Time (s)');
    ylabel('Velocity (m/s)');
    title('Particles Velocity Over Time');
    grid on;  % Add a grid for better visibility

    % Create a figure for the trajectory plot
figure;
hold on;  % Enable overlaying multiple plots on the same figure

% Plot trajectory for each particle
plot(X1_history, Y1_history, 'r-', 'LineWidth', 2, 'DisplayName', 'Particle 1');  % Particle 1 trajectory in red
plot(X2_history, Y2_history, 'b-', 'LineWidth', 2, 'DisplayName', 'Particle 2');  % Particle 2 trajectory in blue
plot(X3_history, Y3_history, 'g-', 'LineWidth', 2, 'DisplayName', 'Particle 3');  % Particle 3 trajectory in green

% Customize the plot
legend();  % Display legend with names of particles
xlabel('X Position');  % Label for x-axis
ylabel('Y Position');  % Label for y-axis
title('Trajectories of Particles');  % Title of the plot
grid on;  % Turn on grid for better visibility
axis equal;  % Use equal scaling on both axes for correct aspect ratio

hold off;  % Release the plot hold


end


