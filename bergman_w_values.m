% Parameters
l1 = 0;
l2 = 0.025;
l3 = 0.000013;
Vl = 12;
m = 5/54;
Gb = 4.5;
Ib = 15;

% Define the input function u(t)
u = @(t) 0; % Constant input, change to any other function as needed

% Define the system of ODEs
odeSystem = @(t, z) [
    -l1*z(1) - z(1)*z(2) - Gb*z(2);  % z1_dot
    -l2*z(2) + l3*z(3);              % z2_dot
    -m*z(3) - m*Ib + (u(t)/Vl)       % z3_dot
];

% Simulation parameters
tspan = [0, 200]; % Time range for simulation
initial_conditions = [ 
    22, 0, 15; % Initial conditions for z1, z2, z3
];

% Solve and plot time responses for z1, z2, and z3
figure;
for i = 1:size(initial_conditions, 1)
    [t, Z] = ode45(odeSystem, tspan, initial_conditions(i, :));
    
    % Plot z1, z2, z3 as functions of time
    subplot(3, 1, 1); % z1 vs time
    hold on;
    plot(t, Z(:, 1), 'DisplayName', ['Init Cond. ', num2str(i)]);
    ylabel('z1');
    title('Time response of z1');
    legend show;
    
    subplot(3, 1, 2); % z2 vs time
    hold on;
    plot(t, Z(:, 2), 'DisplayName', ['Init Cond. ', num2str(i)]);
    ylabel('z2');
    title('Time response of z2');
    legend show;

    subplot(3, 1, 3); % z3 vs time
    hold on;
    plot(t, Z(:, 3), 'DisplayName', ['Init Cond. ', num2str(i)]);
    ylabel('z3');
    xlabel('Time (s)');
    title('Time response of z3');
    legend show;
end

grid on;
