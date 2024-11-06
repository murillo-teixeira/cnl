clear;


% Time span for simulation
tspan = [0, 200];

% Initial conditions (multiple initial conditions for phase portrait)
initial_conditions = [
    100, 10, 0;   % G = 100, I = 10, X = 0
    120, 15, 0;   % G = 120, I = 15, X = 0
    80, 5, 0;     % G = 80, I = 5, X = 0
    110, 8, 0     % G = 110, I = 8, X = 0
];

% Plot phase portrait
figure;
hold on;
for i = 1:size(initial_conditions, 1)
    [t, y] = ode45(@bergman_model, tspan, initial_conditions(i, :));
    plot(y(:, 1), y(:, 2), 'LineWidth', 1.5); % G vs. I
end

xlabel('Glucose (G)');
ylabel('Insulin (I)');
title('Phase Portrait of Bergman Minimal Model');
legend('IC: [100, 10, 0]', 'IC: [120, 15, 0]', 'IC: [80, 5, 0]', 'IC: [110, 8, 0]');
grid on;
hold off;


% Define the Bergman minimal model equations
function dydt = bergman_model(t, y)
    % Parameters
    p1 = 0.0002;  % rate constant for glucose production
    p2 = 0.002;  % rate constant for insulin action decay
    p3 = 0.0001;  % rate constant for insulin effect on glucose uptake
    n = 0.001;    % insulin degradation rate
    gamma = 0.005;% rate of insulin secretion
    Gb = 0.1;    % basal glucose concentration


    G = y(1);
    I = y(2);
    X = y(3);
    
    dGdt = -(X + p1) * G + p1 * Gb;
    dIdt = -n * I + gamma * (G - Gb);
    dXdt = -p2 * X + p3 * I;
    
    dydt = [dGdt; dIdt; dXdt];
end