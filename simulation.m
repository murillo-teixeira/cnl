clear;
% Parameters
l1_normal = 0.028;
l1_diabetes = 0.0;
l2 = 0.025;
l3 = 0.000013;
Vl = 12;
n = 5/54;
Gb = 4.5;
Ib = 15;

u = @(t) n * Vl * Ib; % Basal input

odeSystemNormal = @(t, z) [
    -l1_normal*z(1) - z(1)*z(2) - Gb*z(2);  % z1_dot
    -l2*z(2) + l3*z(3);                    % z2_dot
    -n*z(3) - n*Ib + (u(t)/Vl)             % z3_dot
    ];

odeSystemDiabetes = @(t, z) [
    -l1_diabetes*z(1) - z(1)*z(2) - Gb*z(2);
    -l2*z(2) + l3*z(3);                       % z2_dot
    -n*z(3) - n*Ib + (u(t)/Vl)                % z3_dot
    ];

tspan = [0 300];

z0 = [10; 0.01; 0];

% Solve for normal case
[t_normal, z_normal] = ode45(odeSystemNormal, tspan, z0);

% Solve for diabetic case
[t_diabetes, z_diabetes] = ode45(odeSystemDiabetes, tspan, z0);

% Create subplots
figure('Position', [100, 100, 400, 400]); % Set window size

% Top subplot: Insulin response
subplot(2, 1, 1);
plot(t_normal, z_normal(:,3) + Ib, 'k--', 'LineWidth', 1.5); % Plot I + Ib in blue
hold on;
ylim([Ib-1, Ib+1]); % Set ylim range
text(5, Ib - 0.07 * (Ib+1 - (Ib-1)), 'Basal Insulin Level', 'FontSize', 8, 'Color', 'k'); % Annotation
xlabel('Time (minutes)');
ylabel('Insulin');
title('Insulin');
grid on;

% Bottom subplot: Glucose response
subplot(2, 1, 2);
plot(t_normal, z_normal(:,1) + Gb, 'b', 'LineWidth', 1.5); % Plot G + Gb
hold on;
plot(t_diabetes, z_diabetes(:,1) + Gb, 'r', 'LineWidth', 1.5); % Plot G + Gb for diabetes
plot(tspan, [Gb Gb], 'k--', 'LineWidth', 1.5); % Basal glucose level using plot
ylim([0, max(max(z_normal(:,1) + Gb), max(z_diabetes(:,1) + Gb)) + 1]); % Set ylim range
text(5, Gb - 0.07 * (max(max(z_normal(:,1) + Gb), max(z_diabetes(:,1) + Gb)) + 1), 'Basal Glucose Level', 'FontSize', 8, 'Color', 'k'); % Annotation
xlabel('Time (minutes)');
ylabel('Glucose');
legend('Normal', 'Diabetes');
title('Glucose Response: Normal vs Diabetes');
grid on;

% Save the first figure
exportgraphics(gcf, 'images/time_response.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

% Quiver plot for phase portrait
figure('Position', [100, 100, 700, 250]); % Set window size

% Define grid for quiver plot
[G, I] = meshgrid(4.5:1:20, 15:1:30);

% Initialize derivatives
dG_normal = zeros(size(G));
dI_normal = zeros(size(I));
dG_diabetes = zeros(size(G));
dI_diabetes = zeros(size(I));

% Compute derivatives at each grid point
for i = 1:numel(G)
    z_normal = [G(i) - Gb; 0.01; I(i) - Ib];
    dz_normal = odeSystemNormal(0, z_normal);
    dG_normal(i) = dz_normal(1);
    dI_normal(i) = dz_normal(3);
    
    z_diabetes = [G(i) - Gb; 0.01; I(i) - Ib];
    dz_diabetes = odeSystemDiabetes(0, z_diabetes);
    dG_diabetes(i) = dz_diabetes(1);
    dI_diabetes(i) = dz_diabetes(3);
end

% Plot quiver for normal case
subplot(1, 2, 1);
quiver(G, I, dG_normal, dI_normal, 'b');
hold on;

xlabel('Glucose');
ylabel('Insulin');
title('Normal');
grid on;

% Plot quiver for diabetic case
subplot(1, 2, 2);
quiver(G, I, dG_diabetes, dI_diabetes, 'r');
hold on;

xlabel('Glucose');
ylabel('Insulin');
title('Diabetes');
grid on;

% Save the second figure
exportgraphics(gcf, 'images/phase_portrait.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');