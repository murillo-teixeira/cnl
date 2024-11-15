clear;

% -------------------------------- %
%         Model Parameters         %
% -------------------------------- %

l1_normal = 0.028;
l1_diabetes = 0.0;
l2 = 0.025;
l3 = 0.000013;
Vl = 12;
n = 5/54;
Gb = 4.5;
Ib = 15;

% -------------------------------- %
%         Model Equations          %
% -------------------------------- %

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

% -------------------------------- %
%          Time Response           %
% -------------------------------- %

tspan = [0 600];
z0 = [10; 0.01; 0];

% % Solve for normal case
[t_normal, z_normal] = ode45(odeSystemNormal, tspan, z0);

% % Solve for diabetic case
[t_diabetes, z_diabetes] = ode45(odeSystemDiabetes, tspan, z0);

% Create subplots
figure(1); % Use figure 1
set(gcf, 'Position', [100, 100, 400, 400]); % Set window size

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

% -------------------------------- %
%          Phase Portrait          %
% -------------------------------- %

% Quiver plot for phase portrait
figure(2); % Use figure 2
set(gcf, 'Position', [100, 100, 700, 250]); % Set window size

% Define grid for quiver plot
[G, I] = meshgrid(4.5:1:20, 15:1:30);

% Initialize derivatives
dG_normal = zeros(size(G));
dI_normal = zeros(size(I));
dG_diabetes = zeros(size(G));
dI_diabetes = zeros(size(I));

% Compute derivatives at each grid point
for i = 1:numel(G)
    z_normal_point = [G(i) - Gb; 0.01; I(i) - Ib];
    dz_normal = odeSystemNormal(0, z_normal_point);
    dG_normal(i) = dz_normal(1);
    dI_normal(i) = dz_normal(3);
    
    z_diabetes_point = [G(i) - Gb; 0.01; I(i) - Ib];
    dz_diabetes = odeSystemDiabetes(0, z_diabetes_point);
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

% Add a super title for the entire figure
sgtitle('Phase Portrait of the System');

% Save the second figure
exportgraphics(gcf, 'images/phase_portrait.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

% -------------------------------- %
%       Sliding Mode Control       %
% -------------------------------- %

% Define the sliding surface
% s = z(1) = G - Gb (glucose deviation from basal level)
% Our goal is to make s -> 0

% Control gains
k = 2;    % Proportional gain
eta = 12; % Reaching law parameter

% Modified control input for SMC with eta
u_smc = @(t, z) n * Vl * Ib + k * z(1) + eta * sign(z(1));

% Update the diabetic system with SMC
odeSystemDiabetesSMC = @(t, z) [
    -l1_diabetes*z(1) - z(1)*z(2) - Gb*z(2);         % z1_dot (Glucose dynamics)
    -l2*z(2) + l3*z(3);                              % z2_dot
    -n*z(3) - n*Ib + (u_smc(t, z)/Vl)                % z3_dot (Insulin dynamics with SMC)
    ];

[t_diabetes_smc, z_diabetes_smc] = ode45(odeSystemDiabetesSMC, tspan, z0);

figure(3);
set(gcf, 'Position', [100, 100, 400, 400]); % Set window size

% Top subplot: Insulin response
subplot(2, 1, 1);
plot(t_normal, z_normal(:,3) + Ib, 'k--', 'LineWidth', 1.5);
hold on;
plot(t_diabetes_smc, z_diabetes_smc(:,3) + Ib, 'g', 'LineWidth', 1.5);
ylim([0, 40]);
text(5, Ib - 0.07 * (max([z_normal(:,3) + Ib; z_diabetes(:,3) + Ib; z_diabetes_smc(:,3) + Ib]) + 1), 'Basal Insulin Level', 'FontSize', 8, 'Color', 'k'); % Annotation
xlabel('Time (minutes)');
ylabel('Insulin');
legend('Normal', 'Diabetes with SMC');
title('Insulin Response with Sliding Mode Control');
grid on;

% Bottom subplot: Glucose response
subplot(2, 1, 2);
plot(t_normal, z_normal(:,1) + Gb, 'b', 'LineWidth', 1.5); % Plot G + Gb for normal case
hold on;
plot(t_diabetes, z_diabetes(:,1) + Gb, 'r', 'LineWidth', 1.5); % Plot G + Gb for diabetes
plot(t_diabetes_smc, z_diabetes_smc(:,1) + Gb, 'g', 'LineWidth', 1.5); % Plot G + Gb for diabetes with SMC
plot(tspan, [Gb Gb], 'k--', 'LineWidth', 1.5); % Basal glucose level
ylim([0, max([z_normal(:,1) + Gb; z_diabetes(:,1) + Gb; z_diabetes_smc(:,1) + Gb]) + 1]); % Set ylim range
text(5, Gb - 0.07 * (max([z_normal(:,1) + Gb; z_diabetes(:,1) + Gb; z_diabetes_smc(:,1) + Gb]) + 1), 'Basal Glucose Level', 'FontSize', 8, 'Color', 'k'); % Annotation
xlabel('Time (minutes)');
ylabel('Glucose');
legend('Normal', 'Diabetes', 'Diabetes with SMC');
title('Glucose Response with Sliding Mode Control');
grid on;

% Save the figure
exportgraphics(gcf, 'images/insulin_glucose_response_smc.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');