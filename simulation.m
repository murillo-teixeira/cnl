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
    -n*z(3) - n*Ib + (u(t)/Vl);             % z3_dot
    0];

odeSystemDiabetes = @(t, z) [
    -l1_diabetes*z(1) - z(1)*z(2) - Gb*z(2);
    -l2*z(2) + l3*z(3);                       % z2_dot
    -n*z(3) - n*Ib + (u(t)/Vl);                % z3_dot
    0];

% -------------------------------- %
%          Time Response           %
% -------------------------------- %

tspan = [0 600];
z0 = [10; 0.01; 0; 0];

% % Solve for normal case
[t_normal, z_normal] = ode45(odeSystemNormal, tspan, z0);

% % Solve for diabetic case
[t_diabetes, z_diabetes] = ode45(odeSystemDiabetes, tspan, z0);

% Create subplots
figure(1); % Use figure 1
set(gcf, 'Position', [100, 100, 700, 500], 'Visible', 'off');
% Left subplot: Insulin response
subplot(2, 1, 1);
plot(t_normal, z_normal(:,3) + Ib, 'k', 'LineWidth', 2); % Plot I + Ib in blue
hold on;
set(gca, 'FontSize', 12); % Increase ticks font size
ylim([Ib-1, Ib+1]); % Set ylim range
text(5, Ib - 0.1 * (Ib+1 - (Ib-1)), 'Basal Insulin Level', 'FontSize', 14, 'Color', 'k'); % Annotation
xlabel('Time (minutes)', 'FontSize', 14);
ylabel('I(t) (mU/L)', 'FontSize', 14);
title('Insulin Concentration', 'FontSize', 16);
grid on;

% Right subplot: Glucose response
subplot(2, 1, 2);
plot(t_normal, z_normal(:,1) + Gb, 'b', 'LineWidth', 2); % Plot G + Gb
hold on;
set(gca, 'FontSize', 12); % Increase ticks font size
plot(t_diabetes, z_diabetes(:,1) + Gb, 'r', 'LineWidth', 2); % Plot G + Gb for diabetes
plot(tspan, [Gb Gb], 'k--', 'LineWidth', 2); % Basal glucose level using plot
ylim([0, max(max(z_normal(:,1) + Gb), max(z_diabetes(:,1) + Gb)) + 1]); % Set ylim range
text(5, Gb - 0.1 * (max(max(z_normal(:,1) + Gb), max(z_diabetes(:,1) + Gb)) + 1), 'Basal Glucose Level', 'FontSize', 14, 'Color', 'k'); % Annotation
xlabel('Time (minutes)', 'FontSize', 14);
ylabel('G(t) (mmol/L)', 'FontSize', 14);
% legend('Normal', 'Diabetes', 'FontSize', 12);
title('Glucose Response: Normal vs Diabetes', 'FontSize', 16);
grid on;

% Save the first figure
exportgraphics(gcf, 'images/time_response.png', 'Resolution', 300);

% -------------------------------- %
%          Phase Portrait          %
% -------------------------------- %

% % Quiver plot for phase portrait
% figure(2); % Use figure 2
% set(gcf, 'Position', [100, 100, 1000, 300], 'Visible', 'off');

% % Define grid for quiver plot
% [G, I] = meshgrid(4.5:1:20, 15:1:30);

% % Initialize derivatives
% dG_normal = zeros(size(G));
% dI_normal = zeros(size(I));
% dG_diabetes = zeros(size(G));
% dI_diabetes = zeros(size(I));

% % Compute derivatives at each grid point
% for i = 1:numel(G)
%     z_normal_point = [G(i) - Gb; 0.01; I(i) - Ib;0];
%     dz_normal = odeSystemNormal(0, z_normal_point);
%     dG_normal(i) = dz_normal(1);
%     dI_normal(i) = dz_normal(3);

%     z_diabetes_point = [G(i) - Gb; 0.01; I(i) - Ib;0];
%     dz_diabetes = odeSystemDiabetes(0, z_diabetes_point);
%     dG_diabetes(i) = dz_diabetes(1);
%     dI_diabetes(i) = dz_diabetes(3);
% end

% % Plot quiver for normal case
% subplot(1, 2, 1);
% quiver(G, I, dG_normal, dI_normal, 'b');
% hold on;
% xlabel('Glucose', 'FontSize', 14);
% ylabel('Insulin', 'FontSize', 14);
% title('Normal', 'FontSize', 16);
% set(gca, 'FontSize', 12); % Adjust tick font size
% grid on;

% % Plot quiver for diabetic case
% subplot(1, 2, 2);
% quiver(G, I, dG_diabetes, dI_diabetes, 'r');
% hold on;
% xlabel('Glucose', 'FontSize', 14);
% ylabel('Insulin', 'FontSize', 14);
% title('Diabetes', 'FontSize', 16);
% set(gca, 'FontSize', 12); % Adjust tick font size
% grid on;

% % Save the second figure
% exportgraphics(gcf, 'images/phase_portrait.png', 'Resolution', 300);

% ------------------------------
%       Sliding Mode Control with Chattering      %
% ------------------------------

k = 0.5; eta = 3;
a = -6; b = 1; m = n; phi = 45;

u_smc = @(t, z) n * Vl * Ib + (Vl/b)*(a*l1_diabetes*z(1) + a*z(1)*z(2) + a*Gb*z(2) + b*m*z(3) + b*m*Ib - eta*sign(a*z(1) + b*z(3)) - k*(a*z(1) + b*z(3)));

% Update the diabetic system with SMC
odeSystemDiabetesSMC = @(t, z) [
    -l1_diabetes*z(1) - z(1)*z(2) - Gb*z(2);         % z1_dot (Glucose dynamics)
    -l2*z(2) + l3*z(3);                              % z2_dot
    -n*z(3) - n*Ib + (u_smc(t, z)/Vl);                % z3_dot (Insulin dynamics with SMC)
    u_smc(t,z)];

[t_diabetes_smc, z_diabetes_smc] = ode45(odeSystemDiabetesSMC, tspan, z0);

figure(3);
set(gcf, 'Position', [100, 100, 700, 800]);
% Left subplot: Insulin response
subplot(3, 1, 2);
plot(t_normal, z_normal(:,3) + Ib, 'k--', 'LineWidth', 1.5);
hold on;
plot(t_diabetes_smc, z_diabetes_smc(:,3) + Ib, 'g', 'LineWidth', 1.5);
text(5, Ib - 7, ...
    'Basal Insulin Level', 'FontSize', 14, 'Color', 'k'); % Annotation
xlabel('Time (minutes)', 'FontSize', 14);
ylabel('Insulin', 'FontSize', 14);
ylim([0, 80]);
% legend('Normal', 'Diabetes with SMC', 'FontSize', 12);
title('Insulin Response', 'FontSize', 16);
set(gca, 'FontSize', 12); % Adjust tick font size
grid on;

% Middle subplot: Glucose response
subplot(3, 1, 3);
plot(t_normal, z_normal(:,1) + Gb, 'b', 'LineWidth', 1.5); % Plot G + Gb for normal case
hold on;
plot(t_diabetes, z_diabetes(:,1) + Gb, 'r', 'LineWidth', 1.5); % Plot G + Gb for diabetes
plot(t_diabetes_smc, z_diabetes_smc(:,1) + Gb, 'g', 'LineWidth', 1.5); % Plot G + Gb for diabetes with SMC
plot(tspan, [Gb Gb], 'k--', 'LineWidth', 1.5); % Basal glucose level
ylim([0, max([z_normal(:,1) + Gb; z_diabetes(:,1) + Gb; z_diabetes_smc(:,1) + Gb]) + 1]); % Set ylim range
text(5, Gb - 0.1 * (max([z_normal(:,1) + Gb; z_diabetes(:,1) + Gb; z_diabetes_smc(:,1) + Gb]) + 1), ...
    'Basal Glucose Level', 'FontSize', 14, 'Color', 'k'); % Annotation
xlabel('Time (minutes)', 'FontSize', 14);
ylabel('Glucose', 'FontSize', 14);
% legend('Normal', 'Diabetes', 'Diabetes with SMC', 'FontSize', 12);
title('Glucose Response', 'FontSize', 16);
set(gca, 'FontSize', 12); % Adjust tick font size
grid on;

% Right subplot: Control effort
subplot(3, 1, 1);
plot(t_diabetes_smc(1:end-1), diff(z_diabetes_smc(:,4))./diff(t_diabetes_smc), 'b', 'LineWidth', 1.5); % Plot control effort
xlabel('Time (minutes)', 'FontSize', 14);
ylabel('u(t) (mU/min)', 'FontSize', 14);
title('Control Effort', 'FontSize', 16);
set(gca, 'FontSize', 12); % Adjust tick font size
grid on;

% Save the figure
exportgraphics(gcf, 'images/insulin_glucose_response_smc.png', 'Resolution', 300);

% -------------------------------------------------- %
%       Sliding Mode Control without Chattering      %
% -------------------------------------------------- %

k = 0.5; eta = 3;
a = -6; b = 1; m = n; phi = 45;

u_smc = @(t, z) n * Vl * Ib + (Vl/b)*(a*l1_diabetes*z(1) + a*z(1)*z(2) + a*Gb*z(2) + b*m*z(3) + b*m*Ib - eta*boundary_layer(phi,a*z(1) + b*z(3)) - k*(a*z(1) + b*z(3)));

% Update the diabetic system with SMC
odeSystemDiabetesSMC = @(t, z) [
    -l1_diabetes*z(1) - z(1)*z(2) - Gb*z(2);         % z1_dot (Glucose dynamics)
    -l2*z(2) + l3*z(3);                              % z2_dot
    -n*z(3) - n*Ib + (u_smc(t, z)/Vl);                % z3_dot (Insulin dynamics with SMC)
    u_smc(t,z)];

[t_diabetes_smc, z_diabetes_smc] = ode45(odeSystemDiabetesSMC, tspan, z0);

figure(7);
set(gcf, 'Position', [100, 100, 700, 800]);
% Left subplot: Insulin response
subplot(3, 1, 2);
plot(t_normal, z_normal(:,3) + Ib, 'k--', 'LineWidth', 1.5);
hold on;
plot(t_diabetes_smc, z_diabetes_smc(:,3) + Ib, 'g', 'LineWidth', 1.5);
text(5, Ib - 7, ...
    'Basal Insulin Level', 'FontSize', 14, 'Color', 'k'); % Annotation
xlabel('Time (minutes)', 'FontSize', 14);
ylabel('Insulin', 'FontSize', 14);
ylim([0, 80]);
% legend('Normal', 'Diabetes with SMC', 'FontSize', 12);
title('Insulin Response', 'FontSize', 16);
set(gca, 'FontSize', 12); % Adjust tick font size
grid on;

% Middle subplot: Glucose response
subplot(3, 1, 3);
plot(t_normal, z_normal(:,1) + Gb, 'b', 'LineWidth', 1.5); % Plot G + Gb for normal case
hold on;
plot(t_diabetes, z_diabetes(:,1) + Gb, 'r', 'LineWidth', 1.5); % Plot G + Gb for diabetes
plot(t_diabetes_smc, z_diabetes_smc(:,1) + Gb, 'g', 'LineWidth', 1.5); % Plot G + Gb for diabetes with SMC
plot(tspan, [Gb Gb], 'k--', 'LineWidth', 1.5); % Basal glucose level
text(5, Gb - 0.1 * (max([z_normal(:,1) + Gb; z_diabetes(:,1) + Gb; z_diabetes_smc(:,1) + Gb]) + 1), ...
    'Basal Glucose Level', 'FontSize', 14, 'Color', 'k'); % Annotation
ylim([0, max([z_normal(:,1) + Gb; z_diabetes(:,1) + Gb; z_diabetes_smc(:,1) + Gb]) + 1]); % Set ylim range
xlabel('Time (minutes)', 'FontSize', 14);
ylabel('Glucose', 'FontSize', 14);

% legend('Normal', 'Diabetes', 'Diabetes with SMC', 'FontSize', 12);
title('Glucose Response', 'FontSize', 16);
set(gca, 'FontSize', 12); % Adjust tick font size
grid on;

% Right subplot: Control effort
subplot(3, 1, 1);
plot(t_diabetes_smc(1:end-1), diff(z_diabetes_smc(:,4))./diff(t_diabetes_smc), 'LineWidth', 1.5); % Plot control effort
hold on;
xlabel('Time (minutes)', 'FontSize', 14);
ylabel('u(t) (mU/min)', 'FontSize', 14);
title('Control Effort', 'FontSize', 16);
set(gca, 'FontSize', 12); % Adjust tick font size
grid on;

% Save the figure
exportgraphics(gcf, 'images/insulin_glucose_response_smc_boundary_layer.png', 'Resolution', 300);


% ------------------------------------------------------- %
%          Print the integer of the control effort        %
% ------------------------------------------------------- %

% Compute the integer of the control effort
control_effort = trapz(t_diabetes_smc(1:end-1), diff(z_diabetes_smc(:,4))./diff(t_diabetes_smc));

fprintf('The integer of the control effort is: %.2f\n', control_effort);

% Save parameters, arrays, and integral to a .mat file
save(sprintf('params/params_and_arrays_k%d_eta%d_a%d_b%d_phi%d.mat', k, eta, a, b, phi), 'k', 'eta', 'a', 'b', 'phi', 't_diabetes_smc', 'z_diabetes_smc', 'control_effort');

% -------------------------------- %
%          Phase Portrait          %
% -------------------------------- %

% % Quiver plot for phase portrait
% figure(2); % Use figure 2
% set(gcf, 'Position', [100, 100, 1200, 300]);

% % Define grid for quiver plot
% [G, I] = meshgrid(4.5:1:20, 15:1:30);

% % Initialize derivatives
% dG_normal = zeros(size(G));
% dI_normal = zeros(size(I));
% dG_diabetes = zeros(size(G));
% dI_diabetes = zeros(size(I));
% dG_diabetes_smc = zeros(size(G));
% dI_diabetes_smc = zeros(size(I));

% % Compute derivatives at each grid point
% for i = 1:numel(G)
%     z_normal_point = [G(i) - Gb; 0.01; I(i) - Ib;0];
%     dz_normal = odeSystemNormal(0, z_normal_point);
%     dG_normal(i) = dz_normal(1);
%     dI_normal(i) = dz_normal(3);

%     z_diabetes_point = [G(i) - Gb; 0.01; I(i) - Ib;0];
%     dz_diabetes = odeSystemDiabetes(0, z_diabetes_point);
%     dG_diabetes(i) = dz_diabetes(1);
%     dI_diabetes(i) = dz_diabetes(3);

%     z_diabetes_smc_point = [G(i) - Gb; 0.01; I(i) - Ib; 0];
%     dz_diabetes_smc = odeSystemDiabetesSMC(0, z_diabetes_smc_point);
%     dG_diabetes_smc(i) = dz_diabetes_smc(1);
%     dI_diabetes_smc(i) = dz_diabetes_smc(3);
% end

% % Plot quiver for normal case
% subplot(1, 3, 1);
% quiver(G, I, dG_normal, dI_normal, 'b');
% hold on;
% xlabel('Glucose', 'FontSize', 14);
% ylabel('Insulin', 'FontSize', 14);
% title('Normal', 'FontSize', 16);
% set(gca, 'FontSize', 12); % Adjust tick font size
% grid on;

% % Plot quiver for diabetic case
% subplot(1, 3, 2);
% quiver(G, I, dG_diabetes, dI_diabetes, 'r');
% hold on;
% xlabel('Glucose', 'FontSize', 14);
% ylabel('Insulin', 'FontSize', 14);
% title('Diabetes', 'FontSize', 16);
% set(gca, 'FontSize', 12); % Adjust tick font size
% grid on;

% % Plot quiver for diabetic case
% subplot(1, 3, 3);
% quiver(G, I, dG_diabetes_smc, dI_diabetes_smc, 'r');
% hold on;
% xlabel('Glucose', 'FontSize', 14);
% ylabel('Insulin', 'FontSize', 14);
% title('Diabetes with SMC', 'FontSize', 16);
% set(gca, 'FontSize', 12); % Adjust tick font size
% grid on;

% hold on;
% x = linspace(4.5, 20, 100);
% y = (1 - a * x) / b;
% plot(x, y, 'k--', 'LineWidth', 1.5);
% xlim([4.5, 20]);
% ylim([15, 30]);

% % Save the second figure
% exportgraphics(gcf, 'images/phase_portrait_smc.png', 'Resolution', 300);

