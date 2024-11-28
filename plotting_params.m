Ib = 15;
Gb = 4.5;
tspan = [0 200];

% Define the folder where the .mat files are stored
folder = 'params_k';

% Get a list of all .mat files in the folder
filePattern = fullfile(folder, '*.mat');
matFiles = dir(filePattern);

figure(3);
set(gcf, 'Position', [200, -10, 600, 500]); % Increase the width to accommodate the legend

numFiles=5;
pink = [255, 128, 0]/255;
red = [0, 128, 255]/255;

colors_p = [linspace(red(1),pink(1),numFiles)', linspace(red(2),pink(2),numFiles)', linspace(red(3),pink(3),numFiles)'];

% Initialize a cell array to store legend entries
legendEntries = cell(numel(numFiles), 1);

% Loop through each .mat file
for k = 1:numel(matFiles)
    % Get the full file name
    baseFileName = matFiles(k).name;
    fullFileName = fullfile(folder, baseFileName);
    
    % Load the .mat file
    data = load(fullFileName);
    
    % Left subplot: Insulin response
    subplot(3, 1, 2);
    plot(data.t_diabetes_smc, data.z_diabetes_smc(:,3) + Ib, 'Color', [colors_p(k, :) 0.7], 'LineWidth', 1.1);
    hold on;
    xlabel('Time (minutes)', 'FontSize', 14);
    ylabel('Insulin', 'FontSize', 14);
    ylim([-5, 80]);
    title('Insulin Response', 'FontSize', 16);
    set(gca, 'FontSize', 12); % Adjust tick font size
    grid on;
    
    % Middle subplot: Glucose response
    subplot(3, 1, 3);
    plot(data.t_diabetes_smc, data.z_diabetes_smc(:,1) + Gb, 'Color', [colors_p(k, :) 0.7], 'LineWidth', 1.1); % Plot G + Gb for diabetes with SMC
    hold on;
    plot(tspan, [Gb Gb], 'k--'); % Basal glucose level
    ylim([0, 20]); % Set ylim range
    xlabel('Time (minutes)', 'FontSize', 14);
    ylabel('Glucose', 'FontSize', 14);
    
    % legend('Normal', 'Diabetes', 'Diabetes with SMC', 'FontSize', 12);
    title('Glucose Response', 'FontSize', 16);
    set(gca, 'FontSize', 12); % Adjust tick font size
    grid on;
    
    % Right subplot: Control effort
    subplot(3, 1, 1);
    plot(data.t_diabetes_smc(1:end-1), diff(data.z_diabetes_smc(:,4))./diff(data.t_diabetes_smc), 'Color', [colors_p(k, :) 0.7], 'LineWidth', 1.1); % Plot control effort
    hold on;
    xlabel('Time (minutes)', 'FontSize', 14);
    ylabel('u(t) (mU/min)', 'FontSize', 14);
    title('Control Effort', 'FontSize', 16);
    set(gca, 'FontSize', 12); % Adjust tick font size
    grid on;
    
    % Check if the required variables exist and print them
    if isfield(data, 'k') && isfield(data, 'eta') && isfield(data, 'a') && isfield(data, 'b') && isfield(data, 'phi')
        fprintf('File: %s\n', baseFileName);
        fprintf('k: %f\n', data.k);
        fprintf('eta: %f\n', data.eta);
        fprintf('a: %f\n', data.a);
        fprintf('b: %f\n', data.b);
        fprintf('phi: %f\n', data.phi);
        fprintf('control_effort: %f\n', data.control_effort);
        fprintf('\n');
        
        % Store the legend entry
        legendEntries{k} = sprintf('$k = %2.1f$', data.k);
    else
        fprintf('File: %s does not contain all required variables.\n', baseFileName);
    end
end

% Create the legend outside the loop
subplot(3, 1, 1);
legendHandle = legend(legendEntries, 'FontSize', 12, 'Location', 'northeastoutside', 'Interpreter', 'latex');
subplot(3, 1, 2);
legendHandle = legend(legendEntries, 'FontSize', 12, 'Location', 'northeastoutside', 'Interpreter', 'latex');
subplot(3, 1, 3);
legendHandle = legend(legendEntries, 'FontSize', 12, 'Location', 'northeastoutside', 'Interpreter', 'latex');

exportgraphics(gcf, 'images/varying_params.png', 'Resolution', 300);
