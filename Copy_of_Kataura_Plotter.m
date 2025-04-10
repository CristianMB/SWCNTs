%% Initialize KATAURA structure
clc;
clear;
import UsefulFunctions.*;



%% Choose region
% XAxis = [50, 600];  %RBM Range
XAxis = [0, 4];  %Diameter Range
YAxis = [0, 800];  %WL Range
% YAxis = [450, 520];  %WL Range

% Laserlines = [457.9, 476.50, 496.50, 514.5, 561.0];
% Tolerances = [5, 5, 5, 5, 5];
% Laserlines = [514.5 561.1 496.5 488 476.5 457.9];
% Tolerances = [5 5 5 5 5 5];
Laserlines = [];
Tolerances = [];

hc = 1240.84193;    %h*c value to convert energy to nm

KATAURA.Family = [];
KATAURA.RBM = [];
KATAURA.WL1 = [];
KATAURA.WL2 = [];
KATAURA.WL3 = [];
KATAURA.WL4 = [];
KATAURA.E1 = [];
KATAURA.E2 = [];
KATAURA.E3 = [];
KATAURA.E4 = [];
KATAURA.D = [];
KATAURA.InverseD = [];
KATAURA.Theta = [];
KATAURA.M = [];
KATAURA.N = [];
KATAURA.Chirality = [];
KATAURA.Type = [];

% Kataura plot calculation

% for m = 5:2
%     for n = 0:m
for m = 5:16
    for n = 0:m
        [rbm, wl1, w22, w33, w44, diam, theta, type] = CalculateKataura([n, m]);
        
        rbm_energy = rbm * 1.239841984e-4; % Convert cm^-1 to eV
        % Apply Stokes correction to each wavelength
        wl1_corrected = 1239.841984 / (1239.841984 / wl1 - rbm_energy / 2);
        w22_corrected = 1239.841984 / (1239.841984 / w22 - rbm_energy / 2);
        w33_corrected = 1239.841984 / (1239.841984 / w33 - rbm_energy / 2);
        w44_corrected = 1239.841984 / (1239.841984 / w44 - rbm_energy / 2);
        
        KATAURA.RBM = [KATAURA.RBM, rbm];
        KATAURA.WL1 = [KATAURA.WL1, wl1];
        KATAURA.WL2 = [KATAURA.WL2, w22];
        KATAURA.WL3 = [KATAURA.WL3, w33];
        KATAURA.WL4 = [KATAURA.WL4, w44];
        KATAURA.E1 = [KATAURA.E1, hc./wl1];
        KATAURA.E2 = [KATAURA.E2, hc./w22];
        KATAURA.E3 = [KATAURA.E3, hc./w33];
        KATAURA.E4 = [KATAURA.E4, hc./w44];
        KATAURA.D = [KATAURA.D, diam];
        KATAURA.InverseD = [KATAURA.InverseD, 1/diam];
        KATAURA.Theta = [KATAURA.Theta, theta];
        KATAURA.M = [KATAURA.M, m];
        KATAURA.N = [KATAURA.N, n];
        KATAURA.Type = [KATAURA.Type, {sprintf('%s', type)}];
        KATAURA.Chirality = [KATAURA.Chirality, {sprintf('(%d,%d)', m, n)}]; % Store formatted string in cell array
    end
end



%% Choose region
% XAxis = [50, 600];  %RBM Range
XAxis = [0, 4];  %Diameter Range
YAxis = [0, 800];  %WL Range
% YAxis = [450, 520];  %WL Range

% Laserlines = [457.9, 476.50, 496.50, 514.5, 561.0];
% Tolerances = [5, 5, 5, 5, 5];
% Laserlines = [514.5 561.1 496.5 488 476.5 457.9];
% Tolerances = [5 5 5 5 5 5];
Laserlines = [];
Tolerances = [];
% Plotting
colors = {[0 0.4470 0.7410], [0.4660 0.6740 0.1880], [0.6350 0.0780 0.1840], [0.9290 0.6940 0.1250]};
markers = {'v', 'diamond', 'o', 's'}; % Markers for each energy level

figure;
hold on;
title('KatauraPlot');
xlabel('Diameter (wn)');
ylabel('Wavelength (nm)');

% Plotting the data points
for i = 1:length(KATAURA.D)
    % Determine the color based on the type
    if strcmp(KATAURA.Type{i}, 'M')
        color = colors{1}; % Metallic color
    else
        color = colors{3}; % Semiconducting color
    end
    
    % Loop through each energy level and plot with different markers
    for energy = 1:4
        marker = markers{energy};
        
        % Get the RBM/Diameter and Wavelength for the current energy level
        x_val = KATAURA.D(i);
        y_val = KATAURA.(['WL' num2str(energy)])(i);
        
        % Only plot and label if within the specified limits
        if x_val >= XAxis(1) && x_val <= XAxis(2) && y_val >= YAxis(1) && y_val <= YAxis(2)
            scatter(x_val, y_val, 50, color, marker, 'filled');
            text(x_val, y_val, KATAURA.Chirality{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
        end
    end
end

% Group by (2m+n) family and plot lines connecting points within each energy level
families = unique(2 * [KATAURA.M] + [KATAURA.N]);
for family = families
    for energy = 1:4
        % Find indices of CNTs in the current family and energy level
        family_indices = find((2 * [KATAURA.M] + [KATAURA.N] == family) & ...
                              ~isnan(KATAURA.(['WL' num2str(energy)])) & ...
                              KATAURA.D >= XAxis(1) & KATAURA.D <= XAxis(2) & ...
                              KATAURA.(['WL' num2str(energy)]) >= YAxis(1) & KATAURA.(['WL' num2str(energy)]) <= YAxis(2));
        
        if length(family_indices) > 1
            % Sort by X
            [~, sort_indices] = sort(KATAURA.D(family_indices));
            sorted_indices = family_indices(sort_indices);
            
            % Plot line connecting the family within the current energy level
            plot(KATAURA.D(sorted_indices), KATAURA.(['WL' num2str(energy)])(sorted_indices), 'k--');
        end
    end
end

% Adding laser lines if they exist
if ~isempty(Laserlines) && ~isempty(Tolerances)
    for i = 1:length(Laserlines)
        y1 = Laserlines(i) - Tolerances(i);
        y2 = Laserlines(i) + Tolerances(i);
        patch([min(XAxis) max(XAxis) max(XAxis) min(XAxis)], [y1 y1 y2 y2], ...
              [0.8 0.8 0.8], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end
end

% Adding custom legend
h1 = scatter(nan, nan, 50, colors{1}, 'v', 'filled');
h2 = scatter(nan, nan, 50, colors{3}, 'v', 'filled');
h3 = scatter(nan, nan, 50, colors{1}, 'diamond', 'filled');
h4 = scatter(nan, nan, 50, colors{3}, 'diamond', 'filled');
h5 = scatter(nan, nan, 50, colors{1}, 'o', 'filled');
h6 = scatter(nan, nan, 50, colors{3}, 'o', 'filled');
h7 = scatter(nan, nan, 50, colors{1}, 's', 'filled');
h8 = scatter(nan, nan, 50, colors{3}, 's', 'filled');

legend([h1, h2, h3, h4, h5, h6, h7, h8], ...
    {'M11', 'S11', 'M22', 'S22', 'M33', 'S33', 'M44', 'S44'}, 'Location', 'northeast');

% Set axis limits after plotting all data
xlim(XAxis);
ylim(YAxis);


%%%%%%PLOT RBM


%% Choose region
XAxis = [50, 600];  %RBM Range
% XAxis = [1, 2];  %Diameter Range
YAxis = [400, 600];  %WL Range
% YAxis = [450, 520];  %WL Range

% Laserlines = [457.9, 476.50, 496.50, 514.5, 561.0];
% Tolerances = [5, 5, 5, 5, 5];
% Laserlines = [514.5 561.1 496.5 488 476.5 457.9];
% Tolerances = [5 5 5 5 5 5];
Laserlines = [];
Tolerances = [];
figure;
hold on;
title('KatauraPlot');
xlabel('RBM frequency (cm-1)');
ylabel('Wavelength (nm)');

% Plotting the data points
for i = 1:length(KATAURA.RBM)
    % Determine the color based on the type
    if strcmp(KATAURA.Type{i}, 'M')
        color = colors{1}; % Metallic color
    else
        color = colors{3}; % Semiconducting color
    end
    
    % Loop through each energy level and plot with different markers
    for energy = 1:4
        marker = markers{energy};
        
        % Get the RBM/Diameter and Wavelength for the current energy level
        x_val = KATAURA.RBM(i);
        y_val = KATAURA.(['WL' num2str(energy)])(i);
        
        % Only plot and label if within the specified limits
        if x_val >= XAxis(1) && x_val <= XAxis(2) && y_val >= YAxis(1) && y_val <= YAxis(2)
            scatter(x_val, y_val, 50, color, marker, 'filled');
            text(x_val, y_val, KATAURA.Chirality{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
        end
    end
end

% Group by (2m+n) family and plot lines connecting points within each energy level
families = unique(2 * [KATAURA.M] + [KATAURA.N]);
for family = families
    for energy = 1:4
        % Find indices of CNTs in the current family and energy level
        family_indices = find((2 * [KATAURA.M] + [KATAURA.N] == family) & ...
                              ~isnan(KATAURA.(['WL' num2str(energy)])) & ...
                              KATAURA.RBM >= XAxis(1) & KATAURA.RBM <= XAxis(2) & ...
                              KATAURA.(['WL' num2str(energy)]) >= YAxis(1) & KATAURA.(['WL' num2str(energy)]) <= YAxis(2));
        
        if length(family_indices) > 1
            % Sort by X
            [~, sort_indices] = sort(KATAURA.RBM(family_indices));
            sorted_indices = family_indices(sort_indices);
            
            % Plot line connecting the family within the current energy level
            plot(KATAURA.RBM(sorted_indices), KATAURA.(['WL' num2str(energy)])(sorted_indices), 'k--');
        end
    end
end

% Adding laser lines if they exist
if ~isempty(Laserlines) && ~isempty(Tolerances)
    for i = 1:length(Laserlines)
        y1 = Laserlines(i) - Tolerances(i);
        y2 = Laserlines(i) + Tolerances(i);
        patch([min(XAxis) max(XAxis) max(XAxis) min(XAxis)], [y1 y1 y2 y2], ...
              [0.8 0.8 0.8], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end
end

% Adding custom legend
h1 = scatter(nan, nan, 50, colors{1}, 'v', 'filled');
h2 = scatter(nan, nan, 50, colors{3}, 'v', 'filled');
h3 = scatter(nan, nan, 50, colors{1}, 'diamond', 'filled');
h4 = scatter(nan, nan, 50, colors{3}, 'diamond', 'filled');
h5 = scatter(nan, nan, 50, colors{1}, 'o', 'filled');
h6 = scatter(nan, nan, 50, colors{3}, 'o', 'filled');
h7 = scatter(nan, nan, 50, colors{1}, 's', 'filled');
h8 = scatter(nan, nan, 50, colors{3}, 's', 'filled');

legend([h1, h2, h3, h4, h5, h6, h7, h8], ...
    {'M11', 'S11', 'M22', 'S22', 'M33', 'S33', 'M44', 'S44'}, 'Location', 'northeast');

% Set axis limits after plotting all data
xlim(XAxis);
ylim(YAxis);



