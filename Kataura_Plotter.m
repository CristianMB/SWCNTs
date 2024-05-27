%% Initialize KATAURA structure
clc;
clear;
import UsefulFunctions.*;


KATAURA.RBM = [];
KATAURA.WL1 = [];
KATAURA.WL2 = [];
KATAURA.WL3 = [];
KATAURA.WL4 = [];
KATAURA.D = [];
KATAURA.InverseD = [];
KATAURA.Theta = [];
KATAURA.M = [];
KATAURA.N = [];
KATAURA.Chirality = [];
KATAURA.Type = [];

% Kataura plot calculation
for m = 5:18
    for n = 0:m
        [rbm, wl1, wl2, wl3, wl4, diam, theta, type] = CalculateKataura([n, m]);
        KATAURA.RBM = [KATAURA.RBM, rbm];
        KATAURA.WL1 = [KATAURA.WL1, wl1];
        KATAURA.WL2 = [KATAURA.WL2, wl2];
        KATAURA.WL3 = [KATAURA.WL3, wl3];
        KATAURA.WL4 = [KATAURA.WL4, wl4];
        KATAURA.D = [KATAURA.D, diam];
        KATAURA.InverseD = [KATAURA.InverseD, 1/diam];
        KATAURA.Theta = [KATAURA.Theta, theta];
        KATAURA.M = [KATAURA.M, m];
        KATAURA.N = [KATAURA.N, n];
        KATAURA.Type = [KATAURA.Type, {sprintf('%s', type)}];
        KATAURA.Chirality = [KATAURA.Chirality, {sprintf('(%d,%d)', m, n)}]; % Store formatted string in cell array
    end
end

%% Plotting
colors = {[0 0.4470 0.7410], [0.4660 0.6740 0.1880], [0.6350 0.0780 0.1840], [0.9290 0.6940 0.1250]}; % Elegant colors
markers = {'v', 'diamond', 'o', 's'}; % Markers for each energy level

figure;
hold on;
title('KatauraPlot');
xlabel('Diameter (nm)');
ylabel('Wavelength (nm)');

for i = 1:length(KATAURA.RBM)
    if strcmp(KATAURA.Type(i), 'M')
        color = colors{1}; % Metallic color
    else
        color = colors{3}; % Semiconducting color
    end
    
    % Loop through each energy level and plot with different marker
    for energy = 1:3
        marker = markers{energy};
        scatter(KATAURA.D(i), KATAURA.(['WL' num2str(energy)])(i), 50, color, marker, 'filled');
%         text(KATAURA.RBM(i), KATAURA.(['WL' num2str(energy)])(i), KATAURA.Chirality{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    end
end

% Add custom legend
h1 = scatter(nan, nan, 50, colors{1}, 'v', 'filled');
h2 = scatter(nan, nan, 50, colors{3}, 'v', 'filled');
h3 = scatter(nan, nan, 50, colors{1}, 'diamond', 'filled');
h4 = scatter(nan, nan, 50, colors{3}, 'diamond', 'filled');
h5 = scatter(nan, nan, 50, colors{1}, 'o', 'filled');
h6 = scatter(nan, nan, 50, colors{3}, 'o', 'filled');
% h7 = scatter(nan, nan, 50, colors{1}, 's', 'filled');
% h8 = scatter(nan, nan, 50, colors{3}, 's', 'filled');

legend([h1, h2, h3, h4, h5, h6], {'M11', 'S11','M22', 'S22','M33', 'S33'}, 'Location', 'northeast');
%legend([h1, h2, h3, h4, h5, h6, h7, h8], {'M11', 'S11','M22', 'S22','M33', 'S33','M44', 'S44'}, 'Location', 'northeast');



% %% Plotting
% colors = {[0 0.4470 0.7410], [0.4660 0.6740 0.1880], [0.6350 0.0780 0.1840], [0.9290 0.6940 0.1250]}; % Elegant colors
% markers = {'v', 'diamond', 'o', 's'}; % Markers for each energy level
% 
% figure;
% hold on;
% title('KatauraPlot');
% xlabel('Diameter (nm)');
% ylabel('Wavelength (nm)');
% for i = 1:length(KATAURA.RBM)
%     if strcmp(KATAURA.Type(i), 'M')
%         color = colors{1}; % Metallic color
%     else
%         color = colors{3}; % Semiconducting color
%     end
%     
%     % Loop through each energy level and plot with different marker
%     for energy = 1:2
%         marker = markers{energy};
%         scatter(KATAURA.D(i), KATAURA.(['WL' num2str(energy)])(i), 50, color, marker, 'filled');
%         text(KATAURA.D(i), KATAURA.(['WL' num2str(energy)])(i), KATAURA.Chirality{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
%     end
% end