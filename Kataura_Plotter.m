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
title('Excitation/Emission Wavelengths vs. RBM');
xlabel('RBM frequency (cm^{-1})');
ylabel('Wavelength (nm)');
for i = 1:length(KATAURA.RBM)
    if strcmp(KATAURA.Type(i), 'M')
        color = colors{1}; % Metallic color
    else
        color = colors{3}; % Semiconducting color
    end
    
    % Loop through each energy level and plot with different marker
    for energy = 1:4
        marker = markers{energy};
        scatter(KATAURA.RBM(i), KATAURA.(['WL' num2str(energy)])(i), 50, color, marker, 'filled');
        text(KATAURA.RBM(i), KATAURA.(['WL' num2str(energy)])(i), KATAURA.Chirality{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    end
end