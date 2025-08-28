% MATLAB Script: plot_weight_vs_temperature_multiple.m

files = {'20250826-01 P2 unfilled.txt', '20250826-02 TTF filled.txt', '20250825-01 TCNQ filled.txt'}; % <-- Change to your file names
colors = {'b', 'r', 'g'}; % Colors for each file

figure;
hold on;

for i = 1:length(files)
    % Read file as text and replace commas with dots
    rawText = fileread(files{i});
    rawText = strrep(rawText, ',', '.');
    
    % Save to temporary file
    tempfile = ['temp_' num2str(i) '.txt'];
    fid = fopen(tempfile, 'w');
    fwrite(fid, rawText);
    fclose(fid);
    
    % Read numeric data (skip first 2 header lines)
    data = readmatrix(tempfile, 'NumHeaderLines', 2);
    
    % Extract columns: Weight (5th) vs. Temperature (6th)
    Temperature = data(:,6);
    Weight = data(:,5);
    
    % Plot
    plot(Temperature, Weight, colors{i}, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('File %d', i));
end

xlabel('Temperature [°C]');
ylabel('Weight [%]');
title('Weight vs. Temperature (Multiple Files)');legend('show');
grid on;
hold off;
