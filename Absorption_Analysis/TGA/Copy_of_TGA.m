% MATLAB Script: plot_weight_derivative_multiple.m

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
    
    % Read numeric data (skip header lines)
    data = readmatrix(tempfile, 'NumHeaderLines', 2);
    
    % Extract columns: Temperature (6th) and Weight (5th)
    Temperature = data(:,6);
    Weight = data(:,5);
    
    % Compute numerical derivative d(Weight)/d(Temperature)
    dWdT = gradient(Weight, Temperature);
    
    % Plot derivative vs. Temperature
    plot(Temperature, dWdT, colors{i}, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('File %d', i));
end

xlabel('Temperature [°C]');
ylabel('d(Weight)/d(Temperature) [%/°C]');
title('Derivative of Weight vs. Temperature (Multiple Files)');
legend('show');
grid on;
hold off;
