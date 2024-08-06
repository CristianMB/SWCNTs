
function plot_absorption()

    import UsefulFunctions.*;
    addpath('X:\Measurements Data\');
    % Add the path to the directory containing your functions

    % Step 1: Select the .csv file
    [file_name, file_path] = uigetfile('*.csv', 'Select the CSV file');
    if isequal(file_name, 0)
        disp('No file selected. Exiting.');
        return;
    end

    % Full path to the selected file
    full_file_path = fullfile(file_path, file_name);
    paths = {full_file_path};
    % Step 2: Use ReadAbsorptionFromPaths to read the data
    % Assuming ReadAbsorptionFromPaths returns the required data structure

    % Step 2: Use ReadAbsorptionFromPaths to read the data
    ReadAbsorptionFromPaths(paths);
    
    % Step 3: Extract the date from the file path
    % Assuming the date is embedded in the file name in the format YYMMDD
    [~, name, ~] = fileparts(file_name);
    date_str = regexp(name, '\d{6}', 'match', 'once');
    if isempty(date_str)
        error('No date found in the file name.');
    end
    
    % Construct the variable name based on the date
    data_var_name = ['DATA_', date_str];
    
    % Check if the variable exists in the workspace
    if evalin('base', ['exist(''', data_var_name, ''', ''var'')'])
        % Step 4: Retrieve the data from the workspace
        data_struct = evalin('base', data_var_name);
        
        % Step 5: Plot the absorption spectrum using plotAbsorption
        % Assuming plotAbsorption takes the data structure and an optional offset
        figure;
        plotAbsorption(data_struct, 0);  % No offset provided, modify if needed

        % Customize the plot
        xlabel('Wavelength (nm)');
        ylabel('Absorption');
        title(['Absorption Spectrum: ', file_name]);
        legend('show');
    else
        error(['Variable ', data_var_name, ' does not exist in the workspace.']);
    end
end


