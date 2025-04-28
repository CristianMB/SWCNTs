
import UsefulFunctions.*;
    % Add the path to the directory containing your functions
    defaultpath = 'X:\Measurements Data\Raman';
    addpath('X:\SWCNTs\');
    addpath('X:\SWCNTs\SpecialMatlabFunctions\DrosteEffect-BrewerMap-3.2.5.0')

    % Step 1: Select the .csv file
    [file_name, folder_path] = uigetfile(fullfile(defaultpath, '*.m3d'), 'Select Raman Data File');

    % Check if the user canceled the file selection
    if isequal(file_name, 0)
        disp('No file selected. Exiting.');
        return;
    end

    % Full path to the selected file   
    full_file_path = fullfile(folder_path, file_name);

    % Step 2: Use the full file path to read the file
    S = ReadIndividualRamanFromFile(full_file_path);


    [~, fileName, ~] = fileparts(file_name);  % Remove the extension (.m3d)
    DataName = ['DATA_', upper(fileName)];    % Construct the dynamic variable name
    assignin('caller', DataName, S);


    Datasets = fieldnames(S);
    Spectra = {};

     for i = 1:length(Datasets)
            dataset_name = Datasets{i};
            spectra_fields = fieldnames(S.(dataset_name));
            % Store the spectra data references in the Spectra cell array
            for j = 1:length(spectra_fields)
                Spectra{end+1} = S.(dataset_name).(spectra_fields{j});
            end
     end
%      
% Spectra = NormalizeSample(Spectra, 0, 3000)
% plotRaman(, 0.0)
fields = fieldnames(S);

% Initialize an empty cell array to store the spectra fields
SpectraArray = cell(1, length(fields));

% Loop through each field and store it in the SpectraArray
for i = 1:length(fields)
    SpectraArray{i} = S.(fields{i});
end

plotRaman(SpectraArray, 0.0)

% Step 3: Extract X axis (Raman shift) and first 10 spectra
numSpectraToSave = min(10, length(SpectraArray));  % In case there are fewer than 10 spectra
raman_shift = SpectraArray{1}.X;  % Assume all spectra share the same X axis

% Preallocate matrix: 1 column for X, and 10 for Y values
export_data = zeros(length(raman_shift), numSpectraToSave + 1);
export_data(:,1) = raman_shift;

for i = 1:numSpectraToSave
    y_vals = SpectraArray{i}.Y;

    % Optional cleanup of very small numbers
    y_vals(abs(y_vals) < 1e-10) = 0;

    export_data(:, i+1) = y_vals;
end

% Define export filename
export_filename = fullfile(folder_path, [fileName '_ExportedSpectra.csv']);
fid = fopen(export_filename, 'w');

% Write header
fprintf(fid, 'RamanShift');
for i = 1:numSpectraToSave
    fprintf(fid, ',Spectrum%d', i);
end
fprintf(fid, '\n');

% Write data rows with 4 decimal places
for row = 1:length(raman_shift)
    fprintf(fid, '%.4f', export_data(row, 1));  % Raman shift
    for col = 2:(numSpectraToSave + 1)
        fprintf(fid, ',%.4f', export_data(row, col));  % Spectra
    end
    fprintf(fid, '\n');
end

fclose(fid);
fprintf('✅ Spectra exported to "%s" with 4 decimal precision.\n', export_filename);
