function plot_individualRaman()

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

    plotRaman(SpectraArray, 10.0)
end 