
function plot_absorption()

    import UsefulFunctions.*;
    % Add the path to the directory containing your functions
    defaultpath=('X:\Measurements Data\Absorption');
    addpath(' X:\SWCNTs\SpecialMatlabFunctions\DrosteEffect-BrewerMap-3.2.5.0');
    addpath(' X:\SWCNTs');
    % Step 1: Select the .csv file
    [file_name, file_path] = uigetfile('*.csv', 'Select the CSV file', defaultpath);
    if isequal(file_name, 0)
        disp('No file selected. Exiting.');
        return;
    end
 
    
    % Full path to the selected file
    full_file_path = fullfile(file_path, file_name);
    paths = {full_file_path};
   
    A = ReadAbsorptionFromPaths(paths);
    Datasets = fieldnames(A);
    Spectra = {};

     for i = 1:length(Datasets)
            dataset_name = Datasets{i};
            spectra_fields = fieldnames(A.(dataset_name));
            % Store the spectra data references in the Spectra cell array
            for j = 1:length(spectra_fields)
                 if ~contains(spectra_fields{j}, 'Baseline')
                    Spectra{end+1} = A.(dataset_name).(spectra_fields{j});
                 end
            end
     end
     
    Spectra = NormalizeSample(Spectra, 1000, 1050)
%      Spectra = NormalizeSample(Spectra, 445, 448)

    plotAbsorptionOrdered(Spectra, 0.0)
%     plotAbsorption(Spectra, 0.0)

end


