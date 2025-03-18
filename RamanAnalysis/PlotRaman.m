
function plot_raman()

    import UsefulFunctions.*;
    % Add the path to the directory containing your functions
    addpath('X:\SWCNTs');
    addpath(' X:\SWCNTs\SpecialMatlabFunctions\DrosteEffect-BrewerMap-3.2.5.0');
    defaultpath=('X:\Measurements Data\Raman');
    % Step 1: Select the .csv file
    folder_path  = uigetdir(defaultpath);
    if isequal(folder_path, 0)
        disp('No file selected. Exiting.');
        return;
    end
    paths = {[folder_path, '\']};

    % Full path to the selected file   
    A = UsefulFunctions.ReadRamanFromPaths(paths, 2);
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
     
%     Spectra = NormalizeSample(Spectra, 120, 210)
    plotRaman(Spectra, 0)
    
end


