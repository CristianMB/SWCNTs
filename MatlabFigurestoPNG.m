function exportFiguresAsPNG()
    % Prompt user to select the folder containing .fig files
    inputFolder = uigetdir(pwd, 'Select the folder containing the .fig files');

    % Check if the user canceled the dialog
    if isequal(inputFolder, 0)
        disp('User canceled the operation.');
        return;
    end

    % Create an output folder for PNG figures
    outputFolder = fullfile(inputFolder, 'PNG_Figures');
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);  % Create the subfolder if it doesn't exist
    end

    % Get a list of all .fig files in the selected folder
    figFiles = dir(fullfile(inputFolder, '*.fig'));

    % Loop through each .fig file
    for k = 1:length(figFiles)
        % Construct the full file name and path
        figFileName = fullfile(inputFolder, figFiles(k).name);
        
        % Open the .fig file
        fig = openfig(figFileName, 'invisible');  % Open the figure invisibly

        % Get the name of the file without extension
        [~, name, ~] = fileparts(figFiles(k).name);
        
        % Save the figure as a PNG file with 600 dpi
        exportFileName = fullfile(outputFolder, [name, '.png']);  % Change to '.jpg' for JPG

        % Use 'print' to set resolution (600 dpi in this case)
        print(fig, exportFileName, '-dpng', '-r600');  % Use '-r600' for 600 dpi resolution
        
        % Close the figure after saving
        close(fig);
    end

    disp(['All figures have been exported to ', outputFolder]);
end
