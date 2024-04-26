clc;
clear;
addpath('C:\Users\cborja\OneDrive - Universiteit Antwerpen\SWCNTs\');
%addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\');
import UsefulFunctions.*;

%rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements\Raman\';
rootpath = 'C:\Users\cborja\OneDrive - Universiteit Antwerpen\Measurements Data\Raman\';

%All paths as default
path = [rootpath,'\ShabnamMeasurements\'];

%Select the paths of interest
paths = {
    path
    };




ReadRamanFromPaths(paths);


%%%--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%--------SAMPLE COMPARISION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


GDBand514 = {
                DATA_ShabnamMeasurements.S25424RB,
                DATA_ShabnamMeasurements.S25424S,
                DATA_ShabnamMeasurements.S25424Q,

    };

GDBand514 = {
                    DATA_ShabnamMeasurements.S25424RB,
                DATA_ShabnamMeasurements.S25424S,
                DATA_ShabnamMeasurements.S25424Q,
                DATA_ShabnamMeasurements.S25424M,
                DATA_ShabnamMeasurements.S25424O,
                DATA_ShabnamMeasurements.S25424P,

    };
    
LG = 1560;
HG = 1570;






plotRaman(GDBand514, 0)





%%%%%TODO
%Use an excel file to handle the labeling

function labeledData = assignLabelsFromExcel(data, excelFile)
    % Read the Excel file
    excelData = readtable(excelFile);
    
    % Extract Raman and N columns
    ramanColumn = excelData.Raman;
    nColumn = excelData.N;
    
    % Initialize labeled data structure
    labeledData = struct();
    
    % Loop through each variable in the input data
    variableNames = fieldnames(data);
    for i = 1:numel(variableNames)
        variableName = variableNames{i};
        
        % Find the corresponding N value in the Excel data
        index = strcmp(variableName, ramanColumn);
        if any(index)
            nValue = nColumn(index);
            
            % Assign N value to the labeled data structure
            labeledData.(variableName) = data.(variableName);
            labeledData.(variableName).N = nValue;
        else
            % If the variable name is not found in the Excel data, display a warning
            warning(['Variable "', variableName, '" not found in Excel data.']);
        end
    end
end





