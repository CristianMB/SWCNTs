clc;
clear;
addpath('C:\Users\lucp12322\Desktop\shabnam\');
%addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\');
import UsefulFunctions.*;

%rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements\Raman\';
rootpath = 'C:\Users\lucp12322\Desktop\shabnam\';

%All paths as defaultC:\Users\lucp12322\Desktop\AFM CNT\shabnam\ShabnamM
path = [rootpath, '\ShabnamM\'];

%Select the paths of interest
paths = {
    path
    };




ReadRamanFromPaths(paths);


%%%--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%--------SAMPLE COMPARISION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


GDBand514 = {
                DATA_ShabnamM.S25424rb,
                DATA_ShabnamM.S25424s,
               DATA_ShabnamM.S25424q,

    };

GDBand514 = {
                   DATA_ShabnamM.S25424rb,
                DATA_ShabnamM.S25424s,
                DATA_ShabnamM.S25424q,
                DATA_ShabnamM.S25424m,
                DATA_ShabnamM.S25424o,
                DATA_ShabnamM.S25424p,

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





