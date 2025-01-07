clc;
clear;
addpath('X:\SWCNTs');
% addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\');
import UsefulFunctions.*;

% rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements Data\Absorption\';
rootpath = 'X:\Measurements Data\ExternalData\';

%All paths as default
path_Examples = [rootpath,'Raman_2020412_TestData_LorentzianFitting_EmptyWater\'];


%Select the paths of interest

paths = {
    path_Examples
    };


ReadRamanFromPaths(paths,2);




%%%--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%--------SAMPLE COMPARISION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Explore the data
Examples = {
    DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.SR0H514R_241212
    DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.SR1H514R_241212
%     DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.FFH514R_241212
        
    DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.EAH514R_240517
    DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.WAH514R_240517
%    DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.FF514R_240517
    
    DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.EAH514R_241006
    DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.WAH514R_241006
%     DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.FFH514R_241006
        
    DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.S240111J_240111
    DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.S240111S_240111
%     DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.FLATHD_240111
}


%Order the data
D241212 = {    
        DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.SR0H514R_241212
        DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.SR1H514R_241212
        };
    
D240517 = {
        DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.EAH514R_240517
        DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.WAH514R_240517    
            };
        
D241016 = {
        DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.EAH514R_241006
        DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.WAH514R_241006
            };
        
D240111 = {
        DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.S240111S_240111
        DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.S240111J_240111

            };

%Correct the data
 
D241212 = FlatFieldCorrection(D241212, DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.FFH514R_241212);
D240517 = FlatFieldCorrection(D240517, DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.FF514R_240517);
D241016 = FlatFieldCorrection(D241016, DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.FFH514R_241006);
D240111 = FlatFieldCorrection(D240111, DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.FLATHD_240111); 

%Filter the data (Mathing ranges)

D241212 = clipRangeEdges(D241212, 135, 220);
D240517 = clipRangeEdges(D240517, 135, 220);
D241016 = clipRangeEdges(D241016, 135, 220);
D240111 = clipRangeEdges(D240111, 135, 220);


D241212 = NormalizeSample(D241212, 135, 220)
D240517 = NormalizeSample(D240517, 135, 220)
D241016 = NormalizeSample(D241016, 135, 220)
D240111 = NormalizeSample(D240111, 135, 220)

AllData = {
            D241212{1}
            D241212{2}
            D240517{1}
            D240517{2}
            D241016{1}
            D241016{2}
            D240111{1}
            D240111{2}
            }
AllData = NormalizeSample(AllData, 135, 220)


%Plotting the data
% plotRaman(D241212, 0)
% plotRaman(D240517, 0)
% plotRaman(D241016, 0)
% plotRaman(D240111, 0)

%Export the data

currentDir = pwd; 

%Offset
D241212{1}.Y=D241212{1}.Y+1

AllData{2}.Y = AllData{2}.Y + 1
AllData{3}.Y = AllData{3}.Y + 2
AllData{4}.Y = AllData{4}.Y + 3
AllData{5}.Y = AllData{5}.Y + 4
AllData{6}.Y = AllData{6}.Y + 5
AllData{7}.Y = AllData{7}.Y + 6
AllData{8}.Y = AllData{8}.Y + 7

plotRaman(AllData,0)

ExportMatFile(D241212, [currentDir,'\D241212.mat'])
ExportMatFile(D240517, [currentDir,'\D240517.mat'])
ExportMatFile(D241016, [currentDir,'\D241016.mat'])
ExportMatFile(D240111, [currentDir,'\D240111.mat'])
ExportMatFile(AllData, [currentDir,'\AllData.mat'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ExportMatFile(DS, filename)

    % Initialize variables to store the output data
    numSpectra = length(DS); % Number of spectra (elements in DS)
    xdata = DS{1}.X;         % Use the .X attribute of the first structure as the x-axis
    numPoints = length(xdata); % Number of points in xData
    
    % Preallocate the matrix for efficiency
    ydata = zeros(numPoints, numSpectra); % Rows: xdata, Columns: spectra    
    
    % Loop through each spectrum to collect .Y data
    for i = 1:numSpectra     
        % Assign .Y to the corresponding column
        ydata(:, i) = DS{i}.Y;
    end
    
    % Save the data matrix as a .mat file
    save(filename, 'xdata', 'ydata');
    
    % Notify the user
    fprintf('Data exported successfully to %s\n', filename);
end


function DS = clipRangeEdges(DS, min_value, max_value)
    % Iterate over each element in the array of data structures
    for i = 1:length(DS)
        % Find the indices where X is within the specified range
        idx = DS{i}.X >= min_value & DS{i}.X <= max_value;
        
        % Filter X, Y, and any other fields that have the same number of elements as X
        DS{i}.X = DS{i}.X(idx);
        DS{i}.Y = DS{i}.Y(idx);
    end
end

