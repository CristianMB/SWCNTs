clc;
clear;
addpath('X:\SWCNTs');
addpath('X:\SWCNTs\SpecialMatlabFunctions\DrosteEffect-BrewerMap-3.2.5.0');

% addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\');
import UsefulFunctions.*;

% rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements Data\Absorption\';
rootpath = 'X:\Measurements Data\';

%All paths as default
path_Examples = [rootpath,'ExternalData\Raman_2020412_TestData_LorentzianFitting_EmptyWater\'];
path_FilledA = [rootpath,'Raman\20240517\'];
path_FilledB = [rootpath,'Raman\20241212\'];


%Select the paths of interest

paths = {
    path_Examples
    path_FilledA
    path_FilledB
    };


ReadRamanFromPaths(paths,2);
currentDir = pwd; 




%%%--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DATA_20240517.S2H514R.N = 'PCE@SWCNT';
DATA_20240517.S3H514R.N = 'TCE@SWCNT'; 
DATA_20240517.S4H514R.N = 'TEMED@SWCNT'; 
DATA_20240517.S5H514R.N = 'TDAE@SWCNT'; 
DATA_20240517.S6H514R.N = 'HEXADECANE@SWCNT'; 
DATA_20240517.S7H514R.N = 'DODECANE@SWCNT'; 
DATA_20241212.SR0H514R.N = 'Empty@SWCNT EMPTY'; 
DATA_20241212.SR1H514R.N = 'Water@SWCNT KVD'; 
DATA_20241212.SF6H514R.N = 'TTF@SWCNT SFF6'; 
DATA_20241212.S1BH514R.N = 'TEMED@SWCNT SF1B'; 

%%%--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%--------EXPORTING FILLED SAMPLES----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FilledB = {
            DATA_20241212.SR0H514R
            DATA_20241212.SR1H514R
            DATA_20241212.SF6H514R
            DATA_20241212.S1BH514R
            };
FilledB = FlatFieldCorrection(FilledB, DATA_20241212.FFH514R);
FilledB = clipRangeEdges(FilledB, 135, 220);
FilledB = NormalizeSample(FilledB, 135, 220);


FilledB{2}.Y = FilledB{2}.Y + 1;
FilledB{3}.Y = FilledB{3}.Y + 2;
FilledB{4}.Y = FilledB{4}.Y + 3;
plotRaman(FilledB, 0.0);
ExportMatFile(FilledB, [currentDir,'\FilledB.mat'])


%Explore the data
Filled = {
            DATA_20240517.S2H514R
            DATA_20240517.S3H514R
            DATA_20240517.S4H514R
            DATA_20240517.S5H514R
            DATA_20240517.S6H514R
            DATA_20240517.S7H514R
            };
        
Filled = FlatFieldCorrection(Filled, DATA_20240517.FF514R);
Filled = clipRangeEdges(Filled, 135, 220);
Filled = NormalizeSample(Filled, 135, 220);



%Offset
Filled{2}.Y = Filled{2}.Y + 1;
Filled{3}.Y = Filled{3}.Y + 2;
Filled{4}.Y = Filled{4}.Y + 3;
Filled{5}.Y = Filled{5}.Y + 4;
Filled{6}.Y = Filled{6}.Y + 5;

% plotRaman(Filled, 0.0);
% ExportMatFile(Filled, [currentDir,'\Filled.mat'])



%%%--------EXPORTING REFERENCES--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Order the data
    
D240517 = {
        DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.EAH514R_240517
        DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.WAH514R_240517    
            };      
D240610 = {
        DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.EAH514R_240610
        DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.WAH514R_240610
            };

%Correct the data
D240517 = FlatFieldCorrection(D240517, DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.FF514R_240517);
D240610 = FlatFieldCorrection(D240610, DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.FFH514R_240610);
%Filter the data (Mathing ranges)

D240517 = clipRangeEdges(D240517, 135, 220);
D240610 = clipRangeEdges(D240610, 135, 220);

D240517 = NormalizeSample(D240517, 135, 220);
D240610 = NormalizeSample(D240610, 135, 220);


%Plotting the data
% plotRaman(D240517, 0)
% plotRaman(D240610, 0)

%Export the data


References = {
        D240517{1}
        D240517{2}
        D240610{1}
        D240610{2}
            };

%Offset
References{2}.Y=References{2}.Y+1;
References{3}.Y=References{3}.Y+2;
References{4}.Y=References{4}.Y+3;


% plotRaman(References,0)

% ExportMatFile(References, [currentDir,'\References.mat'])
% ExportMatFile(D240610, [currentDir,'\D240610.mat'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ExportMatFile(DS, filename)

    % Initialize variables to store the output data
    numSpectra = length(DS); % Number of spectra (elements in DS)
    xdata = DS{1}.X;         % Use the .X attribute of the first structure as the x-axis
    numPoints = length(xdata); % Number of points in xData
    
    % Preallocate the matrix for efficiency
    ydata = zeros(numPoints, numSpectra); % Rows: xdata, Columns: spectra    
    labels = cell(numSpectra,1); % Rows: xdata, Columns: spectra    

    % Loop through each spectrum to collect .Y data
    for i = 1:numSpectra     
        % Assign .Y to the corresponding column
        ydata(:, i) = DS{i}.Y;
        labels{i} = DS{i}.N;
    end
    
    % Save the data matrix as a .mat file
    save(filename, 'xdata', 'ydata', 'labels');
    
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

