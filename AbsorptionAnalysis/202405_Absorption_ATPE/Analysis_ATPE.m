clc;
clear;
import UsefulFunctions.*;
rootpath = 'X:\Measurements Data\Absorption\';
addpath('X:\SWCNTs\SpecialMatlabFunctions\DrosteEffect-BrewerMap-3.2.5.0')
addpath('X:\SWCNTs')


%All paths as default

Refs= [rootpath,'References.csv'];

ATPE_2= [rootpath,'20250515\ATPE_Round2_Part2.csv'];
ATPE_1a= [rootpath,'20250425\ATPE_T9_MimicBL.csv'];
ATPE_1b= [rootpath,'20250425\ATPE_T9_TopPhases.csv'];

%Select the paths of interest
paths = {
    Refs
    ATPE_1a
    ATPE_1b
    ATPE_2
        };
    

%Read and structure data from the paths
ReadAbsorptionFromPaths(paths);


DATA_20250515.Parent_30.N = 'Parent SWCNTs dispersion';
DATA_20250515.T5_2MT_1.N = 'Semiconducting enriched sample';

AllSpecs = {
                DATA_20250515.Parent_30
                DATA_20250515.T5_2MT_1

%                 DATA_20250515.T1_1
%                 DATA_20250515.T2_1
%                 DATA_20250515.T3_1
%                 DATA_20250515.T4_MT_20
%                 DATA_20250515.T6_3MT_1
%                 DATA_20250515.T7_MTT_1
%                 DATA_20250515.T8_MT2T_1
%                 DATA_20250515.B6_3MB_30
%                 DATA_20250515.B8_MT2B_1
                
       };
   
AllSpecs = FilterDataByXRange(AllSpecs,200,1270)  
AllSpecs = Normalize(AllSpecs,240, 260, 'M')  

AllSpecs = SubtractInverseBG(AllSpecs, [616, 1270])

% plotAbsorption(AllSpecs, -0.2);
%    
% Abs_303 = get_Y_values(AllSpecs, 303)
% Abs_315 = get_Y_values(AllSpecs, 315)
% Abs_362 = get_Y_values(AllSpecs, 362)
% Abs_447 = get_Y_values(AllSpecs, 447)



plotAbsorption(AllSpecs, 0.15);


function Y_values = get_Y_values(D_list, value)
    N = numel(D_list);
    Y_values = nan(1, N);

    for k = 1:N
        DATA = D_list{k};
        idx = find(DATA.X == value, 1); % Find index of value in X
        if ~isempty(idx)
            Y_values(k) = DATA.Y(idx);  % Get corresponding Y
        end
    end
end

function filteredSamples = FilterDataByXRange(samplesToFilter, xMin, xMax)
    % FilterDataByXRange filters the data of each sample to include only the points within the specified X-range.
    %
    % Inputs:
    %   - samplesToFilter: Cell array of structures, each with fields 'X' and 'Y'.
    %   - xMin: The minimum value of X to include in the filtered data.
    %   - xMax: The maximum value of X to include in the filtered data.
    % Outputs:
    %   - filteredSamples: Cell array of structures with filtered 'X' and 'Y' values within the range [xMin, xMax].

    filteredSamples = cell(size(samplesToFilter));
    
    % Iterate over each sample to filter
    for sampleIdx = 1:length(samplesToFilter)
        currentSample = samplesToFilter{sampleIdx};
        
        % Find the indices of X-values within the specified range
        validIndices = currentSample.X >= xMin & currentSample.X <= xMax;
        
        % Filter the data based on valid indices
        filteredSample = currentSample;
        filteredSample.X = currentSample.X(validIndices);
        filteredSample.Y = currentSample.Y(validIndices);
        
        % Store the filtered sample
        filteredSamples{sampleIdx} = filteredSample;
    end
end
