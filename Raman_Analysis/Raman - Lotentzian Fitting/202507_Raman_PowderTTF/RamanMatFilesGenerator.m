clc;
clear;
addpath('X:\SWCNTs');
% addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\');
import UsefulFunctions.*;

% rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements Data\Absorption\';
rootpath = 'X:\Measurements Data\Raman\';

%All paths as default
path_a = [rootpath,'20250618\'];
path_b = [rootpath,'20250520\'];


%Select the paths of interest

paths = {
    path_a
    path_b
    };


ReadRamanFromPaths(paths,2);




%%%--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %--------SAMPLE COMPARISION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Explore the data
Data_G = {
                    DATA_20250618.PAPL514G
                    DATA_20250618.PP2L514G
%                     DATA_20250520.PP2514G
%                     DATA_20250618.U18L514G
                    DATA_20250618.U20L514G
                    DATA_20250618.R20A514G
                    DATA_20250618.R20B514G
                    DATA_20250618.R20C514G
                    DATA_20250618.R20D514G
                    DATA_20250618.R20E514G
}

%Correct the data (Flat field)
 

%Filter the data (Mathing ranges)

Data_G = clipRangeEdges(Data_G, 1250, 1680);
Data_G = RemovePolyBG(Data_G, 0);
Data_G = SubstractLinearBG(Data_G, 1250, 1680);
Data_G = Normalize(Data_G, 1500, 1680, 'M');

% D241212 = NormalizeSample(D241212, 135, 220)

%Plotting the data
Data_G = OffsetSamples(Data_G,1)
plotRaman(Data_G, 0.0)

%Export the data

currentDir = pwd; 

%Offset



% plotRaman(AllData,0)

ExportMatFile(Data_G, [currentDir,'\Data_G.mat'])

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


function DSList = RemovePolyBG(DSList, degree)
    % Remove baseline from a list of data structures using polynomial fitting
    % DSList: list of structures, each with fields X (Raman shift) and Y (intensity)
    % degree: Degree of the polynomial used for baseline fitting
    
    % Iterate over each structure in the list
    for i = 1:length(DSList)
        DS = DSList{i};  % Extract the current data structure
        
        % Extract the X and Y data
        X = DS.X;  % Raman shift (assumed centered at zero)
        Y = DS.Y;  % Intensity values
        
        % Identify regions to exclude based on peak detection
        % You can implement your own peak detection logic here or use findpeaks
        [pks, locs] = findpeaks(Y, 'MinPeakHeight', 0.05, 'MinPeakDistance', 10);
        
        % Create a mask for excluding the peak regions
        exclude_indices = false(size(Y));
        exclude_indices(locs) = true;

        % Fit a polynomial to the non-peak regions
        p = polyfit(X(~exclude_indices), Y(~exclude_indices), degree);  % Polynomial coefficients
        baseline = polyval(p, X);  % Evaluate the polynomial to get the baseline

        % Subtract the baseline from the original intensity
        Y_corrected = Y - baseline;

        % Find the minimum value of the corrected spectrum
        min_val = min(Y_corrected);

        % Shift the corrected spectrum so that its minimum value is zero
        Y_corrected = Y_corrected - min_val;

        % Update the structure with the corrected Y values
        DS.Y = Y_corrected;
        
        % Optionally, display the polynomial coefficients for debugging
%         disp(['Structure ', num2str(i), ' Polynomial Coefficients: ', num2str(p)]);
        
        % Save the updated structure back to the list
        DSList{i} = DS;
    end
end

function DSList = OffsetSamples(DSList, offset)

        % Get a ColorBrewer colormap (e.g., 'Set1', 'Dark2', etc.)
        numSamples = length(DSList);  % Number of samples to plot
        for sampleIdx = 1:numSamples
            DSList{sampleIdx}.Y = DSList{sampleIdx}.Y -offset * sampleIdx;
        end
end
