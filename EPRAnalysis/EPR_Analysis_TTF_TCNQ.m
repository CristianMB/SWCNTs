clc;
clear;
close all;

addpath('X:\SWCNTs\');
addpath('X:\SWCNTs\SpecialMatlabFunctions\easyspin-6.0.10\');
rootpath = 'X:\Measurements Data\EPR\';

import UsefulFunctions.*;

path_TTF_TCNQ = [rootpath,'20250806\'];

paths = {
            path_TTF_TCNQ
        };

%% Reading Data
ReadEPRfromPaths(paths);

%% Labeling

DATA_20250806.S060825A.N='Annealed P2-SWCNTs (Undoped)@2,5K';
DATA_20250806.S060825B.N='TCNQ@P2-SWCNTs (Rx6)@2,5K';
DATA_20250806.S060825C.N='TCNQ@P2-SWCNTs (Rx6)@2,5K';
DATA_20250806.S060825D.N='TCNQ@P2-SWCNTs (Rx6)@2,3K';
DATA_20250806.S060825E.N='TCNQ@P2-SWCNTs (Rx6)@5K';
DATA_20250806.S060825F.N='TCNQ@P2-SWCNTs (Rx6)@10K';
DATA_20250806.S060825G.N='TCNQ@P2-SWCNTs (Rx6)@25K';
DATA_20250806.S060825H.N='TCNQ@P2-SWCNTs (Rx6)@50K';
DATA_20250806.S060825I.N='TCNQ@P2-SWCNTs (Rx6)@75K';
DATA_20250806.S060825J.N='TCNQ@P2-SWCNTs (Rx6)@100K';
DATA_20250806.S060825K.N='TCNQ@P2-SWCNTs (Rx3)@2,5K';
DATA_20250806.S060825L.N='TTF@P2-SWCNTs (Rx6)@2,5K';
DATA_20250806.S060825M.N='TTF@P2-SWCNTs (Rx6)@5K';
DATA_20250806.S060825N.N='TTF@P2-SWCNTs (Rx6)@10K';
DATA_20250806.S0608250.N='TTF@P2-SWCNTs (Rx6)@20K';
DATA_20250806.S060825O.N='TTF@P2-SWCNTs (Rx6)@20K';
DATA_20250806.S060825P.N='TTF@P2-SWCNTs (Rx6)@50K';
DATA_20250806.S060825Q.N='Pure TCNQ (old)@2,5K';
DATA_20250806.S060825R.N='Pure TCNQ (new)@2,5K';
DATA_20250806.S060825S.N='Pure TTF (old)@2,5K';

% DATA_20250806.S060825A  %REF-AnnP2
% 
% DATA_20250806.S060825B  %R23F@2K broad
% DATA_20250806.S060825C  %R23F@2K broad
% DATA_20250806.S060825D  %R23F@2K
% DATA_20250806.S060825E  %R23F@5K
% DATA_20250806.S060825F  %R23F@10K
% DATA_20250806.S060825G  %R23F@25K
% DATA_20250806.S060825H  %R23F@50K
% DATA_20250806.S060825I  %R23F@75K
% DATA_20250806.S060825J  %R23F@100K
% DATA_20250806.S060825K  %R23C@2K
% 
% DATA_20250806.S060825L  %R13@2K
% DATA_20250806.S060825M  %R13@5K
% DATA_20250806.S060825N  %R13@5K
% DATA_20250806.S060825O  %R13@20K
% DATA_20250806.S060825P  %R13@55K
% 
% DATA_20250806.S060825Q  %TCNQ (old)
% DATA_20250806.S060825R  %TCNQ (new)
% DATA_20250806.S060825S  %TTF (old)

%% Parameters for corrections

% Rescaling for BG substraction

DATA_20250806.S060825B.BS = 1;  %R23F
DATA_20250806.S060825C.BS = 1;  %R23F
DATA_20250806.S060825D.BS = 1;  %R23F
DATA_20250806.S060825E.BS = 1;  %R23F
DATA_20250806.S060825F.BS = 1;  %R23F
DATA_20250806.S060825G.BS = 1;  %R23F
DATA_20250806.S060825H.BS = 1;  %R23F
DATA_20250806.S060825I.BS = 1;  %R23F
DATA_20250806.S060825J.BS = 1;  %R23F
DATA_20250806.S060825K.BS = 1;  %R23C
 
DATA_20250806.S060825L.BS = 1;  %R13
DATA_20250806.S060825M.BS = 1;  %R13
DATA_20250806.S060825N.BS = 1;  %R13
DATA_20250806.S0608250.BS = 1;  %R13
DATA_20250806.S060825O.BS = 1;  %R13
DATA_20250806.S060825P.BS = 1;  %R13

DATA_20250806.S060825Q.BS = 1;  %TCNQ (old)
DATA_20250806.S060825R.BS = 1;  %TCNQ (new)

% Rescaling factors
 
DATA_20250806.S060825B.RS = 1;  %R23F
DATA_20250806.S060825C.RS = 1;  %R23F
DATA_20250806.S060825D.RS = 1;  %R23F
DATA_20250806.S060825E.RS = 1;  %R23F
DATA_20250806.S060825F.RS = 1;  %R23F
DATA_20250806.S060825G.RS = 1;  %R23F
DATA_20250806.S060825H.RS = 1;  %R23F
DATA_20250806.S060825I.RS = 1;  %R23F
DATA_20250806.S060825J.RS = 1;  %R23F
DATA_20250806.S060825K.RS = 1;  %R23C
 
DATA_20250806.S060825L.RS = 1;  %R13
DATA_20250806.S060825M.RS = 1;  %R13
DATA_20250806.S060825N.RS = 1;  %R13
DATA_20250806.S0608250.RS = 1;  %R13
DATA_20250806.S060825O.RS = 1;  %R13
DATA_20250806.S060825P.RS = 1;  %R13

DATA_20250806.S060825Q.RS = 1;  %TCNQ (old)
DATA_20250806.S060825R.RS = 1;  %TCNQ (new)

%% Sample selection and plotting

EPR_TCNQ= {
%             DATA_20250806.S060825A  %REF-AnnP2
            
                DATA_20250806.S060825D  %R23F@2K
                DATA_20250806.S060825E  %R23F@5K
                DATA_20250806.S060825F  %R23F@10K
                DATA_20250806.S060825G  %R23F@25K
                DATA_20250806.S060825H  %R23F@50K
                DATA_20250806.S060825I  %R23F@75K
                DATA_20250806.S060825J  %R23F@100K

%             DATA_20250806.S060825K  %R23C@2K
            
%             DATA_20250806.S060825Q  %TCNQ (old)
%             DATA_20250806.S060825R  %TCNQ (new)
        };

EPR_TTF= {
%             DATA_20250806.S060825A  %REF-AnnP2
                   
            DATA_20250806.S060825L  %R13@2K
            DATA_20250806.S060825M  %R13@5K
            DATA_20250806.S060825N  %R13@5K
            DATA_20250806.S060825O  %R13@20K
            DATA_20250806.S060825P  %R13@55K

%             DATA_20250806.S060825S  %TTF (old)        
        };

%% X Axis correction

freq=9.44*10^9; %Hz

DATA_20250806.S060825A.X=DATA_20250806.S060825A.X.*freq/DATA_20250806.S060825A.PAR.MWFQ;
EPR_TCNQ = CorrectEPRXAxis(EPR_TCNQ,freq);
EPR_TTF = CorrectEPRXAxis(EPR_TTF,freq);

%% Plotting and corrections for plotting

EPR_TTF = FilterDataByXRange(EPR_TTF, 3220, 3520); 
EPR_TTF = EPR_BG_Correction(EPR_TTF, DATA_20250806.S060825A);
EPR_TTF = RemovePolyBG(EPR_TTF, 0);
% EPR_TTF = Normalize(EPR_TTF,3300,3400, 'M');

plotEPR(EPR_TTF, 0)

EPR_TCNQ = FilterDataByXRange(EPR_TCNQ, 3220, 3520); 
EPR_TCNQ = EPR_BG_Correction(EPR_TCNQ, DATA_20250806.S060825A);
EPR_TCNQ = RemovePolyBG(EPR_TCNQ, 0);
% EPR_TCNQ = Normalize(EPR_TCNQ,3300,3400, 'M');

plotEPR(EPR_TCNQ, 0)


% EPR_All = [EPR_TCNQ(:); EPR_TTF(:)];  % forces both to column vectors
% plotEPR(EPR_All, 0)

%% Functions in development for EPR

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

function DSList = EPR_BG_Correction(DSList, BackgroundSample)

    numSamples = length(DSList);

    bgX = BackgroundSample.X;
    bgY = BackgroundSample.Y;

    for i = 1:numSamples
        sample = DSList{i};
        X = sample.X;
        Y = sample.Y;
        BS = sample.BS;
        
        % Match background to sample X via interpolation
        validIndices = bgX >= min(X) & bgX <= max(X);
        bgY_interp = interp1(bgX(validIndices), bgY(validIndices), X, 'linear', 'extrap');
        
        % Subtract background and apply scaling
        DSList{i}.Y = (Y - BS*bgY_interp);       
    end
end

function plotEPR(SamplesToPlot, offset)
% Plot corrected EPR spectra with optional background subtraction and rescaling
%
% Inputs:
%   SamplesToPlot    - Cell array of EPR spectrum structures with .X, .Y, .N
%   BackgroundSample - A single spectrum structure with fields .X and .Y
%   BackgroundScale  - Cell array of scalars to scale the background per sample
%   Rescale          - Cell array of scalars to rescale each spectrum after subtraction
%
% If BackgroundScale or Rescale are empty ({}), defaults to {1, 1, ..., 1}
% If BackgroundSample is empty ([]), background subtraction is skipped.

    numSamples = length(SamplesToPlot);

    figure;
    hold on;

    for i = 1:numSamples
        sample = SamplesToPlot{i};
        X = sample.X;
        Y = sample.Y;
        N = sample.N;
        RS = sample.RS;
        
        % Match background to sample X via interpolation

        % Subtract background and apply scaling
        Y = Y * RS - offset * i;
        
        plot(X, Y, 'DisplayName', N, 'LineWidth', 1.3);
    end
%     plot(bgX, bgY, 'DisplayName', N, 'LineWidth', 1.3);

    xlabel('Magnetic Field (G)', 'FontSize', 14);
    ylabel('Intensity (a.u.)', 'FontSize', 14);
    title('EPR Spectra');
    legend('show', 'FontSize', 11);
    grid on;
    hold off;
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

function correctedEPR = CorrectEPRXAxis(EPRcell, freq)
% Applies X-axis correction for all spectra in the EPRcell array
% X_corrected = X * freq / par.MWFQ
% Input:
%   EPRcell - cell array of EPR spectra (structures with fields X and par)
%   freq    - experimental microwave frequency in Hz
% Output:
%   correctedEPR - corrected cell array (same format)

    correctedEPR = EPRcell; % Initialize output
    
    for i = 1:length(EPRcell)
        sample = EPRcell{i};
        sample.X = sample.X * freq / sample.PAR.MWFQ;
        correctedEPR{i} = sample;
    end
end

function dataStructures = ReadEPRfromPaths(paths)

        dataStructures = struct();  
        for p = 1:length(paths)
            folderPath = paths{p};
            dataset_name = strsplit(paths{p}, "\");
            fieldName = ['DATA_', strrep(dataset_name{end-1}, '.', '')];
            dirInfo = dir(fullfile(paths{p}, '*.DSC'));
            fileList = {dirInfo(~[dirInfo.isdir]).name};
            structure = struct();

            for f = 1:length(fileList)

                fileName = fileList{f};
                fullFilePath = fullfile(folderPath, fileName);
                
                sampleName = upper(strrep(fileList{f}, '.DSC',''));
                
                [structure.(sampleName).X,structure.(sampleName).Y,structure.(sampleName).PAR] = eprload(fullFilePath);
                 
                structure.(sampleName).N = structure.(sampleName).PAR.TITL;
                
                structure.(sampleName).BS = 1;  %BG Rescale factor (for substracting BG)
                structure.(sampleName).RS = 1;  %Rescale factor (for matching)
            end
           dataStructures.(fieldName) = structure;
           assignin('caller', fieldName, structure); % Assign data to a variable in the caller workspac
        end
end
         


