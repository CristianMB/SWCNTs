
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      Setting Folders

date = '20231213\';
path = 'C:\Users\cborja\OneDrive - Universiteit Antwerpen\Measurements\Raman\';
mainPath = [path, date];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      Read and pre-process Data

[Samples, X_values, Y_values] = listFiles(mainPath);
numSamples = numel(Samples);
DATA = struct();

for sampleIdx = 1:numSamples
    %Read X,Y data for each sample
    currentSample = Samples{sampleIdx};
    DATA.(currentSample).X = X_values{sampleIdx};
    DATA.(currentSample).Y = Y_values{sampleIdx};
    
    %Metadata M is Mode, L is the laser WL, T is the Exposure Time
    DATA.(currentSample).M = 'LD';
    DATA.(currentSample).L = '568';
    DATA.(currentSample).T = 10;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      Manually Adjust Metadata
%%% Mode Adjustment
DATA.BS8568C.M = 'HD';
DATA.FL568A.M = 'HD';

%%% Mode Adjustment
DATA.ARC568A.T = 5;
DATA.FL568A.T = 100;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      Normalizations

DATA.BS8568CNORMED = DATA.BS8568C;
DATA.BS8568CNORMED.Y = (DATA.BS8568CNORMED.Y ./ DATA.FL568A.Y);

CCLPeakIntegral = computeIntegral(DATA.CCL4568B, 295, 330);



for sampleIdx = 1:numSamples
    currentSample = Samples{sampleIdx};
    DATA.(currentSample).Y = DATA.(currentSample).Y / DATA.(currentSample).T;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      Plot Preparation

sampleNames = fieldnames(DATA);

samplesToPlot20231204 = {
    'ARC568A',
    'BS8568A',
    'BS8568B',
    %'BS8568C',
    %'BS8568CNORMED',
    %'CCL4568B',
    %'FL568A',
    %'L568A' ,
    %'L568B',
    %'L568C',
    %'L568D'
    };

samplesToPlot20231213 = {
    'ARC568A',
    'BS8568A',
    'BS8568B',
    %'BS8568C',
    %'BS8568CNORMED',
    %'CCL4568B',
    %'FL568A',
    %'L568A' ,
    %'L568B',
    %'L568C',
    %'L568D'
    };


samplesToPlot = samplesToPlot20231204;

plotting = true;


if plotting==true
    % Create a figure for the plot
    figure;

    % Iterate over each sample
    for sampleIdx = 1:length(sampleNames)
        currentSample = sampleNames{sampleIdx};
        % Get the current sample, X values, and Y values
        currentX = DATA.(currentSample).X;
        currentY = DATA.(currentSample).Y;

        if ismember(currentSample, samplesToPlot)
            plot(currentX, currentY, 'DisplayName', currentSample);
            hold on; % Add spectra to the same plot
        end
        hold on; % Add spectra to the same plot
    end

    % Add labels and legend
    
    xlabel('Raman Shift (cm^{-1})');
    ylabel('Intensity (a.u.)');
    title('Raman Spectra');
    legend('show');

    % Optional: Customize the plot further if needed
    grid on;
    % Hold off to stop adding new plots to the current figure
    hold off;
end

function integralValue = computeIntegral(sample, lowerLimit, upperLimit)
    % Check if the sampleName is in DATA
    % Extract X and Y values for the specified sample
    x = sample.X;
    y = sample.Y;
    % Define the function to integrate
    f = @(xi) interp1(x, y, xi, 'pchip');

    % Calculate the integral
    integralValue = integral(f, lowerLimit, upperLimit);
end

function [Samples, X_values, Y_values]  = listFiles(mainPath)
    % Use the dir function to list files
    dirInfo = dir(fullfile(mainPath, '*.*'));
    % Exclude directories from the list
    fileList = {dirInfo(~[dirInfo.isdir]).name};
    
    % Display the list of files
    Samples = cell(1, length(fileList));
    X_values = cell(1, length(fileList));
    Y_values = cell(1, length(fileList));
    
    for f = 1:length(fileList)
        raw_spectrum = RdExp(fileList{f});
        raw_spectrum(:,2)=[];
        
        NumSpec=length(raw_spectrum(1,:))-1;
        NumDel=1;
        
        X=raw_spectrum(:,1);
        for i=1:1024
            spectrum=sort(raw_spectrum(i,2:NumSpec+1));
            Y(:,i)= mean(spectrum(NumDel+1:NumSpec-NumDel));  
        end
        
        % Store X and Y values for the current file
        [~, sample, ~] = fileparts(upper(fileList{f}));
        
        Samples{f} = sample;
        X_values{f} = X;
        Y_values{f} = Y;
        
    end  
end


    
  