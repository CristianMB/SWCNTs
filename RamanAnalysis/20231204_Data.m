
%addpath('C:\Users\cborja\OneDrive - Universiteit Antwerpen\SWCNTs\RamanAnalysis\');
addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\RamanAnalysis\');
mainPath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements\Raman\20240111\';

dirInfo = dir(fullfile(mainPath, '*.*'));
% Exclude directories from the list
fileList = {dirInfo(~[dirInfo.isdir]).name};

DATA = getData(fileList,'LD', 514, 1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      Manually Adjust Metadata

%%% Mode Adjustment
%DATA.FL568A.M = 'HD';

%%% Mode Adjustment
%DATA.FL568A.M = 'HD';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      Normalizations

%DATA.BS8568CNORMED = DATA.BS8568C;
%DATA.BS8568CNORMED.Y = (DATA.BS8568CNORMED.Y ./ DATA.FL568A.Y);

%CCLPeakIntegral = computeIntegral(DATA.CCL4568B, 295, 330);


%for sampleIdx = 1:numSamples
%    currentSample = Samples{sampleIdx};
%    DATA.(currentSample).Y = DATA.(currentSample).Y / DATA.(currentSample).T;
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      Plot Preparation
%sampleNames = fieldnames(DATA);

%samplesToPlot = {
     %'ARC514D',
%     'BS8514A',
%     'BS8514C',
     
      %'ARC514A',       %84nm seems to be a reference instead of a sample
      %'DUMMY',
      %'L514A', 
      %'FLAT514H',      %108nm
%    };

%sampleNames = fieldnames(DATA);
%plotSamples(DATA, samplesToPlot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      Read and pre-process Data

function DATA = getData(fileList,default_mode, laser_wl, time)
    % Use the dir function to list files
    % Display the list of files
    Samples = cell(1, length(fileList));
    X_values = cell(1, length(fileList));
    Y_values = cell(1, length(fileList));

    for f = 1:length(fileList)
        
        raw_spectrum = RdExp([mainPath, fileList{f}]);
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
    
    
    numSamples = numel(Samples);
    DATA = struct();
    
    for sampleIdx = 1:numSamples
        %Read X,Y data for each sample
        currentSample = Samples{sampleIdx};
        DATA.(currentSample).X = X_values{sampleIdx};
        DATA.(currentSample).Y = Y_values{sampleIdx};

        %Metadata M is Mode, L is the laser WL, T is the Exposure Time
        DATA.(currentSample).M = default_mode;
        DATA.(currentSample).L = laser_wl;
        DATA.(currentSample).T = time;
    end

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

function plotSamples(DATA, samplesToPlot)

    % Create a figure for the plot
    figure;

    % Iterate over each sample
    for sampleIdx = 1:length(samplesToPlot)
        currentSample = samplesToPlot{sampleIdx};
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




