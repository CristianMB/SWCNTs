
addpath('C:\Users\cborja\OneDrive - Universiteit Antwerpen\SWCNTs\RamanAnalysis\');
mainPath = 'C:\Users\cborja\OneDrive - Universiteit Antwerpen\Measurements\Raman\20240111\';
dirInfo = dir(fullfile(mainPath, '*.*'));
% Exclude directories from the list
fileList = {dirInfo(~[dirInfo.isdir]).name};


DATA = getData(fileList, 1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      SAMPLES

samplesReferences = {
    %'FLATHD',
    %'LL514',
    %'LL514HD'
    'FLATHDNORMED',
    };

samplesGBand = {
    'S240111A',
    'S240111B',
    'S240111C',
    'S240111D',
    'S240111E',
    'S240111F',
    'S240111G',
    'S240111H',
    'S240111I'
    };

samplesRBM = {
    'S240111J',
    'S240111K',
    'S240111KK'
    'S240111L',
    'S240111M',
    'S240111N',
    'S240111O',
    'S240111P',
    'S240111Q',
    'S240111R',
    'S240111S',
    };


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   DATA ADJUSTMENT
%%% TIME %%%
DATA.S240111A.T = 10;
DATA.S240111B.T = 10;
DATA.S240111BB.T = 10;
DATA.S240111C.T = 20;
DATA.S240111D.T = 15;
DATA.S240111E.T = 10;
DATA.S240111F.T = 10;
DATA.S240111G.T = 10;
DATA.S240111H.T = 10;
DATA.S240111I.T = 15;
DATA.S240111J.T = 90;
DATA.S240111K.T = 50;
DATA.S240111KK.T = 90;
DATA.S240111L.T = 90;
DATA.S240111M.T = 90;
DATA.S240111N.T = 90;
DATA.S240111O.T = 90;
DATA.S240111P.T = 90;
DATA.S240111Q.T = 90;
DATA.S240111R.T = 90;
DATA.S240111S.T = 30;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     DATA NORMALIZATION

FFNORM = computeIntegral(DATA.FLATHD, 120, 230);
DATA.FLATHDNORMED = DATA.FLATHD;
DATA.FLATHDNORMED.Y = DATA.FLATHD.Y/FFNORM;

%Write a function to normalize by flatfield all the samples in a list.
%similar to how we subselect samples in the plotting function. Normalize to
%a flat field

%DATA.BS8568CNORMED = DATA.BS8568C;
%DATA.BS8568CNORMED.Y = (DATA.BS8568CNORMED.Y./ DATA.FL568ANORMED.Y) ;





plotSamples(DATA, samplesReferences)
%plotSamples(DATA, samplesRBM)
%plotSamples(DATA, samplesGBand)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     READ AND PROCESS FUNCTIONS

function DATA = getData(fileList, time)

    % Use the dir function to list files
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
    
    
    numSamples = numel(Samples);
    DATA = struct();
    
    for sampleIdx = 1:numSamples
    %Read X,Y data for each sample
    currentSample = Samples{sampleIdx};
    DATA.(currentSample).X = X_values{sampleIdx};
    DATA.(currentSample).Y = Y_values{sampleIdx};
    
    %Metadata T is the Exposure Time
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