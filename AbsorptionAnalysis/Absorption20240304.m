bl_path = 'C:\Users\cborja\OneDrive - Universiteit Antwerpen\Measurements\Absorption\Baselines.csv';
path = 'C:\Users\cborja\OneDrive - Universiteit Antwerpen\Measurements\Absorption\20240216\DGUC.csv';

% Read data from CSV files


%DATA = readSamplesData(path1)

BL = readSamplesData(bl_path)

DATA = readSamplesData(path); % Skip the first two rows
Samples = labelSamples(DATA)

toPlot = {
    'S2DGUC',
    %'S3DGUC',
    %'S4DDGUC',
    %'S4DDGUC2',
    %'S4LDGUC'
    }
%For normalization 
LL = 1700;
UL = 1900;

%computeIntegral(DATA.S2DGUC, LL, UL)
DATANormal = Normalize(DATA, toPlot, LL, UL)

%CORRECTIONS
DATANormal.S2DGUC.A = DATANormal.S2DGUC.A-BL.H2OinD2O.A*0.003

%plotSamples(DATA, 'All', Samples)
%plotSamples(DATA, toPlot, Samples)
plotSamples(DATANormal, toPlot, Samples)


function samples = readSamplesData(filePath)
    % Read the header
    header = readcell(filePath, 'Range', 'A1:Z1');
    header = cellfun(@(x) strrep(x, ' 100%T', ''), header, 'UniformOutput', false);
    sampleNames = header(1, 1:2:end);
    
    % Read ONLY the datalines
    data = readmatrix(filePath, 'Range', ['A' num2str(3) ':Z' num2str(2328)]);
    
    samples = struct();

    % Iterate through each sample and store wavelength and absorption data
    for i = 1:length(sampleNames)
        % Extract wavelength and absorption data for the current sample
        wavelengths = data(:, 2*i - 1); % Odd columns contain wavelength
        absorption = data(:, 2*i);       % Even columns contain absorption

        % Store the data in the container object with the sample name
        sampleName = sampleNames{i};
        samples.(sampleName).W = wavelengths;
        samples.(sampleName).A = absorption;
    end
end

function labeling = labelSamples(DataStructure)
    labeling = containers.Map;
    sampleNames = fieldnames(DataStructure);
    for i = 1:numel(sampleNames)
        sample = sampleNames{i};
        labeling(sample) = sample;
    end
end

function integralValue = computeIntegral(sample, lowerLimit, upperLimit)
    % Check if the sampleName is in DATA
    % Extract X and Y values for the specified sample
    x = sample.W;
    y = sample.A;
    % Define the function to integrate
    f = @(xi) interp1(x, y, xi, 'pchip');

    % Calculate the integral
    integralValue = integral(f, lowerLimit, upperLimit);
end

function normalizedData = Normalize(DataStructure, samplesToNormalize, lowerLimit, upperLimit)
    % Initialize the container for normalized data
    normalizedData = DataStructure;

    % Iterate over each sample to be normalized
    for sampleIdx = 1:length(samplesToNormalize)
        currentSample = samplesToNormalize{sampleIdx};
        if isfield(DataStructure, currentSample) && isfield(DataStructure.(currentSample), 'A')
            normalizedData.(currentSample).A = normalizedData.(currentSample).A/computeIntegral(normalizedData.(currentSample),lowerLimit, upperLimit);
        end
    end
end

function plotSamples(DATA, samplesToPlot, nameLabelMap)

    if strcmp(samplesToPlot, 'All')
        samplesToPlot=keys(nameLabelMap);
        %Exclude Baseline when selecting All
        samplesToPlot=samplesToPlot(2:end);
    end
    
    % Create a figure for the plot
    figure;
  
    % Iterate over each sample
    for sampleIdx = 1:length(samplesToPlot)
        currentSample = samplesToPlot{sampleIdx};
        if isfield(DATA, currentSample) && isKey(nameLabelMap, currentSample)

            % Get the current sample, X values, and Y values
            currentX = DATA.(currentSample).W;
            currentY = DATA.(currentSample).A;

            if ismember(currentSample, samplesToPlot)
                plot(currentX, currentY, 'DisplayName', nameLabelMap(currentSample));
                hold on; % Add spectra to the same plot
            end
            hold on; % Add spectra to the same plot
        end
    end

    % Add labels and legend

    xlabel('Wavelenght (nm)');
    ylabel('Absorption (a.u.)');
    title('Absorption Spectra');
    legend('show');

    % Optional: Customize the plot further if needed
    grid on;
    % Hold off to stop adding new plots to the current figure
    hold off;
    
end

