bl_path = 'C:\Users\cborja\OneDrive - Universiteit Antwerpen\Measurements\Absorption\Baselines.csv';
path20240216 = 'C:\Users\cborja\OneDrive - Universiteit Antwerpen\Measurements\Absorption\20240216\DGUC.csv';
path20240304 = 'C:\Users\cborja\OneDrive - Universiteit Antwerpen\Measurements\Absorption\20240304\S6S7.csv';

% Read data from CSV files


%DATA = readSamplesData(path1)

BL = readSamplesData(bl_path)

DATA1 = readSamplesData(path20240216); % Skip the first two rows
Samples1 = labelSamples(DATA1)

toPlot1 = {
    'S2DGUC',
    'S3DGUC',
    'S4DDGUC',
    %'S4DDGUC2',
    %'S4LDGUC'
    }

DATA2 = readSamplesData(path20240304); % Skip the first two rows
Samples2 = labelSamples(DATA2)

toPlot2 = {
    'S6',
    'S7',
    }

%For normalization 
LL = 600;
UL = 800;

%computeIntegral(DATA.S2DGUC, LL, UL)
DATANormal1 = Normalize(DATA1, toPlot1, LL, UL)
DATANormal2 = Normalize(DATA2, toPlot2, LL, UL)

%CORRECTIONS
%DATANormal.S2DGUC.A = DATANormal.S2DGUC.A-BL.H2OinD2O.A*0.003

%plotSamples(DATA, 'All', Samples)
%plotSamples(DATA, toPlot, Samples)
%plotSamples(DATANormal1, toPlot1, Samples1)
%plotSamples(DATANormal2, toPlot2, Samples2)


DATA = DATANormal1
DATA.S6 = DATANormal2.S6
DATA.S7 = DATANormal2.S7

Samples = labelSamples(DATA)
Samples('S2DGUC') = 'Tetrachloroethylene@SWCNT DGU C'
Samples('S3DGUC') = 'Trichloroethylene@SWCNT DGU C'
Samples('S4DDGUC') = 'TEMED@SWCNT DGU C'
Samples('S6') = 'Hexadecane@SWCNT'
Samples('S7') = 'Dodecane@SWCNT'

toPlot = {
    'S2DGUC',
    'S3DGUC',
    'S4DDGUC',
    'S6',
    'S7'
    }

plotSamples(DATA, toPlot, Samples)
plotSamples(flattenSpectra(DATA, [370,600,800,1500]), toPlot, Samples)

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

function flattenedData = flattenSpectra(DataStructure, points)
    % Get the fieldnames of the data structure
    sampleNames = fieldnames(DataStructure);

    % Loop through each sample in the data structure
    for i = 1:numel(sampleNames)
        sample = sampleNames{i};

        % Get the wavelength and absorption data for the current sample
        selected_A = DataStructure.(sample).A(ismember(DataStructure.(sample).W, points));
        selected_W = DataStructure.(sample).W(ismember(DataStructure.(sample).W, points));

        % Define the objective function for least-squares optimization
        objective = @(a) sum((selected_A - a ./ selected_W).^2);

        % Choose initial value of 'a'
        initialA = 1; % or any initial value

        % Perform least-squares optimization to find the optimal value of 'a'
        optimalA = fmincon(objective, initialA, [], [], [], [], [], [], []);
        optimalA
        % Subtract the background from the absorption data within the specified range
        A_flattened = DataStructure.(sample).A - optimalA ./ DataStructure.(sample).W;
        
        % Update the absorption data in the data structure
        
        DataStructure.(sample).A = A_flattened;
    end
    
    % Return the modified data structure
    flattenedData = DataStructure;
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

