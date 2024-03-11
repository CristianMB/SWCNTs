
path1 = 'C:\Users\cborja\OneDrive - Universiteit Antwerpen\Measurements\Absorption\20240212\AllDataCentrifuge.csv';
path2 = 'C:\Users\cborja\OneDrive - Universiteit Antwerpen\Measurements\Absorption\20240212\AllDataS5.csv';

% Read data from CSV files


%DATA = readSamplesData(path1)


DATA_CF = readSamplesData(path1); % Skip the first two rows
nameLabelMapCF = containers.Map;
% Add mappings to the dictionary
nameLabelMapCF('Baseline') = 'Baseline';
nameLabelMapCF('S2') = 'S2';
nameLabelMapCF('S3') = 'S3';
nameLabelMapCF('S4') = 'S4';
PlotCF = {
    'S2',
    'S3',
    'S4'
    };

DATA_S5 = readSamplesData(path2); % Skip the first two rows
nameLabelMapS5 = containers.Map;
% Add mappings to the dictionary
nameLabelMapS5('Baseline') = 'Baseline';
nameLabelMapS5('S5R1') = 'S5R1';
nameLabelMapS5('S5R2') = 'S5R2';
nameLabelMapS5('S5R3') = 'S5R3';
nameLabelMapS5('S5R4') = 'S5R4';

PlotS5 = {
    'S5R1',
    'S5R2',
    'S5R3',
    'S5R4'
    };

%plotSamples(DATA_S5, PlotS5, nameLabelMapS5)
plotSamples(DATA_CF, PlotCF, nameLabelMapCF)



function samples = readSamplesData(filePath)
    % Read the data from the CSV file
    data = readmatrix(filePath); % Skip the first two rows
    header = readcell(filePath, 'Range', 'A1:Z1');
    header = cellfun(@(x) strrep(x, ' 100%T', ''), header, 'UniformOutput', false);
    sampleNames = header(1, 1:2:end);

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


function plotSamples(DATA, samplesToPlot, nameLabelMap)

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

