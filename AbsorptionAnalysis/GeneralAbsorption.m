clc;
clear;

path_baselines = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements\Absorption\References.csv';
path_20231206 = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements\Absorption\20231206\alldata20231206.csv';
path_20240117a = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements\Absorption\20240117\S1Data20240117.csv';
path_20240117b = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements\Absorption\20240117\S2Data20240117.csv';
path_20240124 = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements\Absorption\20240124\AllDataS3S4.csv';
path_20240129 = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements\Absorption\20240129\AllData.csv';
path_20240202 = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements\Absorption\20240202\AllData.csv';
path_20240212a = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements\Absorption\20240212\AllDataCentrifuge.csv';
path_20240212b = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements\Absorption\20240212\AllDataS5.csv';
path_20240214 = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements\Absorption\20240214\DGU.csv';
path_20240215 = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements\Absorption\20240215\Dialysis.csv';
path_20240216 = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements\Absorption\20240216\DGUC.csv';
path_20240220 = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements\Absorption\20240220\RinsingS6S7.csv';
path_20240304 = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements\Absorption\20240304\S6S7.csv';
path_20240305 = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements\Absorption\20240305\S5.csv';
path_20240307 = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements\Absorption\20240307\DGUS5S6S7.csv';
path_20240308 = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements\Absorption\20240308\DialS5S6S7.csv';


paths = {
        path_baselines,
        path_20231206,
        %path_20240117a,
        %path_20240117b,
        %path_20240124,
        %path_20240129,
        path_20240202,
        path_20240212a,
        %path_20240212b,
        path_20240214,
        %path_20240215,
        path_20240216,
        %path_20240220,
        path_20240304,
        path_20240305,
        path_20240307,
        path_20240308
        };

%Read DATA
ReadFromPaths(paths);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uncorrected Data Empty Tubes
DATA_20231206.SFF4dil.N = 'SF Methanol@SWCNT A';
DATA_20231206.SFF4dil2.N = 'SF Methanol@SWCNT B';
DATA_20231206.SFF5dil.N = 'SF SFF5 Water@SWCNT';
DATA_20240214.S2DB.N = 'CB Empty@SWCNT DGU B S2';
DATA_20240214.S3DB.N = 'CB Empty@SWCNT DGU B S3';
DATA_20240214.S4DB.N = 'CB Empty@SWCNT DGU B S4D';
DATA_20240214.S4LB.N = 'CB Empty@SWCNT DGU B S4L';
DATA_20240307.S5DGUB.N = 'CB Empty@SWCNT DGU B S5';
DATA_20240307.S6DGUB.N = 'CB Empty@SWCNT DGU B S6';
DATA_20240307.S7DGUB.N = 'CB Empty@SWCNT DGU B S7';
DATA_References.empty_new_P2_0601.N = 'SC Empty@SWCNT A';
DATA_References.empty_P2_0329.N = 'SC Empty@SWCNT B';
DATA_References.empty_P2_0601.N = 'SC Empty@SWCNT C';
DATA_References.empty_P2_0706.N = 'SC Empty@SWCNT D';
DATA_References.empty_P2_APA218.N = 'SC Empty@SWCNT E';
DATA_References.empty_P2_dial_0930.N = 'SC Empty@SWCNT Dial';
DATA_References.AnnoctadecaneP2.N = 'SC Annealed Octadecane@SWCNTP2';
DATA_References.AnnoctadecaneP2650.N = 'SC Annealed Octadecane@SWCNT650';
DATA_References.OctadecaneP2.N = 'SC OctadecaneP2@SWCNTP2';
DATA_References.Octadecene1P2.N = 'SC Octadecne1P2@SWCNTP2';
DATA_References.TriacontaneP2.N = 'SC Triacontane@SWCNTP2';
DATA_References.WaterFilled.N = 'SC Water@SWCNT';

ReferenceSamples = {
    DATA_20231206.SFF4dil,
    DATA_20231206.SFF4dil2,
    DATA_20231206.SFF5dil,
    %DATA_20240214.S2DB,
    %DATA_20240214.S3DB,
    %DATA_20240214.S4DB,
    %DATA_20240214.S4LB,
    %DATA_20240307.S5DGUB,
    %DATA_20240307.S6DGUB,
    %DATA_20240307.S7DGUB.
    DATA_References.AnnoctadecaneP2,
    DATA_References.AnnoctadecaneP2650,
    DATA_References.OctadecaneP2,
    DATA_References.Octadecene1P2,
    DATA_References.TriacontaneP2,
};

% Uncorrected Data for PCE Sample 2

DATA_20231206.SFF2Bdil.N = 'SF PCE@SWCNT dil1';
DATA_20231206.SFF2dil.N = 'SF PCE@SWCNT';
DATA_20240202.S2_1060.N= 'CB PCE@SWCNT 2hCF';
DATA_20240212.S2.N= 'CB PCE@SWCNT 4hCF';
DATA_20240214.S2DC.N= 'CB PCE@SWCNT DGU C (Filled)';
DATA_20240214.S2DD.N = 'CB PCE@SWCNT DGU D (Defective)';
DATA_20240214.S2DE.N = 'CB PCE@SWCNT DGU E (Bundles)';
DATA_20240216.S2DGUC.N = 'CB PCE@SWCNT Dial. DGU C (Filled)';
    
PCESamples = {
    %DATA_20231206.SFF2Bdil,
    DATA_20231206.SFF2dil,
    %DATA_20240202.S2_1060,
    %DATA_20240212.S2,
    %DATA_20240214.S2DC,
    %DATA_20240214.S2DD,
    %DATA_20240214.S2DE,
    DATA_20240216.S2DGUC,
    };


% Uncorrected Data for TCE Sample 3

DATA_20231206.SFF3_3dil.N = 'SF TCE@SWCNT RF';
DATA_20231206.SFF3dil.N = 'SF TCE@SWCNT LP Long Rinsing';
DATA_20231206.SFF3Bdil.N = 'SF TCE@SWCNT LP Short Rinsing';
DATA_20240202.S3_1060.N= 'CB TCE@SWCNT 2hCF';
DATA_20240212.S3.N= 'CB TCE@SWCNT 4hCF';
DATA_20240214.S3DC.N= 'CB TCE@SWCNT DGU C (Filled)';
DATA_20240214.S3DD.N = 'CB TCE@SWCNT DGU D (Defective)';
DATA_20240214.S3DE.N = 'CB TCE@SWCNT DGU E (Bundles)';
DATA_20240216.S3DGUC.N = 'CB TCE@SWCNT Dial. DGU C (Filled)';
    
TCESamples = {
    DATA_20231206.SFF3dil,
    DATA_20231206.SFF3_3dil,
    DATA_20231206.SFF3Bdil,
    %DATA_20240202.S3_1060,
    %DATA_20240212.S3,
    %DATA_20240214.S3DC,
    %DATA_20240214.S3DD,
    %DATA_20240214.S3DE,
    DATA_20240216.S3DGUC
    };



% Uncorrected Data for TEMED Sample 4

DATA_20231206.SFF1Bdil.N = 'SF TEMED@SWCNT';
DATA_20240202.S4_1060.N = 'CB TEMED@SWCNT 2hCF';
DATA_20240212.S4.N = 'CB TEMED@SWCNT 4hCF';
DATA_20240214.S4DC.N = 'CB TEMED@SWCNT DGU C (Filled) D';
DATA_20240214.S4LC.N = 'CB TEMED@SWCNT DGU C (Filled) L';
DATA_20240214.S4DD.N = 'CB TEMED@SWCNT DGU D (Defective)';
DATA_20240214.S4DE.N = 'CB TEMED@SWCNT DGU E (Bundles)';
DATA_20240216.S4DDGUC.N = 'CB TEMED@SWCNT Dial. DGU C (Filled)';
DATA_20240216.S4DDGUC2.N = 'CB TEMED@SWCNT Dial. DGU C (Filled)';
DATA_20240216.S4LDGUC.N = 'CB TEMED@SWCNT Dial. DGU C (Filled)';

TEMEDSamples = {

DATA_20231206.SFF1Bdil,
%DATA_20240202.S4_1060,
%DATA_20240212.S4,
%DATA_20240214.S4DC,
%DATA_20240214.S4LC,
%DATA_20240214.S4DD,
%DATA_20240214.S4DE,
DATA_20240216.S4DDGUC,
%DATA_20240216.S4DDGUC2,
%DATA_20240216.S4LDGUC
    };


% Uncorrected Data for TDAE Sample 5

DATA_20240305.S5.N = 'CB TDAE@SWCNT A';
DATA_20240305.S5_dil.N = 'CB TDAE@SWCNT B';
DATA_20240305.S5_dil2.N = 'CB TDAE@SWCNT C';
DATA_20240305.S5_dil4.N = 'CB TDAE@SWCNT';
DATA_20240307.S5CF.N = 'CB TDAE@SWCNT 4hCF';
DATA_20240307.S5DGUC.N = 'CB TDAE@SWCNT DGU C (Filled)';
DATA_20240308.S5DGUCDial.N = 'CB TDAE@SWCNT Dial. DGU C (Filled)';

TDAESamples = {
    DATA_20240305.S5_dil4,
    DATA_20240307.S5CF,
    DATA_20240307.S5DGUC,
    DATA_20240308.S5DGUCDial
    };

% Uncorrected Data for Hexadecane Sample 6

DATA_20240304.S6.N = 'CB Hexadecane@SWCNT A';
DATA_20240307.S6CF.N = 'CB Hexadecane@SWCNT 4hCF';
DATA_20240307.S6DGUC.N = 'CB Hexadecane@SWCNT DGU C (Filled)';
DATA_20240308.S6DGUCDial.N = 'CB Hexadecane@SWCNT Dial. DGU C (Filled)';

HexadecaneSamples = {
    DATA_20240304.S6,
    DATA_20240307.S6CF,
    DATA_20240307.S6DGUC,
    DATA_20240308.S6DGUCDial,
    };


% Uncorrected Data for Dodecane Sample 7

DATA_20240304.S7.N = 'CB Dodecane@SWCNT A';
DATA_20240307.S7CF.N = 'CB Dodecane@SWCNT 4hCF';
DATA_20240307.S7DGUC.N = 'CB Dodecane@SWCNT DGU C (Filled)';
DATA_20240308.S7DGUCDial.N = 'CB Dodecane@SWCNT Dial. DGU C (Filled)';

DodecaneSamples = {
    DATA_20240304.S7,
    DATA_20240307.S7CF,
    DATA_20240307.S7DGUC,
    DATA_20240308.S7DGUCDial,
    };

%plotSampleList(PCESamples)
%plotSampleList(TCESamples)
%plotSampleList(TEMEDSamples,0)
%plotSampleList(TDAESamples)
%plotSampleList(HexadecaneSamples)
%plotSampleList(DodecaneSamples)
%plotSampleList(ReferenceSamples, 0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Compare all samples after DGU

%DATA_20240308.CORRECTED = DATA_20240308.S7DGUCDial;
%DATA_20240308.CORRECTED.A = DATA_20240308.S7DGUCDial.A - DATA_References.H2OinD2O.A*0.0
%plotSampleList({DATA_20240308.S7DGUCDial, DATA_20240308.CORRECTED, DATA_References.H2OinD2O},0.0)


%H2O in D2O Correction
DATA_20240216.S2DGUC.A = DATA_20240216.S2DGUC.A - DATA_References.H2OinD2O.A*1.22;
DATA_20240216.S3DGUC.A = DATA_20240216.S3DGUC.A - DATA_References.H2OinD2O.A*0.7;
DATA_20240216.S4DDGUC.A = DATA_20240216.S4DDGUC.A - DATA_References.H2OinD2O.A*0.25;
%SAMPLE 5 Doesnt need correction
%SAMPLE 6 Doesnt need correction
%SAMPLE 7 Doesnt need correction

Compare = {
    
    %DATA_References.empty_new_P2_0601,
    %DATA_References.empty_P2_0329,
    %DATA_References.empty_P2_0601,
    %DATA_References.empty_P2_0706,
    %DATA_References.empty_P2_APA218,
    DATA_References.empty_P2_dial_0930,
    
    %Alkanes
    %DATA_References.TriacontaneP2,
    %DATA_References.AnnoctadecaneP2,
    %DATA_References.AnnoctadecaneP2650,  
    DATA_20240308.S7DGUCDial,
    DATA_20240308.S6DGUCDial,
    
    %Dopants
    DATA_20240216.S2DGUC,
    DATA_20240216.S3DGUC,
    DATA_20240216.S4DDGUC,
    DATA_20240308.S5DGUCDial,
    
    DATA_References.WaterFilled
    %DATA_20231206.SFF5dil,

    };

%plotSampleList()
%For normalization 
%plotSampleList(Compare,0)

%LL1= 580;
%UL1= 610;

LL1= 750;
UL1= 860;

LL2= 1100;
UL2= 1400;

plotSampleList(Compare,0.0)

CS = correctSpectra(Compare,LL1, UL1, LL2, UL2);
%plotSampleList(CS,0.0)

LL = 900;
UL = 1100;
NS = Normalize(CS,LL, UL);
plotSampleList(NS,0.0)

NS = Normalize(Compare,LL, UL);
plotSampleList(NS,0.0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dataStructures = ReadFromPaths(paths)
    dataStructures = struct();
    for i = 1:length(paths)
        % Extract the suffix from the path
        try
            % Extract the suffix from the path variable
            dataset_name = strsplit(paths{i}, "\");
            % Create dynamic field name for the structure
            fieldName = ['DATA_', strrep(dataset_name{7}, '.csv', '')];
            % Read the data from the current path
            data = readSamplesData(paths{i});
            assignin('caller', fieldName, data); % Assign data to a variable in the caller workspace

            % Create a structure with dynamic field name and store the data
            structure = struct(fieldName, data);

            % Store the structure in the cell array
            dataStructures.(fieldName) = structure.(fieldName);        
        catch ME
            % Print the path that caused the error
            disp(['Error reading data from path: ' paths{i}]);
            % Re-throw the error
            rethrow(ME);
        end
    end
end

function samples = readSamplesData(filePath)
    % Read the header
    header = readcell(filePath, 'Range', 'A1:AZ1');
    header = cellfun(@(x) strrep(x, ' 100%T', ''), header, 'UniformOutput', false);
    sampleNames = header(1, 1:2:end);
    warning('off', 'MATLAB:strrep:InvalidInputType');
    % Read ONLY the datalines
    data = readmatrix(filePath, 'Range', ['A' num2str(3) ':AZ' num2str(2328)]);
    
    samples = struct();

    % Iterate through each sample and store wavelength and absorption data
    for i = 1:length(sampleNames)
        % Extract wavelength and absorption data for the current sample
        wavelengths = data(:, 2*i - 1); % Odd columns contain wavelength
        absorption = data(:, 2*i);       % Even columns contain absorption
        warning('off', 'MATLAB:strrep:InvalidInputType');
        % Store the data in the container object with the sample name
        sampleName = strrep(sampleNames{i}, '-', '_');
        samples.(sampleName).W = wavelengths;
        samples.(sampleName).A = absorption;
        samples.(sampleName).N = sampleName;
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

function maximumValue = computeMaximum(sample, lowerLimit, upperLimit)
    % Check if the sampleName is in DATA
    % Extract X and Y values for the specified sample
    x = sample.W;
    y = sample.A;
    % Define the function to integrate
    
    indicesInRange = find(x >= lowerLimit & x <= upperLimit);
    yInRange = y(indicesInRange);

    maximumValue = max(yInRange);
end

function NormedSamples = Normalize(samplesToNormalize, lowerLimit, upperLimit)
    NormedSamples = cell(size(samplesToNormalize));
    % Iterate over each sample to be normalized
    for sampleIdx = 1:length(samplesToNormalize)
        currentSample = samplesToNormalize{sampleIdx};
        %currentSample.A = currentSample.A/computeIntegral(currentSample,lowerLimit, upperLimit);
        currentSample.A = currentSample.A/computeMaximum(currentSample,lowerLimit, upperLimit);
        NormedSamples{sampleIdx} = currentSample;
    end
end

function CorrectedSpectra = correctSpectra(samplesToCorrect, LL1, UL1, LL2, UL2)
    CorrectedSpectra = cell(size(samplesToCorrect));
    for i = 1:length(samplesToCorrect)
        currentSample = samplesToCorrect{i};
        % Extract wavelength and absorption data
        W = currentSample.W;
        A = currentSample.A;
        % Find the indices corresponding to the specified wavelength ranges
        idx_range1 = find(W >= LL1 & W <= UL1);
        idx_range2 = find(W >= LL2 & W <= UL2);
        % Find the wavelength where the minimum absorption occurs in each range
        [~, min_idx_range1] = min(A(idx_range1));
        [~, min_idx_range2] = min(A(idx_range2));
        % Get the corresponding wavelengths
        lambda_min_range1 = W(idx_range1(min_idx_range1));
        lambda_min_range2 = W(idx_range2(min_idx_range2));     
        % Extract data around the two minima
        W_range1 = W(idx_range1);
        A_range1 = A(idx_range1);
        W_range2 = W(idx_range2);
        A_range2 = A(idx_range2);

        % Fit a straight line through the two minima
        fitfunc = @(p, x) p(1) ./ x;
        initialGuess = [100]; % Initial guess for the fitting parameter A
        [fitparams, ~] = lsqcurvefit(fitfunc, initialGuess, [W_range1; W_range2], [A_range1; A_range2]);
        background = fitfunc(fitparams, W);

        % Subtract the background
        currentSample.A = currentSample.A - background;
        CorrectedSpectra{i} = currentSample;
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
        % Subtract the background from the absorption data within the specified range
        A_flattened = DataStructure.(sample).A - optimalA ./ DataStructure.(sample).W;
        
        % Update the absorption data in the data structure
        
        DataStructure.(sample).A = A_flattened;
    end
    
    % Return the modified data structure
    flattenedData = DataStructure;
end

function plotSampleList(SamplesToPlot, offset)
    % Create a figure for the plot
    figure;
    % Iterate over each sample
    for sampleIdx = 1:length(SamplesToPlot)
        currentSample = SamplesToPlot{sampleIdx};
            % Get the current sample, X values, and Y values
            currentX = currentSample.W;
            currentY = currentSample.A - offset*sampleIdx;
            currentN = currentSample.N;
            plot(currentX, currentY, 'DisplayName', currentN,'LineWidth', 1.3);
            hold on; % Add spectra to the same plot
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
