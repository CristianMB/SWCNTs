clc;
clear;

rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements\Absorption\';
%rootpath = 'C:\Users\cborja\OneDrive - Universiteit Antwerpen\Measurements\Absorption\';

%All paths as default
path_baselines = [rootpath,'References.csv'];
path_20231206 = [rootpath,'20231206\alldata20231206.csv'];
path_20240202 = [rootpath,'20240202\AllData.csv'];
path_20240212a = [rootpath,'20240212\AllDataCentrifuge.csv'];
path_20240214 = [rootpath,'20240214\DGU.csv'];
path_20240216 = [rootpath,'20240216\DGUC.csv'];
path_20240304 = [rootpath,'20240304\S6S7.csv'];
path_20240305 = [rootpath,'20240305\S5.csv'];
path_20240307 = [rootpath,'20240307\DGUS5S6S7.csv'];
path_20240308 = [rootpath,'20240308\DialS5S6S7.csv'];

%Select the paths of interest
paths = {
        path_baselines,
        path_20231206,
        path_20240202,
        path_20240212a,
        path_20240214,
        path_20240216,
        path_20240304,
        path_20240305,
        path_20240307,
        path_20240308
        };

%Read and structure data from the paths
ReadFromPaths(paths);


%%%--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Data Reference
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

DATA_20240202.BASELINE.N='Baseline (DOC/D2O)';
DATA_20240214.BASELINE.N='Baseline 1%DOCD2O';
DATA_20240304.S7.N='CB Dodecane@SWCNT';
DATA_20240307.S7CF.N='CB Dodecane@SWCNT 4hCF';
DATA_20240307.S7DGUC.N='CB Dodecane@SWCNT DGU C (Filled)';
DATA_20240308.S7DGUCDial.N='CB Dodecane@SWCNT Dial. DGU C (Filled)';
DATA_20240307.S5DGUB.N='CB Empty@SWCNT DGU B S5';
DATA_20240307.S6DGUB.N='CB Empty@SWCNT DGU B S6';
DATA_20240307.S7DGUB.N='CB Empty@SWCNT DGU B S7';
DATA_20240304.S6.N='CB Hexadecane@SWCNT';
DATA_20240307.S6CF.N='CB Hexadecane@SWCNT 4hCF';
DATA_20240307.S6DGUC.N='CB Hexadecane@SWCNT DGU C (Filled)';
DATA_20240308.S6DGUCDial.N='CB Hexadecane@SWCNT Dial. DGU C (Filled)';
DATA_20240124.BASELINE.N='CB Methanol Baseline 240115';
DATA_20240129.BASELINE.N='CB Methanol Baseline 240115';
DATA_20240202.S2_1060.N='CB PCE@SWCNT 2hCF';
DATA_20240212.S2.N='CB PCE@SWCNT 4hCF';
DATA_20240214.S2DB.N='CB PCE@SWCNT DGU C (Empty)';
DATA_20240214.S2DC.N='CB PCE@SWCNT DGU C (Filled)';
DATA_20240214.S2DA.N='CB PCE@SWCNT DGU C (Nicondenz)';
DATA_20240214.S2DD.N='CB PCE@SWCNT DGU D (Defective)';
DATA_20240214.S2DE.N='CB PCE@SWCNT DGU E (Bundles)';
DATA_20240216.S2DGUC.N='CB PCE@SWCNT Dial. DGU C (Filled';
DATA_20240202.S3_1060.N='CB TCE@SWCNT 2hCF';
DATA_20240212.S3.N='CB TCE@SWCNT 4hCF';
DATA_20240214.S3DB.N='CB TCE@SWCNT DGU C (Empty)';
DATA_20240214.S3DC.N='CB TCE@SWCNT DGU C (Filled)';
DATA_20240214.S3DA.N='CB TCE@SWCNT DGU C (Nicondenz)';
DATA_20240214.S3DD.N='CB TCE@SWCNT DGU D (Defective)';
DATA_20240214.S3DE.N='CB TCE@SWCNT DGU E (Bundles)';
DATA_20240216.S3DGUC.N='CB TCE@SWCNT Dial. DGU C (Filled)';
DATA_20240305.S5.N='CB TDAE@SWCNT';
DATA_20240305.S5_DIL.N='CB TDAE@SWCNT';
DATA_20240305.S5_DIL2.N='CB TDAE@SWCNT';
DATA_20240305.S5_DIL4.N='CB TDAE@SWCNT';
DATA_20240307.S5CF.N='CB TDAE@SWCNT 4hCF';
DATA_20240307.S5DGUC.N='CB TDAE@SWCNT DGU C (Filled)';
DATA_20240308.S5DGUCDial.N='CB TDAE@SWCNT Dial. DGU C (Filled)';
DATA_20240202.S4_1060.N='CB TEMED@SWCNT 2hCF';
DATA_20240212.S4.N='CB TEMED@SWCNT 4hCF';
DATA_20240214.S4DB.N='CB TEMED@SWCNT DGU C (Empty)';
DATA_20240214.S4LB.N='CB TEMED@SWCNT DGU C (Empty)';
DATA_20240214.S4DC.N='CB TEMED@SWCNT DGU C (Filled)';
DATA_20240214.S4LC.N='CB TEMED@SWCNT DGU C (Filled)';
DATA_20240214.S4DA.N='CB TEMED@SWCNT DGU C (Nicondenz)';
DATA_20240214.S4LA.N='CB TEMED@SWCNT DGU C (Nicondenz)';
DATA_20240214.S4DD.N='CB TEMED@SWCNT DGU D (Defective)';
DATA_20240214.S4LD.N='CB TEMED@SWCNT DGU D (Defective)';
DATA_20240214.S4DE.N='CB TEMED@SWCNT DGU E (Bundles)';
DATA_20240214.S4LE.N='CB TEMED@SWCNT DGU E (Bundles)';
DATA_20240216.S4DDGUC.N='CB TEMED@SWCNT Dial. DGU C (Filled)';
DATA_20240216.S4DDGUC2.N='CB TEMED@SWCNT Dial. DGU C (Filled)';
DATA_20240216.S4LDGUC.N='CB TEMED@SWCNT Dial. DGU C (Filled)';
DATA_20240216.BASELINE.N='DOCD2O Baseline';
DATA_20240304.BASELINE.N='DOCD2O Baseline';
DATA_20240305.BASELINE.N='DOCD2O Baseline';
DATA_20240307.BASELINE.N='DOCD2O Baseline';
DATA_20240308.BASELINE.N='DOCD2O Baseline';
DATA_20240220.BASELINE.N='EthylAcetate Baseline';
DATA_20240308.S5D1.N='Nicodenz after D(1)';
DATA_20240308.S5D4.N='Nicodenz after D(4)';
DATA_20240308.S5D7.N='Nicodenz after D(7)';
DATA_20240308.S5D9.N='Nicodenz after D(9)';

DATA_20231206.SFF4dil.N='SF Methanol@SWCNT A';
DATA_20231206.SFF4dil2.N='SF Methanol@SWCNT B';
DATA_20231206.SFF2dil.N='SF PCE@SWCNT';
DATA_20231206.SFF2Bdil.N='SF PCE@SWCNT dil1';
DATA_20231206.SFF5dil.N='SF SFF5 Water@SWCNT';
DATA_20231206.SFF3dil.N='SF TCE@SWCNT LP Long Rinsing';
DATA_20231206.SFF3Bdil.N='SF TCE@SWCNT LP Short Rinsing';
DATA_20231206.SFF3_3dil.N='SF TCE@SWCNT RF';
DATA_20231206.SFF1Bdil.N='SF TEMED@SWCNT';
DATA_20231206.SFF6dil.N='SF TTF@SWCNT';

%%%--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Compare before and after H2O in D2O Correction
%SampleToCorrect = DATA_20240307.S7DGUC
%H2OD2OFactor = 0.7
%NicodenzFactor = 1.2

%Corrected = SampleToCorrect;
%Corrected.N = 'Corrected';
%Corrected.A = SampleToCorrect.A - DATA_References.H2OinD2O.A*H2OD2OFactor - DATA_References.Nicodenz.A*NicodenzFactor
%plotSampleList({Corrected, SampleToCorrect},0.0)
%plotSampleList({DATA_References.H2OinD2O, DATA_References.Nicodenz},0.0)

%-----------------------------------------------

%H2O in D2O Correction fixed corrections.
%SAMPLE 2
DATA_20240216.S2DGUC.A = DATA_20240216.S2DGUC.A - DATA_References.H2OinD2O.A*1.22;
%SAMPLE 3
DATA_20240216.S3DGUC.A = DATA_20240216.S3DGUC.A - DATA_References.H2OinD2O.A*0.7;
%SAMPLE 4
DATA_20240216.S4DDGUC.A = DATA_20240216.S4DDGUC.A - DATA_References.H2OinD2O.A*0.25;
%SAMPLE 5 Doesnt need correction
%SAMPLE 6 Doesnt need correction
%SAMPLE 7 Doesnt need correction

%%%--------SAMPLE COMPARISION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AllData = {
    %DATA_References.empty_new_P2_0601,
    %DATA_References.empty_P2_0329,
    %DATA_References.empty_P2_0601,
    %DATA_References.empty_P2_0706,
    %DATA_References.empty_P2_APA218,
    %DATA_References.empty_P2_dial_0930,
    
    %Alkanes
    %DATA_References.TriacontaneP2,
    %DATA_References.AnnoctadecaneP2,
    %DATA_References.AnnoctadecaneP2650,  
    %DATA_20240308.S7DGUCDial,
    %DATA_20240308.S6DGUCDial,
    
    %Dopants
    
    %PCE
    %DATA_20231206.SFF2Bdil,
    %DATA_20231206.SFF2dil,
    %DATA_20240216.S2DGUC,

    %TCE
    %DATA_20231206.SFF3_3dil,
    %DATA_20231206.SFF3dil,
    %DATA_20231206.SFF3Bdil,
    %DATA_20240216.S3DGUC,
   
    %TEMED
    %DATA_20240216.S4DDGUC,
    %DATA_20231206.SFF1Bdil,
    
    %TDAE
    %DATA_20240308.S5DGUCDial,
    
    %DATA_References.WaterFilled
    %DATA_20231206.SFF5dil,
    };


Alkanes = {
    %Empty
    DATA_References.empty_P2_dial_0930,
    
    %Alkanes
    DATA_References.TriacontaneP2,
    DATA_References.AnnoctadecaneP2,
    DATA_References.AnnoctadecaneP2650,  
    DATA_20240308.S7DGUCDial,
    DATA_20240308.S6DGUCDial,
    
    %Water
    DATA_References.WaterFilled
    };

PCE = {
    
    %Empty
    DATA_References.empty_P2_dial_0930,
    
    DATA_20231206.SFF2dil,
    DATA_20231206.SFF2Bdil,
    DATA_20240212.S2,
    DATA_20240216.S2DGUC,
    
    %Water
    DATA_References.WaterFilled
}

TCE = {
    
    %Empty
    DATA_References.empty_P2_dial_0930,
    DATA_20231206.SFF3dil,
    DATA_20231206.SFF3Bdil,
    DATA_20231206.SFF3_3dil,
    DATA_20240202.S3_1060,
    DATA_20240212.S3,
    DATA_20240216.S3DGUC,    
    %Water
    DATA_References.WaterFilled
}

TEMED = {
    
    %Empty
    DATA_References.empty_P2_dial_0930,
    DATA_20231206.SFF1Bdil,
    DATA_20240202.S4_1060,
    DATA_20240212.S4,
    DATA_20240216.S4DDGUC,
    DATA_20240216.S4LDGUC,
    %Water
    DATA_References.WaterFilled
}

TDAE = {
    
    %Empty
    DATA_References.empty_P2_dial_0930,

    DATA_20240305.S5,
    DATA_20240305.S5_dil4,
    DATA_20240307.S5CF,
    DATA_20240308.S5DGUCDial,
    
    %Water
    DATA_References.WaterFilled
}

Hexadecane = {
    
    %Empty
    DATA_References.empty_P2_dial_0930,

    DATA_20240304.S6,
    DATA_20240307.S6CF,
    DATA_20240308.S6DGUCDial,
    
    %Water
    DATA_References.WaterFilled
}

Dodecane = {
    
    %Empty
    DATA_References.empty_P2_dial_0930,

    DATA_20240304.S7,
    DATA_20240307.S7CF,
    DATA_20240308.S7DGUCDial,

    %Water
    DATA_References.WaterFilled
}

Dopants = {
    %Empty
    DATA_References.empty_P2_dial_0930,
    
    %PDopants
    %PCE
    DATA_20231206.SFF2Bdil,
    DATA_20231206.SFF2dil,
    DATA_20240216.S2DGUC,
    %TCE
    DATA_20231206.SFF3_3dil,
    DATA_20231206.SFF3dil,
    DATA_20231206.SFF3Bdil,
    DATA_20240216.S3DGUC,
    
    %NDopants
    %TEMED
    DATA_20240216.S4DDGUC,
    DATA_20231206.SFF1Bdil,
    
    %TDAE
    DATA_20240308.S5DGUCDial,
    
    %Water
    DATA_References.WaterFilled
    };

CFSamples = {
            DATA_References.empty_P2_dial_0930,

            DATA_20240212.S2,
            DATA_20231206.SFF2Bdil,
            DATA_20231206.SFF2dil,
            
            DATA_20240212.S3,
            DATA_20231206.SFF3_3dil,
            DATA_20231206.SFF3dil,
            DATA_20231206.SFF3Bdil,
    
            DATA_20240212.S4,
            DATA_20231206.SFF1Bdil,

            DATA_References.WaterFilled

            };

        
DialSamples = {
            DATA_References.empty_P2_dial_0930,

            DATA_20240308.S6DGUCDial
            DATA_20240308.S7DGUCDial
            DATA_20240216.S2DGUC
            DATA_20240216.S3DGUC
            DATA_20240216.S4DDGUC
            DATA_20240308.S5DGUCDial
            
            DATA_References.WaterFilled

            };



%%%--------BACKGROUND CORRECTION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Correction is made based on local minima (750, 860) and (1100, 1400)
DialSamples = correctSpectra(DialSamples,750, 850, 1150, 1250);
Alkanes = correctSpectra(Alkanes,750, 850, 1150, 1250);
CFSamples = correctSpectra(CFSamples,750, 850, 1150, 1250);
PCE = correctSpectra(PCE,750, 850, 1150, 1250);
TCE = correctSpectra(TCE,750, 850, 1150, 1250);
TEMED = correctSpectra(TEMED,750, 850, 1150, 1250);
TDAE = correctSpectra(TDAE,750, 850, 1150, 1250);
Hexadecane = correctSpectra(Hexadecane,750, 850, 1150, 1250);
Dodecane = correctSpectra(Dodecane,750, 850, 1150, 1250);



%%%--------NORMALIZATION TO S22 PEAK--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Normalization using maximum value of S22 transition, can also use Integral
LS2= 900;  
US2= 1100;

DialSamples = Normalize(DialSamples,LS2, US2);
Alkanes = Normalize(Alkanes,LS2, US2);
CFSamples= Normalize(CFSamples,LS2, US2);
PCE = Normalize(PCE,LS2, US2);
TCE = Normalize(TCE,LS2, US2);
TEMED = Normalize(TEMED,LS2, US2);
TDAE = Normalize(TDAE,LS2, US2);
Hexadecane = Normalize(Hexadecane,LS2, US2);
Dodecane = Normalize(Dodecane,LS2, US2);

%%%--------PEAK CALCULATION AND EXPORT--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Find peaks in these ranges
%S11
LS1= 1600;
US1= 1850;

%Plasmon Right
LPR= 260;
UPR= 285;
%Plasmon Left
LPL= 220;
UPL= 250;

%Normed = PeakCalculation(Normed, LS1, US1, LS2, US2, LPL, UPL, LPR, UPR);
%exportPeaksToCSV(Normed, 'PP.csv')


%%%%PLOTING

%plotSampleList(Alkanes,0.0)
%plotSampleList(Dopants, 0.0)

%plotSampleList(PCE, 2.0)
%plotSampleList(TCE, 2.0)
%plotSampleList(TDAE, 0.0)
%plotSampleList(TEMED, 0.0)
%plotSampleList(Hexadecane, 0.0)
%plotSampleList(Dodecane, 0.0)

plotSampleList(DialSamples,2.0)
plotSampleList(CFSamples, 2.0)
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

function sampleList = PeakCalculation(sampleList, LS1, US1, LS2, US2, LPL, UPL, LPR, UPR)

    % Iterate over each sample to be normalized
    for sampleIdx = 1:length(sampleList)
        currentSample = sampleList{sampleIdx};
        [S11W, S11A] = computePeak(currentSample,LS1, US1);
        currentSample.S11W = S11W;
        currentSample.S11A = S11A;
        
        [S22W, S22A] = computePeak(currentSample,LS2, US2);
        currentSample.S22W = S22W;
        currentSample.S22A = S22A;
        
        [PlasWL, PlasAL] = computePeak(currentSample,LPL, UPL);
        currentSample.PlasWL = PlasWL;
        currentSample.PlasAL = PlasAL;
        
        [PlasWR, PlasAR] = computePeak(currentSample,LPR, UPR);
        currentSample.PlasWR = PlasWR;
        currentSample.PlasAR = PlasAR;
        
        sampleList{sampleIdx} = currentSample;
    end
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

function [peakPosition, peakValue] = computePeak(sample, lowerLimit, upperLimit)
    % Check if the sampleName is in DATA
    % Extract X and Y values for the specified sample
    x = sample.W;
    y = sample.A;
    % Define the function to integrate
    indicesInRange = find(x >= lowerLimit & x <= upperLimit);
    yInRange = y(indicesInRange);
    % Find the maximum value and its index within the range
    [peakValue, maxIndex] = max(yInRange);
    % Convert the index to the corresponding X value
    peakPosition = x(indicesInRange(maxIndex));
end

function exportPeaksToCSV(sampleList, fileName)
    % Define the peak fields to include in the CSV file
    peakFields = {'N', 'S11W', 'S11A', 'S22W', 'S22A', 'PlasWL', 'PlasAL', 'PlasWR', 'PlasAR'};

    % Create a cell array to store the data
    data = cell(length(sampleList), length(peakFields));

    % Fill in the data for each sample
    for sampleIdx = 1:length(sampleList)
        currentSample = sampleList{sampleIdx};
        
        % Fill in the peak information for the current sample
        for peakIdx = 1:length(peakFields)
            peakValue = currentSample.(peakFields{peakIdx});
            data{sampleIdx, peakIdx} = peakValue;
        end
    end

    % Create column names
    columnNames = peakFields;

    % Create a table from the data
    dataTable = cell2table(data, 'VariableNames', columnNames);

    % Write the table to a CSV file
    writetable(dataTable, fileName);
end



