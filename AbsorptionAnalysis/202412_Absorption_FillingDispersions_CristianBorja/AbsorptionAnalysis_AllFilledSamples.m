clc;
clear;
addpath('X:\SWCNTs');
import UsefulFunctions.*;
rootpath = 'X:\Measurements Data\Absorption\';

%All paths as default

References = [rootpath,'References.csv'];
SF_CF = [rootpath,'20231206\alldata20231206.csv'];
CB_CF_S234 = [rootpath,'20240212\AllDataCentrifuge.csv'];
CB_CF_S567 = [rootpath,'20240307\DGUS5S6S7.csv'];
CB_CF_SF6_S5789 = [rootpath,'20241022\CF_FilledSamples.csv'];
KIT_FILM_S234567 = [rootpath,'20241003\FilledSamplesFilms.csv'];
CB_DGU_S234 = [rootpath,'20240216\DGUC.csv'];
CB_DGU_S567 = [rootpath,'20240308\DialS5S6S7.csv'];
CB_DGU_AP = [rootpath,'20241030\DGU_APCNTs.csv'];
KIT_FILM_Round2 = [rootpath,'20241202\Films_KIT_Set2.csv'];

%Select the paths of interest
paths = {   References
            SF_CF
            CB_CF_S234
            CB_CF_S567
            CB_CF_SF6_S5789
            KIT_FILM_S234567
            CB_DGU_S234
            CB_DGU_S567
            CB_DGU_AP
            KIT_FILM_Round2
        };

%Read and structure data from the paths

ReadAbsorptionFromPaths(paths);


%% %--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
DATA_20240212.S2.N='CB PCE@SWCNT 4hCF';
DATA_20240212.S3.N='CB TCE@SWCNT 4hCF';
DATA_20240212.S4.N='CB TEMED@SWCNT 4hCF';
DATA_20240216.S2DGUC.N='CB PCE@SWCNT Dial. DGU C (Filled';
DATA_20240216.S3DGUC.N='CB TCE@SWCNT Dial. DGU C (Filled)';
DATA_20240216.S4DDGUC.N='CB TEMED@SWCNT Dial. DGU C (Filled)';
DATA_20240216.S4DDGUC2.N='CB TEMED@SWCNT Dial. DGU C (Filled)';
DATA_20240216.S4LDGUC.N='CB TEMED@SWCNT Dial. DGU C (Filled)';
DATA_20240216.Baseline.N='DOCD2O Baseline';
DATA_20240307.S7CF.N='CB Dodecane@SWCNT 4hCF';
DATA_20240307.S7DGUC.N='CB Dodecane@SWCNT DGU C (Filled)';
DATA_20240307.S5DGUB.N='CB Empty@SWCNT DGU B S5';
DATA_20240307.S6DGUB.N='CB Empty@SWCNT DGU B S6';
DATA_20240307.S7DGUB.N='CB Empty@SWCNT DGU B S7';
DATA_20240307.S6CF.N='CB Hexadecane@SWCNT 4hCF';
DATA_20240307.S6DGUC.N='CB Hexadecane@SWCNT DGU C (Filled)';
DATA_20240307.S5CF.N='CB TDAE@SWCNT 4hCF';
DATA_20240307.S5DGUC.N='CB TDAE@SWCNT DGU C (Filled)';
DATA_20240307.Baseline.N='DOCD2O Baseline';
DATA_20240308.S7DGUCDial.N='CB Dodecane@SWCNT Dial. DGU C (Filled)';
DATA_20240308.S6DGUCDial.N='CB Hexadecane@SWCNT Dial. DGU C (Filled)';
DATA_20240308.S5DGUCDial.N='CB TDAE@SWCNT Dial. DGU C (Filled)';
DATA_20241003.FS2.N='PCE@P2-ASWCNTs (Film)';
DATA_20241003.FS3.N='TCE@P2-ASWCNTs (Film)';
DATA_20241003.FS4.N='TEMED@P2-ASWCNTs (Film)';
DATA_20241003.FS5.N='TDAE@P2-ASWCNTs (Film)';
DATA_20241003.FS6.N='Hexadecane@P2-ASWCNTs (Film)';
DATA_20241003.FS7.N='Dodecane@P2-ASWCNTs (Film)';
DATA_20241003.Baseline.N='Sapphire Baseline';
DATA_20241022.Baseline.N='1%DOC/D2O';
DATA_20241022.S5CF.N='TDAE@P2-ASWCNTs (24hCF)';
DATA_20241022.S7CF.N='Dodecane@P2-ASWCNTs (24hCF)';
DATA_20241022.S8CF.N='PCE@P2-ASWCNTs (24hCF)';
DATA_20241022.S9CF.N='TEMED@P2-ASWCNTs (24hCF)';
DATA_20241022.SFF6CF.N='TTF@P2-ASWCNTs (24hCF)';
DATA_20241030.Baseline.N='1%DOC/D2O';
DATA_20241030.StockSolution.N='StockSolution';
DATA_20241030.DGU_Empty.N='Empty fraction';
DATA_20241030.DGU_Filled.N='Water filled fraction';
DATA_20241030.CF_APCNTs.N='Centrifuged AP CNTs';
DATA_20241030.DGUDial_Filled.N='Water filled fraction (after dialysis)';
DATA_20241030.DGUDial_Empty.N='Empty fraction (after dialysis)';

DATA_20241202.Baseline.N='Film Holder Baseline';
DATA_20241202.F3TTF.N='KIT Film - TTF@SWCNTs ';
DATA_20241202.F5PCE.N='KIT Film - PCE@SWCNTs ';
DATA_20241202.F1Dodecane.N='KIT Film - Dodecane@SWCNTs ';
DATA_20241202.F2TDAE.N='KIT Film - TDAE@SWCNTs ';
DATA_20241202.F4TEMED.N='KIT Film - TEMED@SWCNTs ';

%% %--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% %--------BACKGROUND CORRECTION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %--------SPECTRA SELECTION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


CF = { 
        DATA_20231206.SFF4dil
        DATA_20231206.SFF2dil
        DATA_20231206.SFF5dil
        DATA_20231206.SFF3dil
        DATA_20231206.SFF3Bdil
        DATA_20231206.SFF3_3dil
        DATA_20231206.SFF1Bdil
        DATA_20231206.SFF6dil
%         DATA_20240212.S2
%         DATA_20240212.S3
%         DATA_20240212.S4
%         DATA_20240307.S7CF
%         DATA_20240307.S6CF
%         DATA_20240307.S5CF
%         DATA_20241003.FS2
%         DATA_20241003.FS3
%         DATA_20241003.FS4
%         DATA_20241003.FS5
%         DATA_20241003.FS6
%         DATA_20241003.FS7
%         DATA_20241022.S5CF
%         DATA_20241022.S7CF
%         DATA_20241022.S8CF
%         DATA_20241022.S9CF
%         DATA_20241022.SFF6CF
%         DATA_20241030.CF_APCNTs
        };
        

CF = NormalizeSample(CF,800, 1100);
% plotAbsorption(CF, 0.0)


DATA_References.empty_P2_dial_0930.N = 'Empty SWCNTs (DGU)'
DATA_References.WaterFilled.N = 'Water@SWCNT (DGU)'
DATA_20240216.S2DGUC.N = 'PCE@SWCNT (DGU)'
DATA_20240216.S3DGUC.N = 'TCE@SWCNT (DGU)'
DATA_20240216.S4DDGUC.N = 'TEMED@SWCNT (DGU)'
DATA_20240308.S5DGUCDial.N = 'TDAE@SWCNT (DGU)'
DATA_20240308.S6DGUCDial.N = 'Hexadecane@SWCNT (DGU)'
DATA_20240308.S7DGUCDial.N = 'Dodecane@SWCNT (DGU)'

DGU = { 
%         DATA_References.empty_P2_dial_0930

        DATA_20240216.S2DGUC
        DATA_20240216.S3DGUC
        DATA_20240216.S4DDGUC
        DATA_20240308.S5DGUCDial        
        DATA_20240308.S6DGUCDial
        DATA_20240308.S7DGUCDial
%         DATA_20241030.DGUDial_Filled
%         DATA_References.WaterFilled
%         DATA_20241030.DGUDial_Empty
        };
        

DGU = NormalizeSample(DGU,815, 1250);
% plotAbsorption(DGU, 0.0)

S11_valsDGU = {};
S22_valsDGU = {};

% for i=1:length(DGU)
%     samp_i = DGU{i}
%     S11_valsDGU{i} = ComputeMaximum(samp_i, 1500, 1865)
%     S22_valsDGU{i} = ComputeMaximum(samp_i, 815, 1250)
% end


DATA_20241003.FS2.N='PCE@SWCNT (Film from S2)';
DATA_20241003.FS3.N='TCE@SWCNT (Film from S3)';
DATA_20241003.FS4.N='TEMED@SWCNT (Film from S4)';
DATA_20241003.FS5.N='TDAE@SWCNT (Film from S5)';
DATA_20241003.FS6.N='Hexadecane@SWCNT (Film from S6)';
DATA_20241003.FS7.N='Dodecane@SWCNT (Film from S7)';

Films = { 
    DATA_20241003.FS2
    DATA_20241003.FS3
    DATA_20241003.FS4
    DATA_20241003.FS5
    DATA_20241003.FS6
    DATA_20241003.FS7
        };
% plotAbsorption(Films, 0)

Films = FilterDataByXRange(Films, 225, 2600)
Films = NormalizeSample(Films,902, 1300); %Normalization to S22 peak
Films = RemovePolyBG(Films, 0)
Films = RemoveBackgroundProfile(Films, [820, 1290, 2500])
Films = NormalizeSample(Films,902, 1300); %Normalization to S22 peak
Films = RemovePolyBG(Films, 0)
Films = NormalizeSample(Films,902, 1300); %Normalization to S22 peak

plotAbsorption(Films, 0)
 
% % Films = NormalizeSample(Films,222, 2500);
% Films = NormalizeSample(Films,902, 1300); %Normalization to S22 peak
% 
% S11_vals = {};
% S22_vals = {};
% 
% % for i=1:length(Films)
% %     samp_i = Films{i}
% %     S11_vals{i} = ComputeIntegral(samp_i, 1600, 2039)
% %     S22_vals{i} = ComputeIntegral(samp_i, 902, 1300)
% % end
% % Films = SubtractInverseBG(Films, [618, 846])

%% KIT Films - Filled Samples Round 2 (2024 12 02)

DATA_20241202.F3TTF.N='TTF@SWCNT (Film from SFF6)';
DATA_20241202.F5PCE.N='PCE@SWCNT (Film from S8)';
DATA_20241202.F1Dodecane.N='Dodecane@SWCNT (Film from S7)';
DATA_20241202.F2TDAE.N='TDAE@SWCNT (Film from S5)';
DATA_20241202.F4TEMED.N='TEMED@SWCNT (Film from S9)';

KITFilms_Set2 = { 
        DATA_20241202.F5PCE
        DATA_20241202.F4TEMED
        DATA_20241202.F2TDAE
        DATA_20241202.F1Dodecane
        
        DATA_20241202.F3TTF
        
        
                };    
            
% plotAbsorption(KITFilms_Set2, 0)

KITFilms_Set2 = FilterDataByXRange(KITFilms_Set2, 225, 2600)
KITFilms_Set2 = NormalizeSample(KITFilms_Set2,902, 1300); %Normalization to S22 peak
KITFilms_Set2 = RemovePolyBG(KITFilms_Set2, 0)
KITFilms_Set2 = RemoveBackgroundProfile(KITFilms_Set2, [820, 1340, 2500])
KITFilms_Set2 = NormalizeSample(KITFilms_Set2,902, 1300); %Normalization to S22 peak
KITFilms_Set2 = RemovePolyBG(KITFilms_Set2, 0)
KITFilms_Set2 = NormalizeSample(KITFilms_Set2,902, 1300); %Normalization to S22 peak

plotAbsorption(KITFilms_Set2, 0)

%%Developing inverse substraction


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

function SpectraList = RemoveBackgroundProfile(SpectraList, Xpoints)
    % Inicializar matrices para almacenar X y Y del fondo
    for i = 1:length(SpectraList)
        sample = SpectraList{i};
        
        % Interpolar los valores Y en los puntos Xpoints que corresponden al fondo
        Ybackground = interp1(sample.X, sample.Y, Xpoints, 'linear', 'extrap');
        Xbackground = Xpoints(:);  % Asegurar que sea un vector columna
        Ybackground = Ybackground(:);

        % Definir la función objetivo para el ajuste de mínimos cuadrados
        % La función modelo es A + B/X, y queremos minimizar la suma de cuadrados
        % con penalización en los valores negativos de la corrección
       objectiveFunc = @(params) sum((Ybackground - (params(1) + params(2)./Xbackground)).^2) + ...
                                10*sum(min(0, Ybackground - (params(1) + params(2)./Xbackground)).^2) + ...
                                10*sum(min(0, Ybackground - (params(1) + params(2)./Xbackground) - (params(1) + params(2)./Xbackground)).^2);

        % Inicializar los parámetros A y B
        initialParams = [0, 0];

        % Usar un optimizador para encontrar los parámetros A y B
        options = optimset('Display', 'off');
        params = fminsearch(objectiveFunc, initialParams, options);

        % Extraer A y B para este espectro
        A = params(1);
        B = params(2);

        % Calcular el fondo ajustado para este espectro
        background = A + B ./ sample.X;

        % Restar el fondo de los valores Y para corregir el espectro
        correctedY = sample.Y - background;

        % Actualizar el espectro corregido en la lista
        sample.Y = correctedY;
        SpectraList{i} = sample;  % Actualizar la lista de espectros
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
        disp(['Structure ', num2str(i), ' Polynomial Coefficients: ', num2str(p)]);
        
        % Save the updated structure back to the list
        DSList{i} = DS;
    end
end

