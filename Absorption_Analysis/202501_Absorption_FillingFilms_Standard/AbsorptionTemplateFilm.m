clc;
clear;
addpath('X:\SWCNTs');
import UsefulFunctions.*;
rootpath = 'X:\Measurements Data\Absorption\';
addpath 'X:\SWCNTs\SpecialMatlabFunctions\DrosteEffect-BrewerMap-3.2.5.0'

%All paths as default

KIT_FILM_S234567 = [rootpath,'20241003\FilledSamplesFilms.csv'];
KIT_FILM_Round2 = [rootpath,'20241202\Films_KIT_Set2.csv'];

%Select the paths of interest
paths = {   
            KIT_FILM_S234567
            KIT_FILM_Round2
        };

%Read and structure data from the paths

ReadAbsorptionFromPaths(paths);


%% %--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


DATA_20241003.FS2.N='PCE@P2-ASWCNTs (Film)';
DATA_20241003.FS3.N='TCE@P2-ASWCNTs (Film)';
DATA_20241003.FS4.N='TEMED@P2-ASWCNTs (Film)';
DATA_20241003.FS5.N='TDAE@P2-ASWCNTs (Film)';
DATA_20241003.FS6.N='Hexadecane@P2-ASWCNTs (Film)';
DATA_20241003.FS7.N='Dodecane@P2-ASWCNTs (Film)';
DATA_20241003.Baseline.N='Sapphire Baseline';

DATA_20241202.Baseline.N='Film Holder Baseline';
DATA_20241202.F3TTF.N='KIT Film - TTF@SWCNTs ';
DATA_20241202.F5PCE.N='KIT Film - PCE@SWCNTs ';
DATA_20241202.F1Dodecane.N='KIT Film - Dodecane@SWCNTs ';
DATA_20241202.F2TDAE.N='KIT Film - TDAE@SWCNTs ';
DATA_20241202.F4TEMED.N='KIT Film - TEMED@SWCNTs ';

%% %--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% %--------BACKGROUND CORRECTION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %--------SPECTRA SELECTION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

% plotAbsorptionOrdered(KITFilms_Set2, 0)

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

