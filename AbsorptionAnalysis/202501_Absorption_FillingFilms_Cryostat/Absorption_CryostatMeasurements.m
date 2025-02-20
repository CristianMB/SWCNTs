clc;
clear all;
addpath('X:\SWCNTs');
import UsefulFunctions.*;
rootpath = 'X:\Measurements Data\Absorption\';

%All paths as default
set(0,'DefaultFigureWindowStyle','normal')


CryoFS4 = [rootpath,'20241216\TDAE_VacuumMeasurements.csv'];
CryoFS7 = [rootpath,'20241216\Dodecane_VacuumMeasurements.csv'];

%Select the paths of interest
paths = {   CryoFS4
            CryoFS7
        };

%Read and structure data from the paths

ReadAbsorptionFromPaths(paths);


%% %--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DATA_20241216.Baseline.N='Baseline';
DATA_20241216.Dodecane_KIT_Air_A.N='Dodecane_KIT_Air_A';
DATA_20241216.Dodecane_KIT_Vac_1.N='Dodecane_KIT_Vac_1';
DATA_20241216.Dodecane_KIT_Vac_20.N='Dodecane_KIT_Vac_20';
DATA_20241216.Dodecane_KIT_Vac_75.N='Dodecane_KIT_Vac_75';
DATA_20241216.Dodecane_KIT_Vac_140.N='Dodecane_KIT_Vac_140';
DATA_20241216.Dodecane_KIT_Vac_350.N='Dodecane_KIT_Vac_350';
DATA_20241216.Dodecane_KIT_Vac_1175.N='Dodecane_KIT_Vac_1175';
DATA_20241216.Dodecane_KIT_N2_B.N='Dodecane_KIT_N2_B';
DATA_20241216.Dodecane_KIT_Vac_1200.N='Dodecane_KIT_Vac_1200';
DATA_20241216.Dodecane_KIT_N2_C.N='Dodecane_KIT_N2_C';
DATA_20241216.Dodecane_KIT_Vac_1300.N='Dodecane_KIT_Vac_1300';
DATA_20241216.Dodecane_KIT_Temp_50.N='Dodecane_KIT_Temp_50';
DATA_20241216.Dodecane_KIT_Temp_50_B.N='Dodecane_KIT_Temp_50_B';
DATA_20241216.Dodecane_KIT_Vac_1420.N='Dodecane_KIT_Vac_1420';
DATA_20241216.Dodecane_KIT_HotWater_A.N='Dodecane_KIT_HotWater_A';
DATA_20241216.Dodecane_KIT_N2_D.N='Dodecane_KIT_N2_D';
DATA_20241216.Dodecane_KIT_Temp_75.N='Dodecane_KIT_Temp_75';
DATA_20241216.Dodecane_KIT_Temp_85.N='Dodecane_KIT_Temp_85';
DATA_20241216.Dodecane_KIT_Temp_98.N='Dodecane_KIT_Temp_98';
DATA_20241216.Dodecane_KIT_Air_B.N='Dodecane_KIT_Air_B';

DATA_20241216.Baseline.N='Baseline';
DATA_20241216.TDAE_KIT_Air.N='TDAE_KIT_Air';
DATA_20241216.TDAE_KIT_Vac_1.N='TDAE_KIT_Vac_1';
DATA_20241216.TDAE_KIT_Vac_10.N='TDAE_KIT_Vac_10';
DATA_20241216.TDAE_KIT_Vac_20.N='TDAE_KIT_Vac_20';
DATA_20241216.TDAE_KIT_Vac_30.N='TDAE_KIT_Vac_30';
DATA_20241216.TDAE_KIT_Vac_60.N='TDAE_KIT_Vac_60';
DATA_20241216.TDAE_KIT_N2.N='TDAE_KIT_N2';
DATA_20241216.TDAE_KIT_Vac_150.N='TDAE_KIT_Vac_150';
DATA_20241216.TDAE_KIT_Vac_300.N='TDAE_KIT_Vac_300';
DATA_20241216.TDAE_KIT_Vac_1254.N='TDAE_KIT_Vac_1254';
DATA_20241216.TDAE_KIT_N2_B.N='TDAE_KIT_N2_B';
DATA_20241216.TDAE_KIT_Vac_1282.N='TDAE_KIT_Vac_1282';
DATA_20241216.TDAE_KIT_Vac_1293.N='TDAE_KIT_Vac_1293';
DATA_20241216.TDAE_KIT_Vac_1300.N='TDAE_KIT_Vac_1300';
DATA_20241216.TDAE_KIT_Vac_1300.N='TDAE_KIT_Vac_1300';
DATA_20241216.TDAE_KIT_Vac_1425.N='TDAE_KIT_Vac_1425';
DATA_20241216.TDAE_KIT_Air_B.N='TDAE_KIT_Air_B';

%% %--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% %--------BACKGROUND CORRECTION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %--------SPECTRA SELECTION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Rinsings


FS4 = {
%         DATA_20241216.Baseline
        DATA_20241216.TDAE_KIT_Air
%         DATA_20241216.TDAE_KIT_Vac_1
%         DATA_20241216.TDAE_KIT_Vac_10
%         DATA_20241216.TDAE_KIT_Vac_20
%         DATA_20241216.TDAE_KIT_Vac_30
%         DATA_20241216.TDAE_KIT_Vac_60
%         DATA_20241216.TDAE_KIT_N2
        DATA_20241216.TDAE_KIT_Vac_150
%         DATA_20241216.TDAE_KIT_Vac_300
%         DATA_20241216.TDAE_KIT_Vac_1254
%         DATA_20241216.TDAE_KIT_N2_B
%         DATA_20241216.TDAE_KIT_Vac_1282
        DATA_20241216.TDAE_KIT_Vac_1293
        DATA_20241216.TDAE_KIT_Vac_1300
%         DATA_20241216.TDAE_KIT_Vac_1300
        DATA_20241216.TDAE_KIT_Vac_1425
        DATA_20241216.TDAE_KIT_Air_B
    };

FS7 = {
%         DATA_20241216.Baseline
        DATA_20241216.Dodecane_KIT_Air_A
%         DATA_20241216.Dodecane_KIT_Vac_1
%         DATA_20241216.Dodecane_KIT_Vac_20
%         DATA_20241216.Dodecane_KIT_Vac_75
        DATA_20241216.Dodecane_KIT_Vac_140
%         DATA_20241216.Dodecane_KIT_Vac_350
%         DATA_20241216.Dodecane_KIT_Vac_1175
%         DATA_20241216.Dodecane_KIT_N2_B
        DATA_20241216.Dodecane_KIT_Vac_1200
%         DATA_20241216.Dodecane_KIT_N2_C
%         DATA_20241216.Dodecane_KIT_Vac_1300
%         DATA_20241216.Dodecane_KIT_Temp_50
        DATA_20241216.Dodecane_KIT_Temp_50_B
%         DATA_20241216.Dodecane_KIT_Vac_1420
%         DATA_20241216.Dodecane_KIT_HotWater_A
%         DATA_20241216.Dodecane_KIT_N2_D
%         DATA_20241216.Dodecane_KIT_Temp_75
%         DATA_20241216.Dodecane_KIT_Temp_85
        DATA_20241216.Dodecane_KIT_Temp_98
        DATA_20241216.Dodecane_KIT_Air_B
    };


% FS4 = FilterDataByXRange(FS4, 0, 2600);
% FS4 = NormalizeSample(FS4,902, 1300); 
% FS4 = RemovePolyBG(FS4, 0);
% FS4 = NormalizeSample(FS4,902, 1300); 
% plotAbsorptionOrdered(FS4, 0);

FS7 = FilterDataByXRange(FS7, 0, 2600);
FS7 = NormalizeSample(FS7,902, 1300); 
FS7 = RemovePolyBG(FS7, 0);
FS7 = NormalizeSample(FS7,902, 1300); 
plotAbsorptionOrdered(FS7, 0);

% close;
% backgr = [330, 610, 840, 1320, 2500];
backgr = [330,610, 1311, 2500];

% FS4 = BackgroundSubtractionExcludeRanges(FS4, [[620,800], [900, 1220], [1578, 2260]]);
% FS4 = BackgroundSubtractionWithSpecifiedPoints(FS4, backgr);
% FS4 = Normalize(FS4,902, 1300, 'I'); 
% plotAbsorptionOrdered(FS4, 0);

% plotAbsorptionOrdered(FS4, 0);
FS4 = FilterDataByXRange(FS4, 0, 2520);
FS4 = BackgroundSubtraction(FS4, [500, 2600]);
FS4 = Normalize(FS4,910, 1290, 'M'); 
% plotAbsorptionOrdered(FS4, 0);

% Define los límites del rango
x_min = 800;
x_max = 860;

% Itera sobre cada dataset en la lista
for i = 1:length(FS4)
    DS = FS4{i};
    
    % Encuentra los índices donde DS.X está entre 800 y 860
    idx_range = (DS.X >= x_min) & (DS.X <= x_max);
    
    % Encuentra el valor mínimo de DS.Y en ese rango
    min_val = min(DS.Y(idx_range));
    
    
    idx_range2 = (DS.X >= 0) & (DS.X <= 900);

    % Sustrae este valor mínimo de DS.Y solo en el rango especificado
    DS.Y(idx_range2) = DS.Y(idx_range2) - min_val;
    
    % Guarda los cambios en el dataset actual
    FS4{i} = DS;
end


plotAbsorptionOrdered(FS4, 0);

% Define the list of datasets

% FS7 = FilterDataByXRange(FS7, 0, 2520);
% FS7 = BackgroundSubtraction(FS7, [500, 2500]);
% FS7 = Normalize(FS7,902, 1300, 'M'); 
% plotAbsorptionOrdered(FS7, 0);


%RIGHT PARAMETERS!! - Dont modify
% plotAbsorptionOrdered(FS7, 0);

% FS7 = FilterDataByXRange(FS7, 0, 2520);
% FS7 = BackgroundSubtraction(FS7, [250, 2500]);
% FS7 = Normalize(FS7,902, 1300, 'M'); 
% plotAbsorptionOrdered(FS7, 0);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DSListOut = BackgroundSubtractionExcludeRanges(DSList, excludeRanges)
    % BackgroundSubtractionExcludeRanges performs background subtraction using the Naumov model,
    % excluding specified ranges from the background fit.
    %
    % Inputs:
    %   DSList: Array of data structures with fields X (wavelength) and Y (absorption).
    %   excludeRanges: Nx2 matrix where each row defines [rmin, rmax] to exclude.
    %
    % Output:
    %   DSListOut: Modified DSList with background-subtracted Y values.

    DSListOut = DSList; % Initialize output as the input DSList

    for i = 1:length(DSList)
        % Extract X and Y values from the current spectrum
        xx = DSList{i}.X; % Wavelength
        yy = DSList{i}.Y; % Absorption

        % Ensure X is in ascending order
        if xx(1) > xx(end)
            xx = flip(xx); 
            yy = flip(yy); 
        end

        % Interpolate the data to ensure even spacing
        xx_interp = round(xx(1)):1:round(xx(end)); % Interpolation range
        yy_interp = interp1(xx, yy, xx_interp, 'linear'); % Interpolated Y values
        xx_interp = xx_interp'; 
        yy_interp = yy_interp';

        % Identify indices to exclude based on excludeRanges
        excludeMask = false(size(xx_interp));
        for j = 1:size(excludeRanges, 1)
            range = excludeRanges(j, :);
            excludeMask = excludeMask | (xx_interp >= range(1) & xx_interp <= range(2));
        end

        % Select only points outside the excluded ranges
        bgPoints = xx_interp(~excludeMask);
        bgY = yy_interp(~excludeMask);

        % Check if sufficient points remain for fitting
        if numel(bgPoints) < 2
            error('Not enough points outside the excluded ranges for fitting.');
        end

        % Write the background data to a temporary text file
        dataToWrite = [bgPoints, bgY];
        if isempty(dataToWrite)
            error('No data to write to the temporary file.');
        else
            disp('Writing background data to temp_data.txt');
            dlmwrite('temp_data.txt', dataToWrite);
        end

        % Optimization options for fmincon
        options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'off');

        % Perform the optimization to find A and b (initial guess: [0.2, 0.00002])
        xn = fmincon(@Naumov, [0.2, 0.00002], [], [], [], [], [], [], @conf_naumov, options); 
        % Naumov optimization using the specified points

        % Define the Naumov background model
        F_naumov = @(x) x(1) * exp(-x(2) * xx_interp);

        % Evaluate the background using the Naumov function
        BKG = F_naumov(xn);

        % Subtract the background from the interpolated Y values
        correctedY = yy_interp - BKG;

        % Update the Y values in the original X range
        DSListOut{i}.Y = yy - interp1(xx_interp, BKG, xx, 'linear', 'extrap');
        DSListOut{i}.X = xx;

        % Delete the temporary file after use
        delete('temp_data.txt');

        % Plot (optional, for debugging/visualization)
%         figure;
%         hold on;
%         plot(xx, yy, 'b', 'DisplayName', 'Original Spectrum');
%         plot(xx_interp, BKG, 'r--', 'DisplayName', 'Background Fit');
%         plot(xx, DSListOut{i}.Y, 'g', 'DisplayName', 'Corrected Spectrum');
%         legend('show');
%         title(['Background Subtraction - Spectrum ', num2str(i)]);
%         xlabel('X (Wavelength)');
%         ylabel('Y (Absorption)');
%         hold off;
    end
end

function DSListOut = BackgroundSubtractionWithSpecifiedPoints(DSList, bgPoints)
    % BackgroundSubtractionWithSpecifiedPoints performs background subtraction using the Naumov model,
    % fitting the background based on user-specified X-values (bgPoints).
    %
    % Inputs:
    %   DSList: Array of data structures with fields X (wavelength) and Y (absorption).
    %   bgPoints: Array of X-values where the background is calculated.
    %
    % Output:
    %   DSListOut: Modified DSList with background-subtracted Y values.

    DSListOut = DSList; % Initialize output as the input DSList

    % Loop over each spectrum in DSList
    for i = 1:length(DSList)
        % Extract X and Y values from the current spectrum
        xx = DSList{i}.X; % Wavelength
        yy = DSList{i}.Y; % Absorption

        % Ensure X is in ascending order
        if xx(1) > xx(end)
            xx = flip(xx); 
            yy = flip(yy); 
        end

        % Interpolate the data to ensure even spacing
        xx_interp = round(xx(1)):1:round(xx(end)); % Interpolation range
        yy_interp = interp1(xx, yy, xx_interp, 'linear'); % Interpolated Y values
        xx_interp = xx_interp'; 
        yy_interp = yy_interp';

        % Ensure bgPoints is a column vector
        bgPoints = bgPoints(:);

        % Interpolate Y-values at the specified background points
        bgY = interp1(xx_interp, yy_interp, bgPoints, 'linear', 'extrap');
        bgY = bgY(:); % Ensure bgY is also a column vector

        % Write the background data to a temporary text file
        dataToWrite = [bgPoints, bgY];
        if isempty(dataToWrite)
            error('No data to write to the temporary file.');
        else
            disp('Writing background data to temp_data.txt');
            dlmwrite('temp_data.txt', dataToWrite);
        end

        % Optimization options for fmincon
        options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'off');

        % Perform the optimization to find A and b (initial guess: [0.2, 0.00002])
        xn = fmincon(@Naumov, [0.2, 0.00002], [], [], [], [], [], [], @conf_naumov, options); 
        % Naumov optimization using the specified points

        % Define the Naumov background model
        F_naumov = @(x) x(1) * exp(-x(2) * xx_interp);

        % Evaluate the background using the Naumov function
        BKG = F_naumov(xn);

        % Subtract the background from the interpolated Y values
        correctedY = yy_interp - BKG;

        % Update the Y values in the original X range
        DSListOut{i}.Y = yy - interp1(xx_interp, BKG, xx, 'linear', 'extrap');
        DSListOut{i}.X = xx;

        % Delete the temporary file after use
        delete('temp_data.txt');

%         % Plot (optional, for debugging/visualization)
%         figure;
%         hold on;
%         plot(xx, yy, 'b', 'DisplayName', 'Original Spectrum');
%         plot(xx_interp, BKG, 'r--', 'DisplayName', 'Background Fit');
%         plot(xx, DSListOut{i}.Y, 'g', 'DisplayName', 'Corrected Spectrum');
% %         scatter(bgPoints, bgY, 'k', 'DisplayName', 'Background Points', 'filled');
%         legend('show');
%         title(['Background Subtraction - Spectrum ', num2str(i)]);
%         xlabel('X (Wavelength)');
%         ylabel('Y (Absorption)');
%         hold off;
    end
end

function DSListOut = BackgroundSubtraction(DSList, range)
    % DSList is the input array of data structures with fields X and Y
    
    DSListOut = DSList; % Initialize output as the input DSList
    
%     figure
%     hold on
    
    % Loop over each structure in DSList
    for i = 1:length(DSList)
        % Extract the X and Y values from the current data structure
        xx = DSList{i}.X; % Wavelength
        yy = DSList{i}.Y; % Absorption
        
        if xx(1) > xx(end)
            xx = flip(xx); % Flip xx to ascending order
            yy = flip(yy); % Flip yy to match xx
        end

        % Interpolate the data to ensure it is evenly spaced
        xx_interp = round(xx(1)):1:round(xx(end)); % Interpolation range
        yy_interp = interp1(xx, yy, xx_interp, 'linear'); % Interpolated Y values
        xx_interp = xx_interp'; 
        yy_interp = yy_interp';
        
               
         % Find the index for range(1) (closest value greater than or equal to range(1))
         [~, start_idx] = min(abs(xx_interp - range(1))); % Closest value to range(1)
         [~, end_idx] = min(abs(xx_interp - range(2))); % Closest value to range(2)
         
       
         % Write the temporary data to a text file
        dataToWrite = [xx_interp(start_idx:end_idx), yy_interp(start_idx:end_idx)];
        if isempty(dataToWrite)
            error('No data to write to the temporary file.');
        else
            disp('Writing data to temp_data.txt');
            dlmwrite('temp_data.txt', dataToWrite);
        end
        
        % Optimization options for fmincon
        options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'off');
        
        % Perform the optimization to find A and b (initial guess: [0.2, 0.00002])
        xn = fmincon(@Naumov, [.2, 0.00002], [], [], [], [], [], [], @conf_naumov, options); % % Approach of Naumov et al. requires two starting values for A (0.2) and b (0.002).
        
        % Define the background subtraction model
        F_naumov = @(x) double(yy_interp(start_idx:end_idx) - x(1) * exp(-x(2) * xx_interp(start_idx:end_idx))); 
        
        BKG = yy_interp(start_idx:end_idx) - F_naumov(xn);
                 
        % Update the Y values in the data structure with the background-subtracted data
        DSListOut{i}.Y = yy - interp1(xx_interp(start_idx:end_idx), BKG, xx, 'linear', 'extrap');
        DSListOut{i}.X = xx;
         
        % Delete the temporary file after use
        delete('temp_data.txt');
    end
end

function err = Naumov(x) % Calculates the difference between the absorption data and the background. The MATLAB function "fmincon" tries to minimize this difference by fitting x(1)=A and x(2)=b

A=dlmread('temp_data.txt');
c = A(:,2)-x(1)*exp(-x(2).*A(:,1));
err = double(sum(c));

end

function [c,ceq] = conf_naumov(x) % Constraint function, that forces the background to be smaller than the absorption data for every single wavelength

A=dlmread('temp_data.txt');
% Nonlinear inequality constraints
c = double(x(1)*exp(-x(2).*A(:,1))-A(:,2));
% Nonlinear equality constraints
ceq = [];
end

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
        % Findpeaks
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
