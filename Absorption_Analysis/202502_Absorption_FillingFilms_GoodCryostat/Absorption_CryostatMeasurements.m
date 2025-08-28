clc;
clear all;
addpath('X:\SWCNTs');
import UsefulFunctions.*;
% rootpath = 'X:\Measurements Data\Absorption\';
rootpath = 'X:\Measurements (RAW)\Absorption\';
% addpath('X:\SWCNTs\SpecialMatlabFunctions\DrosteEffect-BrewerMap-3.2.5.0')

%All paths as default
set(0,'DefaultFigureWindowStyle','normal')


CryoFS6 = [rootpath,'20250220\Cryostat_Hexadecane_Film.csv'];
CryoFS4 = [rootpath,'20250221\Cryostat_TDAE_Film.csv'];
CryoFS4_2 = [rootpath,'20250305\TDAE_Cryostat.csv'];

%Select the paths of interest
paths = {   
            CryoFS6
            CryoFS4
            CryoFS4_2
        };

%Read and structure data from the paths

ReadAbsorptionFromPaths(paths);


%% %--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DATA_20250220.Hexadecane_Initial.N = 'Hexadecane in Air, 22°C'
DATA_20250220.Hexadecane_V1.N = 'Hexadecane in HV, 22°C'
DATA_20250220.Hexadecane_V2.N = 'Hexadecane in Air, 22°C'
DATA_20250220.Hexadecane_V3.N = 'Hexadecane in Air, 22°C'
DATA_20250220.Hexadecane_V4.N = 'Hexadecane in HV, 22°C'
DATA_20250220.Hexadecane_V5_T1.N = 'Hexadecane in HV, 76°C'
DATA_20250220.Hexadecane_V6_T1.N = 'Hexadecane in HV, 76°C'
DATA_20250220.Hexadecane_V6_T2.N = 'Hexadecane in HV, 127°C'
DATA_20250220.Hexadecane_V8_T2.N = 'Hexadecane in HV, 127°C'
DATA_20250220.Hexadecane_V9_T3.N = 'Hexadecane in HV, 177°C'
DATA_20250220.Hexadecane_V10_T3.N = 'Hexadecane in HV, 177°C'
DATA_20250220.Hexadecane_V11_T4.N = 'Hexadecane in HV, 226°C'
DATA_20250220.Hexadecane_V12_T4.N = 'Hexadecane in HV, 226°C'
DATA_20250220.Hexadecane_V13_C1.N = 'Hexadecane in HV, 196°C'
DATA_20250220.Hexadecane_V14_C2.N = 'Hexadecane in HV, 177°C'
DATA_20250220.Hexadecane_V15_C3.N = 'Hexadecane in HV, 137°C'
DATA_20250220.Hexadecane_V15_C4.N = 'Hexadecane in HV, 71°C'
DATA_20250220.Hexadecane_V17_C5.N = 'Hexadecane in HV, 22°C'
DATA_20250220.Hexadecane_V18_A1.N = 'Hexadecane in Air, 22°C'

DATA_20250221.TDAE_V0.N = 'TDAE in Air, 22°C'
DATA_20250221.TDAE_V1.N = 'TDAE in HV, 22°C'
DATA_20250221.TDAE_V2_T1.N = 'TDAE in HV, 77°C'
DATA_20250221.TDAE_V3_T2.N = 'TDAE in HV, 102°C'
DATA_20250221.TDAE_V4_T3.N = 'TDAE in HV, 127°C'
DATA_20250221.TDAE_V5_T4.N = 'TDAE in HV, 137°C'
DATA_20250221.TDAE_V6_T5.N = 'TDAE in HV, 147C'
DATA_20250221.TDAE_V7_T6.N = 'TDAE in HV, 167°C'
DATA_20250221.TDAE_V8_T7.N = 'TDAE in HV, 177°C'
DATA_20250221.TDAE_V9_T8.N = 'TDAE in HV, 88°C'
DATA_20250221.TDAE_V10_T9.N = 'TDAE in HV, 22°C'
DATA_20250221.TDAE_V11_A1.N= 'TDAE in Air, 22°C'
DATA_20250221.TDAE_V12_A2.N = 'TDAE in Air, 22°C'


DATA_20250221.TDAE_V0.N = 'TDAE in Air, 22°C Initial'
DATA_20250220.Hexadecane_Initial.N = 'Hexadecane in Air, 22°C Initial'
DATA_20250220.Hexadecane_V18_A1.N = 'Hexadecane in Air, 22°C Final'
DATA_20250221.TDAE_V12_A2.N = 'TDAE in Air, 22°C Final'

%% %--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %--------BACKGROUND CORRECTION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %--------SPECTRA SELECTION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% #3DBRICKS

FS4 = {
% %     DATA_20250221.Baseline
% %     DATA_20250221.Cryostat_Air
    DATA_20250221.TDAE_V0
%     DATA_20250221.TDAE_V1
%     DATA_20250221.TDAE_V2_T1
%     DATA_20250221.TDAE_V3_T2
%     DATA_20250221.TDAE_V4_T3
%     DATA_20250221.TDAE_V5_T4
%     DATA_20250221.TDAE_V6_T5
%     DATA_20250221.TDAE_V7_T6
%     DATA_20250221.TDAE_V8_T7
%     DATA_20250221.TDAE_V9_T8
%     DATA_20250221.TDAE_V10_T9
%     DATA_20250221.TDAE_V11_A1
    DATA_20250221.TDAE_V12_A2
%     DATA_20250305.TDAE_Air_11days
     };
 
FS4 = FilterDataByXRange(FS4, 230, 2500);
FS4 = matchSpectra(FS4, 900, 5);
FS4 = matchSpectra(FS4, 1199, 5);
% FS4 = Normalize(FS4,950, 1150, 'M');
FS4 = BackgroundSubtraction(FS4, [350, 2500]);
% FS4 = Normalize(FS4,910, 1290, 'M');
% FS4 = RemovePolyBG(FS4, 0);

FS4 = Normalize(FS4,240, 1316, 'M'); 

plotAbsorption(FS4, 0);

% FS6 = {
% %         DATA_20250220.Baseline
% %         DATA_20250220.Cryostat_Air
%         DATA_20250220.Hexadecane_Initial
% %         DATA_20250220.Hexadecane_V1
% %         DATA_20250220.Hexadecane_V2
% %         DATA_20250220.Hexadecane_V3
%         DATA_20250220.Hexadecane_V4
% %         DATA_20250220.Hexadecane_V5_T1
%         DATA_20250220.Hexadecane_V6_T1
% %         DATA_20250220.Hexadecane_V6_T2
%         DATA_20250220.Hexadecane_V8_T2
% %         DATA_20250220.Hexadecane_V9_T3
%         DATA_20250220.Hexadecane_V10_T3
% %         DATA_20250220.Hexadecane_V11_T4
%         DATA_20250220.Hexadecane_V12_T4
%         DATA_20250220.Hexadecane_V13_C1
%         DATA_20250220.Hexadecane_V14_C2
%         DATA_20250220.Hexadecane_V15_C3
%         DATA_20250220.Hexadecane_V15_C4
%         DATA_20250220.Hexadecane_V17_C5
%         DATA_20250220.Hexadecane_V18_A1
%     };
% 
% FS4 = {
% %     DATA_20250221.Baseline
% %     DATA_20250221.Cryostat_Air
%     DATA_20250221.TDAE_V0
% %     DATA_20250221.TDAE_V1
% %     DATA_20250221.TDAE_V2_T1
% %     DATA_20250221.TDAE_V3_T2
% %     DATA_20250221.TDAE_V4_T3
% %     DATA_20250221.TDAE_V5_T4
% %     DATA_20250221.TDAE_V6_T5
% %     DATA_20250221.TDAE_V7_T6
% %     DATA_20250221.TDAE_V8_T7
% %     DATA_20250221.TDAE_V9_T8
% %     DATA_20250221.TDAE_V10_T9
% %     DATA_20250221.TDAE_V11_A1
%     DATA_20250221.TDAE_V12_A2
%     DATA_20250305.TDAE_Air_11days
%     };
% 
% 
% % plotMaxima(FS6, 960, 1060, 0.05);
% 
% % plotAbsorptionOrdered(FS6, 0.00);
% % 
% FS4 = FilterDataByXRange(FS4, 190, 2500);
% FS4 = matchSpectra(FS4, 900, 5);
% FS4 = matchSpectra(FS4, 1199, 5);
% FS4 = Normalize(FS4,950, 1150, 'M');
% 
% FS6 = FilterDataByXRange(FS6, 190, 2500);
% FS6 = matchSpectra(FS6, 900, 5);
% FS6 = matchSpectra(FS6, 1199, 5);
% FS6 = Normalize(FS6,950, 1150, 'M');
% % 
% % plotAbsorption(FS4, 0.00);
% % plotAbsorptionOrdered(FS4, 0.00);
% 
% FS4 = Normalize(FS4,950, 1150, 'M');
% plotAbsorption(FS4, 0.00);
% 
% % FS6 = Normalize(FS6,950, 1150, 'I');
% % plotAbsorption(FS6, 0.00);
% 
% % FS4 = ConvertedEnergy(FS4);
% % plotAbsorptionOrdered(FS4, 0.1);
% % plotMaxima(FS4, 960, 1060, 0.1);

% FS6 = FilterDataByXRange(FS6, 190, 2500);
% FS6 = matchSpectra(FS6, 900, 5);
% FS6 = matchSpectra(FS6, 1199, 5);
% FS6 = RemovePolyBG(FS6, 0);
% FS6 = Normalize(FS6,950, 1150, 'M');
% % FS6 = Normalize(FS6,906, 1380, 'I');
% 
% plotAbsorptionOrdered(FS6, 0.0001);
% plotMaxima(FS6, 1600, 2120, 0.05);
% plotMaxima(FS6, 960, 1060, 0.05);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis on S11/S22 ratio

% FS4 = BackgroundSubtractionWithSpecifiedPoints(FS4, [611, 826, 1300]);
% plotAbsorption(FS4,0);
% calculateS11S22(FS4, 960, 1060, 1600, 2120, 0.1)
% calculateS11S22(FS6, 960, 1060, 1600, 2120, 0.1)


function calculateS11S22(F, x11_min, x11_max, x22_min, x22_max, offset)
    % Store values for plotting
    categories = {};  % Stores qualitative values from DS.N
    ratios = [];      % Stores computed S11/S22 ratios

    % First Figure: Ratio vs DS.N
    figure;
    hold on;

    for i = 1:length(F)
        DS = F{i}; % Use a single dataset

        % Filter data within the specified region
        idx_range11 = (DS.X >= x11_min) & (DS.X <= x11_max);
        X11_fit = DS.X(idx_range11);
        Y11_fit = DS.Y(idx_range11);
        
        idx_range22 = (DS.X >= x22_min) & (DS.X <= x22_max);
        X22_fit = DS.X(idx_range22);
        Y22_fit = DS.Y(idx_range22);
        
        % Fit a Lorentzian model to the selected data
        lorentzEqn = 'a / (1 + ((x-b)/c)^2) + d';
        
        startPoints11 = [max(Y11_fit), mean(X11_fit), std(X11_fit), min(Y11_fit)];
        startPoints22 = [max(Y22_fit), mean(X22_fit), std(X22_fit), min(Y22_fit)];

        fitResult11 = fit(X11_fit, Y11_fit, lorentzEqn, 'Start', startPoints11);
        fitResult22 = fit(X22_fit, Y22_fit, lorentzEqn, 'Start', startPoints22);

        % Extract peak values from fits
        max_x11 = fitResult11.b;
        max_val11 = fitResult11.a + fitResult11.d;
        
        max_x22 = fitResult22.b;
        max_val22 = fitResult22.a + fitResult22.d;

        % Compute S11/S22 ratio
        ratio = max_val11 / max_val22;

        % Store qualitative variable (DS.N) and computed ratio
        categories{end+1} = DS.N;
        ratios(end+1) = ratio;
    end

    % Convert categories to a categorical variable
    categories = categorical(categories);

    % Plot the ratio vs DS.N
    scatter(categories, ratios, 'ro', 'filled');
    xlabel('DS.N (Qualitative)');
    ylabel('S11/S22 Ratio');
    title('S11/S22 Ratio vs. DS.N');
    grid on;
    hold off;

    % Second Figure: Spectra with peak markers
    figure;
    hold on;

    for i = 1:length(F)
        DS = F{i}; % Use a single dataset

        % Apply offset to Y values
        DS.Y = DS.Y - (i - 1) * offset;

        % Filter data within the specified region
        idx_range11 = (DS.X >= x11_min) & (DS.X <= x11_max);
        X11_fit = DS.X(idx_range11);
        Y11_fit = DS.Y(idx_range11);
        
        idx_range22 = (DS.X >= x22_min) & (DS.X <= x22_max);
        X22_fit = DS.X(idx_range22);
        Y22_fit = DS.Y(idx_range22);

        % Fit a Lorentzian model to the selected data
        lorentzEqn = 'a / (1 + ((x-b)/c)^2) + d';
        
        startPoints11 = [max(Y11_fit), mean(X11_fit), std(X11_fit), min(Y11_fit)];
        startPoints22 = [max(Y22_fit), mean(X22_fit), std(X22_fit), min(Y22_fit)];

        fitResult11 = fit(X11_fit, Y11_fit, lorentzEqn, 'Start', startPoints11);
        fitResult22 = fit(X22_fit, Y22_fit, lorentzEqn, 'Start', startPoints22);

        % Extract peak values from fits
        max_x11 = fitResult11.b;
        max_val11 = fitResult11.a + fitResult11.d;
        
        max_x22 = fitResult22.b;
        max_val22 = fitResult22.a + fitResult22.d;

        % Plot the spectra
        plot(DS.X, DS.Y, '-', 'DisplayName', sprintf('Spectrum %d', i));

        % Mark the peaks
        plot(max_x11, max_val11, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
        plot(max_x22, max_val22, 'bo', 'MarkerSize', 8, 'LineWidth', 2);
    end

    xlabel('X');
    ylabel('Y');
    title('Spectra with Peak Markers');
%     legend show;
    grid on;
    hold off;
end


function plotMaxima(FS6, x_min, x_max, offset)
    figure; hold on;
    for i = 1:length(FS6)
        DS = FS6{i};
        
        % Apply offset to Y values
        DS.Y = DS.Y - (i - 1) * offset;
        
        % Filter data within the specified region
        idx_range = (DS.X >= x_min) & (DS.X <= x_max);
        X_fit = DS.X(idx_range);
        Y_fit = DS.Y(idx_range);
        
        % Fit a Lorentzian model to the selected data
        lorentzEqn = 'a / (1 + ((x-b)/c)^2) + d';
        startPoints = [max(Y_fit), mean(X_fit), std(X_fit), min(Y_fit)];
        fitResult = fit(X_fit, Y_fit, lorentzEqn, 'Start', startPoints);
        
        % Find the peak from the fit
        max_x = fitResult.b;
        max_val = fitResult.a / (1 + ((max_x - fitResult.b)/fitResult.c)^2) + fitResult.d;
        
        % Plot the spectrum
        plot(DS.X, DS.Y, '-');
        
        % Mark the maximum point
        plot(max_x, max_val, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    end
    hold off;
end



% FS6 = ConvertedEnergy(FS6);
% plotAbsorptionOrdered(FS6, 0);


% FS7 = NormalizeSample(FS7,902, 1300); 
% FS7 = RemovePolyBG(FS7, 0);
% FS7 = NormalizeSample(FS7,902, 1300); 



% close;
% backgr = [330, 610, 840, 1320, 2500];
% backgr = [330,610, 1311, 2500];

% FS4 = BackgroundSubtractionExcludeRanges(FS4, [[620,800], [900, 1220], [1578, 2260]]);
% FS4 = BackgroundSubtractionWithSpecifiedPoints(FS4, backgr);
% FS4 = Normalize(FS4,902, 1300, 'I'); 
% plotAbsorptionOrdered(FS4, 0);

% plotAbsorptionOrdered(FS4, 0);
% FS4 = FilterDataByXRange(FS4, 0, 2520);
% FS4 = BackgroundSubtraction(FS4, [500, 2600]);
% FS4 = Normalize(FS4,910, 1290, 'M'); 
% % plotAbsorptionOrdered(FS4, 0);
% 
% % Define los límites del rango
% x_min = 800;
% x_max = 860;

% Itera sobre cada dataset en la lista


% plotAbsorptionOrdered(FS4, 0);

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

function FS6 = matchSpectra(FS6, x0, delta)
    for i = 1:length(FS6)
        DS = FS6{i};
        
        % Definir los rangos izquierdo y derecho alrededor de x0
        idx_left = (DS.X >= (x0 - delta)) & (DS.X <= x0);
        idx_right = (DS.X >= x0) & (DS.X <= (x0 + delta));
        
        % Calcular los promedios de Y en cada región
        mean_left = mean(DS.Y(idx_left));
        mean_right = mean(DS.Y(idx_right));
        
        % Calcular el ajuste necesario
        adjustment = mean_right - mean_left;
        
        % Aplicar el ajuste a la parte izquierda o derecha (puedes cambiar esto según necesites)
%         DS.Y(DS.X < x0) = DS.Y(DS.X < x0) + adjustment;
          DS.Y(DS.X > x0) = DS.Y(DS.X > x0) - adjustment;
        
        % Guardar los cambios en el dataset actual
        FS6{i} = DS;
    end
end

function DS_list = ConvertedEnergy(DS_list)
    h = 4.135667696e-15; % Planck's constant in eV*s
    c = 299792458; % Speed of light in m/s
    
    for i = 1:length(DS_list)
        DS = DS_list{i};
        
        % Convert wavelength (nm) to energy (eV)
        DS.X = (h * c) ./ (DS.X * 1e-9);
        
        % Store the converted dataset
        DS_list{i} = DS;
    end
end

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
