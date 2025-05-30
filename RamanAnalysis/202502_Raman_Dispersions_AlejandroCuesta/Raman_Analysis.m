clc;
clear;
addpath('X:\SWCNTs');
import UsefulFunctions.*;
rootpath = 'X:\Measurements (RAW)\Raman\';

% rootpath = 'X:\Measurements Data\Raman\';
addpath('X:\SWCNTs\SpecialMatlabFunctions\DrosteEffect-BrewerMap-3.2.5.0')
%All paths as default
path_powder = [rootpath,'20250211\'];
path_aina = [rootpath,'20250126\'];

%Select the paths of interest

paths = {
        path_powder
        path_aina
        };


ReadRamanFromPaths(paths, 2);

%%%--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



DATA_20250211.S1BH514R.N = 'SWCNTs - DGU B';
DATA_20250211.S1CH514R.N = 'SWCNTs - DGU C';
DATA_20250211.S1DH514R.N = 'SWCNTs - DGU D';
DATA_20250211.S2BH514R.N = 'BiTeI@SWCNTs - DGU B';
DATA_20250211.S2CH514R.N = 'BiTeI@SWCNTs - DGU C';
DATA_20250211.S2DH514R.N = 'BiTeI@SWCNTs - DGU D';
DATA_20250211.S1BH514G.N = 'SWCNTs - DGU B';
DATA_20250211.S1CH514G.N = 'SWCNTs - DGU C';
DATA_20250211.S1DH514G.N = 'SWCNTs - DGU D';
DATA_20250211.S2BH514G.N = 'BiTeI@SWCNTs - DGU B';
DATA_20250211.S2CH514G.N = 'BiTeI@SWCNTs - DGU C';
DATA_20250211.S2DH514G.N = 'BiTeI@SWCNTs - DGU D';

%%%--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%--------SAMPLE COMPARISION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% AINA AND ALEJANDRO

%% CRISTIAN AND ALENADRO
               
R = {
%         DATA_20250126.L496A
%         DATA_20250126.L514A
%         
%         DATA_20250126.DOCRBMA
%         DATA_20250126.DOCRBMB
%         
%         DATA_20250126.MEGPRA
%         DATA_20250126.MERBMA
%         DATA_20250126.MERBMPRA
%         DATA_20250126.MERBMSON
        
%         DATA_20250126.SEGA
%         DATA_20250126.SEGB
        DATA_20250126.SERBMA
%         DATA_20250126.SERBMB
        
%         DATA_20250126.SFGA
%         DATA_20250126.SFGB
        DATA_20250126.SFRBMA
%         DATA_20250126.SFRBMB
            };   
        
               
G = {
        DATA_20250126.L496A
        DATA_20250126.L514A
        
        DATA_20250126.DOCRBMA
        DATA_20250126.DOCRBMB
        
        DATA_20250126.MEGPRA
        DATA_20250126.MERBMA
        DATA_20250126.MERBMPRA
        DATA_20250126.MERBMSON
        
        DATA_20250126.SEGA
        DATA_20250126.SEGB
        DATA_20250126.SERBMA
        DATA_20250126.SERBMB
        
        DATA_20250126.SFGA
        DATA_20250126.SFGB
        DATA_20250126.SFRBMA
        DATA_20250126.SFRBMB
            };   



 
R = FilterDataByXRange(R, 50, 550);
G = FilterDataByXRange(G, 1400, 1800);


% R = RemovePolyBG(R, 0);
% G = RemovePolyBG(G, 0);

R = Normalize(R, 150, 190, 'M');
G = Normalize(G, 1560, 1600, 'M');


plotRaman(R, 0.0, 514);        
plotRaman(G, 0.0, 514);        
% plotRaman(DD, 0, 514);        


%% CRISTIAN AND ALENADRO
               
R = {
        DATA_20250211.S1BH514R
        DATA_20250211.S1CH514R
        DATA_20250211.S1DH514R
        DATA_20250211.S2BH514R
        DATA_20250211.S2CH514R
        DATA_20250211.S2DH514R
        
            };   
        
               
G = {
        DATA_20250211.S1BH514G
        DATA_20250211.S1CH514G
        DATA_20250211.S1DH514G
        DATA_20250211.S2BH514G
        DATA_20250211.S2CH514G
        DATA_20250211.S2DH514G
            };   



R = FlatFieldCorrection(R,DATA_20250211.FFH514R);
G = FlatFieldCorrection(G,DATA_20250211.FFH514G);

% 
R = FilterDataByXRange(R, 130, 220);
G = FilterDataByXRange(G, 1540, 1620);
%         
%         
R = RemovePolyBG(R, 1);
% G = RemovePolyBG(G, 1);
% DD = RemovePolyBG(DD, 0);
% 
R = Normalize(R, 150, 190, 'M');
G = Normalize(G, 1560, 1600, 'M');
% DD = Normalize(DD, 0, 3000, 'M');
% 
plotRaman(R, 0.0, 514);        
% plotRaman(G, 0.0, 514);        
% plotRaman(DD, 0, 514);        


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


