clc;
clear;
addpath('X:\SWCNTs');
import UsefulFunctions.*;
rootpath = 'X:\Measurements Data\Raman\';

% rootpath = 'X:\Measurements Data\Raman\';
addpath('X:\SWCNTs\SpecialMatlabFunctions\DrosteEffect-BrewerMap-3.2.5.0')
%All paths as default
path_powder = [rootpath,'20250131\'];
path_dgua = [rootpath,'20240517\'];
path_dgub = [rootpath,'20241212\'];
path_sf = [rootpath,'20240111\'];
path_ga = [rootpath,'20240610\'];

%Select the paths of interest

paths = {
        path_powder
        path_dgua
        path_dgub
        path_sf
        
        path_ga
        };


ReadRamanFromPaths(paths, 2);

%%%--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DATA_20250131.S10L514R.N='Dispersion S10 - TMG@P2-SWCNTs';
DATA_20250131.S11L514R.N='Dispersion S11 - TTF@P2-SWCNTs';
DATA_20250131.SR0L514R.N='Dispersion SR0 - Empty@P2-SWCNTs';
DATA_20250131.SR1L514R.N='Dispersion SR1 - D2O@P2-SWCNTs - KVD';
DATA_20250131.SR2L514R.N='Dispersion SR2 - MeOH@P2-SWCNTs';
DATA_20250131.SWFL514R.N='Dispersion SWF - D2O@P2-SWCNTs (Salome)';
DATA_20250131.P12L514R.N='Powder S12 - TTF@P2-SWCNTs';
DATA_20250131.P12RBM.N='Powder S12 - TTF@P2-SWCNTs';
DATA_20250131.P2ARBM.N='Powder P2-SWCNTs Annealed';
DATA_20250131.PAPL514R.N='Powder AP-SWCNTs';
DATA_20250131.PP2L514R.N='Powder P2-SWCNTs Annealed';
DATA_20240517.EAH514R.N='Empty Arc SWCNTs';
DATA_20240517.S2H514R.N='CB PCE@SWCNT Dial. DGU C (Filled)';
DATA_20240517.S3H514R.N='CB TCE@SWCNT Dial. DGU C (Filled)';
DATA_20240517.S4H514R.N='CB TEMED@SWCNT Dial. DGU C (Filled)';
DATA_20240517.S5H514R.N='CB TDAE@SWCNT Dial. DGU C (Filled)';
DATA_20240517.S6H514R.N='CB Hexadecane@SWCNT Dial. DGU C (Filled)';
DATA_20240517.S7H514R.N='CB Dodecane@SWCNT Dial. DGU C (Filled)';
DATA_20240517.WAH514R.N='Water Filled Arc SWCNTs';
DATA_20241212.S1BH514R.N='TEMED Filled DGU (SFF1B)';
DATA_20241212.SF6H514R.N='TTF Filled DGU (SFF6)';
DATA_20241212.SR0H514R.N='SF6 Sofie Empty';
DATA_20241212.SR1H514R.N='Water Filled DGU (KVDP2)';
DATA_20240111.S240111J.N='SF D2O@SWCNT 514nm';
DATA_20240111.S240111K.N='SF TCE@SWCNT 514nm';
DATA_20240111.S240111KK.N='SF TCE@SWCNT 514nm';
DATA_20240111.S240111L.N='SF Methanol@SWCNT 514nm';
DATA_20240111.S240111M.N='SF TCE@SWCNT 514nm';
DATA_20240111.S240111N.N='SF TTF@SWCNT 514nm';
DATA_20240111.S240111O.N='SF PCE@SWCNT 514nm';
DATA_20240111.S240111P.N='SF PCE@SWCNT 514nm';
DATA_20240111.S240111Q.N='SF PCE@SWCNT 514nm';
DATA_20240111.S240111R.N='SF TEMED@SWCNT 514nm';
DATA_20240111.S240111S.N='SC Empty@SWCNT 514nm ';


DATA_20240517.EAL514GD.N = 'Empty@SWCNTs'
DATA_20240517.S2L514GD.N = 'PCE@SWCNTs'
DATA_20240517.S3L514GD.N = 'TCE@SWCNTs'
DATA_20240517.S4L514GD.N = 'TEMED@SWCNTs'
DATA_20240517.S5L514GD.N = 'TDAE@SWCNTs'
DATA_20240517.S6L514GD.N = 'Hexadecane@SWCNTs'
DATA_20240517.S7L514GD.N = 'Dodecane@SWCNTs'
DATA_20240517.WAL514GD.N = 'H2O@SWCNTs'

DATA_20240111.S240111A.N='SF D2O@SWCNTs';
DATA_20240111.S240111B.N='SF TCE@SWCNTs';
DATA_20240111.S240111BB.N='SF TCE@SWCNTs';
DATA_20240111.S240111C.N='SF MeOH@SWCNTs';
DATA_20240111.S240111D.N='SF TCE@SWCNTs';
DATA_20240111.S240111E.N='SF TTF@SWCNTs';
DATA_20240111.S240111F.N='SF PCE@SWCNTs';
DATA_20240111.S240111G.N='SF PCE@SWCNTs';
DATA_20240111.S240111H.N='SF PCE@SWCNs';
DATA_20240111.S240111I.N='SF TEMED@SWCNTs';


%%%--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%--------SAMPLE COMPARISION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Centrifuged Samples
               
R = {
%         DATA_20250131.SR0L514R
%         DATA_20250131.SR1L514R
% %         DATA_20250131.SWFL514R
% %         DATA_20250131.SR2L514R
% %         DATA_20250131.S10L514R
%         DATA_20250131.S11L514R

%         DATA_20240517.WAH514R
% %         DATA_20240517.EAH514R
        DATA_20240517.S2H514R
        DATA_20240517.S3H514R
        DATA_20240517.S4H514R
        DATA_20240517.S5H514R
        DATA_20240517.S6H514R
        DATA_20240517.S7H514R

         
%         DATA_20241212.S1BH514R
%         DATA_20241212.SF6H514R
%         DATA_20241212.SR0H514R
%         DATA_20241212.SR1H514R

%         DATA_20240111.S240111J
%         DATA_20240111.S240111K
%         DATA_20240111.S240111KK
%         DATA_20240111.S240111L
%         DATA_20240111.S240111M
%         DATA_20240111.S240111N
%         DATA_20240111.S240111O
%         DATA_20240111.S240111P
%         DATA_20240111.S240111Q
%         DATA_20240111.S240111R
%         DATA_20240610.EAH514R



    };   
     

% R = FlatFieldCorrection(R,DATA_20250131.FFL514R);
% R = FilterDataByXRange(R, 100, 250);           
% R = RemovePolyBG(R, 1);
% R = Normalize(R, 100, 250, 'M');
% plotRaman(R, 0.0, 514);        

%MY SAMPLES AFTER DGU
G = {
        DATA_20240517.EAL514GD
        DATA_20240517.S2L514GD
        DATA_20240517.S3L514GD
        DATA_20240517.S4L514GD
        DATA_20240517.S5L514GD
        DATA_20240517.S6L514GD
        DATA_20240517.S7L514GD
        DATA_20240517.WAL514GD
%         DATA_20240517.WAH514G
%         DATA_20240517.EAH514G
%         DATA_20240517.S2H514G
%         DATA_20240517.S3H514G
%         DATA_20240517.S4H514G
%         DATA_20240517.S5H514G
%         DATA_20240517.S6H514G
%         DATA_20240517.S7H514G
    };   

%SALOME SAMPLES G BAND AFTER CF but not DGU
% G = {
%         DATA_20240111.S240111A
%         DATA_20240111.S240111B
%         DATA_20240111.S240111BB
%         DATA_20240111.S240111C
%         DATA_20240111.S240111D
%         DATA_20240111.S240111E
%         DATA_20240111.S240111F
%         DATA_20240111.S240111G
%         DATA_20240111.S240111H
%         DATA_20240111.S240111I
%               
%     };   

% G = FlatFieldCorrection(G,DATA_20250131.FFL514R);
G = FilterDataByXRange(G, 1260, 1660);           
G = RemovePolyBG(G, 0);
G = Normalize(G, 1260, 1680, 'M');
plotRaman(G, 0.02, 514);    
% plotMaxima(G, 1575, 1610, 0.05)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


