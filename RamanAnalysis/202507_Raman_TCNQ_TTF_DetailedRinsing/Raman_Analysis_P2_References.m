clc;
clear;
import UsefulFunctions.*;
rootpath = 'X:\Measurements Data\Raman\';
addpath('X:\SWCNTs\RamanAnalysis\Raman - Voigt Fitting\Faddeeva_voigt');
addpath('X:\SWCNTs');

% rootpath = 'X:\Measurements Data\Raman\';
%All paths as default

path_a = [rootpath,'20250131\'];
path_b = [rootpath,'20250403\'];
path_c = [rootpath,'20250520\'];
path_d = [rootpath,'20250618\'];

%Select the paths of interest

paths = {
        path_a
        path_b
        path_c
        path_d
        };


ReadRamanFromPaths(paths, 2);

%% %--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% TCNQ ANALYSIS

% G and 2D Band Plots
G = {
        DATA_20250131.PP2L514G 
        DATA_20250403.PP2L514G 
        DATA_20250520.PP2514G
        DATA_20250618.PP2L514G 
        
%         DATA_20250520.U16B514G
%         DATA_20250520.U17L514G
%         DATA_20250618.U18L514G
%         DATA_20250618.U20L514G
     };   
%         
G = FilterDataByXRange(G, 1250, 1680);           
G = RemovePolyBG(G, 0);
G = SubstractLinearBG(G, 1250, 1680);
G = Normalize(G, 1500, 1680, 'M');

plotRaman(G, 0.00, 514);     

DD = {
        DATA_20250131.PP2L514D
        DATA_20250403.PP2L514D
        DATA_20250520.PP2L514D
        DATA_20250618.PP2L514D
        
%         DATA_20250520.U17L514D
%         DATA_20250618.U17L514D
%         DATA_20250618.U20L514D
        
          };   
        
DD = FilterDataByXRange(DD, 2500, 2835);           
DD = RemovePolyBG(DD, 0);
DD = SubstractLinearBG(DD, 2500, 2835);
DD = Normalize(DD, 2500, 2700, 'M');

% close all
plotRaman(DD, 0.2, 514);     


%% DOPING VECTORS


FITTEDD = FitSamples(DD, 2680)

FD = zeros(1, length(FITTEDD));
for i=1:length(FITTEDD)
%    FITTED{i}.N
   FD(i) = FITTEDD{i}.F.Params(2)-FITTEDD{1}.F.Params(2);
   FD(i) = FITTEDD{i}.F.Params(2);
end


% FITTED
FITTEDG = FitSamples(G, 1592)
FG = zeros(1, length(FITTEDG));
NG = cell(1, length(FITTEDG));

for i=1:length(FITTEDG)
%    FITTED{i}.N
   FG(i) = FITTEDG{i}.F.Params(2)-FITTEDG{1}.F.Params(2);
   FG(i) = FITTEDG{i}.F.Params(2);
   NG{i} = num2str(FITTEDG{i}.N);
end

plot(FG,FD)
scatter(FG,FD, 50, 'k','d', 'filled')
text(FG, FD, NG, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 8);
xlim([1590, 1593]);
ylim([2660, 2697]);
ylabel('2D Band (cm-1)')
xlabel('G Band (cm-1)')

%%%% FITTING V2

FITTEDG = FitSamples(G, 1590);
FITTEDD = FitSamples(DD, 2670)

numSamples = length(FITTEDG);
FG = zeros(1, numSamples);
FD = zeros(1, numSamples);
NG = cell(1, numSamples);

% Extract data
for i = 1:numSamples
    FG(i) = FITTEDG{i}.F.Params(2);
    FD(i) = FITTEDD{i}.F.Params(2);  % Assuming you want FD values too
    NG{i} = num2str(FITTEDG{i}.N);   % Sample names
end

% Assign unique colors
colors = lines(numSamples);  % Or use: jet(numSamples), hsv(numSamples), etc.

figure; hold on;
for i = 1:numSamples
    scatter(FG(i), FD(i), 50, colors(i,:), 'd', 'filled', 'DisplayName', NG{i});
end

xlim([1590, 1593]);
ylim([2660, 2697]);
ylabel('2D Band (cm^{-1})');
xlabel('G Band (cm^{-1})');
legend();  % Adjust legend location as needed
title('Doping Vectors');







%% 514nm DATA - 20250131                
% G = {
%           DATA_20250520.PP2514G
% %           DATA_20250520.U16L514G
% %           DATA_20250520.U16B514G
%             DATA_20250520.R12L514G
% 
%             DATA_20250520.R16L514G  
%           
% %           DATA_20250520.U17L514G          
% 
%           DATA_20250520.R17L514G
%           DATA_20250520.R13L514G
%           
% 
% 
% %           DATA_20250520.PTTF514G
% 
%             };   
% DD = {
%           DATA_20250520.PP2L514D
%                     
%           DATA_20250520.U16L514D
% 
%           DATA_20250520.R16L514D
%          DATA_20250520.R12L514D
% 
%           DATA_20250520.U17L514D
% 
%           DATA_20250520.R17L514D   
%           DATA_20250520.R13L514D
% 
%             };  
% 
% R = {
%         DATA_20250520.PP2H514R
%         DATA_20250520.R12H514R
%         DATA_20250520.R13H514R
%         DATA_20250520.R16H514R
%         DATA_20250520.R17H514R
%     };  
% 
% RL = {
%         DATA_20250520.PP2L514R
%         DATA_20250520.R12L514R
%         DATA_20250520.R13L514R
%         DATA_20250520.R16L514R
%         DATA_20250520.R17L514R
% 
%     };  
% 
% R = FlatFieldCorrection(R,DATA_20250520.FFH514R);
% 
% R = FilterDataByXRange(R, 135, 200);
% RL = FilterDataByXRange(RL, 100, 250);
% 
% % G = FilterDataByXRange(G, 1550, 1680);
% G = FilterDataByXRange(G, 1300, 1680);
% 
% DD = FilterDataByXRange(DD, 2500, 2840);
%               
% G = RemovePolyBG(G, 0);
% DD = RemovePolyBG(DD, 0);
% 
% % G = SubstractLinearBG(G, 1250, 1680);
% DD = SubstractLinearBG(DD, 2500,2835);
% RL = SubstractLinearBG(RL, 100, 250);
% 
% R = Normalize(R, 160, 170, 'M');
% RL = Normalize(RL, 160, 170, 'M');
% 
% G = Normalize(G, 1580, 1600, 'M');
% DD = Normalize(DD, 0, 2680, 'M');
% % 
% 
% % 
% % plotRaman(RL, 0.5, 514);  
% % plotRaman(R, 0.2, 514); 
% % G{6}.Y = G{6}.Y/100
% % % % % % % % % % plotRaman(G, 0.00, 514);     
% plotRaman({G{4}}, 0.0, 514);        

% plotRaman(DD, 1, 514);        

%% DOPING VECTORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FITTED
% FITTEDD = FitSamples(DD, 2680)
% 
% FD = zeros(1, length(FITTEDD));
% for i=1:length(FITTEDD)
% %    FITTED{i}.N
%    FD(i) = FITTEDD{i}.F.Params(2)-FITTEDD{1}.F.Params(2);
%    FD(i) = FITTEDD{i}.F.Params(2);
% end
% 
% 
% % FITTED
% FITTEDG = FitSamples(G, 1592)
% FG = zeros(1, length(FITTEDG));
% NG = cell(1, length(FITTEDG));
% 
% for i=1:length(FITTEDG)
% %    FITTED{i}.N
%    FG(i) = FITTEDG{i}.F.Params(2)-FITTEDG{1}.F.Params(2);
%    FG(i) = FITTEDG{i}.F.Params(2);
%    NG{i} = num2str(FITTEDG{i}.N);
% end
% 
% % plot(FG,FD)
% scatter(FG,FD, 50, 'k','d', 'filled')
% text(FG, FD, NG, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 8);
% xlim([1590,1593])
% ylim([2670,2690])
% ylabel('2D Band (cm-1)')
% xlabel('G Band (cm-1)')

%%%%% FITTING V2
% 
% FITTEDG = FitSamples(G, 1592);
% FITTEDD = FitSamples(DD, 2680)
% 
% numSamples = length(FITTEDG);
% FG = zeros(1, numSamples);
% FD = zeros(1, numSamples);
% NG = cell(1, numSamples);
% 
% % Extract data
% for i = 1:numSamples
%     FG(i) = FITTEDG{i}.F.Params(2);
%     FD(i) = FITTEDD{i}.F.Params(2);  % Assuming you want FD values too
%     NG{i} = num2str(FITTEDG{i}.N);   % Sample names
% end
% 
% % Assign unique colors
% colors = lines(numSamples);  % Or use: jet(numSamples), hsv(numSamples), etc.
% 
% figure; hold on;
% for i = 1:numSamples
%     scatter(FG(i), FD(i), 50, colors(i,:), 'd', 'filled', 'DisplayName', NG{i});
% end
% 
% xlim([1590, 1593]);
% ylim([2670, 2690]);
% ylabel('2D Band (cm^{-1})');
% xlabel('G Band (cm^{-1})');
% legend();  % Adjust legend location as needed
% title('Scatter Plot with Unique Colors and Legend');



%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
%         disp(['Structure ', num2str(i), ' Polynomial Coefficients: ', num2str(p)]);
        
        % Save the updated structure back to the list
        DSList{i} = DS;
    end
end


