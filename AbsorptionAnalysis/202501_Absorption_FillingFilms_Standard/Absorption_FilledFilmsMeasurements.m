clc;
clear;
addpath('X:\SWCNTs');
import UsefulFunctions.*;
rootpath = 'X:\Measurements Data\Absorption\';
addpath 'X:\SWCNTs\SpecialMatlabFunctions\DrosteEffect-BrewerMap-3.2.5.0'

%All paths as default
KIT_FILM_Round1= [rootpath,'20241003\FilledSamplesFilms.csv'];
KIT_FILM_Round2 = [rootpath,'20241202\Films_KIT_Set2.csv'];

%Select the paths of interest
paths = {   
            KIT_FILM_Round1
            KIT_FILM_Round2
        };

%Read and structure data from the paths

ReadAbsorptionFromPaths(paths);


%% %--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DATA_20241003.FS2.N = 'PCE@SWCNT p-doping';
DATA_20241003.FS3.N = 'TCE@SWCNT p-doping';
DATA_20241003.FS4.N = 'TEMED@SWCNT n-doping';
DATA_20241003.FS5.N = 'TDAE@SWCNT n-doping';
DATA_20241003.FS6.N = 'Hexadecane@SWCNT neutral';
DATA_20241003.FS7.N = 'Dodecane@SWCNT neutral';

DATA_20241202.F1Dodecane.N = 'Dodecane@SWCNT neutral';
DATA_20241202.F2TDAE.N = 'TDAE@SWCNT n-doping';
DATA_20241202.F3TTF.N = 'TTF@SWCNT n-doping';
DATA_20241202.F5PCE.N = 'PCE@SWCNT p-doping';
DATA_20241202.F4TEMED.N = 'TEMED@SWCNT n-doping';

%% %--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %--------BACKGROUND CORRECTION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %--------SPECTRA SELECTION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


FR1 = {
            DATA_20241003.FS2
            DATA_20241003.FS3
            DATA_20241003.FS4
            DATA_20241003.FS5
            DATA_20241003.FS6
            DATA_20241003.FS7
    };

FR2 = {
            DATA_20241202.F2TDAE
            DATA_20241202.F1Dodecane
            DATA_20241202.F5PCE
            DATA_20241202.F4TEMED
            DATA_20241202.F3TTF
           
        };


FR1 = FilterDataByXRange(FR1, 200, 2600);
% % FR1 = NormalizeSample(FR1,902, 1300); 
% % FR1 = RemovePolyBG(FR1, 0);
% % FR1 = NormalizeSample(FR1,902, 1300); 
FR1 = BackgroundSubtraction(FR1, [500, 2600]);
FR1 = NormalizeSample(FR1,902, 1300); 
plotAbsorption(FR1, 0);

FR2 = FilterDataByXRange(FR2, 200, 2600);
% % % FR2 = NormalizeSample(FR2,902, 1300); 
FR2 = RemovePolyBG(FR2, 0);
% % % FR2 = NormalizeSample(FR2,902, 1300); 
FR2 = BackgroundSubtraction(FR2, [0, 2500]);
FR2 = NormalizeSample(FR2,902, 1300); 
plotAbsorption(FR2, 0);


All = {
            DATA_20241003.FS2
            DATA_20241202.F5PCE

%             DATA_20241003.FS3
            
            DATA_20241003.FS4
            DATA_20241202.F4TEMED
            
            DATA_20241003.FS5
            DATA_20241202.F2TDAE
            
%             DATA_20241003.FS6
            
            DATA_20241003.FS7
            DATA_20241202.F1Dodecane
            
%             DATA_20241202.F3TTF
            
        };

All = FilterDataByXRange(All, 200, 2500);
% % % FR2 = NormalizeSample(FR2,902, 1300); 
% All = RemovePolyBG(All, 0);
% All = Normalize(All,902, 1300,'M'); 
All = BackgroundSubtraction(All, [200, 2500]);
All = Normalize(All,800, 1200,'M'); 


% All{11}.Y = All{11}.Y + 0.01;

% plotAbsorption(All, 0);

DATA_20241003.FS4.N = 'TEMED@SWCNT (Round 1)';
DATA_20241202.F4TEMED.N = 'TEMED@SWCNT (Round 2)';
DATA_20241003.FS7.N = 'Dodecane@SWCNT (Round 1)';
DATA_20241202.F1Dodecane.N = 'Dodecane@SWCNT (Round 2)';
DATA_20241003.FS2.N = 'PCE@SWCNT (Round 1)';
DATA_20241202.F5PCE.N = 'PCE@SWCNT (Round 2)';
DATA_20241003.FS5.N = 'TDAE@SWCNT (Round 1)';
DATA_20241202.F2TDAE.N = 'TDAE@SWCNT (Round 2)';
            
Compare = {
            DATA_20241003.FS2
            DATA_20241202.F5PCE

%             DATA_20241003.FS3
            
            DATA_20241003.FS4
            DATA_20241202.F4TEMED
            
            DATA_20241003.FS5
            DATA_20241202.F2TDAE
            
%             DATA_20241003.FS6
            
            DATA_20241003.FS7
            DATA_20241202.F1Dodecane
            
%             DATA_20241202.F3TTF
            
        };
    

Compare = FilterDataByXRange(Compare, 200, 2500);
% % % FR2 = NormalizeSample(FR2,902, 1300); 
% All = RemovePolyBG(All, 0);
% All = Normalize(All,902, 1300,'M'); 
Compare = BackgroundSubtraction(Compare, [200, 2500]);
Compare = Normalize(Compare,800, 1200,'M'); 

offset = 2.25
Compare{2}.Y = Compare{2}.Y + offset;
Compare{4}.Y = Compare{4}.Y + offset;
Compare{6}.Y = Compare{6}.Y + offset;
Compare{8}.Y = Compare{8}.Y + offset;


plotAbsorption(Compare, 0);

FR2 = FilterDataByXRange(FR2, 200, 2500);
FR2 = RemovePolyBG(FR2, 0);
FR2 = Normalize(FR2,902, 1300,'M'); 

plotAbsorption(FR2,0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
