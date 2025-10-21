clc;
clear;
close all;

import UsefulFunctions.*;
rootpath = 'X:\Measurements Data\Absorption\';

%All paths as default

Dispersions = [rootpath,'20250813\Dispersions_TTF_TCNQ.csv'];
CF_Dispersions = [rootpath,'20250818\24hCF_Dispersions_TTF_TCNQ.csv'];
DGU_Dispersions = [rootpath,'20250819\DGU_Dispersions_TTF_TCNQ.csv'];

References = [rootpath,'References.csv'];
TCNQ_Chloroform = [rootpath,'20250617\TCNQ_Chloroform.csv'];
VacuumInfilSamples = [rootpath,'20251015\SalomeVI_CF_Samples_NIFUB_TCNQinDOCD2O.csv'];



%Select the paths of interest
paths = {
    References
    Dispersions
    CF_Dispersions
    TCNQ_Chloroform
    DGU_Dispersions
    VacuumInfilSamples
};
    

%% Read and structure data from the paths
ReadAbsorptionFromPaths(paths);

%% Labeling

DATA_20250813.Baseline.N='Baseline';
DATA_20250813.B3_10.N='Batch 3 - Annealed P2-SWCNTs';
DATA_20250813.R17_10.N='R17 - TTF@P2-SWCNTs (GP) (Rinsed)';
DATA_20250813.R21D_10.N='R21D - TTF@P2-SWCNTs (Melt) (Rinsed)';
DATA_20250813.R21D_5.N='R21D - TTF@P2-SWCNTs (Melt) (Rinsed)';
DATA_20250813.U22_10.N='U22 - TTF@P2-SWCNTs (Reflux) (Unrinsed)';
DATA_20250813.R18E_10.N='R18E - TCNQ@P2-SWCNTs (GP) (Rinsed)';
DATA_20250813.R23F_10.N='R23F - TCNQ@P2-SWCNTs (GP) (Rinsed)';

DATA_20250818.B3_24hCF_5.N='Batch 3 - 24hCF Annealed P2-SWCNTs';
DATA_20250818.B3_24hCF_10.N='Batch 3 - 24hCF Annealed P2-SWCNTs';
DATA_20250818.R17_24hCF_5.N='R17 - 24hCF TTF@P2-SWCNTs (GP) (Rinsed)';
DATA_20250818.R21D_24hCF_1.N='R21D - 24hCF TTF@P2-SWCNTs (Melt) (Rinsed)';
DATA_20250818.U22_24hCF_5.N='U22 - 24hCF TTF@P2-SWCNTs (Reflux) (Unrinsed)';
DATA_20250818.R18E_24hCF_5.N='R18E - 24hCF TCNQ@P2-SWCNTs (GP) (Rinsed)';
DATA_20250818.R23F_24hCF_5.N='R23F - 24hCF TCNQ@P2-SWCNTs (GP) (Rinsed)';
DATA_20250818.Baseline.N='Baseline';

DATA_20250818.B3_24hCF_5.N='(undoped) P2-SWCNTs';
DATA_20250818.B3_24hCF_10.N='(undoped) P2-SWCNTs';
DATA_20250818.R17_24hCF_5.N='TTF@P2-SWCNTs (Gas Phase)';
DATA_20250818.R21D_24hCF_1.N='TTF@P2-SWCNTs (Melt)';
DATA_20250818.U22_24hCF_5.N='TTF@P2-SWCNTs (Reflux in MeOH)';
DATA_20250818.R18E_24hCF_5.N='TCNQ@P2-SWCNTs (Gas Phase) A';
DATA_20250818.R23F_24hCF_5.N='TCNQ@P2-SWCNTs (Gas Phase) B';
DATA_20250818.Baseline.N='Baseline';


DATA_20250819.R23F_DGU_A_1.N='R23F - TCNQ@P2-SWCNTs (GP) (DGU A)';
DATA_20250819.R23F_DGU_B_1.N='R23F - TCNQ@P2-SWCNTs (GP) (DGU B)';
DATA_20250819.R23F_DGU_C_1.N='R23F - TCNQ@P2-SWCNTs (GP) (DGU C)';
DATA_20250819.R23F_DGU_D_1.N='R23F - TCNQ@P2-SWCNTs (GP) (DGU D)';
DATA_20250819.R23F_DGU_E_1.N='R23F - TCNQ@P2-SWCNTs (GP) (DGU E)';
DATA_20250819.R23F_DGU_F_1.N='R23F - TCNQ@P2-SWCNTs (GP) (DGU F)';
DATA_20250819.R23F_DGU_G_1.N='R23F - TCNQ@P2-SWCNTs (GP) (DGU G)';
DATA_20250819.U22_DGU_A_1.N='U22 - TTF@P2-SWCNTs (R) (DGU A)';
DATA_20250819.U22_DGU_B_1.N='U22 - TTF@P2-SWCNTs (R) (DGU B)';
DATA_20250819.U22_DGU_C_1.N='U22 - TTF@P2-SWCNTs (R) (DGU C)';
DATA_20250819.U22_DGU_D_1.N='U22 - TTF@P2-SWCNTs (R) (DGU D)';
DATA_20250819.U22_DGU_E_1.N='U22 - TTF@P2-SWCNTs (R) (DGU E)';
DATA_20250819.U22_DGU_F_1.N='U22 - TTF@P2-SWCNTs (R) (DGU F)';
DATA_20250819.U22_DGU_G_1.N='U22 - TTF@P2-SWCNTs (R) (DGU G)';
DATA_20250819.R17_DGU_A_1.N='R17 - TTF@P2-SWCNTs (GP) (DGU A)';
DATA_20250819.R17_DGU_B_1.N='R17 - TTF@P2-SWCNTs (GP) (DGU B)';
DATA_20250819.R17_DGU_C_1.N='R17 - TTF@P2-SWCNTs (GP) (DGU C)';
DATA_20250819.R17_DGU_D_1.N='R17 - TTF@P2-SWCNTs (GP) (DGU D)';
DATA_20250819.R17_DGU_E_1.N='R17 - TTF@P2-SWCNTs (GP) (DGU E)';
DATA_20250819.R17_DGU_F_1.N='R17 - TTF@P2-SWCNTs (GP) (DGU F)';
DATA_20250819.R17_DGU_G_1.N='R17 - TTF@P2-SWCNTs (GP) (DGU G)';
DATA_20250819.R21D_DGU_A_1.N='R21D - TTF@P2-SWCNTs (M) (DGU A)';
DATA_20250819.R21D_DGU_B_1.N='R21D - TTF@P2-SWCNTs (M) (DGU B)';
DATA_20250819.R21D_DGU_C_1.N='R21D - TTF@P2-SWCNTs (M) (DGU C)';
DATA_20250819.R21D_DGU_D_1.N='R21D - TTF@P2-SWCNTs (M) (DGU D)';
DATA_20250819.R21D_DGU_E_1.N='R21D - TTF@P2-SWCNTs (M) (DGU E)';
DATA_20250819.R21D_DGU_F_1.N='R21D - TTF@P2-SWCNTs (M) (DGU F)';
DATA_20250819.R21D_DGU_G_1.N='R21D - TTF@P2-SWCNTs (M) (DGU G)';


%% Corrections

% DATA_20250818.B3_24hCF_5.Y = DATA_20250818.B3_24hCF_5.Y - ()*DATA
% DATA_20250818.R17_24hCF_5.Y = DATA_20250818.R17_24hCF_5.Y
% DATA_20250818.R21D_24hCF_1.Y = DATA_20250818.R21D_24hCF_1.Y
% DATA_20250818.U22_24hCF_5.Y = DATA_20250818.U22_24hCF_5.Y
% DATA_20250818.R18E_24hCF_5.Y = DATA_20250818.R18E_24hCF_5.Y
% DATA_20250818.R23F_24hCF_5.Y = DATA_20250818.R23F_24hCF_5.Y

%% Plotting

%% RAW
All = {
%         DATA_20250813.Baseline
        DATA_20250813.B3_10
        
        DATA_20250813.R17_10
        DATA_20250813.R21D_5
        DATA_20250813.U22_10

        DATA_20250813.R18E_10
        DATA_20250813.R23F_10
       };
% plotAbsorption(All, 0.0)

% All = FilterDataByXRange(All, 220, 2500);
% All = BackgroundSubtraction(All, [610, 1330]);

% All = RemovePolyBG(All, 0)
% All = RemoveBackgroundProfile(All, [830, 1350, 2500]);
% All = RemoveBackgroundProfile(All, [625, 860, 1320, 2500]);
% All = RemoveBackgroundProfile(All, [300, 615, 870, 1300]);

% All = matchSpectra(All, 900, 5);
% All = matchSpectra(All, 1199, 5);
% All = Abs_MatchSpectra(All, [1350,1450]);
% All = Normalize(All, 1000, 1200, 'M');
% All = Normalize(All, 200, 1320, 'M');
% 
% plotAbsorption(All, 0.0)

%% CF

DATA_References.WaterInD2O = FilterDataByXRange(DATA_References.WaterInD2O, 240, 2400);
DATA_20250818.B3_24hCF_10 = FilterDataByXRange(DATA_20250818.B3_24hCF_10, 240, 2400);
DATA_20250818.R17_24hCF_5 = FilterDataByXRange(DATA_20250818.R17_24hCF_5, 240, 2400);
DATA_20250818.R21D_24hCF_1 = FilterDataByXRange(DATA_20250818.R21D_24hCF_1, 240, 2400);
DATA_20250818.U22_24hCF_5 = FilterDataByXRange(DATA_20250818.U22_24hCF_5, 240, 2400);
DATA_20250818.R18E_24hCF_5 = FilterDataByXRange(DATA_20250818.R18E_24hCF_5, 240, 2400);
DATA_20250818.R23F_24hCF_5 = FilterDataByXRange(DATA_20250818.R23F_24hCF_5, 240, 2400);
DATA_References.WaterFilled= FilterDataByXRange(DATA_References.WaterFilled, 240, 2400);

DATA_20251015.HS1_22hCF_2= FilterDataByXRange(DATA_20251015.HS1_22hCF_2, 240, 2400);
DATA_20251015.SF1_22hCF_2= FilterDataByXRange(DATA_20251015.SF1_22hCF_2, 240, 2400);
        
%Data correction
DATA_20250818.B3_24hCF_10.Y = DATA_20250818.B3_24hCF_10.Y - (DATA_References.WaterInD2O.Y)*0.12;
DATA_20250818.R17_24hCF_5.Y = DATA_20250818.R17_24hCF_5.Y - (DATA_References.WaterInD2O.Y)*0.12;
DATA_20250818.R21D_24hCF_1.Y = DATA_20250818.R21D_24hCF_1.Y - (DATA_References.WaterInD2O.Y)*0.375;
DATA_20250818.U22_24hCF_5.Y = DATA_20250818.U22_24hCF_5.Y - (DATA_References.WaterInD2O.Y)*0.12;
DATA_20250818.R18E_24hCF_5.Y = DATA_20250818.R18E_24hCF_5.Y - (DATA_References.WaterInD2O.Y)*0.12;
DATA_20250818.R23F_24hCF_5.Y = DATA_20250818.R23F_24hCF_5.Y - (DATA_References.WaterInD2O.Y)*0.12;

DATA_20251015.HS1_22hCF_2.Y = DATA_20251015.HS1_22hCF_2.Y + (DATA_References.WaterInD2O.Y)*0.15;
DATA_20251015.SF1_22hCF_2.Y = DATA_20251015.SF1_22hCF_2.Y + (DATA_References.WaterInD2O.Y)*0.15;


% DATA_20250818.B3_24hCF_10 = Abs_InverseLambdaFixed(DATA_20250818.B3_24hCF_10 , 80);
% DATA_20250818.R17_24hCF_5 = Abs_InverseLambdaFixed(DATA_20250818.R17_24hCF_5, 80);
% DATA_20250818.R21D_24hCF_1 = Abs_InverseLambdaFixed(DATA_20250818.R21D_24hCF_1, 80);
% DATA_20250818.U22_24hCF_5 = Abs_InverseLambdaFixed(DATA_20250818.U22_24hCF_5 , 80);
% DATA_20250818.R18E_24hCF_5 = Abs_InverseLambdaFixed(DATA_20250818.R18E_24hCF_5, 70);
% DATA_20250818.R23F_24hCF_5 = Abs_InverseLambdaFixed(DATA_20250818.R23F_24hCF_5, 70);
        
AllCF = {
%         DATA_20250813.Baseline
%         DATA_References.WaterInD2O
%         DATA_20250818.B3_24hCF_5

        DATA_References.WaterFilled
%         DATA_References.empty_P2_dial_0930
%         DATA_20250818.B3_24hCF_10
%         
        DATA_20250818.R17_24hCF_5
        DATA_20250818.R21D_24hCF_1
        DATA_20250818.U22_24hCF_5
        DATA_20251015.HS1_22hCF_2
        DATA_20251015.SF1_22hCF_2
% 
        DATA_20250818.R18E_24hCF_5
        DATA_20250818.R23F_24hCF_5
       };
% plotAbsorption(AllCF, 0.0)

% close all;
% AllCF = FilterDataByXRange(AllCF, 240, 2300);
% AllCF = BackgroundSubtraction(AllCF, [606, 1276]);
% AllCF = RemoveBackgroundProfile(AllCF, [300, 1200]);

% AllCF = RemovePolyBG(AllCF, 0)
% All = RemoveBackgroundProfile(All, [610, 960, 1255, 2400]);

% AllCF = RemoveBackgroundProfile(AllCF, [300, 615, 870, 1300]);
% AllCF = Abs_MatchSpectra(AllCF, [1350,1450]);

% AllCF = Normalize(AllCF, 250, 1200, 'M');
% AllCF = Normalize(AllCF, 914, 1200, 'M');
AllCF = Normalize(AllCF, 914, 1200, 'M');
% AllCF = RemovePolyBG(AllCF, 0)
% AllCF = Abs_InverseLambdaFixeddd(AllCF, 300)
plotAbsorption(AllCF, 0.0)

% DATA_20250617.TCNQ_Chloroform_30 = FilterDataByXRange(DATA_20250617.TCNQ_Chloroform_30, 240, 500);
% close all
% % hold on
% plot(DATA_20250617.TCNQ_Chloroform_30.X, DATA_20250617.TCNQ_Chloroform_30.Y, 'DisplayName', 'TCNQ', 'LineWidth', 1.0, 'Color', 'black')

%% DGU

DGU = {
%         DATA_20250819.R21D_DGU_A_1
%         DATA_20250819.R21D_DGU_B_1
%         DATA_20250819.R21D_DGU_C_1
%         DATA_20250819.R21D_DGU_D_1
%         DATA_20250819.R21D_DGU_E_1
%         DATA_20250819.R21D_DGU_F_1
%         DATA_20250819.R21D_DGU_G_1

%         DATA_20250819.R17_DGU_A_1
%         DATA_20250819.R17_DGU_B_1
%         DATA_20250819.R17_DGU_C_1
%         DATA_20250819.R17_DGU_D_1
%         DATA_20250819.R17_DGU_E_1
%         DATA_20250819.R17_DGU_F_1
%         DATA_20250819.R17_DGU_G_1
% 
        DATA_20250819.U22_DGU_A_1
        DATA_20250819.U22_DGU_B_1
        DATA_20250819.U22_DGU_C_1
        DATA_20250819.U22_DGU_D_1
        DATA_20250819.U22_DGU_E_1
        DATA_20250819.U22_DGU_F_1
        DATA_20250819.U22_DGU_G_1
% 
%         DATA_20250819.R23F_DGU_A_1
%         DATA_20250819.R23F_DGU_B_1
%         DATA_20250819.R23F_DGU_C_1
%         DATA_20250819.R23F_DGU_D_1
%         DATA_20250819.R23F_DGU_E_1
%         DATA_20250819.R23F_DGU_F_1
%         DATA_20250819.R23F_DGU_G_1
% 
%         DATA_References.WaterFilled
%         DATA_References.empty_P2_dial_0930
       };
   
% DGU = FilterDataByXRange(DGU, 350, 2500);
% DGU = BackgroundSubtraction(DGU, [610, 1330]);
% DGU = RemoveBackgroundProfile(DGU, [625, 860, 1320, 2500]);

% DGU = RemovePolyBG(DGU, 0)
% DGU = RemoveBackgroundProfile(DGU, [830, 1350, 2500]);

% DGU = RemoveBackgroundProfile(DGU, [300, 615, 870, 1300]);
% DGU = Abs_MatchSpectra(DGU, [1350,1450]);
% DGU = Normalize(DGU, 1000, 1200, 'M');
% DGU = Normalize(DGU, 920, 1100, 'M');
% 

% close all
% plotAbsorption(DGU, 0.05)
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SpectraList = Abs_MatchSpectra(SpectraList, range)
% ABS_MATCHSPECTRA Aligns spectra vertically by subtracting
% the average value in a given wavelength range.
%
% Inputs:
%   SpectraList - cell array of structs with fields X and Y
%   range       - [xMin, xMax] specifying averaging range
%
% Output:
%   SpectraList - vertically aligned spectra

    nSpectra = length(SpectraList);
    avgValues = zeros(nSpectra,1);

    % Calculate mean value in the given range for each spectrum
    for i = 1:nSpectra
        spec = SpectraList{i};
        mask = spec.X >= range(1) & spec.X <= range(2);
        avgValues(i) = mean(spec.Y(mask));
    end

    % Subtract each average from the respective spectrum
    for i = 1:nSpectra
        spec = SpectraList{i};
        spec.Y = spec.Y - avgValues(i);
        SpectraList{i} = spec;
    end
end

function SpectraList = Abs_InverseBackground(SpectraList)

    for i = 1:length(SpectraList)
        spec = SpectraList{i};
        X = spec.X(:);
        Y = spec.Y(:);

        % Objective: fit A/x + B to data (least squares)
        fitFunc = @(params) sum((Y - (params(1)./X + params(2))).^2);

        % Initial guess for A, B
        initParams = [1, mean(Y)];

        % Fit using fminsearch
        options = optimset('Display','off');
        params = fminsearch(fitFunc, initParams, options);

        % Subtract background
        background = params(1)./X + params(2);
        spec.Y = Y - background;

        SpectraList{i} = spec;
    end
end

function SpectraList = Abs_InverseLambdaFixeddd(SpectraList, Constant)
    
    for i = 1:length(SpectraList)
            spec = SpectraList{i};
            X = spec.X(:);
            Y = spec.Y(:);


            % Subtract background
            background = Constant./X;
            spec.Y = Y - background;
            
            SpectraList{i} = spec;
    end
end

function Spectrum = Abs_InverseLambdaFixed(Spectrum, Constant)

    X = Spectrum.X(:);
    Y = Spectrum.Y(:);

    % Subtract background
    background = Constant./X;
    Spectrum.Y = Y - background;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    %   - samplesToFilter: Either a single structure with fields 'X' and 'Y',
    %                      or a cell array of such structures.
    %   - xMin: The minimum value of X to include in the filtered data.
    %   - xMax: The maximum value of X to include in the filtered data.
    %
    % Outputs:
    %   - filteredSamples: Same format as input (single structure or cell array),
    %                      with filtered 'X' and 'Y' values within the range [xMin, xMax].

    % Ensure input is in cell array form
    isSingle = ~iscell(samplesToFilter);
    if isSingle
        samplesToFilter = {samplesToFilter};
    end
    
    % Allocate output
    filteredSamples = cell(size(samplesToFilter));
    
    % Iterate over each sample to filter
    for sampleIdx = 1:length(samplesToFilter)
        currentSample = samplesToFilter{sampleIdx};
        
        % Find the indices of X-values within the specified range
        validIndices = currentSample.X >= xMin & currentSample.X <= xMax;
        
        % Filter the data
        filteredSample = currentSample;
        filteredSample.X = currentSample.X(validIndices);
        filteredSample.Y = currentSample.Y(validIndices);
        
        % Store result
        filteredSamples{sampleIdx} = filteredSample;
    end
    
    % Return in the same format as input
    if isSingle
        filteredSamples = filteredSamples{1};
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
%         A = params(1);
%         B = params(2);
        A = 0;
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
