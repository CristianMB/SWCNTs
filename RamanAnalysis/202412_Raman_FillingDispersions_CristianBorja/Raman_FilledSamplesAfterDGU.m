clc;
clear;
addpath('X:\SWCNTs');
import UsefulFunctions.*;
rootpath = 'X:\Measurements Data\Raman\';
addpath('X:\SWCNTs\SpecialMatlabFunctions\DrosteEffect-BrewerMap-3.2.5.0')
%All paths as default
path_20231206 = [rootpath,'20231206\'];

path_20240320 = [rootpath,'20240320\'];
path_20240321 = [rootpath,'20240321\'];
path_20240325 = [rootpath,'20240325\'];
path_20240514 = [rootpath,'20240514\'];
path_20240517 = [rootpath,'20240517\'];
path_20240610 = [rootpath,'20240610\'];
path_20240612 = [rootpath,'20240612\'];

%Select the paths of interest

paths = {
        path_20240320
        path_20240321
        path_20240325
        path_20240514
        path_20240517
        path_20240610
        path_20240612
        path_20231206
        };


ReadRamanFromPaths(paths, 2);

%%%--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DATA_20240517.EAH514R.N='Empty Arc SWCNTs';
DATA_20240517.EAL514GD.N='Empty Arc SWCNTs';
DATA_20240517.S2H514R.N='CB PCE@SWCNT Dial. DGU C (Filled)';
DATA_20240517.S2L514GD.N='CB PCE@SWCNT Dial. DGU C (Filled)';
DATA_20240517.S3H514R.N='CB TCE@SWCNT Dial. DGU C (Filled)';
DATA_20240517.S3L514GD.N='CB TCE@SWCNT Dial. DGU C (Filled)';
DATA_20240517.S4H514R.N='CB TEMED@SWCNT Dial. DGU C (Filled)';
DATA_20240517.S4L514GD.N='CB TEMED@SWCNT Dial. DGU C (Filled)';
DATA_20240517.S5H514R.N='CB TDAE@SWCNT Dial. DGU C (Filled)';
DATA_20240517.S5L514GD.N='CB TDAE@SWCNT Dial. DGU C (Filled)';
DATA_20240517.S6H514R.N='CB Hexadecane@SWCNT Dial. DGU C (Filled)';
DATA_20240517.S6L514GD.N='CB Hexadecane@SWCNT Dial. DGU C (Filled)';
DATA_20240517.S7H514R.N='CB Dodecane@SWCNT Dial. DGU C (Filled)';
DATA_20240517.S7L514GD.N='CB Dodecane@SWCNT Dial. DGU C (Filled)';
DATA_20240517.WAH514R.N='Water Filled Arc SWCNTs';
DATA_20240517.WAL514GD.N='Water Filled Arc SWCNTs';
DATA_20240320.S2H520R.N='CB PCE@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240320.S3H520R.N='CB TCE@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240320.S4H520R.N='CB TEMED@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240320.S5H520R.N='CB TDAE@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240320.S6H520R.N='CB Hexadecane@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240320.S7H520R.N='CB Dodecane@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240321.S2H520G.N='CB PCE@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240321.S2L520G.N='CB PCE@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240321.S3H520G.N='CB TCE@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240321.S3L520G.N='CB TCE@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240321.S4H520G.N='CB TEMED@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240321.S4L520G.N='CB TEMED@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240321.S5H520G.N='CB TDAE@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240321.S5L520G.N='CB TDAE@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240321.S6H520G.N='CB Hexadecane@SWCNT Dial. DGU C (Filled) 520,25nm';
DATA_20240321.S6L520G.N='CB Hexadecane@SWCNT Dial. DGU C (Filled) 520,25nm';
DATA_20240321.S7H520G.N='CB Dodecane@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240321.S7L520G.N='CB Dodecane@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240514.EAL650D.N='Empty Arc SWCNTs';
DATA_20240514.EAL650G.N='Empty Arc SWCNTs';
DATA_20240514.EAL650R.N='Empty Arc SWCNTs';
DATA_20240514.S2L650G.N='CB PCE@SWCNT Dial. DGU C (Filled)';
DATA_20240514.S2L650R.N='CB PCE@SWCNT Dial. DGU C (Filled)';
DATA_20240514.S3L650D.N='CB TCE@SWCNT Dial. DGU C (Filled)';
DATA_20240514.S3L650G.N='CB TCE@SWCNT Dial. DGU C (Filled)';
DATA_20240514.S3L650R.N='CB TCE@SWCNT Dial. DGU C (Filled)';
DATA_20240514.S4L650D.N='CB TEMED@SWCNT Dial. DGU C (Filled)';
DATA_20240514.S4L650G.N='CB TEMED@SWCNT Dial. DGU C (Filled)';
DATA_20240514.S4L650R.N='CB TEMED@SWCNT Dial. DGU C (Filled)';
DATA_20240514.S5L650D.N='CB TDAE@SWCNT Dial. DGU C (Filled)';
DATA_20240514.S5L650G.N='CB TDAE@SWCNT Dial. DGU C (Filled)';
DATA_20240514.S5L650R.N='CB TDAE@SWCNT Dial. DGU C (Filled)';
DATA_20240514.S6L650D.N='CB Hexadecane@SWCNT Dial. DGU C (Filled) 650,06nm';
DATA_20240514.S6L650G.N='CB Hexadecane@SWCNT Dial. DGU C (Filled) 650,06nm';
DATA_20240514.S6L650R.N='CB Hexadecane@SWCNT Dial. DGU C (Filled) 650,06nm';
DATA_20240514.S7L650D.N='CB Dodecane@SWCNT Dial. DGU C (Filled)';
DATA_20240514.S7L650G.N='CB Dodecane@SWCNT Dial. DGU C (Filled)';
DATA_20240514.S7L650R.N='CB Dodecane@SWCNT Dial. DGU C (Filled)';
DATA_20240514.WAL650D.N='Water Filled Arc SWCNTs';
DATA_20240514.WAL650G.N='Water Filled Arc SWCNTs';
DATA_20240514.WAL650R.N='Water Filled Arc SWCNTs';
DATA_20240514.S2L650D.N='CB PCE@SWCNT Dial. DGU C (Filled)';
DATA_20240325.REH680R.N='Empty Reference';
DATA_20240325.RWH680R.N='Water Filled Reference';
DATA_20240325.S2H680R.N='CB PCE@SWCNT Dial. DGU C (Filled) 680nm';
DATA_20240325.S3H680R.N='CB TCE@SWCNT Dial. DGU C (Filled) 680nm';
DATA_20240325.S4H680R.N='CB TEMED@SWCNT Dial. DGU C (Filled) 680nm';
DATA_20240325.S5H680R.N='CB TDAE@SWCNT Dial. DGU C (Filled) 680nm';
DATA_20240325.S6H680R.N='CB Hexadecane@SWCNT Dial. DGU C (Filled) 680nm';
DATA_20240325.S7H680R.N='CB Dodecane@SWCNT Dial. DGU C (Filled) 680nm';

%%%--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%--------SAMPLE COMPARISION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 514nm DATA 
% S2-S6 + References

               
RBMsHD514 = {
                    DATA_20240517.EAH514R
                    DATA_20240517.S2H514R
                    DATA_20240517.S3H514R
                    DATA_20240517.S4H514R
                    DATA_20240517.S5H514R
                    DATA_20240517.S6H514R
                    DATA_20240517.S7H514R
                    DATA_20240517.WAH514R
            };   
        
               
GDBand514 = {
                    DATA_20240517.EAL514GD
                    DATA_20240517.S2L514GD
                    DATA_20240517.S3L514GD
                    DATA_20240517.S4L514GD
                    DATA_20240517.S5L514GD
                    DATA_20240517.S6L514GD
                    DATA_20240517.S7L514GD
                    DATA_20240517.WAL514GD
            };   


RBMsHD514 = FlatFieldCorrection(RBMsHD514, DATA_20240517.FF514R)
RBMsHD514 = ClipSamples(RBMsHD514,140,240);      
RBMsHD514 = RemoveBaselinePolynomial(RBMsHD514, 0);
RBMsHD514 = RemoveLinearBackground(RBMsHD514, 141, 239);
RBMsHD514 = RemoveBaselinePolynomial(RBMsHD514, 0);
RBMsHD514 = NormalizeSample(RBMsHD514,140,200);  

plotRaman(RBMsHD514, 0.95, 514);        


GDBand514 = ClipSamples(GDBand514,1250,1665);  
GDBand514 = RemoveBaselinePolynomial(GDBand514, 0);
GDBand514 = NormalizeSample(GDBand514,1550,1600);  
% plotRaman(GDBand514, 0.0, 514);      


%% 520nm DATA 
% S2-S6
% I forgot to measure the empty and water filled references at this wl.

              
RBMsHD520 = {
                    DATA_20240320.S2H520R
                    DATA_20240320.S3H520R
                    DATA_20240320.S4H520R
                    DATA_20240320.S5H520R
                    DATA_20240320.S6H520R
                    DATA_20240320.S7H520R
            };   
        
               
GBand520HD = {
                    DATA_20240321.S2H520G
                    DATA_20240321.S3H520G
                    DATA_20240321.S4H520G
                    DATA_20240321.S5H520G
                    DATA_20240321.S6H520G
                    DATA_20240321.S7H520G
            };   
        
GBand520LD = {
                    DATA_20240321.S2L520G
                    DATA_20240321.S3L520G
                    DATA_20240321.S4L520G
                    DATA_20240321.S5L520G
                    DATA_20240321.S6L520G
                    DATA_20240321.S7L520G
            };   
        

RBMsHD520 = FlatFieldCorrection(RBMsHD520, DATA_20240320.FFH520R)
RBMsHD520 = ClipSamples(RBMsHD520,124,230);      
RBMsHD520 = RemoveBaselinePolynomial(RBMsHD520, 0);
RBMsHD520 = RemoveLinearBackground(RBMsHD520, 125, 225);
RBMsHD520 = RemoveBaselinePolynomial(RBMsHD520, 0);
RBMsHD520 = NormalizeSample(RBMsHD520,140,200);  
% plotRaman(RBMsHD520, 0.0, 520);        

GBand520HD = FlatFieldCorrection(GBand520HD, DATA_20240321.FFH520G)
GBand520HD = ClipSamples(GBand520HD,1533,1620);      
GBand520HD = RemoveBaselinePolynomial(GBand520HD, 0);
GBand520HD = RemoveLinearBackground(GBand520HD, 1534, 1619);
GBand520HD = RemoveBaselinePolynomial(GBand520HD, 0);
GBand520HD = NormalizeSample(GBand520HD,1580,1600);  
% plotRaman(GBand520HD, 0.0, 520);        

GBand520LD = ClipSamples(GBand520LD,1264,1682);      
GBand520LD = RemoveBaselinePolynomial(GBand520LD, 0);
GBand520LD = RemoveLinearBackground(GBand520LD, 1265, 1681);
GBand520LD = RemoveBaselinePolynomial(GBand520LD, 0);
GBand520LD = NormalizeSample(GBand520LD,1580,1600);  
% plotRaman(GBand520LD, 0.0, 520);        


%% 650nm DATA 
% S2-S6 + References
   
GBands650 =    {
            DATA_20240514.EAL650G
            DATA_20240514.S2L650G
            DATA_20240514.S3L650G
            DATA_20240514.S4L650G
            DATA_20240514.S5L650G
            DATA_20240514.S6L650G
            DATA_20240514.S7L650G
            DATA_20240514.WAL650G

            };
        
DBands650 =    {
            DATA_20240514.EAL650D
            DATA_20240514.S2L650D
            DATA_20240514.S3L650D
            DATA_20240514.S4L650D
            DATA_20240514.S5L650D
            DATA_20240514.S6L650D
            DATA_20240514.S7L650D
            DATA_20240514.WAL650D
            };
        
RBMs650 =    {
            DATA_20240514.EAL650R
            DATA_20240514.S2L650R
            DATA_20240514.S3L650R
            DATA_20240514.S4L650R
            DATA_20240514.S5L650R
            DATA_20240514.S6L650R
            DATA_20240514.S7L650R
            DATA_20240514.WAL650R
            };
        

RBMs650 = ClipSamples(RBMs650,100,250);      
RBMs650 = RemoveLinearBackground(RBMs650, 101, 249);
RBMs650 = RemoveBaselinePolynomial(RBMs650, 0);
RBMs650 = NormalizeSample(RBMs650,140,200);  
% plotRaman(RBMs650, 0.0, 650);        

GBands650 = ClipSamples(GBands650,1439,1655);  
GBands650 = RemoveLinearBackground(GBands650, 1440, 1654);
GBands650 = RemoveBaselinePolynomial(GBands650, 0);
GBands650 = NormalizeSample(GBands650,1450,1650);  
% plotRaman(GBands650, 0.0, 650);      

DBands650 = ClipSamples(DBands650,1233,1460);  
DBands650 = RemoveLinearBackground(DBands650, 1234, 1459);
DBands650 = RemoveBaselinePolynomial(DBands650, 0);
DBands650 = NormalizeSample(DBands650,1400,1410);  
% plotRaman(DBands650, 0.0, 650);      

%Plotting together
% plotRaman([DBands650; GBands650], 0.0, 650);      


%% 680nm DATA 
               
        
RBMs680 =    {
            DATA_20240325.REH680R
            DATA_20240325.S2H680R
            DATA_20240325.S3H680R
            DATA_20240325.S4H680R
            DATA_20240325.S5H680R
            DATA_20240325.S6H680R
            DATA_20240325.S7H680R
            DATA_20240325.RWH680R
            };
        
RBMs680 = ClipSamples(RBMs680,147,200);      
RBMs680 = RemoveLinearBackground(RBMs680, 148, 199);
RBMs680 = RemoveBaselinePolynomial(RBMs680, 0);
RBMs680 = NormalizeSample(RBMs680,150,200);  
% plotRaman(RBMs680, 0.0, 680);    

%%Raman per sample

S2 = {
        RBMsHD514{2}
        RBMsHD520{1}
        RBMs650{2}
        RBMs680{2}
        };
    
S3 = {
    RBMsHD514{3}
    RBMsHD520{2}
    RBMs650{3}
    RBMs680{3}
    };

S4 = {
    RBMsHD514{4}
    RBMsHD520{3}
    RBMs650{4}
    RBMs680{4}
    };

S5 = {
    RBMsHD514{5}
    RBMsHD520{4}
    RBMs650{5}
    RBMs680{5}
    };

S6 = {
    RBMsHD514{6}
    RBMsHD520{5}
    RBMs650{6}
    RBMs680{6}
    };

S7 = {
    RBMsHD514{7}
    RBMsHD520{6}
    RBMs650{7}
    RBMs680{7}
    };

% plotRaman(S2, 0.5)
% plotRaman(S3, 0.5)
% plotRaman(S4, 0.5)
% plotRaman(S5, 0.5)
% plotRaman(S6, 0.5)
% plotRaman(S7, 0.5)

%% FITTING Just for the G Bands

%G 514L
%G 520L
%G 520H
%G 650L

% plotRaman(GDBand514, 0)
% fitLorentzianToSamples(GDBand514)

% plotRaman(GBand520LD, 0)
G520Fit = fitLorentzianToSamplesCurvePeak(GBand520LD, 1560)
% plotRaman([GBand520LD;G520Fit],0)
% plotRaman(GBand520HD, 0)
% plotRaman(GBands650, 0)


%% Raman Mapping???

    
% Assuming each dataset has wavenumber (X) and intensity (Y) data
% Example for 514nm laser line:

wavenumbers514 = RBMsHD514{2}.X;  % Raman shift in cm^-1
intensities514 = RBMsHD514{2}.Y;  % Intensity (arbitrary units)

wavenumbers520 = RBMsHD520{1}.X;  % Raman shift for 520 nm
intensities520 = RBMsHD520{1}.Y;  % Intensity for 520 nm

wavenumbers650 = RBMs650{2}.X;  % Raman shift for 650 nm
intensities650 = RBMs650{2}.X;  % Intensity for 650 nm

wavenumbers680 = RBMs680{2}.X;  % Raman shift for 680 nm
intensities680 = RBMs680{2}.Y;  % Intensity for 680 nm

% % Normalize the intensities (optional)
intensities514 = intensities514 / max(intensities514);
intensities520 = intensities520 / max(intensities520);
intensities650 = intensities650 / max(intensities650);
intensities680 = intensities680 / max(intensities680);

% Find the minimum and maximum wavenumbers across all spectra
minWavenumber = min([min(wavenumbers514), min(wavenumbers520), min(wavenumbers650), min(wavenumbers680)]);
maxWavenumber = max([max(wavenumbers514), max(wavenumbers520), max(wavenumbers650), max(wavenumbers680)]);

% Define the common wavenumber range (you can change the number of points if needed)
commonWavenumbers = linspace(minWavenumber, maxWavenumber, 1000);

% Interpolate each spectrum to the common wavenumber range, using zero for missing data
intensities514_interp = interp1(wavenumbers514, intensities514, commonWavenumbers, 'pchip', 0);
intensities520_interp = interp1(wavenumbers520, intensities520, commonWavenumbers, 'pchip', 0);
intensities650_interp = interp1(wavenumbers650, intensities650, commonWavenumbers, 'pchip', 0);
intensities680_interp = interp1(wavenumbers680, intensities680, commonWavenumbers, 'pchip', 0);

% Create a 2D intensity matrix, where each row corresponds to one laser line
intensityMatrix = [
    intensities514_interp;  % 514 nm
    intensities520_interp;  % 520 nm
    intensities650_interp;  % 650 nm
    intensities680_interp   % 680 nm
];
   

%% TODO: Put this functions into useful functions:


function plotRamanBreakAxis(SamplesToPlot, offset, cutfrom, cutto)
        % Create a figure for the plot
        figure;
        
        for sampleIdx = 1:length(SamplesToPlot)
            currentSample = SamplesToPlot{sampleIdx};
            
            % Get the current sample, X values, and Y values
            currentX = currentSample.X;
            currentY = currentSample.Y - offset*sampleIdx;
            currentN = currentSample.N;
            plot(currentX, currentY, 'DisplayName', currentN,'LineWidth', 1.3);
%              scatter(currentX, currentY, 1.5, 'filled', 'DisplayName', currentN);

            hold on; % Add spectra to the same plot
        end
        breakxaxis([cutfrom cutto]);
        grid on;
        % Add labels and legend
        xlabel('Raman Shift (cm^{-1})', 'FontSize', 14);
        ylabel('Normalized Intesity (a.u.)', 'FontSize', 14)
        title('Raman Spectra');
        legend('show','FontSize', 11);        % Optional: Customize the plot further if needed
        % Hold off to stop adding new plots to the current figure
        hold off;
end

function DS = clipRangeEdges(DS, min_value, max_value)
    % Find the indices where X is within the specified range
    idx = DS.X >= min_value & DS.X <= max_value;
    
    % Filter X and Y based on these indices
    DS.X = DS.X(idx);
    DS.Y = DS.Y(idx);
    DS.P = DS.P(idx);
end

function DS = remove_lorentzian_profileNew(DS)
    % Remove Lorentzian profile and constant offset from Raman spectrum
    % DS: structure with fields X (Raman shift) and Y (intensity)

    % Extract the X and Y data
    X = DS.X;  % Raman shift
    Y = DS.Y;  % Intensity values

    % Initial guesses for the parameters: [A, gamma, C]
    initial_params = [10, 1, min(Y)];  % [Amplitude, HWHM, Constant]

    % Define the Lorentzian function with a constant offset
    lorentzian_func = @(params, x) params(1) * (params(2)^2) ./ (x.^2 + params(2)^2);

    % Set bounds for parameters: [A_min, A_max; gamma_min, gamma_max; C_min, C_max]
    lb = [0, 0.01];  % Lower bounds
    ub = [Inf, Inf];  % Upper bounds

    % Define the objective function for fitting
    objective_func = @(params) Y - lorentzian_func(params, X);

    % Perform the least squares fitting using lsqcurvefit
    options = optimset('Display', 'off');  % Suppress output
    [optimal_params, ~, exitflag] = lsqcurvefit(@(params, x) lorentzian_func(params, x), ...
                                                  initial_params, ...
                                                  X, ...
                                                  Y, ...
                                                  lb, ...
                                                  ub, ...
                                                  options);

    if exitflag <= 0
        warning('Fitting did not converge. Adjust initial parameters or bounds.');
        return; % Exit the function if fitting failed
    end

    % Extract the optimized parameters
    optimal_A = optimal_params(1);
    optimal_gamma = optimal_params(2);

    % Compute the fitted Lorentzian with the constant offset using the optimized parameters
    fitted_lorentzian = lorentzian_func(optimal_params, X);

    % Calculate the corrected spectrum
    Y_corrected = Y - fitted_lorentzian;

    % Check for overcorrection and adjust parameters accordingly
    overcorrection = Y_corrected < 0;

    % Adjust the fitted parameters if necessary
    if any(overcorrection)
        % Decrease amplitude if there are negative values
        optimal_A = optimal_A * 0.8; 
        % Increase gamma to broaden the fit
        optimal_gamma = optimal_gamma * 1.1; 
        % Compute the new fitted Lorentzian with the adjusted parameters
        fitted_lorentzian = lorentzian_func([optimal_A, optimal_gamma], X);
        Y_corrected = Y - fitted_lorentzian; % Recompute the corrected spectrum
    end

    % Store the corrected Y values back into the structure
    DS.Y = Y_corrected;
    DS.X = X;

    % Optionally, display the fitted parameters
    disp(['Fitted Amplitude (A): ', num2str(optimal_A)]);
    disp(['Fitted Gamma (HWHM): ', num2str(optimal_gamma)]);
end

function DS = remove_baseline_polynomial(DS, degree)
    % Remove baseline from Raman spectrum using polynomial fitting
    % DS: structure with fields X (Raman shift) and Y (intensity)
    % degree: Degree of the polynomial used for baseline fitting

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

    % Optionally, display the polynomial coefficients
    disp(['Polynomial Coefficients: ', num2str(p)]);
end

function DS = remove_lorentzian_profile(DS, A, gamma)
    % Remove Lorentzian profile centered at zero from Raman spectrum
    % DS: structure with fields X (Raman shift) and Y (intensity)
    % gamma: HWHM (half width at half maximum) for the Lorentzian profile

    % Extract the X and Y data
    X = DS.X;  % Raman shift (assumed centered at zero)
    Y = DS.Y;  % Intensity values

    % Define the Lorentzian function centered at zero
    lorentzian = @(x, A, gamma) A * (gamma^2) ./ (x.^2 + gamma^2);

    
    fitted_lorentzian = lorentzian(X, A, gamma);

    % Subtract the fitted Lorentzian from the original spectrum
    Y_corrected = Y - fitted_lorentzian;

    % Find the minimum value of the corrected spectrum
    min_val = min(Y_corrected);

    % Shift the corrected spectrum so that its minimum value is zero
    Y_corrected = Y_corrected - min_val;

    % Update the structure with the corrected Y values
    DS.Y = Y_corrected;

    % Optionally, return the fitted Lorentzian parameters
    % if nargout > 1
    %     varargout{1} = fit_params;
    % end
    
     % Optionally, display the fitted parameters
%     disp(['Fitted Amplitude (A): ', num2str(fit_params(1))]);
end

function DS = remove_linear_background(DS, x1, x2)
    % Remove linear background from Raman spectrum
    % DS: structure with fields X (Raman shift) and Y (intensity)
    % x1, x2: X-values (Raman shift) used to estimate the initial linear background

    % Extract the X and Y data
    X = DS.X;  % Raman shift
    Y = DS.Y;  % Intensity values

    % Ensure the selected X-values are within the valid range
    if x1 < min(X) || x2 > max(X) || x1 >= x2
        error('Invalid X values selected for background estimation.');
    end

    % Find the indices of the two X values in the X array
    [~, idx1] = min(abs(X - x1));  % Find the index closest to x1
    [~, idx2] = min(abs(X - x2));  % Find the index closest to x2

    % Define the Y-values corresponding to the selected X-values
    y1 = Y(idx1);
    y2 = Y(idx2);

    % Calculate the initial slope (m) and intercept (b) of the linear background
    slope = (y2 - y1) / (x2 - x1);  % m = (y2 - y1) / (x2 - x1)
    intercept = y1 - slope * x1;    % b = y1 - m * x1

    % Define the linear baseline function
    linear_baseline = @(x) slope * x + intercept;

    % Evaluate the linear baseline across the entire spectrum
    baseline = linear_baseline(X);

    % Subtract the baseline from the original intensity
    Y_corrected = Y - baseline;

    % Find the minimum value of the corrected spectrum
    min_val = min(Y_corrected);

    % If the minimum value is negative, adjust the baseline to make the spectrum non-negative
    if min_val < 0
        % Shift the entire spectrum by the absolute value of the minimum
        Y_corrected = Y_corrected - min_val;
    end

    % Update the structure with the corrected Y values
    DS.Y = Y_corrected;

    % Optionally, display the final slope and intercept values
    disp(['Final Slope: ', num2str(slope)]);
    disp(['Final Intercept: ', num2str(intercept)]);
end

function DS = subtract_constant_baseline(DS, x1, x2)
    % Subtract constant baseline from the spectrum based on a flat region
    % DS: structure with fields X (Raman shift) and Y (intensity)
    % x1: starting X value of the flat region
    % x2: ending X value of the flat region

    % Extract X and Y data
    X = DS.X;  % Raman shift
    Y = DS.Y;  % Intensity values

    % Identify the indices of the flat region between x1 and x2
    flat_region_indices = (X >= x1) & (X <= x2);

    % Calculate the average intensity value in the flat region
    average_baseline = mean(Y(flat_region_indices));

    % Subtract the average baseline from the entire spectrum
    Y_corrected = Y - average_baseline;

    % Update the structure with the corrected Y values
    DS.Y = Y_corrected;

    % Optionally, display the subtracted baseline value
    disp(['Subtracted Baseline: ', num2str(average_baseline)]);
end


%Apply to list of samples
function fittingResults = fitLorentzianToSamplesCurvePeak(samplesToFit, peakPositionGuess)
    % Initialize output cell array to hold the fitting results
    fittingResults = cell(size(samplesToFit));
    
    % Define the Lorentzian function
    lorentzian = @(p, x) p(1) ./ (1 + ((x - p(2)) / p(3)).^2);
    % p(1) = Peak Intensity, p(2) = Peak Position, p(3) = HWHM
    
    % Iterate over each sample to apply Lorentzian fitting
    for sampleIdx = 1:length(samplesToFit)
        % Get the current sample
        currentSample = samplesToFit{sampleIdx};
        
        % Extract X (Raman shift) and Y (Intensity) data
        X = currentSample.X;  % Raman shift (independent variable)
        Y = currentSample.Y;  % Intensity values (dependent variable)
        
        % Find the closest position to the provided peakPositionGuess
        [~, peakIdx] = min(abs(X - peakPositionGuess));
        peakPosition = X(peakIdx);  % Actual position of the peak closest to the guess
        
        % Define a search window around the peak (e.g., +/- 10 cm-1)
        searchWindow = 10;  % 10 cm-1 window around the guessed peak position
        windowIdx = (X >= (peakPosition - searchWindow)) & (X <= (peakPosition + searchWindow));
        
        % Use the data within the window for fitting
        X_fit = X(windowIdx);
        Y_fit = Y(windowIdx);
        
        % Define initial guess for the parameters [Intensity, Position, HWHM]
        initialGuess = [max(Y_fit), peakPosition, 10];  % Guess: peak intensity, peak position, HWHM
        
        % Define bounds for the parameters
        lowerBounds = [0, peakPosition - searchWindow, 0];  % Lower bounds for parameters
        upperBounds = [inf, peakPosition + searchWindow, max(X) - min(X)];  % Upper bounds for parameters
        
        % Perform the fitting using lsqcurvefit (nonlinear least squares)
        options = optimset('Display', 'off');
        [fitParams, ~] = lsqcurvefit(lorentzian, initialGuess, X_fit, Y_fit, lowerBounds, upperBounds, options);
        
        % Extract the fitted parameters: Peak Intensity, Position, HWHM
        peakIntensity = fitParams(1);  % Peak intensity (I_0)
        peakPosition = fitParams(2);   % Peak position (x_0)
        HWHM = fitParams(3);           % Half-width at half maximum (gamma)
        FWHM = 2 * HWHM;               % Full width at half maximum
        
        % Create the fitted Lorentzian curve for the selected peak
        lorentzianCurveY = lorentzian(fitParams, X);  % Y values for the fitted curve
        
        % Store the fitting results and Lorentzian curve for the current sample
        fittingResults{sampleIdx}.PeakIntensity = peakIntensity;
        fittingResults{sampleIdx}.PeakPosition = peakPosition;
        fittingResults{sampleIdx}.HWHM = HWHM;
        fittingResults{sampleIdx}.FWHM = FWHM;
        fittingResults{sampleIdx}.X = X;
        fittingResults{sampleIdx}.Y = lorentzianCurveY;
        fittingResults{sampleIdx}.N = currentSample.N;
        
        % Optionally, display the fitting results
        disp(['Sample ', num2str(sampleIdx), ':']);
        disp(['  Peak Intensity: ', num2str(peakIntensity)]);
        disp(['  Peak Position: ', num2str(peakPosition)]);
        disp(['  HWHM: ', num2str(HWHM)]);
        disp(['  FWHM: ', num2str(FWHM)]);
    end
end
function fittingResults = fitLorentzianToSamplesCurve(samplesToFit)
    % Initialize output cell array to hold the fitting results
    fittingResults = cell(size(samplesToFit));
    
    % Define the Lorentzian function
    lorentzian = @(p, x) p(1) ./ (1 + ((x - p(2)) / p(3)).^2);
    % p(1) = Peak Intensity, p(2) = Peak Position, p(3) = HWHM
    
    % Iterate over each sample to apply Lorentzian fitting
    for sampleIdx = 1:length(samplesToFit)
        % Get the current sample
        currentSample = samplesToFit{sampleIdx};
        
        % Extract X (Raman shift) and Y (Intensity) data
        X = currentSample.X;  % Raman shift (independent variable)
        Y = currentSample.Y;  % Intensity values (dependent variable)
        
        % Define initial guess for the parameters [Intensity, Position, HWHM]
        initialGuess = [max(Y), X(find(Y == max(Y), 1)), 10];  % Guess: peak intensity, peak position, HWHM
        
        % Define bounds for the parameters
        lowerBounds = [0, min(X), 0];  % Lower bounds for parameters
        upperBounds = [inf, max(X), max(X)-min(X)];  % Upper bounds for parameters
        
        % Perform the fitting using lsqcurvefit (nonlinear least squares)
        options = optimset('Display', 'off');
        [fitParams, ~] = lsqcurvefit(lorentzian, initialGuess, X, Y, lowerBounds, upperBounds, options);
        
        % Extract the fitted parameters: Peak Intensity, Position, HWHM
        peakIntensity = fitParams(1);  % Peak intensity (I_0)
        peakPosition = fitParams(2);   % Peak position (x_0)
        HWHM = fitParams(3);           % Half-width at half maximum (gamma)
        FWHM = 2 * HWHM;               % Full width at half maximum
        
        % Create the fitted Lorentzian curve
        lorentzianCurveY = lorentzian(fitParams, X);  % Y values for the fitted curve
        
        % Store the fitting results and Lorentzian curve for the current sample
        fittingResults{sampleIdx}.PeakIntensity = peakIntensity;
        fittingResults{sampleIdx}.PeakPosition = peakPosition;
        fittingResults{sampleIdx}.HWHM = HWHM;
        fittingResults{sampleIdx}.FWHM = FWHM;
        fittingResults{sampleIdx}.X = X;
        fittingResults{sampleIdx}.Y = lorentzianCurveY;
        fittingResults{sampleIdx}.N = currentSample.N;
        
        % Optionally, display the fitting results
        disp(['Sample ', num2str(sampleIdx), ':']);
        disp(['  Peak Intensity: ', num2str(peakIntensity)]);
        disp(['  Peak Position: ', num2str(peakPosition)]);
        disp(['  HWHM: ', num2str(HWHM)]);
        disp(['  FWHM: ', num2str(FWHM)]);
    end
end
function fittingResults = fitLorentzianToSamples(samplesToFit)
    % Initialize output cell array to hold the fitting results
    fittingResults = cell(size(samplesToFit));
    
    % Define the Lorentzian function
    lorentzian = @(p, x) p(1) ./ (1 + ((x - p(2)) / p(3)).^2);
    % p(1) = Peak Intensity, p(2) = Peak Position, p(3) = HWHM
    
    % Iterate over each sample to apply Lorentzian fitting
    for sampleIdx = 1:length(samplesToFit)
        % Get the current sample
        currentSample = samplesToFit{sampleIdx};
        
        % Extract X (Raman shift) and Y (Intensity) data
        X = currentSample.X;  % Raman shift (independent variable)
        Y = currentSample.Y;  % Intensity values (dependent variable)
        
        % Define initial guess for the parameters [Intensity, Position, HWHM]
        initialGuess = [max(Y), X(find(Y == max(Y), 1)), 10];  % Guess: peak intensity, peak position, HWHM
        
        % Define bounds for the parameters
        lowerBounds = [0, min(X), 0];  % Lower bounds for parameters
        upperBounds = [inf, max(X), max(X)-min(X)];  % Upper bounds for parameters
        
        % Perform the fitting using lsqcurvefit (nonlinear least squares)
        options = optimset('Display', 'off');
        [fitParams, ~] = lsqcurvefit(lorentzian, initialGuess, X, Y, lowerBounds, upperBounds, options);
        
        % Extract the fitted parameters: Peak Intensity, Position, HWHM
        peakIntensity = fitParams(1);  % Peak intensity (I_0)
        peakPosition = fitParams(2);   % Peak position (x_0)
        HWHM = fitParams(3);           % Half-width at half maximum (gamma)
        FWHM = 2 * HWHM;               % Full width at half maximum
        
        % Store the fitting results for the current sample
        fittingResults{sampleIdx}.PeakIntensity = peakIntensity;
        fittingResults{sampleIdx}.PeakPosition = peakPosition;
        fittingResults{sampleIdx}.HWHM = HWHM;
        fittingResults{sampleIdx}.FWHM = FWHM;
        
        % Optionally, display the fitting results
        disp(['Sample ', num2str(sampleIdx), ':']);
        disp(['  Peak Intensity: ', num2str(peakIntensity)]);
        disp(['  Peak Position: ', num2str(peakPosition)]);
        disp(['  HWHM: ', num2str(HWHM)]);
        disp(['  FWHM: ', num2str(FWHM)]);
    end
end


function result = concatenateSpectraNew(structArray, structName)
    % Create a new structure to hold the concatenated data
    result = struct();

    % Check if the input is not empty
    if isempty(structArray)
        return; % Return empty if no input structures
    end

    % Initialize empty arrays for concatenating X and Y values
    concatenated_x = [];
    concatenated_y = [];
    
    % Loop through the input structures
    for i = 1:length(structArray)
        current_x = structArray{i}.X;  % Current spectrum X values
        current_y = structArray{i}.Y;  % Current spectrum Y values
        
        % Append current X and Y values
        for j = 1:length(current_x)
            % Check if the current X value already exists in the concatenated array
            overlapIndex = find(concatenated_x == current_x(j), 1);
            if isempty(overlapIndex)
                % No overlap, add the new X and Y values
                concatenated_x = [concatenated_x; current_x(j)];
                concatenated_y = [concatenated_y; current_y(j)];
            else
                % Overlap found: Use the Y value from the current spectrum
                concatenated_y(overlapIndex) = current_y(j);  % Update Y value with the new spectrum value
            end
        end
    end
    
    % Remove duplicate X values and keep the corresponding Y values
    [sorted_x, uniqueIndices] = unique(concatenated_x);  % Unique X values
    sorted_y = concatenated_y(uniqueIndices);  % Corresponding Y values

    % Store the sorted X and Y values in the result structure
    result.X = sorted_x;
    result.Y = sorted_y;
    result.N = structName;
end

function result = concatenateSpectra(structArray, structName)
    % Create a new structure to hold the concatenated data
    result = struct();

    % Check if the input is not empty
    if isempty(structArray)
        return; % Return empty if no input structures
    end
    % Initialize empty arrays for concatenating X and Y values
    concatenated_x = [];
    concatenated_y = [];
    
    % Loop through the input structures and concatenate the X and Y values
    for i = 1:length(structArray)
        concatenated_x = [concatenated_x; structArray{i}.X];  % Concatenate X values
        concatenated_y = [concatenated_y; structArray{i}.Y];  % Concatenate Y values
    end
    
    % Now sort the concatenated X and Y values based on the X-values
    [sorted_x, sort_idx] = sort(concatenated_x);  % Sort X values and get the sorting indices
    sorted_y = concatenated_y(sort_idx);  % Use the sorting indices to reorder the Y values
    
    % Store the sorted X and Y values in the result structure
    result.X = sorted_x;
    result.Y = sorted_y;
    result.N = structName;
end

function baselineCorrectedSamples = SubtractConstantBaseline(samplesToCorrect, x1, x2)
    % Initialize output cell array to hold the baseline-corrected samples
    baselineCorrectedSamples = cell(size(samplesToCorrect));
    
    % Iterate over each sample to apply baseline correction
    for sampleIdx = 1:length(samplesToCorrect)
        % Get the current sample
        currentSample = samplesToCorrect{sampleIdx};
        
        % Extract X and Y data
        X = currentSample.X;  % Raman shift
        Y = currentSample.Y;  % Intensity values
        
        % Identify the indices of the flat region between x1 and x2
        flat_region_indices = (X >= x1) & (X <= x2);
        
        % Calculate the average intensity value in the flat region
        average_baseline = mean(Y(flat_region_indices));
        
        % Subtract the average baseline from the entire spectrum
        Y_corrected = Y - average_baseline;
        
        % Update the structure with the corrected Y values
        currentSample.Y = Y_corrected;
        
        % Store the baseline-corrected sample in the output cell array
        baselineCorrectedSamples{sampleIdx} = currentSample;
        
        % Optionally, display the subtracted baseline value
        disp(['Sample ', num2str(sampleIdx), ' - Subtracted Baseline: ', num2str(average_baseline)]);
    end
end

function baselineCorrectedSamples = RemoveBaselinePolynomial(samplesToCorrect, degree)
    % Initialize output cell array to hold the baseline-corrected samples
    baselineCorrectedSamples = cell(size(samplesToCorrect));
    
    % Iterate over each sample to apply polynomial baseline correction
    for sampleIdx = 1:length(samplesToCorrect)
        % Get the current sample
        currentSample = samplesToCorrect{sampleIdx};
        
        % Extract the X and Y data
        X = currentSample.X;  % Raman shift
        Y = currentSample.Y;  % Intensity values

        % Identify regions to exclude based on peak detection
        % Implement peak detection or use `findpeaks`
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
        currentSample.Y = Y_corrected;

        % Store the corrected sample in the output cell array
        baselineCorrectedSamples{sampleIdx} = currentSample;
        
        % Optionally, display the polynomial coefficients for each sample
        disp(['Sample ', num2str(sampleIdx), ' - Polynomial Coefficients: ', num2str(p)]);
    end
end

function backgroundCorrectedSamples = RemoveLinearBackground(samplesToCorrect, x1, x2)
    % Initialize output cell array to hold the background-corrected samples
    backgroundCorrectedSamples = cell(size(samplesToCorrect));
    
    % Iterate over each sample to apply linear background correction
    for sampleIdx = 1:length(samplesToCorrect)
        % Get the current sample
        currentSample = samplesToCorrect{sampleIdx};
        
        % Extract the X and Y data
        X = currentSample.X;  % Raman shift
        Y = currentSample.Y;  % Intensity values

        % Ensure the selected X-values are within the valid range
        if x1 < min(X) || x2 > max(X) || x1 >= x2
            error('Invalid X values selected for background estimation.');
        end

        % Find the indices of the two X values in the X array
        [~, idx1] = min(abs(X - x1));  % Find the index closest to x1
        [~, idx2] = min(abs(X - x2));  % Find the index closest to x2

        % Define the Y-values corresponding to the selected X-values
        y1 = Y(idx1);
        y2 = Y(idx2);

        % Calculate the initial slope and intercept of the linear background
        slope = (y2 - y1) / (x2 - x1);  % Slope (m)
        intercept = y1 - slope * x1;    % Intercept (b)

        % Define the linear baseline function
        linear_baseline = @(x) slope * x + intercept;

        % Evaluate the linear baseline across the entire spectrum
        baseline = linear_baseline(X);

        % Subtract the baseline from the original intensity
        Y_corrected = Y - baseline;

        % Find the minimum value of the corrected spectrum
        min_val = min(Y_corrected);

        % If the minimum value is negative, adjust the baseline to make the spectrum non-negative
        if min_val < 0
            % Shift the entire spectrum by the absolute value of the minimum
            Y_corrected = Y_corrected - min_val;
        end

        % Update the structure with the corrected Y values
        currentSample.Y = Y_corrected;

        % Store the corrected sample in the output cell array
        backgroundCorrectedSamples{sampleIdx} = currentSample;

        % Optionally, display the final slope and intercept values for each sample
        disp(['Sample ', num2str(sampleIdx), ' - Final Slope: ', num2str(slope)]);
        disp(['Sample ', num2str(sampleIdx), ' - Final Intercept: ', num2str(intercept)]);
    end
end


function clippedSamples = ClipSamples(samplesToClip, min_value, max_value)
    % Initialize output cell array to hold the clipped samples
    clippedSamples = cell(size(samplesToClip));
    
    % Iterate over each sample to apply clipping
    for sampleIdx = 1:length(samplesToClip)
        % Get the current sample
        currentSample = samplesToClip{sampleIdx};
        
        % Find the indices where X is within the specified range
        idx = currentSample.X >= min_value & currentSample.X <= max_value;
        
        % Clip X, Y, and P fields based on these indices
        currentSample.X = currentSample.X(idx);
        currentSample.Y = currentSample.Y(idx);
        currentSample.P = currentSample.P(idx);
        
        % Store the clipped sample in the output cell array
        clippedSamples{sampleIdx} = currentSample;
    end
end