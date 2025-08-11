clc;
clear;
import UsefulFunctions.*;
rootpath = 'X:\Measurements Data\Raman\';
addpath('X:\SWCNTs\RamanAnalysis\Raman - Voigt Fitting\Faddeeva_voigt');
addpath('X:\SWCNTs');

% rootpath = 'X:\Measurements Data\Raman\';
%All paths as default

path_a= [rootpath,'20240111\'];
path_b= [rootpath,'20240517\'];
path_c= [rootpath,'20240610\'];
path_d= [rootpath,'20240612\'];
path_e= [rootpath,'20241007\'];
path_f= [rootpath,'20241008\'];
path_g= [rootpath,'20241129\'];
path_h= [rootpath,'20241212\'];
path_i= [rootpath,'20241213\'];
path_j= [rootpath,'20250131\'];
path_k= [rootpath,'20241212\'];
path_l= [rootpath,'20250411\'];

path_TTF = [rootpath,'20250520\'];
path_TCNQ = [rootpath,'20250618\'];
path_aa = [rootpath,'20250131\'];
path_bbb = [rootpath,'20250403\'];
path_TCNQ_rinsed = [rootpath,'20250718\'];
path_S21_S22_S23_G = [rootpath,'20250731\'];
path_intensity = [rootpath,'20250801\'];

%Select the paths of interest

paths = {   
            path_a
            path_b
            path_c
            path_d
            path_e
            path_f
            path_g
            path_h
            path_i
            path_j
            path_k    
            path_TTF
            path_TCNQ
            path_aa
            path_bbb
            path_TCNQ_rinsed
            path_S21_S22_S23_G
            path_intensity
            path_l
        };


ReadRamanFromPaths(paths, 2);

%% %--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DATA_20240111.S240111A.N='SF D2O@SWCNT 514nm';
DATA_20240111.S240111B.N='SF TCE@SWCNT 514nm';
DATA_20240111.S240111BB.N='SF TCE@SWCNT 514nm';
DATA_20240111.S240111C.N='SF Methanol@SWCNT 514nm';
DATA_20240111.S240111D.N='SF TCE@SWCNT 514nm';
DATA_20240111.S240111E.N='SF TTF@SWCNT 514nm';
DATA_20240111.S240111F.N='SF PCE@SWCNT 514nm';
DATA_20240111.S240111G.N='SF PCE@SWCNT 514nm';
DATA_20240111.S240111H.N='SF PCE@SWCN 514nm';
DATA_20240111.S240111I.N='SF TEMED@SWCNT 514nm';
DATA_20240517.EAL514GD.N='Empty Arc SWCNTs';
DATA_20240517.S2L514GD.N='CB PCE@SWCNT Dial. DGU C (Filled)';
DATA_20240517.S3L514GD.N='CB TCE@SWCNT Dial. DGU C (Filled)';
DATA_20240517.S4L514GD.N='CB TEMED@SWCNT Dial. DGU C (Filled)';
DATA_20240517.S5L514GD.N='CB TDAE@SWCNT Dial. DGU C (Filled)';
DATA_20240517.S6L514GD.N='CB Hexadecane@SWCNT Dial. DGU C (Filled)';
DATA_20240517.S7L514GD.N='CB Dodecane@SWCNT Dial. DGU C (Filled)';
DATA_20240517.WAL514GD.N='Water Filled Arc SWCNTs';
DATA_20240610.EAL514GD.N='Empty Arc SWCNTs';
DATA_20240610.S2L514GD.N='CB PCE@SWCNT After CSA Treatment';
DATA_20240610.S3L514GD.N='CB TCE@SWCNT After CSA Treatment';
DATA_20240610.S4L514GD.N='CB TEMED@SWCNT After CSA Treatment';
DATA_20240610.S5L514GD.N='CB TDAE@SWCNT After CSA Treatment';
DATA_20240610.S6L514GD.N='CB Hexadecane@SWCNT After CSA Treatment';
DATA_20240610.S7L514GD.N='CB Dodecane@SWCNT After CSA Treatment';
DATA_20240610.WAL514GD.N='Water Filled Arc SWCNTs';
DATA_20240612.EAL514GD.N='Empty Arc SWCNTs';
DATA_20240612.S2L514GD.N='CB PCE@SWCNT After CSA Treatment + Dialysis to D2O';
DATA_20240612.S3L514GD.N='CB TCE@SWCNT After CSA Treatment + Dialysis to D2O';
DATA_20240612.S4L514GD.N='CB TEMED@SWCNT After CSA Treatment + Dialysis to D2O';
DATA_20240612.S5L514GD.N='CB TDAE@SWCNT After CSA Treatment + Dialysis to D2O';
DATA_20240612.S6L514GD.N='CB Hexadecane@SWCNT After CSA Treatment + Dialysis to D2O';
DATA_20240612.S7L514GD.N='CB Dodecane@SWCNT After CSA Treatment + Dialysis to D2O';
DATA_20240612.WAL514GD.N='Water Filled Arc SWCNTs';
DATA_20241007.F2L514GD.N='Film PCE@SWCNT';
DATA_20241007.F3L514GD.N='Film TCE@SWCNT';
DATA_20241007.F4L514GD.N='Film TEMED@SWCNT';
DATA_20241007.F5L514GD.N='Film TDEA@SWCNT';
DATA_20241007.F6L514GD.N='Film Hexadecane@SWCNT';
DATA_20241007.F7L514GD.N='Film Dodecane@SWCNT';
DATA_20241008.F0L514C.N='Sapphire SubstrateC';
DATA_20241008.F0L514D.N='Sapphire SubstrateD';
DATA_20241008.F0L514E.N='Sapphire SubstrateE';
DATA_20241008.FFL514C.N='FlatFieldC';
DATA_20241008.FFL514D.N='FlatFieldD';
DATA_20241008.FFL514E.N='FlatFieldE';
DATA_20241213.FF6L514G.N='Film 3 (SFF6 TTF@CNTs)';
DATA_20241213.FS5L514G.N='Film 2 (S5 TDAE@CNTs)';
DATA_20241213.FS7L514G.N='Film 1 (S7 Dodecane@CNTs)';
DATA_20241213.FS8L514G.N='Film 5 (S8 PCE@CNTs)';
DATA_20241213.FS9L514G.N='Film 4 (S9 TEMED@CNTs)';
DATA_20241213.P30.N='Power test to see if sample is heated at 30mW';
DATA_20241213.P50.N='Power test to see if sample is heated at 50mW';
DATA_20241213.P70.N='Power test to see if sample is heated at 70mW';
DATA_20241213.S1BL514G.N='TEMED Filled DGU (SFF1B)';
DATA_20241213.SF6L514G.N='TTF Filled DGU (SFF6)';
DATA_20241213.SR0L514G.N='SF6 Sofie Empty';
DATA_20241213.SR1L514G.N='Water Filled DGU (KVDP2)';
DATA_20250131.KITL514G .N='KIT Converted to DWCNTs (Dispersion)';
DATA_20250131.KT3L514G .N='KIT Converted to DWCNTs (Dispersion)';
DATA_20250131.KT5L514G.N='KIT Converted to DWCNTs (Dispersion)';
DATA_20250131.KTLL514G.N='KIT Converted to DWCNTs (Dispersion)';
DATA_20250131.S10L514G.N='Dispersion S10 - TMG@P2-SWCNTs';
DATA_20250131.S11L514G.N='Dispersion S11 - TTF@P2-SWCNTs';
DATA_20250131.SR0L514G.N='Dispersion SR0 - Empty@P2-SWCNTs';
DATA_20250131.SR1L514G.N='Dispersion SR1 - D2O@P2-SWCNTs - KVD';
DATA_20250131.SR2L514G.N='Dispersion SR2 - MeOH@P2-SWCNTs';
DATA_20250131.SWFL514G.N='Dispersion SWF - D2O@P2-SWCNTs (Salome)';
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
DATA_20250131.SR1L514R.N='Dispersion SR1 - D2O@P2-SWCNTs - KVD';
DATA_20250131.SWFL514R.N='Dispersion SWF - D2O@P2-SWCNTs (Salome)';
DATA_20250131.S10L514R.N='Dispersion S10 - TMG@P2-SWCNTs';
DATA_20250131.S11L514R.N='Dispersion S11 - TTF@P2-SWCNTs';
DATA_20250131.SR0L514R.N='Dispersion SR0 - Empty@P2-SWCNTs';
DATA_20250131.SR2L514R.N='Dispersion SR2 - MeOH@P2-SWCNTs';
DATA_20240111.S240111E.N='SF TTF@SWCNT 514nm';
DATA_20240111.S240111N.N='SF TTF@SWCNT 514nm';
DATA_20250131.S11L514G.N='Dispersion S11 - TTF@P2-SWCNTs';
DATA_20250131.S11L514R.N='Dispersion S11 - TTF@P2-SWCNTs';
DATA_20250411.S11L476D.N='Dispersion S11 - TTF@P2-SWCNTs (Melt)';
DATA_20250411.S12L476D.N='Dispersion S12 - TTF@P2-SWCNTs (Reflux)';
DATA_20250411.S13L476D.N='Dispersion S13 - TTF@P2-SWCNTs (GasPhase)';
DATA_20250411.S11L476G.N='Dispersion S11 - TTF@P2-SWCNTs (Melt)';
DATA_20250411.S12L476G.N='Dispersion S12 - TTF@P2-SWCNTs (Reflux)';
DATA_20250411.S13L476G.N='Dispersion S13 - TTF@P2-SWCNTs (GasPhase)';
DATA_20250411.S11L476R.N='Dispersion S11 - TTF@P2-SWCNTs (Melt)';
DATA_20250411.S12L476R.N='Dispersion S12 - TTF@P2-SWCNTs (Reflux)';
DATA_20250411.S13L476R.N='Dispersion S13 - TTF@P2-SWCNTs (GasPhase)';


%% All G vs Ref

close all;
        
G = {
DATA_20250411.P12L476G
DATA_20250411.P13L476G
DATA_20250411.PP2L476G

DATA_20250411.S11L476G
DATA_20250411.S12L476G
DATA_20250411.S13L476G
DATA_20250411.S14L476G
DATA_20250411.S15L476G


       }; 

% G = FilterDataByXRange(G, 1260, 1660);           
% G = RemovePolyBG(G, 0);
% G = SubstractLinearBG(G, 1260, 1680);
% G = Normalize(G, 1580, 1600, 'M');
% 
% plotRaman(G, 0.0, 476);     


R = {
DATA_20250411.S11L476R
DATA_20250411.S12L476R
DATA_20250411.S13L476R
DATA_20250411.S14L476R
DATA_20250411.S15L476R
    }; 

% R= FlatFieldCorrection(R, DATA_20250411.);           
R = FilterDataByXRange(R, 130, 228);           
R = RemovePolyBG(R, 1);
R = SubstractLinearBG(R, 130, 228);
R = Normalize(R, 140, 220, 'M');
% 
plotRaman(R, 0.0, 476);     

DD= {
    DATA_20241007.F2L514DD
    DATA_20241007.F3L514DD
    DATA_20241007.F4L514DD
    DATA_20241007.F5L514DD
    DATA_20241007.F6L514DD
    DATA_20241007.F7L514DD
    
    DATA_20241213.FF6L514D
    DATA_20241213.FS5L514D
    DATA_20241213.FS7L514D
    DATA_20241213.FS8L514D
    DATA_20241213.FS9L514D
}; 

% DD = FilterDataByXRange(DD, 2500, 2835);           
% DD = RemovePolyBG(DD, 0);
% DD = SubstractLinearBG(DD, 2500, 2835);
% DD = Normalize(DD, 2500, 2700, 'M');
% 
% close all
% plotRaman(DD, 0.20, 514);     




%% DOPING VECTORS


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
% plot(FG,FD)
% scatter(FG,FD, 50, 'k','d', 'filled')
% text(FG, FD, NG, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 8);
% xlim([1590, 1593]);
% ylim([2660, 2697]);
% ylabel('2D Band (cm-1)')
% xlabel('G Band (cm-1)')
% 
% %%%% FITTING V2
% 
% % FITTEDG = FitSamples(G, 1592);
% FITTEDG = FitSamples(G, [1550,1560,1590]);
% 
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
% ylim([2660, 2697]);
% ylabel('2D Band (cm^{-1})');
% xlabel('G Band (cm^{-1})');
% legend();  % Adjust legend location as needed
% title('Doping Vectors');
% 
% 
% 
% %%%%
% 
% plotRamanFits(FITTEDG,0.3)

%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotRamanFits(SamplesToPlot, offset)
    % Create a figure for the plot
    figure;
    
    for sampleIdx = 1:length(SamplesToPlot)
        currentSample = SamplesToPlot{sampleIdx};

        % Get sample data
        X = currentSample.X;
        Y = currentSample.Y - offset * sampleIdx;
        N = currentSample.N;
        fitCurve = currentSample.F.Fit - offset * sampleIdx;
        fitParams = currentSample.F.Params;
        numPeaks = size(fitParams, 2);

        % Plot original spectrum with offset
        plot(X, Y, 'DisplayName', N, 'LineWidth', 1.3);
        hold on;

        % Plot fitted curve
        plot(X, fitCurve, 'k', 'LineWidth', 1.5, 'DisplayName', sprintf('%s - Fit', N));

        % Plot individual Lorentzian peaks with offset
        for j = 1:numPeaks
            amp = fitParams(1, j);
            pos = fitParams(2, j);
            width = fitParams(3, j);
            peakCurve = amp ./ ((X - pos).^2 + width) - offset * sampleIdx;
            plot(X, peakCurve, 'r--', 'LineWidth', 1, ...
                'HandleVisibility', 'off'); % Don't crowd the legend
        end
    end

    % Labels and final formatting
    xlabel('Raman Shift (cm^{-1})', 'FontSize', 14);
    ylabel('Intensity (a.u.)', 'FontSize', 14);
    title('Raman Spectra with Multi-Lorentzian Fits', 'FontSize', 14);
    legend('show', 'FontSize', 11);
    grid on;
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
%         disp(['Structure ', num2str(i), ' Polynomial Coefficients: ', num2str(p)]);
        
        % Save the updated structure back to the list
        DSList{i} = DS;
    end
end


