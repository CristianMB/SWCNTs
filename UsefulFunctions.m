classdef UsefulFunctions
    
   methods (Static)
   
   %% KATAURA PLOT VALUES CALCULATION

    function [nuRBM,w11, w22, w33, w44,diam,theta, type]=CalculateKataura(P)
    % Calculate RBM frequency, energy transitions (in nm scale), diameter, chiral angle and type of nanotube for chirality [n,m]. Based on Bachilo (2002) and Araujo (2010)

    n=P(1);
    m=P(2);

    %Bachilo: https://www.science.org/doi/10.1126/science.1078727
    dcc = 0.144;        %nm Carbon carbon distance
    % Surfacntant environment
    A = 223.5;          %cm-1
    B = 12.5;           %cm-1
    % Free standing environment
%     A = 248;          %cm-1
%     B = 0;           %cm-1    
    diam= dcc*sqrt(3)*(sqrt(n.^2+m.^2+n.*m))/pi;
    theta=atan(sqrt(3)*m./(m+2*n));
    nuRBM=(A./diam)+B ;
    
    %Araujo: https://www.sciencedirect.com/science/article/pii/S1386947710000445
    a=1.049;            %eV nm
    b=0.456;
    c=0.812;            %nm-1
    hc = 1240.84193;    %h*c value to convert energy to nm
    
    
    if mod(n-m,3)==0 % METALLIC TUBES using Araujo Equations
        type = 'M';

        %M11 Araujo, two branches
        e1a=((a*3./diam).*(1+(b*log10(c./(3./diam)))))-0.19*cos(3*theta)./diam^2;    %(in eV)
        e1b=((a*3./diam).*(1+(b*log10(c./(3./diam)))))+0.29*cos(3*theta)./diam^2;    %(in eV)
        w11=hc/e1a;                                                                 %(in nm)
        w22=hc/e1b;                                                                 %(in nm)

        %M22 Araujo, two branches
        e2a=((a*6./diam).*(1+(b*log10(c./(6./diam)))))-0.60*cos(3*theta)./diam^2;    %(in eV)
        e2b=((a*6./diam).*(1+(b*log10(c./(6./diam)))))+0.57*cos(3*theta)./diam^2;    %(in eV)
        w33=hc/e2a;                                                                 %(in nm)
        w44=hc/e2b;                                                                 %(in nm)

    end
       
    if mod(n-m,3)==1 % SEMICONDUCTING TUBES using Bachilo Equations and Araujo Equations
       type = 'S';

       %S11 and S22 Bachilo with (m-n)%3 = 1
        nu11=(1*10^7)./(157.5+1066.9*diam)- 710*cos(3*theta)./diam.^2;              %(in cm-1)
        nu22=(1*10^7)./(145.6+ 575.7*diam)+1375*cos(3*theta)./diam.^2;              %(in cm-1)
        w11=(1*10^7)./nu11;                                                         %(in nm)
        w22=(1*10^7)./nu22;                                                         %(in nm)

        %S33 and S44 Araujo (p=4,5) lower33 (4), upper44 (5)
        e33=((a*4/diam).*(1+(b*log10(c./(4./diam)))))-0.42*cos(3*theta)./diam^2 +(0.0596*4/diam); %(in eV)
        e44=((a*5/diam).*(1+(b*log10(c./(5./diam)))))+0.40*cos(3*theta)./diam^2 +(0.0596*5/diam); %(in eV)
        w33=hc/e33;
        w44=hc/e44;
    end
       
    if mod(n-m,3)==2 % SEMICONDUCTING TUBES using Bachilo Equations and Araujo Equations
       type = 'S';

       %S11 and S22 Bachilo with (m-n)%3 = 2
        nu11=(1*10^7)./(157.5+1066.9*diam)+ 369*cos(3*theta)./diam.^2;      %(in cm-1)
        nu22=(1*10^7)./(145.6+ 575.7*diam)-1475*cos(3*theta)./diam.^2;      %(in cm-1)
        w11=(1*10^7)./nu11;                                                %(in nm)
        w22=(1*10^7)./nu22;                                                %(in nm)

        %S33 and S44 Araujo (p=4,5) upper33 (4), lower44 (5) + term
        %added for all beyond M11
        e33=((a*4/diam).*(1+(b*log10(c./(4./diam)))))+0.42*cos(3*theta)/diam^2+(0.0596*4/diam);  %(in eV)
        e44=((a*5/diam).*(1+(b*log10(c./(5./diam)))))-0.40*cos(3*theta)/diam^2+(0.0596*5/diam);   %(in eV)
        w33=hc/e33;                                                       %(in nm)
        w44=hc/e44;                                                       %(in nm)
    end
    end

   %% GENERAL FUNCTIONS
   

    function maximumValue = ComputeMaximum(sample, lowerLimit, upperLimit)
        % Calculate maximum value between two limit points (X1, X2)
        x = sample.X;
        y = sample.Y;
        indicesInRange = find(x >= lowerLimit & x <= upperLimit);
        maximumValue = max(y(indicesInRange));
    end
    function [peakPosition, peakValue] = ComputePeak(sample, lowerLimit, upperLimit)
        % Calculate peak position and intensity between two limit points (X1, X2)

        x = sample.X;
        y = sample.Y;
        indicesInRange = find(x >= lowerLimit & x <= upperLimit);
        yInRange = y(indicesInRange);
        [peakValue, maxIndex] = max(yInRange);
        peakPosition = x(indicesInRange(maxIndex));
    end
    function integralValue = ComputeIntegral(sample, lowerLimit, upperLimit)
        % Calculate integral between two limit points (X1, X2)
        x = sample.X;
        y = sample.Y;
        f = @(xi) interp1(x, y, xi, 'pchip');
        integralValue = integral(f, lowerLimit, upperLimit);
    end
    function normedSamples = NormalizeSample(samplesToNormalize, lowerLimit, upperLimit)
        % Divide all spectrum by maximum value or integral value between two limit points (X1, X2)
        normedSamples = cell(size(samplesToNormalize));
        % Iterate over each sample to be normalized
        for sampleIdx = 1:length(samplesToNormalize)
            currentSample = samplesToNormalize{sampleIdx};
            currentSample.Y = currentSample.Y/UsefulFunctions.ComputeIntegral(currentSample,lowerLimit, upperLimit);            
%             currentSample.Y = currentSample.Y/UsefulFunctions.ComputeMaximum(currentSample,lowerLimit, upperLimit);
            normedSamples{sampleIdx} = currentSample;
        end
    end
    function correctedSample = remove_baseline_polynomial(sample, degree)
        % Remove baseline using polynomial fitting of a chosen degree

        correctedSample = sample;
        % Extract the X and Y data
        X = correctedSample.X;  % Raman shift (assumed centered at zero)
        Y = correctedSample.Y;  % Intensity values

        % Identify regions to exclude based on peak detection
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
        correctedSample.Y = Y_corrected;

        % Optionally, display the polynomial coefficients
        disp(['Polynomial Coefficients: ', num2str(p)]);
    end
    function correctedSamples = SubstractLinearBG(samplesToCorrect, X1, X2)
     % Take a list of samples from wich you want to substract a linear profile based just on two values (X1,X2)

        CorrectedSpectra = cell(size(samplesToCorrect));
        for i = 1:length(samplesToCorrect)
            currentSample = samplesToCorrect{i};
            X = currentSample.X;
            Y = currentSample.Y;
            % Find the indices of X1 and X2 in the X array
            [~, idx_range1] = min(abs(X - X1));
            [~, idx_range2] = min(abs(X - X2));

            % Calculate the slope m
            m = (Y(idx_range2) - Y(idx_range1)) / (X2 - X1);
            % Calculate the y-intercept b
            b = Y(idx_range1) - m*X1;

            % Calculate the linear background using the equation of the line
            background = m*X + b;
                    
            % Subtract the background
            currentSample.Y = currentSample.Y - background;
            correctedSamples{i} = currentSample;
        end
    end
    function correctedSamples = SubtractInverseBG(samplesToCorrect, zeroPoints)
        % Take a list of samples from wich you want to substract backgroundA/lambdda specifiyng non-signal points to sample the background
 
        CorrectedSpectra = cell(size(samplesToCorrect));

        % Find the minimum and maximum values of the zeroPoints array
        X1 = min(zeroPoints);
        X2 = max(zeroPoints);

        for i = 1:length(samplesToCorrect)
            currentSample = samplesToCorrect{i};
            X = currentSample.X;
            Y = currentSample.Y;

            % Find the indices of X1 and X2 in the X array
            [~, idx_range1] = min(abs(X - X1));
            [~, idx_range2] = min(abs(X - X2));

            % Find the indices of the zeroPoints in the X array
            idx_zeroPoints = arrayfun(@(p) find(abs(X - p) == min(abs(X - p)), 1), zeroPoints);

            % Extract the X and Y values at the zero points
            X_zero = X(idx_zeroPoints);
            Y_zero = Y(idx_zeroPoints);

            % Fit the background using the zero points
            A = sum(Y_zero .* (1 ./ X_zero)) / sum(1 ./ X_zero.^2);

            % Calculate the background using the fitted A
            background = A ./ X;

            % Subtract the background from the entire spectrum
            currentSample.Y = Y - background;

            % Ensure no negative values in the specified range
            while any(currentSample.Y(idx_range1:idx_range2) < 0)
                A = A * 0.99;  % Reduce A slightly if there are negative values
                background = A ./ X;
                currentSample.Y = Y - background;
            end

            correctedSamples{i} = currentSample;
        end
    end
    
    %% RAMAN MEASUREMENTS FUNCTIONS
    
    function dataStructure = ReadIndividualRamanFromFile(path)
    % Read an individual Raman spectrum file and return a structure containing all spectra takes individually for later processing.

    % Read the raw spectrum data from the file
    raw_spectrum = RdExp(path);

    % Prepare the structure to hold the raw data
    dataStructure = struct();

    % Extract the file name (assuming you want to use the file name as a field)
    [~, fileName, ~] = fileparts(path);
    fileName = upper(fileName); % Make sample name uppercase
    DataName = ['DATA_', fileName];

    % Extract X (Raman shift or wavelength) and Y (intensity) values
    X = raw_spectrum(:, 1);  % Raman shift or wavelength (X-axis)
    Y = raw_spectrum(:, 3:end); % Extract only the actual spectra (columns 3 to end)

    % Loop through each spectrum and store it individually in the structure
    for i = 1:size(Y, 2)  % Loop through the number of spectra
        spectrumName = sprintf('S%d', i);
        dataStructure.(spectrumName).X = X;             % Store the X values
        dataStructure.(spectrumName).Y = Y(:, i);       % Store the Y values for each spectrum
        dataStructure.(spectrumName).N = sprintf('S%d', i);
    end
    
    % Optionally assign the data structure to the caller's workspace
    assignin('caller', DataName, dataStructure)
    end            
    function dataStructures = ReadRamanFromPaths(paths, NDel)
    % Read an individual Raman spectrum file, perform multimedian filtering and return a structure containing X (Raman shift),P (Pixel positions),Y values(Intensity)

        dataStructures = struct();  
        for p = 1:length(paths)
            try
                dataset_name = strsplit(paths{p}, "\");
                fieldName = ['DATA_', strrep(dataset_name{end-1}, '.', '')];
                dirInfo = dir(fullfile(paths{p}, '*.m3d'));
                fileList = {dirInfo(~[dirInfo.isdir]).name};
                structure = struct();

                for f = 1:length(fileList)
                    raw_spectrum = RdExp([paths{p},fileList{f}]);
                    raw_spectrum(:,2)=[];
                    NumSpec=length(raw_spectrum(1,:))-1;
                    NumDel=NDel;
                    X=raw_spectrum(:,1);
                    for i=1:1024
                        spectrum=sort(raw_spectrum(i,2:NumSpec+1));
                        Y(:,i)= mean(spectrum(NumDel+1:NumSpec-NumDel));  
                    end
                    sampleName = upper(strrep(fileList{f}, '.m3d',''));

                    structure.(sampleName).X = X;
                    structure.(sampleName).Y = Y';
                    structure.(sampleName).N = sampleName;
                    structure.(sampleName).P = (1:1024)';
                    
                end 
               dataStructures.(fieldName) = structure;
               assignin('caller', fieldName, structure); % Assign data to a variable in the caller workspace

            catch ME
                % Print the path that caused the error
                disp(['Error reading data from path: ' paths{p}]);
                % Re-throw the error
                rethrow(ME);
            end
        end       
    end     
    function correctedSample = FlatFieldCorrectionSingle(sampleToCorrect, FlatField)
    % Takes a single sample to correct by a flat field file.
    NormedFlat = FlatField.Y / UsefulFunctions.ComputeIntegral(FlatField, FlatField.X(end), FlatField.X(1));
    
    % Perform the flat-field correction on the single sample
    correctedSample = sampleToCorrect;  % Initialize the corrected sample
    correctedSample.Y = sampleToCorrect.Y ./ NormedFlat;  % Correct the Y values
    end
    function correctedSamples = FlatFieldCorrection(samplesToCorrect, FlatField)
     % Takes a list of samples to correct by a flat field file and then divide all spectra by the file after normalization of the flat field
        correctedSamples = cell(size(samplesToCorrect));
        NormedFlat = FlatField.Y / UsefulFunctions.ComputeIntegral(FlatField,FlatField.X(end), FlatField.X(1));
        
        for i = 1:length(samplesToCorrect)
            currentSample = samplesToCorrect{i};
            currentSample.Y = currentSample.Y ./ NormedFlat;
            correctedSamples{i} = currentSample;
        end
    end   

    function plotRaman(SamplesToPlot, offset, wl)
        % Create a figure for the plot
        figure;

        % Get a ColorBrewer colormap (e.g., 'Set1', 'Dark2', etc.)
        numSamples = length(SamplesToPlot);  % Number of samples to plot
        cmap = brewermap(numSamples + 1, 'Set1');  % Generate 1 more color than needed to skip the 6th
        cmap(6, :) = [];  % Remove the 6th color from the colormap
        for sampleIdx = 1:numSamples
            currentSample = SamplesToPlot{sampleIdx};

            % Get the current sample, X values, and Y values
            currentX = currentSample.X;
            currentY = currentSample.Y - offset * sampleIdx;
            currentN = currentSample.N;

            % Plot each sample using a different color from the colormap
            plot(currentX, currentY, 'Color', cmap(sampleIdx, :), 'DisplayName', currentN, 'LineWidth', 1.3);
            hold on;  % Add spectra to the same plot
        end

        % Add labels and legend
        xlabel('Raman Shift (cm^{-1})', 'FontSize', 14);
        ylabel('Normalized Intensity (a.u.)', 'FontSize', 14);
        % Conditional title based on wavelength parameter 'wl'
        if nargin < 3 || isempty(wl)
            title('Raman Spectra');
        else
            title(['Raman spectra at ', num2str(wl), ' nm']);
        end
            % Show legend with proper font size
        legend('show', 'FontSize', 11);

        % Optional: Customize the plot further if needed
        grid on;

        % Hold off to stop adding new plots to the current figure
        hold off;
    end
    function plotRamanGroup(SamplesToPlot, offset, groupingIndex, wl)
        % Create a figure for the plot
        figure;

        % Calculate the total number of individual spectra
        totalSpectra = length(SamplesToPlot);

        % Get a ColorBrewer colormap for each individual spectrum
        cmap = brewermap(10 + 1, 'Set1');  % Generate 1 more color than needed to skip the 6th
        cmap(6, :) = [];  % Remove the 6th color from the colormap
        
        % Loop through all samples and spectra
        for spectrumIdx = 1:totalSpectra
            currentSample = SamplesToPlot{spectrumIdx};

            % Get the current sample's X and Y values
            currentX = currentSample.X;
            currentY = currentSample.Y - offset * floor((spectrumIdx - 1) / groupingIndex); % Offset based on group
            currentN = currentSample.N;

            % Plot each spectrum using a unique color from the colormap
            plot(currentX, currentY, 'Color', cmap(spectrumIdx, :), 'DisplayName', currentN, 'LineWidth', 1.3);
            hold on;  % Add spectra to the same plot
        end

        % Add labels and legend
        xlabel('Raman Shift (cm^{-1})', 'FontSize', 14);
        ylabel('Normalized Intensity (a.u.)', 'FontSize', 14);

        % Conditional title based on wavelength parameter 'wl'
        if nargin < 4 || isempty(wl)
            title('Raman Spectra');
        else
            title(['Raman spectra at ', num2str(wl), ' nm']);
        end

        % Show legend with proper font size
        legend('show', 'FontSize', 11);

        % Optional: Customize the plot further if needed
        grid on;

        % Hold off to stop adding new plots to the current figure
        hold off;
    end

    function sampleList = GDBandPeaksCalculation(sampleList, LGp, HGp, LGm, HGm, LD, HD)

        % Iterate over each sample to be normalized
        for sampleIdx = 1:length(sampleList)
            currentSample = sampleList{sampleIdx};
            [GpX, GpY] = UsefulFunctions.ComputePeak(currentSample,LGp, HGp);
            currentSample.GpX = GpX;
            currentSample.GpY = GpY;

            [GmX, GmY] = UsefulFunctions.ComputePeak(currentSample,LGm, HGm);
            currentSample.GmX = GmX;
            currentSample.GmY = GmY;

            [DX, DY] = UsefulFunctions.ComputePeak(currentSample,LD, HD);
            currentSample.DX = DX;
            currentSample.DY = DY;

            sampleList{sampleIdx} = currentSample;
        end
    end
    function exportGDBandPeaks(sampleList, fileName)
        % Define the peak fields to include in the CSV file
        peakFields = {'N', 'GpX', 'GpY', 'GmX', 'GmY', 'DX', 'DY'};

        % Create a cell array to store the data
        data = cell(length(sampleList), length(peakFields));

        % Fill in the data for each sample
        for sampleIdx = 1:length(sampleList)
            currentSample = sampleList{sampleIdx};

            % Fill in the peak information for the current sample
            for peakIdx = 1:length(peakFields)
                peakValue = currentSample.(peakFields{peakIdx});
                data{sampleIdx, peakIdx} = peakValue;
            end
        end

        % Create column names
        columnNames = peakFields;

        % Create a table from the data
        dataTable = cell2table(data, 'VariableNames', columnNames);

        % Write the table to a CSV file
        writetable(dataTable, fileName);
    end
    
    %% Merge and Correct Functions (Adapted from Dmitry)
    function FixedSpectra = ClipSamples(samplesToClip, ClipLeft0, ClipRight0)
        FixedSpectra = cell(size(samplesToClip));

        for i=1:length(samplesToClip)
            currentSample = samplesToClip{i};
            currentSample= UsefulFunctions.clip_spectrum(currentSample, ClipLeft0, ClipRight0);
            FixedSpectra{i} = currentSample;
        end 
    end
    function FixedSpectra = RemoveBackground(samplesToClip)
        FixedSpectra = cell(size(samplesToClip));
        for i=1:length(samplesToClip)
            currentSample = samplesToClip{i};
            currentSample= UsefulFunctions.remove_bg_poly(currentSample);
            FixedSpectra{i} = currentSample;
        end 
    end
    function FixedSpectra = RemoveInclination(samplesToClip, WL)
        FixedSpectra = cell(size(samplesToClip));
        for i=1:length(samplesToClip)
            currentSample = samplesToClip{i};
            currentSample= UsefulFunctions.remove_inclination(currentSample, WL);
            FixedSpectra{i} = currentSample;
        end 
    end
    function FixedSpectra = InstrumentCorrection(samplesToClip, WL)
        FixedSpectra = cell(size(samplesToClip));
        for i=1:length(samplesToClip)
            currentSample = samplesToClip{i};
            currentSample= UsefulFunctions.correct_instrument_response(currentSample, WL);
            FixedSpectra{i} = currentSample;
        end 
    end

    function DS = clip_spectrum(DS,ClipLeft0,ClipRight0)
    %function clips spectrum at the left and right sides by calculating its derivative; 
    %ClipLeft and ClipRight determine number of pixels to remove from the center of left and right derivative peaks,respectively; 
    %By default, ClipLeft = 8, ClipRight = 9;

    load InclinationCoeff.mat IncPoly3;
    LeftLim = 130; %all peaks of the spectrum's derivative curve should be to the left of LeftLim - determined visually from graph with all derivative curves
    RightLim = 900; %all peaks are to the right of RightLim

    ClipRangeL = ClipLeft0; %number of pixels to add/substract to the deriv. peak
    ClipRangeR = ClipRight0;

    arb = 3; %tolerance = PeakHeight/AverageValue; determines how high should be the derivative peak, in respect to average value, to be treated as real; 
    ClipDef = 150; %default number of pixels to clip if no clear peak present in the derivative

    f.x = DS.X;
    f.y = DS.Y;
    f.p = DS.P;

    deriv = diff(f.y); %finding derivative for current spectrum
    deriv = [deriv;nan];

    cond = (f.p >= LeftLim)&(f.p<=RightLim);
    mdn = median(deriv(cond), 'omitnan'); %calculating median of derivative

    condL = (f.p <= LeftLim);
    dleft = max(deriv(condL)) - mdn; %calculating difference between peak and median

    condR = (f.p>=RightLim);
    dright = mdn - min(deriv(condR));
    dev = mad(deriv,1); %calculating average deviation from the median

    rl = dleft/dev; %calculating how different (peak-median) from average deviation from median 
    rr = dright/dev;

    if rl>=arb
       pleft=f.p(deriv == max(deriv(condL)));
       if length(pleft)>=2 %to make sure there is no  such value in the right part condR (happened once)
           pleft = min(pleft);
       end
       ClipLeft = pleft + ClipRangeL;        
    else
       ClipLeft = 1+ClipDef;
    end    

    if rr>=arb
       pright=f.p(deriv == min(deriv(condR)));
       if length(pright)>=2
           pright = max(pright);
       end
       ClipRight = pright - ClipRangeR;
    else
       ClipRight = 1024-ClipDef;
    end

    %here we determine the convergence range (in terms of pixels) of the cubic polynomial IncPoly3 interpolated for the x_mean of the spectrum;
    PixFirst = interp1(IncPoly3(:,9),IncPoly3(:,7),f.x(f.p==512));
    PixLast = interp1(IncPoly3(:,9),IncPoly3(:,6),f.x(f.p==512));

    if ClipLeft<PixFirst %choose the highest starting pixel for clipping
       ClipLeft = PixFirst;
    end 
    if ClipRight >PixLast %choose the smallest ending pixel for clipping
       ClipRight = PixLast;
     end

    condition = (f.p >= ClipLeft)&(f.p <= ClipRight);
    DS.X = f.x(condition);
    DS.Y = f.y(condition);
    DS.P = f.p(condition);
    end
    function DS = remove_bg_poly(DS)
        %removes background according to the polynomial method (Zhao et al. Applied Spectroscopy 61 11 2007 - 10.1366/000370207782597003)
        X0 = DS.X;
        Y0 = DS.Y;
        spectrum0 = [X0,Y0];
        EndCycle = false; %declaring logical variable (= criterion for ending the cycle) to start the first iteration of the while loop
        i=1; %iteration number

        fitfunc = 'poly1'; %order of the polynomial function
        options = fitoptions; % initializing "fitoptions" object for removing points during the fit after the first iteration
        % options = fitoptions(fitfunc, 'Normalize','On','Robust', 'LAR'); %to activate if normalization is needed; however, it works worse than without normalization

        while not(EndCycle)

           %fitting spectrum with high-order polynomial function
           if i == 1 %choose which spectrum to fit for current iteration
               current_spectrum = spectrum0(:,2); %take intensity array of input spectrum
           else
               current_spectrum = spectrum{i-1};
           end

           [cpoly{i},gofs{i}] = fit(spectrum0(:,1), current_spectrum, fitfunc, options);

        %    if options.Normalize == 'on' %code to use when normalization is on for fit function
        %        meanx = mean(spectrum0(:,1));
        %        stdx = std(spectrum0(:,1));
        %        tempx = (spectrum0(:,1)-meanx)./stdx;
        %    else
        %        tempx = spectrum0(:,1);
        %    end
           polyfit{i} = polyval(coeffvalues(cpoly{i}),spectrum0(:,1)); %estimating polynomial function at every x value

           %calculating residual
           res{i} = current_spectrum - polyfit{i}; 

           %calculating standard deviation of residual
           dev{i} = std(res{i},1);
           sum = polyfit{i}+dev{i};

           %reconstructing model input for the next iteration
           spectrum{i} = sum;
           if i == 1    %remove peaks from fitting for first iteration
              PointsToExclude = current_spectrum > sum;
              options.Exclude = excludedata(spectrum0(:,1),current_spectrum,'indices',PointsToExclude);
           end

           cond_mod = current_spectrum < sum; %determine at which frequencies spectrum is less then polynomial function
           spectrum{i}(cond_mod) = current_spectrum(cond_mod);

           if i~= 1
               EndCycle = abs((dev{i}-dev{i-1})/dev{i})<0.05; %calculate condition to end cycle
           end
           i = i+1; %change to next iteration
        end

        X = spectrum0(:,1);
        Y = spectrum0(:,2) - polyfit{i-1}; % spectral intensity without background contribution
        DS.Y = Y;
        DS.X = X;
    end
    function DS = correct_instrument_response(DS, WL)
        load DilorXY_Instrument_Response.mat;
        XL = 10^7./(10^7./WL - DS.X); % x values in nm 
        DS.Y = DS.Y./interp1(instrument_response(:,1),instrument_response(:,2),XL);
    end    
    function DS = remove_inclination(DS,WL)
        %Function removes inclination effect for the spectrum for Dilor XY spectrometer
        
        load InclinationCoeff.mat IncPoly3;
        X0 = DS.X;
        Y0 = DS.Y;
        P0 = DS.P;
        
        XL = 10^7./(10^7./WL - X0); % change x from cm-1 to nm

        center = XL(P0==512);

        p(1) = interp1(IncPoly3(:,8), IncPoly3(:,2),center,'spline');
        p(2) = interp1(IncPoly3(:,8), IncPoly3(:,3),center,'spline');
        p(3) = interp1(IncPoly3(:,8), IncPoly3(:,4),center,'spline');
        p(4) = interp1(IncPoly3(:,8), IncPoly3(:,5),center,'spline');

        y= @(p,x) p(1)*x.^3+p(2)*x.^2+p(3)*x+p(4);

        Y = Y0./y(p,XL);
        DS.Y = Y;
    end
    
   %% ABSORPTION MEASUREMENTS FUNCTIONS
    
    function samples = readSamplesData(filePath)
        % Read the header
        header = readcell(filePath, 'Range', 'A1:AZ1');
        header = cellfun(@(x) strrep(x, ' 100%T', ''), header, 'UniformOutput', false);
        sampleNames = header(1, 1:2:end);
        warning('off', 'MATLAB:strrep:InvalidInputType');
        % Read ONLY the datalines
        data = readmatrix(filePath, 'Range', ['A' num2str(3) ':AZ' num2str(2328)]);

        samples = struct();

        % Iterate through each sample and store wavelength and absorption data
        for i = 1:length(sampleNames)
            % Extract wavelength and absorption data for the current sample
            wavelengths = data(:, 2*i - 1); % Odd columns contain wavelength
            absorption = data(:, 2*i);       % Even columns contain absorption
            warning('off', 'MATLAB:strrep:InvalidInputType');
            % Store the data in the container object with the sample name
            sampleName = strrep(sampleNames{i}, '-', '_');
            samples.(sampleName).X = wavelengths;
            samples.(sampleName).Y = absorption;
            samples.(sampleName).N = sampleName;
        end
    end
    function mergedStruct = mergeStructures(struct1, struct2)
    % Get fieldnames of both structures
    fields2 = fieldnames(struct2);
    % Merge the fields of struct2 into struct1
    for i = 1:length(fields2)
        if ~isfield(struct1, fields2{i})
            struct1.(fields2{i}) = struct2.(fields2{i});
        end
    end
        mergedStruct = struct1;
    end    
    function mergedStruct = mergeStructuresRaman(structs)
            mergedStruct = struct('X', [], 'Y', [], 'N', 'Merged');
        for i = structs
            mergedStruct.X = [mergedStruct.X;i.X];
            mergedStruct.Y = [mergedStruct.Y;i.Y];
        end
    end
    function dataStructures = ReadAbsorptionFromPaths(paths)
        dataStructures = struct();
        for i = 1:length(paths)
            % Extract the suffix from the path
            try
                % Extract the suffix from the path variable
                dataset_name = strsplit(paths{i}, "\");
                % Create dynamic field name for the structure
                fieldName = ['DATA_', strrep(dataset_name{4}, '.csv', '')];
                % Read the data from the current path
                data = UsefulFunctions.readSamplesData(paths{i});
                
                structure = struct(fieldName, data);
                
                % Check if the field already exists
                if isfield(dataStructures, fieldName)
                    dataStructures.(fieldName) = UsefulFunctions.mergeStructures(dataStructures.(fieldName), data);
                else
                    dataStructures.(fieldName) = structure.(fieldName);        
                end
                
                % Assign data to a variable in the caller workspace
                assignin('caller', fieldName, dataStructures.(fieldName)); 
                
            catch ME
                % Print the path that caused the error
                disp(['Error reading data from path: ' paths{i}]);
                % Re-throw the error
                rethrow(ME);
            end
        end
    end
    function plotAbsorption(SamplesToPlot, offset)
        % Create a figure for the plot
        figure;
        % Iterate over each sample
        cmap = brewermap(length(SamplesToPlot) + 1, 'Set1');  % Generate 1 more color than needed to skip the 6th
        cmap(6, :) = [];  % Remove the 6th color from the colormap
        
        for sampleIdx = 1:length(SamplesToPlot)
            currentSample = SamplesToPlot{sampleIdx};
                % Get the current sample, X values, and Y values
                currentX = currentSample.X;
                currentY = currentSample.Y - offset*sampleIdx;
                currentN = currentSample.N;
                plot(currentX, currentY,  'Color', cmap(sampleIdx, :), 'DisplayName', currentN,'LineWidth', 1.3);
                hold on; % Add spectra to the same plot
        end

        % Add labels and legend
        xlabel('Wavelenght (nm)', 'FontSize', 14);
        ylabel('Normalized Absorption (a.u.)', 'FontSize', 14);
        title('Absorption Spectra','FontSize', 14);
        legend('show','FontSize', 11);
        % Optional: Customize the plot further if needed
        grid on;
        % Hold off to stop adding new plots to the current figure
        hold off;

    end
    function sampleList = TransitionPeaksCalculation(sampleList, LS1, US1, LS2, US2)

        % Iterate over each sample to be normalized
        for sampleIdx = 1:length(sampleList)
            currentSample = sampleList{sampleIdx};
            [S1X, S1Y] = UsefulFunctions.ComputePeak(currentSample,LS1, US1);
            currentSample.S11X = S1X;
            currentSample.S11Y = S1Y;

            [S2X, S2Y] = UsefulFunctions.ComputePeak(currentSample,LS2, US2);
            currentSample.S2X = S2X;
            currentSample.S2Y = S2Y;

            sampleList{sampleIdx} = currentSample;
        end
    end
    function exportTransitionPeaks(sampleList, fileName)
        % Define the peak fields to include in the CSV file
        peakFields = {'N', 'S11W', 'S11A', 'S22W', 'S22A'};

        % Create a cell array to store the data
        data = cell(length(sampleList), length(peakFields));

        % Fill in the data for each sample
        for sampleIdx = 1:length(sampleList)
            currentSample = sampleList{sampleIdx};

            % Fill in the peak information for the current sample
            for peakIdx = 1:length(peakFields)
                peakValue = currentSample.(peakFields{peakIdx});
                data{sampleIdx, peakIdx} = peakValue;
            end
        end

        % Create column names
        columnNames = peakFields;

        % Create a table from the data
        dataTable = cell2table(data, 'VariableNames', columnNames);

        % Write the table to a CSV file
        writetable(dataTable, fileName);
    end
    
    
       %% FITTING FUNCTIONS
   
    function result = Lorentzian(params, x)
        %We expect an array of params. 3 params per sample
        numPeaks = numel(params) / 3;
        result = zeros(size(x));
        for i = 1:numPeaks
            amp = params(3*i - 2);
            pos = params(3*i - 1);
            width = params(3*i);
            result = result + amp ./ ((x - pos).^2 + width);
        end
    end
    function [fitParams, fitCurve] = fitMultiLorentzian(DS, peakPositions)
        % Define the multi-Lorentzian function
        X = DS.X;
        Y = DS.Y;
        multiLorentzian = @(params, x) UsefulFunctions.Lorentzian(params, x);
        % Define initial parameter guesses
        initialParams = zeros(1, 3 * length(peakPositions));
        initialParams(1:3:end) = 1;                         % Becase I normalized
        initialParams(2:3:end) = peakPositions;             % peakPositions is the Input
        initialParams(3:3:end) = 1;                        
        % Fit the multi-Lorentzian function to the data
        fitParams = lsqcurvefit(multiLorentzian, initialParams, X, Y);
        fitCurve = multiLorentzian(fitParams, X);
    end   
    function PeakList = PeakID(SampleList, peakThreshold)
        PeakList = {};
        for i = 1:length(SampleList)
            currentSample = SampleList{i};
            X = currentSample.X;
            Y = currentSample.Y;

            %sorting, just because fitting requires it.
            [X, sort_idx] = sort(X);
            Y = Y(sort_idx);

            % Find peaks using the findpeaks function
            [h, x, w, p] = findpeaks(Y, X, 'MinPeakProminence', peakThreshold);
            PeakList{i} = x;
        end
    end
    function FittedSamples = FitSamples(SampleList, peakpos)
        FittedSamples = cell(size(SampleList));
        for i = 1:length(SampleList)
            currentSample = SampleList{i};
            [fitParams, fitCurve] = UsefulFunctions.fitMultiLorentzian(currentSample, peakpos);
            currentSample.F.Fit = fitCurve;
            currentSample.F.Params = reshape(fitParams, 3, length(peakpos)); 
            FittedSamples{i} = currentSample;
         end
    end
    
    
%% FUNCTIONS UNDER DEVELOPMENT
    function plotRamanFit(Sample)
        % Create a figure for the plot
        figure;
        % Iterate over each sample
        % Get the current sample, X values, and Y values
        X = Sample.X;
        Y = Sample.Y;
        N = Sample.N;
        fitCurve = Sample.F.Fit;
        fitParams = Sample.F.Params;
        numPeaks = size(fitParams, 2);

        plot(X, Y, 'DisplayName', N,'LineWidth', 1.3);
        hold on; % Add spectra to the same plot
        plot(X, fitCurve, 'k', 'LineWidth', 1.5, 'DisplayName', 'MultiLorentzFit');
        %Plot each Lorentzian peak individually
        
        for i = 1:numPeaks
             amp = fitParams(1,i);
             pos = fitParams(2,i);  
             width = fitParams(3,i);
             peakCurve = amp ./ ((X - pos).^2 + width);
             plot(X, peakCurve, 'r--', 'LineWidth', 1, 'DisplayName', sprintf('Peak at %.1f', pos));
        end

        xlabel('Raman Shift (cm^{-1})');
        ylabel('Intesity (a.u.)');
        title('Raman Spectra');
        legend('show');
        % Optional: Customize the plot further if needed
        grid on;
        % Hold off to stop adding new plots to the current figure
        hold on;
    end
    function plotRamanPxl(SamplesToPlot, offset, wl)
        % Create a figure for the plot
        figure;

        % Get a ColorBrewer colormap (e.g., 'Set1', 'Dark2', etc.)
        numSamples = length(SamplesToPlot);  % Number of samples to plot
        cmap = brewermap(10 + 1, 'Set1');  % Generate 1 more color than needed to skip the 6th
        cmap(6, :) = [];  % Remove the 6th color from the colormap
        for sampleIdx = 1:numSamples
            currentSample = SamplesToPlot{sampleIdx};

            % Get the current sample, X values, and Y values
            currentX = currentSample.P;
            currentY = currentSample.Y - offset * sampleIdx;
            currentN = currentSample.N;

            % Plot each sample using a different color from the colormap
            plot(currentX, currentY, 'Color', cmap(sampleIdx, :), 'DisplayName', currentN, 'LineWidth', 1.3);
            hold on;  % Add spectra to the same plot
        end

        % Add labels and legend
        xlabel('Raman Shift (cm^{-1})', 'FontSize', 14);
        ylabel('Normalized Intensity (a.u.)', 'FontSize', 14);
        % Conditional title based on wavelength parameter 'wl'
        if nargin < 3 || isempty(wl)
            title('Raman Spectra');
        else
            title(['Raman spectra at ', num2str(wl), ' nm']);
        end
            % Show legend with proper font size
        legend('show', 'FontSize', 11);

        % Optional: Customize the plot further if needed
        grid on;

        % Hold off to stop adding new plots to the current figure
        hold off;
    end
    function plotRamanNewOfsetBetweenPlots(SamplesToPlot, offset)
        % Create a figure for the plot
        figure;
        % Iterate over each sample
        adjustedSampleIdx = 0;
        
        for sampleIdx = 1:length(SamplesToPlot)
            currentSample = SamplesToPlot{sampleIdx};
            
            % Check if the current sample is the 'OffSet' keyword
            if strcmp(currentSample, 'OffSet')
                % If it's 'OffSet', increment the adjustedSampleIdx to add an offset
                adjustedSampleIdx = adjustedSampleIdx + 1; 
                continue; % Skip plotting for this 'OffSet'
            end
            adjustedSampleIdx = adjustedSampleIdx + 1;
            
            % Get the current sample, X values, and Y values
            currentX = currentSample.X;
            currentY = currentSample.Y - offset*adjustedSampleIdx;
            currentN = currentSample.N;
            plot(currentX, currentY, 'DisplayName', currentN,'LineWidth', 1.3);
            hold on; % Add spectra to the same plot
        end
        
        % Add labels and legend
        xlabel('Raman Shift (cm^{-1})', 'FontSize', 14);
        ylabel('Normalized Intesity (a.u.)', 'FontSize', 14)
        title('Raman Spectra');
        legend('show','FontSize', 11);        % Optional: Customize the plot further if needed
        grid on;
        % Hold off to stop adding new plots to the current figure
        hold off;
    end

   end
  
end