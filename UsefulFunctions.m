classdef UsefulFunctions
    
   methods (Static)
   
   %% Kataura Calculation

    function [nuRBM,wl1,wl2,wl3,wl4,diam,theta, type]=CalculateKataura(P)

    n=P(1);
    m=P(2);

    diam=0.144*sqrt(3)*(sqrt(n.^2+m.^2+n.*m))/pi;
    theta=atan(sqrt(3)*m./(m+2*n));
    nuRBM=(223.5./diam)+12.5;

    a=1.049;%eV nm
    b=0.456;
    c=0.812;%nm-1

               if mod(n-m,3)==0 % METALLIC TUBES               
    %                 freq1=(1*10^7)./(150+(370*diam))+(-2000*cos(3*theta))./diam.^2;
    %                 wl1=(1*10^7)./freq1; 
    %                 freq2=10000000./(150+370*diam)+3000*cos(3*theta)./diam.^2;
    %                 wl2=(1*10^7)./freq2;
                    energy1=((a*3./diam).*(1+(b*log10(c./(3./diam)))))-0.18*cos(3*theta)/diam^2;%(in eV)
                    wl1=1240/energy1;
                    energy2=((a*3./diam).*(1+(b*log10(c./(3./diam)))))+0.29*cos(3*theta)./diam^2;%(in eV)
                    wl2=1240/energy2;
                    energy3=((a*6./diam).*(1+(b*log10(c./(6./diam)))))-0.6*cos(3*theta)/diam^2;%(in eV)
                    wl3=1240/energy3;
                    energy4=((a*6./diam).*(1+(b*log10(c./(6./diam)))))+0.87*cos(3*theta)./diam^2;%(in eV)
                    wl4=1240/energy4;
                    type = 'M';
               end
               if mod(n-m,3)==1 % semiconducting tube 
                    freq1=10000000./(157.5+1066.9*diam)-710*cos(3*theta)./diam.^2;
                    wl1=(1*10^7)./freq1; 
                    freq2=10000000./(145.6+575.7*diam)+1375*cos(3*theta)./diam.^2;
                    wl2=(1*10^7)./freq2;
                    energy3=((a*4/diam).*(1+(b*log10(c./(4/diam)))))-0.42*cos(3*theta)./diam^2+(0.0596*4/diam);%(in eV)
                    wl3=1240/energy3;
                    energy4=((a*5/diam).*(1+(b*log10(c./(5/diam)))))+0.4*cos(3*theta)./diam^2+(0.0596*5/diam);%(in eV)
                    wl4=1240/energy4;
                    type = 'S';

               end
               if mod(n-m,3)==2 % semiconducting tube 
                    freq1=10000000./(157.5+1066.9*diam)+369*cos(3*theta)./diam.^2;
                    wl1=(1*10^7)./freq1; 
                    freq2=10000000./(145.6+575.7*diam)-1475*cos(3*theta)./diam.^2;
                    wl2=(1*10^7)./freq2;
                    energy3=((a*4/diam).*(1+(b*log10(c./(4/diam)))))+0.42*cos(3*theta)/diam^2+(0.0596*4/diam);%(in eV)
                    wl3=1240/energy3;
                    energy4=((a*5/diam).*(1+(b*log10(c./(5/diam)))))-0.4*cos(3*theta)/diam^2+(0.0596*5/diam);%(in eV)
                    wl4=1240/energy4;
                    type = 'S';
               end
    end

   %% General Functions

       
    function integralValue = ComputeIntegral(sample, lowerLimit, upperLimit)
        x = sample.X;
        y = sample.Y;
        f = @(xi) interp1(x, y, xi, 'pchip');
        integralValue = integral(f, lowerLimit, upperLimit);
    end

    function maximumValue = ComputeMaximum(sample, lowerLimit, upperLimit)
        x = sample.X;
        y = sample.Y;
        indicesInRange = find(x >= lowerLimit & x <= upperLimit);
        maximumValue = max(y(indicesInRange));
    end
    
    function NormedSamples = NormalizeSample(samplesToNormalize, lowerLimit, upperLimit)
        NormedSamples = cell(size(samplesToNormalize));
        % Iterate over each sample to be normalized
        for sampleIdx = 1:length(samplesToNormalize)
            currentSample = samplesToNormalize{sampleIdx};
            %currentSample.Y = currentSample.Y/ComputeIntegral(currentSample,lowerLimit, upperLimit);            
            currentSample.Y = currentSample.Y/UsefulFunctions.ComputeMaximum(currentSample,lowerLimit, upperLimit);
            NormedSamples{sampleIdx} = currentSample;
        end
    end

    function [peakPosition, peakValue] = ComputePeak(sample, lowerLimit, upperLimit)
        x = sample.X;
        y = sample.Y;
        indicesInRange = find(x >= lowerLimit & x <= upperLimit);
        yInRange = y(indicesInRange);
        [peakValue, maxIndex] = max(yInRange);
        peakPosition = x(indicesInRange(maxIndex));
    end
        
    %% Raman Functions

    function dataStructures = ReadRamanFromPaths(paths)
        dataStructures = struct();  
        for p = 1:length(paths)
            try
                dataset_name = strsplit(paths{p}, "\");
                fieldName = ['DATA_', strrep(dataset_name{7}, '.', '')];
                dirInfo = dir(fullfile(paths{p}, '*.m3d'));
                fileList = {dirInfo(~[dirInfo.isdir]).name};
                structure = struct();

                for f = 1:length(fileList)
                    raw_spectrum = RdExp([paths{p},fileList{f}]);
                    raw_spectrum(:,2)=[];
                    NumSpec=length(raw_spectrum(1,:))-1;
                    NumDel=1;
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
    
    function CorrectedSpectra = FlatFieldCorrection(samplesToCorrect, FlatField)
        CorrectedSpectra = cell(size(samplesToCorrect));
        NormedFlat = FlatField.Y / UsefulFunctions.ComputeIntegral(FlatField,FlatField.X(end), FlatField.X(1));
        
        for i = 1:length(samplesToCorrect)
            currentSample = samplesToCorrect{i};
            currentSample.Y = currentSample.Y ./ NormedFlat;
            CorrectedSpectra{i} = currentSample;
        end
    end
    
    function CorrectedSpectra = SubstractLinearBG(samplesToCorrect, X1, X2)
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
            CorrectedSpectra{i} = currentSample;
        end
    end

    function plotRaman(SamplesToPlot, offset)
        % Create a figure for the plot
        figure;
        % Iterate over each sample
        for sampleIdx = 1:length(SamplesToPlot)
            currentSample = SamplesToPlot{sampleIdx};
                % Get the current sample, X values, and Y values
                currentX = currentSample.X;
                currentY = currentSample.Y - offset*sampleIdx;
                currentN = currentSample.N;
                plot(currentX, currentY, 'DisplayName', currentN,'LineWidth', 1.3);
                hold on; % Add spectra to the same plot
        end

        % Add labels and legend
        xlabel('Raman Shift (cm^{-1})');
        ylabel('Intesity (a.u.)');
        title('Raman Spectra');
        legend('show');
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

    function DS = clip_spectrum(DS,ClipLeft0,ClipRight0)
    %function clips spectrum at the left and right sides by calculating its derivative; 
    %ClipLeft and ClipRight determine number of pixels to remove from the center of left and right derivative peaks,respectively; 
    %By default, ClipLeft = 8, ClipRight = 9;

    load InclinationCoeff.mat IncPoly3;
    LeftLim = 130; %all peaks of the spectrum's derivative curve should be to the left of LeftLim - determined visually from graph with all derivative curves
    RightLim = 900; %all peaks are to the right of RightLim

    ClipRangeL = ClipLeft0; %number of pixels to add/substract to the deriv. peak
    ClipRangeR = ClipRight0;

    arb = 3; %tolerance = PeakHeight/AverageBaselineValue; determines how high should be the derivative peak, in respect to average value, to be treated as real; 
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
    %% Absorption Functions
    
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
                fieldName = ['DATA_', strrep(dataset_name{7}, '.csv', '')];
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

    function CorrectedSpectra = SubstractAbsBG(samplesToCorrect, LL1, UL1, LL2, UL2)
        CorrectedSpectra = cell(size(samplesToCorrect));
        for i = 1:length(samplesToCorrect)
            currentSample = samplesToCorrect{i};
            % Extract wavelength and absorption data
            X = currentSample.X;
            Y = currentSample.Y;
            % Find the indices corresponding to the specified wavelength ranges
            idx_range1 = find(X >= LL1 & X <= UL1);
            idx_range2 = find(X >= LL2 & X <= UL2);
            % Find the wavelength where the minimum absorption occurs in each range
            [~, min_idx_range1] = min(Y(idx_range1));
            [~, min_idx_range2] = min(Y(idx_range2));
            % Get the corresponding wavelengths
            lambda_min_range1 = X(idx_range1(min_idx_range1));
            lambda_min_range2 = X(idx_range2(min_idx_range2));     
            % Extract data around the two minima
            X_range1 = X(idx_range1);
            Y_range1 = Y(idx_range1);
            X_range2 = X(idx_range2);
            Y_range2 = Y(idx_range2);

            % Fit a straight line through the two minima
            fitfunc = @(p, x) p(1) ./ x;
            initialGuess = [100]; % Initial guess for the fitting parameter A
            [fitparams, ~] = lsqcurvefit(fitfunc, initialGuess, [X_range1; X_range2], [Y_range1; Y_range2]);
            background = fitfunc(fitparams, X);

            % Subtract the background
            currentSample.Y = currentSample.Y - background;
            CorrectedSpectra{i} = currentSample;
        end
    end

    function flattenedData = FlattenSpectra(DataStructure, points)
        % Get the fieldnames of the data structure
        sampleNames = fieldnames(DataStructure);

        % Loop through each sample in the data structure
        for i = 1:numel(sampleNames)
            sample = sampleNames{i};

            % Get the wavelength and absorption data for the current sample
            selected_Y = DataStructure.(sample).Y(ismember(DataStructure.(sample).X, points));
            selected_X = DataStructure.(sample).X(ismember(DataStructure.(sample).X, points));

            % Define the objective function for least-squares optimization
            objective = @(a) sum((selected_Y - a ./ selected_X).^2);

            % Choose initial value of 'a'
            initialA = 1; % or any initial value

            % Perform least-squares optimization to find the optimal value of 'a'
            optimalA = fmincon(objective, initialA, [], [], [], [], [], [], []);
            % Subtract the background from the absorption data within the specified range
            A_flattened = DataStructure.(sample).Y - optimalA ./ DataStructure.(sample).X;

            % Update the absorption data in the data structure

            DataStructure.(sample).Y = A_flattened;
        end

        % Return the modified data structure
        flattenedData = DataStructure;
    end

    function plotAbsorption(SamplesToPlot, offset)
        % Create a figure for the plot
        figure;
        % Iterate over each sample
        for sampleIdx = 1:length(SamplesToPlot)
            currentSample = SamplesToPlot{sampleIdx};
                % Get the current sample, X values, and Y values
                currentX = currentSample.X;
                currentY = currentSample.Y - offset*sampleIdx;
                currentN = currentSample.N;
                plot(currentX, currentY, 'DisplayName', currentN,'LineWidth', 1.3);
                hold on; % Add spectra to the same plot
        end

        % Add labels and legend
        xlabel('Wavelenght (nm)');
        ylabel('Absorption (a.u.)');
        title('Absorption Spectra');
        legend('show');
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
       
   end
end