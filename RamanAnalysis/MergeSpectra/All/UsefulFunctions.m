classdef UsefulFunctions
    
   methods (Static)
       
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
        
    %Raman Functions

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
                    structure.(sampleName).Y = Y;
                    structure.(sampleName).N = sampleName;
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
            currentSample.Y = currentSample.Y - background';
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
    
    %Absorption Functions
    
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