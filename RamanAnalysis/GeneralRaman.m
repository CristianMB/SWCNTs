clc;
clear;

rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements\Raman\';
%rootpath = 'C:\Users\cborja\OneDrive - Universiteit Antwerpen\Measurements\Raman\';

%All paths as default
path_20240111 = [rootpath,'20240111\'];
path_20240320 = [rootpath,'20240320\'];
path_20240321 = [rootpath,'20240321\'];
path_20240325 = [rootpath,'20240325\'];

%Select the paths of interest
paths = {
    path_20240111,
    path_20240320,
    path_20240321,
    path_20240325
    };

%Read and structure data from the paths
ReadFromPaths(paths);


%%%--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DATA_20240111.FLATHD.N='FlatField 514.5nm';
DATA_20240111.LL514.N='LaserLine 514.5nm';
DATA_20240111.LL514HD.N='LaserLine 514.5nm';


DATA_20240111.S240111A.N='SF D2O@SWCNT (514.5 nm)';
DATA_20240111.S240111C.N='SF Methanol@SWCNT (514.5 nm)';
DATA_20240111.S240111F.N='SF PCE@SWCNT (514.5 nm)';
DATA_20240111.S240111G.N='SF PCE@SWCNT (514.5 nm)';
DATA_20240111.S240111H.N='SF PCE@SWCN (514.5 nm)';
DATA_20240111.S240111B.N='SF TCE@SWCNT (514.5 nm)';
DATA_20240111.S240111BB.N='SF TCE@SWCNT (514.5 nm)';
DATA_20240111.S240111D.N='SF TCE@SWCNT (514.5 nm)';
DATA_20240111.S240111E.N='SF TTF@SWCNT (514.5 nm)';
DATA_20240111.S240111I.N='SF TEMED@SWCNT (514.5 nm)';

DATA_20240111.S240111S.N='SC Empty@SWCNT (514.5 nm)';
DATA_20240111.S240111J.N='SF D2O@SWCNT (514.5 nm)';
DATA_20240111.S240111L.N='SF Methanol@SWCNT (514.5 nm)';
DATA_20240111.S240111O.N='SF PCE@SWCNT (514.5 nm)';
DATA_20240111.S240111P.N='SF PCE@SWCNT (514.5 nm)';
DATA_20240111.S240111Q.N='SF PCE@SWCN (514.5 nm)';
DATA_20240111.S240111K.N='SF TCE@SWCNT (514.5 nm)';
DATA_20240111.S240111KK.N='SF TCE@SWCNT (514.5 nm)';
DATA_20240111.S240111M.N='SF TCE@SWCNT (514.5 nm)';
DATA_20240111.S240111N.N='SF TTF@SWCNT (514.5 nm)';
DATA_20240111.S240111R.N='SF TEMED@SWCNT (514.5 nm)';

DATA_20240320.FFH520R.N='FlatField (520.25 nm)';
DATA_20240320.L520A.N='LaserLine (520.25 nm)';

DATA_20240320.S2H520R.N='CB PCE@SWCNT Dial. DGU C (Filled) (520.25 nm)';
DATA_20240320.S3H520R.N='CB TCE@SWCNT Dial. DGU C (Filled) (520.25 nm)';
DATA_20240320.S4H520R.N='CB TEMED@SWCNT Dial. DGU C (Filled) (520.25 nm)';
DATA_20240320.S5H520R.N='CB TDAE@SWCNT Dial. DGU C (Filled) (520.25 nm)';
DATA_20240320.S6H520R.N='CB Hexadecane@SWCNT Dial. DGU C (Filled) (520.25 nm)';
DATA_20240320.S7H520R.N='CB Dodecane@SWCNT Dial. DGU C (Filled) (520.25 nm)';

DATA_20240321.FFH520G.N='FlatField (520.25 nm)';
DATA_20240321.FFL520G.N='FlatField (520.25 nm)';
DATA_20240321.L520HG.N='LaserLine HD Before (520.25 nm)';
DATA_20240321.L520L.N='LaserLine LD (520.25 nm)';
DATA_20240321.LLH520GB.N='LaserLine HD After (520.25 nm)';

DATA_20240321.S2H520G.N='CB PCE@SWCNT Dial. DGU C (Filled) (520.25 nm)';
DATA_20240321.S2L520G.N='CB PCE@SWCNT Dial. DGU C (Filled) (520.25 nm)';
DATA_20240321.S3H520G.N='CB TCE@SWCNT Dial. DGU C (Filled) (520.25 nm)';
DATA_20240321.S3L520G.N='CB TCE@SWCNT Dial. DGU C (Filled) (520.25 nm)';
DATA_20240321.S4H520G.N='CB TEMED@SWCNT Dial. DGU C (Filled) (520.25 nm)';
DATA_20240321.S4L520G.N='CB TEMED@SWCNT Dial. DGU C (Filled) (520.25 nm)';
DATA_20240321.S5H520G.N='CB TDAE@SWCNT Dial. DGU C (Filled) (520.25 nm)';
DATA_20240321.S5L520G.N='CB TDAE@SWCNT Dial. DGU C (Filled) (520.25 nm)';
DATA_20240321.S6H520G.N='CB Hexadecane@SWCNT Dial. DGU C (Filled) (520.25 nm)';
DATA_20240321.S6L520G.N='CB Hexadecane@SWCNT Dial. DGU C (Filled) (520.25 nm)';
DATA_20240321.S7H520G.N='CB Dodecane@SWCNT Dial. DGU C (Filled) (520.25 nm)';
DATA_20240321.S7L520G.N='CB Dodecane@SWCNT Dial. DGU C (Filled) (520.25 nm)';

DATA_20240325.REH680R.N='SC Empty@SWCNT (680.02 nm)';
DATA_20240325.RWH680R.N='SC Water@SWCNT (680.02 nm)';
DATA_20240325.S2H680R.N='CB PCE@SWCNT Dial. DGU C (Filled) (680.02 nm)';
DATA_20240325.S3H680R.N='CB TCE@SWCNT Dial. DGU C (Filled) (680.02 nm)';
DATA_20240325.S4H680R.N='CB TEMED@SWCNT Dial. DGU C (Filled) (680.02 nm)';
DATA_20240325.S5H680R.N='CB TDAE@SWCNT Dial. DGU C (Filled) (680.02 nm)';
DATA_20240325.S6H680R.N='CB Hexadecane@SWCNT Dial. DGU C (Filled) (680.02 nm)';
DATA_20240325.S7H680R.N='CB Dodecane@SWCNT Dial. DGU C (Filled) (680.02 nm)';




%%%--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FlatFields Normalization
DATA_20240111.FLATHD = NormFlatField(DATA_20240111.FLATHD);
DATA_20240320.FFH520R = NormFlatField(DATA_20240320.FFH520R);
DATA_20240321.FFH520G = NormFlatField(DATA_20240321.FFH520G);
DATA_20240321.FFL520G = NormFlatField(DATA_20240321.FFL520G);
DATA_20240325.FH680R = NormFlatField(DATA_20240325.FH680R);


%FF Correction to RBMs region
DATA_20240111.S240111J.Y = DATA_20240111.S240111J.Y ./ DATA_20240111.FLATHD.Y;
DATA_20240111.S240111K.Y = DATA_20240111.S240111K.Y ./ DATA_20240111.FLATHD.Y;
DATA_20240111.S240111KK.Y = DATA_20240111.S240111KK.Y ./ DATA_20240111.FLATHD.Y;
DATA_20240111.S240111L.Y = DATA_20240111.S240111L.Y ./ DATA_20240111.FLATHD.Y;
DATA_20240111.S240111M.Y =DATA_20240111.S240111M.Y ./ DATA_20240111.FLATHD.Y;
DATA_20240111.S240111N.Y =DATA_20240111.S240111N.Y ./ DATA_20240111.FLATHD.Y; 
DATA_20240111.S240111O.Y =DATA_20240111.S240111O.Y ./ DATA_20240111.FLATHD.Y;
DATA_20240111.S240111P.Y =DATA_20240111.S240111P.Y ./ DATA_20240111.FLATHD.Y; 
DATA_20240111.S240111Q.Y =DATA_20240111.S240111Q.Y ./ DATA_20240111.FLATHD.Y;
DATA_20240111.S240111R.Y =DATA_20240111.S240111R.Y ./ DATA_20240111.FLATHD.Y;
DATA_20240111.S240111S.Y =DATA_20240111.S240111S.Y ./ DATA_20240111.FLATHD.Y;

%FF Correction to RBMs region
DATA_20240320.S2H520R.Y = DATA_20240320.S2H520R.Y ./ DATA_20240320.FFH520R.Y;
DATA_20240320.S3H520R.Y = DATA_20240320.S3H520R.Y ./ DATA_20240320.FFH520R.Y;
DATA_20240320.S4H520R.Y = DATA_20240320.S4H520R.Y ./ DATA_20240320.FFH520R.Y;
DATA_20240320.S5H520R.Y = DATA_20240320.S5H520R.Y ./ DATA_20240320.FFH520R.Y;
DATA_20240320.S6H520R.Y = DATA_20240320.S6H520R.Y ./ DATA_20240320.FFH520R.Y;
DATA_20240320.S7H520R.Y = DATA_20240320.S7H520R.Y ./ DATA_20240320.FFH520R.Y;

%FF Correction to GBand region in LD Mode
DATA_20240321.S2L520G.Y = DATA_20240321.S2L520G.Y ./ DATA_20240321.FFL520G.Y;
DATA_20240321.S3L520G.Y = DATA_20240321.S3L520G.Y ./ DATA_20240321.FFL520G.Y;
DATA_20240321.S4L520G.Y = DATA_20240321.S4L520G.Y ./ DATA_20240321.FFL520G.Y;
DATA_20240321.S5L520G.Y = DATA_20240321.S5L520G.Y ./ DATA_20240321.FFL520G.Y;
DATA_20240321.S6L520G.Y = DATA_20240321.S6L520G.Y ./ DATA_20240321.FFL520G.Y;
DATA_20240321.S7L520G.Y = DATA_20240321.S7L520G.Y ./ DATA_20240321.FFL520G.Y;

%FF Correction to GBand region in HD Mode
DATA_20240321.S2H520G.Y = DATA_20240321.S2H520G.Y ./ DATA_20240321.FFH520G.Y;
DATA_20240321.S3H520G.Y = DATA_20240321.S3H520G.Y ./ DATA_20240321.FFH520G.Y;
DATA_20240321.S4H520G.Y = DATA_20240321.S4H520G.Y ./ DATA_20240321.FFH520G.Y;
DATA_20240321.S5H520G.Y = DATA_20240321.S5H520G.Y ./ DATA_20240321.FFH520G.Y;
DATA_20240321.S6H520G.Y = DATA_20240321.S6H520G.Y ./ DATA_20240321.FFH520G.Y;
DATA_20240321.S7H520G.Y = DATA_20240321.S7H520G.Y ./ DATA_20240321.FFH520G.Y;

%FF Correction to GBand region in HD Mode
DATA_20240325.REH680R.Y = DATA_20240325.REH680R.Y ./ DATA_20240325.FH680R.Y;
DATA_20240325.RWH680R.Y = DATA_20240325.RWH680R.Y ./ DATA_20240325.FH680R.Y;
DATA_20240325.S2H680R.Y = DATA_20240325.S2H680R.Y ./ DATA_20240325.FH680R.Y;
DATA_20240325.S3H680R.Y = DATA_20240325.S3H680R.Y ./ DATA_20240325.FH680R.Y;
DATA_20240325.S4H680R.Y = DATA_20240325.S4H680R.Y ./ DATA_20240325.FH680R.Y;
DATA_20240325.S5H680R.Y = DATA_20240325.S5H680R.Y ./ DATA_20240325.FH680R.Y;
DATA_20240325.S6H680R.Y = DATA_20240325.S6H680R.Y ./ DATA_20240325.FH680R.Y;
DATA_20240325.S7H680R.Y = DATA_20240325.S7H680R.Y ./ DATA_20240325.FH680R.Y;



%%%--------SAMPLE COMPARISION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


GDBand514 = {
        DATA_20240111.S240111A,
        DATA_20240111.S240111B,
        DATA_20240111.S240111BB,
        DATA_20240111.S240111C,
        DATA_20240111.S240111D,
        DATA_20240111.S240111E,
        DATA_20240111.S240111F,
        DATA_20240111.S240111G,
        DATA_20240111.S240111H,
        DATA_20240111.S240111I
    };
    
GDBand520 = {
        DATA_20240321.S2L520G,
        DATA_20240321.S3L520G,
        DATA_20240321.S4L520G,
        DATA_20240321.S5L520G,
        DATA_20240321.S6L520G,
        DATA_20240321.S7L520G,
                
                
        };
    
GBand520 = {
        DATA_20240321.S2H520G,
        DATA_20240321.S3H520G,
        DATA_20240321.S4H520G,
        DATA_20240321.S5H520G,
        DATA_20240321.S6H520G,
        DATA_20240321.S7H520G,
        };


    
RBM514 = {
        %Empty SC
        DATA_20240111.S240111S
        %D2O Filled
        DATA_20240111.S240111J,    
        %PCE
        DATA_20240111.S240111O,
        DATA_20240111.S240111P,
        DATA_20240111.S240111Q,
        %TCE
        DATA_20240111.S240111K,
        DATA_20240111.S240111M,
        %TEMED
        DATA_20240111.S240111R,
        %TTF
        DATA_20240111.S240111N,       
        };

RBM520 = {  
        %Empty
        DATA_20240111.S240111S
        %D2O Filled
        DATA_20240111.S240111J,

        DATA_20240320.S2H520R,
        DATA_20240320.S3H520R,
        DATA_20240320.S4H520R,
        DATA_20240320.S5H520R,
        DATA_20240320.S6H520R,
        DATA_20240320.S7H520R,
    };
    
RBM680 = {
        DATA_20240325.REH680R,
        DATA_20240325.RWH680R,

        DATA_20240325.S2H680R,
        DATA_20240325.S3H680R,
        DATA_20240325.S4H680R,
        DATA_20240325.S5H680R,
        DATA_20240325.S6H680R,
        DATA_20240325.S7H680R
        };
    
    
LG = 1560;
HG = 1570;

GDBand520 = SubstractBG(GDBand520, 1400, 1500)
GDBand520 = Normalize(GDBand520,LG, HG);

GBand520 = SubstractBG(GBand520, 1540, 1620)
GBand520 = Normalize(GBand520,LG, HG);

GDBand514 = SubstractBG(GDBand514,1400, 1500)
GDBand514 = Normalize(GDBand514,LG,HG)


LR = 150;
HR = 160;

RBM514 = SubstractBG(RBM514, 130, 210)
RBM514 = Normalize(RBM514,LR, HR);

RBM520 = SubstractBG(RBM520, 130, 210)
RBM520 = Normalize(RBM520,LR, HR);

RBM680 = SubstractBG(RBM680,150, 198)
RBM680 = Normalize(RBM680,LR, HR)




%plotSampleList(RBM514, 1)
%plotSampleList(RBM520, 1)
%plotSampleList(RBM680, 1)
%plotSampleList(GDBand514, 1)
%plotSampleList(GDBand520, 1)
%plotSampleList(GBand520, 1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dataStructures = ReadFromPaths(paths)
    dataStructures = struct();  
    for p = 1:length(paths)
        % Extract the suffix from the path
        try
            % Extract the suffix from the path variable
            dataset_name = strsplit(paths{p}, "\");
            fieldName = ['DATA_', strrep(dataset_name{7}, '.', '')];
           
            dirInfo = dir(fullfile(paths{p}, '*.m3d'));
            fileList = {dirInfo(~[dirInfo.isdir]).name};
                        
            structure = struct();

            % Read the data from the current path
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

function integralValue = computeIntegral(sample, lowerLimit, upperLimit)
    % Check if the sampleName is in DATA
    % Extract X and Y values for the specified sample
    x = sample.X;
    y = sample.Y;
    % Define the function to integrate
    f = @(xi) interp1(x, y, xi, 'pchip');

    % Calculate the integral
    integralValue = integral(f, lowerLimit, upperLimit);
end

function maximumValue = computeMaximum(sample, lowerLimit, upperLimit)
    % Check if the sampleName is in DATA
    % Extract X and Y values for the specified sample
    x = sample.X;
    y = sample.Y;
    % Define the function to integrate
    
    indicesInRange = find(x >= lowerLimit & x <= upperLimit);
    yInRange = y(indicesInRange);

    maximumValue = max(yInRange);
end

function Normed = NormFlatField(Flat)
        Flat.Y = Flat.Y / computeIntegral(Flat,Flat.X(end), Flat.X(1));
        Normed = Flat;
end

function NormedSamples = Normalize(samplesToNormalize, lowerLimit, upperLimit)
    NormedSamples = cell(size(samplesToNormalize));
    % Iterate over each sample to be normalized
    for sampleIdx = 1:length(samplesToNormalize)
        currentSample = samplesToNormalize{sampleIdx};
        %currentSample.Y = currentSample.Y/computeIntegral(currentSample,lowerLimit, upperLimit);
        currentSample.Y = currentSample.Y/computeMaximum(currentSample,lowerLimit, upperLimit);
        NormedSamples{sampleIdx} = currentSample;
    end
end

function CorrectedSpectra = SubstractBG(samplesToCorrect, X1, X2)
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

function plotSampleList(SamplesToPlot, offset)
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






