clc;
clear;
addpath('X:\SWCNTs');
% addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\');
import UsefulFunctions.*;

% rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements Data\Absorption\';
rootpath = 'X:\Measurements Data\Raman\';

%All paths as default
path_FS514 = [rootpath,'20241007\'];
path_F0514 = [rootpath,'20241008\'];
path_R2 = [rootpath,'20241213\'];

%Select the paths of interest

paths = {
    path_FS514
    path_F0514
    path_R2
    };


ReadRamanFromPaths(paths, 3);


%%%--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DATA_20241007.LL514A.N='Laser';
DATA_20241007.F2L514R.N='Film PCE@SWCNT';
DATA_20241007.F3L514R.N='Film TCE@SWCNT';
DATA_20241007.F4L514R.N='Film TEMED@SWCNT';
DATA_20241007.F5L514R.N='Film TDEA@SWCNT';
DATA_20241007.F6L514R.N='Film Hexadecane@SWCNT';
DATA_20241007.F7L514R.N='Film Dodecane@SWCNT';
DATA_20241007.LL514B.N='Laser';
DATA_20241007.F2L514GD.N='Film PCE@SWCNT';
DATA_20241007.F3L514GD.N='Film TCE@SWCNT';
DATA_20241007.F4L514GD.N='Film TEMED@SWCNT';
DATA_20241007.F5L514GD.N='Film TDEA@SWCNT';
DATA_20241007.F6L514GD.N='Film Hexadecane@SWCNT';
DATA_20241007.F7L514GD.N='Film Dodecane@SWCNT';
DATA_20241007.LL514C.N='Laser';
DATA_20241007.F2L514DD.N='Film PCE@SWCNT';
DATA_20241007.F3L514DD.N='Film TCE@SWCNT';
DATA_20241007.F4L514DD.N='Film TEMED@SWCNT';
DATA_20241007.F5L514DD.N='Film TDEA@SWCNT';
DATA_20241007.F6L514DD.N='Film Hexadecane@SWCNT';
DATA_20241007.F7L514DD.N='Film Dodecane@SWCNT';
DATA_20241007.F0L514R.N='Sapphire Substrate';


DATA_20241008.F0L514A.N='Sapphire SubstrateA';
DATA_20241008.FFL514A.N='FlatFieldA';
DATA_20241008.F0L514B.N='Sapphire SubstrateB';
DATA_20241008.FFL514B.N='FlatFieldB';
DATA_20241008.F0L514C.N='Sapphire SubstrateC';
DATA_20241008.FFL514C.N='FlatFieldC';
DATA_20241008.F0L514D.N='Sapphire SubstrateD';
DATA_20241008.FFL514D.N='FlatFieldD';
DATA_20241008.F0L514E.N='Sapphire SubstrateE';
DATA_20241008.FFL514E.N='FlatFieldE';
DATA_20241008.F0L514F.N='Sapphire SubstrateF';
DATA_20241008.FFL514F.N='FlatFieldF';
DATA_20241008.F0L514G.N='Sapphire SubstrateG';
DATA_20241008.FFL514G.N='FlatFieldG';

DATA_20241213.FS5L514G.N='Film 2 (S5 TDAE@CNTs)';
DATA_20241213.FS5L514D.N='Film 2 (S5 TDAE@CNTs)';
DATA_20241213.FS7L514G.N='Film 1 (S7 Dodecane@CNTs)';
DATA_20241213.FS7L514D.N='Film 1 (S7 Dodecane@CNTs)';
DATA_20241213.FS8L514G.N='Film 5 (S8 PCE@CNTs)';
DATA_20241213.FS8L514D.N='Film 5 (S8 PCE@CNTs)';
DATA_20241213.FS9L514G.N='Film 4 (S9 TEMED@CNTs)';
DATA_20241213.FS9L514D.N='Film 4 (S9 TEMED@CNTs)';
DATA_20241213.FF6L514G.N='Film 3 (SFF6 TTF@CNTs)';
DATA_20241213.FF6L514D.N='Film 3 (SFF6 TTF@CNTs)';


%%%--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%--------SAMPLE COMPARISION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%RBMs
FS514R = {
    DATA_20241007.F2L514R
    DATA_20241007.F3L514R
    DATA_20241007.F4L514R
    DATA_20241007.F5L514R
    DATA_20241007.F6L514R
    DATA_20241007.F7L514R
%     DATA_20241007.FFL514R
    };
% FS514R = FlatFieldCorrectionPixelWise(FS514R, DATA_20241007.FFL514R);

WL = 514.5;
for i=1:length(FS514R)
    current = FS514R{i};  % Access the cell array element once
%     current = filere(current, 10,10);
    current = remove_inclination(current, WL);
    current = correct_instrument_response(current, WL);
    current = remove_bg_poly(current);
    current = remove_baseline_polynomial(current, 0);

    FS514R{i} = current;  % Save the result back to the cell array
end

% FS514R = NormalizeSample(FS514R, 138, 166);

% plotRaman(FS514R,0,WL)

% %%%%%GBand
% FS514GD = {
%     DATA_20241007.F2L514GD
%     DATA_20241007.F3L514GD
%     DATA_20241007.F4L514GD
%     DATA_20241007.F5L514GD
%     DATA_20241007.F6L514GD
%     DATA_20241007.F7L514GD
% %     DATA_20241007.FFL514GD
%     };
% % FS514GD = FlatFieldCorrectionPixelWise(FS514GD, DATA_20241007.FFL514GD);
% 
% 
% WL = 514.5;
% for i=1:length(FS514GD)
%     current = FS514GD{i};  % Access the cell array element once
%     
%     current = clip_spectrum(current, 10,10);
%     current = remove_inclination(current, WL);
%     current = correct_instrument_response(current, WL);
%     current = remove_bg_poly(current);
%     current = remove_baseline_polynomial(current, 0);
% 
%     FS514GD{i} = current;  % Save the result back to the cell array
% end
% 
% % FS514GD = NormalizeSample(FS514GD, 0, 20000);
% 
% % plotRaman(FS514GD,0,WL)
% 
% 
% %%%%%2DBand
% FS514DD = {
%     DATA_20241007.F2L514DD
%     DATA_20241007.F3L514DD
%     DATA_20241007.F4L514DD
%     DATA_20241007.F5L514DD
%     DATA_20241007.F6L514DD
%     DATA_20241007.F7L514DD
% %     DATA_20241007.FFL514DD
% %     DATA_20241008.FFL514G
%     };
% 
%    %% Flat field correction indetail 
% % figure()
% % 
% % hold on
% % for i=1:length(FS514DD)
% %     current = FS514DD{i}
% %     plot(current.P, current.Y/(current.Y(200)))
% % end
% % hold off
% %     
% % figure()
% % 
% % hold on
% % for i=1:length(FS514DD)
% %     current = FS514DD{i}
% %     plot(current.P, current.Y/DATA_20241008.FFL514G.Y)
% % end
% % hold off
% %%
% 
% % FS514DD = FlatFieldCorrectionPixelWise(FS514DD, DATA_20241007.FFL514DD);
% % 
% % plotRaman(FS514DD,0, WL)
% 
% WL = 514.5;
% for i=1:length(FS514DD)
%     current = FS514DD{i};  % Access the cell array element once
%     current = clip_spectrum(current, 10,10);
%     current = remove_inclination(current, WL);
%     current = correct_instrument_response(current, WL);
%     current = remove_bg_poly(current);
%     current = remove_baseline_polynomial(current, 0);
%     FS514DD{i} = current;  % Save the result back to the cell array
% end
% 
% FS514DD = NormalizeSample(FS514DD, 0, 20000);
% % FS514DD{1}.Y = FS514DD{1}.Y-0.04
% % FS514DD{2}.Y = FS514DD{2}.Y-0.04
% % % FS514DD{3}.Y = FS514DD{3}.Y-0.04
% % FS514DD{4}.Y = FS514DD{4}.Y-0.04
% % FS514DD{5}.Y = FS514DD{5}.Y-0.04
% FS514DD{6}.Y = FS514DD{6}.Y-0.055
% FS514DD = NormalizeSample(FS514DD, 0, 20000);
% 
% plotRaman(FS514DD,0, WL)
% 
%%%%%%%%%%%%%%%%%%%%%%%NEW BUT FOR ROUND 2

%%%%%GBand
GBAND = {
    DATA_20241213.FS5L514G
    DATA_20241213.FS7L514G
    DATA_20241213.FS8L514G
    DATA_20241213.FS9L514G
    DATA_20241213.FF6L514G
    };
GBAND = FlatFieldCorrection(GBAND, DATA_20241213.FFL514GB);


GBAND = FilterDataByXRange(GBAND, 1250, 1660)
WL = 514.5;

for i=1:length(GBAND)
    current = GBAND{i};  % Access the cell array element once
    current = remove_inclination(current, WL);
    current = correct_instrument_response(current, WL);
    current = remove_bg_poly(current);
    current = remove_baseline_polynomial(current, 0);

    GBAND{i} = current;  % Save the result back to the cell array
end

GBAND = NormalizeSample(GBAND, 0, 20000);

% plotRaman(GBAND,0,WL)

% %%%%%2DBand
DDBANDS = {
    DATA_20241213.FS5L514D
    DATA_20241213.FS7L514D
    DATA_20241213.FS8L514D
    DATA_20241213.FS9L514D
    DATA_20241213.FF6L514D

    };

DDBANDS = FlatFieldCorrection(DDBANDS, DATA_20241213.FFL514D);
DDBANDS = FilterDataByXRange(DDBANDS, 2480, 2845)

for i=1:length(DDBANDS)
    current = DDBANDS{i};  % Access the cell array element once
    current = remove_inclination(current, WL);
    current = correct_instrument_response(current, WL);
    current = remove_bg_poly(current);
    current = remove_baseline_polynomial(current, 0);

    DDBANDS{i} = current;  % Save the result back to the cell array
end

DDBANDS = NormalizeSample(DDBANDS, 0, 20000);
plotRaman(DDBANDS,0, WL)
% 
% WL = 514.5;
% for i=1:length(FS514DD)
%     current = FS514DD{i};  % Access the cell array element once
%     current = clip_spectrum(current, 10,10);
%     current = remove_inclination(current, WL);
%     current = correct_instrument_response(current, WL);
%     current = remove_bg_poly(current);
%     current = remove_baseline_polynomial(current, 0);
%     FS514DD{i} = current;  % Save the result back to the cell array
% end
% 
% FS514DD = NormalizeSample(FS514DD, 0, 20000);
% % FS514DD{1}.Y = FS514DD{1}.Y-0.04
% % FS514DD{2}.Y = FS514DD{2}.Y-0.04
% % % FS514DD{3}.Y = FS514DD{3}.Y-0.04
% % FS514DD{4}.Y = FS514DD{4}.Y-0.04
% % FS514DD{5}.Y = FS514DD{5}.Y-0.04
% FS514DD{6}.Y = FS514DD{6}.Y-0.055
% FS514DD = NormalizeSample(FS514DD, 0, 20000);
% 
% plotRaman(FS514DD,0, WL)



%%%--------SUBSTRATE PLOT--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DATA_20241008.F0L514A = FlatFieldCorrectionSingle(DATA_20241008.F0L514A, DATA_20241008.FFL514A);
DATA_20241008.F0L514B = FlatFieldCorrectionSingle(DATA_20241008.F0L514B, DATA_20241008.FFL514B);
DATA_20241008.F0L514C = FlatFieldCorrectionSingle(DATA_20241008.F0L514C, DATA_20241008.FFL514C);
DATA_20241008.F0L514D = FlatFieldCorrectionSingle(DATA_20241008.F0L514D, DATA_20241008.FFL514D);
DATA_20241008.F0L514E = FlatFieldCorrectionSingle(DATA_20241008.F0L514E, DATA_20241008.FFL514E);
DATA_20241008.F0L514F = FlatFieldCorrectionSingle(DATA_20241008.F0L514F, DATA_20241008.FFL514F);
DATA_20241008.F0L514G = FlatFieldCorrectionSingle(DATA_20241008.F0L514G, DATA_20241008.FFL514G);

F0514 = {
    DATA_20241008.F0L514A
    DATA_20241008.F0L514B
    DATA_20241008.F0L514C
    DATA_20241008.F0L514D
    DATA_20241008.F0L514E
    DATA_20241008.F0L514F
    DATA_20241008.F0L514G
    };

WL = 514.5;
for i=1:length(F0514)
    current = F0514{i};  % Access the cell array element once
    current = clip_spectrum(current, 30,30);
    current = remove_inclination(current, WL);
    current = correct_instrument_response(current, WL);
    current = remove_bg_poly(current);
    F0514{i} = current;  % Save the result back to the cell array
end


SAPSUBS = concatenateSpectra({F0514{1},F0514{2},F0514{3},F0514{4},F0514{5},F0514{6},F0514{7}}, 'Sapphire Substrate');


% plotRaman({SAPSUBS}, 0.0, WL);


ALLFF = {
    
DATA_20241008.FFL514A
DATA_20241008.FFL514B
DATA_20241008.FFL514C
DATA_20241008.FFL514D
DATA_20241008.FFL514E
DATA_20241008.FFL514F
DATA_20241008.FFL514G
    }

for i=1:length(ALLFF)
    current = ALLFF{i};  % Access the cell array element once
%     current = clip_spectrum(current, 20,20);
    current = remove_inclination(current, WL);
    current = correct_instrument_response(current, WL);
    current = remove_bg_poly(current);
    current = remove_baseline_polynomial(current, 0);
    ALLFF{i} = current;  % Save the result back to the cell array
end

% plotRamanPxl(ALLFF,0)
% plotRaman(ALLFF,0)


% plotRamanFit(FS514R, 2)


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

