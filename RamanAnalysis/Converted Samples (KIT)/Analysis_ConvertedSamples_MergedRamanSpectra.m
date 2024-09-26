clc;
clear;
addpath('X:\SWCNTs');

addpath('X:\Measurements Data\Raman');
%addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\');
import UsefulFunctions.*;

%rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements Data\Raman\';
rootpath = 'X:\Measurements Data\Raman';

%All paths as default
path_20240426 = [rootpath,'\20240426\'];
path_20240514 = [rootpath,'\20240514\'];
path_20240515 = [rootpath,'\20240515\'];
path_20240614 = [rootpath,'\20240614\'];
path_20240620 = [rootpath,'\20240620\'];
path_20240628 = [rootpath,'\20240628\'];

%Select the paths of interest

paths = {
    path_20240426
    path_20240514
    path_20240515
    path_20240614
    path_20240620
    path_20240628
    };

ReadRamanFromPaths(paths);

%%% --------------------------------------------------------
%%%           Correcting for sensitivity detector and inclination and
%%%           cutting edges away
%%% --------------------------------------------------------

%%  DATA 20240426 FIRST TEST SAMPLES JUST RBMS, D AND G BAND  AT 650nm

%Time normalization
DATA_20240426.BAL650R.Y = DATA_20240426.BAL650R.Y/60
DATA_20240426.BAL650RB.Y = DATA_20240426.BAL650RB.Y/60
DATA_20240426.BAL650D.Y = DATA_20240426.BAL650D.Y/60
DATA_20240426.BAL650G.Y = DATA_20240426.BAL650G.Y/60

DATA_20240426.BBL650RA.Y = DATA_20240426.BBL650RA.Y/30
DATA_20240426.BBL650RB.Y =  DATA_20240426.BBL650RB.Y/30
DATA_20240426.BBL650D.Y = DATA_20240426.BBL650D.Y/30
DATA_20240426.BBL650G.Y = DATA_20240426.BBL650G.Y/30
         
DATA_20240426.BAL650R = clipRangeEdges(DATA_20240426.BAL650R, 0,320)
DATA_20240426.BBL650RA = clipRangeEdges(DATA_20240426.BBL650RA, 120,2000)
% DATA_20240426.BAL650RB = clipRangeEdges(DATA_20240426.BAL650RB, 380,520)
DATA_20240426.BBL650RB = clipRangeEdges(DATA_20240426.BBL650RB, 300,2000)
% DATA_20240426.BBL650D = clipRangeEdges(DATA_20240426.BBL650D, 1230,1444)
% DATA_20240426.BAL650D = clipRangeEdges(DATA_20240426.BAL650D, 1230,1444)
% DATA_20240426.BAL650G = clipRangeEdges(DATA_20240426.BAL650G, 1430,1640)
% DATA_20240426.BBL650G = clipRangeEdges(DATA_20240426.BBL650G, 1430,1640)

DATA_20240426.BAL650R = remove_gaussian_profile(DATA_20240426.BAL650R, 85)

TestSamples650={
         DATA_20240426.BAL650R
         DATA_20240426.BAL650RB
         DATA_20240426.BAL650D
         DATA_20240426.BAL650G
         
         DATA_20240426.BBL650RA
         DATA_20240426.BBL650RB
         DATA_20240426.BBL650D
         DATA_20240426.BBL650G 
        };
WL = 650;

for i=1:length(TestSamples650)
    current = TestSamples650{i};  % Access the cell array element once
    current = clip_spectrum(current, 10,10);
    current = remove_inclination(current, WL);
    current = correct_instrument_response(current, WL);
    current = remove_bg_poly(current);
    TestSamples650{i} = current;  % Save the result back to the cell array
end



TSA650a = concatenateSpectra({TestSamples650{1:4}}, 'Test Film Sample (SWCNT)')
TSB650a = concatenateSpectra({TestSamples650{5:8}}, 'Test Film Sample (Converted to DWCNT)')

GroupedTS650a = {
                    TSA650a
                    TSB650a
                    }

                
             
GroupedTS650a = UsefulFunctions.NormalizeSample(GroupedTS650a,1400, 1900)                                

% plotRaman(TestSamples650, 1.0)
plotRaman(GroupedTS650a, 1)






%% DATA 20240514 FIRST TEST SAMPLES JUST G BAND AND CARBINE REGION AT 650nm
% 
% DATA_20240514.BBL650C1 = clipRangeEdges(DATA_20240514.BBL650C1, 1690,6000)
% DATA_20240514.BAL650C1 = clipRangeEdges(DATA_20240514.BAL650C1, 1690,6000)
% DATA_20240514.BBL650C2 = clipRangeEdges(DATA_20240514.BBL650C2, 1490,1705)
% DATA_20240514.BAL650C2 = clipRangeEdges(DATA_20240514.BAL650C2, 1490,1705)

TestSamplesB650={
    DATA_20240514.BAL650C2    
    DATA_20240514.BAL650C1
    
    DATA_20240514.BBL650C2
    DATA_20240514.BBL650C1
    };

WL = 650;

for i=1:length(TestSamplesB650)
    current = TestSamplesB650{i};  % Access the cell array element once
    current = clip_spectrum(current, 10,10);
    current = remove_inclination(current, WL);
    current = correct_instrument_response(current, WL);
    current = remove_bg_poly(current);
    TestSamplesB650{i} = current;  % Save the result back to the cell array
end

TSA650b = concatenateSpectra({TestSamplesB650{1:2}}, 'Test Film Sample (SWCNT)')
TSB650b = concatenateSpectra({TestSamplesB650{3:4}}, 'Test Film Sample (Converted to DWCNT)')

GroupedTS650b = {
                    TSA650b
                    TSB650b
                    }
             
GroupedTS650b = UsefulFunctions.NormalizeSample(GroupedTS650b,0, 2000)  

% plotRaman(GroupedTS650b, 1)
% plotRaman(TestSamplesB650, 5.0)


%% DATA 20240515 FIRST TEST SAMPLES FULL SPEC AT 570nm

DATA_20240515.BAL570RA = clipRangeEdges(DATA_20240515.BAL570RA, 134, 2000)
DATA_20240515.BBL570RA = clipRangeEdges(DATA_20240515.BBL570RA, 134, 2000)
% 
% DATA_20240515.BAL570RB = clipRangeEdges(DATA_20240515.BAL570RB, 420,780)
% DATA_20240515.BBL570RB = clipRangeEdges(DATA_20240515.BBL570RB, 420,780)
% 
% DATA_20240515.BAL570D = clipRangeEdges(DATA_20240515.BAL570D, 1080,1400)
% DATA_20240515.BBL570D = clipRangeEdges(DATA_20240515.BBL570D, 1080,1400)
% 
% DATA_20240515.BAL570G = clipRangeEdges(DATA_20240515.BAL570G, 1390,1700)
% DATA_20240515.BBL570G = clipRangeEdges(DATA_20240515.BBL570G, 1390,1700)
% 
% DATA_20240515.BAL570C = clipRangeEdges(DATA_20240515.BAL570C, 1690,1990)
% DATA_20240515.BBL570C = clipRangeEdges(DATA_20240515.BBL570C, 1690,1990)

TestSamples570={
        DATA_20240515.BAL570RA
        DATA_20240515.BAL570RB
        DATA_20240515.BAL570D
        DATA_20240515.BAL570G
        DATA_20240515.BAL570C

        DATA_20240515.BBL570RA
        DATA_20240515.BBL570RB
        DATA_20240515.BBL570D
        DATA_20240515.BBL570G
        DATA_20240515.BBL570C
    };

WL = 570;

for i=1:length(TestSamples570)
    current = TestSamples570{i};  % Access the cell array element once
    current = clip_spectrum(current, 10,10);
    current = remove_inclination(current, WL);
    current = correct_instrument_response(current, WL);
    current = remove_bg_poly(current);
    TestSamples570{i} = current;  % Save the result back to the cell array
end


TSA570 = concatenateSpectra({TestSamples570{1:5}}, 'Test Film Sample (SWCNT)')
TSB570 = concatenateSpectra({TestSamples570{6:10}}, 'Test Film Sample (Converted to DWCNT)')

GroupedTS570 = {
                    TSA570
                    TSB570
                    }
                
GroupedTS570 = UsefulFunctions.NormalizeSample(GroupedTS570,0, 2000)                                
% plotRaman(GroupedTS570, 1)
% plotRaman(TestSamples570, 5.0)

%% DATA 20240614 SOLUTION SAMPLES RBMs,G and D Band SPEC AT 514nm

% DATA_20240614.T1SL514R = remove_gaussian_profile(DATA_20240614.T1SL514R, 500)
% DATA_20240614.T4PL514R = remove_gaussian_profile(DATA_20240614.T4PL514R, 500)
% DATA_20240614.T4SL514R = remove_gaussian_profile(DATA_20240614.T4SL514R, 500)

SolutionSamples514={
    DATA_20240614.T1SL514R
    DATA_20240614.T1SL514G
    
    DATA_20240614.T4SL514R
    DATA_20240614.T4SL514G
            
    DATA_20240614.T4PL514R
    DATA_20240614.T4PL514G
    };

WL = 514.5;

for i=1:length(SolutionSamples514)
    current = SolutionSamples514{i};  % Access the cell array element once
    current = clip_spectrum(current, 10,10);
    current = remove_inclination(current, WL);
    current = correct_instrument_response(current, WL);
    current = remove_bg_poly(current);
    SolutionSamples514{i} = current;  % Save the result back to the cell array
end

% plotRaman(SolutionSamples514, 0.0)

SST1S570 = concatenateSpectra({SolutionSamples514{1:2}}, 'Solution Sample - T1S')
SST4S570 = concatenateSpectra({SolutionSamples514{3:4}}, 'Solution Sample - T4S')
SST4P570 = concatenateSpectra({SolutionSamples514{5:6}}, 'Solution Sample - T4P')

GroupedSS514 = {
                    SST1S570
                    SST4S570
                    SST4P570
                    }
                
GroupedSS514 = UsefulFunctions.NormalizeSample(GroupedSS514,0, 2000)                
% plotRaman(GroupedSS514, 1)

%% DATA 20240620 FILM SAMPLES RBMs,G and D Band SPEC AT 514nm


FilmSamples514={
    DATA_20240620.F1SL514R
    DATA_20240620.F1SL514G

    DATA_20240620.F1CL514R
    DATA_20240620.F1CL514G
    
    DATA_20240620.F2SL514R
    DATA_20240620.F2SL514G
    DATA_20240620.F2SL514C

    DATA_20240620.F2CL514R
    DATA_20240620.F2CL514G
    
    DATA_20240620.F9ML514R
    DATA_20240620.F9ML514G
    
    DATA_20240620.F9CL514R
    DATA_20240620.F9CL514G
    };

WL = 514.5;

for i=1:length(FilmSamples514)
    current = FilmSamples514{i};  % Access the cell array element once
    current = clip_spectrum(current, 10,10);
    current = remove_inclination(current, WL);
    current = correct_instrument_response(current, WL);
    current = remove_bg_poly(current);
    FilmSamples514{i} = current;  % Save the result back to the cell array
end

% plotRaman(FilmSamples514, 0.0)
F1S514 = concatenateSpectra({FilmSamples514{1:2}}, 'Film T1S (SWCNT)')
F1C514 = concatenateSpectra({FilmSamples514{3:4}}, 'Film T1S (Converted to DWCNT)')
F2S514 = concatenateSpectra({FilmSamples514{5:7}}, 'Film T2S (SWCNT)')
F2C514 = concatenateSpectra({FilmSamples514{8:9}}, 'Film T2S (Converted to DWCNT)')
F9M514 = concatenateSpectra({FilmSamples514{10:11}}, 'Film T9M (SWCNT)')
F9C514 = concatenateSpectra({FilmSamples514{12:13}}, 'Film T9M (Converted to DWCNT)')

GroupedFS514 = {
                    F1S514
                    F1C514
                    F2S514
                    F2C514
                    F9M514
                    F9C514
                    }
                
GroupedFS514 = UsefulFunctions.NormalizeSample(GroupedFS514,0, 2000)
% plotRaman(GroupedFS514, 1)



%% DATA 20240628 FILM SAMPLES RBMs,G and D Band, Carbine SPEC AT 561nm


FilmSamples561={
    DATA_20240628.F1SL561R
    DATA_20240628.F1SL561G
    DATA_20240628.F1SL561C

    DATA_20240628.F1CL561R
    DATA_20240628.F1CL561G
    DATA_20240628.F1CL561C
    
    DATA_20240628.F2SL561R
    DATA_20240628.F2SL561G
    DATA_20240628.F2SL561C
    
    DATA_20240628.F2CL561R
    DATA_20240628.F2CL561G
    DATA_20240628.F2CL561C
    
    DATA_20240628.F9ML561R
    DATA_20240628.F9ML561G
    DATA_20240628.F9ML561C
    
    DATA_20240628.F9CL561R
    DATA_20240628.F9CL561G
    DATA_20240628.F9CL561C
    };

WL = 561;

for i=1:length(FilmSamples561)
    current = FilmSamples561{i};  % Access the cell array element once
    current = clip_spectrum(current, 10,10);
    current = remove_inclination(current, WL);
    current = correct_instrument_response(current, WL);
    current = remove_bg_poly(current);
    FilmSamples561{i} = current;  % Save the result back to the cell array
end

% plotRaman(FilmSamples561, 0.0)

F1S561 = concatenateSpectra({FilmSamples561{1:3}}, 'Film T1S (SWCNT)')
F1C561 = concatenateSpectra({FilmSamples561{4:6}}, 'Film T1S (Converted to DWCNT)')
F2S561 = concatenateSpectra({FilmSamples561{7:9}}, 'Film T2S (SWCNT)')
F2C561 = concatenateSpectra({FilmSamples561{10:12}}, 'Film T2S (Converted to DWCNT)')
F9M561 = concatenateSpectra({FilmSamples561{13:15}}, 'Film T9M (SWCNT)')
F9C561 = concatenateSpectra({FilmSamples561{16:18}}, 'Film T9M (Converted to DWCNT)')

GroupedFS561 = {
                    F1S561
                    F1C561
                    F2S561
                    F2C561
                    F9M561
                    F9C561
                    }
                
GroupedFS561 = UsefulFunctions.NormalizeSample(GroupedFS561,0, 2000)

% plotRaman(GroupedFS561, 1)


%% Trying to put everything together - same plot, different scales.%%%%%% 


%% Functions written for better corrections

function DS = clipRangeEdges(DS, min_value, max_value)
    % Find the indices where X is within the specified range
    idx = DS.X >= min_value & DS.X <= max_value;
    
    % Filter X and Y based on these indices
    DS.X = DS.X(idx);
    DS.Y = DS.Y(idx);
    DS.P = DS.P(idx);
end

function DS = remove_gaussian_profile(DS, sigma)
    % Remove Gaussian profile from Raman spectrum centered at zero
    % DS: structure with fields X (Raman shift) and Y (intensity)
    % sigma: standard deviation for the Gaussian profile

    % Extract the X and Y data
    X = DS.X;  % Raman shift (assumed centered at zero)
    Y = DS.Y;  % Intensity values

    % Define the Gaussian function centered at zero
    gaussian = @(x, A, sigma) A * exp(-((x).^2) / (2 * sigma^2));

    % Initial guess for parameters [Amplitude]
    A_initial = max(Y);  % Initial amplitude (could be the max value)

    % Fit the Gaussian to the data near the center
    idx_fit = (X > -20) & (X < 20);  % Fit range around zero (adjust as needed)

    % Fit the Gaussian using lsqcurvefit
    options = optimset('Display', 'off');  % Suppress output
    initial_guess = [A_initial];

    % Use lsqcurvefit to find the best-fit parameters
    [fit_params, ~] = lsqcurvefit(@(A, x) gaussian(x, A, sigma), initial_guess, X(idx_fit), Y(idx_fit));

    % Generate the fitted Gaussian profile across the entire spectrum
    fitted_gaussian = gaussian(X, fit_params(1), sigma);

    % Subtract the fitted Gaussian from the original spectrum
    Y_corrected = Y - fitted_gaussian;

    % Update the structure with the corrected Y values
    DS.Y = Y_corrected;

    % Optionally, return the fitted Gaussian parameters
    % if nargout > 1
    %     varargout{1} = fit_params;
    % end
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