clc;
clear;
import UsefulFunctions.*;
rootpath = 'X:\Measurements Data\Absorption\';
addpath('X:\SWCNTs\');
%All paths as default

Rinsing_S18_TCNQ_DCM = [rootpath,'20250714\Rinsing_S18_TCNQ_DCM.csv'];
Rinsing_S23_TCNQ_DCM = [rootpath,'20250724\Rinsing_S23_TCNQ_DCM.csv'];



%Select the paths of interest
paths = {
            Rinsing_S18_TCNQ_DCM
            Rinsing_S23_TCNQ_DCM
};

%% Read and structure data from the paths
ReadAbsorptionFromPaths(paths);

%% Labeling


%% Corrections


DATA_20250724.R2_2.Y = DATA_20250724.R2_2.Y*2
DATA_20250724.R1_30.Y = DATA_20250724.R1_30.Y*30
DATA_20250724.R1_50.Y = DATA_20250724.R1_50.Y*50



%% TCNQ Gas Melted - S18

% S18_TCNQ = {
    
%        };

% S18_TCNQ = Normalize(S18_TCNQ, 280,340,'M')
% plotAbsorption(S18_TCNQ, 0.00)

%% TCNQ Gas - S23

S23_TCNQ = {

%             DATA_20250724.R1_1
%             DATA_20250724.R1_30
            DATA_20250724.R1_50
            DATA_20250724.R2_1
%             DATA_20250724.R2_2
            DATA_20250724.R3_1
            DATA_20250724.R4_1
%             DATA_20250724.R5_1
            DATA_20250724.R5_1_again
            DATA_20250724.R6_1

                
       };

% S23_TCNQ = Normalize(S23_TCNQ, 280,340,'M')
plotAbsorption(S23_TCNQ, 0.00)

