clc;
clear;
import UsefulFunctions.*;
rootpath = 'X:\Measurements Data\Absorption\';

%All paths as default

TTF_Rinsing = [rootpath,'20250616\Rinsing_S20_TTF_MeOH.csv'];

%Select the paths of interest
paths = {
    TTF_Rinsing
};

%% Read and structure data from the paths
ReadAbsorptionFromPaths(paths);

%% Labeling

DATA_20250616.Baseline.N='Methanol Baseline';
DATA_20250616.R1_1.N='R1_1';
DATA_20250616.R1_30.N='R1_30';
DATA_20250616.R2_1.N='R2_1';
DATA_20250616.R2_2.N='R2_2';
DATA_20250616.R3_1.N='R3_1';
DATA_20250616.R4_1.N='R4_1';
DATA_20250616.R5_1.N='R5_1';



%% Corrections

DATA_20250616.R1_30.Y = DATA_20250616.R1_30.Y*30
DATA_20250616.R2_2.Y = DATA_20250616.R2_2.Y*2

%% Plotting

TTF = {
        DATA_20250616.R1_1
%         DATA_20250616.R1_30
        DATA_20250616.R2_1
%         DATA_20250616.R2_2
        DATA_20250616.R3_1
        DATA_20250616.R4_1
        DATA_20250616.R5_1
       };

plotAbsorption(TTF, 0.0)

   
