clc;
clear;
import UsefulFunctions.*;
rootpath = 'X:\Measurements Data\Absorption\';

%All paths as default

TEMED_Rinsing = [rootpath,'20250611\Rinsing_S19.csv'];

%Select the paths of interest
paths = {
    TEMED_Rinsing
};

%% Read and structure data from the paths
ReadAbsorptionFromPaths(paths);

%% Labeling

DATA_20250611.Baseline.N='UltraPure Water H2O';
DATA_20250611.S19_R1_2.N='S19 TEMED@P2-SWCNTs Rinsing 1';
DATA_20250611.S19_R2_1.N='S19 TEMED@P2-SWCNTs Rinsing 2';
DATA_20250611.S19_R3_1.N='S19 TEMED@P2-SWCNTs Rinsing 3';
DATA_20250611.R19_R1_2.N='R19 TEMED + AP-SWCNTs Rinsing 1';
DATA_20250611.R19_R2_1.N='R19 TEMED + AP-SWCNTs Rinsing 2';
DATA_20250611.R19_R3_1.N='R19 TEMED + AP-SWCNTs Rinsing 3';


%% Corrections

        
%% Plotting

TEMED = {
%             DATA_20250611.R19_R1_1
            DATA_20250611.R19_R1_2
            DATA_20250611.R19_R2_1
            DATA_20250611.R19_R3_1
            
            DATA_20250611.S19_R1_2
            DATA_20250611.S19_R2_1
            DATA_20250611.S19_R3_1
       };

plotAbsorption(TEMED, 0.0)

   
