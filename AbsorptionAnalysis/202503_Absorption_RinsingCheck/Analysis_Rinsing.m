clc;
clear;
import UsefulFunctions.*;
rootpath = 'X:\Measurements Data\Absorption\';
addpath('X:\SWCNTs\SpecialMatlabFunctions\DrosteEffect-BrewerMap-3.2.5.0')
addpath('X:\SWCNTs')


%All paths as default

S11_TTF_Melt_A = [rootpath,'20241115\Rinsing_S11_TTF.csv'];
S11_TTF_Melt_B = [rootpath,'20241128\TTF_MeOH_S11.csv'];
S11_TTF_Melt_C = [rootpath,'20241211\TTF_MeOH_S11.csv'];

S12_TTF_Refl_MeOH = [rootpath,'20250128\Rinsing_TTF_MeOH_S12.csv'];

S13_TTF_GasPhase = [rootpath,'20250327\S13_TTF_Rinsings.csv'];



Refs= [rootpath,'References.csv'];

%Select the paths of interest
paths = {
    S11_TTF_Melt_A
    S11_TTF_Melt_B
    S11_TTF_Melt_C
    
    S12_TTF_Refl_MeOH
    
    S13_TTF_GasPhase
            };
    

%Read and structure data from the paths
ReadAbsorptionFromPaths(paths);

% Labeling

DATA_20241115.TTF_S11_R1.N = 'S11-TTF@SWCNTs Rinsing #1 (+45ml MeOH) Day1';
DATA_20241115.TTF_S11_R2.N = 'S11-TTF@SWCNTs Rinsing #2 (+25ml MeOH)Day1';
DATA_20241115.TTF_S11_R3.N = 'S11-TTF@SWCNTs Rinsing #3 (+30ml MeOH)Day1';
DATA_20241115.TTF_S11_R4.N = 'S11-TTF@SWCNTs Rinsing #4 (+45ml MeOH)Day1';
DATA_20241115.TTF_S11_R5.N = 'S11-TTF@SWCNTs Rinsing #5 (+35ml MeOH)Day1';
DATA_20241128.TTF_S11_R6_1.N = 'S11-TTF@SWCNTs Rinsing #6 (+50ml MeOH)Day2';
DATA_20241211.TTF_S11_R7_1.N = 'S11-TTF@SWCNTs Rinsing #7 (+25ml MeOH)Day2';
DATA_20241211.TTF_S11_R8_1.N = 'S11-TTF@SWCNTs Rinsing #8 (+30ml MeOH)Day3';
DATA_20241211.TTF_S11_R9_1.N = 'S11-TTF@SWCNTs Rinsing #9 (+30ml MeOH)Day3';
DATA_20241211.TTF_S11_R10_1.N = 'S11-TTF@SWCNTs Rinsing #10 (+30ml MeOH)Day3';

DATA_20250128.TTF_S12_R1_1.N = 'S12-TTF@SWCNTs Rinsing #1 (+40ml MeOH)';
DATA_20250128.TTF_S12_R2_1.N = 'S12-TTF@SWCNTs Rinsing #2 (+40ml MeOH)';
DATA_20250128.TTF_S12_R3_1.N = 'S12-TTF@SWCNTs Rinsing #3 (+40ml MeOH)';
DATA_20250128.TTF_S12_R4_1.N = 'S12-TTF@SWCNTs Rinsing #4 (+40ml MeOH)';
DATA_20250128.TTF_S12_R5_1.N = 'S12-TTF@SWCNTs Rinsing #5 (+45ml MeOH)';
DATA_20250128.TTF_S12_R6_1.N = 'S12-TTF@SWCNTs Rinsing #6 (+40ml MeOH)';

DATA_20250327.S13_TTF_MeOH_R1.N = 'S13-TTF@SWCNTs Rinsing #1 (+40ml MeOH)';
DATA_20250327.S13_TTF_MeOH_R2.N = 'S13-TTF@SWCNTs Rinsing #2 (+40ml MeOH)';
DATA_20250327.S13_TTF_MeOH_R3.N = 'S13-TTF@SWCNTs Rinsing #3 (+40ml MeOH)';
DATA_20250327.S13_TTF_MeOH_R4.N = 'S13-TTF@SWCNTs Rinsing #4 (+40ml MeOH)';
DATA_20250327.S13_TTF_MeOH_R5.N = 'S13-TTF@SWCNTs Rinsing #5 (+40ml MeOH)';
            
% Plotting
S11 = {
            DATA_20241115.TTF_S11_R1
            DATA_20241115.TTF_S11_R2
            DATA_20241115.TTF_S11_R3
            DATA_20241115.TTF_S11_R4
            DATA_20241115.TTF_S11_R5
            DATA_20241128.TTF_S11_R6_1
            DATA_20241211.TTF_S11_R7_1
            DATA_20241211.TTF_S11_R8_1
            DATA_20241211.TTF_S11_R9_1
            DATA_20241211.TTF_S11_R10_1
            };
        
% plotAbsorption(S11, 0.0);

S12 = {
            DATA_20250128.TTF_S12_R1_1
            DATA_20250128.TTF_S12_R2_1
            DATA_20250128.TTF_S12_R3_1
            DATA_20250128.TTF_S12_R4_1
            DATA_20250128.TTF_S12_R5_1
            DATA_20250128.TTF_S12_R6_1
            };        
plotAbsorption(S12, 0.0);


S13 = {
            DATA_20250327.S13_TTF_MeOH_R1
            DATA_20250327.S13_TTF_MeOH_R2
            DATA_20250327.S13_TTF_MeOH_R3
            DATA_20250327.S13_TTF_MeOH_R4
            DATA_20250327.S13_TTF_MeOH_R5
            
            };
% plotAbsorption(S13, 0.0);
