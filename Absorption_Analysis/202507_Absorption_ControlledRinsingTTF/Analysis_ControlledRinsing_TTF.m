clc;
clear;
import UsefulFunctions.*;
rootpath = 'X:\Measurements Data\Absorption\';

%All paths as default

TTF_MeOH = [rootpath,'20250509\TTF_MeOH_Solubility.csv'];

S11_TTF_Rinsing1 = [rootpath,'20241211\TTF_MeOH_S11.csv'];
S11_TTF_Rinsing2 = [rootpath,'20241128\TTF_MeOH_S11.csv'];
S11_TTF_Rinsing3 = [rootpath,'20241115\Rinsing_S11_TTF.csv'];

S12_TTF_Rinsing = [rootpath,'20250128\Rinsing_TTF_MeOH_S12.csv'];
S13_TTF_Rinsing = [rootpath,'20250327\S13_TTF_Rinsings.csv'];
S16_TTF_Rinsing = [rootpath,'20250430\Rinsing_S16_TTF_MeOH.csv'];
S17_TTF_Rinsing = [rootpath,'20250514\S17_Rinsing_TTF_MeOH.csv'];
S20_TTF_Rinsing = [rootpath,'20250616\Rinsing_S20_TTF_MeOH.csv'];
S21_TTF_Rinsing = [rootpath,'20250716\Rinsing_S21_TTF_MeOH.csv'];

%Select the paths of interest
paths = {
    TTF_MeOH
    S11_TTF_Rinsing1
    S11_TTF_Rinsing2
    S11_TTF_Rinsing3

    S12_TTF_Rinsing
    S13_TTF_Rinsing
    S16_TTF_Rinsing
    S17_TTF_Rinsing
    S20_TTF_Rinsing
    S21_TTF_Rinsing
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

DATA_20250716.Baseline.N='Baseline Methanol';
DATA_20250716.R5_1.N='Rinsing 5 - TTF MeOH';
DATA_20250716.R4_1.N='Rinsing 4 - TTF MeOH';
DATA_20250716.R3_1.N='Rinsing 3 - TTF MeOH';
DATA_20250716.R2_1.N='Rinsing 2 - TTF MeOH';
DATA_20250716.R1_5.N='Rinsing 1 - TTF MeOH';
DATA_20250716.R1_30.N='Rinsing 1 - TTF MeOH';



%% Corrections

DATA_20250616.R1_30.Y = DATA_20250616.R1_30.Y*30
DATA_20250616.R2_2.Y = DATA_20250616.R2_2.Y*2

DATA_20250716.R1_30.Y = DATA_20250716.R1_30.Y*30
DATA_20250716.R1_5.Y = DATA_20250716.R1_5.Y*5

%% TTF in MeOH

TTFinMeOH = {
                DATA_20250509.TTF_MeOH_LC_5
            };
        
% plotAbsorption(TTFinMeOH, 0.00)
%% TTF filled by melt

S11_TTF = {
            DATA_20241211.TTF_S11_R9_1
            DATA_20241211.TTF_S11_R8_1
            DATA_20241211.TTF_S11_R7_1
            DATA_20241211.TTF_S11_R6_1
            
            DATA_20241128.TTF_S11_R6_1
            DATA_20241128.TTF_S11_R5_1
            
            DATA_20241115.TTF_S11_R5
            DATA_20241115.TTF_S11_R4
            
            DATA_20241115.TTF_S11_R3
%             DATA_20241115.TTF_S11_R3B
%             DATA_20241115.TTF_S11_R3C
%             DATA_20241115.TTF_S11_R3D
            
            DATA_20241115.TTF_S11_R2
            DATA_20241115.TTF_S11_R1       };

% S11_TTF = Normalize(S11_TTF, 280,340,'M')
% plotAbsorption(S11_TTF, 0.00)
%% TTF filled by reflux in meoh

S12_TTF = {
                DATA_20250128.TTF_S12_R1_1
                DATA_20250128.TTF_S12_R2_1
                DATA_20250128.TTF_S12_R3_1
                DATA_20250128.TTF_S12_R4_1
                DATA_20250128.TTF_S12_R5_1
                DATA_20250128.TTF_S12_R6_1
                DATA_20250128.TTF_S12_R6_1_remeasured
             };

% S12_TTF = Normalize(S11_TTF, 280,340,'M')
% plotAbsorption(S12_TTF, 0.00)

%% TTF filled by gas phase, not totally in vacuum

S13_TTF = {
            DATA_20250327.S13_TTF_MeOH_R1
            DATA_20250327.S13_TTF_MeOH_R2
            DATA_20250327.S13_TTF_MeOH_R3
            DATA_20250327.S13_TTF_MeOH_R4
            DATA_20250327.S13_TTF_MeOH_R5
            DATA_20250327.S13_TTF_MeOH_R5_remeasured
             };

% S13_TTF = Normalize(S13_TTF, 280,340,'M')
% plotAbsorption(S13_TTF, 0.00)


%% TTF filled by gas phase, not totally in vacuum

S16_TTF = {
            DATA_20250430.S16_TTF_MeOH_R1_8
            DATA_20250430.S16_TTF_MeOH_R2_8
            DATA_20250430.S16_TTF_MeOH_R3_8
            DATA_20250430.S16_TTF_MeOH_R4_8
            DATA_20250430.S16_TTF_MeOH_R5_8
            DATA_20250430.S16_TTF_MeOH_R6_4
            };

% S16_TTF = Normalize(S16_TTF, 280,340,'M')
% plotAbsorption(S16_TTF, 0.00)

%% TTF filled by gas phase 

S17_TTF = {
           DATA_20250514.R1_10
           DATA_20250514.R2_1
           DATA_20250514.R3_1
           DATA_20250514.R4_1
           DATA_20250514.R5_1
           DATA_20250514.R5_1_repeat
           DATA_20250514.R6_1
           DATA_20250514.R7_1
       };

% S17_TTF = Normalize(S17_TTF, 280,340,'M')
% plotAbsorption(S17_TTF, 0.00)


%% TTF filled by gas phase 

S20_TTF = {
%         DATA_20250616.R1_1
%         DATA_20250616.R1_30
        DATA_20250616.R2_1
%         DATA_20250616.R2_2
        DATA_20250616.R3_1
        DATA_20250616.R4_1
        DATA_20250616.R5_1
       };

% S20_TTF = Normalize(S20_TTF, 280,340,'M')
plotAbsorption(S20_TTF, 0.00)

%% TTF filled by melt 

S21_TTF = {
%             DATA_20250716.R1_5
%             DATA_20250716.R1_30
            DATA_20250716.R2_1
            DATA_20250716.R3_1
            DATA_20250716.R4_1
            DATA_20250716.R5_1         
       };

% S21_TTF = Normalize(S21_TTF, 280,340,'M')
% plotAbsorption(S21_TTF, 0.00)

