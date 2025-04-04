clc;
clear;
import UsefulFunctions.*;
rootpath = 'X:\Measurements Data\Absorption\';
addpath('X:\SWCNTs\SpecialMatlabFunctions\DrosteEffect-BrewerMap-3.2.5.0')
addpath('X:\SWCNTs')


%All paths as default

S6S7_Hexa_Dode = [rootpath,'20240220\RinsingS6S7.csv'];

S11_TTF_Melt_A = [rootpath,'20241115\Rinsing_S11_TTF.csv'];
S11_TTF_Melt_B = [rootpath,'20241128\TTF_MeOH_S11.csv'];
S11_TTF_Melt_C = [rootpath,'20241211\TTF_MeOH_S11.csv'];

S12_TTF_Refl_MeOH = [rootpath,'20250128\Rinsing_TTF_MeOH_S12.csv'];

S13_TTF_GasPhase = [rootpath,'20250327\S13_TTF_Rinsings.csv'];

S14_Hexa_LP = [rootpath,'20250331\Rinsing_S14_Hexadecane_EA.csv'];
S15_Dode_LP = [rootpath,'20250331\Rinsing_S15_Dodecane_EA.csv'];

Solvents = [rootpath,'20250331\PureSolvents.csv'];


Refs= [rootpath,'References.csv'];

%Select the paths of interest
paths = {
    S6S7_Hexa_Dode
    
    S11_TTF_Melt_A
    S11_TTF_Melt_B
    S11_TTF_Melt_C
    
    S12_TTF_Refl_MeOH
    
    S13_TTF_GasPhase
    
    S14_Hexa_LP
    
    S15_Dode_LP
    
    Solvents
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
            
DATA_20250331.S14_Hexa_EA_R1_1.N = 'S14-C_{16}H_{34}@SWCNTs Rinsing #1 (+30ml EA)';
DATA_20250331.S14_Hexa_EA_R2_1.N = 'S14-C_{16}H_{34}@SWCNTs Rinsing #2 (+45ml EA)';

DATA_20250331.S15_Dode_EA_R1_1.N = 'S15-C_{12}H_{26}@SWCNTs Rinsing #1 (+30ml EA)';
DATA_20250331.S15_Dode_EA_R2_1.N = 'S15-C_{12}H_{26}@SWCNTs Rinsing #2 (+45ml EA)';
            
% Plotting

S6 = {
            DATA_20240220.S6R1
            DATA_20240220.S6R2
            DATA_20240220.S6R3
            };
        
% plotAbsorption(S6, 0.0);

S7 = {
            DATA_20240220.S7R1
            DATA_20240220.S7R2
            DATA_20240220.S7R3
            };
        
% plotAbsorption(S7, 0.0);


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

S14 = {
            DATA_20250331.S14_Hexa_EA_R1_1
            DATA_20250331.S14_Hexa_EA_R2_1
            };
        
% plotAbsorption(S14, 0.0);

S15 = {     
            DATA_20250331.S15_Dode_EA_R1_1
            DATA_20250331.S15_Dode_EA_R2_1
            };
        
% plotAbsorption(S15, 0.0);

% DATA_20250331.DodeEA_10pc.Y = DATA_20250331.DodeEA_10pc.Y - DATA_20250331.EthylAcetate.Y

Solvents = {    
            DATA_20250331.EthylAcetate
            DATA_20250331.Dodecane
            DATA_20250331.DodeEA_10pc
            };
        

% Solvents = Normalize(Solvents, 0, 3000, 'M');
% plotAbsorption(Solvents, 0.00);
