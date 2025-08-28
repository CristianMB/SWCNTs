clc;
clear;
addpath('X:\SWCNTs');
import UsefulFunctions.*;
rootpath = 'X:\Measurements Data\Absorption\';

%All paths as default


P1 = [rootpath,'20241211\TTF_MeOH_S11.csv'];
P2 = [rootpath,'20241128\TTF_MeOH_S11.csv'];
P3 = [rootpath,'20241115\Rinsing_S11_TTF.csv'];

R1 = [rootpath,'20241126\TMG_MeOH_RinsingS10.csv'];

%Select the paths of interest
paths = {  P1
        P2
        P3
        R1
        };

%Read and structure data from the paths

ReadAbsorptionFromPaths(paths);


%% %--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% %--------BACKGROUND CORRECTION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %--------SPECTRA SELECTION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Rinsings

S11Rinsings = {
            DATA_20241115.TTF_S11_R1
            DATA_20241115.TTF_S11_R2
            DATA_20241115.TTF_S11_R3
            DATA_20241115.TTF_S11_R4
            DATA_20241128.TTF_S11_R5_1
            DATA_20241128.TTF_S11_R6_1
            DATA_20241211.TTF_S11_R6_1
            DATA_20241211.TTF_S11_R7_1
            DATA_20241211.TTF_S11_R8_1
            DATA_20241211.TTF_S11_R9_1
            DATA_20241211.TTF_S11_R10_1
%             DATA_20241211.TTF_S11_R10_pip
                }
% plotAbsorptionOrdered(S11Rinsings, 0)

% S11Rinsings=UsefulFunctions.NormalizeSample(S11Rinsings, 299, 301)
% plotAbsorption(S11Rinsings, 0)


S10Rinsings = {
            DATA_20241126.TMG_S10_R1_8
            DATA_20241126.TMG_S10_R2_8
            DATA_20241126.TMG_S10_R3_8
            DATA_20241126.TMG_S10_R4_8
                }

% S11Rinsings=UsefulFunctions.NormalizeSample(S11Rinsings, 299, 301)
% plotAbsorption(S10Rinsings, 0)




M2411 = {
        DATA_20241211.Baseline
        DATA_20241211.Methanol
        DATA_20241211.Methanol_nocap
        DATA_20241211.Methanol_nocap2
        DATA_20241211.Methanol_withCap
                DATA_20241211.Methanol_withCap1
    
                }
plotAbsorption(M2411, 0)

% S11Rinsings=UsefulFunctions.NormalizeSample(S11Rinsings, 299, 301)
% plotAbsorption(S11Rinsings, 0)