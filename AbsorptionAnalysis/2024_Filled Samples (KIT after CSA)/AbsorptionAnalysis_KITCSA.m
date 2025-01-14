clc;
clear;
addpath('X:\SWCNTs');
% addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\');
import UsefulFunctions.*;

% rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements Data\Absorption\';
rootpath = 'X:\Measurements Data\Absorption\';

%All paths as default
path_OriginalCSA = [rootpath,'20240531\KIT_Samples.csv'];
path_ShortDial = [rootpath,'20240611\KIT_Samples.csv'];
path_SD_DGUA = [rootpath,'20240723\DGU_KIT_Filled.csv'];
path_SD_DGUB = [rootpath,'20240724\DGU_KIT_Filled.csv'];
path_LongDial = [rootpath,'20240807\DGU_KIT_Filled_Dial.csv'];
path_LD_DGUA = [rootpath,'20240809\KIT_CSA_1WDial_DGU.csv'];
path_References = [rootpath,'References.csv'];




%Select the paths of interest
paths = {
        path_References
        path_OriginalCSA
        path_ShortDial
        path_SD_DGUA
        path_SD_DGUB
        path_LongDial
        path_LD_DGUA
        };


ReadAbsorptionFromPaths(paths);


%% %--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% Water Correction For Long Dial

DATA_20240807.CSA_Dial_S2.Y = DATA_20240807.CSA_Dial_S2.Y - 0.240*DATA_References.WaterInD2O.Y
DATA_20240807.CSA_Dial_S3.Y = DATA_20240807.CSA_Dial_S3.Y - 0.375*DATA_References.WaterInD2O.Y
DATA_20240807.CSA_Dial_S4.Y = DATA_20240807.CSA_Dial_S4.Y - 0.360*DATA_References.WaterInD2O.Y
DATA_20240807.CSA_Dial_S5.Y = DATA_20240807.CSA_Dial_S5.Y - 0.360*DATA_References.WaterInD2O.Y
DATA_20240807.CSA_Dial_S6.Y = DATA_20240807.CSA_Dial_S6.Y - 0.240*DATA_References.WaterInD2O.Y
DATA_20240807.CSA_Dial_S7.Y = DATA_20240807.CSA_Dial_S7.Y - 0.291*DATA_References.WaterInD2O.Y

%% %--------BACKGROUND CORRECTION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %--------SPECTRA SELECTION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        

Compare = {
            DATA_20240531.S2
            DATA_20240531.S2B
            DATA_20240531.S3
            DATA_20240531.S4
            DATA_20240531.S5
            DATA_20240531.S6
            DATA_20240531.S7
%             
            DATA_20240611.S2Dial
            DATA_20240611.S3Dial
            DATA_20240611.S4Dial
            DATA_20240611.S5Dial
            DATA_20240611.S6Dial
            DATA_20240611.S7Dial
            
%             DATA_20240723.S2A
%             DATA_20240723.S3A
%             DATA_20240723.S4A
%             DATA_20240723.S5A
%             DATA_20240724.S6A
%             DATA_20240724.S7A

%             DATA_20240723.S2B
%             DATA_20240723.S3B
%             DATA_20240723.S4B
%             DATA_20240723.S5B
%             DATA_20240724.S6B
%             DATA_20240724.S7B

%             DATA_20240807.CSA_Dial_S7
%             DATA_20240807.CSA_Dial_S6
%             DATA_20240807.CSA_Dial_S2
%             DATA_20240807.CSA_Dial_S3
%             DATA_20240807.CSA_Dial_S5
%             DATA_20240807.CSA_Dial_S4
% 
%             DATA_20240809.DGU_F_S3
%             DATA_20240809.DGU_F_S4
%             DATA_20240809.DGU_F_S5
%             DATA_20240809.DGU_B_S2
%             DATA_20240809.DGU_B_S3
%             DATA_20240809.DGU_B_S4
%             DATA_20240809.DGU_B_S5
            };
        
Compare = NormalizeSample(Compare,650, 750);
plotAbsorption(Compare, 0.05)
% xlim([192 2500])
% ylim([0 4])
