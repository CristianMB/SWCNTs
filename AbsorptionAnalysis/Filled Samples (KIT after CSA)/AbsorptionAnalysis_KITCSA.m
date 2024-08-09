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
path_LongDial = [rootpath,'20240807\DGU_KIT_Filled_Dial.csv'];
path_References = [rootpath,'References.csv'];




%Select the paths of interest
paths = {
        path_References
        path_OriginalCSA
        path_ShortDial
        path_LongDial
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
%             DATA_20240531.S2
%             DATA_20240531.S2B
%             DATA_20240531.S3
%             DATA_20240531.S4
%             DATA_20240531.S5
%             DATA_20240531.S6
%             DATA_20240531.S7
            
%             DATA_20240611.S2Dial
%             DATA_20240611.S3Dial
%             DATA_20240611.S4Dial
%             DATA_20240611.S5Dial
%             DATA_20240611.S6Dial
%             DATA_20240611.S7Dial  
            DATA_20240807.CSA_Dial_S2
            DATA_20240807.CSA_Dial_S3
            DATA_20240807.CSA_Dial_S4
            DATA_20240807.CSA_Dial_S5
            DATA_20240807.CSA_Dial_S6
            DATA_20240807.CSA_Dial_S7
            };
        
Compare = NormalizeSample(Compare,950, 1050);
plotAbsorption(Compare, 0.0)
% xlim([192 2500])
% ylim([0 4])
