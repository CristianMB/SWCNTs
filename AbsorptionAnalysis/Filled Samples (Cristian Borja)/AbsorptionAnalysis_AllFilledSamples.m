clc;
clear;
addpath('X:\SWCNTs');
import UsefulFunctions.*;
rootpath = 'X:\Measurements Data\Absorption\';

%All paths as default

References = [rootpath,'References.csv'];
SF_CF = [rootpath,'20231206\alldata20231206.csv'];
CB_CF_S234 = [rootpath,'20240212\AllDataCentrifuge.csv'];
CB_CF_S567 = [rootpath,'20240307\DGUS5S6S7.csv'];
CB_CF_SF6_S5789 = [rootpath,'20241022\CF_FilledSamples.csv'];
KIT_FILM_S234567 = [rootpath,'20241003\FilledSamplesFilms.csv'];
CB_DGU_S234 = [rootpath,'20240216\DGUC.csv'];
CB_DGU_S567 = [rootpath,'20240308\DialS5S6S7.csv'];
CB_DGU_AP = [rootpath,'20241030\DGU_APCNTs.csv'];

%Select the paths of interest
paths = {   References
            SF_CF
            CB_CF_S234
            CB_CF_S567
            CB_CF_SF6_S5789
            KIT_FILM_S234567
            CB_DGU_S234
            CB_DGU_S567
            CB_DGU_AP
        };

%Read and structure data from the paths

ReadAbsorptionFromPaths(paths);


%% %--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


DATA_20231206.SFF4dil.N='SF Methanol@SWCNT A';
DATA_20231206.SFF4dil2.N='SF Methanol@SWCNT B';
DATA_20231206.SFF2dil.N='SF PCE@SWCNT';
DATA_20231206.SFF2Bdil.N='SF PCE@SWCNT dil1';
DATA_20231206.SFF5dil.N='SF SFF5 Water@SWCNT';
DATA_20231206.SFF3dil.N='SF TCE@SWCNT LP Long Rinsing';
DATA_20231206.SFF3Bdil.N='SF TCE@SWCNT LP Short Rinsing';
DATA_20231206.SFF3_3dil.N='SF TCE@SWCNT RF';
DATA_20231206.SFF1Bdil.N='SF TEMED@SWCNT';
DATA_20231206.SFF6dil.N='SF TTF@SWCNT';
DATA_20240212.S2.N='CB PCE@SWCNT 4hCF';
DATA_20240212.S3.N='CB TCE@SWCNT 4hCF';
DATA_20240212.S4.N='CB TEMED@SWCNT 4hCF';
DATA_20240216.S2DGUC.N='CB PCE@SWCNT Dial. DGU C (Filled';
DATA_20240216.S3DGUC.N='CB TCE@SWCNT Dial. DGU C (Filled)';
DATA_20240216.S4DDGUC.N='CB TEMED@SWCNT Dial. DGU C (Filled)';
DATA_20240216.S4DDGUC2.N='CB TEMED@SWCNT Dial. DGU C (Filled)';
DATA_20240216.S4LDGUC.N='CB TEMED@SWCNT Dial. DGU C (Filled)';
DATA_20240216.Baseline.N='DOCD2O Baseline';
DATA_20240307.S7CF.N='CB Dodecane@SWCNT 4hCF';
DATA_20240307.S7DGUC.N='CB Dodecane@SWCNT DGU C (Filled)';
DATA_20240307.S5DGUB.N='CB Empty@SWCNT DGU B S5';
DATA_20240307.S6DGUB.N='CB Empty@SWCNT DGU B S6';
DATA_20240307.S7DGUB.N='CB Empty@SWCNT DGU B S7';
DATA_20240307.S6CF.N='CB Hexadecane@SWCNT 4hCF';
DATA_20240307.S6DGUC.N='CB Hexadecane@SWCNT DGU C (Filled)';
DATA_20240307.S5CF.N='CB TDAE@SWCNT 4hCF';
DATA_20240307.S5DGUC.N='CB TDAE@SWCNT DGU C (Filled)';
DATA_20240307.Baseline.N='DOCD2O Baseline';
DATA_20240308.S7DGUCDial.N='CB Dodecane@SWCNT Dial. DGU C (Filled)';
DATA_20240308.S6DGUCDial.N='CB Hexadecane@SWCNT Dial. DGU C (Filled)';
DATA_20240308.S5DGUCDial.N='CB TDAE@SWCNT Dial. DGU C (Filled)';
DATA_20241003.FS2.N='PCE@P2-ASWCNTs (Film)';
DATA_20241003.FS3.N='TCE@P2-ASWCNTs (Film)';
DATA_20241003.FS4.N='TEMED@P2-ASWCNTs (Film)';
DATA_20241003.FS5.N='TDAE@P2-ASWCNTs (Film)';
DATA_20241003.FS6.N='Hexadecane@P2-ASWCNTs (Film)';
DATA_20241003.FS7.N='Dodecane@P2-ASWCNTs (Film)';
DATA_20241003.Baseline.N='Sapphire Baseline';
DATA_20241022.Baseline.N='1%DOC/D2O';
DATA_20241022.S5CF.N='TDAE@P2-ASWCNTs (24hCF)';
DATA_20241022.S7CF.N='Dodecane@P2-ASWCNTs (24hCF)';
DATA_20241022.S8CF.N='PCE@P2-ASWCNTs (24hCF)';
DATA_20241022.S9CF.N='TEMED@P2-ASWCNTs (24hCF)';
DATA_20241022.SFF6CF.N='TTF@P2-ASWCNTs (24hCF)';
DATA_20241030.Baseline.N='1%DOC/D2O';
DATA_20241030.StockSolution.N='StockSolution';
DATA_20241030.DGU_Empty.N='Empty fraction';
DATA_20241030.DGU_Filled.N='Water filled fraction';
DATA_20241030.CF_APCNTs.N='Centrifuged AP CNTs';
DATA_20241030.DGUDial_Filled.N='Water filled fraction (after dialysis)';
DATA_20241030.DGUDial_Empty.N='Empty fraction (after dialysis)';


%% %--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% %--------BACKGROUND CORRECTION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %--------SPECTRA SELECTION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


CF = { 
        DATA_20231206.SFF4dil
        DATA_20231206.SFF2dil
        DATA_20231206.SFF5dil
        DATA_20231206.SFF3dil
        DATA_20231206.SFF3Bdil
        DATA_20231206.SFF3_3dil
        DATA_20231206.SFF1Bdil
        DATA_20231206.SFF6dil
%         DATA_20240212.S2
%         DATA_20240212.S3
%         DATA_20240212.S4
%         DATA_20240307.S7CF
%         DATA_20240307.S6CF
%         DATA_20240307.S5CF
%         DATA_20241003.FS2
%         DATA_20241003.FS3
%         DATA_20241003.FS4
%         DATA_20241003.FS5
%         DATA_20241003.FS6
%         DATA_20241003.FS7
%         DATA_20241022.S5CF
%         DATA_20241022.S7CF
%         DATA_20241022.S8CF
%         DATA_20241022.S9CF
%         DATA_20241022.SFF6CF
%         DATA_20241030.CF_APCNTs
        };
        

CF = NormalizeSample(CF,800, 1100);
plotAbsorption(CF, 0.0)




DGU = { 
        DATA_20240216.S2DGUC
        DATA_20240216.S3DGUC
        DATA_20240216.S4DDGUC
        DATA_20240308.S5DGUCDial        
        DATA_20240308.S6DGUCDial
        DATA_20240308.S7DGUCDial
        DATA_20241030.DGUDial_Filled
        DATA_References.WaterFilled
        DATA_20241030.DGUDial_Empty
        DATA_References.empty_P2_dial_0930
        };
        

% DGU = NormalizeSample(DGU,1500, 1850);
% plotAbsorption(DGU, 0.0)