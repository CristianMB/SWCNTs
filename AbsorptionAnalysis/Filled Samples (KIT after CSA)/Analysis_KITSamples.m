clc;
clear;
addpath('X:\SWCNTs');
% addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\');
import UsefulFunctions.*;

% rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements Data\Absorption\';
rootpath = 'X:\Measurements Data\Absorption\';

%All paths as default
%path_baselines = [rootpath,'References.csv'];
path_KIT = [rootpath,'20240531\KIT_Samples.csv'];
path_REF = [rootpath,'References.csv'];
path_KITDial = [rootpath,'20240611\KIT_Samples.csv'];
path_KIT_DialWastes = [rootpath,'20240730\DialysisKITCSA_D1&D2.csv'];
path_KIT_FinalDialWastes = [rootpath,'20240802\DialysisKITCSA.csv'];

%Select the paths of interest
paths = {
        path_KIT
        path_REF
        path_KITDial
        path_KIT_DialWastes
        path_KIT_FinalDialWastes
        };

%Read and structure data from the paths
ReadAbsorptionFromPaths(paths);


%%%--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Compare before and after H2O in D2O Correction

%SampleToCorrect = DATA_20240307.S7DGUC
%H2OD2OFactor = 0.7
%NicodenzFactor = 1.2


%%%%%%D2O Correction for Wastes:

DATA_20240730.S2.Y = DATA_20240730.S2.Y
DATA_20240730.S3.Y = DATA_20240730.S3.Y
DATA_20240730.S4.Y = DATA_20240730.S4.Y
DATA_20240730.S5.Y = DATA_20240730.S5.Y
DATA_20240730.S6.Y = DATA_20240730.S6.Y
DATA_20240730.S7.Y = DATA_20240730.S7.Y
DATA_20240730.D2S2.Y = DATA_20240730.D2S2.Y
DATA_20240730.D2S3.Y = DATA_20240730.D2S3.Y
DATA_20240730.D2S4.Y = DATA_20240730.D2S4.Y
DATA_20240730.D2S5.Y = DATA_20240730.D2S5.Y
DATA_20240730.D2S6.Y = DATA_20240730.D2S6.Y
DATA_20240730.D2S7.Y = DATA_20240730.D2S7.Y


%%%--------SAMPLE COMPARISION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DATA_20240531 - Original Samples from KIT in DOC/H2O
DATA_20240531.Baseline.N =  '1% DOC/D2O Baseline';
DATA_20240531.EA.N =  'Empty ArcSWCNTs';
DATA_20240531.WA.N =  'Water ArcSWCNTs';
DATA_20240531.T1Converted.N =  'T1 Converted';
DATA_20240531.T4Converted.N =  'T4 Converted';
DATA_20240531.T4ConvertedPellet.N =  'T4 Converted Pellet';
DATA_20240531.T4ConvertedPelletB.N =  'T4 Converted Pellet';
DATA_20240531.S2.N =  'PCE@P2-ASWCNTs KIT after CSA in 1%DOC/H2O';
DATA_20240531.S2B.N =  'PCE@P2-ASWCNTs KIT after CSA in 1%DOC/H2O';
DATA_20240531.S3.N =  'TCE@P2-ASWCNTs KIT after CSA in 1%DOC/H2O';
DATA_20240531.S4.N =  'TEMED@P2-ASWCNTs KIT after CSA in 1%DOC/H2O';
DATA_20240531.S5.N =  'TDAE@P2-ASWCNTs KIT after CSA in 1%DOC/H2O';
DATA_20240531.S6.N =  'Hexadecane@P2-ASWCNTs KIT after CSA in 1%DOC/H2O';
DATA_20240531.S7.N =  'Dodecane@P2-ASWCNTs KIT after CSA in 1%DOC/H2O';

DATA_20240611.S3Dial.N='TCE@P2-ASWCNTs (KIT after exchange)';
DATA_20240611.S4Dial.N='TEMED@P2-ASWCNTs (KIT after exchange)';
DATA_20240611.S5Dial.N='TDAE@P2-ASWCNTs (KIT after exchange)';
DATA_20240611.S6Dial.N='Hexadecane@P2-ASWCNTs (KIT after exchange)';
DATA_20240611.S7Dial.N='Dodecane@P2-ASWCNTs (KIT after exchange)';
DATA_20240611.S2Dial.N='PCE@P2-ASWCNTs (KIT after exchange)';
       



%DATA_20240611 - Samples from KIT in DOC/D2O after exchange in dialysis

    

KIT = {     
%     DATA_20240531.Baseline
%      DATA_20240531.EA
%      DATA_20240531.WA
     
%      DATA_20240531.T1Converted
%      DATA_20240531.T4Converted
% %      DATA_20240531.T4ConvertedPellet
%      DATA_20240531.T4ConvertedPelletB

%       DATA_References.H2OinD2O
%     DATA_20240531.S2  
      DATA_20240531.S2B
      DATA_20240531.S3
      DATA_20240531.S4
      DATA_20240531.S5
      DATA_20240531.S6
      DATA_20240531.S7
             };
%  
% DATA_20240611.S2Dial.Y = DATA_20240611.S2Dial.Y - 0.2*DATA_References.H2OinD2O.Y- 0.3*DATA_References.DOC.Y
% DATA_20240611.S3Dial.Y = DATA_20240611.S3Dial.Y - 0.2*DATA_References.H2OinD2O.Y
% DATA_20240611.S4Dial.Y = DATA_20240611.S4Dial.Y - 0.5*DATA_References.H2OinD2O.Y
% DATA_20240611.S5Dial.Y = DATA_20240611.S5Dial.Y - 0.5*DATA_References.H2OinD2O.Y
% DATA_20240611.S6Dial.Y = DATA_20240611.S6Dial.Y - 0.2*DATA_References.H2OinD2O.Y
% DATA_20240611.S7Dial.Y = DATA_20240611.S7Dial.Y - 0.2*DATA_References.H2OinD2O.Y    
         
KITDial = {
      DATA_20240611.S2Dial
      DATA_20240611.S3Dial
      DATA_20240611.S4Dial
      DATA_20240611.S5Dial
      DATA_20240611.S6Dial
      DATA_20240611.S7Dial
%        DATA_References.H2OinD2O
            };
   
KITDialWastes = {
    
                DATA_20240730.S2
                DATA_20240730.S3
                DATA_20240730.S4
                DATA_20240730.S5
                DATA_20240730.S6
                DATA_20240730.S7
%                 DATA_20240730.D2S2
%                 DATA_20240730.D2S3
%                 DATA_20240730.D2S4
%                 DATA_20240730.D2S5
%                 DATA_20240730.D2S6
%                 DATA_20240730.D2S7
                DATA_20240802.S2D8
                DATA_20240802.S3D8
                DATA_20240802.S4D8
                DATA_20240802.S5D8
                DATA_20240802.S6D8
                DATA_20240802.S7D8
                
                }

        
%%%--------BACKGROUND CORRECTION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LS2= 1000-20;  
US2= 1000+20;


% KITDialWastes = NormalizeSample(KITDialWastes,800, 1000);
plotAbsorption(KITDialWastes, 0.1)




% KITDial = SubtractInverseBG(KITDial,[418 610 813]);
KITDial = NormalizeSample(KITDial,900, 1100);
% KITDial = NormalizeSample(KITDial,1450, 1500);
% plotAbsorption(KITDial,0.7)
% 
% KIT = SubtractInverseBG(KIT,[418 610 813]);
% KIT = NormalizeSample(KIT,620, 820);

% plotAbsorption(KIT,0.0)
 
% for i=1:length(KIT)
%     current = KIT{i} 
%     current.Y = current.Y - KIT{4}.Y
%     AKIT{i} = current
% end




% peaks = [463 489 551 586];
% KIT = FitSamples(KIT,peaks);

% plotAbsorption(KIT,0.0)
% 
% for i=1:length(KIT)
%     plotRamanFit(KIT{i})
% end

% Blines = {
%     DATA_20240531.Baseline
%     DATA_20240308.Baseline
%             }

% 
% KITDial = SubtractInverseBG(KITDial,[418 610 813]);
% KITDial = NormalizeSample(KITDial,900, 1100);
% KITDial = NormalizeSample(KITDial,1450, 1500);
% plotAbsorption(KITDial,0.0)



