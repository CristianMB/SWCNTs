clc;
clear;
addpath('C:\Users\cborja\OneDrive - Universiteit Antwerpen\SWCNTs\');
% addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\');
import UsefulFunctions.*;

% rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements Data\Absorption\';
rootpath = 'C:\Users\cborja\OneDrive - Universiteit Antwerpen\Measurements Data\Absorption\';

%All paths as default
%path_baselines = [rootpath,'References.csv'];
path_KIT = [rootpath,'20240531\KIT_Samples.csv'];
path_REF = [rootpath,'References.csv'];

%Select the paths of interest
paths = {
        path_KIT
        path_REF
        };

%Read and structure data from the paths
ReadAbsorptionFromPaths(paths);


%%%--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Compare before and after H2O in D2O Correction
%SampleToCorrect = DATA_20240307.S7DGUC
%H2OD2OFactor = 0.7
%NicodenzFactor = 1.2


%%%--------SAMPLE COMPARISION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DATA_20240531.Baseline.N =  '1% DOC/D2O Baseline';
DATA_20240531.EA.N =  'Empty ArcSWCNTs';
DATA_20240531.WA.N =  'Water ArcSWCNTs';
DATA_20240531.T1Converted.N =  'T1 Converted';
DATA_20240531.T4Converted.N =  'T4 Converted';
DATA_20240531.T4ConvertedPellet.N =  'T4 Converted Pellet';
DATA_20240531.T4ConvertedPelletB.N =  'T4 Converted Pellet';
DATA_20240531.S2.N =  'PCE@P2-ASWCNTs';
DATA_20240531.S2B.N =  'PCE@P2-ASWCNTs';
DATA_20240531.S3.N =  'TCE@P2-ASWCNTs';
DATA_20240531.S4.N =  'TEMED@P2-ASWCNTs';
DATA_20240531.S5.N =  'TDAE@P2-ASWCNTs';
DATA_20240531.S6.N =  'Hexadecane@P2-ASWCNTs';
DATA_20240531.S7.N =  'Dodecane@P2-ASWCNTs';
       
          

KIT = {     
%     DATA_20240531.Baseline
%      DATA_20240531.EA
%      DATA_20240531.WA
     
     DATA_20240531.T1Converted
     DATA_20240531.T4Converted
%      DATA_20240531.T4ConvertedPellet
     DATA_20240531.T4ConvertedPelletB

%       DATA_References.H2OinD2O
%     DATA_20240531.S2  

%       DATA_20240531.S2B
%       DATA_20240531.S3
%       DATA_20240531.S4
%       DATA_20240531.S5
%       DATA_20240531.S6
%       DATA_20240531.S7
             };
         
%%%--------BACKGROUND CORRECTION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LS2= 1000-20;  
US2= 1000+20;

KIT = SubtractInverseBG(KIT,[418 610 813]);
KIT = NormalizeSample(KIT,620, 820);

plotAbsorption(KIT,1.0)
 
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




