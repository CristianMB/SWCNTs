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
path_KITDial = [rootpath,'20240611\KIT_Samples.csv'];
path_KITConverted = [rootpath,'20240613\KIT_ConvertedSamples.csv'];
path_KITFilms = [rootpath,'20240614\KIT_FilmConvertedSamples.csv'];

%Select the paths of interest
paths = {
        path_KITConverted
        path_KIT
        path_REF
        path_KITDial
        path_KITFilms
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
       

DATA_20240614.Baseline.N='Baseline';
DATA_20240614.FilmT1S.N='T1SFilm';
DATA_20240614.FilmT1SConverted.N='T1SFilmConverted';
DATA_20240614.FilmT2S.N='T2SFilm';
DATA_20240614.FilmT2SConverted.N='T2SFilmConverted';
DATA_20240614.FilmT9M.N='T9MFilm';
DATA_20240614.FilmT9MConverted.N='T9MFilmConverted';




%DATA_20240611 - Samples from KIT in DOC/D2O after exchange in dialysis

KITConverted = {
    DATA_20240613.T1S
    DATA_20240613.T4S
    DATA_20240613.T4P
    };                 

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
   
KITFilms = {
%     DATA_20240614.FilmT1S
%     DATA_20240614.FilmT1SConverted
%     DATA_20240614.FilmT2S
%     DATA_20240614.FilmT2SConverted
    DATA_20240614.FilmT9M
    DATA_20240614.FilmT9MConverted
            };
         
        
%%%--------BACKGROUND CORRECTION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LS2= 1000-20;  
US2= 1000+20;

% KITDial = SubtractInverseBG(KITDial,[418 610 813]);
KITDial = NormalizeSample(KITDial,900, 1100);
% KITDial = NormalizeSample(KITDial,1450, 1500);
plotAbsorption(KITDial,0.7)
% 
% KIT = SubtractInverseBG(KIT,[418 610 813]);
% KIT = NormalizeSample(KIT,620, 820);

% plotAbsorption(KIT,0.0)
 
% for i=1:length(KIT)
%     current = KIT{i} 
%     current.Y = current.Y - KIT{4}.Y
%     AKIT{i} = current
% end


% KITFilms = SubtractInverseBG(KITFilms,[418 610 1260]);
% KITFilms = SubtractInverseBG(KITFilms,[688 1380 2480]);
% KITFilms = NormalizeSample(KITFilms,900, 1100);
% plotAbsorption(KITFilms,0)

% KITConverted = FitSamples(KITConverted,peaks);
% for i=1:length(KITConverted)
%     plotRamanFit(KITConverted{i})
% end

plotAbsorption(KITConverted, 0.7)

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



