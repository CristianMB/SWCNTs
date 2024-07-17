clc;
clear;
addpath('X:\SWCNTs');
%addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\');
import UsefulFunctions.*;

%rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements Data\Absorption\';
rootpath = 'X:\Measurements Data\Absorption\';

%All paths as default
%path_baselines = [rootpath,'References.csv'];
path_20240220 = [rootpath,'20240220\RinsingS6S7.csv'];


%Select the paths of interest
paths = {
        path_20240220,
        };
    

%Read and structure data from the paths
ReadAbsorptionFromPaths(paths);


%%%--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Compare before and after H2O in D2O Correction
%SampleToCorrect = DATA_20240307.S7DGUC
%H2OD2OFactor = 0.7
%NicodenzFactor = 1.2

%Corrected = SampleToCorrect;
%Corrected.N = 'Corrected';
%Corrected.A = SampleToCorrect.A - DATA_References.H2OinD2O.A*H2OD2OFactor - DATA_References.Nicodenz.A*NicodenzFactor
%plotSampleList({Corrected, SampleToCorrect},0.0)
%plotSampleList({DATA_References.H2OinD2O, DATA_References.Nicodenz},0.0)

%-----------------------------------------------

%%%--------SAMPLE COMPARISION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        

          
Samples = {
               DATA_20240220.Baseline
                DATA_20240220.S6R1
                DATA_20240220.S6R2
                DATA_20240220.S6R3
            };
          
       
        
%%%--------BACKGROUND CORRECTION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Correction is made based on local minima (750, 860) and (1100, 1400)
% DialSamples = SubstractAbsBG(DialSamples,750, 850, 1150, 1250);

%%%--------NORMALIZATION TO S22 PEAK--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Normalization using maximum value of S22 transition, can also use Integral
% LS2= 900;  
% US2= 1100;
% %600a800
% 
% DialSamples = NormalizeSample(DialSamples,LS2, US2);
% BenSamp = NormalizeSample(BenSamp,LS2, US2);
% plotAbsorption(BenSamp,0.0)


%%%--------PEAK CALCULATION AND EXPORT--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Find peaks in these ranges
% %S11
% LS1= 1600;
% US1= 1850;
% 
% %Plasmon Right
% LPR= 260;
% UPR= 285;
% %Plasmon Left
% LPL= 220;
% UPL= 250;


plotAbsorption(Samples, 0.2)

%plotAbsorption(Alkanes,1.0)
%plotAbsorption(Dopants, 0.0)

%plotAbsorption(PCE, 0.0)
%plotAbsorption(TCE, 0.0)
%plotAbsorption(TDAE, 0.0)
%plotAbsorption(TEMED, 0.0)
%plotAbsorption(Hexadecane, 0.0)
%plotAbsorption(Dodecane, 0.0)

%plotAbsorption(DialSamples,1.0)


%plotAbsorption(DialSamples, 1.5)

%plotAbsorption(DialSamples, 0.0)

%plotAbsorption([DialSamples;DATABEN], 0.0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



