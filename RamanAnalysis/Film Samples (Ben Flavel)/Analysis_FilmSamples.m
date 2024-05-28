clc;
clear;
addpath('C:\Users\cborja\OneDrive - Universiteit Antwerpen\SWCNTs\');
%addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\');
import UsefulFunctions.*;

%rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements Data\Raman\';
rootpath = 'C:\Users\cborja\OneDrive - Universiteit Antwerpen\Measurements Data\Raman\';

%All paths as default
path_20240426 = [rootpath,'20240426\'];
path_20240514 = [rootpath,'20240514\'];
path_20240515 = [rootpath,'20240515\'];

%Select the paths of interest

paths = {
    path_20240426
    path_20240514
    path_20240515
    };


ReadRamanFromPaths(paths);


%%%--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DATA_20240426.BAL650D.N = 'S-SWCNTs Film DBand 650nm';
DATA_20240426.BAL650G.N = 'S-SWCNTs Film GBand 650nm';
DATA_20240426.BAL650R.N = 'S-SWCNTs Film RBMa 650nm';
DATA_20240426.BAL650RB.N = 'S-SWCNTs Film RBMb 650nm';
DATA_20240426.BBL650G.N = 'S-SWCNTs Converted Film GBand 650nm';
DATA_20240426.BBL650D.N = 'S-SWCNTs Converted Film DBand 650nm';
DATA_20240426.BBL650RA.N = 'S-SWCNTs Converted Film RBMa 650nm';
DATA_20240426.BBL650RB.N = 'S-SWCNTs Converted Film RBMb 650nm';

DATA_20240514.BBL650C1.N = 'S-SWCNTs Converted Film Carbide 650nm';
DATA_20240514.BBL650C2.N = 'S-SWCNTs Converted Film GBand 650nm';  
DATA_20240514.BAL650C1.N = 'S-SWCNTs Film Carbide 650nm';
DATA_20240514.BAL650C2.N = 'S-SWCNTs Film GBand 650nm';

DATA_20240515.BAL570C.N = 'S-SWCNTs Film Carbide 570nm';
DATA_20240515.BAL570D.N = 'S-SWCNTs Film DBand 570nm';
DATA_20240515.BAL570G.N = 'S-SWCNTs Film GBand 570nm';
DATA_20240515.BAL570RA.N = 'S-SWCNTs Film RBMs1 570nm';
DATA_20240515.BAL570RB.N = 'S-SWCNTs Film RBMs2 570nm';
DATA_20240515.BBL570C.N = 'S-SWCNTs Converted Film Carbide 570nm';
DATA_20240515.BBL570D.N = 'S-SWCNTs Converted Film DBand 570nm';
DATA_20240515.BBL570G.N = 'S-SWCNTs Converted Film GBand 570nm';
DATA_20240515.BBL570RA.N = 'S-SWCNTs Converted Film RBMs1 570nm';
DATA_20240515.BBL570RB.N = 'S-SWCNTs Converted Film RBMs2 570nm';


%%%--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DATA_20240426.BAL650R.Y = DATA_20240426.BAL650R.Y/2;
%DATA_20240426.BAL650RB.Y = DATA_20240426.BAL650RB.Y/2;

%%%--------SAMPLE COMPARISION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

% LG = 1425;
% HG = 1645;
% GBand650 = {DATA_20240426.BAL650G, DATA_20240426.BBL650G}
% %GBand650 = SubstractLinearBG(GBand650, LG, HG);
% %GBand650 = NormalizeSample(GBand650,1435, 1460)
% 
% LG = 1220;
% HG = 1430;
% DBand650 = {DATA_20240426.BAL650D, DATA_20240426.BBL650D}
% DBand650 = SubstractLinearBG(DBand650, LG, HG);
% DBand650 = NormalizeSample(DBand650,1435, 1460)
% 
% RBMBs650 = {DATA_20240426.BAL650RB, DATA_20240426.BBL650RB}
% RBMBs650 = NormalizeSample(RBMBs650,480, 500)
% RBMBs650 = SubstractLinearBG(RBMBs650,340, 520)
% 
% RBMAs650 = {DATA_20240426.BAL650R, DATA_20240426.BBL650RA}
% RBMAs650 = NormalizeSample(RBMAs650,330, 350)
% RBMAs650 = SubstractLinearBG(RBMAs650,340, 520)
% 


%GDBand520 = NormalizeSample(GDBand520,LG, HG);

%GBand520 = SubstractLinearBG(GBand520, 1540, 1620);
%GBand520 = NormalizeSample(GBand520,LG, HG);

%GDBand514 = SubstractLinearBG(GDBand514,1400, 1500);
%GDBand514 = NormalizeSample(GDBand514,LG,HG);


%LR = 150;
%HR = 160;

%RBM514 = SubstractLinearBG(RBM514, 130, 210);
%RBM514 = NormalizeSample(RBM514,LR, HR);

%RBM520 = SubstractLinearBG(RBM520, 130, 210);
%RBM520 = NormalizeSample(RBM520,LR, HR);

%RBM680 = SubstractLinearBG(RBM680,147, 200);
%RBM680 = NormalizeSample(RBM680,LR, HR);




FullSpectra650 = {  DATA_20240426.BAL650D,
                    DATA_20240426.BAL650G,
                    DATA_20240426.BAL650R,
                    DATA_20240426.BAL650RB,
                    DATA_20240426.BBL650G,
                    DATA_20240426.BBL650D,
                    DATA_20240426.BBL650RA,
                    DATA_20240426.BBL650RB
                };

GBandCarbide650 = {
       DATA_20240514.BAL650C1
       DATA_20240514.BAL650C2
       DATA_20240514.BBL650C1
       DATA_20240514.BBL650C2  
       };

FullSpectra570 = {
                  DATA_20240515.BAL570C
                  DATA_20240515.BAL570D
                  DATA_20240515.BAL570G
                  DATA_20240515.BAL570RA
                  DATA_20240515.BAL570RB
                  DATA_20240515.BBL570C
                  DATA_20240515.BBL570D
                  DATA_20240515.BBL570G
                  DATA_20240515.BBL570RA
                  DATA_20240515.BBL570RB
                  };

              
WL = 650;
for i=1:length(GBandCarbide650)
    current = GBandCarbide650{i};  % Access the cell array element once
    current = clip_spectrum(current, 40, 40);
    current = remove_inclination(current, WL);
    current = correct_instrument_response(current, WL);
    current = remove_bg_poly(current);
    GBandCarbide650{i} = current;  % Save the result back to the cell array
end        

WL = 650;
for i=1:length(FullSpectra650)
    current = FullSpectra650{i};  % Access the cell array element once
    current = clip_spectrum(current, 40, 40);
    current = remove_inclination(current, WL);
    current = correct_instrument_response(current, WL);
    current = remove_bg_poly(current);
    FullSpectra650{i} = current;  % Save the result back to the cell array
end        


WL = 570;
for i=1:length(FullSpectra570)
    current = FullSpectra570{i};  % Access the cell array element once
    current = clip_spectrum(current, 40, 40);
    current = remove_inclination(current, WL);
    current = correct_instrument_response(current, WL);
    current = remove_bg_poly(current);
    FullSpectra570{i} = current;  % Save the result back to the cell array
end      



DATA_20240515.FullCarbide = mergeStructuresRaman([DATA_20240514.BAL650C1, DATA_20240514.BAL650C2])

plotRaman({DATA_20240515.FullCarbide}, 0)

%plotRaman(GBandCarbide650, 0)
%plotRaman(FullSpectra570, 0)
%plotRaman(FullSpectra650, 0)

%plotRaman([GBand650,DBand650], 0)
%plotRaman([RBMBs650, RBMAs650], 0)

%GDBand520 = GDBandPeaksCalculation(GDBand520, 1580,1600,1560,1570,1300,1400);
%exportGDBandPeaks(GDBand520, 'dsfsd.csv');