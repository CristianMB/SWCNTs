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


%% %--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %--------SAMPLE COMPARISION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

LG = 1425;
HG = 1645;
NP = 1590;
Tol = 20;

GBand650 = {DATA_20240426.BAL650G, DATA_20240426.BBL650G};
% 
% GBand650 = SubstractLinearBG(GBand650, LG, HG);
% GBand650 = ClipSamples(GBand650,10,0); 
% GBand650 = NormalizeSample(GBand650,NP-Tol, NP+Tol); 
% peaks = [1580 1600]; 
% GBand650 = FitSamples(GBand650,peaks); 
% 
% for i=1:length(GBand650)
%    plotRamanFit(GBand650{i})
% end



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
FullSpectra650 = ClipSamples(FullSpectra650, 40, 40);
FullSpectra650 = RemoveInclination(FullSpectra650, WL);
FullSpectra650 = InstrumentCorrection(FullSpectra650, WL);
FullSpectra650 = RemoveBackground(FullSpectra650);


WL = 570;
FullSpectra570 = ClipSamples(FullSpectra570, 40, 40);
FullSpectra570 = RemoveInclination(FullSpectra570, WL);
FullSpectra570 = InstrumentCorrection(FullSpectra570, WL);
FullSpectra570 = RemoveBackground(FullSpectra570);

plotRaman(FullSpectra570,0)

SampleA.X = [FullSpectra650{9}.X; FullSpectra650{10}.X]
SampleA.Y = [FullSpectra650{9}.Y; FullSpectra650{10}.Y]
SampleA.N = 'S-SWCNT'

SampleB.X = [FullSpectra650{11}.X; FullSpectra650{12}.X]
SampleB.Y = [FullSpectra650{11}.Y; FullSpectra650{12}.Y]
SampleB.N = 'Covnerted S-SWCNT'

[SampleA.X, sort_idx] = sort(SampleA.X);
SampleA.Y = SampleA.Y(sort_idx);


plotRaman({SampleA},0)



%plotRaman(GBandCarbide650, 0)
%plotRaman(FullSpectra570, 0)
%plotRaman(FullSpectra650, 0)

%plotRaman([GBand650,DBand650], 0)
%plotRaman([RBMBs650, RBMAs650], 0)

%GDBand520 = GDBandPeaksCalculation(GDBand520, 1580,1600,1560,1570,1300,1400);
%exportGDBandPeaks(GDBand520, 'dsfsd.csv');