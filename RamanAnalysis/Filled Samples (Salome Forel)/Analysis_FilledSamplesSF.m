clc;
clear;
addpath('C:\Users\cborja\OneDrive - Universiteit Antwerpen\SWCNTs\');
%addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\');
import UsefulFunctions.*;

%rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements Data\Raman\';
rootpath = 'C:\Users\cborja\OneDrive - Universiteit Antwerpen\Measurements Data\Raman\';

%All paths as default
path_20240111 = [rootpath,'20240111\'];

%Select the paths of interest

paths = {
    path_20240111
    };


ReadRamanFromPaths(paths);


%%%--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DATA_20240514.BBL570RB.N = 'S-SWCNTs Converted Film RBMs2 570nm';

%FILLED SAMPLES AFTER DIALISYS AT 650nm + WaterFilled/Empty References


DATA_20240111.FLATHD.N ='FlatField 514nm';
DATA_20240111.LL514HD.N ='Laser 514nm';
DATA_20240111.LL514.N ='Laser 514nm';
DATA_20240111.S240111J.N ='SF D2O@SWCNT 514nm';
DATA_20240111.S240111K.N ='SF TCE@SWCNT 514nm';
DATA_20240111.S240111KK.N ='SF TCE@SWCNT 514nm';
DATA_20240111.S240111L.N ='SF Methanol@SWCNT 514nm';
DATA_20240111.S240111M.N ='SF TCE@SWCNT 514nm';
DATA_20240111.S240111N.N ='SF TTF@SWCNT 514nm';
DATA_20240111.S240111O.N ='SF PCE@SWCNT 514nm';
DATA_20240111.S240111P.N ='SF PCE@SWCNT 514nm';
DATA_20240111.S240111Q.N ='SF PCE@SWCNT 514nm';
DATA_20240111.S240111R.N ='SF TEMED@SWCNT 514nm';
DATA_20240111.S240111S.N ='SC Empty@SWCNT 514nm';
DATA_20240111.S240111A.N ='SF D2O@SWCNT 514nm';
DATA_20240111.S240111B.N ='SF TCE@SWCNT 514nm';
DATA_20240111.S240111BB.N ='SF TCE@SWCNT 514nm';
DATA_20240111.S240111C.N ='SF Methanol@SWCNT 514nm';
DATA_20240111.S240111D.N ='SF TCE@SWCNT 514nm';
DATA_20240111.S240111E.N ='SF TTF@SWCNT 514nm';
DATA_20240111.S240111F.N ='SF PCE@SWCNT 514nm';
DATA_20240111.S240111G.N ='SF PCE@SWCNT 514nm';
DATA_20240111.S240111H.N ='SF PCE@SWCN 514nm';
DATA_20240111.S240111I.N ='SF TEMED@SWCNT 514nm';


%%%--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DATA_20240426.BAL650R.Y = DATA_20240426.BAL650R.Y/2;
%DATA_20240426.BAL650RB.Y = DATA_20240426.BAL650RB.Y/2;

%%%--------SAMPLE COMPARISION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
FilledSamples514 = {
                    DATA_20240111.S240111J
                    DATA_20240111.S240111K
                    DATA_20240111.S240111KK
                    DATA_20240111.S240111L
                    DATA_20240111.S240111M
                    DATA_20240111.S240111N
                    DATA_20240111.S240111O
                    DATA_20240111.S240111P
                    DATA_20240111.S240111Q
                    DATA_20240111.S240111R
                    DATA_20240111.S240111S
                    DATA_20240111.S240111A
                    DATA_20240111.S240111B
                    DATA_20240111.S240111BB
                    DATA_20240111.S240111C
                    DATA_20240111.S240111D
                    DATA_20240111.S240111E
                    DATA_20240111.S240111F
                    DATA_20240111.S240111G
                    DATA_20240111.S240111H
                    DATA_20240111.S240111I
                   };
               
GDBands514 =    {
            DATA_20240111.S240111A
            DATA_20240111.S240111B
            DATA_20240111.S240111BB
            DATA_20240111.S240111C
            DATA_20240111.S240111D
            DATA_20240111.S240111E
            DATA_20240111.S240111F
            DATA_20240111.S240111G
            DATA_20240111.S240111H
            DATA_20240111.S240111I
            };

        
RBMs514 =    {
            DATA_20240111.S240111J
            DATA_20240111.S240111K
            DATA_20240111.S240111KK
            DATA_20240111.S240111L
            DATA_20240111.S240111M
            DATA_20240111.S240111N
            DATA_20240111.S240111O
            DATA_20240111.S240111P
            DATA_20240111.S240111Q
            DATA_20240111.S240111R
            DATA_20240111.S240111S
            };
        
plotRaman(RBMs514, 0)

%% 
%Normalization
LG = 1400;
HG = 1470;
NP = 1565
Tol = 5
GDBands514 = SubstractLinearBG(GDBands514, LG, HG);
GDBands514 = NormalizeSample(GDBands514,NP-Tol, NP+Tol)       

%Normalization
LG = 130;
HG = 220;
NP = 175
Tol = 25
RBMs514 = SubstractLinearBG(RBMs514, LG, HG);
RBMs514 = NormalizeSample(RBMs514,NP-Tol, NP+Tol)       
 
%% 
%plotRaman(FilledSamples514, 0)
%plotRaman(GDBands514, 0)
%plotRaman(RBMs514, 0.2)

%GDBand520 = GDBandPeaksCalculation(GDBand520, 1580,1600,1560,1570,1300,1400);
%exportGDBandPeaks(GDBand520, 'dsfsd.csv');