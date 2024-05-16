clc;
clear;
addpath('C:\Users\cborja\OneDrive - Universiteit Antwerpen\SWCNTs\');
%addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\');
import UsefulFunctions.*;

%rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements Data\Raman\';
rootpath = 'C:\Users\cborja\OneDrive - Universiteit Antwerpen\Measurements Data\Raman\';

%All paths as default
path_20240514 = [rootpath,'20240514\'];

%Select the paths of interest

paths = {
    path_20240514
    };


ReadRamanFromPaths(paths);


%%%--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DATA_20240514.BBL570RB.N = 'S-SWCNTs Converted Film RBMs2 570nm';

%FILLED SAMPLES AFTER DIALISYS AT 650nm + WaterFilled/Empty References

DATA_20240514.EAL650D.N = 'EmptyArc@SWCNTs DBand 650nm';
DATA_20240514.EAL650G.N = 'EmptyArc@SWCNTs GBand 650nm';
DATA_20240514.EAL650R.N = 'EmptyArc@SWCNTs RBMs 650nm';
DATA_20240514.WAL650D.N = 'Water@ArcSWCNTs DBand 650nm';
DATA_20240514.WAL650G.N = 'Water@ArcSWCNTs GBand 650nm';
DATA_20240514.WAL650R.N = 'Water@ArcSWCNTs RBMs 650nm';
DATA_20240514.S2L650D.N = 'S2 PCE@SWCNTs DBand 650nm';
DATA_20240514.S2L650G.N = 'S2 PCE@SWCNTs GBand 650nm';
DATA_20240514.S2L650R.N = 'S2 PCE@SWCNTs RBMs 650nm';
DATA_20240514.S3L650D.N = 'S3 TCE@SWCNTs DBand 650nm';
DATA_20240514.S3L650G.N = 'S3 TCE@SWCNTs GBand 650nm';
DATA_20240514.S3L650R.N = 'S3 TCE@SWCNTs RBMs 650nm';
DATA_20240514.S4L650D.N = 'S4 TEMED@SWCNTs DBand 650nm';
DATA_20240514.S4L650G.N = 'S4 TEMED@SWCNTs GBand 650nm';
DATA_20240514.S4L650R.N = 'S4 TEMED@SWCNTs RBMs 650nm';
DATA_20240514.S5L650D.N = 'S5 TDAE@SWCNTs DBand 650nm';
DATA_20240514.S5L650G.N = 'S5 TDAE@SWCNTs GBand 650nm';
DATA_20240514.S5L650R.N = 'S5 TDAE@SWCNTs RBMs 650nm';
DATA_20240514.S6L650D.N = 'S6 Hexadecane@SWCNTs DBand 650nm';
DATA_20240514.S6L650G.N = 'S6 Hexadecane@SWCNTs GBand 650nm';
DATA_20240514.S6L650R.N = 'S6 Hexadecane@SWCNTs RBMs 650nm';
DATA_20240514.S7L650D.N = 'S7 Dodecane@SWCNTs DBand 650nm';
DATA_20240514.S7L650G.N = 'S7 Dodecane@SWCNTs GBand 650nm';
DATA_20240514.S7L650R.N = 'S7 Dodecane@SWCNTs RBMs 650nm';


%%%--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DATA_20240426.BAL650R.Y = DATA_20240426.BAL650R.Y/2;
%DATA_20240426.BAL650RB.Y = DATA_20240426.BAL650RB.Y/2;

%%%--------SAMPLE COMPARISION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

% LG = 1425;
% HG = 1645;
% GBand650 = {DATA_20240426.BAL650G, DATA_20240426.BBL650G}
% %GBand650 = SubstractLinearBG(GBand650, LG, HG);
% %GBand650 = NormalizeSample(GBand650,1435, 1460)


FilledSamples650 = {
                    DATA_20240514.EAL650D
                    DATA_20240514.EAL650G
                    DATA_20240514.EAL650R
                    DATA_20240514.WAL650D
                    DATA_20240514.WAL650G
                    DATA_20240514.WAL650R
                    DATA_20240514.S2L650D
                    DATA_20240514.S2L650G
                    DATA_20240514.S2L650R
                    DATA_20240514.S3L650D
                    DATA_20240514.S3L650G
                    DATA_20240514.S3L650R
                    DATA_20240514.S4L650D
                    DATA_20240514.S4L650G
                    DATA_20240514.S4L650R
                    DATA_20240514.S5L650D
                    DATA_20240514.S5L650G
                    DATA_20240514.S5L650R
                    DATA_20240514.S6L650D
                    DATA_20240514.S6L650G
                    DATA_20240514.S6L650R
                    DATA_20240514.S7L650D
                    DATA_20240514.S7L650G
                    DATA_20240514.S7L650R
                   };
               
GBands650 =    {
            DATA_20240514.EAL650G
            DATA_20240514.WAL650G
            DATA_20240514.S2L650G
            DATA_20240514.S3L650G
            DATA_20240514.S4L650G
            DATA_20240514.S5L650G
            DATA_20240514.S6L650G
            DATA_20240514.S7L650G
            };
        
DBands650 =    {
            DATA_20240514.EAL650D
            DATA_20240514.WAL650D
            DATA_20240514.S2L650D
            DATA_20240514.S3L650D
            DATA_20240514.S4L650D
            DATA_20240514.S5L650D
            DATA_20240514.S6L650D
            DATA_20240514.S7L650D
            };
        
RBMs650 =    {
            DATA_20240514.EAL650R
            DATA_20240514.WAL650R
            DATA_20240514.S2L650R
            DATA_20240514.S3L650R
            DATA_20240514.S4L650R
            DATA_20240514.S5L650R
            DATA_20240514.S6L650R
            DATA_20240514.S7L650R
            };
        
        
plotRaman(FilledSamples650, 0)
plotRaman(GBands, 0)
plotRaman(DBands, 0)
plotRaman(RBMs, 0)

%GDBand520 = GDBandPeaksCalculation(GDBand520, 1580,1600,1560,1570,1300,1400);
%exportGDBandPeaks(GDBand520, 'dsfsd.csv');