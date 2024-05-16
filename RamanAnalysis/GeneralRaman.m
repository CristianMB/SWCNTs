clc;
clear;

%% 
%%%--------IMPORT DATA--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('C:\Users\cborja\OneDrive - Universiteit Antwerpen\SWCNTs\');
%addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\');
import UsefulFunctions.*;

%rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements Data\Raman\';
rootpath = 'C:\Users\cborja\OneDrive - Universiteit Antwerpen\Measurements Data\Raman\';

%All paths as default
path_20240111 = [rootpath,'20240111\'];
path_20240320 = [rootpath,'20240320\'];
path_20240321 = [rootpath,'20240321\'];
path_20240325 = [rootpath,'20240325\'];
path_20240426 = [rootpath,'20240426\'];
path_20240514 = [rootpath,'20240514\'];
path_20240515 = [rootpath,'20240515\'];

%Select the paths of interest

paths = {
    path_20240514,
    path_20240515
    };



ReadRamanFromPaths(paths);

%% 
%%%--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DATA_20240325.REH680R.N='SC Empty@SWCNT (680.02 nm)';
%DATA_20240325.RWH680R.N='SC Water@SWCNT (680.02 nm)';
%DATA_20240325.S2H680R.N='CB PCE@SWCNT Dial. DGU C (Filled) (680.02 nm)';

%% 
%%%--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FlatFields Normalization

%GBand520 = {
%        DATA_20240321.S2H520G,
%        DATA_20240321.S3H520G,
%        DATA_20240321.S4H520G,
%        DATA_20240321.S5H520G,
%        DATA_20240321.S6H520G,
%        DATA_20240321.S7H520G,
%        };


%SList = FlatFieldCorrection(GBand520, DATA_20240111.FLATHD)

%FF Correction to GBand region in LD Mode
%DATA_20240321.S2L520G.Y = DATA_20240321.S2L520G.Y ./ DATA_20240321.FFL520G.Y;
%DATA_20240321.S3L520G.Y = DATA_20240321.S3L520G.Y ./ DATA_20240321.FFL520G.Y;
%DATA_20240321.S4L520G.Y = DATA_20240321.S4L520G.Y ./ DATA_20240321.FFL520G.Y;
%DATA_20240321.S5L520G.Y = DATA_20240321.S5L520G.Y ./ DATA_20240321.FFL520G.Y;
%DATA_20240321.S6L520G.Y = DATA_20240321.S6L520G.Y ./ DATA_20240321.FFL520G.Y;
%DATA_20240321.S7L520G.Y = DATA_20240321.S7L520G.Y ./ DATA_20240321.FFL520G.Y;

%FF Correction to GBand region in HD Mode
%DATA_20240321.S2H520G.Y = DATA_20240321.S2H520G.Y ./ DATA_20240321.FFH520G.Y;
%DATA_20240321.S3H520G.Y = DATA_20240321.S3H520G.Y ./ DATA_20240321.FFH520G.Y;
%DATA_20240321.S4H520G.Y = DATA_20240321.S4H520G.Y ./ DATA_20240321.FFH520G.Y;
%DATA_20240321.S5H520G.Y = DATA_20240321.S5H520G.Y ./ DATA_20240321.FFH520G.Y;
%DATA_20240321.S6H520G.Y = DATA_20240321.S6H520G.Y ./ DATA_20240321.FFH520G.Y;
%DATA_20240321.S7H520G.Y = DATA_20240321.S7H520G.Y ./ DATA_20240321.FFH520G.Y;

%% 
%%%--------SAMPLE COMPARISION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
%LG = 1560;
%HG = 1570;
%GDBand520 = SubstractLinearBG(GDBand520, 1400, 1500);
%GDBand520 = NormalizeSample(GDBand520,LG, HG);

plotRaman(GDBand520, 0)
%% 
%%%--------PEAKS EXPORT --------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%GDBand520 = GDBandPeaksCalculation(GDBand520, 1580,1600,1560,1570,1300,1400);
%exportGDBandPeaks(GDBand520, 'dsfsd.csv');




