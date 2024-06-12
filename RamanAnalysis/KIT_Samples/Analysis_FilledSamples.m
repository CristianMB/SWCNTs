clc;
clear;
addpath('C:\Users\cborja\OneDrive - Universiteit Antwerpen\SWCNTs\');
%addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\');
import UsefulFunctions.*;

%rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements Data\Raman\';
rootpath = 'C:\Users\cborja\OneDrive - Universiteit Antwerpen\Measurements Data\Raman\';

%All paths as default
path_20240610 = [rootpath,'20240610\'];
path_20240612 = [rootpath,'20240612\'];

%Select the paths of interest

paths = {
    path_20240610
    path_20240612
    };


ReadRamanFromPaths(paths);


%%%--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DATA_20240610.EAH514R.N='Empty ArcSWCNTs';
DATA_20240610.EAL514GD.N='Empty ArcSWCNTs';
DATA_20240610.FFH514R.N='Flat';
DATA_20240610.LL514A.N='Laser';
DATA_20240610.LL514B.N='Laser';
DATA_20240610.S2H514R.N='PCE@P2-ASWCNTs';
DATA_20240610.S2L514GD.N='PCE@P2-ASWCNTs';
DATA_20240610.S3H514R.N='TCE@P2-ASWCNTs';
DATA_20240610.S3L514GD.N='TCE@P2-ASWCNTs';
DATA_20240610.S4H514R.N='TEMED@P2-ASWCNTs';
DATA_20240610.S4L514GD.N='TEMED@P2-ASWCNTs';
DATA_20240610.S5H514R.N='TDAE@P2-ASWCNTs';
DATA_20240610.S5L514GD.N='TDAE@P2-ASWCNTs';
DATA_20240610.S6H514R.N= 'Hexadecane@P2-ASWCNTs';
DATA_20240610.S6L514GD.N= 'Hexadecane@P2-ASWCNTs';
DATA_20240610.S7H514R.N='Dodecane@P2-ASWCNTs';
DATA_20240610.S7L514GD.N='Dodecane@P2-ASWCNTs';
DATA_20240610.WAH514R.N='Water ArcSWCNTs';
DATA_20240610.WAL514GD.N='Water ArcSWCNTs';

DATA_20240612.S2L514R.N='KIT PCE@SWCNT Dial';
DATA_20240612.S3L514R.N='KIT TCE@SWCNT Dial';
DATA_20240612.S4L514R.N='KIT TEMED@SWCNT Dial';
DATA_20240612.S5L514R.N='KIT TDAE@SWCNT Dial';
DATA_20240612.S6L514R.N='KIT Hexadecane@SWCNT Dial';
DATA_20240612.S7L514R.N='KIT Dodecane@SWCNT Dial';
DATA_20240612.LL514A.N='Laser';
DATA_20240612.WAL514R.N='Water Filled Arc SWCNTs';
DATA_20240612.EAL514R.N='Empty Arc SWCNTs';

DATA_20240612.S2L514GD.N='KIT PCE@SWCNT Dial';
DATA_20240612.S3L514GD.N='KIT TCE@SWCNT Dial';
DATA_20240612.S4L514GD.N='KIT TEMED@SWCNT Dial';
DATA_20240612.S5L514GD.N='KIT TDAE@SWCNT Dial';
DATA_20240612.S6L514GD.N='KIT Hexadecane@SWCNT Dial';
DATA_20240612.S7L514GD.N='KIT Dodecane@SWCNT Dial';
DATA_20240612.WAL514GD.N='Water Filled Arc SWCNTs';
DATA_20240612.EAL514GD.N='Empty Arc SWCNTs';


%%%--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%--------SAMPLE COMPARISION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SamplesG = {
% DATA_20240610.EAH514R
DATA_20240610.EAL514GD
% DATA_20240610.S2H514R
DATA_20240610.S2L514GD
% DATA_20240610.S3H514R
DATA_20240610.S3L514GD
% DATA_20240610.S4H514R
DATA_20240610.S4L514GD
% DATA_20240610.S5H514R
DATA_20240610.S5L514GD
% DATA_20240610.S6H514R
DATA_20240610.S6L514GD
% DATA_20240610.S7H514R
DATA_20240610.S7L514GD
% DATA_20240610.WAH514R
DATA_20240610.WAL514GD    
};

SamplesR = {
        DATA_20240610.WAH514R
        DATA_20240610.EAH514R
        DATA_20240610.S2H514R
        DATA_20240610.S3H514R
        DATA_20240610.S4H514R
        DATA_20240610.S5H514R
        DATA_20240610.S6H514R
        DATA_20240610.S7H514R
        };
    
SamplesDialR = {
        DATA_20240612.WAL514R
        DATA_20240612.EAL514R
        DATA_20240612.S2L514R
        DATA_20240612.S3L514R
        DATA_20240612.S4L514R
        DATA_20240612.S5L514R
        DATA_20240612.S6L514R
        DATA_20240612.S7L514R
        };    
SamplesDialG = {
        DATA_20240612.WAL514GD
        DATA_20240612.EAL514GD
        DATA_20240612.S2L514GD
        DATA_20240612.S3L514GD
        DATA_20240612.S4L514GD
        DATA_20240612.S5L514GD
        DATA_20240612.S6L514GD
        DATA_20240612.S7L514GD
        };
    
%for RBMs
% 
SamplesR = FlatFieldCorrection(SamplesR, DATA_20240610.FFH514R)
SamplesR = UsefulFunctions.SubstractLinearBG(SamplesR, 137, 200)
SamplesR = UsefulFunctions.NormalizeSample(SamplesR,140, 160)

%for G/D Bands 
SamplesG = UsefulFunctions.SubstractLinearBG(SamplesG, 1250, 1650)
SamplesG = UsefulFunctions.NormalizeSample(SamplesG,1585, 1595)


%For DialSamples
SamplesDialR = UsefulFunctions.SubstractLinearBG(SamplesDialR, 137, 200)
SamplesDialR = UsefulFunctions.NormalizeSample(SamplesDialR,140, 160)
% 
% SamplesDialG = UsefulFunctions.SubstractLinearBG(SamplesDialG, 1250, 1650)
% SamplesDialG = UsefulFunctions.NormalizeSample(SamplesDialG,1585, 1595)

plotRaman(SamplesR, 0.7)
plotRaman(SamplesG, 0.0)
plotRaman(SamplesDialR, 0.7)
plotRaman(SamplesDialG, 0.0)
