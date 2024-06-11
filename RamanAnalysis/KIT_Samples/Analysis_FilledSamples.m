clc;
clear;
addpath('C:\Users\cborja\OneDrive - Universiteit Antwerpen\SWCNTs\');
%addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\');
import UsefulFunctions.*;

%rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements Data\Raman\';
rootpath = 'C:\Users\cborja\OneDrive - Universiteit Antwerpen\Measurements Data\Raman\';

%All paths as default
path_20240610 = [rootpath,'20240610\'];

%Select the paths of interest

paths = {
    path_20240610
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

%%%--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%--------SAMPLE COMPARISION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Samples = {
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

Samples = {
        DATA_20240610.WAH514R
        DATA_20240610.EAH514R
        DATA_20240610.S2H514R
        DATA_20240610.S3H514R
        DATA_20240610.S4H514R
        DATA_20240610.S5H514R
        DATA_20240610.S6H514R
        DATA_20240610.S7H514R
        };

%for RBMs

Samples = FlatFieldCorrection(Samples, DATA_20240610.FFH514R)
Samples = UsefulFunctions.SubstractLinearBG(Samples, 137, 200)
Samples = UsefulFunctions.NormalizeSample(Samples,140, 160)

%for G/D Bands
% Samples = UsefulFunctions.SubstractLinearBG(Samples, 1250, 1650)
% Samples = UsefulFunctions.NormalizeSample(Samples,1585, 1595)

plotRaman(Samples, 0.0)
