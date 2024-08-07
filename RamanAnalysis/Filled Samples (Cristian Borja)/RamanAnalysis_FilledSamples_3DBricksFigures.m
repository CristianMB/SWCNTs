clc;
clear;
addpath('X:\SWCNTs');
% addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\');
import UsefulFunctions.*;

% rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements Data\Absorption\';
rootpath = 'X:\Measurements Data\Raman\';

%All paths as default
path_SF514 = [rootpath,'20240111\'];
path_CB514 = [rootpath,'20240517\'];

%Select the paths of interest

paths = {
    path_SF514
    path_CB514
    };


ReadRamanFromPaths(paths);


%%%--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DATA_20240111.S240111N.N='TTF';
DATA_20240517.S7H514R.N='Dodecane';
DATA_20240517.S6H514R.N='Hexadecane';
DATA_20240517.S2H514R.N='PCE';
DATA_20240517.S3H514R.N='TCE';
DATA_20240517.S5H514R.N='TDAE';
DATA_20240517.S4H514R.N='TEMED';
DATA_20240517.EAH514R.N='Empty';
DATA_20240517.WAH514R.N='Water';

DATA_20240111.S240111E.N='TTF';
DATA_20240517.S7L514GD.N='Dodecane';
DATA_20240517.S6L514GD.N='Hexadecane';
DATA_20240517.S2L514GD.N='PCE';
DATA_20240517.S3L514GD.N='TCE';
DATA_20240517.S5L514GD.N='TDAE';
DATA_20240517.S4L514GD.N='TEMED';
DATA_20240517.EAL514GD.N='Empty';
DATA_20240517.WAL514GD.N='Water';

%%%--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DATA_20240426.BAL650R.Y = DATA_20240426.BAL650R.Y/2;
%DATA_20240426.BAL650RB.Y = DATA_20240426.BAL650RB.Y/2;

%%%--------SAMPLE COMPARISION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%GBand
FilledSamples514G = {
        DATA_20240517.EAL514GD
        DATA_20240111.S240111E
        DATA_20240517.S2L514GD
        DATA_20240517.S3L514GD
        DATA_20240517.S4L514GD
        DATA_20240517.S5L514GD
        DATA_20240517.S6L514GD
        DATA_20240517.S7L514GD
        DATA_20240517.WAL514GD
    }
% %Normalization
LG = 1400;
HG = 1500;

NP = 1560;
Tol = 20;

plotRaman(FilledSamples514G,0)

% FilledSamples514G = SubstractLinearBG(FilledSamples514G,LG, HG);       
% FilledSamples514G = NormalizeSample(FilledSamples514G,NP-Tol, NP+Tol);       
% plotRaman(FilledSamples514G, 1.5)
% hold on
% xline(1343, '--','LineWidth', 1.0,'HandleVisibility','off')
% xline(1565, '--','LineWidth', 1.0,'HandleVisibility','off')
% xline(1589, '--','LineWidth', 1.0,'HandleVisibility','off')
% xlim([1300 1650])
% % ylim([-2 1.0])



%%%%%RBMs
FilledSamples514R = {
        DATA_20240517.EAH514R
        DATA_20240111.S240111N
        DATA_20240517.S2H514R
        DATA_20240517.S3H514R
        DATA_20240517.S4H514R
        DATA_20240517.S5H514R
        DATA_20240517.S6H514R
        DATA_20240517.S7H514R
        DATA_20240517.WAH514R
    }


% %Normalization
LG = 1450;
HG = 1640;

NP = 170;
Tol = 20;

% FilledSamples514R = FlatFieldCorrection(FilledSamples514R, DATA_20240517.FF514R);
% FilledSamples514R = NormalizeSample(FilledSamples514R,NP-Tol, NP+Tol);       
% plotRaman(FilledSamples514R, 0.65)
% hold on
% xline(167.5, '--','LineWidth', 1.0,'HandleVisibility','off')
% xline(151.5, '--','LineWidth', 1.0,'HandleVisibility','off')
% xline(173.5, '--','LineWidth', 1.0,'HandleVisibility','off')
% xline(179.5, '--','LineWidth', 1.0,'HandleVisibility','off')
% xlim([140 230])
% % ylim([-2 1.0])

