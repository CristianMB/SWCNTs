%% DATA LOAD
clc;
clear;

%Local path where general functions are stored
addpath('C:\Users\cborja\OneDrive - Universiteit Antwerpen\SWCNTs\');
% addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\');
import UsefulFunctions.*;

%Root path for data
%rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements Data\Raman\';
rootpath = 'C:\Users\cborja\OneDrive - Universiteit Antwerpen\Measurements Data\Raman\';

%Specific path
path_20240515 = [rootpath,'20240515\'];

%Select the paths of interest
paths = {
        path_20240515
        };

%Read and structure data from selected paths.
%This create a Data Structure whose fields are the spectra.
ReadRamanFromPaths(paths);


%% DATA LABELING

%DATA_20240515.Sample has 3 fields: X, Y, P (pixels) and N
DATA_20240515.BAL570C.N = 'Ben Sample A at 570nm - Carbide Peak';
 
%% BUILD A SAMPLE LIST

SL = {    
        DATA_20240515.BBL570G
        };
          
%% DATA CORRECTION
%Some of the data corrections are wavelenght dependent
WL = 570;
ClipLeft = 50;
ClipRight = 50;
PP = 1600;
Tol = 50;

%Flat Field correction if Low Dispersion Mode
%SL = FlatFieldCorrection(SL, DATA.FlatSpectrum);
%Remove Linear Background
SL = RemoveBackground(SL);
%Clip spectrum tails
SL = ClipSamples(SL, ClipLeft, ClipRight);
%Correct instrument response
SL = RemoveInclination(SL, WL);
SL = InstrumentCorrection(SL, WL);
%Custom Normalization
SL = NormalizeSample(SL, PP-Tol, PP+Tol);

plotRaman(SL, 0);