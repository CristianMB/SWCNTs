clc;
clear;
import UsefulFunctions.*;
rootpath = 'X:\Measurements Data\Absorption\';
addpath('X:\SWCNTs\SpecialMatlabFunctions\DrosteEffect-BrewerMap-3.2.5.0')
addpath('X:\SWCNTs')


%All paths as default

OriginalSamp = [rootpath,'20250305\DWCNT_Megasonicated_KIT.csv'];
Dial = [rootpath,'20250318\DialysisCheckDay1.csv'];
Refs= [rootpath,'References.csv'];

%Select the paths of interest
paths = {
    OriginalSamp
    Dial
    Refs
        };
    

%Read and structure data from the paths
ReadAbsorptionFromPaths(paths);


Spec = {
            DATA_References.WaterInD2O
            DATA_20250305.Megasonicated_KIT
            DATA_20250318.DialWaste_D1
       };

% Spec = Normalize(Spec, 200, 1900, 'M');

plotAbsorption(Spec, 0.0);

   
