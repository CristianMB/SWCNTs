clc;
clear;
import UsefulFunctions.*;

rootpath = 'X:\Measurements Data\Absorption\';

%All paths as default

TTF_A = [rootpath,'20241122\TTF_MeOH.csv'];
TTF_B = [rootpath,'20241114\Solubility_TTF_MeOH.csv'];

%Select the paths of interest
paths = {
    TTF_A
    TTF_B
        };
    

%Read and structure data from the paths
ReadAbsorptionFromPaths(paths);


TTF = {
%         DATA_20241114.TTF_MeOH_LowC
%         DATA_20241114.TTF_MeOH_LowC_DilA
        DATA_20241122.TTF_MeOH_LC_1
        DATA_20241122.TTF_MeOH_LC_14
        DATA_20241122.TTF_MeOH_S_35
       };

plotAbsorption(TTF, 0.0)

   
