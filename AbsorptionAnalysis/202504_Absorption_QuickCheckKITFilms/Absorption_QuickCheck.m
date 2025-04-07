clc;
clear;
import UsefulFunctions.*;
rootpath = 'X:\Measurements Data\Absorption\';
addpath('X:\SWCNTs\SpecialMatlabFunctions\DrosteEffect-BrewerMap-3.2.5.0')
addpath('X:\SWCNTs')


%All paths as default

Films = [rootpath,'20250407\KIT_ME_SC_Films.csv'];


Refs= [rootpath,'References.csv'];

%Select the paths of interest
paths = {
        Refs   
        Films
        };
    

%Read and structure data from the paths
ReadAbsorptionFromPaths(paths);

% Labeling
            
% Plotting

F = {
        DATA_20250407.FilmA
        DATA_20250407.FilmB
%         DATA_20250407.Substrate
%         DATA_20250407.Glass_4mm_vial
            };
        
F = Normalize(F, 960, 1060, 'M')
plotAbsorption(F, 0.00);
