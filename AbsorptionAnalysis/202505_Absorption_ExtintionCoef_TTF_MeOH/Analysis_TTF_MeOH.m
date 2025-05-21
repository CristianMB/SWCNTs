clc;
clear;
import UsefulFunctions.*;
rootpath = 'X:\Measurements Data\Absorption\';
addpath('X:\SWCNTs\SpecialMatlabFunctions\DrosteEffect-BrewerMap-3.2.5.0')
addpath('X:\SWCNTs')


%All paths as default

OriginalSamp = [rootpath,'20250305\DWCNT_Megasonicated_KIT.csv'];
Dial = [rootpath,'20250318\DialysisCheckDay1.csv'];
DialFinal = [rootpath,'20250321\LastDialisysAndSamplePostDial.csv'];

Refs= [rootpath,'References.csv'];

TTFSOLUB= [rootpath,'20250509\TTF_MeOH_Solubility.csv'];

%Select the paths of interest
paths = {
    OriginalSamp
    Dial
    Refs
    DialFinal
    
    TTFSOLUB
        };
    

%Read and structure data from the paths
ReadAbsorptionFromPaths(paths);

% % % 
% % % Spec = {
% % % %             DATA_References.WaterInD2O
% % %             DATA_20250305.Megasonicated_KIT
% % % %             DATA_20250318.DialWaste_D1
% % %             DATA_20250321.DWCNT_Megasonicated_KIT_Dial
% % %        };
% % % 
% % % Spec = Normalize(Spec, 485, 490, 'M');
% % % 
% % % plotAbsorption(Spec, 0.0);


% DATA_20250509.TTF_MeOH_LC_3.Y = DATA_20250509.TTF_MeOH_LC_3.Y*3.3333
% DATA_20250509.TTF_MeOH_LC_5.Y = DATA_20250509.TTF_MeOH_LC_5.Y*5
% DATA_20250509.TTF_MeOH_LC_10.Y = DATA_20250509.TTF_MeOH_LC_10.Y*10
% DATA_20250509.TTF_MeOH_LC_70.Y = DATA_20250509.TTF_MeOH_LC_70.Y*70
% 
% DATA_20250509.TTF_MeOH_SAT_10.Y = DATA_20250509.TTF_MeOH_SAT_10.Y*10
% DATA_20250509.TTF_MeOH_SAT_70.Y = DATA_20250509.TTF_MeOH_SAT_70.Y*70
% DATA_20250509.TTF_MeOH_SAT_100.Y = DATA_20250509.TTF_MeOH_SAT_100.Y*100

Spec = {
%             DATA_20250509.TTF_MeOH_LC_3
%             DATA_20250509.TTF_MeOH_LC_5
%             DATA_20250509.TTF_MeOH_LC_10
%             DATA_20250509.TTF_MeOH_LC_70
            
            DATA_20250509.TTF_MeOH_SAT_10
            DATA_20250509.TTF_MeOH_SAT_70
            DATA_20250509.TTF_MeOH_SAT_100
       };
% plotAbsorption(Spec, 0.0);

Spec = {
            DATA_20250509.TTF_MeOH_LC_3
            DATA_20250509.TTF_MeOH_LC_5
            DATA_20250509.TTF_MeOH_LC_10
            DATA_20250509.TTF_MeOH_LC_70
            
%             DATA_20250509.TTF_MeOH_SAT_10
%             DATA_20250509.TTF_MeOH_SAT_70
%             DATA_20250509.TTF_MeOH_SAT_100
       };
% plotAbsorption(Spec, 0.0);

AllSpecs = {
            DATA_20250509.TTF_MeOH_LC_3
            DATA_20250509.TTF_MeOH_LC_5
            DATA_20250509.TTF_MeOH_LC_10
            DATA_20250509.TTF_MeOH_LC_70
            
            DATA_20250509.TTF_MeOH_SAT_10
            DATA_20250509.TTF_MeOH_SAT_70
            DATA_20250509.TTF_MeOH_SAT_100
       };
   
plotAbsorption(AllSpecs, 0.0);
   
Abs_303 = get_Y_values(AllSpecs, 303)
Abs_315 = get_Y_values(AllSpecs, 315)
Abs_362 = get_Y_values(AllSpecs, 362)
Abs_447 = get_Y_values(AllSpecs, 447)


function Y_values = get_Y_values(D_list, value)
    N = numel(D_list);
    Y_values = nan(1, N);

    for k = 1:N
        DATA = D_list{k};
        idx = find(DATA.X == value, 1); % Find index of value in X
        if ~isempty(idx)
            Y_values(k) = DATA.Y(idx);  % Get corresponding Y
        end
    end
end