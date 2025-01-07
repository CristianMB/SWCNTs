clc;
clear;
addpath('X:\SWCNTs');
% addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\');
import UsefulFunctions.*;

% rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements Data\Absorption\';
rootpath = 'X:\Measurements Data\ExternalData\';

%All paths as default
path_Examples = [rootpath,'Raman_2020412_TestData_LorentzianFitting_EmptyWater\'];


%Select the paths of interest

paths = {
    path_Examples
    };


ReadRamanFromPaths(paths,2);




%%%--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%--------SAMPLE COMPARISION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Examples = {
    DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.SR0H514R_241212
    DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.SR1H514R_241212
%     DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.FFH514R_241212
        
    DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.EAH514R_240517
    DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.WAH514R_240517
%    DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.FF514R_240517
    
    DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.EAH514R_241006
    DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.WAH514R_241006
%     DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.FFH514R_241006
        
    DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.S240111J_240111
    DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.S240111S_240111
%     DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.FLATHD_240111
}


D241212 = {    
        DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.SR0H514R_241212
        DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.SR1H514R_241212
        };
    
D240517 = {
        DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.EAH514R_240517
        DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.WAH514R_240517    
            };
        
D241016 = {
        DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.EAH514R_241006
        DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.WAH514R_241006
            };
        
D240111 = {
        DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.S240111J_240111
        DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.S240111S_240111
            };

        
D241212 = FlatFieldCorrection(D241212, DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.FFH514R_241212);
D240517 = FlatFieldCorrection(D240517, DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.FF514R_240517);
D241016 = FlatFieldCorrection(D241016, DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.FFH514R_241006);
D240111 = FlatFieldCorrection(D240111, DATA_Raman_2020412_TestData_LorentzianFitting_EmptyWater.FLATHD_240111); 


%FilledSamples514R = FlatFieldCorrection(FilledSamples514R, DATA_20240517.FF514R);
plotRaman(Examples, 0)
