clc;
clear;
import UsefulFunctions.*;
rootpath = 'X:\Measurements Data\Absorption\';

%All paths as default

TCNQ_Acetone = [rootpath,'20250617\TCNQ_Acetone.csv'];
TCNQ_Acetonitrile = [rootpath,'20250617\TCNQ_Acetonitrile.csv'];
TCNQ_Chloroform = [rootpath,'20250617\TCNQ_Chloroform.csv'];
TCNQ_DCM = [rootpath,'20250617\TCNQ_DCM.csv'];
TCNQ_EthylAcetate = [rootpath,'20250617\TCNQ_EthylAcetate.csv'];
TCNQ_MeOH = [rootpath,'20250617\TCNQ_MeOH.csv'];
TCNQ_THF = [rootpath,'20250617\TCNQ_THF.csv'];
TCNQ_Toluene = [rootpath,'20250714\TCNQ_Toluene.csv'];

TCNQ_DCM_LC = [rootpath,'20250714\Rinsing_S18_TCNQ_DCM.csv'];


%Select the paths of interest
paths = {
    TCNQ_Acetone
    TCNQ_Acetonitrile
    TCNQ_Chloroform
    TCNQ_DCM
    TCNQ_DCM_LC
    TCNQ_EthylAcetate
    TCNQ_MeOH
    TCNQ_THF
    TCNQ_Toluene
};
    

%% Read and structure data from the paths
ReadAbsorptionFromPaths(paths);

%% Labeling
DATA_20250617.Acetone.N='Acetone';
DATA_20250617.TCNQ_Acetone_1.N='TCNQ_Acetone_1';
DATA_20250617.TCNQ_Acetone_30.N='TCNQ_Acetone_30';
DATA_20250617.Acetonitrile.N='Acetonitrile';
DATA_20250617.TCNQ_Acetonitrile_1.N='TCNQ_Acetonitrile_1';
DATA_20250617.TCNQ_Acetonitrile_30.N='TCNQ_Acetonitrile_30';
DATA_20250617.Chloroform.N='Chloroform';
DATA_20250617.TCNQ_Chloroform_1.N='TCNQ_Chloroform_1';
DATA_20250617.TCNQ_Chloroform_15.N='TCNQ_Chloroform_15';
DATA_20250617.TCNQ_Chloroform_30.N='TCNQ_Chloroform_30';
DATA_20250617.DCM.N='DCM';
DATA_20250617.TCNQ_DCM_1.N='TCNQ_DCM_1';
DATA_20250617.TCNQ_DCM_30.N='TCNQ_DCM_30';
DATA_20250617.EthylAcetate.N='EthylAcetate';
DATA_20250617.TCNQ_EthylAcetate_1.N='TCNQ_EthylAcetate_1';
DATA_20250617.TCNQ_EthylAcetate_30.N='TCNQ_EthylAcetate_30';
DATA_20250617.MeOH.N='MeOH';
DATA_20250617.TCNQ_MeOH_1.N='TCNQ_MeOH_1';
DATA_20250617.TCNQ_MeOH_30.N='TCNQ_MeOH_30';
DATA_20250617.THF.N='THF';
DATA_20250617.TCNQ_THF_1.N='TCNQ_THF_1';
DATA_20250617.TCNQ_THF_30.N='TCNQ_THF_30';

%% Corrections

Con_TCNQ_Acetone = 0.6807;
Con_TCNQ_Acetonitrile = 0.9974;
Con_TCNQ_Chloroform = 0.7480;
Con_TCNQ_DCM = 0.6108;
Con_TCNQ_DCM_LC = 0.1203;

Con_TCNQ_EthylAcetate = 0.5680;
Con_TCNQ_MeOH = 0.6998;
Con_TCNQ_THF = 0.7702;
Con_TCNQ_Toluene = 0.8990;


% DATA_20250617.TCNQ_Acetone_1.Y = DATA_20250617.TCNQ_Acetone_1.Y/Con_TCNQ_Acetone;
% DATA_20250617.TCNQ_Acetone_30.Y = DATA_20250617.TCNQ_Acetone_30.Y*30/Con_TCNQ_Acetone;
% 
% DATA_20250617.TCNQ_Acetonitrile_1.Y = DATA_20250617.TCNQ_Acetonitrile_1.Y/Con_TCNQ_Acetonitrile;
% DATA_20250617.TCNQ_Acetonitrile_30.Y = DATA_20250617.TCNQ_Acetonitrile_30.Y*30/Con_TCNQ_Acetonitrile;
% 
% DATA_20250617.TCNQ_Chloroform_1.Y = DATA_20250617.TCNQ_Chloroform_1.Y/Con_TCNQ_Chloroform;
% DATA_20250617.TCNQ_Chloroform_15.Y = DATA_20250617.TCNQ_Chloroform_15.Y*15/Con_TCNQ_Chloroform;
% DATA_20250617.TCNQ_Chloroform_30.Y = DATA_20250617.TCNQ_Chloroform_30.Y*30/Con_TCNQ_Chloroform;
% 
% DATA_20250617.TCNQ_DCM_1.Y = DATA_20250617.TCNQ_DCM_1.Y/Con_TCNQ_DCM;
% DATA_20250617.TCNQ_DCM_30.Y = DATA_20250617.TCNQ_DCM_30.Y*30/Con_TCNQ_DCM;
% 
% DATA_20250617.TCNQ_EthylAcetate_1.Y = DATA_20250617.TCNQ_EthylAcetate_1.Y/Con_TCNQ_EthylAcetate;
% DATA_20250617.TCNQ_EthylAcetate_30.Y = DATA_20250617.TCNQ_EthylAcetate_30.Y*30/Con_TCNQ_EthylAcetate;
% 
% DATA_20250617.TCNQ_MeOH_1.Y = DATA_20250617.TCNQ_MeOH_1.Y/Con_TCNQ_MeOH;
% DATA_20250617.TCNQ_MeOH_30.Y = DATA_20250617.TCNQ_MeOH_30.Y*30/Con_TCNQ_MeOH;
% 
% DATA_20250617.TCNQ_THF_1.Y = DATA_20250617.TCNQ_THF_1.Y/Con_TCNQ_THF;
% DATA_20250617.TCNQ_THF_30.Y = DATA_20250617.TCNQ_THF_30.Y*30/Con_TCNQ_THF;
% 
% DATA_20250714.TCNQ_Toluene_30.Y = DATA_20250714.TCNQ_Toluene_30.Y*30/Con_TCNQ_Toluene;
% 
% DATA_20250714.TCNQ_DCM_15.Y = DATA_20250714.TCNQ_DCM_15.Y*15/Con_TCNQ_DCM_LC;


%% Plotting

TCNQ = {
%         DATA_20250617.TCNQ_Acetone_1
        DATA_20250617.TCNQ_Acetone_30
        
%         DATA_20250617.TCNQ_Acetonitrile_1
        DATA_20250617.TCNQ_Acetonitrile_30
        
%         DATA_20250617.TCNQ_Chloroform_1
%         DATA_20250617.TCNQ_Chloroform_15
        DATA_20250617.TCNQ_Chloroform_30
        
%         DATA_20250617.TCNQ_DCM_1
        DATA_20250617.TCNQ_DCM_30
        
%         DATA_20250617.TCNQ_EthylAcetate_1
        DATA_20250617.TCNQ_EthylAcetate_30
        
%         DATA_20250617.TCNQ_MeOH_1
        DATA_20250617.TCNQ_MeOH_30

%         DATA_20250617.TCNQ_THF_1
        DATA_20250617.TCNQ_THF_30
        
        DATA_20250714.TCNQ_Toluene_30

       };

% plotAbsorption(TCNQ, 0.0)

%%

TCNQ = {
        DATA_20250617.TCNQ_DCM_1
        DATA_20250617.TCNQ_DCM_30
        
        DATA_20250714.TCNQ_DCM_2
        DATA_20250714.TCNQ_DCM_15
        DATA_20250714.TCNQ_DCM_30
       };

plotAbsorption(TCNQ, 0.0)
