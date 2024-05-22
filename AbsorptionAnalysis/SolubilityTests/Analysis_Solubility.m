clc;
clear;
addpath('C:\Users\cborja\OneDrive - Universiteit Antwerpen\SWCNTs\');
%addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\');
import UsefulFunctions.*;

%rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements Data\Absorption\';
rootpath = 'C:\Users\cborja\OneDrive - Universiteit Antwerpen\Measurements Data\Absorption\';

%All paths as default

path_20240415a = [rootpath,'20240415\Acetone.csv'];
path_20240415b = [rootpath,'20240415\Ethanol.csv'];
path_20240415c = [rootpath,'20240415\Methanol.csv'];
path_20240415d = [rootpath,'20240415\THF.csv'];

path_20240422a = [rootpath,'20240422\Ethanol.csv'];
path_20240422b = [rootpath,'20240422\EthylAcetate.csv'];
path_20240422c = [rootpath,'20240422\Methanol.csv'];
path_20240422d = [rootpath,'20240422\THF.csv'];

path_20240506a = [rootpath,'20240506\Acetone.csv'];
path_20240506b = [rootpath,'20240506\Acetonitrile.csv'];
path_20240506c = [rootpath,'20240506\DMF.csv'];
path_20240506d = [rootpath,'20240506\THF.csv'];
path_20240506e = [rootpath,'20240506\EthylAcetate.csv'];
path_20240506f = [rootpath,'20240506\Methanol.csv'];

%Select the paths of interest
paths = {
           path_20240415a
           path_20240415b
           path_20240415c
           path_20240415d
           
           path_20240422a
           path_20240422b
           path_20240422c
           path_20240422d
           
           path_20240506a
           path_20240506b
           path_20240506c
           path_20240506d
           path_20240506e
           path_20240506f
        };
    

%Read and structure data from the paths
ReadAbsorptionFromPaths(paths);

%%
% TO DO > Include Methanol on folder 240506

BLines = {
            %DATA_20240415.Acetone_Baseline
            %DATA_20240506.Acetone_Baseline
            %DATA_20240506.Acetonitrile_Baseline
            %DATA_20240506.DMF_Baseline
            %DATA_20240422.Ethanol_Baseline
            %DATA_20240422.EthylAcetate_Baseline
            %DATA_20240415.Methanol_Baseline     %Error en Baselines
            %DATA_20240422.Methanol_Baseline
            %DATA_20240506.Methanol_Baseline

            %DATA_20240415.THF_Baseline
            %DATA_20240422.THF_Baseline
            %DATA_20240506.THF_Baseline
        };
    
%plotAbsorption(BLines, 3.0)

TCNQ = {
        DATA_20240415.TCNQ_THF
       
        DATA_20240422.TCNQ_THF
        DATA_20240422.TCNQ_EthylAcetate
        
        DATA_20240506.TCNQ_Acetone
        DATA_20240506.TCNQ_Acetonitrile
        DATA_20240506.TCNQ_DMF
        DATA_20240506.TCNQ_EthylAcetate
        DATA_20240506.TCNQ_EthylAcetate_B
        DATA_20240506.TCNQ_THF
       };
   
   
TCBQ = {
        DATA_20240415.TCBQ_Acetone
        DATA_20240415.TCBQ_Ethanol
        
        DATA_20240422.TCBQ_Ethanol
        DATA_20240422.TCBQ_EthylAcetate
        
        DATA_20240506.TCBQ_Acetonitrile
        DATA_20240506.TCBQ_Acetonitrile_B
        DATA_20240506.TCBQ_DMF
        DATA_20240506.TCBQ_EthylAcetate
        DATA_20240506.TCBQ_EthylAcetate_B
        DATA_20240506.TCBQ_THF
        DATA_20240506.TCBQ_Methanol
       };
   
FN = {
        DATA_20240415.FN_Methanol
        DATA_20240415.FN_Ethanol
        
        DATA_20240422.FN_Ethanol
        DATA_20240422.FN_EthylAcetate
        DATA_20240422.FN_Methanol
                
        DATA_20240506.FN_Acetone
        DATA_20240506.FN_Acetonitrile
        DATA_20240506.FN_Acetonitrile_B
        DATA_20240506.FN_DMF
        DATA_20240506.FN_THF
        DATA_20240506.FN_Methanol
        DATA_20240506.FN_Methanol_B

        };   
 
DCB = {
        DATA_20240415.DCB_Acetone
        DATA_20240415.DCB_THF
        
        DATA_20240422.DCB_THF
        DATA_20240422.DCB_EthylAcetate

        DATA_20240506.DCB_EthylAcetate
        DATA_20240506.DCB_THF
        DATA_20240506.DCB_THF_B
        DATA_20240506.DCB_Acetonitrile
        DATA_20240506.DCB_Acetonitrile_B
        DATA_20240506.DCB_DMF
        DATA_20240506.DCB_Methanol
        };  
   
CuCl = {
        DATA_20240506.CuCl_Acetone
        DATA_20240506.CuCl_Acetonitrile
        DATA_20240506.CuCl_Acetonitrile_B
        DATA_20240506.CuCl_DMF
        DATA_20240506.CuCl_EthylAcetate
        DATA_20240506.CuCl_THF
        DATA_20240506.CuCl_Methanol

        };
    
  
%plotAbsorption(TCNQ, 1.0)
%plotAbsorption(TCBQ, 1.0)
%plotAbsorption(FN, 1.0)
%plotAbsorption(DCB, 1.0)
%plotAbsorption(CuCl, 1.0)


