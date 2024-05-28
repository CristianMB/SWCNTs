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

DATA_20240506.CuCl_Acetone.N='CuCl_Acetone 2,18mg in 20ml - conc 0,109 - Original';
DATA_20240506.CuCl_Acetonitrile.N='CuCl_Acetonitrile 1,86mg in 20ml - conc 0,093 - Original';
DATA_20240506.CuCl_Acetonitrile_B.N='CuCl_Acetonitrile_B 1,86mg in 20ml - conc 0,0465 - Diluited';
DATA_20240506.CuCl_DMF.N='CuCl_DMF 2,06mg in 20ml - conc 0,103 - Original';
DATA_.CuCl_Ethanol.N='CuCl_Ethanol 1,89mg in 20ml - conc 0,0945 - Original';
DATA_20240506.CuCl_EthylAcetate.N='CuCl_EthylAcetate 1,45mg in 20ml - conc 0,0725 - Original';
DATA_20240506.CuCl_Methanol.N='CuCl_Methanol 2,29mg in 20ml - conc 0,1145 - Original';
DATA_20240506.CuCl_THF.N='CuCl_THF 3,16mg in 20ml - conc 0,158 - Original';
DATA_20240415.DCB_THF.N='DCB in THF 21,35mg in 2ml - conc 10,675 - Original';
DATA_20240415.DCB_Acetone.N='DCB_Acetone 22,15mg in 2ml - conc 11,075 - Original';
DATA_20240506.DCB_Acetonitrile.N='DCB_Acetonitrile 1,38mg in 20ml - conc 0,069 - Original';
DATA_20240506.DCB_DMF.N='DCB_DMF 1,1mg in 20ml - conc 0,055 - Original';
DATA_.DCB_Ethanol.N='DCB_Ethanol 3,6mg in 20ml - conc 0,18 - Original';
DATA_20240506.DCB_EthylAcetate.N='DCB_EthylAcetate 1,42mg in 20ml - conc 0,071 - Original';
DATA_20240422.DCB_EthylAcetate.N='DCB_EthylAcetate 7,51mg in 10,2ml - conc 0,736274509803922 - Original';
DATA_20240506.DCB_Methanol.N='DCB_Methanol 1,23mg in 20ml - conc 0,0615 - Original';
DATA_20240506.DCB_THF.N='DCB_THF 21,35mg in 21ml - conc 1,01666666666667 - Diluited';
DATA_20240422.DCB_THF.N='DCB_THF 21,35mg in 8ml - conc 2,66875 - Diluited';
DATA_20240506.FN_Acetone.N='FN_Acetone 3,1mg in 20ml - conc 0,155 - Original';
DATA_20240506.FN_Acetonitrile.N='FN_Acetonitrile 1,45mg in 20ml - conc 0,0725 - Original';
DATA_20240506.FN_DMF.N='FN_DMF 1,29mg in 20ml - conc 0,0645 - Original';
DATA_20240415.FN_Ethanol.N='FN_Ethanol 21,3mg in 1,5ml - conc 14,2 - Original';
DATA_.FN_Ethanol.N='FN_Ethanol 21,3mg in 23,5ml - conc 0,906382978723404 - Diluited';
DATA_20240422.FN_Ethanol.N='FN_Ethanol 21,3mg in 8,5ml - conc 2,50588235294118 - Diluited';
DATA_20240422.FN_EthylAcetate.N='FN_EthylAcetate 13,27mg in 10,2ml - conc 1,30098039215686 - Original';
DATA_20240415.FN_Methanol.N='FN_Methanol 22,13mg in 0,5ml - conc 44,26 - Original';
DATA_20240506.FN_Methanol.N='FN_Methanol 22,13mg in 22,5ml - conc 0,983555555555556 - Diluited';
DATA_20240422.FN_Methanol.N='FN_Methanol 22,13mg in 6,5ml - conc 3,40461538461538 - Diluited';
DATA_20240422.FN_Methanol_B.N='FN_Methanol_B 22,13mg in 6,5ml - conc 0,486373626373626 - Diluited';
DATA_20240506.FN_THF.N='FN_THF 1,25mg in 20ml - conc 0,0625 - Original';
DATA_20240415.TCBQ_Acetone.N='TCBQ_Acetone 30,8mg in 5ml - conc 6,16 - Original';
DATA_20240506.TCBQ_Acetonitrile.N='TCBQ_Acetonitrile 1,89mg in 20ml - conc 0,0945 - Original';
DATA_20240506.TCBQ_DMF.N='TCBQ_DMF 1,87mg in 20ml - conc 0,0935 - Original';
DATA_20240422.TCBQ_Ethanol.N='TCBQ_Ethanol 22,06mg in 16ml - conc 1,37875 - Diluited';
DATA_.TCBQ_Ethanol.N='TCBQ_Ethanol 22,06mg in 24ml - conc 0,919166666666667 - Diluited';
DATA_20240415.TCBQ_Ethanol.N='TCBQ_Ethanol 22,06mg in 6ml - conc 3,67666666666667 - Original';
DATA_20240422.TCBQ_EthylAcetate.N='TCBQ_EthylAcetate 6,86mg in 10,2ml - conc 0,672549019607843 - Original';
DATA_20240506.TCBQ_EthylAcetate.N='TCBQ_EthylAcetate 6,86mg in 23,2ml - conc 0,295689655172414 - Diluited';
DATA_20240506.TCBQ_Methanol.N='TCBQ_Methanol 2,15mg in 20ml - conc 0,1075 - Original';
DATA_20240506.TCBQ_THF.N='TCBQ_THF 1,31mg in 20ml - conc 0,0655 - Original';
DATA_20240506.TCNQ_Acetone.N='TCNQ_Acetone 1,8mg in 20ml - conc 0,09 - Original';
DATA_20240506.TCNQ_Acetonitrile.N='TCNQ_Acetonitrile 1,3mg in 20ml - conc 0,065 - Original';
DATA_20240506.TCNQ_DMF.N='TCNQ_DMF 1,93mg in 10ml - conc 0,193 - Original';
DATA_.TCNQ_Ethanol.N='TCNQ_Ethanol 1,05mg in 20ml - conc 0,0525 - Original';
DATA_20240422.TCNQ_EthylAcetate.N='TCNQ_EthylAcetate 4,02mg in 10,2ml - conc 0,394117647058824 - Original';
DATA_20240506.TCNQ_EthylAcetate.N='TCNQ_EthylAcetate 4,02mg in 23,2ml - conc 0,173275862068966 - Diluited';
DATA_20240506.TCNQ_EthylAcetate_B.N='TCNQ_EthylAcetate_B 4,02mg in 23,2ml - conc 0,0433189655172414 - Diluited';
DATA_20240506.TCNQ_Methanol.N='TCNQ_Methanol 1,27mg in 20ml - conc 0,0635 - Original';
DATA_20240415.TCNQ_THF.N='TCNQ_THF 10,36mg in 15ml - conc 0,690666666666667 - Original';
DATA_20240422.TCNQ_THF.N='TCNQ_THF 10,36mg in 23ml - conc 0,450434782608696 - Diluited';
DATA_20240506.TCNQ_THF.N='TCNQ_THF 10,36mg in 40ml - conc 0,259 - Diluited';

%%

BLines = {
            DATA_20240506.Acetonitrile_Baseline
            DATA_20240506.Methanol_Baseline
            DATA_20240422.Ethanol_Baseline
            DATA_20240422.EthylAcetate_Baseline
            DATA_20240506.DMF_Baseline
            DATA_20240506.THF_Baseline
            DATA_20240506.Acetone_Baseline

%             DATA_20240415.Acetone_Baseline
%             DATA_20240415.Methanol_Baseline     %Error in Methanol
%             Baselines, should be flat
%             DATA_20240422.Methanol_Baseline
%             DATA_20240415.THF_Baseline
%             DATA_20240422.THF_Baseline
            %
        };
    
% plotAbsorption(BLines, 1.0)

TCNQ = {
        DATA_20240506.TCNQ_Acetonitrile
        DATA_20240506.TCNQ_Methanol
        DATA_20240422.TCNQ_EthylAcetate
        DATA_20240506.TCNQ_EthylAcetate
        DATA_20240506.TCNQ_EthylAcetate_B
        DATA_20240506.TCNQ_DMF
        DATA_20240415.TCNQ_THF
        DATA_20240422.TCNQ_THF
        DATA_20240506.TCNQ_THF
        DATA_20240506.TCNQ_Acetone
       };
TCNQ_OK = {
        DATA_20240506.TCNQ_Acetonitrile
        DATA_20240506.TCNQ_Methanol
        DATA_20240506.TCNQ_EthylAcetate_B
        DATA_20240506.TCNQ_DMF
        DATA_20240506.TCNQ_THF
        DATA_20240506.TCNQ_Acetone
       };
%plotAbsorption(TCNQ_OK, 0.0)
   
   
TCBQ = {
        DATA_20240506.TCBQ_Acetonitrile
        DATA_20240506.TCBQ_Acetonitrile_B
        DATA_20240506.TCBQ_Methanol
        DATA_20240415.TCBQ_Ethanol
        DATA_20240422.TCBQ_Ethanol
        DATA_20240422.TCBQ_EthylAcetate
        DATA_20240506.TCBQ_EthylAcetate
        DATA_20240506.TCBQ_EthylAcetate_B
        DATA_20240506.TCBQ_DMF
        DATA_20240506.TCBQ_THF
        DATA_20240415.TCBQ_Acetone
       };

TCBQOK = {
        DATA_20240506.TCBQ_Acetonitrile
        DATA_20240506.TCBQ_Acetonitrile_B
        DATA_20240506.TCBQ_Methanol
        %DATA_20240415.TCBQ_Ethanol
        %DATA_20240422.TCBQ_Ethanol
        %DATA_20240422.TCBQ_EthylAcetate
        %DATA_20240506.TCBQ_EthylAcetate
        DATA_20240506.TCBQ_EthylAcetate_B
        DATA_20240506.TCBQ_DMF
        DATA_20240506.TCBQ_THF
        %DATA_20240415.TCBQ_Acetone
       };

%plotAbsorption(TCBQOK, 0.0)

   
FN = {  
        DATA_20240506.FN_Acetonitrile
        DATA_20240506.FN_Acetonitrile_B
        DATA_20240415.FN_Methanol
        DATA_20240422.FN_Methanol
        DATA_20240506.FN_Methanol
        DATA_20240506.FN_Methanol_B
        DATA_20240415.FN_Ethanol
        DATA_20240422.FN_Ethanol
        DATA_20240422.FN_EthylAcetate
        DATA_20240506.FN_DMF
        DATA_20240506.FN_THF                    
        DATA_20240506.FN_Acetone
        }; 
    
FNOK = {  
         DATA_20240506.FN_Acetonitrile
         DATA_20240506.FN_Methanol_B  

        }; 
    
plotAbsorption(FNOK, 0.0)

  
 
DCB = {
        DATA_20240506.DCB_Acetonitrile
        DATA_20240506.DCB_Acetonitrile_B         
        DATA_20240506.DCB_Methanol
        DATA_20240422.DCB_EthylAcetate
        DATA_20240506.DCB_EthylAcetate
        DATA_20240506.DCB_DMF
        DATA_20240415.DCB_THF
        DATA_20240422.DCB_THF
        DATA_20240506.DCB_THF
        DATA_20240506.DCB_THF_B
        DATA_20240415.DCB_Acetone
        };  
DCBOK = {
        DATA_20240506.DCB_Acetonitrile
        DATA_20240506.DCB_Acetonitrile_B         
        DATA_20240506.DCB_Methanol
        DATA_20240422.DCB_EthylAcetate
        DATA_20240506.DCB_EthylAcetate
        DATA_20240506.DCB_DMF
        DATA_20240506.DCB_THF
        DATA_20240506.DCB_THF_B
        };    
%  plotAbsorption(DCBOK, 0.0)


CuCl = {
        DATA_20240506.CuCl_Acetonitrile
        DATA_20240506.CuCl_Acetonitrile_B
        DATA_20240506.CuCl_Methanol
        DATA_20240506.CuCl_EthylAcetate
        DATA_20240506.CuCl_DMF
        DATA_20240506.CuCl_THF
        DATA_20240506.CuCl_Acetone
        }; 
    
CuClOK = {
        DATA_20240506.CuCl_Acetonitrile
        DATA_20240506.CuCl_Acetonitrile_B
        DATA_20240506.CuCl_Methanol
        };     
    
% plotAbsorption(CuClOK, 0.0)


