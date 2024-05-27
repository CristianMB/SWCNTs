clc;
clear;
addpath('C:\Users\cborja\OneDrive - Universiteit Antwerpen\SWCNTs\');
%addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\');
import UsefulFunctions.*;

%rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements Data\Absorption\';
rootpath = 'C:\Users\cborja\OneDrive - Universiteit Antwerpen\Measurements Data\Absorption\';

%All paths as default
path_baselines = [rootpath,'References.csv'];
path_20231206 = [rootpath,'20231206\alldata20231206.csv'];
path_20240202 = [rootpath,'20240202\AllData.csv'];
path_20240212a = [rootpath,'20240212\AllDataCentrifuge.csv'];
path_20240214 = [rootpath,'20240214\DGU.csv'];
path_20240216 = [rootpath,'20240216\DGUC.csv'];
path_20240304 = [rootpath,'20240304\S6S7.csv'];
path_20240305 = [rootpath,'20240305\S5.csv'];
path_20240307 = [rootpath,'20240307\DGUS5S6S7.csv'];
path_20240308 = [rootpath,'20240308\DialS5S6S7.csv'];
path_20240308 = [rootpath,'20240308\DialS5S6S7.csv'];
path_KIT = [rootpath,'KIT_FilledSamples\AbsorptionFilled_KIT.csv'];

%Select the paths of interest
paths = {
        path_baselines,
        path_20231206,
        path_20240202,
        path_20240212a,
        path_20240214,
        path_20240216,
        path_20240304,
        path_20240305,
        path_20240307,
        path_20240308,
        path_KIT
        };
    

%Read and structure data from the paths
ReadAbsorptionFromPaths(paths);


%%%--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Data Reference
DATA_References.empty_new_P2_0601.N = 'SC Empty@SWCNT A';
DATA_References.empty_P2_0329.N = 'SC Empty@SWCNT B';
DATA_References.empty_P2_0601.N = 'SC Empty@SWCNT C';
DATA_References.empty_P2_0706.N = 'SC Empty@SWCNT D';
DATA_References.empty_P2_APA218.N = 'SC Empty@SWCNT E';
DATA_References.empty_P2_dial_0930.N = 'SC Empty@SWCNT Dial';
DATA_References.AnnoctadecaneP2.N = 'SC Annealed Octadecane@SWCNTP2';
DATA_References.AnnoctadecaneP2650.N = 'SC Annealed Octadecane@SWCNT650';
DATA_References.OctadecaneP2.N = 'SC OctadecaneP2@SWCNTP2';
DATA_References.Octadecene1P2.N = 'SC Octadecne1P2@SWCNTP2';
DATA_References.TriacontaneP2.N = 'SC Triacontane@SWCNTP2';
DATA_References.WaterFilled.N = 'SC Water@SWCNT';

DATA_20240202.BASELINE.N='Baseline (DOC/D2O)';
DATA_20240214.BASELINE.N='Baseline 1%DOCD2O';
DATA_20240304.S7.N='CB Dodecane@SWCNT';
DATA_20240307.S7CF.N='CB Dodecane@SWCNT 4hCF';
DATA_20240307.S7DGUC.N='CB Dodecane@SWCNT DGU C (Filled)';
DATA_20240308.S7DGUCDial.N='CB Dodecane@SWCNT Dial. DGU C (Filled)';
DATA_20240307.S5DGUB.N='CB Empty@SWCNT DGU B S5';
DATA_20240307.S6DGUB.N='CB Empty@SWCNT DGU B S6';
DATA_20240307.S7DGUB.N='CB Empty@SWCNT DGU B S7';
DATA_20240304.S6.N='CB Hexadecane@SWCNT';
DATA_20240307.S6CF.N='CB Hexadecane@SWCNT 4hCF';
DATA_20240307.S6DGUC.N='CB Hexadecane@SWCNT DGU C (Filled)';
DATA_20240308.S6DGUCDial.N='CB Hexadecane@SWCNT Dial. DGU C (Filled)';
DATA_20240124.BASELINE.N='CB Methanol Baseline 240115';
DATA_20240129.BASELINE.N='CB Methanol Baseline 240115';
DATA_20240202.S2_1060.N='CB PCE@SWCNT 2hCF';
DATA_20240212.S2.N='CB PCE@SWCNT 4hCF';
DATA_20240214.S2DB.N='CB PCE@SWCNT DGU C (Empty)';
DATA_20240214.S2DC.N='CB PCE@SWCNT DGU C (Filled)';
DATA_20240214.S2DA.N='CB PCE@SWCNT DGU C (Nicondenz)';
DATA_20240214.S2DD.N='CB PCE@SWCNT DGU D (Defective)';
DATA_20240214.S2DE.N='CB PCE@SWCNT DGU E (Bundles)';
DATA_20240216.S2DGUC.N='CB PCE@SWCNT Dial. DGU C (Filled';
DATA_20240202.S3_1060.N='CB TCE@SWCNT 2hCF';
DATA_20240212.S3.N='CB TCE@SWCNT 4hCF';
DATA_20240214.S3DB.N='CB TCE@SWCNT DGU C (Empty)';
DATA_20240214.S3DC.N='CB TCE@SWCNT DGU C (Filled)';
DATA_20240214.S3DA.N='CB TCE@SWCNT DGU C (Nicondenz)';
DATA_20240214.S3DD.N='CB TCE@SWCNT DGU D (Defective)';
DATA_20240214.S3DE.N='CB TCE@SWCNT DGU E (Bundles)';
DATA_20240216.S3DGUC.N='CB TCE@SWCNT Dial. DGU C (Filled)';
DATA_20240305.S5.N='CB TDAE@SWCNT';
DATA_20240305.S5_DIL.N='CB TDAE@SWCNT';
DATA_20240305.S5_DIL2.N='CB TDAE@SWCNT';
DATA_20240305.S5_DIL4.N='CB TDAE@SWCNT';
DATA_20240307.S5CF.N='CB TDAE@SWCNT 4hCF';
DATA_20240307.S5DGUC.N='CB TDAE@SWCNT DGU C (Filled)';
DATA_20240308.S5DGUCDial.N='CB TDAE@SWCNT Dial. DGU C (Filled)';
DATA_20240202.S4_1060.N='CB TEMED@SWCNT 2hCF';
DATA_20240212.S4.N='CB TEMED@SWCNT 4hCF';
DATA_20240214.S4DB.N='CB TEMED@SWCNT DGU C (Empty)';
DATA_20240214.S4LB.N='CB TEMED@SWCNT DGU C (Empty)';
DATA_20240214.S4DC.N='CB TEMED@SWCNT DGU C (Filled)';
DATA_20240214.S4LC.N='CB TEMED@SWCNT DGU C (Filled)';
DATA_20240214.S4DA.N='CB TEMED@SWCNT DGU C (Nicondenz)';
DATA_20240214.S4LA.N='CB TEMED@SWCNT DGU C (Nicondenz)';
DATA_20240214.S4DD.N='CB TEMED@SWCNT DGU D (Defective)';
DATA_20240214.S4LD.N='CB TEMED@SWCNT DGU D (Defective)';
DATA_20240214.S4DE.N='CB TEMED@SWCNT DGU E (Bundles)';
DATA_20240214.S4LE.N='CB TEMED@SWCNT DGU E (Bundles)';
DATA_20240216.S4DDGUC.N='CB TEMED@SWCNT Dial. DGU C (Filled)';
DATA_20240216.S4DDGUC2.N='CB TEMED@SWCNT Dial. DGU C (Filled)';
DATA_20240216.S4LDGUC.N='CB TEMED@SWCNT Dial. DGU C (Filled)';
DATA_20240216.BASELINE.N='DOCD2O Baseline';
DATA_20240304.BASELINE.N='DOCD2O Baseline';
DATA_20240305.BASELINE.N='DOCD2O Baseline';
DATA_20240307.BASELINE.N='DOCD2O Baseline';
DATA_20240308.BASELINE.N='DOCD2O Baseline';
DATA_20240220.BASELINE.N='EthylAcetate Baseline';
DATA_20240308.S5D1.N='Nicodenz after D(1)';
DATA_20240308.S5D4.N='Nicodenz after D(4)';
DATA_20240308.S5D7.N='Nicodenz after D(7)';
DATA_20240308.S5D9.N='Nicodenz after D(9)';

DATA_20231206.SFF4dil.N='SF Methanol@SWCNT A';
DATA_20231206.SFF4dil2.N='SF Methanol@SWCNT B';
DATA_20231206.SFF2dil.N='SF PCE@SWCNT';
DATA_20231206.SFF2Bdil.N='SF PCE@SWCNT dil1';
DATA_20231206.SFF5dil.N='SF SFF5 Water@SWCNT';
DATA_20231206.SFF3dil.N='SF TCE@SWCNT LP Long Rinsing';
DATA_20231206.SFF3Bdil.N='SF TCE@SWCNT LP Short Rinsing';
DATA_20231206.SFF3_3dil.N='SF TCE@SWCNT RF';
DATA_20231206.SFF1Bdil.N='SF TEMED@SWCNT';
DATA_20231206.SFF6dil.N='SF TTF@SWCNT';

%%%--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Compare before and after H2O in D2O Correction
%SampleToCorrect = DATA_20240307.S7DGUC
%H2OD2OFactor = 0.7
%NicodenzFactor = 1.2

%Corrected = SampleToCorrect;
%Corrected.N = 'Corrected';
%Corrected.A = SampleToCorrect.A - DATA_References.H2OinD2O.A*H2OD2OFactor - DATA_References.Nicodenz.A*NicodenzFactor
%plotSampleList({Corrected, SampleToCorrect},0.0)
%plotSampleList({DATA_References.H2OinD2O, DATA_References.Nicodenz},0.0)

%-----------------------------------------------

%%%--------SAMPLE COMPARISION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
DialSamples = {
            DATA_References.empty_P2_dial_0930,
            DATA_References.WaterFilled

            DATA_20240308.S6DGUCDial,
            DATA_20240308.S7DGUCDial,
            DATA_20240216.S2DGUC,
            DATA_20240216.S4DDGUC,
            DATA_20240216.S3DGUC,
            DATA_20240308.S5DGUCDial
            };
          
BenSamp = {
           DATA_KIT_FilledSamples.C12H26SWCNT_stage1_UF1DOC2
           DATA_20240308.S7DGUCDial

           DATA_KIT_FilledSamples.C16H38SWCNT_stage1_UF1DOC
           DATA_20240308.S6DGUCDial

           DATA_KIT_FilledSamples.PCESWCNT_stage1_UF1DOC
           DATA_20240216.S2DGUC
           
           DATA_KIT_FilledSamples.TcESWCNT_stage1_UF1DOC
           DATA_20240216.S3DGUC
           
           DATA_KIT_FilledSamples.TDAESWCNT_stage1_UF1DOC1
           DATA_20240308.S5DGUCDial
           
           DATA_KIT_FilledSamples.TEMEDSWCNT_stage1_UF1DOC1
           DATA_20240216.S4DDGUC
            };
          
       
        
%%%--------BACKGROUND CORRECTION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Correction is made based on local minima (750, 860) and (1100, 1400)
DialSamples = SubstractAbsBG(DialSamples,750, 850, 1150, 1250);

%%%--------NORMALIZATION TO S22 PEAK--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Normalization using maximum value of S22 transition, can also use Integral
LS2= 900;  
US2= 1100;
%600a800

DialSamples = NormalizeSample(DialSamples,LS2, US2);
BenSamp = NormalizeSample(BenSamp,LS2, US2);
plotAbsorption(BenSamp,0.0)


%%%--------PEAK CALCULATION AND EXPORT--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Find peaks in these ranges
%S11
LS1= 1600;
US1= 1850;

%Plasmon Right
LPR= 260;
UPR= 285;
%Plasmon Left
LPL= 220;
UPL= 250;



%plotAbsorption(Alkanes,1.0)
%plotAbsorption(Dopants, 0.0)

%plotAbsorption(PCE, 0.0)
%plotAbsorption(TCE, 0.0)
%plotAbsorption(TDAE, 0.0)
%plotAbsorption(TEMED, 0.0)
%plotAbsorption(Hexadecane, 0.0)
%plotAbsorption(Dodecane, 0.0)

%plotAbsorption(DialSamples,1.0)


%plotAbsorption(DialSamples, 1.5)

%plotAbsorption(DialSamples, 0.0)

%plotAbsorption([DialSamples;DATABEN], 0.0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



