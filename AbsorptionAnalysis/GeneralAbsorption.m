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

path_20240415a = [rootpath,'20240415\Acetone.csv'];
path_20240415b = [rootpath,'20240415\Ethanol.csv'];
path_20240415c = [rootpath,'20240415\Methanol.csv'];
path_20240415d = [rootpath,'20240415\THF.csv'];

path_20240422a = [rootpath,'20240422\EthylAcetate.csv'];
path_20240422b = [rootpath,'20240422\Ethanol.csv'];
path_20240422c = [rootpath,'20240422\Methanol.csv'];
path_20240422d = [rootpath,'20240422\THF.csv'];

path_20240506a = [rootpath,'20240506\EthylAcetate.csv'];
path_20240506b = [rootpath,'20240506\Acetone.csv'];
path_20240506c = [rootpath,'20240506\Acetonitrile.csv'];
path_20240506d = [rootpath,'20240506\THF.csv'];
path_20240506e = [rootpath,'20240506\DMF.csv'];

path_benupdate = [rootpath,'Data.csv']

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
        path_20240415a,
        path_20240415b,
        path_20240415c,
        path_20240415d,
        path_20240422a,
        path_20240422b,
        path_20240422c,
        path_20240422d,
        path_benupdate,
        path_20240506a,
        path_20240506b,
        path_20240506c,
        path_20240506d,
        path_20240506e
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

%H2O in D2O Correction fixed corrections.
%SAMPLE 2
DATA_20240216.S2DGUC.Y = DATA_20240216.S2DGUC.Y - DATA_References.H2OinD2O.Y*1.22;
%SAMPLE 3
DATA_20240216.S3DGUC.Y = DATA_20240216.S3DGUC.Y - DATA_References.H2OinD2O.Y*0.7;
%SAMPLE 4
DATA_20240216.S4DDGUC.Y = DATA_20240216.S4DDGUC.Y - DATA_References.H2OinD2O.Y*0.25;
%SAMPLE 5 Doesnt need correction
%SAMPLE 6 Doesnt need correction
%SAMPLE 7 Doesnt need correction

%%%--------SAMPLE COMPARISION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AllData = {
    %DATA_References.empty_new_P2_0601,
    %DATA_References.empty_P2_0329,
    %DATA_References.empty_P2_0601,
    %DATA_References.empty_P2_0706,
    %DATA_References.empty_P2_APA218,
    %DATA_References.empty_P2_dial_0930,
    
    %Alkanes
    %DATA_References.TriacontaneP2,
    %DATA_References.AnnoctadecaneP2,
    %DATA_References.AnnoctadecaneP2650,  
    %DATA_20240308.S7DGUCDial,
    %DATA_20240308.S6DGUCDial,
    
    %Dopants
    
    %PCE
    %DATA_20231206.SFF2Bdil,
    %DATA_20231206.SFF2dil,
    %DATA_20240216.S2DGUC,

    %TCE
    %DATA_20231206.SFF3_3dil,
    %DATA_20231206.SFF3dil,
    %DATA_20231206.SFF3Bdil,
    %DATA_20240216.S3DGUC,
   
    %TEMED
    %DATA_20240216.S4DDGUC,
    %DATA_20231206.SFF1Bdil,
    
    %TDAE
    %DATA_20240308.S5DGUCDial,
    
    %DATA_References.WaterFilled
    %DATA_20231206.SFF5dil,
    };


Alkanes = {
    %Empty
    DATA_References.empty_P2_dial_0930,
    
    %Alkanes
    DATA_References.TriacontaneP2,
    DATA_References.AnnoctadecaneP2,
    DATA_References.AnnoctadecaneP2650,  
    DATA_20240308.S7DGUCDial,
    DATA_20240308.S6DGUCDial,
    
    %Water
    DATA_References.WaterFilled
    };

PCE = {
    
    %Empty
    %DATA_References.empty_P2_dial_0930,
    
    %DATA_20231206.SFF2dil,
    %DATA_20231206.SFF2Bdil,
    DATA_20240202.S2_1060,
    DATA_20240212.S2,
    DATA_20240216.S2DGUC,
    
    %Water
    %DATA_References.WaterFilled
}

TCE = {
    
    %Empty
    %DATA_References.empty_P2_dial_0930,
    
    %DATA_20231206.SFF3dil,
    %DATA_20231206.SFF3Bdil,
    %DATA_20231206.SFF3_3dil,
    DATA_20240202.S3_1060,
    DATA_20240212.S3,
    DATA_20240216.S3DGUC,    
    %Water
    %DATA_References.WaterFilled
}

TEMED = {
    
    %Empty
    %DATA_References.empty_P2_dial_0930,
    
    %DATA_20231206.SFF1Bdil,
    DATA_20240202.S4_1060,
    DATA_20240212.S4,
    DATA_20240216.S4DDGUC,
    %DATA_20240216.S4LDGUC,
    %Water
    %DATA_References.WaterFilled
}

TDAE = {
    
    %Empty
    %DATA_References.empty_P2_dial_0930,

    DATA_20240305.S5,
    %DATA_20240305.S5_dil4,
    DATA_20240307.S5CF,
    DATA_20240308.S5DGUCDial,
    
    %Water
    %DATA_References.WaterFilled
}

Hexadecane = {
    
    %Empty
    %DATA_References.empty_P2_dial_0930,

    DATA_20240304.S6,
    DATA_20240307.S6CF,
    DATA_20240308.S6DGUCDial,
    
    %Water
    %DATA_References.WaterFilled
}

Dodecane = {
    
    %Empty
    %DATA_References.empty_P2_dial_0930,

    DATA_20240304.S7,
    DATA_20240307.S7CF,
    DATA_20240308.S7DGUCDial,

    %Water
    %DATA_References.WaterFilled
}

Dopants = {
    %Empty
    DATA_References.empty_P2_dial_0930,
    
    %PDopants
    %PCE
    DATA_20231206.SFF2Bdil,
    DATA_20231206.SFF2dil,
    DATA_20240216.S2DGUC,
    %TCE
    DATA_20231206.SFF3_3dil,
    DATA_20231206.SFF3dil,
    DATA_20231206.SFF3Bdil,
    DATA_20240216.S3DGUC,
    
    %NDopants
    %TEMED
    DATA_20240216.S4DDGUC,
    DATA_20231206.SFF1Bdil,
    
    %TDAE
    DATA_20240308.S5DGUCDial,
    
    %Water
    DATA_References.WaterFilled
    };

CFSamples = {
            %DATA_References.empty_P2_dial_0930,

            DATA_20240212.S2,
            %DATA_20231206.SFF2Bdil,
            %DATA_20231206.SFF2dil,
            
            DATA_20240212.S3,
            %DATA_20231206.SFF3_3dil,
            %DATA_20231206.SFF3dil,
            %DATA_20231206.SFF3Bdil,
    
            DATA_20240212.S4,
            %DATA_20231206.SFF1Bdil,

            %DATA_References.WaterFilled

            };

        
DialSamples = {
            %DATA_References.empty_P2_dial_0930,
            %DATA_References.WaterFilled

            DATA_20240308.S6DGUCDial,
            DATA_20240308.S7DGUCDial,
            DATA_20240216.S2DGUC,
            DATA_20240216.S4DDGUC,
            DATA_20240216.S3DGUC,
            DATA_20240308.S5DGUCDial
            

            };

DATABEN = { DATA_Data.C16H38SWCNT_stage1_UF1DOC,
            DATA_Data.C12H26SWCNT_stage1_UF1DOC2,
            DATA_Data.PCESWCNT_stage1_UF1DOC,
            DATA_Data.TEMEDSWCNT_stage1_UF1DOC1,
            DATA_Data.TcESWCNT_stage1_UF1DOC,
            DATA_Data.TDAESWCNT_stage1_UF1DOC1,
            
            };

        
        
%%%--------BACKGROUND CORRECTION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Correction is made based on local minima (750, 860) and (1100, 1400)
DialSamples = SubstractAbsBG(DialSamples,750, 850, 1150, 1250);
DATABEN = SubstractAbsBG(DATABEN,750, 850, 1150, 1250);

Alkanes = SubstractAbsBG(Alkanes,750, 850, 1150, 1250);
%CFSamples = SubstractAbsBG(CFSamples,750, 850, 1150, 1250);
%PCE = SubstractAbsBG(PCE,750, 850, 1150, 1250);
%TCE = SubstractAbsBG(TCE,750, 850, 1150, 1250);
%TEMED = SubstractAbsBG(TEMED,750, 850, 1150, 1250);
%TDAE = SubstractAbsBG(TDAE,750, 850, 1150, 1250);
%Hexadecane = SubstractAbsBG(Hexadecane,750, 850, 1150, 1250);
%Dodecane = SubstractAbsBG(Dodecane,750, 850, 1150, 1250);



%%%--------NORMALIZATION TO S22 PEAK--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Normalization using maximum value of S22 transition, can also use Integral
LS2= 900;  
US2= 1100;

DialSamples = NormalizeSample(DialSamples,LS2, US2);
DATABEN = NormalizeSample(DATABEN,LS2, US2);

Alkanes = NormalizeSample(Alkanes,LS2, US2);
%CFSamples= NormalizeSample(CFSamples,LS2, US2);
%PCE = NormalizeSample(PCE,LS2, US2);
%TCE = NormalizeSample(TCE,LS2, US2);
%TEMED = NormalizeSample(TEMED,LS2, US2);
%TDAE = NormalizeSample(TDAE,LS2, US2);
%Hexadecane = NormalizeSample(Hexadecane,LS2, US2);
%Dodecane = NormalizeSample(Dodecane,LS2, US2);

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

DialSamples = TransitionPeaksCalculation(DialSamples, LS1, US1, LS2, US2);
%exportPeaksToCSV(DialSamples, 'DialSamples.csv')

CFSamples = TransitionPeaksCalculation(CFSamples, LS1, US1, LS2, US2);
%exportPeaksToCSV(CFSamples, 'CFSamples.csv')

%%%%PLOTING

%plotAbsorption(Alkanes,1.0)
%plotAbsorption(Dopants, 0.0)

%plotAbsorption(PCE, 0.0)
%plotAbsorption(TCE, 0.0)
%plotAbsorption(TDAE, 0.0)
%plotAbsorption(TEMED, 0.0)
%plotAbsorption(Hexadecane, 0.0)
%plotAbsorption(Dodecane, 0.0)

%plotAbsorption(DialSamples,1.0)

Baselines = {   DATA_20240415.Acetone_Baseline,
                DATA_20240415.Ethanol_Baseline,
                DATA_20240422.Ethanol_Baseline,
                DATA_20240415.Methanol_Baseline,
                DATA_20240422.Methanol_Baseline,
                DATA_20240415.THF_Baseline,
                DATA_20240422.THF_Baseline,
                DATA_20240422.EthylAcetate_Baseline
                };
            


SolubilitySolvents = {
            %DATA_20240415.Acetone_Baseline
            
            DATA_20240506.AcetoneBaseline
            DATA_20240415.TCBQ_Acetone,
            DATA_20240415.DCB_Acetone,  
            DATA_20240506.CuCl_Acetone,
            DATA_20240506.FN_Acetone,
            DATA_20240506.TCNQ_Acetone,
            
            %DATA_20240415.Ethanol_Baseline,
            %DATA_20240422.Ethanol_Baseline
            %DATA_20240415.TCBQ_Ethanol,
            %DATA_20240422.TCBQ_Ethanol
            %DATA_20240415.FN_Ethanol,
            %DATA_20240422.FN_Ethanol
           
            %DATA_20240415.Methanol_Baseline,
            %DATA_20240422.Methanol_Baseline,
            %DATA_20240415.FN_Methanol,
            %DATA_20240422.FN_Methanol,
            
            %DATA_20240415.THF_Baseline,
            %DATA_20240422.THF_Baseline
            %DATA_20240415.DCB_THF,
            %DATA_20240422.DCB_THF,
            %DATA_20240415.TCNQ_THF
            %DATA_20240422.TCNQ_THF,
%             DATA_20240506.DCB_THF,
%             DATA_20240506.DCB_THF_B,
%             DATA_20240506.FN_THF,
%             DATA_20240506.CuCl_THF,
%             DATA_20240506.TCBQ_THF,
%             DATA_20240506.TCNQ_THF
            
            %DATA_20240422.EthylAcetate_Baseline,
            %DATA_20240422.TCNQ_EthylAcetate,
            %DATA_20240422.TCBQ_EthylAcetate,
            %DATA_20240422.FN_EthylAcetate,
            %DATA_20240422.DCB_EthylAcetate,
            
%             DATA_20240506.DCB_DMF,
%             DATA_20240506.FN_DMF,
%             DATA_20240506.TCBQ_DMF,
%             DATA_20240506.TCNQ_DMF,
%             DATA_20240506.CuCl_DMF,
            
            };
            
Solubility = {              
%               DATA_20240506.THFBaseline
%              DATA_20240506.DCB_Acetonitrile,
%              DATA_20240506.DCB_Acetonitrile_B,
%              DATA_20240506.DCB_DMF,
%              DATA_20240422.DCB_EthylAcetate,
%              DATA_20240506.DCB_EthylAcetate,
%               DATA_20240422.DCB_THF
%               DATA_20240506.DCB_THF,
%               DATA_20240506.DCB_THF_B,
            

%                DATA_20240506.FN_Acetone,
%               DATA_20240506.FN_Acetonitrile,
               %DATA_20240506.FN_Acetonitrile_B,
%             DATA_20240506.FN_DMF,
              %DATA_20240506.FN_THF,
%              DATA_20240422.FN_EthylAcetate
             
%             DATA_20240506.AcetonitrileBaseline
                  %DATA_20240506.TCBQ_Acetonitrile,
                  DATA_20240506.TCBQ_Acetonitrile_B,
                  DATA_20240506.TCBQ_DMF,
                 %DATA_20240506.TCBQ_EthylAcetate_A,
                 DATA_20240506.TCBQ_EthylAcetate_B,
                 DATA_20240506.TCBQ_THF,
                 DATA_20240415.TCBQ_Acetone

%             DATA_20240506.DMFBaseline
%             DATA_20240422.EthylAcetate_Baseline
%              DATA_20240506.TCNQ_Acetone,
%              DATA_20240506.TCNQ_Acetonitrile,
%              DATA_20240506.TCNQ_DMF,
%              DATA_20240506.TCNQ_EthylAcetate,
%              DATA_20240506.TCNQ_EthylAcetate_B,
%              DATA_20240506.TCNQ_THF,
%             
          %  DATA_20240506.THFBaseline
%             DATA_20240506.CuCl_Acetone,
%              DATA_20240506.CuCl_Acetonitrile,
%               DATA_20240506.CuCl_Acetonitrile_B,
%               DATA_20240506.CuCl_DMF,
%              DATA_20240506.CuCl_EthylAcetate,
%              DATA_20240506.CuCl_THF
            };    
        
%plotAbsorption(SolubilitySolvents, 2.0)

%Solubility = NormalizeSample(Solubility,290, 296);
plotAbsorption(Solubility, 0.0)

%plotAbsorption(DialSamples, 0.0)

%plotAbsorption([DialSamples;DATABEN], 0.0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



