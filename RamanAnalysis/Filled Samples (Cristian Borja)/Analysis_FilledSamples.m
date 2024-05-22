clc;
clear;
addpath('C:\Users\cborja\OneDrive - Universiteit Antwerpen\SWCNTs\');
%addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\');
import UsefulFunctions.*;

%rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements Data\Raman\';
rootpath = 'C:\Users\cborja\OneDrive - Universiteit Antwerpen\Measurements Data\Raman\';

%All paths as default
path_20240320 = [rootpath,'20240320\'];
path_20240321 = [rootpath,'20240321\'];
path_20240514 = [rootpath,'20240514\'];
path_20240517 = [rootpath,'20240517\'];

%Select the paths of interest

paths = {
    path_20240320
    path_20240321
    path_20240514
    path_20240517
    };


ReadRamanFromPaths(paths);


%%%--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DATA_20240514.BBL570RB.N = 'S-SWCNTs Converted Film RBMs2 570nm';

%FILLED SAMPLES AFTER DIALISYS AT 650nm + WaterFilled/Empty References

DATA_20240517.EAH514R.N='Empty Arc SWCNTs';
DATA_20240517.S2H514R.N='CB PCE@SWCNT Dial. DGU C (Filled)';
DATA_20240517.S3H514R.N='CB TCE@SWCNT Dial. DGU C (Filled)';
DATA_20240517.S4H514R.N='CB TEMED@SWCNT Dial. DGU C (Filled)';
DATA_20240517.S5H514R.N='CB TDAE@SWCNT Dial. DGU C (Filled)';
DATA_20240517.S6H514R.N='CB Hexadecane@SWCNT Dial. DGU C (Filled)';
DATA_20240517.S7H514R.N='CB Dodecane@SWCNT Dial. DGU C (Filled)';
DATA_20240517.WAH514R.N='Water Filled Arc SWCNTs';
DATA_20240517.EAL514GD.N='Empty Arc SWCNTs';
DATA_20240517.S2L514GD.N='CB PCE@SWCNT Dial. DGU C (Filled)';
DATA_20240517.S3L514GD.N='CB TCE@SWCNT Dial. DGU C (Filled)';
DATA_20240517.S4L514GD.N='CB TEMED@SWCNT Dial. DGU C (Filled)';
DATA_20240517.S5L514GD.N='CB TDAE@SWCNT Dial. DGU C (Filled)';
DATA_20240517.S6L514GD.N='CB Hexadecane@SWCNT Dial. DGU C (Filled)';
DATA_20240517.S7L514GD.N='CB Dodecane@SWCNT Dial. DGU C (Filled)';
DATA_20240517.WAL514GD.N='Water Filled Arc SWCNTs';
DATA_20240515.LL514A.N='LaserLine';
DATA_20240515.LL514B.N='LaserLine';

DATA_20240514.EAL650D.N = 'EmptyArc@SWCNTs DBand 650nm';
DATA_20240514.EAL650G.N = 'EmptyArc@SWCNTs GBand 650nm';
DATA_20240514.EAL650R.N = 'EmptyArc@SWCNTs RBMs 650nm';
DATA_20240514.WAL650D.N = 'Water@ArcSWCNTs DBand 650nm';
DATA_20240514.WAL650G.N = 'Water@ArcSWCNTs GBand 650nm';
DATA_20240514.WAL650R.N = 'Water@ArcSWCNTs RBMs 650nm';
DATA_20240514.S2L650D.N = 'S2 PCE@SWCNTs DBand 650nm';
DATA_20240514.S2L650G.N = 'S2 PCE@SWCNTs GBand 650nm';
DATA_20240514.S2L650R.N = 'S2 PCE@SWCNTs RBMs 650nm';
DATA_20240514.S3L650D.N = 'S3 TCE@SWCNTs DBand 650nm';
DATA_20240514.S3L650G.N = 'S3 TCE@SWCNTs GBand 650nm';
DATA_20240514.S3L650R.N = 'S3 TCE@SWCNTs RBMs 650nm';
DATA_20240514.S4L650D.N = 'S4 TEMED@SWCNTs DBand 650nm';
DATA_20240514.S4L650G.N = 'S4 TEMED@SWCNTs GBand 650nm';
DATA_20240514.S4L650R.N = 'S4 TEMED@SWCNTs RBMs 650nm';
DATA_20240514.S5L650D.N = 'S5 TDAE@SWCNTs DBand 650nm';
DATA_20240514.S5L650G.N = 'S5 TDAE@SWCNTs GBand 650nm';
DATA_20240514.S5L650R.N = 'S5 TDAE@SWCNTs RBMs 650nm';
DATA_20240514.S6L650D.N = 'S6 Hexadecane@SWCNTs DBand 650nm';
DATA_20240514.S6L650G.N = 'S6 Hexadecane@SWCNTs GBand 650nm';
DATA_20240514.S6L650R.N = 'S6 Hexadecane@SWCNTs RBMs 650nm';
DATA_20240514.S7L650D.N = 'S7 Dodecane@SWCNTs DBand 650nm';
DATA_20240514.S7L650G.N = 'S7 Dodecane@SWCNTs GBand 650nm';
DATA_20240514.S7L650R.N = 'S7 Dodecane@SWCNTs RBMs 650nm';


DATA_20240320.L520A.N='LaserLine';
DATA_20240320.FFH520R.N='FlatField';
DATA_20240320.S2H520R.N='CB PCE@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240320.S3H520R.N='CB TCE@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240320.S4H520R.N='CB TEMED@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240320.S5H520R.N='CB TDAE@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240320.S6H520R.N='CB Hexadecane@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240320.S7H520R.N='CB Dodecane@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240321.L520HG.N='LaserLine 520nm';
DATA_20240321.L520L.N='LaserLine 520nm';
DATA_20240321.FFL520G.N='FlatField 520nm';
DATA_20240321.S2L520G.N='CB PCE@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240321.S3L520G.N='CB TCE@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240321.S4L520G.N='CB TEMED@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240321.S5L520G.N='CB TDAE@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240321.S6L520G.N='CB Hexadecane@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240321.S7L520G.N='CB Dodecane@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240321.FFH520G.N='FlatField 520nm';
DATA_20240321.LLH520GB.N='LaserLine 520nm';
DATA_20240321.S2H520G.N='CB PCE@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240321.S3H520G.N='CB TCE@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240321.S4H520G.N='CB TEMED@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240321.S5H520G.N='CB TDAE@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240321.S6H520G.N='CB Hexadecane@SWCNT Dial. DGU C (Filled) 520nm';
DATA_20240321.S7H520G.N='CB Dodecane@SWCNT Dial. DGU C (Filled) 520nm';

%%%--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DATA_20240426.BAL650R.Y = DATA_20240426.BAL650R.Y/2;
%DATA_20240426.BAL650RB.Y = DATA_20240426.BAL650RB.Y/2;

%%%--------SAMPLE COMPARISION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 650nm DATA 
FilledSamples650 = {
                    DATA_20240514.EAL650D
                    DATA_20240514.EAL650G
                    DATA_20240514.EAL650R
                    DATA_20240514.WAL650D
                    DATA_20240514.WAL650G
                    DATA_20240514.WAL650R
                    DATA_20240514.S2L650D
                    DATA_20240514.S2L650G
                    DATA_20240514.S2L650R
                    DATA_20240514.S3L650D
                    DATA_20240514.S3L650G
                    DATA_20240514.S3L650R
                    DATA_20240514.S4L650D
                    DATA_20240514.S4L650G
                    DATA_20240514.S4L650R
                    DATA_20240514.S5L650D
                    DATA_20240514.S5L650G
                    DATA_20240514.S5L650R
                    DATA_20240514.S6L650D
                    DATA_20240514.S6L650G
                    DATA_20240514.S6L650R
                    DATA_20240514.S7L650D
                    DATA_20240514.S7L650G
                    DATA_20240514.S7L650R
                   };
               
GBands650 =    {
            DATA_20240514.EAL650G
            DATA_20240514.WAL650G
            DATA_20240514.S2L650G
            DATA_20240514.S3L650G
            DATA_20240514.S4L650G
            DATA_20240514.S5L650G
            DATA_20240514.S6L650G
            DATA_20240514.S7L650G
            };
        
DBands650 =    {
            DATA_20240514.EAL650D
            DATA_20240514.WAL650D
            DATA_20240514.S2L650D
            DATA_20240514.S3L650D
            DATA_20240514.S4L650D
            DATA_20240514.S5L650D
            DATA_20240514.S6L650D
            DATA_20240514.S7L650D
            };
        
RBMs650 =    {
            DATA_20240514.EAL650R
            DATA_20240514.WAL650R
            DATA_20240514.S2L650R
            DATA_20240514.S3L650R
            DATA_20240514.S4L650R
            DATA_20240514.S5L650R
            DATA_20240514.S6L650R
            DATA_20240514.S7L650R
            };
        

% %Normalization
LG = 1450;
HG = 1640;
NP = 1588;
Tol = 5;
GBands650 = SubstractLinearBG(GBands650, LG, HG);
GBands650 = NormalizeSample(GBands650,NP-Tol, NP+Tol);       

% %Normalization
LG = 1240;
HG = 1400;
NP = 1315;
Tol = 5;
DBands650 = SubstractLinearBG(DBands650, LG, HG);
DBands650 = NormalizeSample(DBands650,NP-Tol, NP+Tol);    

%Normalization
LG = 120;
HG = 300;
NP = 175;
Tol = 25;
RBMs650 = SubstractLinearBG(RBMs650, LG, HG);
RBMs650 = NormalizeSample(RBMs650,NP-Tol, NP+Tol);       
%         
%plotRaman([GBands650; DBands650], 0)
% plotRaman(GBands650, 0.0)
% plotRaman(DBands650, 0.0)
% plotRaman(RBMs650, 0.25)


%% 520nm DATA 

FilledSamples520 = {
                    DATA_20240320.S2H520R
                    DATA_20240320.S3H520R
                    DATA_20240320.S4H520R
                    DATA_20240320.S5H520R
                    DATA_20240320.S6H520R
                    DATA_20240320.S7H520R
                    DATA_20240321.S2L520G
                    DATA_20240321.S3L520G
                    DATA_20240321.S4L520G
                    DATA_20240321.S5L520G
                    DATA_20240321.S6L520G
                    DATA_20240321.S7L520G
                    DATA_20240321.S2H520G
                    DATA_20240321.S3H520G
                    DATA_20240321.S4H520G
                    DATA_20240321.S5H520G
                    DATA_20240321.S6H520G
                    DATA_20240321.S7H520G
                   };
               
RBMs520 = {
                    DATA_20240320.S2H520R
                    DATA_20240320.S3H520R
                    DATA_20240320.S4H520R
                    DATA_20240320.S5H520R
                    DATA_20240320.S6H520R
                    DATA_20240320.S7H520R
            };   
        
               
GBand520HD = {
                    DATA_20240321.S2H520G
                    DATA_20240321.S3H520G
                    DATA_20240321.S4H520G
                    DATA_20240321.S5H520G
                    DATA_20240321.S6H520G
                    DATA_20240321.S7H520G
            };   
        
GBand520LD = {
                    DATA_20240321.S2L520G
                    DATA_20240321.S3L520G
                    DATA_20240321.S4L520G
                    DATA_20240321.S5L520G
                    DATA_20240321.S6L520G
                    DATA_20240321.S7L520G
            };   
        

% %Normalization
LG = 1280;
HG = 1450;
NP = 1340;
Tol = 10;
GBand520LD = SubstractLinearBG(GBand520LD, LG, HG);
GBand520LD = NormalizeSample(GBand520LD,NP-Tol, NP+Tol);       

% %Normalization
LG = 1536;
HG = 1618;
NP = 1565;
Tol = 5;
GBand520HD = SubstractLinearBG(GBand520HD, LG, HG);
GBand520HD = NormalizeSample(GBand520HD,NP-Tol, NP+Tol);       

%Normalization
LG = 130;
HG = 220;
NP = 175;
Tol = 5;
RBMs520 = SubstractLinearBG(RBMs520, LG, HG);
RBMs520 = NormalizeSample(RBMs520,NP-Tol, NP+Tol);       

%plotRaman(FilledSamples520, 0)

% plotRaman(RBMs520, 0.25)
% plotRaman(GBand520HD, 0.0)
% plotRaman(GBand520LD, 0.0)

%% 514nm DATA 

FilledSamples514 = {
%                     DATA_20240517.EAH514R
%                     DATA_20240517.EAL514GD
%                      DATA_20240517.FF514R
%                      DATA_20240517.LL514A
%                      DATA_20240517.LL514B
%                     DATA_20240517.S2H514R
%                     DATA_20240517.S2L514GD
%                     DATA_20240517.S3H514R
%                     DATA_20240517.S3L514GD
%                     DATA_20240517.S4H514R
%                     DATA_20240517.S4L514GD
%                     DATA_20240517.S5H514R
%                     DATA_20240517.S5L514GD
%                     DATA_20240517.S6H514R
%                     DATA_20240517.S7L514GD
%                     DATA_20240517.S7H514R
%                     DATA_20240517.WAH514R
%                     DATA_20240517.WAL514GD
                   };
               
RBMsHD514 = {
                    DATA_20240517.EAH514R
                    DATA_20240517.S2H514R
                    DATA_20240517.S3H514R
                    DATA_20240517.S4H514R
                    DATA_20240517.S5H514R
                    DATA_20240517.S6H514R
                    DATA_20240517.S7H514R
                    DATA_20240517.WAH514R
            };   
        
               
GDBand514 = {
                    DATA_20240517.EAL514GD
                    DATA_20240517.S2L514GD
                    DATA_20240517.S3L514GD
                    DATA_20240517.S4L514GD
                    DATA_20240517.S5L514GD
                    DATA_20240517.S7L514GD
                    DATA_20240517.WAL514GD
            };   

% %Normalization
LG = 1280;
HG = 1490;
NP = 1590;
Tol = 5;
GDBand514 = SubstractLinearBG(GDBand514, LG, HG);
GDBand514 = NormalizeSample(GDBand514,NP-Tol, NP+Tol);       

%Normalization
LG = 140;
HG = 250;
NP = 170;
Tol = 5;
RBMsHD514 = FlatFieldCorrection(RBMsHD514, DATA_20240517.FF514R);
RBMsHD514 = SubstractLinearBG(RBMsHD514, LG, HG);
RBMsHD514 = NormalizeSample(RBMsHD514,NP-Tol, NP+Tol);       


%plotRaman(FilledSamples514, 0.25)
% plotRaman(RBMsHD514, 0.25)
% plotRaman(GDBand514, 0.25)


%% Peak Calculation
%GDBand520 = GDBandPeaksCalculation(GDBand520, 1580,1600,1560,1570,1300,1400);
%exportGDBandPeaks(GDBand520, 'dsfsd.csv');