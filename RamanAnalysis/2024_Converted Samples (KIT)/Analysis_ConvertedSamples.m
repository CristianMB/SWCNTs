clc;
clear;
% addpath('C:\Users\cborja\OneDrive - Universiteit Antwerpen\SWCNTs\');
%addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\');
addpath('X:\SWCNTs\');

import UsefulFunctions.*;

%rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements Data\Raman\';
%rootpath = 'C:\Users\cborja\OneDrive - Universiteit Antwerpen\Measurements Data\Raman\';
rootpath = 'X:\Measurements Data\Raman\';

%All paths as default
path_20240610 = [rootpath,'20240610\'];
path_20240612 = [rootpath,'20240612\'];
path_20240614 = [rootpath,'20240614\'];
path_20240620 = [rootpath,'20240620\'];
path_20240628 = [rootpath,'20240628\'];

%Select the paths of interest

paths = {
    path_20240610
    path_20240612
    path_20240614
    path_20240620
    path_20240628
    };


ReadRamanFromPaths(paths);


%%%--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DATA_20240610.EAH514R.N='Empty ArcSWCNTs';
DATA_20240610.EAL514GD.N='Empty ArcSWCNTs';
DATA_20240610.FFH514R.N='Flat';
DATA_20240610.LL514A.N='Laser';
DATA_20240610.LL514B.N='Laser';
DATA_20240610.S2H514R.N='PCE@P2-ASWCNTs';
DATA_20240610.S2L514GD.N='PCE@P2-ASWCNTs';
DATA_20240610.S3H514R.N='TCE@P2-ASWCNTs';
DATA_20240610.S3L514GD.N='TCE@P2-ASWCNTs';
DATA_20240610.S4H514R.N='TEMED@P2-ASWCNTs';
DATA_20240610.S4L514GD.N='TEMED@P2-ASWCNTs';
DATA_20240610.S5H514R.N='TDAE@P2-ASWCNTs';
DATA_20240610.S5L514GD.N='TDAE@P2-ASWCNTs';
DATA_20240610.S6H514R.N= 'Hexadecane@P2-ASWCNTs';
DATA_20240610.S6L514GD.N= 'Hexadecane@P2-ASWCNTs';
DATA_20240610.S7H514R.N='Dodecane@P2-ASWCNTs';
DATA_20240610.S7L514GD.N='Dodecane@P2-ASWCNTs';
DATA_20240610.WAH514R.N='Water ArcSWCNTs';
DATA_20240610.WAL514GD.N='Water ArcSWCNTs';

DATA_20240612.S2L514R.N='KIT PCE@SWCNT Dial';
DATA_20240612.S3L514R.N='KIT TCE@SWCNT Dial';
DATA_20240612.S4L514R.N='KIT TEMED@SWCNT Dial';
DATA_20240612.S5L514R.N='KIT TDAE@SWCNT Dial';
DATA_20240612.S6L514R.N='KIT Hexadecane@SWCNT Dial';
DATA_20240612.S7L514R.N='KIT Dodecane@SWCNT Dial';
DATA_20240612.LL514A.N='Laser';
DATA_20240612.WAL514R.N='Water Filled Arc SWCNTs';
DATA_20240612.EAL514R.N='Empty Arc SWCNTs';

DATA_20240612.S2L514GD.N='KIT PCE@SWCNT Dial';
DATA_20240612.S3L514GD.N='KIT TCE@SWCNT Dial';
DATA_20240612.S4L514GD.N='KIT TEMED@SWCNT Dial';
DATA_20240612.S5L514GD.N='KIT TDAE@SWCNT Dial';
DATA_20240612.S6L514GD.N='KIT Hexadecane@SWCNT Dial';
DATA_20240612.S7L514GD.N='KIT Dodecane@SWCNT Dial';
DATA_20240612.WAL514GD.N='Water Filled Arc SWCNTs';
DATA_20240612.EAL514GD.N='Empty Arc SWCNTs';

%Exchange for filled samples was made on 2024/06/11, exchange for converted
%solutions was made on 2024/06/14. So I only measured Raman before and
%after exchange just for the filled samples. I only measured the converted
%solutions after exchange 

DATA_20240614.LL514A.N='Laser';
DATA_20240614.T4SL514R.N='T4 S-SWCNTs Converted';
DATA_20240614.T1SL514R.N='T1 S-SWCNTs Converted';
DATA_20240614.T4PL514R.N='T4 P-SWCNTs Converted';
DATA_20240614.T4SL514G.N='T4 S-SWCNTs Converted';
DATA_20240614.T1SL514G.N='T1 S-SWCNTs Converted';
DATA_20240614.T4PL514G.N='T4 P-SWCNTs Converted';



DATA_20240620.LL514A.N='Laser';
DATA_20240620.F1SL514R.N='KIT Film T1S';
DATA_20240620.F1CL514R.N='KIT Film T1S Converted';
DATA_20240620.F2SL514R.N='KIT Film T2S';
DATA_20240620.F2CL514R.N='KIT Film T2S Converted';
DATA_20240620.F9ML514R.N='KIT Film T9M';
DATA_20240620.F9CL514R.N='KIT Film T9M Converted';
DATA_20240620.F1SL514G.N='KIT Film T1S';
DATA_20240620.F1CL514G.N='KIT Film T1S Converted';
DATA_20240620.F2SL514G.N='KIT Film T2S';
DATA_20240620.F2CL514G.N='KIT Film T2S Converted';
DATA_20240620.F9ML514G.N='KIT Film T9M';
DATA_20240620.F9CL514G.N='KIT Film T9M Converted';
DATA_20240620.F2SL514C.N='KIT Film T2S';



DATA_20240628.LL561A.N='Laser';
DATA_20240628.F1SL561R.N='KIT Film T1S';
DATA_20240628.F1SL561G.N='KIT Film T1S';
DATA_20240628.F1SL561C.N='KIT Film T1S';
DATA_20240628.F1CL561R.N='KIT Film T1S Converted';
DATA_20240628.F1CL561G.N='KIT Film T1S Converted';
DATA_20240628.F1CL561C.N='KIT Film T1S Converted';
DATA_20240628.F2SL561R.N='KIT Film T2S';
DATA_20240628.F2SL561G.N='KIT Film T2S';
DATA_20240628.F2SL561C.N='KIT Film T2S';
DATA_20240628.F2CL561R.N='KIT Film T2S Converted';
DATA_20240628.F2CL561G.N='KIT Film T2S Converted';
DATA_20240628.F2CL561C.N='KIT Film T2S Converted';
DATA_20240628.F9ML561R.N='KIT Film T9M';
DATA_20240628.F9ML561G.N='KIT Film T9M';
DATA_20240628.F9ML561C.N='KIT Film T9M';
DATA_20240628.F9CL561R.N='KIT Film T9M Converted';
DATA_20240628.F9CL561G.N='KIT Film T9M Converted';
DATA_20240628.F9CL561C.N='KIT Film T9M Converted';
DATA_20240628.LL561B.N='Laser';
DATA_20240628.LL561C.N='Laser';

%%%--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%--------SAMPLE COMPARISION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SamplesG = {
            DATA_20240610.EAL514GD
            DATA_20240610.S2L514GD
            DATA_20240610.S3L514GD
            DATA_20240610.S4L514GD
            DATA_20240610.S5L514GD
            DATA_20240610.S6L514GD
            DATA_20240610.S7L514GD
            DATA_20240610.WAL514GD    
            };

SamplesR = {
            DATA_20240610.WAH514R
            DATA_20240610.EAH514R
            DATA_20240610.S2H514R
            DATA_20240610.S3H514R
            DATA_20240610.S4H514R
            DATA_20240610.S5H514R
            DATA_20240610.S6H514R
            DATA_20240610.S7H514R
            };
    
SamplesDialR = {
            DATA_20240612.WAL514R
            DATA_20240612.EAL514R
            DATA_20240612.S2L514R
            DATA_20240612.S3L514R
            DATA_20240612.S4L514R
            DATA_20240612.S5L514R
            DATA_20240612.S6L514R
            DATA_20240612.S7L514R
            };    
        
SamplesDialG = {
            DATA_20240612.WAL514GD
            DATA_20240612.EAL514GD
            DATA_20240612.S2L514GD
            DATA_20240612.S3L514GD
            DATA_20240612.S4L514GD
            DATA_20240612.S5L514GD
            DATA_20240612.S6L514GD
            DATA_20240612.S7L514GD
            };
    
KITConvertedR = {
            DATA_20240614.T1SL514R
            DATA_20240614.T4SL514R
            DATA_20240614.T4PL514R
            };

KITConvertedG = {
            DATA_20240614.T1SL514G
            DATA_20240614.T4SL514G
            DATA_20240614.T4PL514G
            };
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
Lines561= {
            DATA_20240628.LL561A  
            DATA_20240628.LL561B  
            DATA_20240628.LL561C  
            DATA_20240628.LL561D  
       };
%Gband 561
KITF1G561= {
            DATA_20240628.F1SL561G
            DATA_20240628.F1CL561G    
       };   
KITF2G561= {
            DATA_20240628.F2SL561G
            DATA_20240628.F2CL561G  
       };    
KITF9G561= {
            DATA_20240628.F9ML561G
            DATA_20240628.F9CL561G  
       };  
   
%RBM 561
KITF1R561= {
            DATA_20240628.F1SL561R
            DATA_20240628.F1CL561R
       };   
KITF2R561= {
            DATA_20240628.F2SL561R
            DATA_20240628.F2CL561R
       };    
KITF9R561= {
            DATA_20240628.F9ML561R
            DATA_20240628.F9CL561R
       };  

%Carbyde 561
KITF1C561= {
            DATA_20240628.F1SL561C
            DATA_20240628.F1CL561C
       };   
KITF2C561= {
            DATA_20240628.F2SL561C
            DATA_20240628.F2CL561C
       };    
KITF9C561= {
            DATA_20240628.F9ML561C
            DATA_20240628.F9CL561C
       };  
   
   
%for RBMs
% SamplesR = FlatFieldCorrection(SamplesR, DATA_20240610.FFH514R)
% SamplesR = UsefulFunctions.SubstractLinearBG(SamplesR, 137, 200)
% SamplesR = UsefulFunctions.NormalizeSample(SamplesR,140, 160)

%for G/D Bands 
% SamplesG = UsefulFunctions.SubstractLinearBG(SamplesG, 1250, 1650)
% SamplesG = UsefulFunctions.NormalizeSample(SamplesG,1585, 1595)
% 

% %For DialSamples
% SamplesDialR = UsefulFunctions.SubstractLinearBG(SamplesDialR, 137, 200)
% SamplesDialR = UsefulFunctions.NormalizeSample(SamplesDialR,140, 160)
% %
% SamplesDialG = UsefulFunctions.SubstractLinearBG(SamplesDialG, 1250, 1650)
% SamplesDialG = UsefulFunctions.NormalizeSample(SamplesDialG,1585, 1595)
% 
% 
% %For KITConvertedR
% KITConvertedR = UsefulFunctions.SubstractLinearBG(KITConvertedR, 137, 200)
% KITConvertedR = UsefulFunctions.NormalizeSample(KITConvertedR,140, 160)
% % 
% KITConvertedG = UsefulFunctions.SubstractLinearBG(KITConvertedG, 1250, 1650)
% KITConvertedG = UsefulFunctions.NormalizeSample(KITConvertedG,1585, 1595)



%For KITFilms
% KITFilmsR = UsefulFunctions.SubstractLinearBG(KITFilmsR, 137, 200)
% KITFilmsR = UsefulFunctions.NormalizeSample(KITFilmsR,100, 200)
% % 
% KITFilmsG = UsefulFunctions.SubstractLinearBG(KITFilmsG, 1250, 1650)
% KITFilmsG = UsefulFunctions.NormalizeSample(KITFilmsG,1585, 1595)
% 



% 
% plotRaman(KITConvertedR, 0.7)
% plotRaman(KITConvertedG, 0.0)

% plotRaman(SamplesR, 0.7)
% plotRaman(SamplesG, 0.0)
% 
% plotRaman(SamplesDialR, 0.0)
% plotRaman(SamplesDialG, 0.0)

% plotRaman(KITFilmsR, 3.0)
% plotRaman(KITFilmsG, 3.0)



%Individual Film Analysis KIT at 514
% KITF1R = UsefulFunctions.SubstractLinearBG(KITF1R, 137, 200);
% KITF1R = UsefulFunctions.NormalizeSample(KITF1R,100, 200);
% KITF2R = UsefulFunctions.SubstractLinearBG(KITF2R, 137, 200);
% KITF2R = UsefulFunctions.NormalizeSample(KITF2R,100, 200);
% KITF9R = UsefulFunctions.SubstractLinearBG(KITF9R, 137, 200);
% KITF9R = UsefulFunctions.NormalizeSample(KITF9R,100, 200);
% 
% KITF1G = UsefulFunctions.SubstractLinearBG(KITF1G, 1250, 1650);
% KITF1G = UsefulFunctions.NormalizeSample(KITF1G,1550, 1650);
% KITF2G = UsefulFunctions.SubstractLinearBG(KITF2G, 1250, 1650);
% KITF2G = UsefulFunctions.NormalizeSample(KITF2G,1550, 1650);
% KITF9G = UsefulFunctions.SubstractLinearBG(KITF9G, 1250, 1650);
% KITF9G = UsefulFunctions.NormalizeSample(KITF9G,1550, 1650);

XAxis = [140, 500];


%Individual Film Analysis KIT at 561
% KITF1R561 = UsefulFunctions.SubstractLinearBG(KITF1R561, 100, 500);
% KITF1R561 = UsefulFunctions.NormalizeSample(KITF1R561,100, 200);
% KITF2R561 = UsefulFunctions.SubstractLinearBG(KITF2R561, 100, 500);
% KITF2R561 = UsefulFunctions.NormalizeSample(KITF2R561,100, 200);
% KITF9R561 = UsefulFunctions.SubstractLinearBG(KITF9R561, 100, 500);
% KITF9R561 = UsefulFunctions.NormalizeSample(KITF9R561,100, 500);
% % 
% KITF1G561 = UsefulFunctions.SubstractLinearBG(KITF1G561, 1280, 1624);
% KITF1G561 = UsefulFunctions.NormalizeSample(KITF1G561,1550, 1650);
% KITF2G561 = UsefulFunctions.SubstractLinearBG(KITF2G561, 1400, 1450);
% KITF2G561 = UsefulFunctions.NormalizeSample(KITF2G561,1550, 1650);
% KITF9G561 = UsefulFunctions.SubstractLinearBG(KITF9G561, 1400, 1450);
% KITF9G561 = UsefulFunctions.NormalizeSample(KITF9G561,1550, 1650);
% % 
KITF1C561 = UsefulFunctions.SubstractLinearBG(KITF1C561, 1800, 1900);
KITF1C561 = UsefulFunctions.NormalizeSample(KITF1C561,1600, 1610);
KITF2C561 = UsefulFunctions.SubstractLinearBG(KITF2C561, 1800, 1900);
KITF2C561 = UsefulFunctions.NormalizeSample(KITF2C561,1600, 1610);
KITF9C561 = UsefulFunctions.SubstractLinearBG(KITF9C561, 1800, 1900);
KITF9C561 = UsefulFunctions.NormalizeSample(KITF9C561,1600, 1610);

% plotRaman(Lines561, 0.0)


% plotRaman(KITF1R561, 0.7)
% plotRaman(KITF2R561, 0.7)
% plotRaman(KITF9R561, 0.7)
% 
% plotRaman(KITF1G561, 0.0)
% plotRaman(KITF2G561, 0.0)
% plotRaman(KITF9G561, 0.0)

plotRaman(KITF1C561, 0.0)
plotRaman(KITF2C561, 0.0)
plotRaman(KITF9C561, 0.0)