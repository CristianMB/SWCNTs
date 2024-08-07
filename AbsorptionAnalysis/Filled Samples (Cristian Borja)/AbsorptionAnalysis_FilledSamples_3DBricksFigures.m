clc;
clear;
addpath('X:\SWCNTs');
% addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\');
import UsefulFunctions.*;

% rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements Data\Absorption\';
rootpath = 'X:\Measurements Data\Absorption\';

%All paths as default
path_REF = [rootpath,'References.csv'];
path_SF = [rootpath,'20231206\alldata20231206.csv'];
path_CBS234 = [rootpath,'20240216\DGUC.csv'];
path_CBS567 = [rootpath,'20240308\DialS5S6S7.csv'];
path_Films = [rootpath,'20240614\KIT_FilmConvertedSamples.csv'];
path_DGUCSAAlkane = [rootpath,'20240724\DGU_KIT_Filled.csv'];



%Select the paths of interest
paths = {
        path_REF
        path_SF
        path_CBS234
        path_CBS567
        path_Films
        path_DGUCSAAlkane
        };




%% %--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DATA_References.empty_P2_dial_0930.N = 'Empty';
DATA_References.WaterFilled.N = 'Water';

DATA_20231206.SFF6dil.N='TTF';

DATA_20240216.S2DGUC.N='PCE';
DATA_20240216.S3DGUC.N='TCE';
DATA_20240216.S4DDGUC.N='TEMED';
DATA_20240308.S7DGUCDial.N='Dodecane';
DATA_20240308.S6DGUCDial.N='Hexadecane';
DATA_20240308.S5DGUCDial.N='TDAE';

DATA_20240614.FilmT1S.N='Film A';
DATA_20240614.FilmT1SConverted.N='Film B';
DATA_20240614.FilmT2S.N='Film C';
DATA_20240614.FilmT2SConverted.N='Film D';
DATA_20240614.FilmT9M.N='Film E';
DATA_20240614.FilmT9MConverted.N='Film F';

DATA_20240724.S6A.N='Hexadecane (After CSA)';
DATA_20240724.S7A.N='Dodecane (After CSA)';

%% %--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DATA_20240216.S2DGUC.Y = DATA_20240216.S2DGUC.Y - 0.42*3*DATA_References.WaterInD2O.Y
DATA_20240216.S3DGUC.Y = DATA_20240216.S3DGUC.Y - 0.23*3*DATA_References.WaterInD2O.Y
DATA_20240216.S4DDGUC.Y = DATA_20240216.S4DDGUC.Y- 0.079*3*DATA_References.WaterInD2O.Y
DATA_20240308.S7DGUCDial.Y = DATA_20240308.S7DGUCDial.Y - 0.01*3*DATA_References.WaterInD2O.Y
DATA_20240308.S6DGUCDial.Y = DATA_20240308.S6DGUCDial.Y- 0.01*3*DATA_References.WaterInD2O.Y
DATA_20240308.S5DGUCDial.Y = DATA_20240308.S5DGUCDial.Y- 0.01*3*DATA_References.WaterInD2O.Y
DATA_20231206.SFF6dil.Y = DATA_20231206.SFF6dil.Y - 0.087*3*DATA_References.WaterInD2O.Y

DATA_20240724.S6A.Y = DATA_20240724.S6A.Y - 3*0.525*DATA_References.Nicodenz.Y - 0.16*3*DATA_References.WaterInD2O.Y 
DATA_20240724.S7A.Y = DATA_20240724.S7A.Y - 0.5*3*DATA_References.Nicodenz.Y - 0.165*3*DATA_References.WaterInD2O.Y 
DATA_References.WaterFilled.Y = DATA_References.WaterFilled.Y - 3*0.0670*DATA_References.WaterInD2O.Y
DATA_References.empty_P2_dial_0930.Y = DATA_References.empty_P2_dial_0930.Y - 3*0.13*DATA_References.WaterInD2O.Y



%% %--------BACKGROUND CORRECTION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %--------SPECTRA SELECTION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


FilledCNTs = {
               
                DATA_References.empty_P2_dial_0930
                DATA_20231206.SFF6dil

                DATA_20240216.S2DGUC
                DATA_20240216.S3DGUC
                DATA_20240216.S4DDGUC
                DATA_20240308.S5DGUCDial
                DATA_20240308.S6DGUCDial
                DATA_20240308.S7DGUCDial

                DATA_References.WaterFilled
            };
        

% 
% FilledCNTs = NormalizeSample(FilledCNTs,800, 1100);
% 
% plotAbsorption(FilledCNTs, 0.75)
% hold on
% xline(1000, '--','LineWidth', 1.0,'HandleVisibility','off')
% xline(1764, '--','LineWidth', 1.0,'HandleVisibility','off')
% xlim([205 1880])
% % ylim([-0.5 8.0])





SFilms = {
    DATA_20240614.FilmT1S
    DATA_20240614.FilmT1SConverted
    DATA_20240614.FilmT2S
    DATA_20240614.FilmT2SConverted
    DATA_20240614.FilmT9M
    DATA_20240614.FilmT9MConverted
            };
        
% SFilms = NormalizeSample(SFilms,950, 1050);
% plotAbsorption(SFilms, 0.4)
% xlim([192 2500])
% ylim([0 4])

MFilms = {
%     DATA_20240614.FilmT1S
%     DATA_20240614.FilmT1SConverted
%     DATA_20240614.FilmT2S
%     DATA_20240614.FilmT2SConverted
%     DATA_20240614.FilmT9M
%     DATA_20240614.FilmT9MConverted
            };
        
% MFilms = NormalizeSample(MFilms,950, 1050);
% plotAbsorption(MFilms, 0.25)
% xlim([192 2500])
% ylim([0 4])

% hold on
% xline(1000, '--','LineWidth', 1.0,'HandleVisibility','off')
% xline(1764, '--','LineWidth', 1.0,'HandleVisibility','off')
% xlim([205 1880])
% ylim([-0.5 8.0])


CSACompareDode = {
                DATA_References.empty_P2_dial_0930

                DATA_20240308.S6DGUCDial
                DATA_20240724.S6A
                
                DATA_20240308.S7DGUCDial
                DATA_20240724.S7A
                
                DATA_References.WaterFilled
            };
        
% CSACompareDode = NormalizeSample(CSACompareDode,950, 1050);
% plotAbsorption(CSACompareDode, 0.8)
% xline(1774, '--','LineWidth', 1.0,'HandleVisibility','off')
% xline(1002, '--','LineWidth', 1.0,'HandleVisibility','off')
% xlim([400 2100])
% ylim([-5 1.5])
% 
% 
% CSACompareHexa = {
%                 DATA_References.empty_P2_dial_0930
%                 
%                 DATA_20240308.S6DGUCDial
%                 DATA_20240724.S6A
%                 
%                 DATA_References.WaterFilled
%             };
%         
% CSACompareHexa = NormalizeSample(CSACompareHexa,950, 1050);
% plotAbsorption(CSACompareHexa, 0.8)
% xline(1777, '--','LineWidth', 1.0,'HandleVisibility','off')
% xline(1003, '--','LineWidth', 1.0,'HandleVisibility','off')
% 
% xlim([450 2200])
% ylim([-3.5 1.5])
