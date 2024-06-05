%% DATA LOAD
clc;
clear;

%Local path where general functions are stored
addpath('C:\Users\cborja\OneDrive - Universiteit Antwerpen\SWCNTs\');
% addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\');
import UsefulFunctions.*;

%Root path for data
% rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements Data\Absorption\';
rootpath = 'C:\Users\cborja\OneDrive - Universiteit Antwerpen\Measurements Data\Absorption\';

%Specific path
path_references = [rootpath,'References.csv'];

%Select the paths of interest
paths = {
        path_references
        };

%Read and structure data from selected paths.
%This create a Data Structure whose fields are the spectra.
ReadAbsorptionFromPaths(paths);


%% DATA LABELING

%DATA_References.Sample has 3 fields: X, Y and N

DATA_References.TriacontaneP2.N = 'Triacontane@SWCNT';
DATA_References.AnnoctadecaneP2.N = 'Annoctadecane@SWCNT';
DATA_References.WaterFilled.N = 'H2O@SWCNT';       
 
%% BUILD A SAMPLE LIST

SL = {    
        DATA_References.TriacontaneP2
        DATA_References.AnnoctadecaneP2
        DATA_References.WaterFilled
        };
          
%% DATA CORRECTION

%[203 604 803 1264] are points that should go to zero after 1/wavelenght subst.
SL = SubtractInverseBG(SL,[203 604 803 1264]);  

%Normalization is made to the peak found between 620 and 820 nm
SL = NormalizeSample(SL,1700, 1800);

%Plot all the spectra in SL list with an offset of 1
plotAbsorption(SL,0);

%This function tries to identify peaks with a prominence criterium
%peaksss = PeakID(SL,10);
%You can also specify the peaks
% peaks = [1000];
% SL = FitSamples(SL, peaks); 

%Plot all fittings
% for i=1:length(SL)
%     plotRamanFit(SL{i})
% end
