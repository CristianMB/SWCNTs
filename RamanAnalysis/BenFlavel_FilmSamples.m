clc;
clear;
addpath('C:\Users\cborja\OneDrive - Universiteit Antwerpen\SWCNTs\');
%addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\');
import UsefulFunctions.*;

%rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements Data\Raman\';
rootpath = 'C:\Users\cborja\OneDrive - Universiteit Antwerpen\Measurements Data\Raman\';

%All paths as default
path_20240426 = [rootpath,'20240426\'];

%Select the paths of interest

paths = {
    path_20240426
    };


ReadRamanFromPaths(paths);


%%%--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       
DATA_20240426.BAL650D.N = 'S-SWCNTs Film DBand';
DATA_20240426.BAL650G.N = 'S-SWCNTs Film GBand';
DATA_20240426.BAL650R.N = 'S-SWCNTs Film RBMa';
DATA_20240426.BAL650RB.N = 'S-SWCNTs Film RBMb';

DATA_20240426.BBL650G.N = 'S-SWCNTs Converted Film GBand';
DATA_20240426.BBL650D.N = 'S-SWCNTs Converted Film DBand';
DATA_20240426.BBL650RA.N = 'S-SWCNTs Converted Film RBMa';
DATA_20240426.BBL650RB.N = 'S-SWCNTs Converted Film RBMb';

%%%--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DATA_20240426.BAL650R.Y = DATA_20240426.BAL650R.Y/2;
DATA_20240426.BAL650RB.Y = DATA_20240426.BAL650RB.Y/2;

%%%--------SAMPLE COMPARISION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

LG = 1425;
HG = 1645;
GBand650 = {DATA_20240426.BAL650G, DATA_20240426.BBL650G}
GBand650 = SubstractLinearBG(GBand650, LG, HG);
GBand650 = NormalizeSample(GBand650,1435, 1460)

LG = 1220;
HG = 1430;
DBand650 = {DATA_20240426.BAL650D, DATA_20240426.BBL650D}
DBand650 = SubstractLinearBG(DBand650, LG, HG);
DBand650 = NormalizeSample(DBand650,1435, 1460)

RBMBs650 = {DATA_20240426.BAL650RB, DATA_20240426.BBL650RB}
%RBMBs650 = NormalizeSample(RBMBs650,480, 500)
%RBMBs650 = SubstractLinearBG(RBMBs650,340, 520)

RBMAs650 = {DATA_20240426.BAL650R, DATA_20240426.BBL650RA}
%RBMAs650 = NormalizeSample(RBMAs650,330, 350)
%RBMAs650 = SubstractLinearBG(RBMAs650,340, 520)

plotRaman([GBand650,DBand650], 0)
plotRaman([RBMBs650, RBMAs650], 0)

%GDBand520 = NormalizeSample(GDBand520,LG, HG);

%GBand520 = SubstractLinearBG(GBand520, 1540, 1620);
%GBand520 = NormalizeSample(GBand520,LG, HG);

%GDBand514 = SubstractLinearBG(GDBand514,1400, 1500);
%GDBand514 = NormalizeSample(GDBand514,LG,HG);


%LR = 150;
%HR = 160;

%RBM514 = SubstractLinearBG(RBM514, 130, 210);
%RBM514 = NormalizeSample(RBM514,LR, HR);

%RBM520 = SubstractLinearBG(RBM520, 130, 210);
%RBM520 = NormalizeSample(RBM520,LR, HR);

%RBM680 = SubstractLinearBG(RBM680,147, 200);
%RBM680 = NormalizeSample(RBM680,LR, HR);





%plotRaman(RBM680,3)
%plotRaman(RBM520, 3)

%plotRaman(RBM514, 1)
%plotRaman(RBM520, 1)
%plotRaman(RBM680, 1)

%plotRaman(GDBand514, 0)
%plotRaman(SList, 0.1);
%SList = SubstractLinearBG(SList,130, 210);
%SList = NormalizeSample(SList,LR, HR);
%plotRaman(SList, 0.0);

All = {%DATA_20240426.LL650,
       %DATA_20240426.LL650B,
       %DATA_20240426.LL650C,
       
       DATA_20240426.BAL650D,
       DATA_20240426.BAL650G,
       DATA_20240426.BAL650R,
       DATA_20240426.BAL650RB,
       
       DATA_20240426.BBL650G,
       DATA_20240426.BBL650D,
       DATA_20240426.BBL650RA,
       DATA_20240426.BBL650RB
    }

plotRaman(All, 0)

%GDBand520 = GDBandPeaksCalculation(GDBand520, 1580,1600,1560,1570,1300,1400);
%exportGDBandPeaks(GDBand520, 'dsfsd.csv');



%%%%%TODO
%Use an excel file to handle the labeling

function labeledData = assignLabelsFromExcel(data, excelFile)
    % Read the Excel file
    excelData = readtable(excelFile);
    
    % Extract Raman and N columns
    ramanColumn = excelData.Raman;
    nColumn = excelData.N;
    
    % Initialize labeled data structure
    labeledData = struct();
    
    % Loop through each variable in the input data
    variableNames = fieldnames(data);
    for i = 1:numel(variableNames)
        variableName = variableNames{i};
        
        % Find the corresponding N value in the Excel data
        index = strcmp(variableName, ramanColumn);
        if any(index)
            nValue = nColumn(index);
            
            % Assign N value to the labeled data structure
            labeledData.(variableName) = data.(variableName);
            labeledData.(variableName).N = nValue;
        else
            % If the variable name is not found in the Excel data, display a warning
            warning(['Variable "', variableName, '" not found in Excel data.']);
        end
    end
end





