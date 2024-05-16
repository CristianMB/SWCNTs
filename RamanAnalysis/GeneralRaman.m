clc;
clear;
addpath('C:\Users\cborja\OneDrive - Universiteit Antwerpen\SWCNTs\');
%addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\');
import UsefulFunctions.*;

%rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements Data\Raman\';
rootpath = 'C:\Users\cborja\OneDrive - Universiteit Antwerpen\Measurements Data\Raman\';

%All paths as default
path_20240111 = [rootpath,'20240111\'];
path_20240320 = [rootpath,'20240320\'];
path_20240321 = [rootpath,'20240321\'];
path_20240325 = [rootpath,'20240325\'];
path_20240426 = [rootpath,'20240426\'];

path_20240514 = [rootpath,'20240514\'];
path_20240515 = [rootpath,'20240515\'];

%Select the paths of interest
paths = {
    path_20240111,
    path_20240320,
    path_20240321,
    path_20240325
    };

paths = {
    path_20240514,
    path_20240515
    };



ReadRamanFromPaths(paths);


%%%--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DATA_20240111.FLATHD.N='FlatField 514.5nm';
%DATA_20240111.LL514.N='LaserLine 514.5nm';
%DATA_20240111.LL514HD.N='LaserLine 514.5nm';

%DATA_20240111.S240111A.N='SF D2O@SWCNT (514.5 nm)';
%DATA_20240111.S240111C.N='SF Methanol@SWCNT (514.5 nm)';
%DATA_20240111.S240111F.N='SF PCE@SWCNT (514.5 nm)';
%DATA_20240111.S240111G.N='SF PCE@SWCNT (514.5 nm)';
%DATA_20240111.S240111H.N='SF PCE@SWCN (514.5 nm)';
%DATA_20240111.S240111B.N='SF TCE@SWCNT (514.5 nm)';
%DATA_20240111.S240111BB.N='SF TCE@SWCNT (514.5 nm)';
%DATA_20240111.S240111D.N='SF TCE@SWCNT (514.5 nm)';
%DATA_20240111.S240111E.N='SF TTF@SWCNT (514.5 nm)';
%DATA_20240111.S240111I.N='SF TEMED@SWCNT (514.5 nm)';

%DATA_20240111.S240111S.N='SC Empty@SWCNT (514.5 nm)';
%DATA_20240111.S240111J.N='SF D2O@SWCNT (514.5 nm)';
%DATA_20240111.S240111L.N='SF Methanol@SWCNT (514.5 nm)';
%DATA_20240111.S240111O.N='SF PCE@SWCNT (514.5 nm)';
%DATA_20240111.S240111P.N='SF PCE@SWCNT (514.5 nm)';
%DATA_20240111.S240111Q.N='SF PCE@SWCN (514.5 nm)';
%DATA_20240111.S240111K.N='SF TCE@SWCNT (514.5 nm)';
%DATA_20240111.S240111KK.N='SF TCE@SWCNT (514.5 nm)';
%DATA_20240111.S240111M.N='SF TCE@SWCNT (514.5 nm)';
%DATA_20240111.S240111N.N='SF TTF@SWCNT (514.5 nm)';
%DATA_20240111.S240111R.N='SF TEMED@SWCNT (514.5 nm)';

%DATA_20240320.FFH520R.N='FlatField (520.25 nm)';
%DATA_20240320.L520A.N='LaserLine (520.25 nm)';

%DATA_20240320.S2H520R.N='CB PCE@SWCNT Dial. DGU C (Filled) (520.25 nm)';
%DATA_20240320.S3H520R.N='CB TCE@SWCNT Dial. DGU C (Filled) (520.25 nm)';
%DATA_20240320.S4H520R.N='CB TEMED@SWCNT Dial. DGU C (Filled) (520.25 nm)';
%DATA_20240320.S5H520R.N='CB TDAE@SWCNT Dial. DGU C (Filled) (520.25 nm)';
%DATA_20240320.S6H520R.N='CB Hexadecane@SWCNT Dial. DGU C (Filled) (520.25 nm)';
%DATA_20240320.S7H520R.N='CB Dodecane@SWCNT Dial. DGU C (Filled) (520.25 nm)';

%DATA_20240321.FFH520G.N='FlatField (520.25 nm)';
%DATA_20240321.FFL520G.N='FlatField (520.25 nm)';
%DATA_20240321.L520HG.N='LaserLine HD Before (520.25 nm)';
%DATA_20240321.L520L.N='LaserLine LD (520.25 nm)';
%DATA_20240321.LLH520GB.N='LaserLine HD After (520.25 nm)';

%DATA_20240321.S2H520G.N='CB PCE@SWCNT Dial. DGU C (Filled) (520.25 nm)';
%DATA_20240321.S2L520G.N='CB PCE@SWCNT Dial. DGU C (Filled) (520.25 nm)';
%DATA_20240321.S3H520G.N='CB TCE@SWCNT Dial. DGU C (Filled) (520.25 nm)';
%DATA_20240321.S3L520G.N='CB TCE@SWCNT Dial. DGU C (Filled) (520.25 nm)';
%DATA_20240321.S4H520G.N='CB TEMED@SWCNT Dial. DGU C (Filled) (520.25 nm)';
%DATA_20240321.S4L520G.N='CB TEMED@SWCNT Dial. DGU C (Filled) (520.25 nm)';
%DATA_20240321.S5H520G.N='CB TDAE@SWCNT Dial. DGU C (Filled) (520.25 nm)';
%DATA_20240321.S5L520G.N='CB TDAE@SWCNT Dial. DGU C (Filled) (520.25 nm)';
%DATA_20240321.S6H520G.N='CB Hexadecane@SWCNT Dial. DGU C (Filled) (520.25 nm)';
%DATA_20240321.S6L520G.N='CB Hexadecane@SWCNT Dial. DGU C (Filled) (520.25 nm)';
%DATA_20240321.S7H520G.N='CB Dodecane@SWCNT Dial. DGU C (Filled) (520.25 nm)';
%DATA_20240321.S7L520G.N='CB Dodecane@SWCNT Dial. DGU C (Filled) (520.25 nm)';

%DATA_20240325.REH680R.N='SC Empty@SWCNT (680.02 nm)';
%DATA_20240325.RWH680R.N='SC Water@SWCNT (680.02 nm)';
%DATA_20240325.S2H680R.N='CB PCE@SWCNT Dial. DGU C (Filled) (680.02 nm)';
%DATA_20240325.S3H680R.N='CB TCE@SWCNT Dial. DGU C (Filled) (680.02 nm)';
%DATA_20240325.S4H680R.N='CB TEMED@SWCNT Dial. DGU C (Filled) (680.02 nm)';
%DATA_20240325.S5H680R.N='CB TDAE@SWCNT Dial. DGU C (Filled) (680.02 nm)';
%DATA_20240325.S6H680R.N='CB Hexadecane@SWCNT Dial. DGU C (Filled) (680.02 nm)';
%DATA_20240325.S7H680R.N='CB Dodecane@SWCNT Dial. DGU C (Filled) (680.02 nm)';




%%%--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FlatFields Normalization


%FF Correction to RBMs region
%SList = {
%        DATA_20240111.S240111J,
 %       DATA_20240111.S240111K,
 %       DATA_20240111.S240111KK,
 %       DATA_20240111.S240111L,
 %       DATA_20240111.S240111M,
 %       DATA_20240111.S240111N,
 %       DATA_20240111.S240111O,
 %       DATA_20240111.S240111P, 
 %       DATA_20240111.S240111Q,
 %       DATA_20240111.S240111R,
 %       DATA_20240111.S240111S
 %       };

%SList = FlatFieldCorrection(SList, DATA_20240111.FLATHD)


%FF Correction to RBMs region
%DATA_20240320.S2H520R.Y = DATA_20240320.S2H520R.Y ./ DATA_20240320.FFH520R.Y;
%DATA_20240320.S3H520R.Y = DATA_20240320.S3H520R.Y ./ DATA_20240320.FFH520R.Y;
%DATA_20240320.S4H520R.Y = DATA_20240320.S4H520R.Y ./ DATA_20240320.FFH520R.Y;
%DATA_20240320.S5H520R.Y = DATA_20240320.S5H520R.Y ./ DATA_20240320.FFH520R.Y;
%DATA_20240320.S6H520R.Y = DATA_20240320.S6H520R.Y ./ DATA_20240320.FFH520R.Y;
%DATA_20240320.S7H520R.Y = DATA_20240320.S7H520R.Y ./ DATA_20240320.FFH520R.Y;

%FF Correction to GBand region in LD Mode
%DATA_20240321.S2L520G.Y = DATA_20240321.S2L520G.Y ./ DATA_20240321.FFL520G.Y;
%DATA_20240321.S3L520G.Y = DATA_20240321.S3L520G.Y ./ DATA_20240321.FFL520G.Y;
%DATA_20240321.S4L520G.Y = DATA_20240321.S4L520G.Y ./ DATA_20240321.FFL520G.Y;
%DATA_20240321.S5L520G.Y = DATA_20240321.S5L520G.Y ./ DATA_20240321.FFL520G.Y;
%DATA_20240321.S6L520G.Y = DATA_20240321.S6L520G.Y ./ DATA_20240321.FFL520G.Y;
%DATA_20240321.S7L520G.Y = DATA_20240321.S7L520G.Y ./ DATA_20240321.FFL520G.Y;

%FF Correction to GBand region in HD Mode
%DATA_20240321.S2H520G.Y = DATA_20240321.S2H520G.Y ./ DATA_20240321.FFH520G.Y;
%DATA_20240321.S3H520G.Y = DATA_20240321.S3H520G.Y ./ DATA_20240321.FFH520G.Y;
%DATA_20240321.S4H520G.Y = DATA_20240321.S4H520G.Y ./ DATA_20240321.FFH520G.Y;
%DATA_20240321.S5H520G.Y = DATA_20240321.S5H520G.Y ./ DATA_20240321.FFH520G.Y;
%DATA_20240321.S6H520G.Y = DATA_20240321.S6H520G.Y ./ DATA_20240321.FFH520G.Y;
%DATA_20240321.S7H520G.Y = DATA_20240321.S7H520G.Y ./ DATA_20240321.FFH520G.Y;

%FF Correction to GBand region in HD Mode
%%DATA_20240325.REH680R.Y = DATA_20240325.REH680R.Y ./ DATA_20240325.FH680R.Y;
%DATA_20240325.RWH680R.Y = DATA_20240325.RWH680R.Y ./ DATA_20240325.FH680R.Y;
%DATA_20240325.S2H680R.Y = DATA_20240325.S2H680R.Y ./ DATA_20240325.FH680R.Y;
%DATA_20240325.S3H680R.Y = DATA_20240325.S3H680R.Y ./ DATA_20240325.FH680R.Y;
%DATA_20240325.S4H680R.Y = DATA_20240325.S4H680R.Y ./ DATA_20240325.FH680R.Y;
%DATA_20240325.S5H680R.Y = DATA_20240325.S5H680R.Y ./ DATA_20240325.FH680R.Y;
%DATA_20240325.S6H680R.Y = DATA_20240325.S6H680R.Y ./ DATA_20240325.FH680R.Y;
%DATA_20240325.S7H680R.Y = DATA_20240325.S7H680R.Y ./ DATA_20240325.FH680R.Y;



%%%--------SAMPLE COMPARISION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%GDBand514 = {
%        DATA_20240111.S240111A,
%        DATA_20240111.S240111B,
%        DATA_20240111.S240111BB,
%        DATA_20240111.S240111C,
%        DATA_20240111.S240111D,
%        DATA_20240111.S240111E,
%        DATA_20240111.S240111F,
%        DATA_20240111.S240111G,
%        DATA_20240111.S240111H,
%        DATA_20240111.S240111I
%    };
    
%GDBand520 = {
%        DATA_20240321.S2L520G,
%        DATA_20240321.S3L520G,
%        DATA_20240321.S4L520G,
%        DATA_20240321.S5L520G,
%        DATA_20240321.S6L520G,
%        DATA_20240321.S7L520G,
%        };
    
%GBand520 = {
%        DATA_20240321.S2H520G,
%        DATA_20240321.S3H520G,
%        DATA_20240321.S4H520G,
%        DATA_20240321.S5H520G,
%        DATA_20240321.S6H520G,
%        DATA_20240321.S7H520G,
%        };

    
%RBM514 = {
        %Empty SC
%        DATA_20240111.S240111S
        %D2O Filled
%        DATA_20240111.S240111J,    
        %PCE
 %       DATA_20240111.S240111O,
%        DATA_20240111.S240111P,
 %       DATA_20240111.S240111Q,
        %TCE
 %       DATA_20240111.S240111K,
 %       DATA_20240111.S240111M,
        %TEMED
  %      DATA_20240111.S240111R,
        %TTF
   %     DATA_20240111.S240111N,       
  %      };

%RBM520 = {  
        %Empty
 %       DATA_20240111.S240111S
  %      %D2O Filled
   %     DATA_20240111.S240111J,

    %    DATA_20240320.S2H520R,
     %   DATA_20240320.S3H520R,
      %  DATA_20240320.S4H520R,
       % DATA_20240320.S5H520R,
        %DATA_20240320.S6H520R,
        %DATA_20240320.S7H520R,
    %};
    
%RBM680 = {
%        DATA_20240325.REH680R,
%        DATA_20240325.RWH680R,

%        DATA_20240325.S2H680R,
 %       DATA_20240325.S3H680R,
 %       DATA_20240325.S4H680R,
%        DATA_20240325.S5H680R,
%        DATA_20240325.S6H680R,
%        DATA_20240325.S7H680R
 %       };
    
    
%LG = 1560;
%HG = 1570;

%GDBand520 = SubstractLinearBG(GDBand520, 1400, 1500);
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
% 
% All = {DATA_20240426.LL650,
%         DATA_20240426.LL650B,
%         DATA_20240426.BAL650D,
%         DATA_20240426.BAL650G,
%         DATA_20240426.BAL650R,
%         DATA_20240426.BAL650RB
%     }

BENSAM = {DATA_20240515.BAL570C
          DATA_20240515.BAL570D
          DATA_20240515.BAL570G
          DATA_20240515.BAL570RA
          DATA_20240515.BAL570RB
          DATA_20240515.BBL570C
          DATA_20240515.BBL570D
          DATA_20240515.BBL570G
          DATA_20240515.BBL570RA
          DATA_20240515.BBL570RB
          }

plotRaman(BENSAM, 0)

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





