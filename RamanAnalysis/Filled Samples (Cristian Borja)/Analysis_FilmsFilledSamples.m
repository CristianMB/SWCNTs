clc;
clear;
addpath('X:\SWCNTs');
% addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\');
import UsefulFunctions.*;

% rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements Data\Absorption\';
rootpath = 'X:\Measurements Data\Raman\';

%All paths as default
path_FS514 = [rootpath,'20241007\'];
path_F0514 = [rootpath,'20241008\'];

%Select the paths of interest

paths = {
    path_FS514
    path_F0514
    };


ReadRamanFromPaths(paths, 2);


%%%--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DATA_20241007.LL514A.N='Laser';
DATA_20241007.F2L514R.N='Film PCE@SWCNT';
DATA_20241007.F3L514R.N='Film TCE@SWCNT';
DATA_20241007.F4L514R.N='Film TEMED@SWCNT';
DATA_20241007.F5L514R.N='Film TDEA@SWCNT';
DATA_20241007.F6L514R.N='Film Hexadecane@SWCNT';
DATA_20241007.F7L514R.N='Film Dodecane@SWCNT';
DATA_20241007.LL514B.N='Laser';
DATA_20241007.F2L514GD.N='Film PCE@SWCNT';
DATA_20241007.F3L514GD.N='Film TCE@SWCNT';
DATA_20241007.F4L514GD.N='Film TEMED@SWCNT';
DATA_20241007.F5L514GD.N='Film TDEA@SWCNT';
DATA_20241007.F6L514GD.N='Film Hexadecane@SWCNT';
DATA_20241007.F7L514GD.N='Film Dodecane@SWCNT';
DATA_20241007.LL514C.N='Laser';
DATA_20241007.F2L514DD.N='Film PCE@SWCNT';
DATA_20241007.F3L514DD.N='Film TCE@SWCNT';
DATA_20241007.F4L514DD.N='Film TEMED@SWCNT';
DATA_20241007.F5L514DD.N='Film TDEA@SWCNT';
DATA_20241007.F6L514DD.N='Film Hexadecane@SWCNT';
DATA_20241007.F7L514DD.N='Film Dodecane@SWCNT';
DATA_20241007.F0L514R.N='Sapphire Substrate';


DATA_20241008.F0L514A.N='Sapphire SubstrateA';
DATA_20241008.FFL514A.N='FlatFieldA';
DATA_20241008.F0L514B.N='Sapphire SubstrateB';
DATA_20241008.FFL514B.N='FlatFieldB';
DATA_20241008.F0L514C.N='Sapphire SubstrateC';
DATA_20241008.FFL514C.N='FlatFieldC';
DATA_20241008.F0L514D.N='Sapphire SubstrateD';
DATA_20241008.FFL514D.N='FlatFieldD';
DATA_20241008.F0L514E.N='Sapphire SubstrateE';
DATA_20241008.FFL514E.N='FlatFieldE';
DATA_20241008.F0L514F.N='Sapphire SubstrateF';
DATA_20241008.FFL514F.N='FlatFieldF';
DATA_20241008.F0L514G.N='Sapphire SubstrateG';
DATA_20241008.FFL514G.N='FlatFieldG';

%%%--------MANUAL CORRECTIONS--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%--------SAMPLE COMPARISION--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%RBMs
FS514R = {
    DATA_20241007.F2L514R
    DATA_20241007.F3L514R
    DATA_20241007.F4L514R
    DATA_20241007.F5L514R
    DATA_20241007.F6L514R
    DATA_20241007.F7L514R
%     DATA_20241007.FFL514R
    };
% 
% plotRaman(FS514R,0)
FS514R = FlatFieldCorrection(FS514R, DATA_20241007.FFL514R);
% plotRaman(FS514R,0)

%%%%%GBand
FS514GD = {
    DATA_20241007.F2L514GD
    DATA_20241007.F3L514GD
    DATA_20241007.F4L514GD
    DATA_20241007.F5L514GD
    DATA_20241007.F6L514GD
    DATA_20241007.F7L514GD
%     DATA_20241007.FFL514GD
    };

% plotRaman(FS514GD,0)
FS514GD = FlatFieldCorrection(FS514GD, DATA_20241007.FFL514GD);
% plotRaman(FS514GD,0)


%%%%%2DBand
FS514DD = {
    DATA_20241007.F2L514DD
    DATA_20241007.F3L514DD
    DATA_20241007.F4L514DD
    DATA_20241007.F5L514DD
    DATA_20241007.F6L514DD
    DATA_20241007.F7L514DD
%     DATA_20241007.FFL514DD
    };

% plotRaman(FS514DD,0)
FS514DD = FlatFieldCorrection(FS514DD, DATA_20241007.FFL514DD);
% plotRaman(FS514DD,0)



%%%--------SUBSTRATE PLOT--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DATA_20241008.F0L514A = FlatFieldCorrectionSingle(DATA_20241008.F0L514A, DATA_20241008.FFL514A);
DATA_20241008.F0L514B = FlatFieldCorrectionSingle(DATA_20241008.F0L514B, DATA_20241008.FFL514B);
DATA_20241008.F0L514C = FlatFieldCorrectionSingle(DATA_20241008.F0L514C, DATA_20241008.FFL514C);
DATA_20241008.F0L514D = FlatFieldCorrectionSingle(DATA_20241008.F0L514D, DATA_20241008.FFL514D);
DATA_20241008.F0L514E = FlatFieldCorrectionSingle(DATA_20241008.F0L514E, DATA_20241008.FFL514E);
DATA_20241008.F0L514F = FlatFieldCorrectionSingle(DATA_20241008.F0L514F, DATA_20241008.FFL514F);
DATA_20241008.F0L514G = FlatFieldCorrectionSingle(DATA_20241008.F0L514G, DATA_20241008.FFL514G);

F0514 = {
    DATA_20241008.F0L514A
    DATA_20241008.F0L514B
    DATA_20241008.F0L514C
    DATA_20241008.F0L514D
    DATA_20241008.F0L514E
    DATA_20241008.F0L514F
    DATA_20241008.F0L514G
    };

WL = 514.5;
for i=1:length(F0514)
    current = F0514{i};  % Access the cell array element once
    current = clip_spectrum(current, 30,30);
    current = remove_inclination(current, WL);
    current = correct_instrument_response(current, WL);
    current = remove_bg_poly(current);
    F0514{i} = current;  % Save the result back to the cell array
end



plotRaman(F0514, 0.0);




