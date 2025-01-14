clc;
clear;
addpath('X:\SWCNTs');
addpath('X:\Measurements Data\Raman');
addpath('X:\SWCNTs\SpecialMatlabFunctions\DrosteEffect-BrewerMap-3.2.5.0');
import UsefulFunctions.*;
rootpath = 'X:\Measurements Data\Raman';


%All paths as default
path_20240426 = [rootpath,'\20240426\'];
path_20240514 = [rootpath,'\20240514\'];
path_20240515 = [rootpath,'\20240515\'];
path_20240614 = [rootpath,'\20240614\'];
path_20240620 = [rootpath,'\20240620\'];
path_20240628 = [rootpath,'\20240628\'];

%Select the paths of interest

paths = {
    path_20240426
    path_20240514
    path_20240515
    path_20240614
    path_20240620
    path_20240628
    };
ReadRamanFromPaths(paths, 2);


%%  DATA 20240426 FIRST TEST SAMPLES RBMS AT 650nm
       
TestSamples650_RBMs={
     DATA_20240426.BAL650R
     DATA_20240426.BBL650RA
     DATA_20240426.BAL650RB
     DATA_20240426.BBL650RB
    };

WL = 650;

for i=1:length(TestSamples650_RBMs)
    current = TestSamples650_RBMs{i};  % Access the cell array element once
    current = clip_spectrum(current, 5 ,0);
    current = remove_inclination(current, WL);
    current = correct_instrument_response(current, WL);
%     current = remove_baseline_polynomial(current, 1);
    TestSamples650_RBMs{i} = current;  % Save the result back to the cell array
end

%Forcing matching of spectra
TestSamples650_RBMs{1} = clipRangeEdges(TestSamples650_RBMs{1}, 120, 378);
TestSamples650_RBMs{3} = clipRangeEdges(TestSamples650_RBMs{3}, 378, 520);
TestSamples650_RBMs{2} = clipRangeEdges(TestSamples650_RBMs{2}, 120, 289);
TestSamples650_RBMs{4} = clipRangeEdges(TestSamples650_RBMs{4}, 289, 520);

% plotRaman(TestSamples650_RBMs, 0);

TSA650RBMs = concatenateSpectra({TestSamples650_RBMs{1} TestSamples650_RBMs{3}}, 'Test Film Sample (SWCNT)');
TSB650RBMs = concatenateSpectra({TestSamples650_RBMs{2} TestSamples650_RBMs{4}}, 'Test Film Sample (Converted to DWCNT)');
GroupedTS650RBMs = {
                    TSA650RBMs
                    TSB650RBMs
                    }; 
                

for i=1:length(GroupedTS650RBMs)
     current = GroupedTS650RBMs{i};  % Access the cell array element once
% %      current = remove_lorentzian_profileNew(current);
% %      current = remove_baseline_polynomial(current, 1);
% %      current = remove_linear_background(current, 312, 518);

    current = remove_linear_background(current, 312, 518);
    current = remove_lorentzian_profileNew(current);
    current = remove_baseline_polynomial(current, 1);
    current = remove_linear_background(current, 312, 518);

    GroupedTS650RBMs{i} = current;  % Save the result back to the cell array
end

GroupedTS650RBMs = NormalizeSample(GroupedTS650RBMs, 450, 500); 

% plotRaman(GroupedTS650RBMs, 0, WL);

% exportRamanData(GroupedTS650RBMs, 'RamanTestSamples_RBM_650nm.csv')


%%  DATA 20240426 FIRST TEST SAMPLES G and D BAND AT 650nm
% 
% DATA_20240426.BAL650D = clipRangeEdges(DATA_20240426.BAL650D, 0,1500);
% DATA_20240426.BBL650D = clipRangeEdges(DATA_20240426.BBL650D, 0,1500);
% DATA_20240426.BAL650G = clipRangeEdges(DATA_20240426.BAL650G, 1300,3000);
% DATA_20240426.BBL650G = clipRangeEdges(DATA_20240426.BBL650G, 1300,3000);


TestSamples650_DG={
     DATA_20240426.BAL650D
     DATA_20240426.BAL650G
     DATA_20240426.BBL650D
     DATA_20240426.BBL650G
    };

% clipedges = 1;
WL = 650;

for i=1:length(TestSamples650_DG)
    current = TestSamples650_DG{i};  % Access the cell array element once
    current = clip_spectrum(current, 10,10);
    current = remove_inclination(current, WL);
    current = correct_instrument_response(current, WL);
    current = remove_bg_poly(current);
    TestSamples650_DG{i} = current;  % Save the result back to the cell array
end


%Forcing matching of spectra
TestSamples650_DG{1} = clipRangeEdges(TestSamples650_DG{1}, 0, 1450);
TestSamples650_DG{2} = clipRangeEdges(TestSamples650_DG{2}, 1450, 2000);
TestSamples650_DG{3} = clipRangeEdges(TestSamples650_DG{3}, 0, 1450);
TestSamples650_DG{4} = clipRangeEdges(TestSamples650_DG{4}, 1450, 2000);


TSA650GD = concatenateSpectraNew({TestSamples650_DG{1} TestSamples650_DG{2}}, 'Test Film Sample (SWCNT)');
TSB650GD = concatenateSpectraNew({TestSamples650_DG{3} TestSamples650_DG{4}}, 'Test Film Sample (Converted to DWCNT)');
GroupedTS650GD = {
                    TSA650GD
                    TSB650GD
                    };
 
                
for i=1:length(GroupedTS650GD)
    current = GroupedTS650GD{i};  % Access the cell array element once
    current = remove_baseline_polynomial(current, 1);
    current = remove_linear_background(current, 1224, 1642);

%     current = remove_lorentzian_profileNew(current);
    GroupedTS650GD{i} = current;  % Save the result back to the cell array
end                

GroupedTS650GD = NormalizeSample(GroupedTS650GD, 0, 2000);         


% plotRaman(GroupedTS650GD, 0.0, WL);
% exportRamanData(GroupedTS650GD, 'RamanTestSamples_GDBand_650nm.csv')




%% DATA 20240514 FIRST TEST SAMPLES JUST G BAND AND CARBINE REGION AT 650nm
% 
% DATA_20240514.BBL650C1 = clipRangeEdges(DATA_20240514.BBL650C1, 1700,6000);
% DATA_20240514.BAL650C1 = clipRangeEdges(DATA_20240514.BAL650C1, 1700,6000);
% DATA_20240514.BBL650C2 = clipRangeEdges(DATA_20240514.BBL650C2, 1490,1705)
% DATA_20240514.BAL650C2 = clipRangeEdges(DATA_20240514.BAL650C2, 1490,1705)

TestSamplesB650={
    DATA_20240514.BAL650C2    
    DATA_20240514.BAL650C1
    
    DATA_20240514.BBL650C2
    DATA_20240514.BBL650C1
    };

WL = 650;

for i=1:length(TestSamplesB650)
    current = TestSamplesB650{i};  % Access the cell array element once
    current = clip_spectrum(current, 40,40);
    current = remove_inclination(current, WL);
    current = correct_instrument_response(current, WL);
    %current = remove_bg_poly(current);
    current = remove_linear_background(current, min(current.X), max(current.X));
    TestSamplesB650{i} = current;  % Save the result back to the cell array
end


TSA650b = concatenateSpectraNew({TestSamplesB650{1:2}}, 'Test Film Sample (SWCNT)');
TSB650b = concatenateSpectraNew({TestSamplesB650{3:4}}, 'Test Film Sample (Converted to DWCNT)');
GroupedTS650b = {
                    TSA650b
                    TSB650b
                    };
             
GroupedTS650b = UsefulFunctions.NormalizeSample(GroupedTS650b,0, 2000)  ;
% exportRamanData(GroupedTS650b, 'RamanTestSamples_GBandCarbyneRegion_650nm.csv')

% plotRaman(GroupedTS650b, 0.0, WL)


%% DATA 20240515 FIRST TEST SAMPLES RBMs AT 570nm

% DATA_20240515.BAL570RA = clipRangeEdges(DATA_20240515.BAL570RA, 0, 2000);
% DATA_20240515.BBL570RA = clipRangeEdges(DATA_20240515.BBL570RA, 0, 2000);
% DATA_20240515.BAL570RB = clipRangeEdges(DATA_20240515.BAL570RB, 420,780);
% DATA_20240515.BBL570RB = clipRangeEdges(DATA_20240515.BBL570RB, 420,780);


TestSamples570R={
        DATA_20240515.BAL570RA
        DATA_20240515.BAL570RB

        DATA_20240515.BBL570RA
        DATA_20240515.BBL570RB
    };

WL = 570;

for i=1:length(TestSamples570R)
    current = TestSamples570R{i};  % Access the cell array element once
    current = clip_spectrum(current, 10,10);
    current = remove_inclination(current, WL);
    current = correct_instrument_response(current, WL);
    current = subtract_constant_baseline(current, 429, 431);

    TestSamples570R{i} = current;  % Save the result back to the cell array
end

%Forcing matching of spectra
TestSamples570R{1} = clipRangeEdges(TestSamples570R{1}, 0, 435);
TestSamples570R{2} = clipRangeEdges(TestSamples570R{2}, 435, 2000);
TestSamples570R{3} = clipRangeEdges(TestSamples570R{3}, 0, 435);
TestSamples570R{4} = clipRangeEdges(TestSamples570R{4}, 435, 2000);


TSA570R = concatenateSpectraNew({TestSamples570R{1:2}}, 'Test Film Sample (SWCNT)');
TSB570R = concatenateSpectraNew({TestSamples570R{3:4}}, 'Test Film Sample (Converted to DWCNT)');

GroupedTS570R = {
                    TSA570R
                    TSB570R
                    };
                
for i=1:length(GroupedTS570R)
    current = GroupedTS570R{i};  % Access the cell array element once
    current = remove_lorentzian_profileNew(current);
    current = remove_linear_background(current, 123, 779);
    current = remove_baseline_polynomial(current, 0)

    GroupedTS570R{i} = current;  % Save the result back to the cell array
end
                
GroupedTS570R = NormalizeSample(GroupedTS570R,390, 410);                              
% plotRaman(GroupedTS570R, 0, WL)
% exportRamanData(GroupedTS570R, 'RamanTestSamples_RBM_570nm.csv')


%% DATA 20240515 FIRST TEST SAMPLES D,G,C Band AT 570nm

% 
% DATA_20240515.BAL570D = clipRangeEdges(DATA_20240515.BAL570D, 1080,1400)
% DATA_20240515.BBL570D = clipRangeEdges(DATA_20240515.BBL570D, 1080,1400)
% 
% DATA_20240515.BAL570G = clipRangeEdges(DATA_20240515.BAL570G, 1390,1700)
% DATA_20240515.BBL570G = clipRangeEdges(DATA_20240515.BBL570G, 1390,1700)
% 
% DATA_20240515.BAL570C = clipRangeEdges(DATA_20240515.BAL570C, 1690,1990)
% DATA_20240515.BBL570C = clipRangeEdges(DATA_20240515.BBL570C, 1690,1990)


TestSamples570DGC={
        DATA_20240515.BAL570D
        DATA_20240515.BAL570G
        DATA_20240515.BAL570C

        DATA_20240515.BBL570D
        DATA_20240515.BBL570G
        DATA_20240515.BBL570C
    };

WL = 570;

for i=1:length(TestSamples570DGC)
    current = TestSamples570DGC{i};  % Access the cell array element once
    current = clip_spectrum(current, 10,10);
    current = remove_inclination(current, WL);
    current = correct_instrument_response(current, WL);
    current = remove_bg_poly(current);
    TestSamples570DGC{i} = current;  % Save the result back to the cell array
end

TestSamples570DGC{1} = subtract_constant_baseline(TestSamples570DGC{1}, 1405, 1410);
TestSamples570DGC{4} = subtract_constant_baseline(TestSamples570DGC{4}, 1405, 1410);
TestSamples570DGC{2} = subtract_constant_baseline(TestSamples570DGC{2}, 1405, 1410);
TestSamples570DGC{5} = subtract_constant_baseline(TestSamples570DGC{5}, 1405, 1410);

TestSamples570DGC{1} = clipRangeEdges(TestSamples570DGC{1}, 0, 1405);
TestSamples570DGC{2} = clipRangeEdges(TestSamples570DGC{2}, 1405, 3000);
TestSamples570DGC{4} = clipRangeEdges(TestSamples570DGC{4}, 0, 1405);
TestSamples570DGC{5} = clipRangeEdges(TestSamples570DGC{5}, 1405, 3000);

TSA570DG = concatenateSpectra({TestSamples570DGC{1:2}}, 'Test Film Sample (SWCNT)');
TSB570DG = concatenateSpectra({TestSamples570DGC{4:5}}, 'Test Film Sample (Converted to DWCNT)');


TestSamples570DGC{3} = subtract_constant_baseline(TestSamples570DGC{3}, 1695, 1705);
TestSamples570DGC{6} = subtract_constant_baseline(TestSamples570DGC{6}, 1695, 1705);
TSA570DG = subtract_constant_baseline(TSA570DG, 1695, 1705);
TSB570DG = subtract_constant_baseline(TSB570DG, 1695, 1705);


TestSamples570DGC{3} = clipRangeEdges(TestSamples570DGC{3}, 1690, 3000);
TestSamples570DGC{6} = clipRangeEdges(TestSamples570DGC{6}, 1690, 3000);

TSA570DGC = concatenateSpectra({TSA570DG, TestSamples570DGC{3}}, 'Test Film Sample (SWCNT)');
TSB570DGC = concatenateSpectra({TSB570DG, TestSamples570DGC{6}}, 'Test Film Sample (Converted to DWCNT)');

GroupedTS570DGC = {
                    TSA570DGC
                    TSB570DGC
                    };
                
for i=1:length(GroupedTS570DGC)
    current = GroupedTS570DGC{i};  % Access the cell array element once
    current = remove_baseline_polynomial(current, 1);
    GroupedTS570DGC{i} = current;  % Save the result back to the cell array
end
       
                
GroupedTS570DGC = UsefulFunctions.NormalizeSample(GroupedTS570DGC,0, 2000)  ;                              
% plotRaman(GroupedTS570DGC,0, WL)
% exportRamanData(GroupedTS570DGC, 'RamanTestSamples_GDBand_570nm.csv')

%% DATA 20240614 SOLUTION SAMPLES RBMs AT 514nm

DATA_20240614.T1SL514G.N='Sol. T1 S-SWCNTs Converted';
DATA_20240614.T1SL514R.N='Sol.T1 S-SWCNTs Converted';
DATA_20240614.T4PL514G.N='Sol.T4P S-SWCNTs Converted';
DATA_20240614.T4SL514G.N='Sol.T4 S-SWCNTs Converted';
DATA_20240614.T4PL514R.N='Sol.T4P S-SWCNTs Converted';
DATA_20240614.T4SL514R.N='Sol.T4 S-SWCNTs Converted';

SolutionSamples514R={
    DATA_20240614.T1SL514R    
    DATA_20240614.T4SL514R
    DATA_20240614.T4PL514R
    };

WL = 514.5;

for i=1:length(SolutionSamples514R)
    current = SolutionSamples514R{i};  % Access the cell array element once
    current = clip_spectrum(current, 10,10);
%     current = remove_inclination(current, WL);
%     current = correct_instrument_response(current, WL);
%     current = remove_bg_poly(current);
    SolutionSamples514R{i} = current;  % Save the result back to the cell array
end

% plotRaman(SolutionSamples514R, 0.0, WL)
% exportRamanData(SolutionSamples514R, 'RamanSolutionSamples_RBM_514nm.csv')


%% DATA 20240614 SOLUTION SAMPLES GDBand AT 514nm

SolutionSamples514GD={
    DATA_20240614.T1SL514G
    DATA_20240614.T4SL514G
    DATA_20240614.T4PL514G
    };

WL = 514.5;

for i=1:length(SolutionSamples514GD)
    current = SolutionSamples514GD{i};  % Access the cell array element once
%     current = clip_spectrum(current, 10,10);
%     current = remove_inclination(current, WL);
%     current = correct_instrument_response(current, WL);
%     current = remove_bg_poly(current);
%     current = remove_baseline_polynomial(current, 0);
% 
    SolutionSamples514GD{i} = current;  % Save the result back to the cell array
end

% SolutionSamples514GD = NormalizeSample(SolutionSamples514GD, 1500, 1600);
% plotRaman(SolutionSamples514GD, 0.0, WL)
% exportRamanData(SolutionSamples514GD, 'RamanSolutionSamples_GDBand_514nm.csv')

                
%% DATA 20240620 FILM SAMPLES RBMs AT 514nm


DATA_20240620.F1SL514G.N='Film T1 S-SWCNTs';
DATA_20240620.F1CL514G.N='Film T1 S-SWCNTs Converted';
DATA_20240620.F1SL514R.N='Film T1 S-SWCNTs';
DATA_20240620.F2SL514G.N='Film T2 S-SWCNTs';
DATA_20240620.F2SL514C.N='Film T2 S-SWCNTs';
DATA_20240620.F2CL514G.N='Film T2 S-SWCNTs Converted';
DATA_20240620.F1CL514R.N='Film T1 S-SWCNTs Converted';
DATA_20240620.F9ML514G.N='Film T9 M-SWCNTs';
DATA_20240620.F9CL514G.N='Film T9 M-SWCNTs Converted';
DATA_20240620.F2SL514R.N='Film T2 S-SWCNTs';
DATA_20240620.F2CL514R.N='Film T2 S-SWCNTs Converted';
DATA_20240620.F9ML514R.N='Film T9 M-SWCNTs';
DATA_20240620.F9CL514R.N='Film T9 M-SWCNTs Converted';

FilmSamples514R={
    DATA_20240620.F1SL514R
    DATA_20240620.F1CL514R    
    DATA_20240620.F2SL514R
    DATA_20240620.F2CL514R
    DATA_20240620.F9ML514R
    DATA_20240620.F9CL514R
    };

WL = 514.5;

for i=1:length(FilmSamples514R)
    current = FilmSamples514R{i};  % Access the cell array element once
    current = clip_spectrum(current, 5,5);
    current = remove_inclination(current, WL);
    current = correct_instrument_response(current, WL);
    current = remove_lorentzian_profileNew(current);
    current = remove_linear_background(current, 240, 528);
    current = remove_baseline_polynomial(current, 0);
    FilmSamples514R{i} = current;  % Save the result back to the cell array
end
                
FilmSamples514R = UsefulFunctions.NormalizeSample(FilmSamples514R,100, 200);
% plotRamanGroup(FilmSamples514R, 2, 2, WL)
% exportRamanData(FilmSamples514R, 'RamanFilmSamples_RBM_514nm.csv')


%% DATA 20240620 FILM SAMPLES G and D Band AT 514nm


FilmSamples514GD={
    DATA_20240620.F1SL514G
    DATA_20240620.F1CL514G
    DATA_20240620.F2SL514G
    DATA_20240620.F2CL514G
    DATA_20240620.F9ML514G
    DATA_20240620.F9CL514G
    };

WL = 514.5;

for i=1:length(FilmSamples514GD)
    current = FilmSamples514GD{i};  % Access the cell array element once
    current = clip_spectrum(current, 10,10);
    current = remove_inclination(current, WL);
    current = correct_instrument_response(current, WL);
    current = remove_bg_poly(current);
    current = remove_baseline_polynomial(current, 0);
%     current = remove_linear_background(current,1243, 1679);
    FilmSamples514GD{i} = current;  % Save the result back to the cell array
end

               
FilmSamples514GD = UsefulFunctions.NormalizeSample(FilmSamples514GD,0, 2000);

% plotRamanGroup(FilmSamples514GD, 0.8, 2, WL)
% exportRamanData(FilmSamples514GD, 'RamanFilmSamples_GDBand_514nm.csv')


%% DATA 20240628 FILM SAMPLES RBMs SPEC AT 561nm

DATA_20240628.F1SL561C.N='Film T1 S-SWCNTs';
DATA_20240628.F1CL561C.N='Film T1 S-SWCNTs Converted';
DATA_20240628.F2SL561C.N='Film T2 S-SWCNTs';
DATA_20240628.F1SL561G.N='Film T1 S-SWCNTs';
DATA_20240628.F2CL561C.N='Film T2 S-SWCNTs Converted';
DATA_20240628.F1CL561G.N='Film T1 S-SWCNTs Converted';
DATA_20240628.F9ML561C.N='Film T9 M-SWCNTs';
DATA_20240628.F2SL561G.N='Film T2 S-SWCNTs';
DATA_20240628.F1SL561R.N='Film T1 S-SWCNTs';
DATA_20240628.F1CL561R.N='Film T1 S-SWCNTs Converted';
DATA_20240628.F2CL561G.N='Film T2 S-SWCNTs Converted';
DATA_20240628.F9CL561C.N='Film T9 M-SWCNTs Converted';
DATA_20240628.F2SL561R.N='Film T2 S-SWCNTs';
DATA_20240628.F9ML561G.N='Film T9 M-SWCNTs';
DATA_20240628.F2CL561R.N='Film T2 S-SWCNTs Converted';
DATA_20240628.F9ML561R.N='Film T9 M-SWCNTs';
DATA_20240628.F9CL561G.N='Film T9 M-SWCNTs Converted';
DATA_20240628.F9CL561R.N='Film T9 M-SWCNTs Converted';


FilmSamples561R={
    DATA_20240628.F1SL561R
    DATA_20240628.F1CL561R
    DATA_20240628.F2SL561R
    DATA_20240628.F2CL561R
%     DATA_20240628.F9ML561R
%     DATA_20240628.F9CL561R
    };

WL = 561;

for i=1:length(FilmSamples561R)
    current = FilmSamples561R{i};  % Access the cell array element once
    current = clip_spectrum(current, 10,10);
    current = remove_inclination(current, WL);
    current = correct_instrument_response(current, WL);
    current = remove_bg_poly(current);
    current = remove_baseline_polynomial(current, 0);
    FilmSamples561R{i} = current;  % Save the result back to the cell array
end

FilmSamples561R = UsefulFunctions.NormalizeSample(FilmSamples561R,0, 2000);

% plotRamanGroup(FilmSamples561R, 1, 2, WL)

% exportRamanData(FilmSamples561R, 'RamanSemicondFilmSamples_RBM_561nm.csv')

%% Just for metallic 
FilmSamples561RM={
    DATA_20240628.F9ML561R
    DATA_20240628.F9CL561R
    };

WL = 561;

for i=1:length(FilmSamples561RM)
    current = FilmSamples561RM{i};  % Access the cell array element once
%     current = clip_spectrum(current, 10,10);
%     current = remove_inclination(current, WL);
%     current = correct_instrument_response(current, WL);
%     current = remove_bg_poly(current);
    current = remove_linear_background(current,100,503);
    current = remove_baseline_polynomial(current, 0);

%     current = remove_baseline_polynomial(current, 0);
    FilmSamples561RM{i} = current;  % Save the result back to the cell array
end

FilmSamples561RM = UsefulFunctions.NormalizeSample(FilmSamples561RM,0, 2000);

% plotRamanPxl(FilmSamples561RM, 0, WL)
% exportRamanData(FilmSamples561RM, 'RamanMetallicFilmSamples_RBM_561nm.csv')


%% DATA 20240628 FILM SAMPLES G and D Band, Carbine SPEC AT 561nm


FilmSamples561GCD={
    DATA_20240628.F1SL561G
    DATA_20240628.F1SL561C

    DATA_20240628.F1CL561G
    DATA_20240628.F1CL561C
    
    DATA_20240628.F2SL561G
    DATA_20240628.F2SL561C
    
    DATA_20240628.F2CL561G
    DATA_20240628.F2CL561C
    
    DATA_20240628.F9ML561G
    DATA_20240628.F9ML561C
    
    DATA_20240628.F9CL561G
    DATA_20240628.F9CL561C
    };

WL = 561;

for i=1:length(FilmSamples561GCD)
    current = FilmSamples561GCD{i};  % Access the cell array element once
    current = clip_spectrum(current, 6, 6);
    current = remove_inclination(current, WL);
    current = correct_instrument_response(current, WL);
    current = remove_bg_poly(current);
    FilmSamples561GCD{i} = current;  % Save the result back to the cell array
end

FilmSamples561GCD{1} = clipRangeEdges(FilmSamples561GCD{1}, 0, 1595);
FilmSamples561GCD{3} = clipRangeEdges(FilmSamples561GCD{3}, 0, 1595);
FilmSamples561GCD{5} = clipRangeEdges(FilmSamples561GCD{5}, 0, 1595);
FilmSamples561GCD{7} = clipRangeEdges(FilmSamples561GCD{7}, 0, 1595);
FilmSamples561GCD{9} = clipRangeEdges(FilmSamples561GCD{9}, 0, 1595);
FilmSamples561GCD{11} = clipRangeEdges(FilmSamples561GCD{11}, 0, 1595);

FilmSamples561GCD{2} = clipRangeEdges(FilmSamples561GCD{2}, 1595, 3000);
FilmSamples561GCD{4} = clipRangeEdges(FilmSamples561GCD{4}, 1595, 3000);
FilmSamples561GCD{6} = clipRangeEdges(FilmSamples561GCD{6}, 1595, 3000);
FilmSamples561GCD{8} = clipRangeEdges(FilmSamples561GCD{8}, 1595, 3000);
FilmSamples561GCD{10} = clipRangeEdges(FilmSamples561GCD{10}, 1595, 3000);
FilmSamples561GCD{12} = clipRangeEdges(FilmSamples561GCD{12}, 1595, 3000);

F1S561 = concatenateSpectra({FilmSamples561GCD{1:2}}, 'Film T1S (SWCNT)');
F1C561 = concatenateSpectra({FilmSamples561GCD{3:4}}, 'Film T1S (Converted to DWCNT)');
F2S561 = concatenateSpectra({FilmSamples561GCD{5:6}}, 'Film T2S (SWCNT)');
F2C561 = concatenateSpectra({FilmSamples561GCD{7:8}}, 'Film T2S (Converted to DWCNT)');
F9M561 = concatenateSpectra({FilmSamples561GCD{9:10}}, 'Film T9M (SWCNT)');
F9C561 = concatenateSpectra({FilmSamples561GCD{11:12}}, 'Film T9M (Converted to DWCNT)');

GroupedFS561 = {
                    F1S561
                    F1C561
                    F2S561
                    F2C561
                    F9M561
                    F9C561
                    };

for i=1:length(GroupedFS561)
    current = GroupedFS561{i};  % Access the cell array element once
    current = remove_baseline_polynomial(current, 0);
    GroupedFS561{i} = current;  % Save the result back to the cell array
end
                
GroupedFS561 = UsefulFunctions.NormalizeSample(GroupedFS561,0, 2000);
% plotRamanGroup(GroupedFS561, 1, 2, WL)
% exportRamanData(GroupedFS561, 'RamanFilmSamples_GDBand_561nm.csv')

plotRamanData('RamanFilmSamples_GDBand_561nm')
plotRamanData('RamanMetallicFilmSamples_RBM_561nm')
plotRamanData('RamanSemicondFilmSamples_RBM_561nm')
plotRamanData('RamanFilmSamples_GDBand_514nm')
plotRamanData('RamanFilmSamples_RBM_514nm')
plotRamanData('RamanSolutionSamples_GDBand_514nm')
plotRamanData('RamanSolutionSamples_RBM_514nm')
plotRamanData('RamanTestSamples_GDBand_570nm')
plotRamanData('RamanTestSamples_GBandCarbyneRegion_650nm')
plotRamanData('RamanTestSamples_RBM_570nm')
plotRamanData('RamanTestSamples_GDBand_650nm')
plotRamanData('RamanTestSamples_RBM_650nm')

%% PLOTTING PREVIOUSLY EXPORTED DATA

function plotRamanData(filename)
    % filename: the CSV file generated by exportRamanData, containing repeated columns for each sample
    
    % Read data from the CSV file
    data = readcell(filename);
    
    % Extract headers and data
    headers = data(1, :);  % First row are headers
    data = data(2:end, :); % Remaining rows are data
    
    % Identify the number of samples by counting unique sets of "Name", "X", "Y"
    numSamples = sum(strcmp(headers, 'Name'));
    
    % Loop over each sample and plot X and Y
    figure;
    hold on;
    
    for i = 1:numSamples
        % Get columns for the i-th sample
        nameCol = (i - 1) * 3 + 1;
        xCol = nameCol + 1;
        yCol = nameCol + 2;
        
        % Extract sample data
        sampleNames = data(:, nameCol);
        xValues = cell2mat(data(:, xCol));
        yValues = cell2mat(data(:, yCol));
        
        % Remove NaN values for plotting
        validIdx = ~isnan(xValues) & ~isnan(yValues);
        xValues = xValues(validIdx);
        yValues = yValues(validIdx);
        
        % Get the name of the sample for legend purposes
        sampleName = unique(sampleNames(validIdx));
        if ~isempty(sampleName)
            sampleName = sampleName{1};
        else
            sampleName = ['Sample ', num2str(i)];
        end
        
        % Plot the data
        plot(xValues, yValues, 'DisplayName', sampleName);
    end
    
    % Customize the plot
    xlabel('Raman Shift');
    ylabel('Intensity');
    title('Raman Spectra');
    legend;
    hold off;
end


%% EXPORTING CORRECTED DATA

function exportRamanData(groupedData, filename)
    % groupedData: a cell array of structures with fields N (name), X, and Y.
    % filename: name of the file to save the table, e.g., 'output.csv'.
    
    % Determine the maximum length of X and Y values across all samples
    maxLength = 0;
    for i = 1:length(groupedData)
        sample = groupedData{i};
        maxLength = max(maxLength, length(sample.X));
    end
    
    % Initialize a cell array to hold all columns for the table
    tableData = {};
    
    % Loop through each sample in groupedData
    for i = 1:length(groupedData)
        % Get the current data structure
        sample = groupedData{i};
        
        % Retrieve name, X, and Y
        name = sample.N;
        xValues = sample.X(:);  % Ensure column format
        yValues = sample.Y(:);  % Ensure column format
        
        % Pad X and Y values with NaN to match maxLength
        xValuesPadded = NaN(maxLength, 1);
        yValuesPadded = NaN(maxLength, 1);
        
        xValuesPadded(1:length(xValues)) = xValues;
        yValuesPadded(1:length(yValues)) = yValues;
        
        % Create columns: name, X (padded), and Y (padded)
        nameColumn = repmat({name}, maxLength, 1); % Column with sample name
        tableData = [tableData, nameColumn, num2cell(xValuesPadded), num2cell(yValuesPadded)]; %#ok<AGROW>
    end
    
    % Convert to cell array and write to CSV file
    headers = repmat({'Name', 'X', 'Y'}, 1, length(groupedData));
    tableCell = [headers; tableData];  % Add headers as the first row
    
    % Write to CSV file without enforcing unique headers
    writecell(tableCell, filename);
    
    disp(['Table saved to ', filename]);
end



%% Functions written for better corrections

function plotRamanBreakAxis(SamplesToPlot, offset, cutfrom, cutto)
        % Create a figure for the plot
        figure;
        
        for sampleIdx = 1:length(SamplesToPlot)
            currentSample = SamplesToPlot{sampleIdx};
            
            % Get the current sample, X values, and Y values
            currentX = currentSample.X;
            currentY = currentSample.Y - offset*sampleIdx;
            currentN = currentSample.N;
            plot(currentX, currentY, 'DisplayName', currentN,'LineWidth', 1.3);
%              scatter(currentX, currentY, 1.5, 'filled', 'DisplayName', currentN);

            hold on; % Add spectra to the same plot
        end
        breakxaxis([cutfrom cutto]);
        grid on;
        % Add labels and legend
        xlabel('Raman Shift (cm^{-1})', 'FontSize', 14);
        ylabel('Normalized Intesity (a.u.)', 'FontSize', 14)
        title('Raman Spectra');
        legend('show','FontSize', 11);        % Optional: Customize the plot further if needed
        % Hold off to stop adding new plots to the current figure
        hold off;
end

function DS = clipRangeEdges(DS, min_value, max_value)
    % Find the indices where X is within the specified range
    idx = DS.X >= min_value & DS.X <= max_value;
    
    % Filter X and Y based on these indices
    DS.X = DS.X(idx);
    DS.Y = DS.Y(idx);
    DS.P = DS.P(idx);
end

function DS = remove_lorentzian_profileNew(DS)
    % Remove Lorentzian profile and constant offset from Raman spectrum
    % DS: structure with fields X (Raman shift) and Y (intensity)

    % Extract the X and Y data
    X = DS.X;  % Raman shift
    Y = DS.Y;  % Intensity values

    % Initial guesses for the parameters: [A, gamma, C]
    initial_params = [10, 1, min(Y)];  % [Amplitude, HWHM, Constant]

    % Define the Lorentzian function with a constant offset
    lorentzian_func = @(params, x) params(1) * (params(2)^2) ./ (x.^2 + params(2)^2);

    % Set bounds for parameters: [A_min, A_max; gamma_min, gamma_max; C_min, C_max]
    lb = [0, 0.01];  % Lower bounds
    ub = [Inf, Inf];  % Upper bounds

    % Define the objective function for fitting
    objective_func = @(params) Y - lorentzian_func(params, X);

    % Perform the least squares fitting using lsqcurvefit
    options = optimset('Display', 'off');  % Suppress output
    [optimal_params, ~, exitflag] = lsqcurvefit(@(params, x) lorentzian_func(params, x), ...
                                                  initial_params, ...
                                                  X, ...
                                                  Y, ...
                                                  lb, ...
                                                  ub, ...
                                                  options);

    if exitflag <= 0
        warning('Fitting did not converge. Adjust initial parameters or bounds.');
        return; % Exit the function if fitting failed
    end

    % Extract the optimized parameters
    optimal_A = optimal_params(1);
    optimal_gamma = optimal_params(2);

    % Compute the fitted Lorentzian with the constant offset using the optimized parameters
    fitted_lorentzian = lorentzian_func(optimal_params, X);

    % Calculate the corrected spectrum
    Y_corrected = Y - fitted_lorentzian;

    % Check for overcorrection and adjust parameters accordingly
    overcorrection = Y_corrected < 0;

    % Adjust the fitted parameters if necessary
    if any(overcorrection)
        % Decrease amplitude if there are negative values
        optimal_A = optimal_A * 0.8; 
        % Increase gamma to broaden the fit
        optimal_gamma = optimal_gamma * 1.1; 
        % Compute the new fitted Lorentzian with the adjusted parameters
        fitted_lorentzian = lorentzian_func([optimal_A, optimal_gamma], X);
        Y_corrected = Y - fitted_lorentzian; % Recompute the corrected spectrum
    end

    % Store the corrected Y values back into the structure
    DS.Y = Y_corrected;
    DS.X = X;

    % Optionally, display the fitted parameters
    disp(['Fitted Amplitude (A): ', num2str(optimal_A)]);
    disp(['Fitted Gamma (HWHM): ', num2str(optimal_gamma)]);
end

function DS = remove_baseline_polynomial(DS, degree)
    % Remove baseline from Raman spectrum using polynomial fitting
    % DS: structure with fields X (Raman shift) and Y (intensity)
    % degree: Degree of the polynomial used for baseline fitting

    % Extract the X and Y data
    X = DS.X;  % Raman shift (assumed centered at zero)
    Y = DS.Y;  % Intensity values

    % Identify regions to exclude based on peak detection
    % You can implement your own peak detection logic here or use findpeaks
    [pks, locs] = findpeaks(Y, 'MinPeakHeight', 0.05, 'MinPeakDistance', 10);
    
    % Create a mask for excluding the peak regions
    exclude_indices = false(size(Y));
    exclude_indices(locs) = true;

    % Fit a polynomial to the non-peak regions
    p = polyfit(X(~exclude_indices), Y(~exclude_indices), degree);  % Polynomial coefficients
    baseline = polyval(p, X);  % Evaluate the polynomial to get the baseline

    % Subtract the baseline from the original intensity
    Y_corrected = Y - baseline;

    % Find the minimum value of the corrected spectrum
    min_val = min(Y_corrected);

    % Shift the corrected spectrum so that its minimum value is zero
    Y_corrected = Y_corrected - min_val;

    % Update the structure with the corrected Y values
    DS.Y = Y_corrected;

    % Optionally, display the polynomial coefficients
    disp(['Polynomial Coefficients: ', num2str(p)]);
end


function DS = remove_lorentzian_profile(DS, A, gamma)
    % Remove Lorentzian profile centered at zero from Raman spectrum
    % DS: structure with fields X (Raman shift) and Y (intensity)
    % gamma: HWHM (half width at half maximum) for the Lorentzian profile

    % Extract the X and Y data
    X = DS.X;  % Raman shift (assumed centered at zero)
    Y = DS.Y;  % Intensity values

    % Define the Lorentzian function centered at zero
    lorentzian = @(x, A, gamma) A * (gamma^2) ./ (x.^2 + gamma^2);

    
    fitted_lorentzian = lorentzian(X, A, gamma);

    % Subtract the fitted Lorentzian from the original spectrum
    Y_corrected = Y - fitted_lorentzian;

    % Find the minimum value of the corrected spectrum
    min_val = min(Y_corrected);

    % Shift the corrected spectrum so that its minimum value is zero
    Y_corrected = Y_corrected - min_val;

    % Update the structure with the corrected Y values
    DS.Y = Y_corrected;

    % Optionally, return the fitted Lorentzian parameters
    % if nargout > 1
    %     varargout{1} = fit_params;
    % end
    
     % Optionally, display the fitted parameters
%     disp(['Fitted Amplitude (A): ', num2str(fit_params(1))]);
end

function DS = remove_linear_background(DS, x1, x2)
    % Remove linear background from Raman spectrum
    % DS: structure with fields X (Raman shift) and Y (intensity)
    % x1, x2: X-values (Raman shift) used to estimate the initial linear background

    % Extract the X and Y data
    X = DS.X;  % Raman shift
    Y = DS.Y;  % Intensity values

    % Ensure the selected X-values are within the valid range
    if x1 < min(X) || x2 > max(X) || x1 >= x2
        error('Invalid X values selected for background estimation.');
    end

    % Find the indices of the two X values in the X array
    [~, idx1] = min(abs(X - x1));  % Find the index closest to x1
    [~, idx2] = min(abs(X - x2));  % Find the index closest to x2

    % Define the Y-values corresponding to the selected X-values
    y1 = Y(idx1);
    y2 = Y(idx2);

    % Calculate the initial slope (m) and intercept (b) of the linear background
    slope = (y2 - y1) / (x2 - x1);  % m = (y2 - y1) / (x2 - x1)
    intercept = y1 - slope * x1;    % b = y1 - m * x1

    % Define the linear baseline function
    linear_baseline = @(x) slope * x + intercept;

    % Evaluate the linear baseline across the entire spectrum
    baseline = linear_baseline(X);

    % Subtract the baseline from the original intensity
    Y_corrected = Y - baseline;

    % Find the minimum value of the corrected spectrum
    min_val = min(Y_corrected);

    % If the minimum value is negative, adjust the baseline to make the spectrum non-negative
    if min_val < 0
        % Shift the entire spectrum by the absolute value of the minimum
        Y_corrected = Y_corrected - min_val;
    end

    % Update the structure with the corrected Y values
    DS.Y = Y_corrected;

    % Optionally, display the final slope and intercept values
    disp(['Final Slope: ', num2str(slope)]);
    disp(['Final Intercept: ', num2str(intercept)]);
end

function DS = subtract_constant_baseline(DS, x1, x2)
    % Subtract constant baseline from the spectrum based on a flat region
    % DS: structure with fields X (Raman shift) and Y (intensity)
    % x1: starting X value of the flat region
    % x2: ending X value of the flat region

    % Extract X and Y data
    X = DS.X;  % Raman shift
    Y = DS.Y;  % Intensity values

    % Identify the indices of the flat region between x1 and x2
    flat_region_indices = (X >= x1) & (X <= x2);

    % Calculate the average intensity value in the flat region
    average_baseline = mean(Y(flat_region_indices));

    % Subtract the average baseline from the entire spectrum
    Y_corrected = Y - average_baseline;

    % Update the structure with the corrected Y values
    DS.Y = Y_corrected;

    % Optionally, display the subtracted baseline value
    disp(['Subtracted Baseline: ', num2str(average_baseline)]);
end


function result = concatenateSpectraNew(structArray, structName)
    % Create a new structure to hold the concatenated data
    result = struct();

    % Check if the input is not empty
    if isempty(structArray)
        return; % Return empty if no input structures
    end

    % Initialize empty arrays for concatenating X and Y values
    concatenated_x = [];
    concatenated_y = [];
    
    % Loop through the input structures
    for i = 1:length(structArray)
        current_x = structArray{i}.X;  % Current spectrum X values
        current_y = structArray{i}.Y;  % Current spectrum Y values
        
        % Append current X and Y values
        for j = 1:length(current_x)
            % Check if the current X value already exists in the concatenated array
            overlapIndex = find(concatenated_x == current_x(j), 1);
            if isempty(overlapIndex)
                % No overlap, add the new X and Y values
                concatenated_x = [concatenated_x; current_x(j)];
                concatenated_y = [concatenated_y; current_y(j)];
            else
                % Overlap found: Use the Y value from the current spectrum
                concatenated_y(overlapIndex) = current_y(j);  % Update Y value with the new spectrum value
            end
        end
    end
    
    % Remove duplicate X values and keep the corresponding Y values
    [sorted_x, uniqueIndices] = unique(concatenated_x);  % Unique X values
    sorted_y = concatenated_y(uniqueIndices);  % Corresponding Y values

    % Store the sorted X and Y values in the result structure
    result.X = sorted_x;
    result.Y = sorted_y;
    result.N = structName;
end


function result = concatenateSpectra(structArray, structName)
    % Create a new structure to hold the concatenated data
    result = struct();

    % Check if the input is not empty
    if isempty(structArray)
        return; % Return empty if no input structures
    end
    % Initialize empty arrays for concatenating X and Y values
    concatenated_x = [];
    concatenated_y = [];
    
    % Loop through the input structures and concatenate the X and Y values
    for i = 1:length(structArray)
        concatenated_x = [concatenated_x; structArray{i}.X];  % Concatenate X values
        concatenated_y = [concatenated_y; structArray{i}.Y];  % Concatenate Y values
    end
    
    % Now sort the concatenated X and Y values based on the X-values
    [sorted_x, sort_idx] = sort(concatenated_x);  % Sort X values and get the sorting indices
    sorted_y = concatenated_y(sort_idx);  % Use the sorting indices to reorder the Y values
    
    % Store the sorted X and Y values in the result structure
    result.X = sorted_x;
    result.Y = sorted_y;
    result.N = structName;
end