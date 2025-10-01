clc;
clear;
import UsefulFunctions.*;
addpath('X:\SWCNTs\SpecialMatlabFunctions\DrosteEffect-BrewerMap-3.2.5.0')

rootpath = 'X:\Measurements Data\Absorption\';

%All paths as default

DGU_path = [rootpath,'20250819\DGU_Dispersions_TTF_TCNQ.csv'];
References = [rootpath,'References.csv'];

%Select the paths of interest
paths = {
    DGU_path
    References
        };
    

%Read and structure data from the paths
ReadAbsorptionFromPaths(paths);

%% Labelling

DATA_20250819.R23F_DGU_A_1.N='R23F - TCNQ@P2-SWCNTs (GP) (DGU A)';
DATA_20250819.R23F_DGU_B_1.N='R23F - TCNQ@P2-SWCNTs (GP) (DGU B)';
DATA_20250819.R23F_DGU_C_1.N='R23F - TCNQ@P2-SWCNTs (GP) (DGU C)';
DATA_20250819.R23F_DGU_D_1.N='R23F - TCNQ@P2-SWCNTs (GP) (DGU D)';
DATA_20250819.R23F_DGU_E_1.N='R23F - TCNQ@P2-SWCNTs (GP) (DGU E)';
DATA_20250819.R23F_DGU_F_1.N='R23F - TCNQ@P2-SWCNTs (GP) (DGU F)';
DATA_20250819.R23F_DGU_G_1.N='R23F - TCNQ@P2-SWCNTs (GP) (DGU G)';
DATA_20250819.U22_DGU_A_1.N='U22 - TTF@P2-SWCNTs (R) (DGU A)';
DATA_20250819.U22_DGU_B_1.N='U22 - TTF@P2-SWCNTs (R) (DGU B)';
DATA_20250819.U22_DGU_C_1.N='U22 - TTF@P2-SWCNTs (R) (DGU C)';
DATA_20250819.U22_DGU_D_1.N='U22 - TTF@P2-SWCNTs (R) (DGU D)';
DATA_20250819.U22_DGU_E_1.N='U22 - TTF@P2-SWCNTs (R) (DGU E)';
DATA_20250819.U22_DGU_F_1.N='U22 - TTF@P2-SWCNTs (R) (DGU F)';
DATA_20250819.U22_DGU_G_1.N='U22 - TTF@P2-SWCNTs (R) (DGU G)';
DATA_20250819.R17_DGU_A_1.N='R17 - TTF@P2-SWCNTs (GP) (DGU A)';
DATA_20250819.R17_DGU_B_1.N='R17 - TTF@P2-SWCNTs (GP) (DGU B)';
DATA_20250819.R17_DGU_C_1.N='R17 - TTF@P2-SWCNTs (GP) (DGU C)';
DATA_20250819.R17_DGU_D_1.N='R17 - TTF@P2-SWCNTs (GP) (DGU D)';
DATA_20250819.R17_DGU_E_1.N='R17 - TTF@P2-SWCNTs (GP) (DGU E)';
DATA_20250819.R17_DGU_F_1.N='R17 - TTF@P2-SWCNTs (GP) (DGU F)';
DATA_20250819.R17_DGU_G_1.N='R17 - TTF@P2-SWCNTs (GP) (DGU G)';
DATA_20250819.R21D_DGU_A_1.N='R21D - TTF@P2-SWCNTs (M) (DGU A)';
DATA_20250819.R21D_DGU_B_1.N='R21D - TTF@P2-SWCNTs (M) (DGU B)';
DATA_20250819.R21D_DGU_C_1.N='R21D - TTF@P2-SWCNTs (M) (DGU C)';
DATA_20250819.R21D_DGU_D_1.N='R21D - TTF@P2-SWCNTs (M) (DGU D)';
DATA_20250819.R21D_DGU_E_1.N='R21D - TTF@P2-SWCNTs (M) (DGU E)';
DATA_20250819.R21D_DGU_F_1.N='R21D - TTF@P2-SWCNTs (M) (DGU F)';
DATA_20250819.R21D_DGU_G_1.N='R21D - TTF@P2-SWCNTs (M) (DGU G)';

%% DGU Fraction Correction  by Water and Nycodenz

DATA_20250819.R23F_DGU_A_1.Y
DATA_20250819.R23F_DGU_B_1.Y
DATA_20250819.R23F_DGU_C_1.Y
DATA_20250819.R23F_DGU_D_1.Y
DATA_20250819.R23F_DGU_E_1.Y
DATA_20250819.R23F_DGU_F_1.Y
DATA_20250819.R23F_DGU_G_1.Y
DATA_20250819.U22_DGU_A_1.Y
DATA_20250819.U22_DGU_B_1.Y
DATA_20250819.U22_DGU_C_1.Y
DATA_20250819.U22_DGU_D_1.Y
DATA_20250819.U22_DGU_E_1.Y
DATA_20250819.U22_DGU_F_1.Y
DATA_20250819.U22_DGU_G_1.Y
DATA_20250819.R17_DGU_A_1.Y
DATA_20250819.R17_DGU_B_1.Y
DATA_20250819.R17_DGU_C_1.Y
DATA_20250819.R17_DGU_D_1.Y
DATA_20250819.R17_DGU_E_1.Y
DATA_20250819.R17_DGU_F_1.Y
DATA_20250819.R17_DGU_G_1.Y
DATA_20250819.R21D_DGU_A_1.Y
DATA_20250819.R21D_DGU_B_1.Y
DATA_20250819.R21D_DGU_C_1.Y
DATA_20250819.R21D_DGU_D_1.Y
DATA_20250819.R21D_DGU_E_1.Y
DATA_20250819.R21D_DGU_F_1.Y
DATA_20250819.R21D_DGU_G_1.Y

%% Plotting

Plot = {
%             DATA_20250819.R23F_DGU_A_1
%             DATA_20250819.R23F_DGU_B_1
%             DATA_20250819.R23F_DGU_C_1
%             DATA_20250819.R23F_DGU_D_1
%             DATA_20250819.R23F_DGU_E_1
%             DATA_20250819.R23F_DGU_F_1
%             DATA_20250819.R23F_DGU_G_1
%             DATA_20250819.U22_DGU_A_1
%             DATA_20250819.U22_DGU_B_1
%             DATA_20250819.U22_DGU_C_1
%             DATA_20250819.U22_DGU_D_1
%             DATA_20250819.U22_DGU_E_1
%             DATA_20250819.U22_DGU_F_1
%             DATA_20250819.U22_DGU_G_1
%             DATA_20250819.R17_DGU_A_1
%             DATA_20250819.R17_DGU_B_1
%             DATA_20250819.R17_DGU_C_1
%             DATA_20250819.R17_DGU_D_1
%             DATA_20250819.R17_DGU_E_1
%             DATA_20250819.R17_DGU_F_1
%             DATA_20250819.R17_DGU_G_1
%             DATA_20250819.R21D_DGU_A_1
%             DATA_20250819.R21D_DGU_B_1
%             DATA_20250819.R21D_DGU_C_1
%             DATA_20250819.R21D_DGU_D_1
%             DATA_20250819.R21D_DGU_E_1
%             DATA_20250819.R21D_DGU_F_1
%             DATA_20250819.R21D_DGU_G_1
        DATA_References.Nicodenz
        DATA_References.WaterInD2O
        };
 
plotAbsorption(Plot, 0.0)

   


function corrected_list = generate_corrections(filenames, A_values, B_values)

    corrected_list = {};
    idx = 1;
    for i = 1:length(filenames)
        for j = 1:length(A_values)
            A = A_values(j);
            B = B_values(j);
            corrected_str = sprintf('%s - %.3f*(AbsSpectWater.Y) - %.3f*(AbsSpectNycodenz.Y)', ...
                                     filenames{i}, A, B);
            corrected_list{idx} = corrected_str; 
            idx = idx + 1;
        end
    end
end