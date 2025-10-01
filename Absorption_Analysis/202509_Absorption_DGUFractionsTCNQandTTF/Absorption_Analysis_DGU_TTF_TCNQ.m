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

samplesToCorrect = {
                    DATA_20250819.R23F_DGU_A_1
                    DATA_20250819.R23F_DGU_B_1
                    DATA_20250819.R23F_DGU_C_1
                    DATA_20250819.R23F_DGU_D_1
                    DATA_20250819.R23F_DGU_E_1
                    DATA_20250819.R23F_DGU_F_1
                    DATA_20250819.R23F_DGU_G_1
                    DATA_20250819.U22_DGU_A_1
                    DATA_20250819.U22_DGU_B_1
                    DATA_20250819.U22_DGU_C_1
                    DATA_20250819.U22_DGU_D_1
                    DATA_20250819.U22_DGU_E_1
                    DATA_20250819.U22_DGU_F_1
                    DATA_20250819.U22_DGU_G_1
                    DATA_20250819.R17_DGU_A_1
                    DATA_20250819.R17_DGU_B_1
                    DATA_20250819.R17_DGU_C_1
                    DATA_20250819.R17_DGU_D_1
                    DATA_20250819.R17_DGU_E_1
                    DATA_20250819.R17_DGU_F_1
                    DATA_20250819.R17_DGU_G_1
                    DATA_20250819.R21D_DGU_A_1
                    DATA_20250819.R21D_DGU_B_1
                    DATA_20250819.R21D_DGU_C_1
                    DATA_20250819.R21D_DGU_D_1
                    DATA_20250819.R21D_DGU_E_1
                    DATA_20250819.R21D_DGU_F_1
                    DATA_20250819.R21D_DGU_G_1
};

A_values = {
    0.16
    0.27
    0.38
    0.435
    0.495
    0.58
    0.88
    0.16
    0.27
    0.41
    0.435
    0.48
    0.54
    0.81
    0.1
    0.15
    0.36
    0.42
    0.47
    0.54
    0.77
    0.16
    0.28
    0.34
    0.43
    0.5
    0.64
    0.97
};

B_values = {
    0.49
    0.32
    0.22
    0.18
    0.137
    0.066
    -0.15
    0.49
    0.32
    0.2
    0.185
    0.14
    0.11
    -0.09
    0.41
    0.5
    0.25
    0.2
    0.135
    0.12
    -0.07
    0.51
    0.31
    0.25
    0.22
    0.125
    0.06
    -0.16
};

correctedSamples = CorrectSamples(samplesToCorrect, A_values, DATA_References.Nicodenz, B_values, DATA_References.WaterInD2O)

DATA_20250819.R23F_DGU_A_1 = correctedSamples{1};
DATA_20250819.R23F_DGU_B_1 = correctedSamples{2};
DATA_20250819.R23F_DGU_C_1 = correctedSamples{3};
DATA_20250819.R23F_DGU_D_1 = correctedSamples{4};
DATA_20250819.R23F_DGU_E_1 = correctedSamples{5};
DATA_20250819.R23F_DGU_F_1 = correctedSamples{6};
DATA_20250819.R23F_DGU_G_1 = correctedSamples{7};
DATA_20250819.U22_DGU_A_1 = correctedSamples{8};
DATA_20250819.U22_DGU_B_1 = correctedSamples{9};
DATA_20250819.U22_DGU_C_1 = correctedSamples{10};
DATA_20250819.U22_DGU_D_1 = correctedSamples{11};
DATA_20250819.U22_DGU_E_1 = correctedSamples{12};
DATA_20250819.U22_DGU_F_1 = correctedSamples{13};
DATA_20250819.U22_DGU_G_1 = correctedSamples{14};
DATA_20250819.R17_DGU_A_1 = correctedSamples{15};
DATA_20250819.R17_DGU_B_1 = correctedSamples{16};
DATA_20250819.R17_DGU_C_1 = correctedSamples{17};
DATA_20250819.R17_DGU_D_1 = correctedSamples{18};
DATA_20250819.R17_DGU_E_1 = correctedSamples{19};
DATA_20250819.R17_DGU_F_1 = correctedSamples{20};
DATA_20250819.R17_DGU_G_1 = correctedSamples{21};
DATA_20250819.R21D_DGU_A_1 = correctedSamples{22};
DATA_20250819.R21D_DGU_B_1 = correctedSamples{23};
DATA_20250819.R21D_DGU_C_1 = correctedSamples{24};
DATA_20250819.R21D_DGU_D_1 = correctedSamples{25};
DATA_20250819.R21D_DGU_E_1 = correctedSamples{26};
DATA_20250819.R21D_DGU_F_1 = correctedSamples{27};
DATA_20250819.R21D_DGU_G_1 = correctedSamples{28};

%% Plotting

 
DGUFractions = {
        DATA_20250819.R23F_DGU_A_1
        DATA_20250819.R23F_DGU_B_1
        DATA_20250819.R23F_DGU_C_1
        DATA_20250819.R23F_DGU_D_1
        DATA_20250819.R23F_DGU_E_1
        DATA_20250819.R23F_DGU_F_1
        DATA_20250819.R23F_DGU_G_1
        
%         DATA_20250819.U22_DGU_A_1
%         DATA_20250819.U22_DGU_B_1
%         DATA_20250819.U22_DGU_C_1
%         DATA_20250819.U22_DGU_D_1
%         DATA_20250819.U22_DGU_E_1
%         DATA_20250819.U22_DGU_F_1
%         DATA_20250819.U22_DGU_G_1
%         
%         DATA_20250819.R17_DGU_A_1
%         DATA_20250819.R17_DGU_B_1
%         DATA_20250819.R17_DGU_C_1
%         DATA_20250819.R17_DGU_D_1
%         DATA_20250819.R17_DGU_E_1
%         DATA_20250819.R17_DGU_F_1
%         DATA_20250819.R17_DGU_G_1
% %         
%         DATA_20250819.R21D_DGU_A_1
%         DATA_20250819.R21D_DGU_B_1
%         DATA_20250819.R21D_DGU_C_1
%         DATA_20250819.R21D_DGU_D_1
%         DATA_20250819.R21D_DGU_E_1
%         DATA_20250819.R21D_DGU_F_1
%         DATA_20250819.R21D_DGU_G_1
};

DGUFractions = Normalize(DGUFractions, 920, 1120, 'M');
plotAbsorption(DGUFractions, 0.0)


   


function correctedSamples = CorrectSamples(samplesToCorrect, A_values, A_spec, B_values, B_spec)
    % Correct the samples by subtracting A*A_spec.Y + B*B_spec.Y
    % Works even if A_spec / B_spec have different lengths by interpolating.
    %
    % samplesToCorrect : cell array of structs with fields .X and .Y
    % A_values, B_values : arrays (same length as samplesToCorrect)
    % A_spec, B_spec : structs with fields .X and .Y (reference spectra)

    correctedSamples = cell(size(samplesToCorrect));

    % Iterate over each sample
    for sampleIdx = 1:length(samplesToCorrect)
        currentSample = samplesToCorrect{sampleIdx};

        A = A_values{sampleIdx};
        B = B_values{sampleIdx};

        % Interpolate reference spectra onto current sample's X-grid
        A_interp = interp1(A_spec.X, A_spec.Y, currentSample.X, 'linear', 'extrap');
        B_interp = interp1(B_spec.X, B_spec.Y, currentSample.X, 'linear', 'extrap');

        % Apply correction
        correctedSample = currentSample;
        % I have put by hand a factor of 3 to consider the 3mm cell
        correctedSample.Y = currentSample.Y - 3*A*A_interp - 3*B*B_interp;

        correctedSamples{sampleIdx} = correctedSample;
    end
end
