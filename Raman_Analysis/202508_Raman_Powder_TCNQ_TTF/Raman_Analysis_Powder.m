clc;
clear;
import UsefulFunctions.*;
rootpath = 'X:\Measurements Data\Raman\';
addpath('X:\SWCNTs\RamanAnalysis\Raman - Voigt Fitting\Faddeeva_voigt');
addpath('X:\SWCNTs');

% rootpath = 'X:\Measurements Data\Raman\';
%All paths as default

path_TTF = [rootpath,'20250520\'];
path_TCNQ = [rootpath,'20250618\'];
path_a = [rootpath,'20250131\'];
path_b = [rootpath,'20250403\'];

path_TCNQ_rinsed = [rootpath,'20250718\'];
path_S21_S22_S23_G = [rootpath,'20250731\'];
path_intensity = [rootpath,'20250801\'];


%Select the paths of interest

paths = {
    path_TTF
    path_TCNQ
    path_a
    path_b
    path_TCNQ_rinsed
    path_S21_S22_S23_G
    path_intensity
        };


ReadRamanFromPaths(paths, 2);

%% %--------LABELING--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


DATA_20250131.FFL514D.N='FlatField';
DATA_20250131.FFL514G.N='FlatField';
DATA_20250131.FFL514R.N='FlatField';
DATA_20250131.KITL514G .N='KIT Converted to DWCNTs (Dispersion)';
DATA_20250131.KITL514R  .N='KIT Converted to DWCNTs (Dispersion)';
DATA_20250131.KT3L514G .N='KIT Converted to DWCNTs (Dispersion)';
DATA_20250131.KT5L514G.N='KIT Converted to DWCNTs (Dispersion)';
DATA_20250131.KTLL514G.N='KIT Converted to DWCNTs (Dispersion)';
DATA_20250131.KTLL514R.N='KIT Converted to DWCNTs (Dispersion)';
DATA_20250131.LL514A.N='Laser';
DATA_20250131.LL514B.N='Laser';
DATA_20250131.LL514C .N='Laser';
DATA_20250131.LL514D.N='Laser';
DATA_20250131.P12DD.N='Powder S12 - TTF@P2-SWCNTs (Rinsed)';
DATA_20250131.P12G.N='Powder S12 - TTF@P2-SWCNTs (Rinsed)';
DATA_20250131.P12L514D.N='Powder S12 - TTF@P2-SWCNTs (Rinsed)';
DATA_20250131.P12L514G.N='Powder S12 - TTF@P2-SWCNTs (Rinsed)';
DATA_20250131.P12L514R.N='Powder S12 - TTF@P2-SWCNTs (Rinsed)';
DATA_20250131.P12RBM.N='Powder S12 - TTF@P2-SWCNTs (Rinsed)';
DATA_20250131.P2ADD.N='Powder P2-SWCNTs Annealed';
DATA_20250131.P2AG.N='Powder P2-SWCNTs Annealed';
DATA_20250131.P2ARBM.N='Powder P2-SWCNTs Annealed';
DATA_20250131.PAPL514D.N='Powder AP-SWCNTs';
DATA_20250131.PAPL514G.N='Powder AP-SWCNTs';
DATA_20250131.PAPL514R.N='Powder AP-SWCNTs';
DATA_20250131.PP2L514D.N='Powder P2-SWCNTs Annealed';
DATA_20250131.PP2L514G .N='Powder P2-SWCNTs Annealed';
DATA_20250131.PP2L514R.N='Powder P2-SWCNTs Annealed';
DATA_20250131.S10L514G.N='Dispersion S10 - TMG@P2-SWCNTs';
DATA_20250131.S10L514R.N='Dispersion S10 - TMG@P2-SWCNTs';
DATA_20250131.S11L514G.N='Dispersion S11 - TTF@P2-SWCNTs';
DATA_20250131.S11L514R.N='Dispersion S11 - TTF@P2-SWCNTs';
DATA_20250131.SR0L514G.N='Dispersion SR0 - Empty@P2-SWCNTs';
DATA_20250131.SR0L514R.N='Dispersion SR0 - Empty@P2-SWCNTs';
DATA_20250131.SR1L514G.N='Dispersion SR1 - D2O@P2-SWCNTs - KVD';
DATA_20250131.SR1L514R.N='Dispersion SR1 - D2O@P2-SWCNTs - KVD';
DATA_20250131.SR2L514G.N='Dispersion SR2 - MeOH@P2-SWCNTs';
DATA_20250131.SR2L514R.N='Dispersion SR2 - MeOH@P2-SWCNTs';
DATA_20250131.SWFL514G.N='Dispersion SWF - D2O@P2-SWCNTs (Salome)';
DATA_20250131.SWFL514R.N='Dispersion SWF - D2O@P2-SWCNTs (Salome)';
DATA_20250403.FFL514D.N='FlatField';
DATA_20250403.FFL514G.N='FlatField';
DATA_20250403.LL514A.N='Laser';
DATA_20250403.P12L514D.N='Powder S12 - TTF@P2-SWCNTs (Reflux in MeOH)';
DATA_20250403.P12L514G.N='Powder S12 - TTF@P2-SWCNTs (Reflux in MeOH)';
DATA_20250403.P13L514D.N='Powder S13 - TTF@P2-SWCNTs (Gas phase)';
DATA_20250403.P13L514G.N='Powder S13 - TTF@P2-SWCNTs (Gas phase)';
DATA_20250403.P14L514D.N='Powder S14 - C_{16}H_{34}@P2-SWCNTs (Liquid phase)';
DATA_20250403.P14L514G.N='Powder S14 - C_{16}H_{34}@P2-SWCNTs (Liquid phase)';
DATA_20250403.P15L514D.N='Powder S15 - C_{12}H_{26}@P2-SWCNTs (Liquid phase)';
DATA_20250403.P15L514G.N='Powder S15 - C_{12}H_{26}@P2-SWCNTs (Liquid phase)';
DATA_20250403.PP2L514D.N='Powder P2-SWCNTs Annealed';
DATA_20250403.PP2L514G .N='Powder P2-SWCNTs Annealed';
DATA_20250404.FAL458D.N='KIT Empty CNT film (A)';
DATA_20250404.FAL458G.N='KIT Empty CNT film (A)';
DATA_20250404.FAL458R.N='KIT Empty CNT film (A)';
DATA_20250404.FBL458D .N='KIT Empty CNT film (B)';
DATA_20250404.FBL458G  .N='KIT Empty CNT film (B)';
DATA_20250404.FBL458R.N='KIT Empty CNT film (B)';
DATA_20250404.LL458A .N='Laser';
DATA_20250404.P12L458D.N='Powder S12 - TTF@P2-SWCNTs (Reflux in MeOH)';
DATA_20250404.P12L458G .N='Powder S12 - TTF@P2-SWCNTs (Reflux in MeOH)';
DATA_20250404.P12L458R .N='Powder S12 - TTF@P2-SWCNTs (Reflux in MeOH)';
DATA_20250404.P13L458D .N='Powder S13 - TTF@P2-SWCNTs (Gas phase)';
DATA_20250404.P13L458G.N='Powder S13 - TTF@P2-SWCNTs (Gas phase)';
DATA_20250404.P13L458R.N='Powder S13 - TTF@P2-SWCNTs (Gas phase)';
DATA_20250404.P14L458D  .N='Powder S14 - C_{16}H_{34}@P2-SWCNTs (Liquid phase)';
DATA_20250404.P14L458G  .N='Powder S14 - C_{16}H_{34}@P2-SWCNTs (Liquid phase)';
DATA_20250404.P14L458R.N='Powder S14 - C_{16}H_{34}@P2-SWCNTs (Liquid phase)';
DATA_20250404.P15L458D.N='Powder S15 - C_{12}H_{26}@P2-SWCNTs (Liquid phase)';
DATA_20250404.P15L458G .N='Powder S15 - C_{12}H_{26}@P2-SWCNTs (Liquid phase)';
DATA_20250404.P15L458R.N='Powder S15 - C_{12}H_{26}@P2-SWCNTs (Liquid phase)';
DATA_20250404.PP2L458D.N='Powder P2-SWCNTs Annealed';
DATA_20250404.PP2L458G.N='Powder P2-SWCNTs Annealed';
DATA_20250404.PP2L458R.N='Powder P2-SWCNTs Annealed';
DATA_20250404.SR0L458G.N='Empty Arc SWCNTs (SF6)';
DATA_20250404.SR0L458R.N='Empty Arc SWCNTs (SF6)';
DATA_20250411.LL476A.N='Laser';
DATA_20250411.LL476B.N='Laser';
DATA_20250411.P12L476D.N='Powder S12 - TTF@P2-SWCNTs (Reflux in MeOH)';
DATA_20250411.P12L476G.N='Powder S13 - TTF@P2-SWCNTs (Gas phase)';
DATA_20250411.P13L476D.N='Powder S13 - TTF@P2-SWCNTs (Gas phase)';
DATA_20250411.P13L476G.N='Powder S13 - TTF@P2-SWCNTs (Gas phase)';
DATA_20250411.PP2L476D.N='Powder P2-SWCNTs Annealed';
DATA_20250411.PP2L476G.N='Powder P2-SWCNTs Annealed';
DATA_20250411.S11L476D.N='Dispersion S11 - TTF@P2-SWCNTs (Melt)';
DATA_20250411.S11L476G.N='Dispersion S11 - TTF@P2-SWCNTs (Melt)';
DATA_20250411.S11L476R.N='Dispersion S11 - TTF@P2-SWCNTs (Melt)';
DATA_20250411.S12L476D.N='Dispersion S12 - TTF@P2-SWCNTs (Reflux)';
DATA_20250411.S12L476G.N='Dispersion S12 - TTF@P2-SWCNTs (Reflux)';
DATA_20250411.S12L476R.N='Dispersion S12 - TTF@P2-SWCNTs (Reflux)';
DATA_20250411.S13L476D.N='Dispersion S13 - TTF@P2-SWCNTs (GasPhase)';
DATA_20250411.S13L476G.N='Dispersion S13 - TTF@P2-SWCNTs (GasPhase)';
DATA_20250411.S13L476R.N='Dispersion S13 - TTF@P2-SWCNTs (GasPhase)';
DATA_20250411.S14L476D.N='Dispersion S14 - C_{16}H_{34}@P2-SWCNTs (Liquid)';
DATA_20250411.S14L476G.N='Dispersion S14 - C_{16}H_{34}@P2-SWCNTs (Liquid)';
DATA_20250411.S14L476R.N='Dispersion S14 - C_{16}H_{34}@P2-SWCNTs (Liquid)';
DATA_20250411.S15L476D.N='Dispersion S15 - C_{12}H_{26}@P2-SWCNTs (Liquid)';
DATA_20250411.S15L476G.N='Dispersion S15 - C_{12}H_{26}@P2-SWCNTs (Liquid)';
DATA_20250411.S15L476R.N='Dispersion S15 - C_{12}H_{26}@P2-SWCNTs (Liquid)';
DATA_20250411.SR0L476D.N='Dispersion SR0 - Empty@P2-SWCNTs';
DATA_20250411.SR0L476G.N='Dispersion SR0 - Empty@P2-SWCNTs';
DATA_20250411.SR0L476R.N='Dispersion SR0 - Empty@P2-SWCNTs';
DATA_20250411.SWFL476D.N='Dispersion SWF - D2O@P2-SWCNTs (SF6)';
DATA_20250411.SWFL476G.N='Dispersion SWF - D2O@P2-SWCNTs (SF6)';
DATA_20250411.SWFL476R.N='Dispersion SWF - D2O@P2-SWCNTs (SF6)';
DATA_20250520.FFH514R.N='FlatField';
DATA_20250520.LL514A.N='Laser';
DATA_20250520.LL514B.N='Laser';
DATA_20250520.LLH514A.N='Laser';
DATA_20250520.PP2514G.N='P2 Annealed Powder';
DATA_20250520.PP2H514R.N='P2 Annealed Powder';
DATA_20250520.PP2L514D.N='P2 Annealed Powder';
DATA_20250520.PP2L514R.N='P2 Annealed Powder';
DATA_20250520.PTTF514G.N='TTF Powder';
DATA_20250520.R12H514R.N='Powder S12 TTF@P2-SWCNTs (Reflux) Rinsed';
DATA_20250520.R12L514D.N='Powder S12 TTF@P2-SWCNTs (Reflux) Rinsed';
DATA_20250520.R12L514G.N='Powder S12 TTF@P2-SWCNTs (Reflux) Rinsed';
DATA_20250520.R12L514R.N='Powder S12 TTF@P2-SWCNTs (Reflux) Rinsed';
DATA_20250520.R13H514R.N='Powder S13 TTF@P2-SWCNTs (GP) Rinsed';
DATA_20250520.R13L514D.N='Powder S13 TTF@P2-SWCNTs (GP) Rinsed';
DATA_20250520.R13L514G.N='Powder S13 TTF@P2-SWCNTs (GP) Rinsed';
DATA_20250520.R13L514R.N='Powder S13 TTF@P2-SWCNTs (GP) Rinsed';
DATA_20250520.R16H514R.N='Powder S16 TTF@P2-SWCNTs (GP/Melt) Rinsed';
DATA_20250520.R16L514D.N='Powder S16 TTF@P2-SWCNTs (GP/Melt) Rinsed';
DATA_20250520.R16L514G.N='Powder S16 TTF@P2-SWCNTs (GP/Melt) Rinsed';
DATA_20250520.R16L514R.N='Powder S16 TTF@P2-SWCNTs (GP/Melt) Rinsed';
DATA_20250520.R17H514R.N='Powder S17 TTF@P2-SWCNTs (GP) Rinsed';
DATA_20250520.R17L514D.N='Powder S17 TTF@P2-SWCNTs (GP) Rinsed';
DATA_20250520.R17L514G.N='Powder S17 TTF@P2-SWCNTs (GP) Rinsed';
DATA_20250520.R17L514R.N='Powder S17 TTF@P2-SWCNTs (GP) Rinsed';
DATA_20250520.U16B514G.N='Powder S16 TTF@P2-SWCNTs (GP/Melt) Unrinsed Verify';
DATA_20250520.U16L514D.N='Powder S16 TTF@P2-SWCNTs (GP/Melt) Unrinsed';
DATA_20250520.U16L514G.N='Powder S16 TTF@P2-SWCNTs (GP/Melt) Unrinsed';
DATA_20250520.U17L514D.N='Powder S17 TTF@P2-SWCNTs (GP) Unrinsed';
DATA_20250520.U17L514G.N='Powder S17 TTF@P2-SWCNTs (GP) Unrinsed';
DATA_20250618.LL514A.N='Laser';
DATA_20250618.LL514B.N='Laser';
DATA_20250618.PAPL514D.N='Powder AP-SWCNTs';
DATA_20250618.PAPL514G.N='Powder AP-SWCNTs';
DATA_20250618.PP2L514D.N='P2 Annealed Powder (Same as Salome)';
DATA_20250618.PP2L514G .N='P2 Annealed Powder (Same as Salome)';
DATA_20250618.R19A514D.N='Powder S19 - TEMED + AP Closed CNTs (Rinsed x1)';
DATA_20250618.R19A514G.N='Powder S19 - TEMED + AP Closed CNTs (Rinsed x1)';
DATA_20250618.R20A514D.N='Powder S20 - TTF@P2-SWCNTs (GP) (Rinsed x1)';
DATA_20250618.R20A514G.N='Powder S20 - TTF@P2-SWCNTs (GP) (Rinsed x1)';
DATA_20250618.R20B514D.N='Powder S20 - TTF@P2-SWCNTs (GP) (Rinsed x2)';
DATA_20250618.R20B514G.N='Powder S20 - TTF@P2-SWCNTs (GP) (Rinsed x2)';
DATA_20250618.R20C514D.N='Powder S20 - TTF@P2-SWCNTs (GP) (Rinsed x3)';
DATA_20250618.R20C514G.N='Powder S20 - TTF@P2-SWCNTs (GP) (Rinsed x3)';
DATA_20250618.R20D514D.N='Powder S20 - TTF@P2-SWCNTs (GP) (Rinsed x4)';
DATA_20250618.R20D514G.N='Powder S20 - TTF@P2-SWCNTs (GP) (Rinsed x4)';
DATA_20250618.R20E514D.N='Powder S20 - TTF@P2-SWCNTs (GP) (Rinsed x5)';
DATA_20250618.R20E514G.N='Powder S20 - TTF@P2-SWCNTs (GP) (Rinsed x5)';


% DATA_20250618.R20A514D.N='Powd.TTF@P2-SWCNTs (Rinsed x1)';
% DATA_20250618.R20A514G.N='Powd.TTF@P2-SWCNTs (Rinsed x1)';
% DATA_20250618.R20B514D.N='Powd.TTF@P2-SWCNTs (Rinsed x2)';
% DATA_20250618.R20B514G.N='Powd.TTF@P2-SWCNTs (Rinsed x2)';
% DATA_20250618.R20C514D.N='Powd.TTF@P2-SWCNTs (Rinsed x3)';
% DATA_20250618.R20C514G.N='Powd.TTF@P2-SWCNTs (Rinsed x3)';
% DATA_20250618.R20D514D.N='Powd.TTF@P2-SWCNTs (Rinsed x4)';
% DATA_20250618.R20D514G.N='Powd.TTF@P2-SWCNTs (Rinsed x4)';
% DATA_20250618.R20E514D.N='Powd.TTF@P2-SWCNTs (Rinsed x5)';
% DATA_20250618.R20E514G.N='Powd.TTF@P2-SWCNTs (Rinsed x5)';

DATA_20250618.S19A514D.N='Powder S19 - TEMED@SWCNTs (Rinsed x1)';
DATA_20250618.S19A514G.N='Powder S19 - TEMED@SWCNTs (Rinsed x1)';

DATA_20250618.TCNQ514D.N='TCNQ Powder';
DATA_20250618.TCNQ514G.N='TCNQ Powder';
DATA_20250618.U18L514D.N='Powder S18 - TCNQ@P2-SWCNTs (GP) (Unrinsed)';
DATA_20250618.U18L514G.N='Powd.TTF@P2-SWCNTs (Unrinsed)';
DATA_20250618.U20L514D.N='Powder S20 - TTF@P2-SWCNTs (GP) (Unrinsed)';
DATA_20250618.U20L514G.N='Powd.TTF@P2-SWCNTs (Unrinsed)';
DATA_20250718.B1L514D.N='Batch 1';
DATA_20250718.B1L514G.N='Batch 1';
DATA_20250718.B2L514D.N='Batch 2';
DATA_20250718.B2L514G.N='Batch 2';
DATA_20250718.B3L514D.N='Batch 3';
DATA_20250718.B3L514G.N='Batch 3';
DATA_20250718.FFL514D.N='FlatField';
DATA_20250718.FFL514G.N='FlatField';
DATA_20250718.LL514A.N='Laser';
DATA_20250718.R18A514G.N='Powder S18 - TCNQ@P2-SWCNTs (GP) (Rinsed x1)';
DATA_20250718.R18B514G.N='Powder S18 - TCNQ@P2-SWCNTs (GP) (Rinsed x1)';
DATA_20250718.R18C514D.N='Powder S18 - TCNQ@P2-SWCNTs (GP) (Rinsed x3)';
DATA_20250718.R18C514G.N='Powder S18 - TCNQ@P2-SWCNTs (GP) (Rinsed x2)';
DATA_20250718.R18D514D.N='Powder S18 - TCNQ@P2-SWCNTs (GP) (Rinsed x4)';
DATA_20250718.R18D514G.N='Powder S18 - TCNQ@P2-SWCNTs (GP) (Rinsed x2)';
DATA_20250718.R18E514D.N='Powder S18 - TCNQ@P2-SWCNTs (GP) (Rinsed x5)';
DATA_20250718.R18E514G.N='Powder S18 - TCNQ@P2-SWCNTs (GP) (Rinsed x3)';
DATA_20250718.R18F514G.N='Powder S18 - TCNQ@P2-SWCNTs (GP) (Rinsed x3)';
DATA_20250718.R18G514G.N='Powder S18 - TCNQ@P2-SWCNTs (GP) (Rinsed x4)';
DATA_20250718.R18H514G.N='Powder S18 - TCNQ@P2-SWCNTs (GP) (Rinsed x5)';
DATA_20250718.R18I514G.N='Powder S18 - TCNQ@P2-SWCNTs (GP) (Rinsed x5)';
DATA_20250718.TCNQ514G.N='TCNQ Powder';
DATA_20250718.U18L514G.N='Powder U18 - TCNQ@P2-SWCNTs (GP) (Unrinsed)';

DATA_20250731.R21A514G.N='Powder S21 - TTF@P2-SWCNTs (Melt) (Rinsed x1)';
DATA_20250731.R21E514G.N='Powder S21 - TTF@P2-SWCNTs (Melt) (Rinsed x5)';
DATA_20250731.U22L514G.N='Powder S22 - TTF@P2-SWCNTs (Reflux) (R0 Filtered)';
DATA_20250731.R22C514G.N='Powder S22 - TTF@P2-SWCNTs (Reflux) (Rinsed x3)';
DATA_20250731.U23L514G.N='Powder S23 - TCNQ@P2-SWCNTs (GP) (Unrinsed)';
DATA_20250731.R23F514G.N='Powder S23 - TCNQ@P2-SWCNTs (GP) (Rinsed x6)';
DATA_20250731.R23G514G.N='Powder S23 - TCNQ@P2-SWCNTs (GP) (Rinsed x6)';
DATA_20250801.B3L514D.N='Batch 3';
DATA_20250801.B3L514G.N='Empty (Undoped) P2-SWCNTs';
DATA_20250801.B3L514R.N='Empty (Undoped) P2-SWCNTs';
DATA_20250801.LL514A.N='Laser';
DATA_20250801.R21A514D.N='Powder S21 - TTF@P2-SWCNTs (Melt) (Rinsed x1)';
DATA_20250801.R21A514G.N='Powder S21 - TTF@P2-SWCNTs (Melt) (Rinsed x1)';
DATA_20250801.R21A514R.N='Powder S21 - TTF@P2-SWCNTs (Melt) (Rinsed x1)';
DATA_20250801.R23F514D.N='Powder S23 - TCNQ@P2-SWCNTs (GP) (Rinsed x6)';
DATA_20250801.R23F514G.N='Powder S23 - TCNQ@P2-SWCNTs (GP) (Rinsed x6)';
DATA_20250801.R23F514H.N='Powder S23 - TCNQ@P2-SWCNTs (GP) (Rinsed x6)';
DATA_20250801.R23F514N.N='Powder S23 - TCNQ@P2-SWCNTs (GP) (Rinsed x6)';
DATA_20250801.R23F514R.N='Powder S23 - TCNQ@P2-SWCNTs (GP) (Rinsed x6)';
DATA_20250801.TCNQ514D.N='TCNQ Powder';
DATA_20250801.TCNQ514G.N='TCNQ Powder';
DATA_20250801.TCNQ514H.N='TCNQ Powder';
DATA_20250801.TCNQ514N.N='TCNQ Powder';
DATA_20250801.TCNQ514R.N='TCNQ Powder';
DATA_20250801.TTFL514D.N='TTF Powder';
DATA_20250801.TTFL514G.N='TTF Powder';
DATA_20250801.TTFL514R.N='TTF Powder';
DATA_20250801.U22L514D.N='Powder S22 - TTF@P2-SWCNTs (Reflux) (R0 Filtered)';
DATA_20250801.U22L514G.N='Powder S22 - TTF@P2-SWCNTs (Reflux) (R0 Filtered)';
DATA_20250801.U22L514R.N='Powder S22 - TTF@P2-SWCNTs (Reflux) (R0 Filtered)';

%% Full overview of powdered samples measured at 514.5nm

close all;
% G and 2D Band Plots          
        
G = {
%     DATA_20250131.PAPL514G
%     DATA_20250618.PAPL514G


%     DATA_20250718.B1L514G
%     DATA_20250718.B2L514G
%     DATA_20250801.B3L514G
%     DATA_20250718.B3L514G
%     DATA_20250731.B3L514G
%     DATA_20250131.P2AG
%     DATA_20250520.PP2514G
%     DATA_20250131.PP2L514G 
%     DATA_20250403.PP2L514G 
%     DATA_20250618.PP2L514G 


    DATA_20250801.B3L514G
    
    DATA_20250131.P12L514G
    DATA_20250131.P12G
    DATA_20250403.P12L514G
    DATA_20250520.R12L514G

    DATA_20250403.P13L514G
    DATA_20250520.R13L514G
%     
%     DATA_20250403.P14L514G
%     DATA_20250403.P15L514G
% 
    DATA_20250520.U16L514G
    DATA_20250520.U16B514G
    DATA_20250520.R16L514G
%     
    DATA_20250520.U17L514G
    DATA_20250520.R17L514G
% 
%     DATA_20250618.U18L514G
%     DATA_20250718.U18L514G
%     DATA_20250718.R18A514G
%     DATA_20250718.R18B514G
%     DATA_20250718.R18C514G
%     DATA_20250718.R18D514G
%     DATA_20250718.R18E514G
%     DATA_20250718.R18F514G
%     DATA_20250718.R18G514G
%     DATA_20250718.R18H514G
%     DATA_20250718.R18I514G
% 
%     DATA_20250618.S19A514G
% 
    DATA_20250618.U20L514G
    DATA_20250618.R20A514G
    DATA_20250618.R20B514G
    DATA_20250618.R20C514G
    DATA_20250618.R20D514G
    DATA_20250618.R20E514G
    
    DATA_20250731.R21A514G
    DATA_20250801.R21A514G
    DATA_20250731.R21E514G

    DATA_20250731.U22L514G
    DATA_20250801.U22L514G
    DATA_20250731.R22C514G
    
%     DATA_20250731.U23L514G
%     DATA_20250731.R23F514G
%     DATA_20250801.R23F514G

}; 



% G = FlatFieldCorrection(G, DATA_20250718.FFL514G);           
% G = FilterDataByXRange(G, 1540, 1620);     
G = FilterDataByXRange(G, 1300, 1680);           

G = RemovePolyBG(G, 0);
G = SubstractLinearBG(G, 1250, 1680);
G = Normalize(G, 1580, 1600, 'M');
% G = Normalize(G, 1577, 1578, 'M');

% DATA_20250801.TTFL514G = FilterDataByXRange(DATA_20250801.TTFL514G, 1500, 1640);           
% DATA_20250801.TTFL514G.Y = DATA_20250801.TTFL514G.Y/max(DATA_20250801.TTFL514G.Y) -2.2;

plotRaman(G, 0.00, 514);     

% hold on
% plot(DATA_20250801.TTFL514G.X, DATA_20250801.TTFL514G.Y, 'LineWidth', 1.3, 'DisplayName', 'TTF')
% xline(1589.5 , '--', 'LineWidth', 1.0,'HandleVisibility', 'off');

% mask = (DATA_20250801.B3L514G.X >= 1580) & (DATA_20250801.B3L514G.X<= 1600);
% X_fit = DATA_20250801.B3L514G.X(mask);
% Y_fit = DATA_20250801.B3L514G.Y(mask);
% lorentzian = @(b,x) b(1) * (b(2)^2 ./ ((x - b(3)).^2 + b(2)^2)) + b(4);
% [~, idxMax] = max(Y_fit);
% b0 = [max(Y_fit)-min(Y_fit), 2, (1580+1600)/2, min(Y_fit)]; 
% opts = optimset('Display','off');
% bFit = lsqcurvefit(lorentzian, b0, X_fit, Y_fit, [], [], opts);
% xPeak = bFit(3);


DD= {

DATA_20250718.B1L514D
DATA_20250718.B2L514D
DATA_20250718.B3L514D
DATA_20250801.B3L514D
DATA_20250131.P2ADD
DATA_20250131.PP2L514D
DATA_20250403.PP2L514D
DATA_20250520.PP2L514D
DATA_20250618.PP2L514D
DATA_20250131.P12DD
DATA_20250131.P12L514D
DATA_20250403.P12L514D
DATA_20250520.R12L514D
DATA_20250403.P13L514D
DATA_20250520.R13L514D
DATA_20250403.P14L514D
DATA_20250403.P15L514D
DATA_20250520.R16L514D
DATA_20250520.U16L514D
DATA_20250520.R17L514D
DATA_20250520.U17L514D
DATA_20250618.U18L514D
DATA_20250718.R18C514D
DATA_20250718.R18D514D
DATA_20250718.R18E514D
DATA_20250618.R19A514D
DATA_20250618.S19A514D
DATA_20250618.R20A514D
DATA_20250618.R20B514D
DATA_20250618.R20C514D
DATA_20250618.R20D514D
DATA_20250618.R20E514D
DATA_20250618.U20L514D
DATA_20250801.R21A514D
DATA_20250801.U22L514D
DATA_20250801.R23F514D

    
}; 



DD = FilterDataByXRange(DD, 2500, 2835);           
DD = RemovePolyBG(DD, 0);
DD = SubstractLinearBG(DD, 2500, 2835);
DD = Normalize(DD, 2500, 2700, 'M');

% close all
% plotRaman(DD, 0.25, 514);     

RBM = {
        DATA_20250131.P2ARBM
        DATA_20250131.PP2L514R
        DATA_20250520.PP2L514R
        DATA_20250801.B3L514R
        DATA_20250131.P12L514R
        DATA_20250131.P12RBM
        DATA_20250520.R12L514R
        DATA_20250520.R13L514R
        DATA_20250520.R16L514R
        DATA_20250520.R17L514R
        DATA_20250801.R21A514R
        DATA_20250801.U22L514R
        DATA_20250801.R23F514R
    };

RBM = FilterDataByXRange(RBM, 100, 530);           
RBM = RemovePolyBG(RBM, 0);
% RBM = SubstractLinearBG(RBM, 2500, 2835);
RBM = Normalize(RBM, 140, 200, 'M');

% close all
% plotRaman(RBM, 0.25, 514);     

%% DOPING VECTORS
% 
% G = {
%     DATA_20250801.B3L514G
%     DATA_20250618.U20L514G
%     DATA_20250618.R20A514G
%     DATA_20250618.R20B514G	
%     DATA_20250618.R20C514G
%     DATA_20250618.R20D514G
%     DATA_20250618.R20E514G    
%     
%     DATA_20250718.R18F514G
%     DATA_20250718.R18G514G   
%     DATA_20250718.R18H514G  
%     
%     DATA_20250131.P12L514G
%     DATA_20250403.P12L514G
%     DATA_20250403.P13L514G
%  
%     DATA_20250520.R12L514G
%     DATA_20250520.R13L514G
%     DATA_20250520.R16L514G
%     DATA_20250520.R17L514G
%     DATA_20250520.U16L514G
%     DATA_20250520.U17L514G
%     
%     DATA_20250801.R23F514G
%     DATA_20250801.U22L514G
%     
% }; 
% 
% G = FilterDataByXRange(G, 1500, 1640);           
% G = RemovePolyBG(G, 0);
% G = SubstractLinearBG(G, 1250, 1680);
% G = Normalize(G, 1580, 1600, 'M');
% plotRaman(G, 0.10, 514);    
% 
% DD= {
%     DATA_20250801.B3L514D  
%     DATA_20250618.U20L514D
%     DATA_20250618.R20A514D
%     DATA_20250618.R20B514D	
%     DATA_20250618.R20C514D
%     DATA_20250618.R20D514D
%     DATA_20250618.R20E514D
%     
%     DATA_20250718.R18C514D
%     DATA_20250718.R18D514D
%     DATA_20250718.R18E514D
%     
%     DATA_20250131.P12L514D
%     DATA_20250403.P12L514D
%     DATA_20250403.P13L514D
%     
%     DATA_20250520.R12L514D
%     DATA_20250520.R13L514D
%     DATA_20250520.R16L514D
%     DATA_20250520.R17L514D
%     DATA_20250520.U16L514D
%     DATA_20250520.U17L514D
%     
%     DATA_20250801.R23F514D  
%     DATA_20250801.U22L514D
% }; 
% 
% DD = FilterDataByXRange(DD, 2500, 2835);           
% DD = RemovePolyBG(DD, 0);
% DD = SubstractLinearBG(DD, 2500, 2835);
% DD = Normalize(DD, 2500, 2700, 'M');
% plotRaman(DD, 0.10, 514);    



% % % % FITTEDD = FitSamples(DD, 2680)
% % % % 
% % % % FD = zeros(1, length(FITTEDD));
% % % % for i=1:length(FITTEDD)
% % % % %    FITTED{i}.N
% % % %    FD(i) = FITTEDD{i}.F.Params(2)-FITTEDD{1}.F.Params(2);
% % % %    FD(i) = FITTEDD{i}.F.Params(2);
% % % % end
% % % % 
% % % % 
% % % % % FITTED
% % % % FITTEDG = FitSamples(G, 1592)
% % % % FG = zeros(1, length(FITTEDG));
% % % % NG = cell(1, length(FITTEDG));
% % % % 
% % % % for i=1:length(FITTEDG)
% % % % %    FITTED{i}.N
% % % %    FG(i) = FITTEDG{i}.F.Params(2)-FITTEDG{1}.F.Params(2);
% % % %    FG(i) = FITTEDG{i}.F.Params(2);
% % % %    NG{i} = num2str(FITTEDG{i}.N);
% % % % end
% % % % 
% % % % plot(FG,FD)
% % % % scatter(FG,FD, 50, 'k','d', 'filled')
% % % % text(FG, FD, NG, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 8);
% % % % % xlim([1590, 1593]);
% % % % % ylim([2660, 2697]);
% % % % ylabel('2D Band (cm-1)')
% % % % xlabel('G Band (cm-1)')


%% 
% 
% %%%% FITTING V2
% 
% G = {
%     DATA_20250801.B3L514G
%     DATA_20250618.U20L514G
%     DATA_20250618.R20A514G
%     DATA_20250618.R20B514G	
%     DATA_20250618.R20C514G
%     DATA_20250618.R20D514G
%     DATA_20250618.R20E514G
%     
% }; 
% 
% DD= {
%     DATA_20250801.B3L514D  
%     DATA_20250618.U20L514D
%     DATA_20250618.R20A514D
%     DATA_20250618.R20B514D	
%     DATA_20250618.R20C514D
%     DATA_20250618.R20D514D
%     DATA_20250618.R20E514D
%    
% }; 

% FITTEDG = FitSamples(G, 1592);
% % FITTEDG = FitSamples(G, [1550,1560,1590]);
% 
% FITTEDD = FitSamples(DD, 2680)
% 
% numSamples = length(FITTEDG);
% FG = zeros(1, numSamples);
% FD = zeros(1, numSamples);
% NG = cell(1, numSamples);
% 
% % % % Extract data
% for i = 1:numSamples
%     FG(i) = FITTEDG{i}.F.Params(2);
%     FD(i) = FITTEDD{i}.F.Params(2);  % Assuming you want FD values too
%     NG{i} = num2str(FITTEDG{i}.N);   % Sample names
% end
% 
% % % % Assign unique colors
% colors = lines(numSamples);  % Or use: jet(numSamples), hsv(numSamples), etc.
% 
% figure; hold on;
% for i = 1:numSamples
%     scatter(FG(i), FD(i), 50, colors(i,:), 'd', 'filled', 'DisplayName', NG{i});
% end
% 
% xlim([1590, 1593]);
% ylim([2660, 2697]);
% ylabel('2D Band (cm^{-1})');
% xlabel('G Band (cm^{-1})');
% legend();  % Adjust legend location as needed
% title('Doping Vectors');
% 
% 
% 
% %%%
% 
% plotRamanFits(FITTEDG,0.3)

%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotRamanFits(SamplesToPlot, offset)
    % Create a figure for the plot
    figure;
    
    for sampleIdx = 1:length(SamplesToPlot)
        currentSample = SamplesToPlot{sampleIdx};

        % Get sample data
        X = currentSample.X;
        Y = currentSample.Y - offset * sampleIdx;
        N = currentSample.N;
        fitCurve = currentSample.F.Fit - offset * sampleIdx;
        fitParams = currentSample.F.Params;
        numPeaks = size(fitParams, 2);

        % Plot original spectrum with offset
        plot(X, Y, 'DisplayName', N, 'LineWidth', 1.3);
        hold on;

        % Plot fitted curve
        plot(X, fitCurve, 'k', 'LineWidth', 1.5, 'DisplayName', sprintf('%s - Fit', N));

        % Plot individual Lorentzian peaks with offset
        for j = 1:numPeaks
            amp = fitParams(1, j);
            pos = fitParams(2, j);
            width = fitParams(3, j);
            peakCurve = amp ./ ((X - pos).^2 + width) - offset * sampleIdx;
            plot(X, peakCurve, 'r--', 'LineWidth', 1, ...
                'HandleVisibility', 'off'); % Don't crowd the legend
        end
    end

    % Labels and final formatting
    xlabel('Raman Shift (cm^{-1})', 'FontSize', 14);
    ylabel('Intensity (a.u.)', 'FontSize', 14);
    title('Raman Spectra with Multi-Lorentzian Fits', 'FontSize', 14);
    legend('show', 'FontSize', 11);
    grid on;
    hold off;
end

function DSListOut = BackgroundSubtractionExcludeRanges(DSList, excludeRanges)
    % BackgroundSubtractionExcludeRanges performs background subtraction using the Naumov model,
    % excluding specified ranges from the background fit.
    %
    % Inputs:
    %   DSList: Array of data structures with fields X (wavelength) and Y (absorption).
    %   excludeRanges: Nx2 matrix where each row defines [rmin, rmax] to exclude.
    %
    % Output:
    %   DSListOut: Modified DSList with background-subtracted Y values.

    DSListOut = DSList; % Initialize output as the input DSList

    for i = 1:length(DSList)
        % Extract X and Y values from the current spectrum
        xx = DSList{i}.X; % Wavelength
        yy = DSList{i}.Y; % Absorption

        % Ensure X is in ascending order
        if xx(1) > xx(end)
            xx = flip(xx); 
            yy = flip(yy); 
        end

        % Interpolate the data to ensure even spacing
        xx_interp = round(xx(1)):1:round(xx(end)); % Interpolation range
        yy_interp = interp1(xx, yy, xx_interp, 'linear'); % Interpolated Y values
        xx_interp = xx_interp'; 
        yy_interp = yy_interp';

        % Identify indices to exclude based on excludeRanges
        excludeMask = false(size(xx_interp));
        for j = 1:size(excludeRanges, 1)
            range = excludeRanges(j, :);
            excludeMask = excludeMask | (xx_interp >= range(1) & xx_interp <= range(2));
        end

        % Select only points outside the excluded ranges
        bgPoints = xx_interp(~excludeMask);
        bgY = yy_interp(~excludeMask);

        % Check if sufficient points remain for fitting
        if numel(bgPoints) < 2
            error('Not enough points outside the excluded ranges for fitting.');
        end

        % Write the background data to a temporary text file
        dataToWrite = [bgPoints, bgY];
        if isempty(dataToWrite)
            error('No data to write to the temporary file.');
        else
            disp('Writing background data to temp_data.txt');
            dlmwrite('temp_data.txt', dataToWrite);
        end

        % Optimization options for fmincon
        options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'off');

        % Perform the optimization to find A and b (initial guess: [0.2, 0.00002])
        xn = fmincon(@Naumov, [0.2, 0.00002], [], [], [], [], [], [], @conf_naumov, options); 
        % Naumov optimization using the specified points

        % Define the Naumov background model
        F_naumov = @(x) x(1) * exp(-x(2) * xx_interp);

        % Evaluate the background using the Naumov function
        BKG = F_naumov(xn);

        % Subtract the background from the interpolated Y values
        correctedY = yy_interp - BKG;

        % Update the Y values in the original X range
        DSListOut{i}.Y = yy - interp1(xx_interp, BKG, xx, 'linear', 'extrap');
        DSListOut{i}.X = xx;

        % Delete the temporary file after use
        delete('temp_data.txt');

        % Plot (optional, for debugging/visualization)
%         figure;
%         hold on;
%         plot(xx, yy, 'b', 'DisplayName', 'Original Spectrum');
%         plot(xx_interp, BKG, 'r--', 'DisplayName', 'Background Fit');
%         plot(xx, DSListOut{i}.Y, 'g', 'DisplayName', 'Corrected Spectrum');
%         legend('show');
%         title(['Background Subtraction - Spectrum ', num2str(i)]);
%         xlabel('X (Wavelength)');
%         ylabel('Y (Absorption)');
%         hold off;
    end
end

function DSListOut = BackgroundSubtractionWithSpecifiedPoints(DSList, bgPoints)
    % BackgroundSubtractionWithSpecifiedPoints performs background subtraction using the Naumov model,
    % fitting the background based on user-specified X-values (bgPoints).
    %
    % Inputs:
    %   DSList: Array of data structures with fields X (wavelength) and Y (absorption).
    %   bgPoints: Array of X-values where the background is calculated.
    %
    % Output:
    %   DSListOut: Modified DSList with background-subtracted Y values.

    DSListOut = DSList; % Initialize output as the input DSList

    % Loop over each spectrum in DSList
    for i = 1:length(DSList)
        % Extract X and Y values from the current spectrum
        xx = DSList{i}.X; % Wavelength
        yy = DSList{i}.Y; % Absorption

        % Ensure X is in ascending order
        if xx(1) > xx(end)
            xx = flip(xx); 
            yy = flip(yy); 
        end

        % Interpolate the data to ensure even spacing
        xx_interp = round(xx(1)):1:round(xx(end)); % Interpolation range
        yy_interp = interp1(xx, yy, xx_interp, 'linear'); % Interpolated Y values
        xx_interp = xx_interp'; 
        yy_interp = yy_interp';

        % Ensure bgPoints is a column vector
        bgPoints = bgPoints(:);

        % Interpolate Y-values at the specified background points
        bgY = interp1(xx_interp, yy_interp, bgPoints, 'linear', 'extrap');
        bgY = bgY(:); % Ensure bgY is also a column vector

        % Write the background data to a temporary text file
        dataToWrite = [bgPoints, bgY];
        if isempty(dataToWrite)
            error('No data to write to the temporary file.');
        else
            disp('Writing background data to temp_data.txt');
            dlmwrite('temp_data.txt', dataToWrite);
        end

        % Optimization options for fmincon
        options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'off');

        % Perform the optimization to find A and b (initial guess: [0.2, 0.00002])
        xn = fmincon(@Naumov, [0.2, 0.00002], [], [], [], [], [], [], @conf_naumov, options); 
        % Naumov optimization using the specified points

        % Define the Naumov background model
        F_naumov = @(x) x(1) * exp(-x(2) * xx_interp);

        % Evaluate the background using the Naumov function
        BKG = F_naumov(xn);

        % Subtract the background from the interpolated Y values
        correctedY = yy_interp - BKG;

        % Update the Y values in the original X range
        DSListOut{i}.Y = yy - interp1(xx_interp, BKG, xx, 'linear', 'extrap');
        DSListOut{i}.X = xx;

        % Delete the temporary file after use
        delete('temp_data.txt');

%         % Plot (optional, for debugging/visualization)
%         figure;
%         hold on;
%         plot(xx, yy, 'b', 'DisplayName', 'Original Spectrum');
%         plot(xx_interp, BKG, 'r--', 'DisplayName', 'Background Fit');
%         plot(xx, DSListOut{i}.Y, 'g', 'DisplayName', 'Corrected Spectrum');
% %         scatter(bgPoints, bgY, 'k', 'DisplayName', 'Background Points', 'filled');
%         legend('show');
%         title(['Background Subtraction - Spectrum ', num2str(i)]);
%         xlabel('X (Wavelength)');
%         ylabel('Y (Absorption)');
%         hold off;
    end
end

function DSListOut = BackgroundSubtraction(DSList, range)
    % DSList is the input array of data structures with fields X and Y
    
    DSListOut = DSList; % Initialize output as the input DSList
    
%     figure
%     hold on
    
    % Loop over each structure in DSList
    for i = 1:length(DSList)
        % Extract the X and Y values from the current data structure
        xx = DSList{i}.X; % Wavelength
        yy = DSList{i}.Y; % Absorption
        
        if xx(1) > xx(end)
            xx = flip(xx); % Flip xx to ascending order
            yy = flip(yy); % Flip yy to match xx
        end

        % Interpolate the data to ensure it is evenly spaced
        xx_interp = round(xx(1)):1:round(xx(end)); % Interpolation range
        yy_interp = interp1(xx, yy, xx_interp, 'linear'); % Interpolated Y values
        xx_interp = xx_interp'; 
        yy_interp = yy_interp';
        
               
         % Find the index for range(1) (closest value greater than or equal to range(1))
         [~, start_idx] = min(abs(xx_interp - range(1))); % Closest value to range(1)
         [~, end_idx] = min(abs(xx_interp - range(2))); % Closest value to range(2)
         
       
         % Write the temporary data to a text file
        dataToWrite = [xx_interp(start_idx:end_idx), yy_interp(start_idx:end_idx)];
        if isempty(dataToWrite)
            error('No data to write to the temporary file.');
        else
            disp('Writing data to temp_data.txt');
            dlmwrite('temp_data.txt', dataToWrite);
        end
        
        % Optimization options for fmincon
        options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'off');
        
        % Perform the optimization to find A and b (initial guess: [0.2, 0.00002])
        xn = fmincon(@Naumov, [.2, 0.00002], [], [], [], [], [], [], @conf_naumov, options); % % Approach of Naumov et al. requires two starting values for A (0.2) and b (0.002).
        
        % Define the background subtraction model
        F_naumov = @(x) double(yy_interp(start_idx:end_idx) - x(1) * exp(-x(2) * xx_interp(start_idx:end_idx))); 
        
        BKG = yy_interp(start_idx:end_idx) - F_naumov(xn);
                 
        % Update the Y values in the data structure with the background-subtracted data
        DSListOut{i}.Y = yy - interp1(xx_interp(start_idx:end_idx), BKG, xx, 'linear', 'extrap');
        DSListOut{i}.X = xx;
         
        % Delete the temporary file after use
        delete('temp_data.txt');
    end
end

function err = Naumov(x) % Calculates the difference between the absorption data and the background. The MATLAB function "fmincon" tries to minimize this difference by fitting x(1)=A and x(2)=b

A=dlmread('temp_data.txt');
c = A(:,2)-x(1)*exp(-x(2).*A(:,1));
err = double(sum(c));

end

function [c,ceq] = conf_naumov(x) % Constraint function, that forces the background to be smaller than the absorption data for every single wavelength

A=dlmread('temp_data.txt');
% Nonlinear inequality constraints
c = double(x(1)*exp(-x(2).*A(:,1))-A(:,2));
% Nonlinear equality constraints
ceq = [];
end

function filteredSamples = FilterDataByXRange(samplesToFilter, xMin, xMax)
    % FilterDataByXRange filters the data of each sample to include only the points within the specified X-range.
    %
    % Inputs:
    %   - samplesToFilter: Either a single structure with fields 'X' and 'Y',
    %                      or a cell array of such structures.
    %   - xMin: The minimum value of X to include in the filtered data.
    %   - xMax: The maximum value of X to include in the filtered data.
    %
    % Outputs:
    %   - filteredSamples: Same format as input (single structure or cell array),
    %                      with filtered 'X' and 'Y' values within the range [xMin, xMax].

    % Ensure input is in cell array form
    isSingle = ~iscell(samplesToFilter);
    if isSingle
        samplesToFilter = {samplesToFilter};
    end
    
    % Allocate output
    filteredSamples = cell(size(samplesToFilter));
    
    % Iterate over each sample to filter
    for sampleIdx = 1:length(samplesToFilter)
        currentSample = samplesToFilter{sampleIdx};
        
        % Find the indices of X-values within the specified range
        validIndices = currentSample.X >= xMin & currentSample.X <= xMax;
        
        % Filter the data
        filteredSample = currentSample;
        filteredSample.X = currentSample.X(validIndices);
        filteredSample.Y = currentSample.Y(validIndices);
        
        % Store result
        filteredSamples{sampleIdx} = filteredSample;
    end
    
    % Return in the same format as input
    if isSingle
        filteredSamples = filteredSamples{1};
    end
end

function SpectraList = RemoveBackgroundProfile(SpectraList, Xpoints)
    % Inicializar matrices para almacenar X y Y del fondo
    for i = 1:length(SpectraList)
        sample = SpectraList{i};
        
        % Interpolar los valores Y en los puntos Xpoints que corresponden al fondo
        Ybackground = interp1(sample.X, sample.Y, Xpoints, 'linear', 'extrap');
        Xbackground = Xpoints(:);  % Asegurar que sea un vector columna
        Ybackground = Ybackground(:);

        % Definir la función objetivo para el ajuste de mínimos cuadrados
        % La función modelo es A + B/X, y queremos minimizar la suma de cuadrados
        % con penalización en los valores negativos de la corrección
       objectiveFunc = @(params) sum((Ybackground - (params(1) + params(2)./Xbackground)).^2) + ...
                                10*sum(min(0, Ybackground - (params(1) + params(2)./Xbackground)).^2) + ...
                                10*sum(min(0, Ybackground - (params(1) + params(2)./Xbackground) - (params(1) + params(2)./Xbackground)).^2);

        % Inicializar los parámetros A y B
        initialParams = [0, 0];

        % Usar un optimizador para encontrar los parámetros A y B
        options = optimset('Display', 'off');
        params = fminsearch(objectiveFunc, initialParams, options);

        % Extraer A y B para este espectro
        A = params(1);
        B = params(2);

        % Calcular el fondo ajustado para este espectro
        background = A + B ./ sample.X;

        % Restar el fondo de los valores Y para corregir el espectro
        correctedY = sample.Y - background;

        % Actualizar el espectro corregido en la lista
        sample.Y = correctedY;
        SpectraList{i} = sample;  % Actualizar la lista de espectros
    end
end

function DSList = RemovePolyBG(DSList, degree)
    % Remove baseline from a list of data structures using polynomial fitting
    % DSList: list of structures, each with fields X (Raman shift) and Y (intensity)
    % degree: Degree of the polynomial used for baseline fitting
    
    % Iterate over each structure in the list
    for i = 1:length(DSList)
        DS = DSList{i};  % Extract the current data structure
        
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
        
        % Optionally, display the polynomial coefficients for debugging
%         disp(['Structure ', num2str(i), ' Polynomial Coefficients: ', num2str(p)]);
        
        % Save the updated structure back to the list
        DSList{i} = DS;
    end
end


