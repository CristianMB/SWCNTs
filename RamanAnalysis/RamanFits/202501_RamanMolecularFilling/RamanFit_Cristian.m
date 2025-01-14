clear all;
clc
addpath('X:\SWCNTs');
import UsefulFunctions.*;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%       PROGRAM TO FIT MULTIPLE RAMAN DATA      %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%       INPUTS!      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Import the experimental data files
% xdata = Raman shifts (all data should be on the same x-axis)
% ydata = Raman spectra (can be matrix with in each column another
% spectrum)
% usually I create a .mat file and read in this .mat file

load FilledB.mat

figure;
plot(xdata,ydata)
title('Raw data');
xlabel('Raman shift (cm-1)'); 
ylabel('Intensity (a.u.)');

%% %%%%%%%%%%%%%%%%%%%%%%%%       INPUTS fit paramters     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create a new .xlsx file. Enter number of peaks and starting parameters of each peak, the parameters 
%of the Lorentz function x0=[xc1,w1,xc2,w2,...], where xci is the center 
%point of the i peak and wi its width. Set the lower & upper boundaries (LB & UB)
%for each parameter. NOTE! If no boundaries are desired to be set, input 
% '-inf' on the corresponding position, if you want to "fix" a parameter,
% put the LB and UB the same as your parameter

% Input = xlsread('InputParameters.xlsx','References','A2:D100');
Input = xlsread('InputParameters.xlsx','ReferencesFixed','A2:D100');
NumLor = Input(1,1); % number of lorentzians
x0=Input(:,2); 
LB=Input(:,3); 
UB=Input(:,4);
NumSpec=size(ydata,2);% number of spectra to fit simultaneously

% Check if the initial parameters are good (adjust if needed)

FITinital=RamanFit(xdata,ydata,x0,NumSpec,NumLor);
       
%We plot it to test our guess
figure;
plot(xdata,ydata,'-',xdata,FITinital,'k');
title('Initial guess');
xlabel('Raman shift (cm-1)'); 
ylabel('Intensity (a.u.)');
 


%% FITTING SECTION. we use lsqcurvefit
%experimental data. fminsearchbnd(fun,x0,LB,UB,options,varargin)
% options=optimset('MaxFunEvals',2e15,'MaxIter',2e15,'TolX',1e-20);
options = optimset('MaxFunEvals',1e4,'MaxIter',100,'Display','iter','TolX',1e-10,'TolFun',1e-10);
[FittedParams,chi2,residual,~,~,~,Jacobian] =  lsqcurvefit(@(x0,xdata)RamanFit(xdata,ydata,x0,NumSpec,NumLor),x0,xdata,ydata,LB,UB,options);
[FIT,L,A] = RamanFit(xdata,ydata,FittedParams,NumSpec,NumLor);

figure;
plot(xdata,ydata,'k',xdata,FIT,'r')
title('SecondPlot FittingSection');
xlabel('Raman shift (cm-1)'); 
ylabel('Intensity (a.u.)');

% Combine FittedParams with their errors for clarity
FittedDataWithErrors = FittedParams(:); 
% Define where to write te data in Excel
outputFile = 'FittingResults.xlsx';
sheetName = 'References';
outputRange = 'A2'; % Starting from column H, row 2
% Write the data into the Excel sheet
writematrix(FittedDataWithErrors, outputFile, 'Sheet', sheetName, 'Range', outputRange);
disp('Fitted parameters and errors have been written to Excel.');


%% PLOT SECTION
%Create the fitted function and plot it
figure;
title('FinalPlot');
xlabel('Raman shift (cm-1)');
ylabel('Intensity (a.u.)')

hold on

for i = 1:NumSpec
%Plot the data and the fit
    plot(xdata,ydata(:,i),'.k',xdata,FIT(:,i),'g');
    xlim([min(xdata) max(xdata)]);
%Plot the individual components of the sum of lorentzians to check the fits  
    for k = 1:2:NumLor-1 % these are the empty peaks
        IndivComp = L(:,k)*A(k,i)+L(:,end-1)*A(end-1,i)+L(:,end)*A(end,i);
        plot(xdata,IndivComp,'color','r');
        IndCompSave(:,i,k) = IndivComp;
        
    end
    for k = 2:2:NumLor % these are the water-filled peaks
        IndivComp = L(:,k)*A(k,i)+L(:,end-1)*A(end-1,i)+L(:,end)*A(end,i);
        plot(xdata,IndivComp,'color','b');
        IndCompSave(:,i,k) = IndivComp;
        
    end
end


%% Make error bars on the fitted parameters
ci=nlparci(FittedParams,residual,'Jacobian',Jacobian); % gives two-sigma upper and lower boundaries of the fit parameters
error=max(abs(FittedParams-ci(:,1)),abs(FittedParams-ci(:,2)))/2; % take the maximum of the difference, and then divide by 2 to have one sigma errors.

