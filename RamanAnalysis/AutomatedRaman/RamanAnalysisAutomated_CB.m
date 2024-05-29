clearvars; 
close all; 

%% Input data:
% The program expects multiple (up to 10) files in .txt or .dat format. In each file, the Raman shift
% should be included in the first column, followed by the Intensity data in
% consecutive columns (as many spectra available, and not neccesarily the same number of spectra per file).

% The location and names of the files should be indicated in path and file_name
% variables, respectively. The total number of files should be included in
% total.
addpath = 'C:\Users\cborja\OneDrive - Universiteit Antwerpen\Measurements Data\Raman\Automated\';
file_name={ [addpath 'DATA']    };
total=1; % Total number of files
type = '.txt';
name={'A'}; % Indicate here the names of the data for legend, titles and data output. No spaces allowed.


%% Model: Choose the method to analyse spectral features 

rbm=1;% Set to 1 if RBM analysis desired
lorentz=1; % lorentz fits for peaks 1, 2 and 3 (=0 for simple identification; =1 for lorentzian peak fitting)
nt=0; % Splitting peak 1 in 1+ and 1-? (=0 for no; =1 for yes), note only valid if lorentz=1
map=0; % Set to 1 if spatial mapping desired

%% Normalization (mandatory): 
% Choose normalization spectral range. The spectra will be normalised to
% the maximun within the specified spectral range: 
normLow=100; %Lower limit in cm-1
normHigh=220; %Upper limit in cm-1

%% Peak identification: Intensity and Shifts (mandatory)
% Select appropriate spectral range where the full peaks are resolved. 

% Peak 1:  * Note that peak 1 can be optionally fitted into 2 Lorentzian
band1Name='G';
band1Low=1470; %Lower limit in cm-1
band1High=1650; %Upper limit in cm-1
% Peak 2: 
band2Name='D';
band2Low=1200; %Lower limit in cm-1
band2High=1400; %Upper limit in cm-1
% Peak 3: 
band3Name='2D';
band3Low=2400; %Lower limit in cm-1
band3High=2900; %Upper limit in cm-1

% RBM modes (only if RBM=1)
RBMregion_Low=90; % Lower limit in cm-1
RBMregion_High=280; % Upper limit in cm-1
Prom=[0.01]; % This value sets the max limit at which peaks will be considered. 
%Local maxima with a maximum height lower than Prom will be discarded.

%% Fitting parameters (only if lorentz=1)

% If lorentzian fits are selected, the user is required to indicate
% initial guess for the peak position (cm-1), FWHM (cm-1) and intensity

%Peak 1: 
Init_peak1min=[1530 25 0.3; 1530 20 0.20]; %[center FWHM intensity_max] for peak 1- (Only if nt=1)
Init_peak1plus=[1578 30 1; 1580 30 1]; %[center FWHM intensity_max] for peak 1+ (or 1 if nt=0)
%Peak 2: 
Init_peak2=[1300 30 0.08; 1305 40 0.06]; %[center FWHM intensity_max] for peak 2
%Peak 3: 
Init_peak3=[2590 70 0.4; 2640 50 0.3]; %[center FWHM intensity_max] for peak 3

%% Mapping options (only if map=1)
% Select here if mapping is wanted, and indicate the details of the map

col=32 ; % number of columns / pixels in x
raws=32 ; % number of raws / pixels in y

colmum=32.5; % micrometers of the map in the x axis
rawsmum=24; % micrometers of the map in the y axis

%% Ploting options/Output plots:
% Choose the desired output plots (1 is yes and 0 is no). 
raw=1; % Make figure af all raw spectra.
norm=1; % Plot all normalised spectra.
range=1; % Plot spectral regions chosen for peaks 1, 2 and 3 in intensity calculation
peaks=1; % plot found peaks in the RBM region, Note that if the
% number of spectra is very high, the computing is going to slow down
% significantly.
% Plot correlations between spectral features
correlation1=0; % position 1 vs I21
correlation2=0; % position 2 vs I21
correlation3=0; % position 3 vs I21
correlation4=0; % position 3 vs position 1
correlation5=0; % position vs FWHM for all bands

% Plot Raman maps. By default, is map=1 (line 27), the program will plot
% maps of all peaks positions and intensity ratio I21. In addition, the
% user might select:
map1=0; % FWHM peak 1
map2=0; % FWHM peak 2
map3=0; % FWHM peak 3
map4=0; % I3/I1

% Width of histograms
width=3; % width of Raman shift histograms 
width_fw=5; % width of FWHM histograms 
width_int=0.03; %width of intensity ratio histograms 

% Saving figures;
saveFigs=0; % All figures will be saved in .fig format in the current datafolder
nameFigs='Figures'; % Basename used for the saving of figures. Additional labels will be attached to the basename to identify the figure

%% Code Starts here:

if lorentz==0
    nt=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%% Set distribution of subfigures %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if total==1
    div=[1 1];
elseif total==2
    div=[1 2];
elseif total==3
    div=[1 3];
elseif total==4
    div=[2 2];
elseif total==5
    div=[2 3];
elseif total==6
    div=[2 3];
elseif total==10
    div=[2 5];
else total==7
    div=[3 3];
end

for z=1:total

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Import data %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IN = importdata([file_name{z} type]);
Shift(:) = IN(:,1); % Raman shift is the first column of the file
Intensity(:,:)=IN(:,2:end); % Raman intensities are all data from column 2 til end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%% Set color scales for plots %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmap = {gray(length(Intensity(1,:))+20);
winter(length(Intensity(1,:))+5);
summer(length(Intensity(1,:))+5);
autumn(length(Intensity(1,:))+5);
parula(length(Intensity(1,:))+5);
cool(length(Intensity(1,:))+5);
spring(length(Intensity(1,:))+5);
copper(length(Intensity(1,:))+5);
pink(length(Intensity(1,:))+5);
hot(length(Intensity(1,:))+5)}; % Sets color maps for different files

color=cmap{z};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%% Set distribution of subfigures %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if total==1
    div=[1 1];
elseif total==2
    div=[1 2];
elseif total==3
    div=[1 3];
elseif total==4
    div=[2 2];
elseif total==5
    div=[2 3];
elseif total==6
    div=[2 3];
elseif total==10
    div=[2 5];
else total==7
    div=[3 3];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Normalization %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, indLow]=min(abs(Shift-normLow));
[~, indHigh]=min(abs(Shift-normHigh));

Intensity_norm=zeros(size(Intensity));

for n=1:length(Intensity(1,:))
Intensity_norm(:,n)=Intensity(:,n)-min(Intensity(:,n)); % substract min
Intensity_norm(:,n)=Intensity_norm(:,n)./max(Intensity_norm(:,n)); % divide by norm
end

Intensity_av=mean(Intensity_norm'); % Calculate average spectra

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Spectral features calculation %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define spectral range for G, D and 2D
[~ ,ind1Low]=min(abs(Shift-band1Low));
[~ ,ind1High]=min(abs(Shift-band1High));

[~, ind2Low]=min(abs(Shift-band2Low));
[~, ind2High]=min(abs(Shift-band2High));

[~, ind3Low]=min(abs(Shift-band3Low));
[~, ind3High]=min(abs(Shift-band3High));

% Define Intensity and shift vectors

Intensity_1range=zeros(ind1High-ind1Low+1,length(Intensity_norm(1,:)));
Shift_range1=Shift(ind1Low:ind1High);
Intensity_2range=zeros(ind2High-ind2Low+1,length(Intensity_norm(1,:)));
Shift_range2=Shift(ind2Low:ind2High);
Intensity_3range=zeros(ind3High-ind3Low+1,length(Intensity_norm(1,:)));
Shift_range3=Shift(ind3Low:ind3High);
Int_1=zeros(1,length(Intensity_norm(1,:)));
Int_1min=zeros(1,length(Intensity_norm(1,:)));
Int_1plus=zeros(1,length(Intensity_norm(1,:)));
Int_2=zeros(1,length(Intensity_norm(1,:)));
Int_3=zeros(1,length(Intensity_norm(1,:)));
center_1=zeros(1,length(Intensity_norm(1,:)));
center_2=zeros(1,length(Intensity_norm(1,:)));
center_3=zeros(1,length(Intensity_norm(1,:)));
center_1min=zeros(1,length(Intensity_norm(1,:)));
center_1plus=zeros(1,length(Intensity_norm(1,:)));
FWHM_1=zeros(1,length(Intensity_norm(1,:)));
FWHM_2=zeros(1,length(Intensity_norm(1,:)));
FWHM_3=zeros(1,length(Intensity_norm(1,:)));
FWHM_1min=zeros(1,length(Intensity_norm(1,:)));
FWHM_1plus=zeros(1,length(Intensity_norm(1,:)));


for n=1:length(Intensity_norm(1,:))   
   
   Intensity_1range(:,n)=Intensity_norm(ind1Low:ind1High,n);
   Intensity_2range(:,n)=Intensity_norm(ind2Low:ind2High,n);
   Intensity_3range(:,n)=Intensity_norm(ind3Low:ind3High,n);
    
   %%
   %%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%% Method 1 %%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%
   if lorentz==0;
   % peak 1/G: 
   Int_1(n)=max(Intensity_1range(:,n))-min(Intensity_1range(:,n)); % Calculate intensity
   [~, indMax1]=max(Intensity_1range(:,n));
   center_1(n)=Shift_range1(indMax1); % Calculate position
   
   % peak 2/D: 
   Int_2(n)=max(Intensity_2range(:,n))-min(Intensity_2range(:,n));
   [~, indMax2]=max(Intensity_2range(:,n));
   center_2(n)=Shift_range2(indMax2);
   
   % peak 3/2D: 
   Int_3(n)=max(Intensity_3range(:,n))-min(Intensity_3range(:,n));
   [~, indMax3]=max(Intensity_3range(:,n));
   center_3(n)=Shift_range3(indMax3);
   
   end
  %% 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%% Method 2: Lorentzian fitting %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
   if lorentz==1;
    
    fit_func=@(gamma,x)(gamma(1)./((x-gamma(2)).^2+gamma(3))+gamma(4));% Single lorentz fitting function
    
    % Peak 2: D band
    lb=[0 band2Low 0 0]; 
    ub=[4e8 band2High 4e4 1000];
    InitGuess_2=[Init_peak2(z,3)*(Init_peak2(z,2)/2)^2 Init_peak2(z,1) (Init_peak2(z,2)/2)^2 0];
    Intensity_2range(:,n)=Intensity_norm(ind2Low:ind2High,n);
    gamma01=InitGuess_2(:);
    gamma_2 = zeros(length(gamma01), 1);
    %[gamma_2(:),R_2] = nlinfit(Shift_range2',Intensity_2range(:,n),fit_func,gamma01);  
    [gamma_2(:),R_2] = lsqcurvefit(fit_func,gamma01,Shift_range2',Intensity_2range(:,n),lb,ub); 
    R2_2(n)= 1 - sum(R_2.^2)/sum(((Intensity_2range(:,n)-mean(Intensity_2range(:,n))).^2)); 

    
    I_2_fit(:,n)=(gamma_2(1)./((Shift_range2(:)'-gamma_2(2)).^2+gamma_2(3))+gamma_2(4));
    center_2(n)=gamma_2(2);
    FWHM_2(n)=real(2*sqrt(gamma_2(3)));
    Int_2(n)=gamma_2(1)/gamma_2(3);  
       
   % Peak 3: 2D band
    lb=[0 band3Low 0 0]; 
    ub=[4e8 band3High 4e4 1000];
    InitGuess_3=[Init_peak3(z,3)*(Init_peak3(z,2)/2)^2 Init_peak3(z,1) (Init_peak3(z,2)/2)^2 0];
    Intensity_3range(:,n)=Intensity_norm(ind3Low:ind3High,n);
    gamma01=InitGuess_3(:);
    gamma_3 = zeros(length(gamma01), 1);
    %[gamma_3(:),R_3] = nlinfit(Shift_range3',Intensity_3range(:,n),fit_func,gamma01);  
    [gamma_3(:),R_3] =lsqcurvefit(fit_func,gamma01,Shift_range3',Intensity_3range(:,n),lb,ub); 
    R2_3(n)= 1 - sum(R_3.^2)/sum(((Intensity_3range(:,n)-mean(Intensity_3range(:,n))).^2)); 

    
    I_3_fit(:,n)=(gamma_3(1)./((Shift_range3(:)'-gamma_3(2)).^2+gamma_3(3))+gamma_3(4));
    center_3(n)=gamma_3(2);
    FWHM_3(n)=2*sqrt(gamma_3(3));
    Int_3(n)=gamma_3(1)/gamma_3(3);
    
    if nt==0
    % Peak 1: G band, single peak fitting
    lb=[0 band1Low 0 0]; 
    ub=[4e8 band1High 4e4 1000];
    InitGuess_1=[Init_peak1plus(z,3)*(Init_peak1plus(z,2)/2)^2 Init_peak1plus(z,1) (Init_peak1plus(z,2)/2)^2 0];    
    Intensity_1range(:,n)=Intensity_norm(ind1Low:ind1High,n);
    gamma01=InitGuess_1(:);
    gamma_1 = zeros(length(gamma01), 1);
    [gamma_1(:),R_1] = lsqcurvefit(fit_func,gamma01,Shift_range1',Intensity_1range(:,n),lb,ub); 
    R2_1(n)= 1 - sum(R_1.^2)/sum(((Intensity_1range(:,n)-mean(Intensity_1range(:,n))).^2)); 
    
    I_1_fit(:,n)=(gamma_1(1)./((Shift_range1(:)'-gamma_1(2)).^2+gamma_1(3))+gamma_1(4));
    center_1(n)=gamma_1(2);
    FWHM_1(n)=2*sqrt(gamma_1(3));
    Int_1(n)=gamma_1(1)/gamma_1(3); 
    
    else
    % Peak 1: G band split G+ and G-, double peak fitting
    lb=[0 band1Low 0 0 band1Low 0 0]; 
    ub=[4e8 band1High 4e4 4e8 band1High 4e4 1000];
    InitGuess_1=[Init_peak1min(z,3)*(Init_peak1min(z,2)/2)^2 Init_peak1min(z,1) (Init_peak1min(z,2)/2)^2 Init_peak1plus(z,3)*(Init_peak1plus(z,2)/2)^2 Init_peak1plus(z,1) (Init_peak1plus(z,2)/2)^2 0.1];
    Intensity_1range(:,n)=Intensity_norm(ind1Low:ind1High,n);
    fit_func=@(gamma,x)(gamma(1)./((x-gamma(2)).^2+gamma(3))+gamma(4)./((x-gamma(5)).^2+gamma(6))+gamma(7));% lorentz
    
    gamma01=InitGuess_1(:);
    gamma_1 = zeros(length(gamma01), 1);
    [gamma_1(:),R_1] = lsqcurvefit(fit_func,gamma01,Shift_range1',Intensity_1range(:,n),lb,ub); 
    R2_1(n)= 1 - sum(R_1.^2)/sum(((Intensity_1range(:,n)-mean(Intensity_1range(:,n))).^2)); 

    I_1_fit(n,:)=(gamma_1(1)./((Shift_range1(:)'-gamma_1(2)).^2+gamma_1(3))+gamma_1(4)./((Shift_range1(:)'-gamma_1(5)).^2+gamma_1(6))+gamma_1(7));
 
        if gamma_1(2)<gamma_1(5)
        center_1min(n)=gamma_1(2);
        FWHM_1min(n)=2*sqrt(gamma_1(3));
        Int_1min(n)=gamma_1(1)/gamma_1(3);
        center_1plus(n)=gamma_1(5);
        FWHM_1plus(n)=2*sqrt(gamma_1(6));
        Int_1plus(n)=gamma_1(4)/gamma_1(6);
        else
        center_1plus(n)=gamma_1(2);
        FWHM_1plus(n)=2*sqrt(gamma_1(3));
        Int_1plus(n)=gamma_1(1)/gamma_1(3);
        center_1min(n)=gamma_1(5);
        FWHM_1min(n)=2*sqrt(gamma_1(6));
        Int_1min(n)=gamma_1(4)/gamma_1(6);
        end
 
    end
   end
end

%% Calculate mean values and standard deviation

if nt==1
    
        if Int_1plus>Int_1min
    I21=Int_2./Int_1plus; % Calculate intensity ratio
    I21_av=mean(I21);
    I21_error=std(I21);

    I31=Int_3./Int_1plus; % Calculate intensity ratio
    I31_av=mean(I31);
    I31_error=std(I31);
        else
    I21=Int_2./Int_1min; % Calculate intensity ratio
    I21_av=mean(I21);
    I21_error=std(I21);

    I31=Int_3./Int_1min; % Calculate intensity ratio
    I31_av=mean(I31);
    I31_error=std(I31);       
        end

pos_1plus_av=mean(center_1plus);
pos_1plus_err=std(center_1plus);
pos_1min_av=mean(center_1min);
pos_1min_err=std(center_1min);
pos_3_av=mean(center_3);
pos_3_err=std(center_3);
pos_2_av=mean(center_2);
pos_2_err=std(center_2);

else
    
I21=Int_2./Int_1; % Calculate intensity ratio
I21_av=mean(I21);
I21_error=std(I21);

I31=Int_3./Int_1; % Calculate intensity ratio
I31_av=mean(I31);
I31_error=std(I31);

pos_1_av=mean(center_1);
pos_1_err=std(center_1);
pos_3_av=mean(center_3);
pos_3_err=std(center_3);
pos_2_av=mean(center_2);
pos_2_err=std(center_2);
end

if lorentz==1 & nt==0
FWHM_1_av=mean(FWHM_1);
FWHM_1_err=std(FWHM_1);
FWHM_3_av=mean(FWHM_3);
FWHM_3_err=std(FWHM_3);
FWHM_2_av=mean(FWHM_2);
FWHM_2_err=std(FWHM_2);   
end

if lorentz==1 & nt==1
FWHM_1plus_av=mean(FWHM_1plus);
FWHM_1plus_err=std(FWHM_1plus);   
FWHM_1min_av=mean(FWHM_1min);
FWHM_1min_err=std(FWHM_1min); 
FWHM_3_av=mean(FWHM_3);
FWHM_3_err=std(FWHM_3);
FWHM_2_av=mean(FWHM_2);
FWHM_2_err=std(FWHM_2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% RBM modes Shifts  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if rbm
[~, indRBMLow]=min(abs(Shift-RBMregion_Low));
[~, indRBMHigh]=min(abs(Shift-RBMregion_High));
PeaksInt=[];
PeaksLoc=[];
IntRBM={};
LocRBM={};

for n=1:length(Intensity_norm(1,:))

if peaks
figure(401),
subplot(div(1),div(2),z);hold on

[Shift_sorted, sort_idx] = sort(Shift);
Intensity_norm_sorted = Intensity_norm(sort_idx, n);

findpeaks(Intensity_norm_sorted, Shift_sorted, 'MinPeakProminence', Prom(z))
xlabel('Raman shift / cm^{-1}','FontSize',12);
ylabel('Intensity/a.u.','FontSize',12);
title([name{z} ' :RBM region peaks']);
set(gca,'Box','on');
%plotbrowser('on');
end

Intensity_norm_no_nan = Intensity_norm;
Intensity_norm_no_nan(isnan(Intensity_norm_no_nan)) = 0;

[pksRBM,locsRBM]=findpeaks(Intensity_norm_no_nan(:,n),Shift(:),'MinPeakProminence',Prom(z));
PeaksInt=[PeaksInt [pksRBM]'];
PeaksLoc=[PeaksLoc [locsRBM]];

IntRBM(n,:)={[pksRBM]};
LocRBM(n,:)={[locsRBM]};
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Output data %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save data in txt file

if rbm

if lorentz==0
T=table(Int_1',center_1', Int_2',center_2', Int_3',center_3',IntRBM,LocRBM,'VariableNames',{['Intensity_' band1Name],['Shift_' band1Name],['Intensity_' band2Name],['Shift_' band2Name],['Intensity_' band3Name],['Shift_' band3Name],'Intensity_RBM','Shift_RBM'});
end

if lorentz==1 & nt==1
T=table(Int_1min',center_1min',FWHM_1min',Int_1plus',center_1plus',FWHM_1plus',R2_1', Int_2',center_2',FWHM_2',R2_2', Int_3',center_3',FWHM_3',R2_3',IntRBM,LocRBM,'VariableNames',{['Intensity_' band1Name '-'],['Shift_' band1Name '-'],['FWHM_' band1Name '-'],['Intensity_' band1Name '+'],['Shift_' band1Name '+'],['FWHM_' band1Name '+'],['R^2_' band1Name],['Intensity_' band2Name],['Shift_' band2Name],['FWHM_' band2Name],['R^2_' band2Name],['Intensity_' band3Name],['Shift_' band3Name],['FWHM_' band3Name],['R^2_' band3Name],'Intensity_RBM','Shift_RBM'});
end

if lorentz==1 & nt==0
T=table(Int_1',center_1',FWHM_1',R2_1', Int_2',center_2',FWHM_2',R2_2', Int_3',center_3',FWHM_3',R2_3',IntRBM,LocRBM,'VariableNames',{['Intensity_' band1Name],['Shift_' band1Name],['FWHM_' band1Name],['R^2_' band1Name],['Intensity_' band2Name],['Shift_' band2Name],['FWHM_' band2Name],['R^2_' band2Name],['Intensity_' band3Name],['Shift_' band3Name],['FWHM_' band3Name],['R^2_' band3Name],'Intensity_RBM','Shift_RBM'});
end

else

if lorentz==0
T=table(Int_1',center_1', Int_2',center_2', Int_3',center_3','VariableNames',{['Intensity_' band1Name],['Shift_' band1Name],['Intensity_' band2Name],['Shift_' band2Name],['Intensity_' band3Name],['Shift_' band3Name]});
end

if lorentz==1 & nt==1
T=table(Int_1min',center_1min',FWHM_1min',Int_1plus',center_1plus',FWHM_1plus',R2_1', Int_2',center_2',FWHM_2',R2_2', Int_3',center_3',FWHM_3',R2_3','VariableNames',{['Intensity_' band1Name '-'],['Shift_' band1Name '-'],['FWHM_' band1Name '-'],['Intensity_' band1Name '+'],['Shift_' band1Name '+'],['FWHM_' band1Name '+'],['R^2_' band1Name],['Intensity_' band2Name],['Shift_' band2Name],['FWHM_' band2Name],['R^2_' band2Name],['Intensity_' band3Name],['Shift_' band3Name],['FWHM_' band3Name],['R^2_' band3Name]});
end

if lorentz==1 & nt==0
T=table(Int_1',center_1',FWHM_1',R2_1', Int_2',center_2',FWHM_2',R2_2', Int_3',center_3',FWHM_3',R2_3','VariableNames',{['Intensity_' band1Name],['Shift_' band1Name],['FWHM_' band1Name],['R^2_' band1Name],['Intensity_' band2Name],['Shift_' band2Name],['FWHM_' band2Name],['R^2_' band2Name],['Intensity_' band3Name],['Shift_' band3Name],['FWHM_' band3Name],['R^2_' band3Name]});
end
    
end

nameT=strcat(name(z),'_results','.csv');
namefinal=nameT{1};
writetable(T,namefinal);


%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Figures %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

%% General
% Raw spectra
if raw
%    figure(1), hold on;
%    subplot(div(1),div(2),z);
%    set(gca, 'ColorOrder', color, 'NextPlot', 'replacechildren');
%    plot(Shift,Intensity(:,:));
%    xlabel('Raman shift / cm^{-1}','FontSize',12);
%    ylabel('Intensity / a.u.','FontSize',12);
%    title([name{z} ' :Raw data']);
%    set(gca,'Box','on');
   
end

% Normalised spectra
if norm
figure(2), hold on;
subplot(div(1),div(2),z);
set(gca, 'ColorOrder', color, 'NextPlot', 'replacechildren');
plot(Shift,Intensity_norm(:,1:end));
xlabel('Raman shift / cm^{-1}','FontSize',12);
ylabel('Intensity / a.u.','FontSize',12);
title([name{z} ' :Normalised spectra']);
set(gca,'Box','on');
end

% Average spectra
figure(3), hold on;
plot(Shift,Intensity_av(),'color',color(round(end/2),:),'DisplayName',name{z});
xlabel('Raman shift / cm^{-1}','FontSize',12);
ylabel('Intensity / a.u.','FontSize',12);
title('Average spectra');
set(gca,'Box','on');
legend show

% Spectral range for Intensity ratio
if range
   figure(10),
   subplot(div(1),div(2),z);hold on
   set(gca, 'ColorOrder', color, 'NextPlot', 'replacechildren');
   plot(Shift(ind1Low:ind1High),Intensity_1range(:,:));hold on
   plot(Shift(ind2Low:ind2High),Intensity_2range(:,:));hold on
   plot(Shift(ind3Low:ind3High),Intensity_3range(:,:));
   xlabel('Raman shift / cm^{-1}','FontSize',12);
   ylabel('Intensity / a.u.','FontSize',12);
   title([name{z} ': Range used for Intensity ratio calculation'])
   set(gca,'Box','on');
end

%% Lorentz fit results: plots the spectra together with the best fit

if lorentz
    
figure(101);
    subplot(div(1),div(2),z); hold on
    set(gca, 'ColorOrder', color, 'NextPlot', 'replacechildren');
    plot(Shift_range1(:),Intensity_1range(:,:));hold on,
    plot(Shift_range1(:),I_1_fit(:,:),'r');
    xlabel('Raman shift / cm^{-1}','FontSize',12)
    ylabel('Intensity / a.u.','FontSize',12)
    set(gca,'Box','on');     
    title([name{z} ': Fitting results ' band1Name],'FontSize',12);
    

figure(201);
    subplot(div(1),div(2),z); hold on
    set(gca, 'ColorOrder', color, 'NextPlot', 'replacechildren');
    plot(Shift_range2(:),Intensity_2range(:,:));hold on,
    plot(Shift_range2(:),I_2_fit(:,:),'r');
    xlabel('Raman shift / cm^{-1}','FontSize',12)
    ylabel('Intensity / a.u.','FontSize',12)
    set(gca,'Box','on'); 
    title([name{z} ': Fitting results ' band2Name],'FontSize',12);
    
figure(301);
    subplot(div(1),div(2),z); hold on
    set(gca, 'ColorOrder', color, 'NextPlot', 'replacechildren');
    plot(Shift_range3(:),Intensity_3range(:,:));hold on,
    plot(Shift_range3(:),I_3_fit(:,:),'r');
    xlabel('Raman shift / cm^{-1}','FontSize',12)
    ylabel('Intensity / a.u.','FontSize',12)
    set(gca,'Box','on'); 
    title([name{z} ': Fitting results ' band3Name],'FontSize',12);

end

%% Histograms

% Histogram intensity ratio 21 and 31

figure(11), hold on;
g0=histogram(I21,'FaceColor',color(round(end/2),:),'DisplayName',[name{z} ': ' 'I_' band2Name '/I_' band1Name '=' num2str(I21_av) '\pm' num2str(I21_error)]);
g0.BinWidth = width_int; 
xlabel(['I(' band2Name ')/I(' band1Name ')'],'FontSize',12);
ylabel('counts','FontSize',12)
title(['Intensity ratio: I(' band2Name ')/I' band1Name ')'],'FontSize',12);
set(gca,'Box','on');
legend show

figure(13), hold on;
g01=histogram(I31,'FaceColor',color(round(end/2),:),'DisplayName',[name{z} ': ' 'I_' band3Name '/I_' band1Name '=' num2str(I31_av) '\pm' num2str(I31_error)]);
g01.BinWidth = width_int; 
xlabel(['I(' band3Name ')/I(' band1Name ')'],'FontSize',12);
ylabel('counts','FontSize',12)
title(['Intensity ratio: I(' band3Name ')/I(' band1Name ')'],'FontSize',12);
set(gca,'Box','on');
legend show


% Histograms peak 1

if lorentz

    if nt
figure(100), hold on;
g1=histogram(center_1plus,'FaceColor',color(round(end/3),:),'DisplayName',[name{z} ': shift ' band1Name '^+=' num2str(pos_1plus_av) '\pm' num2str(pos_1plus_err)]);
g1.BinWidth = width;  
g2=histogram(center_1min,'FaceColor',color(round(end/2),:),'DisplayName',[name{z} ': shift ' band1Name '^-=' num2str(pos_1min_av) '\pm' num2str(pos_1min_err)]);
g2.BinWidth = width; 
xlabel('Raman shift / cm^{-1}','FontSize',12);
ylabel('counts','FontSize',12);
title(['Raman shift ' band1Name ' modes'],'FontSize',12);
set(gca,'Box','on');
xlim([band1Low-10 band1High+10])
legend show

figure(99), hold on;
Iplus_minus=Int_1plus./Int_1min;
%binsInt=round((max(Iplus_minus)-min(Iplus_minus))/0.1);
g3=histogram(Iplus_minus,'FaceColor',color(round(end/2),:),'DisplayName',name{z}) ;
g3.BinWidth = width_int;
xlabel(['I_{' band1Name '^+}/I_{' band1Name '^-}'],'FontSize',12);
ylabel('counts','FontSize',12);
title(['Intensity ratio ' band1Name ' modes'],'FontSize',12);
set(gca,'Box','on');
legend show

figure(98), hold on;
g4=histogram(FWHM_1min,'FaceColor',color(round(end/2),:),'DisplayName',[name{z} ' FWHM ' band1Name '^-=' num2str(FWHM_1min_av) '\pm' num2str(FWHM_1min_err)]);
g4.BinWidth = width_fw;
g5=histogram(FWHM_1plus,'FaceColor',color(round(end/3),:),'DisplayName',[name{z} ' FWHM ' band1Name '^+=' num2str(FWHM_1plus_av) '\pm' num2str(FWHM_1plus_err)]);
g5.BinWidth = width_fw;
xlabel('FWHM / cm^{-1}','FontSize',12);
ylabel('counts','FontSize',12);
title(['FWHM ' band1Name ' band'],'FontSize',12);
set(gca,'Box','on');
legend show

else
figure(100), hold on
h1=histogram(center_1,'FaceColor',color(round(end/2),:),'DisplayName',[name{z} '; Shift ' band1Name '=' num2str(pos_1_av) '\pm' num2str(pos_1_err)]);
h1.BinWidth = width;
xlabel('Raman shift / cm^{-1}','FontSize',12);
ylabel('counts','FontSize',12);
title(['Raman shift ' band1Name ' band'],'FontSize',12);
legend show;
set(gca,'Box','on');

figure(98), hold on;
g4=histogram(FWHM_1,'FaceColor',color(round(end/2),:),'DisplayName',[name{z} ' FWHM ' band1Name '=' num2str(FWHM_1_av) '\pm' num2str(FWHM_1_err)]);
g4.BinWidth = width_fw;
xlabel('FWHM / cm^{-1}','FontSize',12);
ylabel('counts','FontSize',12);
title(['FWHM ' band1Name ' band'],'FontSize',12);
set(gca,'Box','on');
legend show
end

end

% Spectral features peak 2
binsShift=round(max(center_2)-min(center_2));
figure(200), hold on;
h2=histogram(center_2,'FaceColor',color(round(end/2),:),'DisplayName',[name{z} '; Shift ' band2Name '=' num2str(pos_2_av) '\pm' num2str(pos_2_err)]);
h2.BinWidth = width;
xlabel('Raman shift / cm^{-1}','FontSize',12);
ylabel('counts','FontSize',12);
title(['Raman shift ' band2Name ' band'],'FontSize',12);
set(gca,'Box','on');
legend show

if lorentz
figure(199), hold on;
h3=histogram(FWHM_2,'FaceColor',color(round(end/2),:),'DisplayName',[name{z} '; FWHM ' band2Name '=' num2str(FWHM_2_av) '\pm' num2str(FWHM_2_err)]);
h3.BinWidth = width_fw;
xlabel('FWHM / cm^{-1}','FontSize',12);
ylabel('counts','FontSize',12);
title(['FWHM ' band2Name ' band'],'FontSize',12);
set(gca,'Box','on');
legend show
end

% Spectral features peak 3
binsShift=round(max(center_3)-min(center_3));
figure(300), hold on;
h4=histogram(center_3,'FaceColor',color(round(end/2),:),'DisplayName',[name{z} '; Shift ' band2Name '=' num2str(pos_3_av) '\pm' num2str(pos_3_err)]);
h4.BinWidth = width;
xlabel('Raman shift / cm^{-1}','FontSize',12);
ylabel('counts','FontSize',12);
title(['Raman shift ' band3Name ' band'],'FontSize',12);
set(gca,'Box','on');
legend show

if lorentz
figure(299), hold on;
h5=histogram(FWHM_3,'FaceColor',color(round(end/2),:),'DisplayName',[name{z} '; FWHM ' band3Name '=' num2str(FWHM_3_av) '\pm' num2str(FWHM_3_err)]);
h5.BinWidth = width_fw;
xlabel('FWHM / cm^{-1}','FontSize',12);
ylabel('counts','FontSize',12);
title(['FWHM ' band3Name ' band'],'FontSize',12);
set(gca,'Box','on');
legend show
end

% Raman shift RBM
if rbm
binsShift=round((max(PeaksLoc)-min(PeaksLoc))/2);
figure(400), hold on
h4=histogram(PeaksLoc,'FaceColor',color(round(end/2),:),'DisplayName',name{z});
h4.BinWidth = width;
xlabel('Raman shift / cm^{-1}','FontSize',12);
ylabel('counts','FontSize',12);
title('Raman shift RBM modes','FontSize',12);
set(gca,'Box','on');
legend show
end

%% Correlations

% Position Band 1 vs I21
if correlation1

    if nt 

        if Int_1plus>Int_1min
        p1=fit(I21',center_1plus','poly1');
        p1_coeff = coeffvalues(p1);
        p1_confint = confint(p1);

        figure(1000), hold on;
        plot(I21(:),center_1plus,'o','color',color(round(end/2),:),'MarkerFaceColor',color(round(end/2),:),'DisplayName',name{z});
        plot(linspace(min(I21),max(I21)),polyval(p1_coeff,linspace(min(I21),max(I21))),'-','DisplayName',['Linear fit, m=' num2str(p1_coeff(1)) '\pm' num2str((p1_confint(2,1) - p1_confint(1,1))/2)],'color',color(round(end/2),:));
        xlabel(['I(' band2Name ')/I(' band1Name ')'],'FontSize',12);
        ylabel(['Raman shift ' band1Name '^+ band / cm^{-1}'],'FontSize',12);
        title([band1Name '^+ band vs I(' band2Name ')/I_' band1Name ')'],'FontSize',12);
        set(gca,'Box','on');
        legend show;  

        else
        p1=fit(I21',center_1min','poly1');
        p1_coeff = coeffvalues(p1);
        p1_confint = confint(p1);

        figure(1000), hold on;
        plot(I21(:),center_1min,'o','color',color(round(end/2),:),'MarkerFaceColor',color(round(end/2),:),'DisplayName',name{z});
        plot(linspace(min(I21),max(I21)),polyval(p1_coeff,linspace(min(I21),max(I21))),'-','DisplayName',['Linear fit, m=' num2str(p1_coeff(1)) '\pm' num2str((p1_confint(2,1) - p1_confint(1,1))/2)],'color',color(round(end/2),:));
        xlabel(['I(' band2Name ')/I(' band1Name ')'],'FontSize',12);
        ylabel(['Raman shift ' band1Name '^- band / cm^{-1}'],'FontSize',12);
        title([band1Name '^- band vs I(' band2Name ')/I(' band1Name ')'],'FontSize',12);
        set(gca,'Box','on');
        legend show;   
        end

    else
    
    p1=fit(I21',center_1','poly1');
    p1_coeff = coeffvalues(p1);
    p1_confint = confint(p1);

    figure(1000), hold on;
    plot(I21(:),center_1,'o','color',color(round(end/2),:),'MarkerFaceColor',color(round(end/2),:),'DisplayName',name{z});
    plot(linspace(min(I21),max(I21)),polyval(p1_coeff,linspace(min(I21),max(I21))),'-','DisplayName',['Linear fit, m=' num2str(p1_coeff(1)) '\pm' num2str((p1_confint(2,1) - p1_confint(1,1))/2)],'color',color(round(end/2),:));
    xlabel(['I(' band2Name ')/I(' band1Name ')'],'FontSize',12);
    ylabel(['Raman shift ' band1Name ' band / cm^{-1}'],'FontSize',12);
    title([band1Name ' band vs I(' band2Name ')/I(' band1Name ')'],'FontSize',12);
    set(gca,'Box','on');
    legend show;

    end

end


% Position Band 2 vs I21
if correlation2
p2=fit(I21',center_2','poly1');
p2_coeff = coeffvalues(p2);
p2_confint = confint(p2);

figure(1001), hold on;
plot(I21(:),center_2(:),'o','DisplayName',name{z},'color',color(round(end/2),:),'MarkerFaceColor',color(round(end/2),:))
plot(linspace(min(I21),max(I21)),polyval(p2_coeff,linspace(min(I21),max(I21))),'-','DisplayName',['Linear fit, m=' num2str(p2_coeff(1)) '\pm' num2str((p2_confint(2,1) - p2_confint(1,1))/2)],'color',color(round(end/2),:));
xlabel(['I_' band2Name '/I_' band1Name],'FontSize',12);
ylabel(['Raman shift ' band2Name ' band / cm^{-1}'],'FontSize',12);
title([band2Name ' band vs I(' band2Name ')/I(' band1Name ')'],'FontSize',12);
set(gca,'Box','on');
legend show;
end

% Position Band 3 vs I21
if correlation3
p3=fit(I21',center_3','poly1');
p3_coeff = coeffvalues(p3);
p3_confint = confint(p3);

figure(1002), hold on;
plot(I21(:),center_3(:),'o','DisplayName',name{z},'color',color(round(end/2),:),'MarkerFaceColor',color(round(end/2),:))
plot(linspace(min(I21),max(I21)),polyval(p3_coeff,linspace(min(I21),max(I21))),'-','DisplayName',['Linear fit, m=' num2str(p3_coeff(1)) '\pm' num2str((p3_confint(2,1) - p3_confint(1,1))/2)],'color',color(round(end/2),:));
xlabel(['I(' band2Name ')/I(' band1Name ')'],'FontSize',12);
ylabel(['Raman shift ' band3Name ' band / cm^{-1}'],'FontSize',12);
title([band3Name ' band vs I(' band2Name ')/I(' band1Name ')'],'FontSize',12);
set(gca,'Box','on');
legend show;
end


% Position Band 3 vs position band 1
if correlation4
    if nt
        if Int_1plus>Int_1min    
        p4=fit(center_1plus',center_3','poly1');
        p4_coeff = coeffvalues(p4);
        p4_confint = confint(p4);

        figure(1003), hold on;
        plot(center_1plus(:),center_3(:),'o','DisplayName',name{z},'color',color(round(end/2),:),'MarkerFaceColor',color(round(end/2),:))
        plot(linspace(min(center_1plus),max(center_1plus)),polyval(p4_coeff,linspace(min(center_1plus),max(center_1plus))),'-','DisplayName',['Linear fit, m=' num2str(p4_coeff(1)) '\pm' num2str((p4_confint(2,1) - p4_confint(1,1))/2)],'color',color(round(end/2),:));
        xlabel(['Raman shift ' band1Name '^+ band / cm^{-1}'],'FontSize',12);
        ylabel(['Raman shift ' band3Name ' band / cm^{-1}'],'FontSize',12);
        title(['Raman shift ' band3Name ' vs ' band1Name '^+ band'],'FontSize',12);
        set(gca,'Box','on');
        legend show;
        else
        p4=fit(center_1min',center_3','poly1');
        p4_coeff = coeffvalues(p4);
        p4_confint = confint(p4);

        figure(1003), hold on;
        plot(center_1min(:),center_3(:),'o','DisplayName',name{z},'color',color(round(end/2),:),'MarkerFaceColor',color(round(end/2),:))
        plot(linspace(min(center_1min),max(center_1min)),polyval(p4_coeff,linspace(min(center_1min),max(center_1min))),'-','DisplayName',['Linear fit, m=' num2str(p4_coeff(1)) '\pm' num2str((p4_confint(2,1) - p4_confint(1,1))/2)],'color',color(round(end/2),:));
        xlabel(['Raman shift ' band1Name '^- band / cm^{-1}'],'FontSize',12);
        ylabel(['Raman shift ' band3Name ' band / cm^{-1}'],'FontSize',12);
        title(['Raman shift ' band3Name ' vs ' band1Name '^- band'],'FontSize',12);
        set(gca,'Box','on');
        legend show;
        end

    else
    p4=fit(center_1',center_3','poly1');
    p4_coeff = coeffvalues(p4);
    p4_confint = confint(p4);

    figure(1003), hold on;
    plot(center_1(:),center_3(:),'o','DisplayName',name{z},'color',color(round(end/2),:),'MarkerFaceColor',color(round(end/2),:))
    plot(linspace(min(center_1),max(center_1)),polyval(p4_coeff,linspace(min(center_1),max(center_1))),'-','DisplayName',['Linear fit, m=' num2str(p4_coeff(1)) '\pm' num2str((p4_confint(2,1) - p4_confint(1,1))/2)],'color',color(round(end/2),:));
    xlabel(['Raman shift ' band1Name ' band / cm^{-1}'],'FontSize',12);
    ylabel(['Raman shift ' band3Name ' band / cm^{-1}'],'FontSize',12);
    title(['Raman shift ' band3Name ' vs ' band1Name ' band'],'FontSize',12);
    set(gca,'Box','on');
    legend show;
    end
end

% Position vs FWHM

if correlation5
    
    if lorentz
    % Position band 2 vs FWHM band 2
    figure(1005), hold on;
    plot(center_2(:),FWHM_2(:),'o','DisplayName',name{z},'color',color(round(end/2),:),'MarkerFaceColor',color(round(end/2),:))
    xlabel(['Raman shift ' band2Name ' band / cm^{-1}'],'FontSize',12);
    ylabel(['FWHM ' band2Name ' band / cm^{-1}'],'FontSize',12);
    title(['FWHM vs Shift ' band2Name ' band'],'FontSize',12);
    set(gca,'Box','on');
    legend show;

    % Position band 3 vs FWHM band 3
    figure(1006), hold on;
    plot(center_3(:),FWHM_3(:),'o','DisplayName',name{z},'color',color(round(end/2),:),'MarkerFaceColor',color(round(end/2),:))
    xlabel(['Raman shift ' band3Name ' band / cm^{-1}'],'FontSize',12);
    ylabel(['FWHM ' band3Name ' band / cm^{-1}'],'FontSize',12);
    title(['FWHM vs Shift ' band3Name ' band'],'FontSize',12);
    set(gca,'Box','on');
    legend show;

    % Position band 1 vs FWHM band 1
        if nt
            if Int_1plus>Int_1min    
            figure(1004), hold on;
            plot(center_1plus(:),FWHM_1plus(:),'o','DisplayName',name{z},'color',color(round(end/2),:),'MarkerFaceColor',color(round(end/2),:))
            xlabel(['Raman shift ' band1Name '^+ band / cm^{-1}'],'FontSize',12);
            ylabel(['FWHM ' band1Name '^+ band / cm^{-1}'],'FontSize',12);
            title(['FWHM vs Shift ' band1Name '^+ band'],'FontSize',12);
            set(gca,'Box','on');
            legend show;
            else
            figure(1004), hold on;
            plot(center_1min(:),FWHM_1min(:),'o','DisplayName',name{z},'color',color(round(end/2),:),'MarkerFaceColor',color(round(end/2),:))
            xlabel(['Raman shift ' band1Name '^- band / cm^{-1}'],'FontSize',12);
            ylabel(['FWHM ' band1Name '^- band / cm^{-1}'],'FontSize',12);
            title(['FWHM vs Shift ' band1Name '^- band'],'FontSize',12);
            set(gca,'Box','on');
            legend show;    
            end
        else
        figure(1004), hold on;
        plot(center_1(:),FWHM_1(:),'o','DisplayName',name{z},'color',color(round(end/2),:),'MarkerFaceColor',color(round(end/2),:))
        xlabel(['Raman shift ' band1Name ' band / cm^{-1}'],'FontSize',12);
        ylabel(['FWHM ' band1Name ' band / cm^{-1}'],'FontSize',12);
        title(['FWHM vs Shift ' band1Name ' band'],'FontSize',12);
        set(gca,'Box','on');
        legend show;
        end
    end
end

%% MAPS 

if map
raws_real=linspace(0,rawsmum,raws);
columns_real=linspace(0,colmum,col);

Imap=zeros(raws,col);
pos1map=zeros(raws,col);
pos2map=zeros(raws,col);
pos3map=zeros(raws,col);
I_2d_map=zeros(raws,col);
FWHM1map=zeros(raws,col);
FWHM2map=zeros(raws,col);
FWHM3map=zeros(raws,col);

for i=1:raws
   for j=1:col
      Imap(i,j)=I21(col*(i-1)+j);%I21
      I_2d_map(i,j)=I31(col*(i-1)+j);%I31
      pos2map(i,j)=center_2(col*(i-1)+j);
      pos3map(i,j)=center_3(col*(i-1)+j);
            
      if lorentz
       FWHM2map(i,j)=FWHM_2(col*(i-1)+j);   
       FWHM3map(i,j)=FWHM_3(col*(i-1)+j);
      end
      
      if nt
        if Int_1plus>Int_1min        
        pos1map(i,j)=center_1plus(col*(i-1)+j); 
        FWHM1map(i,j)=FWHM_1plus(col*(i-1)+j);
        else
        pos1map(i,j)=center_1min(col*(i-1)+j);
        FWHM1map(i,j)=FWHM_1min(col*(i-1)+j);
        end
      else
      pos1map(i,j)=center_1(col*(i-1)+j);   
          if lorentz
           FWHM1map(i,j)=FWHM_1(col*(i-1)+j);   
          end
      end
   end
end

% Intensity ratio 1/2
figure(40), hold on
subplot(div(1),div(2),z); hold on
levels=10;
[M,c]=contourf(columns_real,raws_real,Imap,levels);
set(c,'LineColor','none')
xlabel('x / \mum','FontSize',12);
ylabel('y / \mum','FontSize',12);
%caxis([0,1]); % Set custom defined color range
h = colorbar;
title([name(z) ': I(' band2Name ')/I(' band1Name ')']);

% Position band 1
figure(41), hold on
subplot(div(1),div(2),z); hold on
levels=10;
[M1,c1]=contourf(columns_real,raws_real,pos1map,levels);
set(c1,'LineColor','none')
xlabel('x / \mum','FontSize',12);
ylabel('y / \mum','FontSize',12);
%caxis([1583,1602]); % Set custom defined color range
colorbar
if nt
    if Int_1plus>Int_1min  
    title([name(z) ': ' band1Name '+ band position']);
    else
    title([name(z) ': ' band1Name '- band position']);
    end   
else
title([name(z) ': ' band1Name ' band position']);
end

% Position band 2
figure(42), hold on
subplot(div(1),div(2),z); hold on
levels=10;
[M2,c2]=contourf(columns_real,raws_real,pos2map,levels);
set(c2,'LineColor','none')
xlabel('x / \mum','FontSize',12);
ylabel('y / \mum','FontSize',12);
%caxis([1339,1356]);
colorbar
title([name(z) ': ' band2Name 'band position']);

% Postion band 3
figure(43), hold on
subplot(div(1),div(2),z); hold on
levels=10;
[M3,c3]=contourf(columns_real,raws_real,pos3map,levels);
set(c3,'LineColor','none')
xlabel('x / \mum','FontSize',12);
ylabel('y / \mum','FontSize',12);
%caxis([2670,2700]); % Set custom defined color range
colorbar
set(gca,'ColorScale','log')
title([name(z) ': ' band3Name 'band position']);

% Intensity ratio 1/3
if map4
figure(44), hold on
subplot(div(1),div(2),z); hold on
levels=10;
[M,c]=contourf(columns_real,raws_real,I_2d_map,levels);
set(c,'LineColor','none')
xlabel('x / \mum','FontSize',12);
ylabel('y / \mum','FontSize',12);
%caxis([0,1]); % Set custom defined color range
h = colorbar;
title([name(z) ': I(' band3Name ')/I(' band1Name ')']);
end

if lorentz
% FWHM band 1
if map1
figure(45), hold on
subplot(div(1),div(2),z); hold on
levels=10;
[M1,c1]=contourf(columns_real,raws_real,FWHM1map,levels);
set(c1,'LineColor','none')
xlabel('x / \mum','FontSize',12);
ylabel('y / \mum','FontSize',12);
%caxis([1583,1602]); % Set custom defined color range
colorbar
if nt
    if Int_1plus>Int_1min  
    title([name(z) ': FWHM ' band1Name '+']);
    else
    title([name(z) ': FWHM ' band1Name '-']);
    end   
else
title([name(z) ': FWHM ' band1Name]);
end
end

% FWHM band 2
if map2
figure(46), hold on
subplot(div(1),div(2),z); hold on
levels=10;
[M2,c2]=contourf(columns_real,raws_real,FWHM2map,levels);
set(c2,'LineColor','none')
xlabel('x / \mum','FontSize',12);
ylabel('y / \mum','FontSize',12);
%caxis([1339,1356]);
colorbar
title([name(z) ': FWHM ' band2Name]);
end

% FWHM band 3
if map3
figure(47), hold on
subplot(div(1),div(2),z); hold on
levels=10;
[M3,c3]=contourf(columns_real,raws_real,FWHM3map,levels);
set(c3,'LineColor','none')
xlabel('x / \mum','FontSize',12);
ylabel('y / \mum','FontSize',12);
%caxis([2670,2700]); % Set custom defined color range
colorbar
set(gca,'ColorScale','log')
title([name(z) ': FWHM ' band3Name]);
end
end
end

clearvars Shift Intensity IN

end


%%Save figures
if saveFigs
savefig(figure(1),[nameFigs '_RawSpectra.fig']); 
savefig(figure(2),[nameFigs '_NormalizedSpectra.fig']);    
savefig(figure(3),[nameFigs '_AverageSpectra.fig']);
%Fittings
savefig(figure(101),[nameFigs '_Fitting_' band1Name '.fig']);
savefig(figure(201),[nameFigs '_Fitting_' band2Name '.fig']);
savefig(figure(301),[nameFigs '_Fitting_' band3Name '.fig']);
%Histograms
savefig(figure(11),[nameFigs '_I' band2Name band1Name '.fig']);
savefig(figure(13),[nameFigs '_I' band3Name band1Name '.fig']);
savefig(figure(100),[nameFigs '_Shift ' band1Name '.fig']);
savefig(figure(99),[nameFigs '_Int Ratio' band1Name '.fig']);
savefig(figure(98),[nameFigs '_FWHM ' band1Name '.fig']);
savefig(figure(199),[nameFigs '_FWHM ' band2Name '.fig']);
savefig(figure(200),[nameFigs '_Shift ' band2Name '.fig']);
savefig(figure(299),[nameFigs '_FWHM ' band3Name '.fig']);
savefig(figure(300),[nameFigs '_Shift ' band3Name '.fig']);
% Correlations
savefig(figure(1000),[nameFigs '_I' band2Name band1Name 'vs shift_' band1Name '.fig']);
savefig(figure(1001),[nameFigs '_I' band2Name band1Name 'vs shift_' band2Name '.fig']);
savefig(figure(1002),[nameFigs '_I' band2Name band1Name 'vs shift_' band3Name '.fig']);
savefig(figure(1003),[nameFigs '_Shift' band1Name 'vs shift' band3Name '.fig']);
savefig(figure(1004),[nameFigs '_Center vs FWHM ' band1Name '.fig']);
savefig(figure(1005),[nameFigs '_Center vs FWHM ' band2Name '.fig']);
savefig(figure(1006),[nameFigs '_Center vs FWHM ' band3Name '.fig']);
%Maps
savefig(figure(40),[nameFigs '_Map_I' band2Name band1Name '.fig']);
savefig(figure(41),[nameFigs '_Map_Shift ' band1Name '.fig']);
savefig(figure(42),[nameFigs '_Map_Shift ' band2Name '.fig']);
savefig(figure(43),[nameFigs '_Map_Shift ' band3Name '.fig']);
savefig(figure(44),[nameFigs '_Map_I' band3Name band1Name '.fig']);
savefig(figure(45),[nameFigs '_Map_FWHM ' band1Name '.fig']);
savefig(figure(46),[nameFigs '_Map_FWHM ' band2Name '.fig']);
savefig(figure(47),[nameFigs '_Map_FWHM ' band3Name '.fig']);

end
