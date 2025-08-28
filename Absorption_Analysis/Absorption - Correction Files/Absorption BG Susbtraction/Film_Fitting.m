%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% Function to fit single walled carbon nanotube film absorption spectra
%
% Written by Moritz Pfohl, last modified on March 19th 2017.
%
% If you use this code to analyze your spectra, please cite:
%
% Pfohl, Tune, Graf, Zaumseil, Krupke, Flavel, "Fitting Single Walled
% Carbon Nanotube Optical Spectra"
% ACS OMEGA
% DOI: 10.1021/acsomega.6b00468
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function Film_Fitting
global center FWHM phonon_pos phonon_pos2 diam_n FWHM_g h c FWHM_g_22 center_22 w_solution s_ud FWHM_22 fit_region


%--------------------------------------------------------------------------
%
% List of adjustable variables:
%
% PEAK fit
%
% min_shift = Lower limit of peak shift compared to solution fit.
%
% max_shift = Upper limit of peak shift compared to solution fit.
%
% start_shift = Initial guess of peak shift comapred to solution fit.
%
% conc_change = Plus/Minus change in relative concentration compared to
%               solution fit (in %).
%
% broad_start = Initial guess of FWHM broadening.
%
% broad_min = Lower limit of FWHM broadening.
%
% broad_max = Upper limit of FWHM broadening.
%
% EXCITON PHONON SIDEBAND
%
% area_part_low = Lower limit of correction factor f1.
%
% area_part_up = Upper limit of correction factor f1.
%
% eps_broad_min = Lower limit of broadening of FWHM of exciton phonon
%                 sideband.
%
% eps_broad_start = Initial guess of broadening of FWHM of exciton phonon
%                   sideband.
%
% eps_broad_max = Upper limit of broadening of FWHM of exciton phonon
%                 sideband.
%
% GENERAL CONSTANTS
%
% h = Planck's constant in eV/s.
%
% c = Speed of light in m/s.
%--------------------------------------------------------------------------

% Peak fit

min_shift=0;
max_shift=40;
start_shift=30;
conc_change=100; % Plus/Minus change in relative concentration compared to solution fit (in %)
broad_start=2;
broad_min=0.5;
broad_max=2.5;
doping=1; % If doping=1: Area of S11 is decreased and area ratio of S11:S22 is allowed to change

% Exciton phonon sideband

area_part_low = -0.05;
area_part_up = 0.1;
eps_broad_min = 1;
eps_broad_start = 2;
eps_broad_max = 3;

% General constants

h=4.135667662e-15; % Planck's constant
c=299792458;       % Speed of light


%--------------------------------------------------------------------------
% Sub-function "get_data" located in Section 1 performs:
% - user selection of film absorption data
% - user defined background subtraction
% - user selection of wavelength range to be fitted (S11 or S22)
% - loads the results of the solution absorption fit
%--------------------------------------------------------------------------


[lambda,Film,center,center_22,f1,f1_22,FWHM,FWHM_22,w_solution,FWHM_g,FWHM_g_22,diam_n,legstr,legstr_22,phonon_pos,phonon_pos2,colors,w,area,area_22,colors_m,folder_name,parallel_fit] = get_data(h,c);


%--------------------------------------------------------------------------
% Sub-function "fit_film" located in Section 2 performs:
% - fit of absorption spectrum based on the line-shape chosen for the
%   solution absorption fit
% - calculation of chi-square value to evaluate the quality of the fit
%
%
% Sub-function "plot_data" located in Section 3 performs:
% - plot the results from "fit_film"
%
%
% Sub-function "check_fit" located in Section 4 performs:
% - user is asked to check whether they are satisfied with the quality of
%   the fit
%--------------------------------------------------------------------------


s_ud=0;
wl=[];

while 1>0
    
    [x_fit,sse,lambda_n,Film_n,w_solution,wl]=fit_film(lambda,Film+max(Film)*s_ud,FWHM,FWHM_22,FWHM_g,FWHM_g_22,center,center_22,f1,f1_22,w_solution,phonon_pos,phonon_pos2,diam_n,w,min_shift,max_shift,start_shift,h,c,conc_change,broad_start,broad_min,broad_max,area,area_22,area_part_low,area_part_up,eps_broad_min,eps_broad_start,eps_broad_max,wl,doping,parallel_fit);
    
    [area_f,y,legstrs]=plot_data(w,x_fit,FWHM,FWHM_22,lambda_n,Film_n,w_solution,phonon_pos,phonon_pos2,diam_n,FWHM_g,legstr,colors,sse,h,c,FWHM_g_22,center,center_22,legstr_22,fit_region,colors_m,doping);
    
    [tt,max_shift,start_shift,conc_change,wl]=check_fit(max_shift,start_shift,conc_change,wl);
    
    if tt=='y'
        break
    end
end


%--------------------------------------------------------------------------
% Sub-function "save_data" located in Section 5 performs:
% - save data to a cell "Results_Film" and save it as a MATLAB variable
%--------------------------------------------------------------------------


save_data(w,legstrs,x_fit,center,center_22,FWHM,FWHM_22,area_f,phonon_pos,phonon_pos2,FWHM_g,FWHM_g_22,diam_n,w_solution,h,c,fit_region,folder_name)


%--------------------------------------------------------------------------
% Sub-function "export_data" located in Section 6 performs:
% - save all individual peaks in one .txt file
%--------------------------------------------------------------------------


export_data(legstrs,lambda_n,Film_n,y,folder_name)


%--------------------------------------------------------------------------
%
%............................Sub-Functions.................................
%
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Section 1 - "get_data"
%--------------------------------------------------------------------------


function [lambda,Film,center,center_22,f1,f1_22,FWHM,FWHM_22,w_solution,FWHM_g,FWHM_g_22,diam_n,legstr,legstr_22,phonon_pos,phonon_pos2,colors,w,area,area_22,colors_m,folder_name,parallel_fit] = get_data(h,c)

% Prompt user to select their absorption data. An absorption measurement
% should have two columns (1 - wavelength, 2 - absorption data), separated
% by a comma and must be saved as a .txt file. The wavelength data can
% either be from high-low or low-high, because it will be corrected to
% display the data from low-high wavelength values.

vers=version('-release'); % Matlab version
rel=str2double(vers(1:4)); % Release year

if rel>=2015
    parallel_fit=input('Do you want to fit your data on parallel cores (y/n)?: ','s');
    parallel_fit=sscanf(parallel_fit,'%s');
    
    v=ver;
    parallel_box=any(strcmp('Parallel Computing Toolbox', {v.Name}));
    
    if parallel_box==1 && parallel_fit=='y'
        if isempty(gcp('nocreate'))==1
            nc=feature('numCores');
            parpool(nc-1);
        end
        vers=version('-release'); % Matlab version
        rel=str2double(vers(1:4)); % Release year
    end
else
    parallel_fit='n';
end

disp('Select the file containing the absorption data');
filename = uigetfile('*.txt','Select the file containing the absorption data');

string='%f32 %f32';

fid=fopen(filename);
D = textscan(fid,string,'Delimiter',',','HeaderLines',2);
fclose(fid);

Data(:,1)=D{1};
Data(:,2)=D{2};

if Data(2,1)<Data(1,1)
    Data=flipud(Data);
end

tt=figure; hold on

plot(Data(:,1),Data(:,2),'k','linewidth',2)
xlabel('Wavelength (nm)','FontSize',33)
ylabel('Absorption (a.u.)','FontSize',33)
title('Measurement Data','FontSize',36)
set(gca,'LineWidth',2)
set(gca,'FontSize',30)

% Prompt the user to select a method of background subtraction. Method (1)
% is adopted from Tian et al. (RSC Adv. 2015, 5, 102974). Method (2)
% is based on the work of Naumov et al. (ACS Nano, 2011, 5, 3, 1639–1648).
% Method (3) is based on the work of Nair et al. (Anal. Chem., 2006, 78,
% 7689-7696). Method (4) performs no background subtraction. For method 1,
% 2 and 3, the background subtraction is performed and constrained by
% functions defined in Section 1.1.
%
% Definition of Fano function:
% Fano = x(7)*(x(1)+(E-x(2))/(x(3)/2)).^2./(1+((E-x(2))/(x(3)/2)).^2)
% with: x(7)  =cross section far away from the resonance
%       x(1)^2=ratio of the strength of the excitonic transition to free pi
%              -> pi* transition
%       x(2)  = peak position
%       x(3)  = FWHM of the peak
%
% Definiton of Lorentzian:
% Lorentzian = x(4)./(1+((E-x(5))/(0.5*x(6))).^2)
% with: x(4)=height of the Lorentzian peak
%       x(5)=Lorentzian peak position
%       x(6)=FWHM of the Lorentzian peak

while 1>0
    q=input('Enter the starting and ending wavelengths for background subtraction, e.g. 400 1400. \nBe aware that the S11 and S22 transitions to be fitted \nshould be similiar to the one defined during the solution absorption fit: ','s');
    q=sscanf(q,'%i');
    ck=0;
    
    if q(1)<Data(1,1)
        warning('Starting wavelength is outside the measured range. Enter a larger starting wavelength.')
        ck=1;
    end
    if q(2)>Data(end,1)
        warning('Ending wavelength is outside the measured range. Enter a smaller ending wavelength.')
        ck=1;
    end
    if ck==0
        close(tt)
        break
    end
end

rr=1;
bkg=1;

% Fano + Lorentzian Background
advanced=0;

while 1>0
    
    if bkg==1
        w=input('Please enter number of background subtraction:\n (1) - Fano + Lorentzian\n (2) - A*exp(-b*lambda)\n (3) - k/lambda^b\n (4) - No background\n ','s');
        w=sscanf(w,'%i');
        bkg=2;
    end
    
    switch w
        case 1 % Based on Tian et al., RSC Adv. 2015, 5, 102974
            
            if rr==2
                while 1>0
                    q=input('Enter the starting and ending wavelengths for the region of interest, e.g. 400 1400: ','s');
                    q=sscanf(q,'%i');
                    ck=0;
                    
                    if q(1)<Data(1,1)
                        warning('Starting wavelength is outside the measured range. Enter a larger starting wavelength.')
                        ck=1;
                    end
                    if q(2)>Data(end,1)
                        warning('Ending wavelength is outside the measured range. Enter a smaller ending wavelength.')
                        ck=1;
                    end
                    if ck==0
                        break
                    end
                end
            end
            
            xx=round(Data(1,1)):1:round(Data(end,1));
            yy=interp1(Data(:,1),Data(:,2),xx);
            xx=xx'; yy=yy';
            
            start=find(xx==q(1));
            End=find(xx==q(2));
            
            if advanced==1
                cnames = {'Lower Limit','Initial Guess','Upper Limit'};
                rnames = {'Ratio of strength of excitonic to pi* transition','Center of Fano peak (eV)','FWHM of Fano peak (eV)', 'Height of Fano peak','Center of Lorentzian peak (eV)','FWHM of Lorentzian peak (eV)', 'Height of Lorentzian peak'};
                vals=[true true true true true true true true true true true true true true true true true true true true true];
                
                while 1>0
                    f=figure('Position',[100 500 800 150]); t = uitable('Position',[0 0 800 150],'Data', data_FL,'ColumnName',cnames,'RowName',rnames,'ColumnEditable', vals);
                    vv=version('-release');
                    year=str2double(vv(1:4));
                    if year>=2014
                        t.Position(3) = t.Extent(3);
                        t.Position(4) = t.Extent(4);
                    end
                    opt=input('Please change the variables in the table and press d if done: ','s');
                    if opt=='d'
                        data_FL=t.Data;
                        close(f)
                        break
                    end
                end
                advanced=0;
            else
                if exist('data_FL')==0 % if there is no variable called data_FL, create it
                    lb=[-4,4.4,0.4,-0.03,4.8,1.5,0];
                    ub=[-2,4.7,0.8,max(yy(1:200))/9,6,3.5,max(yy(1:200))]; % define upper and lower boundary conditions
                    x_start=[-3,4.4,0.6,0.1,4.8,1,0.02];                  
                    
                    data_FL(:,1)=lb;
                    data_FL(:,2)=x_start;
                    data_FL(:,3)=ub;
                end
            end
            
            E=h*c./(xx(start:End,1)*1e-9); % transform nm in eV
            
            dlmwrite('test.txt',[E yy(start:End,1)]); % save data in test.txt so that the fit functions can access it easily
            
            options = optimoptions('fmincon','Algorithm','interior-point','Display','off');
            options_l = optimoptions(@lsqnonlin,'Display','off');
            
            F_fano=@(x) double(yy(start:End,1)-x(7)*(x(1)+(E-x(2))/(x(3)/2)).^2./(1+((E-x(2))/(x(3)/2)).^2)-x(4)./(1+((E-x(5))/(0.5*x(6))).^2)); % Residual function
            x=lsqnonlin(F_fano,double(data_FL(:,2)'),[],[],options_l); % unrestricted fit with initial guess of starting values
            
            [x_bkg,~] = fmincon(@Fano,x,[],[],[],[],double(data_FL(:,1)'),double(data_FL(:,3)'),@conf_fano,options); % constrained background subtraction with lower and upper boundary conditions
            
            BKG=yy(start:End,1)-F_fano(x_bkg);
            
            delete('test.txt');
            
            tt=figure; hold on
            plot(xx(start:End,1),yy(start:End,1),'k','linewidth',2)
            plot(xx(start:End,1),BKG,'r','linewidth',2)
            xlabel('Wavelength (nm)','FontSize',33)
            ylabel('Absorption (a.u.)','FontSize',33)
            legend('Measured Data','Background')
            set(gca,'LineWidth',2)
            set(gca,'FontSize',30)
            xlim(q)
            ylim([min(BKG)*0.9 max(BKG)*1.1])
            
            opt=input('Are you satisfied with the range of your background subtraction (y/n)? ','s');
            
            if opt=='y'
                close(tt)
                lambda=Data(start:End,1);
                Film=smooth(F_fano(x_bkg)); % Smooth data: MATLAB function "smooth" is used to avoid the incorrect assignment of peaks to nanotubes
                break
            else
                ask2=input('Would you like to change the wavelength range (1), the type of background subtraction (2) or\nmodify the parameters of the Fano+Lorentian fit (3) ?: ','s');
                ask2=sscanf(ask2,'%i');
                switch ask2
                    case 1
                        rr=2;
                        bkg=2;
                    case 2
                        rr=1;
                        bkg=1;
                    case 3
                        rr=1;
                        bkg=2;
                        advanced=1;
                end
            end
            
        case 2
            
            if rr==2
                while 1>0
                    q=input('Enter the starting and ending wavelengths for the region of interest, e.g. 400 1400: ','s');
                    q=sscanf(q,'%i');
                    ck=0;
                    
                    if q(1)<Data(1,1)
                        warning('Starting wavelength is outside the measured range. Enter a larger starting wavelength.')
                        ck=1;
                    end
                    if q(2)>Data(end,1)
                        warning('Ending wavelength is outside the measured range. Enter a smaller ending wavelength.')
                        ck=1;
                    end
                    if ck==0
                        break
                    end
                end
            end
            
            xx=round(Data(1,1)):1:round(Data(end,1));
            yy=interp1(Data(:,1),Data(:,2),xx);
            xx=xx'; yy=yy';
            
            start=find(xx==q(1));
            End=find(xx==q(2));
            
            dlmwrite('test.txt',[xx(start:End,1) yy(start:End,1)]);
            
            options = optimoptions('fmincon','Algorithm','interior-point','Display','off');
            
            xn = fmincon(@Naumov,[.2 0.00002],[],[],[],[],[],[],@conf_naumov,options); % % Approach of Naumov et al. requires two starting values for A (0.2) and b (0.002).
            
            F_naumov=@(x) double(yy(start:End,1)-x(1)*exp(-x(2).*xx(start:End,1))); % Background subtracted data!
            
            BKG=yy(start:End,1)-F_naumov(xn);
            
            delete('test.txt');
            
            tt=figure; hold on
            
            plot(xx(start:End,1),yy(start:End,1),'k','linewidth',2)
            plot(xx(start:End,1),BKG,'r','linewidth',2)
            xlabel('Wavelength (nm)','FontSize',33)
            ylabel('Absorption (a.u.)','FontSize',33)
            legend('Measured Data','Background')
            set(gca,'LineWidth',2)
            set(gca,'FontSize',30)
            xlim(q)
            ylim([min(BKG)*0.9 max(BKG)*1.1])
            
            opt=input('Are you satisfied with the range of your background subtraction (y/n)? ','s');
            if opt=='y'
                close(tt)
                lambda=xx(start:End,1);
                Film=smooth(F_naumov(xn)); % Smooth data: MATLAB function "smooth" is used to avoid the incorrect assignment of peaks to nanotubes
                break
            else
                ask2=input('Would you like to change the wavelength range (1) or the type of background subtraction (2) ?: ','s');
                ask2=sscanf(ask2,'%i');
                if ask2==1
                    rr=2;
                    bkg=2;
                else
                    rr=1;
                    bkg=1;
                end
                close(tt)
            end
            
        case 3
            
            if rr==2
                while 1>0
                    q=input('Enter the starting and ending wavelengths for the region of interest, e.g. 400 1400: ','s');
                    q=sscanf(q,'%i');
                    ck=0;
                    
                    if q(1)<Data(1,1)
                        warning('Starting wavelength is outside the measured range. Enter a larger starting wavelength.')
                        ck=1;
                    end
                    if q(2)>Data(end,1)
                        warning('Ending wavelength is outside the measured range. Enter a smaller ending wavelength.')
                        ck=1;
                    end
                    if ck==0
                        break
                    end
                end
            end
            
            xx=round(Data(1,1)):1:round(Data(end,1));
            yy=interp1(Data(:,1),Data(:,2),xx);
            xx=xx'; yy=yy';
            
            start=find(xx==q(1));
            End=find(xx==q(2));
            
            dlmwrite('test.txt',[xx(start:End,1) yy(start:End,1)]);
            
            options = optimoptions('fmincon','Algorithm','interior-point','Display','off');
            
            xn = fmincon(@Nair,[.2 0.002],[],[],[],[],[],[],@conf_nair,options); % Approach of Nair et al. requires two starting values for k (0.2) and b (0.002).
            
            F_nair=@(x) double(yy(start:End,1)-x(1)./xx(start:End,1).^x(2)); % Background subtracted data!
            
            BKG=yy(start:End,1)-F_nair(xn);
            
            delete('test.txt');
            
            tt=figure; hold on
            
            plot(xx(start:End,1),yy(start:End,1),'k','linewidth',2)
            plot(xx(start:End,1),BKG,'r','linewidth',2)
            xlabel('Wavelength (nm)','FontSize',33)
            ylabel('Absorption (a.u.)','FontSize',33)
            legend('Measured Data','Background')
            set(gca,'LineWidth',2)
            set(gca,'FontSize',30)
            xlim(q)
            ylim([min(BKG)*0.9 max(BKG)*1.1])
            
            opt=input('Are you satisfied with the range of your background subtraction (y/n)? ','s');
            if opt=='y'
                close(tt)
                lambda=Data(start:End,1);
                Film=smooth(F_nair(xn)); % Smooth data: MATLAB function "smooth" is used to avoid the incorrect assignment of peaks to nanotubes
                break
            else
                ask2=input('Would you like to change the wavelength range (1) or the type of background subtraction (2) ?: ','s');
                ask2=sscanf(ask2,'%i');
                if ask2==1
                    rr=2;
                    bkg=2;
                else
                    rr=1;
                    bkg=1;
                end
                close(tt)
            end
            
        case 4 % No Background correction
            
            xx=round(Data(1,1)):1:round(Data(end,1));
            yy=interp1(Data(:,1),Data(:,2),xx);
            xx=xx'; yy=yy';
            
            start=find(xx==q(1));
            End=find(xx==q(2));
            
            lambda=xx(start:End,1);
            Film=smooth(yy(start:End,1)); % MATLAB function "smooth" is used to avoid the incorrect assignment of peaks to nanotubes
            
            break
    end
end

% Prompt the user to select the folder containing the results of the
% solution absorption fit.

disp('Please select the folder containing the variables for the film fit, that was automatically created during the solution fit.')
folder_name=uigetdir(pwd,'Please select the folder containing the variables for the film fit, that was automatically created during the solution fit.');

oldD=cd(folder_name);

D=dir;

Results_M11=[];

for i=1:size(D,1)
    a=strfind(D(i,1).name,'Results');
    b=strfind(D(i,1).name,'Result_S22andM11'); % Check if any metallic nanotubes were assigned to the S22 region
    if isempty(a)==0
        file=D(i,1).name;
        r=load(file);
        Results=r.Results; % All relevant results from the solution fit
    end
    if isempty(b)==0
        file=D(i,1).name;
        r=load(file);
        Results_M11=r.Results;
    end
end

ph=load('cc');
phonon_pos=ph.cc; % If no exciton phonon sidebands (EPS) were considered, cc is an empty matrix; otherwise it contains the number of the CNT peak(s) that have an EPS

ph2=load('cc2');
phonon_pos2=ph2.cc2;

ll=load('legstr');
legstr_sol=ll.legstr;

col=load('colors_n');
colors=col.colors_n;

profile=load('w');
w=profile.w;

colors_m=[];

cd(oldD);

% Load results of the solution absorption fit.

switch w
    case {1,2} % Lorentzian or Gaussian Fit
        if size(Results,2)==9 % Either S11 or S22
            nn=1;
        else % Entire region was fitted
            nn=2;
        end
        
        numb_CNT=size(Results,1)-length(phonon_pos)-length(phonon_pos2); % number of CNTs that were used to fit the solution data
        
        center=zeros(numb_CNT-1,1);
        FWHM=zeros(numb_CNT-1,1);
        area=zeros(numb_CNT-1,1);
        diam_n=zeros(numb_CNT-1,1);
        
        center_22=[];
        FWHM_22=[];
        area_22=[];
        legstr_22=[];
        
        for i=2:numb_CNT % S11 or S22
            center(i-1)=Results{i,3};
            FWHM(i-1)=Results{i,4};
            area(i-1)=Results{i,7};
            diam_n(i-1)=Results{i,2};
            legstr{i-1}=legstr_sol{i-1};
        end
        
        if isempty(Results_M11)==0 && nn==1 % Metallic nanotubes were assigned to S22 region
            for i=2:numb_CNT
                if isempty(Results_M11{i,3})==0
                    center(i-1)=Results_M11{i,3};
                    FWHM(i-1)=Results_M11{i,4};
                    area(i-1)=Results_M11{i,7};
                end
            end
            for i=numb_CNT+1:size(Results_M11,1)-length(phonon_pos)
                center=[center;Results_M11{i,3}];
                FWHM=[FWHM;Results_M11{i,4}];
                area=[area;Results_M11{i,7}];
                legstr=[legstr{:} Results_M11(i,1)];
                colors_m(i-numb_CNT,:)=[1,(Results_M11{i,3}-lambda(1))/(lambda(end)-lambda(1)),1-(Results_M11{i,3}-lambda(1))/(lambda(end)-lambda(1))];
                diam_n=[diam_n;Results_M11{i,2}];
            end
        end
        
        if nn==2 % S11 and S22
            
            n=0;
            for i=2:numb_CNT
                if isempty(Results{i,8})==0
                    center_22(i-1)=Results{i,8};
                    FWHM_22(i-1)=Results{i,9};
                    area_22(i-1)=Results{i,12};
                    legstr_22{i-1-n}=legstr_sol{length(center)+i-n-1};
                else
                    n=n+1;
                end
            end
            
            if isempty(Results_M11)==0 % Metallic nanotubes were assigned to S22 region
                for i=2:numb_CNT
                    if isempty(Results_M11{i,3})==0
                        center_22(i-1)=Results_M11{i,3};
                        FWHM_22(i-1)=Results_M11{i,4};
                        area_22(i-1)=Results_M11{i,7};
                    end
                end
                for i=numb_CNT+1:size(Results_M11,1)-length(phonon_pos2)
                    center_22=[center_22,Results_M11{i,3}];
                    FWHM_22=[FWHM_22,Results_M11{i,4}];
                    area_22=[area_22,Results_M11{i,7}];
                    legstr_22=[legstr_22{:} Results_M11(i,1)];
                    colors_m(i-numb_CNT,:)=[1,(Results_M11{i,3}-lambda(1))/(lambda(end)-lambda(1)),1-(Results_M11{i,3}-lambda(1))/(lambda(end)-lambda(1))];
                    diam_n=[diam_n;Results_M11{i,2}];
                end
            end
        end
        
        FWHM_g=[]; FWHM_g_22=[];
        f1=[]; f1_22=[];
        
        for i=1:length(phonon_pos)
            FWHM_g(i)=Results{numb_CNT+i,4};
            f1=Results{numb_CNT+i,6};
        end
        
        for i=1:length(phonon_pos2)
            if isempty(Results_M11)==0
                FWHM_g_22(i)=Results_M11{size(Results_M11,1)-length(phonon_pos2)+i,4};
                f1_22=Results_M11{size(Results_M11,1)-length(phonon_pos2)+i,6};
            else
                FWHM_g_22(i)=Results{numb_CNT+length(phonon_pos)+i,9};
                f1_22=Results{numb_CNT+length(phonon_pos)+i,11};
            end
        end
        
    case 3 % Voigt
        
        if size(Results,2)==12 % Either S11 or S22
            nn=1;
        else % Entire region was fitted
            nn=2;
        end
        
        numb_CNT=size(Results,1)-length(phonon_pos)-length(phonon_pos2); % number of CNTs that were used to fit the solution data
        
        center=zeros(numb_CNT-1,1);
        FWHM=zeros(numb_CNT-1,1);
        area=zeros(numb_CNT-1,1);
        diam_n=zeros(numb_CNT-1,1);
        
        center_22=[];
        FWHM_22=[];
        area_22=[];
        legstr_22=[];
        
        for i=2:numb_CNT % S11 or S22
            center(i-1)=Results{i,3};
            FWHM(i-1,1)=Results{i,4};
            FWHM(i-1,2)=Results{i,5};
            FWHM(i-1,3)=Results{i,6};
            area(i-1)=Results{i,10};
            diam_n(i-1)=Results{i,2};
            legstr{i-1}=legstr_sol{i-1};
        end
        
        if isempty(Results_M11)==0 && nn==1 % Metallic nanotubes were assigned to S22 region
            for i=2:numb_CNT
                if isempty(Results_M11{i,3})==0
                    center(i-1)=Results_M11{i,3};
                    FWHM(i-1,1)=Results_M11{i,4};
                    FWHM(i-1,2)=Results_M11{i,5};
                    FWHM(i-1,3)=Results_M11{i,6};
                    area(i-1)=Results_M11{i,10};
                end
            end
            for i=numb_CNT+1:size(Results_M11,1)-length(phonon_pos)
                center=[center;Results_M11{i,3}];
                FWHM=[FWHM;Results_M11{i,4} Results_M11{i,5} Results_M11{i,6}];
                area=[area;Results_M11{i,10}];
                legstr=[legstr{:} Results_M11(i,1)];
                colors_m(i-numb_CNT,:)=[1,(Results_M11{i,3}-lambda(1))/(lambda(end)-lambda(1)),1-(Results_M11{i,3}-lambda(1))/(lambda(end)-lambda(1))];
                diam_n=[diam_n;Results_M11{i,2}];
            end
        end
        
        if nn==2 % S11 and S22
            n=0;
            for i=2:numb_CNT
                if isempty(Results{i,11})==0
                    center_22(i-1)=Results{i,11};
                    FWHM_22(i-1,1)=Results{i,12};
                    FWHM_22(i-1,2)=Results{i,13};
                    FWHM_22(i-1,3)=Results{i,14};
                    area_22(i-1)=Results{i,18};
                    legstr_22{i-1-n}=legstr_sol{length(center)+i-n-1};
                else
                    n=n+1;
                end
            end
            
            
            if isempty(Results_M11)==0 % Metallic nanotubes were assigned to S22 region
                for i=2:numb_CNT
                    if isempty(Results_M11{i,3})==0
                        center_22(i-1)=Results_M11{i,3};
                        FWHM_22(i-1,1)=Results_M11{i,4};
                        FWHM_22(i-1,2)=Results_M11{i,5};
                        FWHM_22(i-1,3)=Results_M11{i,6};
                        area_22(i-1)=Results_M11{i,10};
                    end
                end
                for i=numb_CNT+1:size(Results_M11,1)-length(phonon_pos2)
                    center_22=[center_22,Results_M11{i,3}];
                    FWHM_22=[FWHM_22;Results_M11{i,4} Results_M11{i,5} Results_M11{i,6}];
                    area_22=[area_22,Results_M11{i,10}];
                    legstr_22=[legstr_22{:} Results_M11(i,1)];
                    colors_m(i-numb_CNT,:)=[1,(Results_M11{i,3}-lambda(1))/(lambda(end)-lambda(1)),1-(Results_M11{i,3}-lambda(1))/(lambda(end)-lambda(1))];
                    diam_n=[diam_n;Results_M11{i,2}];
                end
            end
        end
        
        FWHM_g=[]; FWHM_g_22=[];
        f1=[]; f1_22=[];
        
        for i=1:length(phonon_pos)
            FWHM_g(i)=Results{numb_CNT+i,5};
            f1(i)=Results{numb_CNT+i,9};
        end
        
        for i=1:length(phonon_pos2)
            if isempty(Results_M11)==0
                FWHM_g_22(i)=Results_M11{size(Results_M11,1)-length(phonon_pos2)+i,5};
                f1_22=Results_M11{size(Results_M11,1)-length(phonon_pos2)+i,9};
            else
                FWHM_g_22(i)=Results{numb_CNT+length(phonon_pos)+i,13};
                f1_22(i)=Results{numb_CNT+length(phonon_pos)+i,17};
            end
        end
end

w_solution=area/sum(area);


%--------------------------------------------------------------------------
% Section 1.1 - Background Subtraction
%--------------------------------------------------------------------------

% Method (1) based on Tian et al.

function err = Fano(x)

A=dlmread('test.txt');
c = A(:,2)-x(7)*(x(1)+(A(:,1)-x(2))/(x(3)/2)).^2./(1+((A(:,1)-x(2))/(x(3)/2)).^2)-x(4)./(1+((A(:,1)-x(5))/(0.5*x(6))).^2);
err = double(c'*c);

function [c,ceq] = conf_fano(x) % Constraint function, that forces the background to be smaller than the absorption data for every single wavelength

A=dlmread('test.txt');
% Nonlinear inequality constraints
c = double(x(7)*(x(1)+(A(:,1)-x(2))/(x(3)/2)).^2./(1+((A(:,1)-x(2))/(x(3)/2)).^2)+x(4)./(1+((A(:,1)-x(5))/(0.5*x(6))).^2) - A(:,2));
% Nonlinear equality constraints
ceq = [];

% Method (2) based on Naumov et al.

function err = Naumov(x) % Calculates the difference between the absorption data and the background. The MATLAB function "fmincon" tries to minimize this difference by fitting x(1)=A and x(2)=b

A=dlmread('test.txt');
c = A(:,2)-x(1)*exp(-x(2).*A(:,1));
err = double(sum(c));

function [c,ceq] = conf_naumov(x) % Constraint function, that forces the background to be smaller than the absorption data for every single wavelength

A=dlmread('test.txt');
% Nonlinear inequality constraints
c = double(x(1)*exp(-x(2).*A(:,1))-A(:,2));
% Nonlinear equality constraints
ceq = [];

% Method (3) based on Nair et al.

function err = Nair(x) % Calculates the difference between the absorption data and the background. The MATLAB function "fmincon" tries to minimize this difference by fitting x(1)=k and x(2)=b

A=dlmread('test.txt');
c = A(:,2)-x(1)./A(:,1).^x(2);
err = double(sum(c));

function [c,ceq] = conf_nair(x) % Constraint function, that forces the background to be smaller than the absorption data for every single wavelength

A=dlmread('test.txt');
% Nonlinear inequality constraints
c = double(x(1)./A(:,1).^x(2)-A(:,2));
% Nonlinear equality constraints
ceq = [];


%..........................................................................


%--------------------------------------------------------------------------
% Section 2 - "fit_film"
%--------------------------------------------------------------------------


function [x_fit,sse,lambda_n,Film_n,w_solution,wl]=fit_film(lambda,Film,FWHM,FWHM_22,FWHM_g,FWHM_g_22,center,center_22,f1,f1_22,w_solution,phonon_pos,phonon_pos2,diam_n,w,min_shift,max_shift,start_shift,h,c,conc_change,broad_start,broad_min,broad_max,area,area_22,area_part_low,area_part_up,eps_broad_min,eps_broad_start,eps_broad_max,wl,doping,parallel_fit)
global fit_region

% In the event that the entire region of the solution absorption
% measurement was fitted, the user is asked to choose S11, S22 or the
% entire region to fit the film absorption measurement. If a Lorentzian or
% Gaussian line-profile was chosen for the solution fit, the subfunction
% "generate_string" located in Section 2.1 will be used to create a
% function handle to calculate the residuals of the fit. If a Voigtian
% line-profile was chosen, the subfunction "generate_start_values_voigt"
% located in Section 2.2 will be used to create the starting values and
% constraints for the Voigt fit performed in Section 2.3 with the
% subfunctions "Voigt" and "complexErrorFunction".

if isempty(area_22)==0 && isempty(wl)==1 % Entire region was fitted for solution absorption measurements
    fit_region=input('Would you like to fit S11 (1), S22 (2) or the entire region of your spectrum (3)?\n','s');
    fit_region=sscanf(fit_region,'%i');
end

if isempty(area_22)==1
    fit_region=[];
end

if isempty(wl)==1 % Select wavelength for the first evaluation of "fit_film"
    wl=input('Please enter the starting and ending wavelength of your fit regime.\nBe aware that the starting and ending wavelengths should be inside the wavelength regime defined during background subtraction:\n ','s');
    wl=sscanf(wl,'%i');
end

start=find(lambda==wl(1));
End=find(lambda==wl(2));

lambda_n=lambda(start:End)';
Film_n=Film(start:End)';

switch w
    case {1,2}
        if parallel_fit=='y'
            try
                options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','TolX',1e-20,'Display','off','UseParallel',true);
            catch
                options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','TolX',1e-20,'Display','off');
            end
        else
            options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','TolX',1e-20,'Display','off');
        end
        
        [plotString,x0,lb,ub,w_solution]=generate_string(lambda_n,Film_n,FWHM,FWHM_g,center,f1,w_solution,phonon_pos,diam_n,w,FWHM_22,FWHM_g_22,center_22,f1_22,phonon_pos2,min_shift,max_shift,start_shift,h,c,conc_change,broad_start,broad_min,broad_max,area,area_22,area_part_low,area_part_up,eps_broad_min,eps_broad_start,eps_broad_max,fit_region,doping);
        
        Fit_sum=str2func(plotString);
        
        [x_fit,~]=lsqnonlin(Fit_sum,x0,lb,ub,options);
        
        sse=sum(Fit_sum(x_fit).^2)/sum((Film_n-mean(Film_n)).^2);   
    case 3
        
        if parallel_fit=='y'
            try
                options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','FinDiffRelStep',1e-6,'Display','off','UseParallel',true);
            catch
                options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','FinDiffRelStep',1e-6,'Display','off');
            end
        else
            options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','FinDiffRelStep',1e-6,'Display','off');
        end
        
        [x0,lb,ub,w_solution]=generate_start_values_voigt(center,center_22,lambda_n,Film_n,f1,w_solution,phonon_pos,area,area_22,min_shift,start_shift,max_shift,conc_change,phonon_pos2,f1_22,area_part_low,area_part_up,eps_broad_min,eps_broad_start,eps_broad_max,broad_start,broad_min,broad_max,doping);
        
        x_fit=lsqnonlin(@(x) Voigt(x,lambda_n,Film_n,w_solution,center,center_22,FWHM,phonon_pos,phonon_pos2,diam_n,FWHM_g,h,c,FWHM_22,FWHM_g_22),double(x0),double(lb),double(ub),options);
        
        sse=sum(Voigt(x_fit,lambda_n,Film_n,w_solution,center,center_22,FWHM,phonon_pos,phonon_pos2,diam_n,FWHM_g,h,c,FWHM_22,FWHM_g_22).^2)/sum((Film_n-mean(Film_n)).^2);        
end


%--------------------------------------------------------------------------
% Section 2.1 - "generate_string"
%--------------------------------------------------------------------------


function [plotString,x0,lb,ub,w_sol]=generate_string(lambda,Film,fwhm,FWHM_g,Cent,f1,w_solution,phonon_pos,diam_n,w,FWHM_22,FWHM_g_22,Cent_22,f1_22,phonon_pos2,min_shift,max_shift,start_shift,h,c,conc_change,broad_start,broad_min,broad_max,area,area_22,area_part_low,area_part_up,eps_broad_min,eps_broad_start,eps_broad_max,fit_region,doping)

% Lorentzian and Gaussian functions are saved in a string "plotString"
% to create a function handle later on. The function handle has the form:
%
% y(x) = Absoprtion - Lorentzian/Gaussian(x) - Phonon sideband(x).
% The code then tries to minimize y(x) to obtain the best fit.
%
% Lorentzian line shapes are defined as:
% y_L = Height / (1 + ((lambda-lambda_c)/(0.5*FWHM))^2)
% with lambda_c = peak position
%
% Gaussian line shapes are defined as:
% y_G = Height * exp(-log(2) * ((lambda-lambda_c )/(FWHM/2)).^2);
%
% Exciton phonon line shapes are defined as:
% y_G = area_L / area_G * spectral_weight * exp(-2 * (lambda - EPS position)^2 / (FWHM_G / sqrt(ln(4)) )^2) / (FWHM_G * sqrt(pi / (2 * ln(4))))
% with area_L = Height * pi * FWHM/ 2 or
% with area_G = Height * FWHM * sqrt(pi/(2*log(4)))
% and spectral weight transfer = 0.017 + 0.1 / diam_S11 + correction_factor

if isempty(fit_region)==0 % Fit S11 and S22
    switch fit_region
        case 1
            areas=area;
            FWHM=fwhm;
            center=Cent;
            center_22=[];
            phonon_pos2=[];
            variables=2*(length(center)-sum(center==0))+1;
        case 2
            areas=area_22';
            FWHM=FWHM_22;
            FWHM_g=FWHM_g_22;
            center=Cent_22;
            center_22=[];
            phonon_pos=phonon_pos2;
            phonon_pos2=[];
            f1=f1_22;
            variables=2*(length(center)-sum(center==0))+1;
        case 3
            areas=[area;area_22'];
            FWHM=fwhm;
            center=Cent;
            center_22=Cent_22;
            variables=2*(length(center)+length(center_22)-sum(center_22==0))+1;
            if isempty(phonon_pos)==0
                l1=1;
            else
                l1=0;
            end
            if isempty(phonon_pos2)==0
                l2=1;
            else
                l2=0;
            end
            if doping==1
                var=2*(length(center)+length(center_22)-sum(center_22==0))+2+2*(l1+l2);
                w1_old=sum(area/sum(areas));
                w2_old=sum(area_22/sum(areas));
            end
    end
    w_sol=areas/sum(areas);
else
    w_sol=w_solution;
    FWHM=fwhm;
    center=Cent;
    center_22=Cent_22;
    variables=2*(length(center)+length(center_22)-sum(center_22==0))+1;
end

area_film=trapz(lambda,Film); % initial guess of entire film area to have a starting value for the area of the most intense peak

lb=[];
ub=[];
x0=[];
plotString='';

n=0;

for j=1:length(center)
    if center(j)>0
        if j==1
            if w==1 % Lorentzian
                plotString = strcat(plotString, ['(2*x(',num2str(1),')/pi*',num2str(FWHM(j)),'*x(',num2str(variables),')./(4*([',num2str(lambda,' %f32'),']-',num2str(center(j)),'-x(',num2str(2),')).^2+(',num2str(FWHM(j)),'*x(',num2str(variables),')).^2)']);
            else % Gaussian
                plotString = strcat(plotString,['(x(',num2str(1),')*exp(-4*log(2)*([',num2str(lambda,' %f32'),']-x(',num2str(2),')-',num2str(center(j)),').^2/(',num2str(FWHM(j)),'*x(',num2str(variables),'))^2)/((',num2str(FWHM(j)),'*x(',num2str(variables),'))*sqrt(pi/(4*log(2))))']);
            end
            
            % x(1)=area of first peak, x(2)=shift of peak position,
            % x(2*length(center)+1)=overall broadening of FWHM
            
            x0=[double(area_film*w_sol(1)*0.8),start_shift];
            lb=[0,min_shift];
            ub=[double(area_film),max_shift];
            
        else
            if doping==0
                if w==1 % Lorentzian
                    plotString = strcat(plotString, ['+ 2*x(',num2str(1),')/(',num2str(w_sol(1)),')*(',num2str(w_sol(j)),')*x(',num2str(2*(j-n)-1),')/pi*',num2str(FWHM(j)),'*x(',num2str(variables),')./(4*([',num2str(lambda,' %f32'),']-',num2str(center(j)),'-x(',num2str(2*(j-n)),')).^2+(',num2str(FWHM(j)),'*x(',num2str(variables),')).^2)']);
                else % Gaussian
                    plotString = strcat(plotString,['+ x(',num2str(1),')/(',num2str(w_sol(1)),')*(',num2str(w_sol(j)),')*x(',num2str(2*(j-n)-1),')*exp(-4*log(2)*([',num2str(lambda,' %f32'),']-x(',num2str(2*(j-n)),')-',num2str(center(j)),').^2/(',num2str(FWHM(j)),'*x(',num2str(variables),'))^2)/((',num2str(FWHM(j)),'*x(',num2str(variables),'))*sqrt(pi/(4*log(2))))']);
                end
            else
                if w==1 % Lorentzian
                    plotString = strcat(plotString, ['+ x(',num2str(var),')/(',num2str(w1_old),')*2*x(',num2str(1),')/(',num2str(w_sol(1)),')*(',num2str(w_sol(j)),')*x(',num2str(2*(j-n)-1),')/pi*',num2str(FWHM(j)),'*x(',num2str(variables),')./(4*([',num2str(lambda,' %f32'),']-',num2str(center(j)),'-x(',num2str(2*(j-n)),')).^2+(',num2str(FWHM(j)),'*x(',num2str(variables),')).^2)']);
                else % Gaussian
                    plotString = strcat(plotString,['+ x(',num2str(var),')/(',num2str(w1_old),')*x(',num2str(1),')/(',num2str(w_sol(1)),')*(',num2str(w_sol(j)),')*x(',num2str(2*(j-n)-1),')*exp(-4*log(2)*([',num2str(lambda,' %f32'),']-x(',num2str(2*(j-n)),')-',num2str(center(j)),').^2/(',num2str(FWHM(j)),'*x(',num2str(variables),'))^2)/((',num2str(FWHM(j)),'*x(',num2str(variables),'))*sqrt(pi/(4*log(2))))']);
                end
            end
            
            % x(2*j-1)=change in relative concentration compared to
            % solution fit (e.g., +- 10 %), x(2*j)=shift of peak
            
            x0=[x0,1,start_shift];
            lb=[lb,1-conc_change/100,min_shift];
            ub=[ub,1+conc_change/100,max_shift];
            
        end
    else
        n=n+1;
    end
end

if fit_region==3
    n=0;
    for j=1:length(center_22)
        if center_22(j)>0
            if doping==0
                if w==1 % Lorentzian
                    plotString = strcat(plotString, ['+ 2*x(',num2str(1),')/(',num2str(w_sol(1)),')*(',num2str(w_sol(j+length(center))),')*x(',num2str(2*(j+length(center)-n)-1),')/pi*',num2str(FWHM_22(j)),'*x(',num2str(variables),')./(4*([',num2str(lambda,' %f32'),']-',num2str(center_22(j)),'-x(',num2str(2*(j+length(center)-n)),')).^2+(',num2str(FWHM_22(j)),'*x(',num2str(variables),')).^2)']);
                else % Gaussian
                    plotString = strcat(plotString,['+ x(',num2str(1),')/(',num2str(w_sol(1)),')*(',num2str(w_sol(j+length(center))),')*x(',num2str(2*(j+length(center)-n)-1),')*exp(-4*log(2)*([',num2str(lambda,' %f32'),']-x(',num2str(2*(j+length(center)-n)),')-',num2str(center_22(j)),').^2/(',num2str(FWHM_22(j)),'*x(',num2str(variables),'))^2)/((',num2str(FWHM_22(j)),'*x(',num2str(variables),'))*sqrt(pi/(4*log(2))))']);
                end
            else
                if w==1 % Lorentzian
                    plotString = strcat(plotString, ['+ (1-x(',num2str(var),'))/(',num2str(w2_old),')*2*x(',num2str(1),')/(',num2str(w_sol(1)),')*(',num2str(w_sol(j+length(center))),')*x(',num2str(2*(j+length(center)-n)-1),')/pi*',num2str(FWHM_22(j)),'*x(',num2str(variables),')./(4*([',num2str(lambda,' %f32'),']-',num2str(center_22(j)),'-x(',num2str(2*(j+length(center)-n)),')).^2+(',num2str(FWHM_22(j)),'*x(',num2str(variables),')).^2)']);
                else % Gaussian
                    plotString = strcat(plotString,['+ (1-x(',num2str(var),'))/(',num2str(w2_old),')*x(',num2str(1),')/(',num2str(w_sol(1)),')*(',num2str(w_sol(j+length(center))),')*x(',num2str(2*(j+length(center)-n)-1),')*exp(-4*log(2)*([',num2str(lambda,' %f32'),']-x(',num2str(2*(j+length(center)-n)),')-',num2str(center_22(j)),').^2/(',num2str(FWHM_22(j)),'*x(',num2str(variables),'))^2)/((',num2str(FWHM_22(j)),'*x(',num2str(variables),'))*sqrt(pi/(4*log(2))))']);
                end
            end
            
            % x(2*j-1)=change in relative concentration compared to
            % solution fit (e.g., +- 10 %), x(2*j)=shift of peak
            
            x0=[x0,1,start_shift];
            lb=[lb,1-conc_change/100,min_shift];
            ub=[ub,1+conc_change/100,max_shift];
            
        else
            n=n+1;
        end
    end
end

x0=[x0,broad_start];   % Initial guess of peak broadening
lb=[lb,broad_min]; % Lower limit
ub=[ub,broad_max]; % Upper limit

count=1;
kk=0;

for j=1:length(phonon_pos)
    
    if doping==0
        if phonon_pos(j)==1
            plotString = strcat(plotString,['+ x(',num2str(1),')/(',num2str(w_sol(1)),')*(',num2str(w_sol(phonon_pos(j))),')*(0.017+0.1/',num2str(diam_n(phonon_pos(j))),' + x(',num2str(variables+1),'))*exp(-4*log(2)*([',num2str(lambda,' %f32'),']-',num2str(h*c),'/(1e-9*(',num2str(h*c),'/((',num2str(center(phonon_pos(j))),'+x(',num2str(2*phonon_pos(j)),'))*1e-9)+0.2))).^2/(',num2str(FWHM_g(j)),'*x(',num2str(variables+2),'))^2)/((',num2str(FWHM_g(j)),'*x(',num2str(variables+2),'))*sqrt(pi/(4*log(2))))']);
        else
            plotString = strcat(plotString,['+ x(',num2str(1),')/(',num2str(w_sol(1)),')*(',num2str(w_sol(phonon_pos(j))),')*x(',num2str(2*phonon_pos(j)-1),')*(0.017+0.1/',num2str(diam_n(phonon_pos(j))),' + x(',num2str(variables+1),'))*exp(-4*log(2)*([',num2str(lambda,' %f32'),']-',num2str(h*c),'/(1e-9*(',num2str(h*c),'/((',num2str(center(phonon_pos(j))),'+x(',num2str(2*phonon_pos(j)),'))*1e-9)+0.2))).^2/(',num2str(FWHM_g(j)),'*x(',num2str(variables+2),'))^2)/((',num2str(FWHM_g(j)),'*x(',num2str(variables+2),'))*sqrt(pi/(4*log(2))))']);
        end
    else
        if phonon_pos(j)==1
            plotString = strcat(plotString,['+ x(',num2str(1),')/(',num2str(w_sol(1)),')*(',num2str(w_sol(phonon_pos(j))),')*(0.017+0.1/',num2str(diam_n(phonon_pos(j))),' + x(',num2str(variables+1),'))*exp(-4*log(2)*([',num2str(lambda,' %f32'),']-',num2str(h*c),'/(1e-9*(',num2str(h*c),'/((',num2str(center(phonon_pos(j))),'+x(',num2str(2*phonon_pos(j)),'))*1e-9)+0.2))).^2/(',num2str(FWHM_g(j)),'*x(',num2str(variables+2),'))^2)/((',num2str(FWHM_g(j)),'*x(',num2str(variables+2),'))*sqrt(pi/(4*log(2))))']);
        else
            plotString = strcat(plotString,['+ x(',num2str(var),')/(',num2str(w1_old),')*x(',num2str(1),')/(',num2str(w_sol(1)),')*(',num2str(w_sol(phonon_pos(j))),')*x(',num2str(2*phonon_pos(j)-1),')*(0.017+0.1/',num2str(diam_n(phonon_pos(j))),' + x(',num2str(variables+1),'))*exp(-4*log(2)*([',num2str(lambda,' %f32'),']-',num2str(h*c),'/(1e-9*(',num2str(h*c),'/((',num2str(center(phonon_pos(j))),'+x(',num2str(2*phonon_pos(j)),'))*1e-9)+0.2))).^2/(',num2str(FWHM_g(j)),'*x(',num2str(variables+2),'))^2)/((',num2str(FWHM_g(j)),'*x(',num2str(variables+2),'))*sqrt(pi/(4*log(2))))']);
        end
    end
    
    count=count+1;
    if count==length(phonon_pos)+1
        x0=[x0,f1(1),eps_broad_start]; % Gauss broadening
        lb=[lb,area_part_low,eps_broad_min];
        ub=[ub,area_part_up,eps_broad_max];
        kk=2;
    end
end

count=1;

for j=1:length(phonon_pos2)
    
    if doping==0
        plotString = strcat(plotString,['+ x(',num2str(1),')/(',num2str(w_sol(1)),')*(',num2str(w_sol(length(center)+phonon_pos2(j))),')*x(',num2str(2*(length(center)+phonon_pos2(j))-1),')*(0.017+0.1/',num2str(diam_n(phonon_pos2(j))),' + x(',num2str(variables+1+kk),'))*exp(-4*log(2)*([',num2str(lambda,' %f32'),']-',num2str(h*c),'/(1e-9*(',num2str(h*c),'/((',num2str(center_22(phonon_pos2(j))),'+x(',num2str(2*(length(center)+phonon_pos2(j))),'))*1e-9)+0.2))).^2/(',num2str(FWHM_g_22(j)),'*x(',num2str(variables+2+kk),'))^2)/((',num2str(FWHM_g_22(j)),'*x(',num2str(variables+2+kk),'))*sqrt(pi/(4*log(2))))']);
    else
        plotString = strcat(plotString,['+ (1-x(',num2str(var),'))/(',num2str(w2_old),')*x(',num2str(1),')/(',num2str(w_sol(1)),')*(',num2str(w_sol(length(center)+phonon_pos2(j))),')*x(',num2str(2*(length(center)+phonon_pos2(j))-1),')*(0.017+0.1/',num2str(diam_n(phonon_pos2(j))),' + x(',num2str(variables+1+kk),'))*exp(-4*log(2)*([',num2str(lambda,' %f32'),']-',num2str(h*c),'/(1e-9*(',num2str(h*c),'/((',num2str(center_22(phonon_pos2(j))),'+x(',num2str(2*(length(center)+phonon_pos2(j))),'))*1e-9)+0.2))).^2/(',num2str(FWHM_g_22(j)),'*x(',num2str(variables+2+kk),'))^2)/((',num2str(FWHM_g_22(j)),'*x(',num2str(variables+2+kk),'))*sqrt(pi/(4*log(2))))']);
    end
    
    count=count+1;
    if count==length(phonon_pos2)+1
        x0=[x0,f1_22(1),eps_broad_start]; % Gauss broadening
        lb=[lb,area_part_low,eps_broad_min];
        ub=[ub,area_part_up,eps_broad_max];
    end
    
end

if doping==1
    x0=[x0,0.5]; 
    lb=[lb,0];
    ub=[ub,0.85];
end

plotString=strcat(plotString,['));']);

plotString=strcat(['@(x) double([',num2str(Film,' %f32'),'] -'],plotString);


%--------------------------------------------------------------------------
% Section 2.2 - "generate_start_values_voigt"
%--------------------------------------------------------------------------


function [x0,lb,ub,w_sol]=generate_start_values_voigt(Cent,Cent_22,lambda,Film,f1,w_solution,phonon_pos,area,area_22,min_shift,start_shift,max_shift,conc_change,phonon_pos2,f1_22,area_part_low,area_part_up,eps_broad_min,eps_broad_start,eps_broad_max,broad_start,broad_min,broad_max,doping)
global fit_region

% Unlike in the solution absorption measurement fit, the Voigtian functions
% are defined based on their area, not their height to account for relative
% concentration changes of spectral weight.

area_film=trapz(lambda,Film); % initial guess of entire film area to have a starting value for the area of the most intense peak

if isempty(fit_region)==0 % Fit S11 and S22
    switch fit_region
        case 1
            areas=area;
            center=Cent;
            phonon_pos2=[];
        case 2
            areas=area_22';
            center=Cent_22;
            phonon_pos=phonon_pos2;
            phonon_pos2=[];
            f1=f1_22;
        case 3
            center=Cent;
            center_22=Cent_22;
            areas=[area;area_22'];
    end
    w_sol=areas/sum(areas);
else
    w_sol=w_solution;
    center=Cent;
    center_22=Cent_22;
end

lb=[];
ub=[];
x0=[];

for j=1:length(center)
    if center(j)>0
        if j==1
            % x(1)=area of first peak, x(2)=shift of peak position,
            % x(2*length(center)+1)=overall broadening of FWHM
            
            x0=[double(area_film*w_sol(1)*0.8),start_shift];
            lb=[0,min_shift];
            ub=[double(area_film),max_shift];
        else
            % x(2*j-1)=change in relative concentration compared to
            % solution fit (e.g., +- 10 %), x(2*j)=shift of peak
            
            x0=[x0,1,start_shift];
            lb=[lb,(1-conc_change/100),min_shift];
            ub=[ub,(1+conc_change/100),max_shift];
        end
    end
end

if fit_region==3
    for j=1:length(center_22)
        if center_22(j)>0
            % x(2*j-1)=change in relative concentration compared to
            % solution fit (e.g., +- 10 %), x(2*j)=shift of peak
            
            x0=[x0,area_film*w_sol(j+length(center)),start_shift];
            lb=[lb,area_film*w_sol(j+length(center))*(1-conc_change/100),min_shift];
            ub=[ub,area_film*w_sol(j+length(center))*(1+conc_change/100),max_shift];
        end
    end
end

x0=[x0,broad_start];% Initial guess of peak broadening
lb=[lb,broad_min];  % Lower limit
ub=[ub,broad_max];  % Upper limit

if isempty(phonon_pos)==0
    x0=[x0,f1(1),eps_broad_start];
    lb=[lb,area_part_low,eps_broad_min];
    ub=[ub,area_part_up,eps_broad_max];
end

if isempty(phonon_pos2)==0
    x0=[x0,f1_22(1),eps_broad_start];
    lb=[lb,area_part_low,eps_broad_min];
    ub=[ub,area_part_up,eps_broad_max];
end

if doping==1
    x0=[x0,sum(w_sol(1:length(center)))];
    lb=[lb,0];
    ub=[ub,.85];
end


%--------------------------------------------------------------------------
% Section 2.3 - "Voigt" and "complexErrorFunction"
%--------------------------------------------------------------------------


function diff=Voigt(x0,lambda_n,Film_n,w_solution,center,center_22,FWHM,phonon_pos,phonon_pos2,diam_n,FWHM_g,h,c,FWHM_22,FWHM_g_22)
global fit_region

if isempty(fit_region)==0 % Fit S11 and S22
    switch fit_region
        case 1
            variables=2*(length(center)-sum(center==0))+1;
            Cent=center;
            fwhm=FWHM/2;
            fwhm_g=FWHM_g;
            pp=phonon_pos;
            pp2=[];
            y=zeros(length(lambda_n),length(center)+length(phonon_pos)); % Initialization of empty vector "y" to increase calculation speed.
            dope=0;
        case 2
            variables=2*(length(center_22)-sum(center_22==0))+1;
            Cent=center_22;
            fwhm=FWHM_22;
            fwhm_g=FWHM_g_22;
            pp=phonon_pos2;
            pp2=[];
            y=zeros(length(lambda_n),length(center_22)-sum(center_22==0)+length(phonon_pos2));
            dope=0;
        case 3
            variables=2*(length(center)+length(center_22)-sum(center_22==0))+1;            
            Cent=center;
            Cent_22=center_22;
            fwhm=FWHM/2;
            fwhm_22=FWHM_22;
            fwhm_g=FWHM_g;
            fwhm_g_22=FWHM_g_22;
            pp=phonon_pos;
            pp2=phonon_pos2;
            y=zeros(length(lambda_n),length(center)+length(center_22)-sum(center_22==0)+length(phonon_pos)+length(phonon_pos2));
            if isempty(phonon_pos)==0
                l1=1;
            else
                l1=0;
            end
            if isempty(phonon_pos2)==0
                l2=1;
            else
                l2=0;
            end
            if length(x0)>variables+2*(l1+l2)
                dope=1;
                var=length(x0);
                w1_old=sum(w_solution(1:length(center)));
                w2_old=sum(w_solution(length(center)+1:end));
            else
                dope=0;
            end
    end
else
    variables=2*(length(center)+length(center_22)-sum(center_22==0))+1;
    dope=0;
    Cent=center;
    Cent_22=center_22;
    fwhm=FWHM/2;
    fwhm_22=FWHM_22;
    fwhm_g=FWHM_g;
    fwhm_g_22=FWHM_g_22;
    pp=phonon_pos;
    pp2=phonon_pos2;
    y=zeros(length(lambda_n),length(center)+length(center_22)-sum(center_22==0)+length(phonon_pos)+length(phonon_pos2));
end

kk=0;

n=0;

for j=1:length(Cent) % Either S11 or S22
    if Cent(j)>0
        wL(j)=fwhm(j,3)*x0(variables)/(0.5436+sqrt((fwhm(j,2)/fwhm(j,1))^2+.2166));
        wg(j)=fwhm(j,2)/fwhm(j,1)*wL(j);
        XX=sqrt(log(2))*(lambda_n-Cent(j)-x0(2*(j-n)))/(wg(j));
        YY=sqrt(log(2))*fwhm(j,1)/fwhm(j,2);
        if dope==0
            if j==1
                y(:,j-n)=x0(1)*sqrt(log(2)/pi)/(wg(j)).*real(complexErrorFunction(XX,YY));
            else
                y(:,j-n)=x0(2*(j-n)-1)*x0(1)/w_solution(1)*w_solution(j)*sqrt(log(2)/pi)/(wg(j)).*real(complexErrorFunction(XX,YY));
            end
        else
            if j==1
                y(:,j-n)=x0(1)*sqrt(log(2)/pi)/(wg(j)).*real(complexErrorFunction(XX,YY));
            else
                y(:,j-n)=x0(var)/w1_old*x0(2*(j-n)-1)*x0(1)/w_solution(1)*w_solution(j)*sqrt(log(2)/pi)/(wg(j)).*real(complexErrorFunction(XX,YY));
            end
        end
        kk=j;
    else
        n=n+1;
    end
end

kk2=kk;

n=0;
if fit_region==3
    for j=1:length(Cent_22) % S11 and S22
        if Cent_22(j)>0
            wL(j+kk2-n)=fwhm_22(j,3)*x0(variables)/(0.5436+sqrt((fwhm_22(j,2)/fwhm_22(j,1))^2+.2166));
            wg(j+kk2-n)=fwhm_22(j,2)/fwhm_22(j,1)*wL(j+kk2-n);
            XX=sqrt(log(2))*(lambda_n-Cent_22(j)-x0(2*(j+kk2-n)))/(wg(j+kk2-n));
            YY=sqrt(log(2))*fwhm_22(j,1)/fwhm_22(j,2);
            if dope==0
                y(:,j+kk2-n)=x0(2*(j+kk2-n)-1)*sqrt(log(2)/pi)/(wg(j+kk2-n)).*real(complexErrorFunction(XX,YY));
            else
                y(:,j+kk2-n)=(1-x0(var))/w2_old*x0(2*(j+kk2-n)-1)*sqrt(log(2)/pi)/(wg(j+kk2-n)).*real(complexErrorFunction(XX,YY));
            end
            kk=kk2+j;
        else
            n=n+1;
        end
    end
end

kk2=kk;
ph=0;

for j=1:length(pp)
    if dope==0
        if pp(j)==1
            y(:,kk2+j-n)=x0(1)*(0.017+0.1/diam_n(pp(j))+x0(variables+1))*exp(-2*(lambda_n-h*c/(1e-9*(h*c/((Cent(pp(j))+x0(pp(j)*2))*1e-9)+0.2))).^2/(fwhm_g(j)*x0(variables+2)/sqrt(log(4)))^2)/(fwhm_g(j)*x0(variables+2)*sqrt(pi/(2*log(4))));
        else
            y(:,kk2+j-n)=x0(2*pp(j)-1)*x0(1)/w_solution(1)*w_solution(pp(j))*(0.017+0.1/diam_n(pp(j))+x0(variables+1))*exp(-2*(lambda_n-h*c/(1e-9*(h*c/((Cent(pp(j))+x0(pp(j)*2))*1e-9)+0.2))).^2/(fwhm_g(j)*x0(variables+2)/sqrt(log(4)))^2)/(fwhm_g(j)*x0(variables+2)*sqrt(pi/(2*log(4))));
        end
    else
        if pp(j)==1
            y(:,kk2+j-n)=x0(1)*(0.017+0.1/diam_n(pp(j))+x0(variables+1))*exp(-2*(lambda_n-h*c/(1e-9*(h*c/((Cent(pp(j))+x0(pp(j)*2))*1e-9)+0.2))).^2/(fwhm_g(j)*x0(variables+2)/sqrt(log(4)))^2)/(fwhm_g(j)*x0(variables+2)*sqrt(pi/(2*log(4))));
        else
            y(:,kk2+j-n)=x0(var)/w1_old*x0(2*pp(j)-1)*x0(1)/w_solution(1)*w_solution(pp(j))*(0.017+0.1/diam_n(pp(j))+x0(variables+1))*exp(-2*(lambda_n-h*c/(1e-9*(h*c/((Cent(pp(j))+x0(pp(j)*2))*1e-9)+0.2))).^2/(fwhm_g(j)*x0(variables+2)/sqrt(log(4)))^2)/(fwhm_g(j)*x0(variables+2)*sqrt(pi/(2*log(4))));
        end
    end
    
    kk=kk2+j;
    ph=2;
end

kk2=kk;

for j=1:length(pp2)
    if dope==0
        y(:,kk2+j-n)=x0(2*(pp2(j)+length(center))-1)*(0.017+0.1/diam_n(pp2(j))+x0(variables+1+ph))*exp(-2*(lambda_n-h*c/(1e-9*(h*c/((Cent_22(pp2(j))+x0((pp2(j)+length(center))*2))*1e-9)+0.2))).^2/(fwhm_g_22(j)*x0(variables+2+ph)/sqrt(log(4)))^2)/(fwhm_g_22(j)*x0(variables+2+ph)*sqrt(pi/(2*log(4))));
    else
        y(:,kk2+j-n)=(1-x0(var))/w2_old*x0(2*(pp2(j)+length(center))-1)*(0.017+0.1/diam_n(pp2(j))+x0(variables+1+ph))*exp(-2*(lambda_n-h*c/(1e-9*(h*c/((Cent_22(pp2(j))+x0((pp2(j)+length(center))*2))*1e-9)+0.2))).^2/(fwhm_g_22(j)*x0(variables+2+ph)/sqrt(log(4)))^2)/(fwhm_g_22(j)*x0(variables+2+ph)*sqrt(pi/(2*log(4))));
    end
end

diff=double(Film_n'-sum(y,2));


function [w] = complexErrorFunction(x,y)
% complexErrorFunction - Calculation of complex error function using
% dimensionless coordinates
%
% [w] = complexErrorFunction(x,y)   Computes the complex error function
%   using the algorithm developed by Dr. F. Schreier and kindly presented
%   in Fortran. The function was rewriten to MATLAB by Dr. N. Cherkasov
%   For more details on algorithm see the publication:
%   F. Schreier: Optimized Implementations of Rational Approximations for
%                the Voigt ane Complex Error Function.
%   J. Quant. Spectrosc. & Radiat. Transfer, 112(6), 10101025,
%   doi 10.1016/j.jqsrt.2010.12.010, 2011.
%
%   Briefly, the algorithm is compiled from two:
%       for    large x+y     J  Humlicek, JQSRT 27, 437, 1982
%       for small x+y:    J.A.C. Weideman,  SIAM J. Numer. Anal. 31 (1994)
%       pp. 1497-1518,  equation (38.I) and table I
%
% INPUT ARGUMENTS are dimensionless coordinates x and y
%   x - array 1*N, and y - single variable
%
% OUTPUT
%   w - complex array 1*N
%
% The function was used for the deconvolution of IR spectra
% see the publication
%
% 27-December-2013 N. Cherkasov
% Comments and questions to: n.b.cherkasov@gmail.com


half=0.5;
one=1;
two=2;
recSqrtPi=1/sqrt(pi);

lenX=length(x);
w=zeros(lenX,1);

%   "Weideman" constants
% n=24;
l=4.1195342878142354; % l=sqrt(n/sqrt(2.))  ! L = 2**(-1/4) * N**(1/2)
a=[-1.5137461654527820e-10,  4.9048215867870488e-09,  1.3310461806370372e-09, -3.0082822811202271e-08, ...
    -1.9122258522976932e-08,  1.8738343486619108e-07,  2.5682641346701115e-07, -1.0856475790698251e-06, ...
    -3.0388931839840047e-06,  4.1394617248575527e-06,  3.0471066083243790e-05,  2.4331415462641969e-05, ...
    -2.0748431511424456e-04, -7.8166429956142650e-04, -4.9364269012806686e-04,  6.2150063629501763e-03, ...
    3.3723366855316413e-02,  1.0838723484566792e-01,  2.6549639598807689e-01,  5.3611395357291292e-01, ...
    9.2570871385886788e-01,  1.3948196733791203e+00,  1.8562864992055408e+00,  2.1978589365315417e+00];
%   humlicek prbFct region I bounds
s15=15e0;
% -------------------------------------------------------------------------
x12 = y - s15; %        left wing -- center
x21 = -x12;    % 15-y   center -- right wing

if (y>s15 || x(1)>x21 || x(lenX)<x12)
    %       all points are in Humlicek's region I
    for ii=1:lenX
        t= y - x(ii)*1i;
        w(ii) = (recSqrtPi*t) / (half + t*t);
    end
else
    for ii=1:lenX
        s  = abs(x(ii)) + y;
        if (s>s15)
            t     = y-x(ii)*1i;
            w(ii) = (recSqrtPi*t) / (half + t*t);
        else
            recLmZ  = one / (l+y-x(ii)*1i);
            t       = (l-y+x(ii)*1i) * recLmZ;
            w(ii) =  recLmZ  *  (recSqrtPi + two*recLmZ*...
                (a(24)+(a(23)+(a(22)+(a(21)+(a(20)+(a(19)+(a(18)+(a(17)+(a(16)+(a(15)+(a(14)+(a(13)+(a(12)+(a(11)+(a(10)+(a(9)+...
                (a(8)+(a(7)+(a(6)+(a(5)+(a(4)+(a(3)+(a(2)+a(1)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t));
        end
    end
end


%..........................................................................


%--------------------------------------------------------------------------
% Section 3 - "plot_data"
%--------------------------------------------------------------------------


function [area_f,y,legstrs]=plot_data(w,x_fit,FWHM,FWHM_22,lambda_n,Film,w_solution,phonon_pos,phonon_pos2,diam_n,FWHM_g,legstr,colors,sse,h,c,FWHM_g_22,center,center_22,legstr_22,fit_region,colors_m,doping)

% Based on the line-shape used for solution fit and the results of the
% fit performed in Section 2 (stored in "x_fit"), the individual (n,m)
% species are plotted along with the calculated spectrum (sum of individual
% (n,m) fits) and the absorption spectrum initially loaded by the user. If
% Lorentzian or Gaussian line-shapes were defined, they are plotted with
% "plot_fit" defined in Section 3.1. If Voigtian line-shapes were chosen,
% they are plotted with "plot_fit_voigt" defined in Section 3.2.

switch w
    case {1,2}
        [area_f,y,legstrs]=plot_fit(x_fit,FWHM,lambda_n,center,center_22,Film,w_solution,phonon_pos,phonon_pos2,diam_n,w,FWHM_g,legstr,colors,sse,h,c,FWHM_g_22,legstr_22,FWHM_22,fit_region,colors_m,doping);
    case 3
        [area_f,y,legstrs]=plot_fit_voigt(x_fit,FWHM,lambda_n,center,center_22,Film,phonon_pos,phonon_pos2,diam_n,FWHM_g,legstr,colors,sse,h,c,FWHM_g_22,legstr_22,FWHM_22,w_solution,colors_m);
end


%--------------------------------------------------------------------------
% Section 3.1 - "plot_fit"
%--------------------------------------------------------------------------


function [area_f,y,legstrs]=plot_fit(x_fit,fwhm,lambda,Cent,Cent_22,Film,w_sol,pp,pp2,diam_n,w,FWHM_g,legstr,colors,sse,h,c,FWHM_g_22,legstr_22,FWHM_22,fit_region,colors_m,doping)

% This subfunction is used to plot the result of a Gaussian or Lorentzian
% fit. The structure of the code follows the rountine outlined in Section
% 2.1 "generate_string".

if isempty(fit_region)==0 % Fit S11 and S22
    switch fit_region
        case 1
            FWHM=fwhm;
            center=Cent;
            center_22=[];
            phonon_pos=pp;
            phonon_pos2=[];
            y=zeros(length(lambda),length(center)+length(phonon_pos));
            area_f=zeros(length(center),1);
            variables=2*(length(center))+1;
            legstr_22=[];
        case 2
            FWHM=FWHM_22;
            FWHM_g=FWHM_g_22;
            center=Cent_22;
            center_22=[];
            phonon_pos=pp2;
            phonon_pos2=[];
            y=zeros(length(lambda),length(center)-sum(center==0)+length(phonon_pos));
            area_f=zeros(length(center),1);
            variables=2*(length(center)-sum(center==0))+1;
            legstr=legstr_22;
            legstr_22=[];
        case 3
            FWHM=fwhm;
            center=Cent;
            center_22=Cent_22;
            phonon_pos=pp;
            phonon_pos2=pp2;
            y=zeros(length(lambda),length(center)+length(center_22)-sum(center_22==0)+length(phonon_pos)+length(phonon_pos2));
            area_f=zeros(length(center)+length(center_22),1);
            variables=2*(length(center)+length(center_22)-sum(center_22==0))+1;
            w1_old=sum(w_sol(1:length(center)));
            w2_old=sum(w_sol(length(center)+1:end));
    end
else
    FWHM=fwhm;
    center=Cent;
    center_22=Cent_22;
    phonon_pos=pp;
    phonon_pos2=pp2;
    area_f=zeros(length(center)+length(center_22),1);
    variables=2*(length(center)+length(center_22)-sum(center_22==0))+1;
    y=zeros(length(lambda),length(center)+length(center_22)-sum(center_22==0)+length(phonon_pos)+length(phonon_pos2));
end

figure, hold on

n=0;

for j=1:length(center)
    if center(j)>0
        if j==1
            if w==1 % Lorentzian
                y(:,j-n)=2*x_fit(1)/pi*FWHM(j)*x_fit(variables)./(4*(lambda-center(j)-x_fit(2*(j-n))).^2+(FWHM(j)*x_fit(variables)).^2);
                title_name='Lorentzian';
            else % Gaussian
                y(:,j-n)=x_fit(1)*exp(-4*log(2)*(lambda-x_fit(2)-center(j)).^2/(FWHM(j)*x_fit(variables))^2)/(FWHM(j)*x_fit(variables)*sqrt(pi/(4*log(2))));
                title_name='Gaussian';
            end
            area_f(j)=x_fit(1);
        else
            if doping==0
                if w==1
                    y(:,j-n)=2*x_fit(1)/w_sol(1)*w_sol(j)*x_fit(2*(j-n)-1)/pi*FWHM(j)*x_fit(variables)./(4*(lambda-center(j)-x_fit(2*(j-n))).^2+(FWHM(j)*x_fit(variables)).^2);
                else
                    y(:,j-n)=x_fit(1)/w_sol(1)*w_sol(j)*x_fit(2*(j-n)-1)*exp(-4*log(2)*(lambda-x_fit(2*(j-n))-center(j)).^2/(FWHM(j)*x_fit(variables))^2)/(FWHM(j)*x_fit(variables)*sqrt(pi/(4*log(2))));
                end
                area_f(j-n)=x_fit(1)/w_sol(1)*w_sol(j)*x_fit(2*(j-n)-1);
            else
                if w==1
                    y(:,j-n)=x_fit(end)/w1_old*2*x_fit(1)/w_sol(1)*w_sol(j)*x_fit(2*(j-n)-1)/pi*FWHM(j)*x_fit(variables)./(4*(lambda-center(j)-x_fit(2*(j-n))).^2+(FWHM(j)*x_fit(variables)).^2);
                else
                    y(:,j-n)=x_fit(end)/w1_old*x_fit(1)/w_sol(1)*w_sol(j)*x_fit(2*(j-n)-1)*exp(-4*log(2)*(lambda-x_fit(2*(j-n))-center(j)).^2/(FWHM(j)*x_fit(variables))^2)/(FWHM(j)*x_fit(variables)*sqrt(pi/(4*log(2))));
                end
                area_f(j-n)=x_fit(end)/w1_old*x_fit(1)/w_sol(1)*w_sol(j)*x_fit(2*(j-n)-1);
            end
        end
        
        if j>length(colors)
            plot(lambda,y(:,j-n),'color',colors_m(j-length(colors),:),'linewidth',1)
        else
            plot(lambda,y(:,j-n),'color',colors(j,:),'linewidth',1)
        end
    else
        n=n+1;
    end
end

if fit_region==3
    n=0;
    for j=1:length(center_22)
        if center_22(j)>0
            if doping==0
                if w==1
                    y(:,j+length(center)-n)=2*x_fit(1)/w_sol(1)*w_sol(j+length(center))*x_fit(2*(j+length(center)-n)-1)/pi*FWHM_22(j)*x_fit(variables)./(4*(lambda-center_22(j)-x_fit(2*(j+length(center)-n))).^2+(FWHM_22(j)*x_fit(variables)).^2);
                else
                    y(:,j+length(center)-n)=x_fit(1)/w_sol(1)*w_sol(j+length(center))*x_fit(2*(j+length(center)-n)-1)*exp(-4*log(2)*(lambda-x_fit(2*(j+length(center)-n))-center_22(j)).^2/(FWHM_22(j)*x_fit(variables))^2)/(FWHM_22(j)*x_fit(variables)*sqrt(pi/(4*log(2))));
                end
                area_f(j+length(center)-n)=x_fit(1)/w_sol(1)*w_sol(j+length(center))*x_fit(2*(j+length(center)-n)-1);
            else 
                if w==1
                    y(:,j+length(center)-n)=(1-x_fit(end))/w2_old*2*x_fit(1)/w_sol(1)*w_sol(j+length(center))*x_fit(2*(j+length(center)-n)-1)/pi*FWHM_22(j)*x_fit(variables)./(4*(lambda-center_22(j)-x_fit(2*(j+length(center)-n))).^2+(FWHM_22(j)*x_fit(variables)).^2);
                else
                    y(:,j+length(center)-n)=(1-x_fit(end))/w2_old*x_fit(1)/w_sol(1)*w_sol(j+length(center))*x_fit(2*(j+length(center)-n)-1)*exp(-4*log(2)*(lambda-x_fit(2*(j+length(center)-n))-center_22(j)).^2/(FWHM_22(j)*x_fit(variables))^2)/(FWHM_22(j)*x_fit(variables)*sqrt(pi/(4*log(2))));
                end
                area_f(j+length(center)-n)=(1-x_fit(end))/w2_old*x_fit(1)/w_sol(1)*w_sol(j+length(center))*x_fit(2*(j+length(center)-n)-1);
            end
            
            if j>length(colors)
                plot(lambda,y(:,j+length(center)-n),'color',colors_m(j-length(colors),:),'linewidth',1)
            else
                plot(lambda,y(:,j+length(center)-n),'color',colors(j,:),'linewidth',1)
            end
            
        else
            n=n+1;
        end
    end
end

% Exciton phonon sidebands

ph=0;
legstr_ph=[];

for j=1:length(phonon_pos)
    if j==1
        y(:,length(center)+length(center_22)-n+j)=x_fit(1)/w_sol(1)*w_sol(phonon_pos(j))*(0.017+0.1/diam_n(phonon_pos(j))+x_fit(variables+1))*exp(-4*log(2)*(lambda-h*c/(1e-9*(h*c/((center(phonon_pos(j))+x_fit(2*phonon_pos(j)))*1e-9)+0.2))).^2/(FWHM_g(j)*x_fit(variables+2))^2)/(FWHM_g(j)*x_fit(variables+2)*sqrt(pi/(4*log(2))));
    else
        y(:,length(center)+length(center_22)-n+j)=x_fit(1)/w_sol(1)*w_sol(phonon_pos(j))*x_fit(2*phonon_pos(j)-1)*(0.017+0.1/diam_n(phonon_pos(j))+x_fit(variables+1))*exp(-4*log(2)*(lambda-h*c/(1e-9*(h*c/((center(phonon_pos(j))+x_fit(2*phonon_pos(j)))*1e-9)+0.2))).^2/(FWHM_g(j)*x_fit(variables+2))^2)/(FWHM_g(j)*x_fit(variables+2)*sqrt(pi/(4*log(2))));
    end
    ph=2;
    
    plot(lambda,y(:,j+length(center)+length(center_22)-n),'color',[colors(phonon_pos(j),1),.8,colors(phonon_pos(j),3)],'linewidth',1)
    
    legstr_ph{j}=['EPS-',legstr{phonon_pos(j)}];
end

for j=1:length(phonon_pos2)
    if j==1
        y(:,length(center)+length(center_22)-n+length(phonon_pos)+j)=area_f(length(center)+phonon_pos2(j))*(0.017+0.1/diam_n(phonon_pos2(j))+x_fit(variables+1+ph))*exp(-4*log(2)*(lambda-h*c/(1e-9*(h*c/((center_22(phonon_pos2(j))+x_fit(2*(length(center)+phonon_pos2(j))))*1e-9)+0.2))).^2/(FWHM_g_22(j)*x_fit(variables+2+ph))^2)/(FWHM_g_22(j)*x_fit(variables+2+ph)*sqrt(pi/(4*log(2))));
    else
        y(:,length(center)+length(center_22)-n+length(phonon_pos)+j)=area_f(length(center)+phonon_pos2(j))*(0.017+0.1/diam_n(phonon_pos2(j))+x_fit(variables+1+ph))*exp(-4*log(2)*(lambda-h*c/(1e-9*(h*c/((center_22(phonon_pos2(j))+x_fit(2*(length(center)+phonon_pos2(j))))*1e-9)+0.2))).^2/(FWHM_g_22(j)*x_fit(variables+2+ph))^2)/(FWHM_g_22(j)*x_fit(variables+2+ph)*sqrt(pi/(4*log(2))));
    end
    
    plot(lambda,y(:,j+length(center)+length(center_22)-n+length(phonon_pos)),'color',[colors(phonon_pos2(j),1),.8,colors(phonon_pos2(j),3)],'linewidth',1)
    legstr_ph{j+length(phonon_pos)}=['EPS-',legstr_22{phonon_pos2(j)}];
end

plot(lambda,Film,'k','linewidth',2)

plot(lambda,sum(y,2),'r','linewidth',2)

legstrs=[legstr legstr_22 legstr_ph];

legend([legstrs {'Measured Spectrum'} {'Calculated Spectrum'}])
title([title_name,' Fit with: nSSE = ',num2str(sse)],'FontSize',26)

xlabel('Wavelength (nm)','FontSize',24)
ylabel('Absorption (a.u.)','FontSize',24)
set(gca,'LineWidth',2)
set(gca,'FontSize',22)


%--------------------------------------------------------------------------
% Section 3.2 - "plot_fit_voigt"
%--------------------------------------------------------------------------


function [area_f,y,legstrs]=plot_fit_voigt(x0,fwhm,lambda_n,Cent,center_22,Film,phonon_pos,phonon_pos2,diam_n,fwhm_g,legstr,colors,sse,h,c,fwhm_g_22,legstr_22,FWHM_22,w_solution,colors_m)
global fit_region

% This subfunction is used to plot the result of a Voigtian fit. The
% structure of the code follows the rountine outlined in Section 2.3
% "Voigt".

dope=0;

if isempty(fit_region)==0 % Fit S11 and S22
    switch fit_region
        case 1
            center=Cent;
            variables=2*length(center)+1;
            FWHM=fwhm;
            FWHM_g=fwhm_g;
            pp=phonon_pos;
            pp2=[];
            y=zeros(length(lambda_n),length(center)+length(phonon_pos));
            area_f=zeros(length(center),1);
            legstr_22=[];
        case 2
            center=center_22;
            variables=2*(length(center)-sum(center==0))+1;
            FWHM=FWHM_22;
            FWHM_g=fwhm_g_22;
            pp=phonon_pos2;
            pp2=[];
            y=zeros(length(lambda_n),length(center_22)-sum(center_22==0)+length(phonon_pos2));
            area_f=zeros(length(center_22)-sum(center_22==0),1);
            legstr=legstr_22;
            legstr_22=[];
        case 3
            center=Cent;
            variables=2*(length(center)+length(center_22)-sum(center_22==0))+1;
            FWHM=fwhm;
            FWHM_g=fwhm_g;
            FWHM_g_22=fwhm_g_22;
            pp=phonon_pos;
            pp2=phonon_pos2;
            y=zeros(length(lambda_n),length(center)+length(center_22)-sum(center_22==0)+length(phonon_pos)+length(phonon_pos2));
            area_f=zeros(length(center)+length(center_22)-sum(center_22==0),1);
            
            if isempty(phonon_pos)==0
                l1=1;
            else
                l1=0;
            end
            if isempty(phonon_pos2)==0
                l2=1;
            else
                l2=0;
            end
            
            if length(x0)>variables+2*(l1+l2)
                dope=1;
                var=length(x0);
                w1_old=sum(w_solution(1:length(center)));
                w2_old=sum(w_solution(length(center)+1:end));
            end
    end
else
    center=Cent;
    variables=2*(length(center)+length(center_22)-sum(center_22==0))+1;
    FWHM=fwhm;
    FWHM_g=fwhm_g;
    FWHM_g_22=fwhm_g_22;
    pp=phonon_pos;
    pp2=phonon_pos2;
    y=zeros(length(lambda_n),length(center)+length(center_22)-sum(center_22==0)+length(phonon_pos)+length(phonon_pos2));
    area_f=zeros(length(center)+length(center_22)-sum(center_22==0),1);
end

kk=0;

figure, hold on

n=0;
for j=1:length(center) % Fit either S11 or S22
    if center(j)>0
        XX=sqrt(log(2))*(lambda_n-center(j)-x0(2*(j-n)))/(0.5*FWHM(j,2)*x0(variables));
        YY=sqrt(log(2))*FWHM(j,1)/FWHM(j,2);
        if dope==0
            if j==1
                y(:,j-n)=x0(1)*sqrt(log(2)/pi)/(0.5*FWHM(j,2)*x0(variables)).*real(complexErrorFunction(XX,YY));
                area_f(1)=x0(1);
            else
                y(:,j-n)=x0(2*(j-n)-1)*x0(1)/w_solution(1)*w_solution(j)*sqrt(log(2)/pi)/(0.5*FWHM(j,2)*x0(variables)).*real(complexErrorFunction(XX,YY));
                area_f(j)=x0(2*(j-n)-1)*x0(1)/w_solution(1)*w_solution(j);
            end
        else
            if j==1
                y(:,j-n)=x0(1)*sqrt(log(2)/pi)/(0.5*FWHM(j,2)*x0(variables)).*real(complexErrorFunction(XX,YY));
                area_f(1)=x0(1);
            else
                y(:,j-n)=x0(var)/w1_old*x0(2*(j-n)-1)*x0(1)/w_solution(1)*w_solution(j)*sqrt(log(2)/pi)/(0.5*FWHM(j,2)*x0(variables)).*real(complexErrorFunction(XX,YY));
                area_f(j)=x0(var)/w1_old*x0(2*(j-n)-1)*x0(1)/w_solution(1)*w_solution(j);
            end
        end
        kk=j;
        
        if j>length(colors)
            plot(lambda_n,y(:,j-n),'color',colors_m(j-length(colors),:),'linewidth',1)
        else
            plot(lambda_n,y(:,j-n),'color',colors(j,:),'linewidth',1)
        end
    else
        n=n+1;
    end
end

kk2=kk;

n=0;
if fit_region==3
    for j=1:length(center_22) % Fit S11 and S22
        if center_22(j)>0
            XX=sqrt(log(2))*(lambda_n-center_22(j)-x0(2*(j+kk2-n)))/(0.5*FWHM_22(j,2)*x0(variables));
            YY=sqrt(log(2))*FWHM_22(j,1)/FWHM_22(j,2);
            if dope==0
                y(:,j+kk2-n)=x0(2*(j+kk2-n)-1)*sqrt(log(2)/pi)/(0.5*FWHM_22(j,2)*x0(variables)).*real(complexErrorFunction(XX,YY));
                area_f(j+kk2)=x0(2*(j+kk2-n)-1);
            else
                y(:,j+kk2-n)=(1-x0(var))/w2_old*x0(2*(j+kk2-n)-1)*sqrt(log(2)/pi)/(0.5*FWHM_22(j,2)*x0(variables)).*real(complexErrorFunction(XX,YY));
                area_f(j+kk2)=(1-x0(var))/w2_old*x0(2*(j+kk2-n)-1);
            end
            kk=kk2+j;
            
            
            if j>length(colors)
                plot(lambda_n,y(:,j+kk2-n),'color',colors_m(j-length(colors),:),'linewidth',1)
            else
                plot(lambda_n,y(:,j+kk2-n),'color',colors(j,:),'linewidth',1)
            end
        else
            n=n+1;
        end
    end
end

kk2=kk;
ph=0;
legstr_ph=[];

for j=1:length(pp)
    y(:,kk2+j-n)=area_f(pp(j))*(0.017+0.1/diam_n(pp(j))+x0(variables+1))*exp(-2*(lambda_n-h*c/(1e-9*(h*c/((center(pp(j))+x0(pp(j)*2))*1e-9)+0.2))).^2/(FWHM_g(j)*x0(variables+2)/sqrt(log(4)))^2)/(FWHM_g(j)*x0(variables+2)*sqrt(pi/(4*log(2))));
    kk=kk2+j;
    ph=2;
    
    plot(lambda_n,y(:,kk2+j-n),'color',[colors(pp(j),1),.8,colors(pp(j),3)],'linewidth',1)
    
    legstr_ph{j}=['EPS-',legstr{pp(j)}];
end

kk2=kk;

for j=1:length(pp2)
    y(:,kk2+j-n)=area_f(pp2(j)+length(center))*(0.017+0.1/diam_n(pp2(j))+x0(variables+1+ph))*exp(-2*(lambda_n-h*c/(1e-9*(h*c/((center_22(pp2(j))+x0((pp2(j)+length(center))*2))*1e-9)+0.2))).^2/(FWHM_g_22(j)*x0(variables+2+ph)/sqrt(log(4)))^2)/(FWHM_g_22(j)*x0(variables+2+ph)*sqrt(pi/(2*log(4))));
    
    plot(lambda_n,y(:,kk2+j-n),'color',[colors(pp2(j),1),.8,colors(pp2(j),3)],'linewidth',1)
    legstr_ph{j+length(pp)}=['EPS-',legstr_22{pp2(j)}];
end

plot(lambda_n,Film,'k','linewidth',2)

plot(lambda_n,sum(y,2),'r','linewidth',2)

legstrs=[legstr legstr_22 legstr_ph];

legend([legstrs {'Measured Spectrum'} {'Calculated Spectrum'}])
title(['Voigtian Fit with: nSSE = ',num2str(sse)],'FontSize',26)

xlabel('Wavelength (nm)','FontSize',24)
ylabel('Absorption (a.u.)','FontSize',24)
set(gca,'LineWidth',2)
set(gca,'FontSize',22)


%..........................................................................


%--------------------------------------------------------------------------
% Section 4 - "check_fit"
%--------------------------------------------------------------------------


function [tt,max_s,start_s,conc_c,w_l]=check_fit(max_shift,start_shift,conc_change,wl)
global s_ud fit_region

% The user is presented with the fit and asked whether they are satisfied
% with its quality. If not, the user is given various options to improve
% the quality of the fit.

tt=input('Are you satisfied with the fit (y/n)? ','s');

max_s=max_shift;
start_s=start_shift;
conc_c=conc_change;
w_l=wl;

if tt=='n'
    if isempty(fit_region)==1
        opt=input('Please enter the number of the variable you want to adjust:\n(1) - Shift absorption data up or down in percentage of the maximum peak intensity (in %)\n(2) - Maximum possible red-shift of a nanotube transition (in nm)\n(3) - Initial red shift of each nanotube transition (in nm)\n(4) - Change the relative concentration change of (n,m) species in solution and film (in %)\n(5) - Change of wavelength range\n ','s');
    else
        opt=input('Please enter the number of the variable you want to adjust:\n(1) - Shift absorption data up or down in percentage of the maximum peak intensity (in %)\n(2) - Maximum possible red-shift of a nanotube transition (in nm)\n(3) - Initial red shift of each nanotube transition (in nm)\n(4) - Change the relative concentration change of (n,m) species in solution and film (in %)\n(5) - Change of wavelength range\n(6) - Change of fit region\n ','s');
    end
    
    opt=sscanf(opt,'%i');
    
    switch opt
        case 1
            s_ud=input('Please enter shift in percentage of the maximum value,\n e.g. a down shift of -1% or an up-shift of 3%: ','s');
            s_ud=sscanf(s_ud,'%i');
            s_ud=s_ud/100;
        case {2,3,4}
            cnames = {};
            rnames = {'Maximum Red-Shift (nm)','Initial Red-Shift (nm)','Relative Concentration Change (%)'};
            vals=[true true true];
            
            while 1>0
                f=figure('Position',[100 100 800 150]); t = uitable('Data', [max_shift;start_shift;conc_change],'ColumnName',cnames,'RowName',rnames,'ColumnEditable', vals);
                t.Position(3) = t.Extent(3);
                t.Position(4) = t.Extent(4);
                
                opt=input('Please change the variables in the table and press d if done: ','s');
                if opt=='d'
                    max_s=t.Data(1);
                    start_s=t.Data(2);
                    conc_c=t.Data(3);
                    close(f)
                    break
                end
            end
        case 5
            wl=input('Please enter the starting and ending wavelength for the fit.\nBe aware that the starting and ending wavelengths should be inside the wavelength regime defined during background subtraction:\n ','s');
            w_l=sscanf(wl,'%i');
        case 6
            fr=input('Would you like to fit S11 (1), S22 (2) or the entire region of your spectrum (3)?\n','s');
            fit_region=sscanf(fr,'%i');
            wl=input('Please enter the starting and ending wavelength of your fit regime.\nBe aware that the starting and ending wavelengths should be inside the wavelength regime defined during background subtraction:\n ','s');
            w_l=sscanf(wl,'%i');
    end
end


%..........................................................................


%--------------------------------------------------------------------------
% Section 5 - "save_data"
%--------------------------------------------------------------------------


function save_data(w,legstrs,x_fit,Cent,Cent_22,fwhm,fwhm_22,area_f,pp,pp2,FWHM_g,FWHM_g_22,diam_n,w_solution,h,c,fit_region,folder_name)

switch w
    case {1,2}
        if isempty(fit_region)==1 || fit_region~=3
            Results_Film{1,1}='CNTs';
            Results_Film{1,2}='Diameter (nm)';
            Results_Film{1,3}='Center (nm)';
            Results_Film{1,4}='FWHM (nm)';
            Results_Film{1,5}='FWHM Broadening';
            Results_Film{1,6}='Height / Area-Frac';
            Results_Film{1,7}='Area';
            Results_Film{1,8}='Rel. Concentration Solution (%)';
            Results_Film{1,9}='Rel. Concentration Film (%)';
            x_len=0;
            
            if isempty(fit_region)==0
                switch fit_region
                    case 1
                        center=Cent;
                        FWHM=fwhm;
                        phonon_pos=pp;
                        phonon_pos2=[];
                        variables=2*length(center);
                    case 2
                        center=Cent_22;
                        FWHM=fwhm_22;
                        FWHM_g=FWHM_g_22;
                        phonon_pos=pp2;
                        phonon_pos2=[];
                        variables=2*(length(center)-sum(center==0));
                end
                center_22=[];
            else
                center=Cent;
                center_22=[];
                FWHM=fwhm;
                phonon_pos=pp;
                phonon_pos2=[];
                variables=2*(length(center)-sum(center==0));
            end
            
        else
            Results_Film{1,1}='CNTs';
            Results_Film{1,2}='Diameter (nm)';
            Results_Film{1,3}='Center S11 (nm)';
            Results_Film{1,4}='FWHM S11 (nm)';
            Results_Film{1,5}='FWHM Broadening S11';
            Results_Film{1,6}='Height / Area-Frac S11';
            Results_Film{1,7}='Area S11';
            Results_Film{1,8}='Center (nm) S22';
            Results_Film{1,9}='FWHM (nm) S22';
            Results_Film{1,10}='FWHM Broadening S22';
            Results_Film{1,11}='Height / Area-Frac S22';
            Results_Film{1,12}='Area S22';
            Results_Film{1,13}='Rel. Concentration (%) S11 Solution';
            Results_Film{1,14}='Rel. Concentration (%) S22 Solution';
            Results_Film{1,15}='Rel. Concentration (%) S11 Film';
            Results_Film{1,16}='Rel. Concentration (%) S22 Film';
            
            center=Cent;
            center_22=Cent_22;
            FWHM=fwhm;
            FWHM_22=fwhm_22;
            
            if length(center_22)>length(center) % Metallic nanotubes were assigned to the S22 region
                len_s1=length(center);
                len_s2=length(center_22)-sum(center_22==0);
                len_m1=length(center_22)-length(center);
                center=[center;center_22(1,len_s1+1:end)'];
                center_22(len_s1+1:end)=[];
                FWHM=[FWHM;FWHM_22(1,len_s1+1:end)'];
                FWHM_22(len_s1+1:end)=[];
                x_fit=[x_fit(1:2*len_s1),x_fit(2*(len_s1+len_s2-len_m1)+1:2*(len_s1+len_s2)),x_fit(2*len_s1+1:2*(len_s1+len_s2-len_m1)),x_fit(2*(len_s1+len_s2)+1:end)];
                area_f=[area_f(1:len_s1);area_f(len_s1+len_s2-len_m1+1:len_s1+len_s2);area_f(len_s1+1:len_s1+len_s2-len_m1)];
                legstrs=[legstrs(1:len_s1) legstrs(len_s1+len_s2-len_m1+1:len_s1+len_s2) legstrs(len_s1+1:len_s1+len_s2-len_m1) legstrs(len_s1+len_s2+1:end)];
                w_solution=[w_solution(1:len_s1);w_solution(len_s1+length(center_22)+1:len_s1+length(center_22)+len_m1);w_solution(len_s1+1:len_s1+length(center_22))];
            end
            
            phonon_pos=pp;
            phonon_pos2=pp2;
            variables=2*(length(center)+length(center_22)-sum(center_22==0));
            x_len=length(center_22)-sum(center_22==0);
        end
        
        n=0;
        nn=0;
        for i=1:length(center)+length(phonon_pos)+length(phonon_pos2)
            if i<=length(center)
                if center(i)>0
                    Results_Film{i-nn+1,1}=legstrs{i-nn};
                    Results_Film{i-nn+1,2}=diam_n(i);
                    Results_Film{i-nn+1,3}=center(i)+x_fit((i-nn)*2);
                    Results_Film{i-nn+1,4}=FWHM(i)*x_fit(variables+1);
                    Results_Film{i-nn+1,5}=x_fit(variables+1);
                    if w==1 % Lorentzian
                        if i==1
                            Results_Film{i-nn+1,6}=2*x_fit(1)/(pi*Results_Film{i-nn+1,4});
                            Results_Film{i-nn+1,7}=x_fit(1);
                        else
                            Results_Film{i-nn+1,6}=2*x_fit(1)/w_solution(1)*w_solution(i)*x_fit(2*(i-nn)-1)/(pi*Results_Film{i-nn+1,4});
                            Results_Film{i-nn+1,7}=x_fit(1)/w_solution(1)*w_solution(i)*x_fit(2*(i-nn)-1);
                        end
                    else % Gaussian
                        if i==1
                            Results_Film{i-nn+1,6}=x_fit(1)/Results_Film{i-nn+1,4}*sqrt(2*log(4)/pi);
                            Results_Film{i-nn+1,7}=x_fit(1);
                        else
                            Results_Film{i-nn+1,6}=x_fit(1)/w_solution(1)*w_solution(i)*x_fit(2*(i-nn)-1)/Results_Film{i-nn+1,4}*sqrt(2*log(4)/pi);
                            Results_Film{i-nn+1,7}=x_fit(1)/w_solution(1)*w_solution(i)*x_fit(2*(i-nn)-1);
                        end
                    end
                else
                    nn=nn+1;
                end
                
                if fit_region==3
                    if i<=length(center_22)
                        if center_22(i)>0
                            
                            Results_Film{i+1,8}=center_22(i)+x_fit((i+length(center)-n)*2);
                            Results_Film{i+1,9}=FWHM_22(i)*x_fit(variables+1);
                            Results_Film{i+1,10}=x_fit(variables+1);
                            
                            if w==1 % Lorentzian
                                Results_Film{i+1,11}=2*x_fit(1)/w_solution(1)*w_solution(i+length(center))*x_fit(2*(i+length(center)-n)-1)/(pi*Results_Film{i+1,9});
                                Results_Film{i+1,12}=x_fit(1)/w_solution(1)*w_solution(i+length(center))*x_fit(2*(i+length(center)-n)-1);
                            else % Gaussian
                                Results_Film{i+1,11}=x_fit(1)/w_solution(1)*w_solution(i+length(center))*x_fit(2*(i+length(center)-n)-1)/Results_Film{i+1,9}*sqrt(2*log(4)/pi);
                                Results_Film{i+1,12}=x_fit(1)/w_solution(1)*w_solution(i+length(center))*x_fit(2*(i+length(center)-n)-1);
                            end
                            
                            Results_Film{i+1,14}=w_solution(i+length(center))/sum(w_solution(1+length(center):end))*100;
                            Results_Film{i+1,16}=area_f(i+length(center)-n)/sum(area_f((length(center)+1):end))*100;
                        else
                            n=n+1;
                        end
                    end
                    Results_Film{i+1,13}=w_solution(i)/sum(w_solution(1:length(center)))*100;
                    Results_Film{i+1,15}=area_f(i)/sum(area_f(1:length(center)))*100;
                else
                    if center(i)>0
                        Results_Film{i-nn+1,8}=w_solution(i)*100;
                        Results_Film{i-nn+1,9}=area_f(i-nn)/sum(area_f)*100;
                    end
                end
                
            else
                Results_Film{i-nn+1,1}=legstrs{i+x_len-nn};
            end
        end
        
        xx_start=2;
        
        for j=1:length(phonon_pos) % Only if there was a phonon sideband
            Results_Film{length(center)-nn+1+j,3}=h*c/(1e-9*(h*c/((center(phonon_pos(j))+x_fit((phonon_pos(j)*2)))*1e-9)+0.2));
            Results_Film{length(center)-nn+1+j,4}=FWHM_g(j)*x_fit(variables+xx_start+1);
            Results_Film{length(center)-nn+1+j,5}=x_fit(variables+xx_start+1);
            Results_Film{length(center)-nn+1+j,6}=x_fit(variables+xx_start);
            Results_Film{length(center)-nn+1+j,7}=area_f(phonon_pos(j))*(0.017+0.1/diam_n(phonon_pos(j))+x_fit(variables+xx_start));
            xx_start=xx_start+2;
        end
        for j=1:length(phonon_pos2) % Only if there was a phonon sideband
            Results_Film{length(center)+length(phonon_pos)+1+j,8}=h*c/(1e-9*(h*c/((center_22(phonon_pos2(j))+x_fit(((length(center)+phonon_pos2(j))*2)))*1e-9)+0.2));
            Results_Film{length(center)+length(phonon_pos)+1+j,9}=FWHM_g_22(j)*x_fit(variables+xx_start+1);
            Results_Film{length(center)+length(phonon_pos)+1+j,10}=x_fit(variables+xx_start+1);
            Results_Film{length(center)+length(phonon_pos)+1+j,11}=x_fit(variables+xx_start);
            Results_Film{length(center)+length(phonon_pos)+1+j,12}=area_f(phonon_pos2(j)+length(center))*(0.017+0.1/diam_n(phonon_pos2(j))+x_fit(variables+xx_start));
        end
        
    case 3
                
        if isempty(fit_region)==1 || fit_region~=3
            Results_Film{1,1}='CNTs';
            Results_Film{1,2}='Diameter (nm)';
            Results_Film{1,3}='Center (nm)';
            Results_Film{1,4}='FWHM_l (nm)';
            Results_Film{1,5}='FWHM_g (nm)';
            Results_Film{1,6}='FWHM_v (nm)';
            Results_Film{1,7}='FWHM Broadening';
            Results_Film{1,8}='Height / Area-Frac';
            Results_Film{1,9}='Area';
            Results_Film{1,10}='Rel. Concentration Solution (%)';
            Results_Film{1,11}='Rel. Concentration Film (%)';
            x_len=0;
            
            if isempty(fit_region)==0
                switch fit_region
                    case 1
                        center=Cent;
                        FWHM=fwhm;
                        phonon_pos=pp;
                        phonon_pos2=[];
                        variables=2*length(center);
                    case 2
                        center=Cent_22;
                        FWHM=fwhm_22;
                        FWHM_g=FWHM_g_22;
                        phonon_pos=pp2;
                        phonon_pos2=[];
                        variables=2*(length(center)-sum(center==0));
                end
                center_22=[];
            else
                center=Cent;
                center_22=[];
                FWHM=fwhm;
                phonon_pos=pp;
                phonon_pos2=[];
                variables=2*(length(center)-sum(center==0));
            end
            
        else
            Results_Film{1,1}='CNTs';
            Results_Film{1,2}='Diameter (nm)';
            Results_Film{1,3}='Center S11 (nm)';
            Results_Film{1,4}='FWHM_l S11 (nm)';
            Results_Film{1,5}='FWHM_g S11 (nm)';
            Results_Film{1,6}='FWHM_v S11 (nm)';
            Results_Film{1,7}='FWHM Broadening S11';
            Results_Film{1,8}='Height / Area-Frac S11';
            Results_Film{1,9}='Area S11';
            Results_Film{1,10}='Center S22 (nm)';
            Results_Film{1,11}='FWHM_l  S22(nm)';
            Results_Film{1,12}='FWHM_g S22 (nm)';
            Results_Film{1,13}='FWHM_v S22 (nm)';
            Results_Film{1,14}='FWHM Broadening S22';
            Results_Film{1,15}='Height / Area-Frac S22';
            Results_Film{1,16}='Area S22';
            Results_Film{1,17}='Rel. Concentration S11 Solution (%)';
            Results_Film{1,18}='Rel. Concentration S22 Solution (%)';
            Results_Film{1,19}='Rel. Concentration S11 Film (%)';
            Results_Film{1,20}='Rel. Concentration S22 Film (%)';
            
            center=Cent;
            center_22=Cent_22;
            FWHM=fwhm;
            FWHM_22=fwhm_22;
            
            if length(center_22)>length(center) % Metallic nanotubes were assigned to the S22 region
                len_s1=length(center);
                len_s2=length(center_22)-sum(center_22==0);
                len_m1=length(center_22)-length(center);
                center=[center;center_22(1,len_s1+1:end)'];
                center_22(len_s1+1:end)=[];
                FWHM=[FWHM;FWHM_22(len_s1+1:end,1:3)];
                FWHM_22=FWHM_22(1:len_s1,1:3);
                x_fit=[x_fit(1:2*len_s1),x_fit(2*(len_s1+len_s2-len_m1)+1:2*(len_s1+len_s2)),x_fit(2*len_s1+1:2*(len_s1+len_s2-len_m1)),x_fit(2*(len_s1+len_s2)+1:end)];
                area_f=[area_f(1:len_s1);area_f(len_s1+length(center)-len_m1+1:len_s1+length(center));area_f(len_s1+1:len_s1+length(center)-len_m1)];
                legstrs=[legstrs(1:len_s1) legstrs(len_s1+len_s2-len_m1+1:len_s1+len_s2) legstrs(len_s1+1:len_s1+len_s2-len_m1) legstrs(len_s1+len_s2+1:end)];
                w_solution=[w_solution(1:len_s1);w_solution(len_s1+length(center_22)+1:len_s1+length(center_22)+len_m1);w_solution(len_s1+1:len_s1+length(center_22))];
            end
            
            phonon_pos=pp;
            phonon_pos2=pp2;
            variables=2*(length(center)+length(center_22)-sum(center_22==0));
            x_len=length(center_22)-sum(center_22==0);
        end
        
        n=0;
        nn=0;
        for i=1:length(center)+length(phonon_pos)+length(phonon_pos2)
            if i<=length(center)
                if center(i)>0
                    Results_Film{i-nn+1,1}=legstrs{i-nn};
                    Results_Film{i-nn+1,2}=diam_n(i);
                    Results_Film{i-nn+1,3}=center(i)+x_fit(2*(i-nn));
                    Results_Film{i-nn+1,4}=FWHM(i,3)*x_fit(variables+1)/(.5436+sqrt(0.2166+(FWHM(i,2)/(FWHM(i,1)))^2));
                    Results_Film{i-nn+1,5}=Results_Film{i-nn+1,4}*FWHM(i,2)/(FWHM(i,1));
                    Results_Film{i-nn+1,6}=FWHM(i,3)*x_fit(variables+1);
                    Results_Film{i-nn+1,7}=x_fit(variables+1);
                    Results_Film{i-nn+1,8}=area_f(i)*real(complexErrorFunction(0,sqrt(log(2))*Results_Film{i-nn+1,4}/Results_Film{i+1,5}))*sqrt(log(2)/pi)/(0.5*Results_Film{i-nn+1,5});
                    Results_Film{i-nn+1,9}=area_f(i);
                else
                    nn=nn+1;
                end
                
                if fit_region==3
                    if i<=length(center_22)
                        if center_22(i)>0
                            
                            Results_Film{i+1,10}=center_22(i)+x_fit((i+length(center)-n)*2);
                            Results_Film{i+1,11}=FWHM_22(i,3)*x_fit(variables+1)/(.5436+sqrt(0.2166+(FWHM_22(i,2)/(FWHM_22(i,1)))^2));
                            Results_Film{i+1,12}=Results_Film{i+1,12}*FWHM_22(i,2)/(FWHM_22(i,1));
                            Results_Film{i+1,13}=FWHM_22(i,3)*x_fit(variables+1);
                            Results_Film{i+1,14}=x_fit(variables+1);
                            Results_Film{i+1,15}=area_f(i)*real(complexErrorFunction(0,sqrt(log(2))*Results_Film{i+1,11}/Results_Film{i+1,12}))*sqrt(log(2)/pi)/(0.5*Results_Film{i+1,12});
                            Results_Film{i+1,16}=area_f(i+length(center));
                            
                            Results_Film{i+1,18}=w_solution(i+length(center))/sum(w_solution(1+length(center):end))*100;
                            Results_Film{i+1,20}=area_f(i+length(center))/sum(area_f((length(center)+1):end))*100;
                        else
                            n=n+1;
                        end
                    end
                    
                    Results_Film{i+1,17}=w_solution(i)/sum(w_solution(1:length(center)))*100;
                    Results_Film{i+1,19}=area_f(i)/sum(area_f(1:length(center)))*100;
                else
                    if center(i)>0
                        Results_Film{i-nn+1,10}=w_solution(i)*100;
                        Results_Film{i-nn+1,11}=area_f(i)/sum(area_f)*100;
                    end
                end
                
            else
                Results_Film{i-nn+1,1}=legstrs{i+x_len-nn};
            end
        end
        xx_start=2;
        for j=1:length(phonon_pos) % Only if there was a phonon sideband
            Results_Film{length(center)-nn+1+j,3}=h*c/(1e-9*(h*c/((center(phonon_pos(j))+x_fit(2*phonon_pos(j)))*1e-9)+0.2));
            Results_Film{length(center)-nn+1+j,5}=FWHM_g(j)*x_fit(variables+3);
            Results_Film{length(center)-nn+1+j,7}=x_fit(variables+3);
            Results_Film{length(center)-nn+1+j,8}=x_fit(variables+2);
            Results_Film{length(center)-nn+1+j,9}=area_f(phonon_pos(j))*(0.017+0.1/diam_n(phonon_pos(j))+x_fit(variables+2));
            xx_start=xx_start+2;
        end
        
        for j=1:length(phonon_pos2) % Only if there was a phonon sideband
            Results_Film{length(center)+length(phonon_pos)+1+j,10}=h*c/(1e-9*(h*c/((center_22(phonon_pos2(j))+x_fit(((length(center)+phonon_pos2(j))*2)))*1e-9)+0.2));
            Results_Film{length(center)+length(phonon_pos)+1+j,12}=FWHM_g_22(j)*x_fit(variables+xx_start+1);
            Results_Film{length(center)+length(phonon_pos)+1+j,14}=x_fit(variables+xx_start+1);
            Results_Film{length(center)+length(phonon_pos)+1+j,15}=x_fit(variables+xx_start);
            Results_Film{length(center)+length(phonon_pos)+1+j,16}=area_f(phonon_pos2(j)+length(center))*(0.017+0.1/diam_n(phonon_pos2(j))+x_fit(variables+xx_start));
        end
end

oldD=cd(folder_name);
save('Results_Film','Results_Film')
cd(oldD)


%..........................................................................


function export_data(legstr,lambda_n,Film_n,y,folder_name)

opt=input('Do you want to export the fitted data as a .txt file (y/n)? ','s');

if opt=='y'
    
    oldD=cd(folder_name);
    
    legstr=[{'Wavelength' 'Measurement Data' 'Calculated Data'},legstr];
    
    if size(lambda_n,1)<size(lambda_n,2)
        lambda_n=lambda_n';
        Film_n=Film_n';
    end
    
    y=[lambda_n,Film_n,sum(y,2),y];
    
    str=[];
    fileID = fopen('Individual_fits_Film.txt','w');
    for i=1:length(legstr)
        fprintf(fileID,['%',num2str(6*i),'s '],legstr{1,i});
        str=[str,'%',num2str(6*i),'f '];
    end
    fprintf(fileID,'\n');
    fprintf(fileID,[str,'\n'],y);
    fclose(fileID);
    
    cd(oldD);
end