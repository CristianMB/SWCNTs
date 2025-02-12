%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% Function to fit single walled carbon nanotube solution absorption spectra
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
%
% List of adjustable variables:
%
% PEAK FIT
%
% upper_lim = Upper boundary limit for peak shift during the fit.
%
% upper_lim_peak = If an (n,m) species was assigned to a peak in the
%                 sub-function "check_peaks" the upper boundary limit
%                 should be smaller than for those species that were not
%                 detected.
%
% lower_lim = Lower boundary limit for peak shift during the fit. This
%             value will be added to 0. E.g. lower_lim = -2 -> Lower
%             boundary limit in shifting the peak position during the fit
%             = -2.
%
% lower_lim_peak = Lower boundary limit for the shift of peak positions
%                  that were assigned in the sub-function "check_peaks".
%                  This value will be added to 0.
%
% broad_down = Lower limit for peak broadening during the fit in percentage
%              of the initial full width at half maximum (FWHM).
%
% broad_up = Upper limit for peak broadening during the fit.
%
% broad_up_turnover = Upper limit for peak broadening during the fit for
%                     peak positions greater than a certain "turnover"
%                     value. The initial broadening of FWHM might be
%                     insufficient beyond a certain wavelength.
%
% turnover = Wavelength above which the FWHM is allowed to broaden up to
%            broad_up_turnover, in nm.
%
% height_guess = It is possible to obtain the best fit with two different
%                approaches and this can be defined by the value assigned
%                to "height_guess". Approach (1): all SWCNTs are assigned
%                with 90 % of their initial intensity, taken from the
%                background corrected data. Approach (2): each peak that
%                was assigned is subtracted from the background corrected
%                data to yield a "new" absorption profile. This new profile
%                is then used to determine the starting height of the
%                peaks. As outlined in the paper the fit can be performed
%                as a combination of Approach (1) and (2) by defining an
%                approriate turnover wavelength.
%
% height_init_low = Lower limit of most intense peak; given as percentage
%                   of peak height.
%
% height_init = Starting value of most intense peak; given as percentage of
%               peak height.
%
% height_init_up = Upper limit of most intense peak; given as percentage of
%                  peak height.
%
% height_low_hg2 = Lower limit of peak height for heigt assignment method
%                  2.
%
% height_start_hg2 = Starting value for estimating peak height for height
%                    assignment method 2.
%
% height_up_hg2 = Upper limit of peak height for heigt assignment method 2.
%
% height_low_hg1 = Lower limit of peak height for heigt assignment method
%                  1.
%
% height_start_hg1 = Starting value for estimating peak height for height
%                    assignment method 1.
%
% height_up_hg1 = Upper limit of peak height for heigt assignment method 1.
%
% lorentz_up = Maximum HWHM of Lorentzian part of the Voigtian line profile
%              (in nm).
%
% ratio = Initial guess of ratio of Gaussian to Lorentzian HWHM for Voigt
%         fit.
%
% ratio_low = Lower limit of ratio of Gaussian to Lorentzian HWHM for Voigt
%             fit.
%
% ratio_up = Upper limit of ratio of Gaussian to Lorentzian HWHM for Voigt
%            fit.
%
% ratio_change = Change of the Gaussian / Lorentzian ratio for Voigtian fit
%                as percentage of initial ratio.
%
% max_shift = Maximal shift that is considered to be reasonable for the
%             red-shift of a nanotube transition compared to reference
%             data.
%
% path_length = Pathlength of light through the cuvette during absorption
%               measurement. Given in cm.
%
% EXCITON PHONON SIDEBAND
%
% fwhm_gauss = Intitial guess of full width at half maximum of the exction
%              phonon sideband (EPS) that is modelled to be a Gaussian line
%              shape.
%
% area_part = Diameter independent factor, that accounts for a greater or
%             smaller spectral weight transfer from the nanotube peak to
%             its associated exciton phonon sideband as proposed by
%             Perebeinos et al. (PRL, 2005 94, 027402).
%
% area_part_low = Lower boundary limit for area_part during the fit.
%
% area_part_up = Upper boundary limit for area_part during the fit.
%
% fwhm_gauss_b_up = Upper limit of broadening of Gaussian FWHM of the
%                   exciton phonon sideband.
%
% fwhm_gauss_b_down = Lower limit of broadening of Gaussian FWHM of the
%                     exciton phonon sideband.
%
% eps_shift_down = Lower limit of the shift for the exciton phonon sideband
%                  away from 0.2 eV (in eV).
%
% eps_shift_start = Initial guess of the shift for the exciton phonon
%                   sideband away from 0.2 eV (in eV).
%
% eps_shift_up = Upper limit of the shift for the exciton phonon sideband
%                away from 0.2 eV (in eV).
%
% ph2_broad_d = Lower limit of the change in FWHM for the first exciton
%               phonon sideband.
%
% ph2_broad_s = Initial guess of the change in FWHM for the first exciton
%               phonon sideband. If the value is equal to one, the FWHM of
%               following phonon sidebands will be the same as the very
%               first.
%
% ph2_broad_up = Upper limit of the change in FWHM for the first exciton
%                phonon sideband.
%
% FIT OF ENTIRE REGION
%
% height_ratio = Ratio of S11 / S22 peak intensity.
%
% height_ratio_low = Lower limit of height_ratio during the fit.
%
% height_ratio_up = Upper limit of height_ratio during the fit.
%
% height_change = Factor that is multiplied with height_ratio to keep the
%                 intensity ratio of S11 / S22 similar.
%
% height_change_low = Lower limit of height_change during fit.
%
% height_change_up = Upper limit of height_change during fit.
%
% FIT OF S22 AND M11
%
% height_low = Lower limit of height that was calculated for S22 based on a
%              previous fit of the entire region.
%
% height_up = Upper limit of height that was calculated for S22 based on a
%             previous fit of the entire region.
%
% height_start = Initial percentage of height that was calculated for S22
%                based on a previous fit of the entire region.
%
% up_shift = Upper limit of the shift of the center position that was
%            calculated for S22 based on a previous fit of the entire
%            region.
%
% low_shift = Lower limit of the shift of the center position that was
%             calculated for S22 based on a previous fit of the entire
%             region.
%
% fwhm_broad_down = Lower limit of the FWHM of S22 (n,m) species from
%                   previous fit of entire region.
%
% fwhm_broad_up = Upper limit of the FWHM of S22 (n,m) species from
%                   previous fit of entire region.
%
% fwhm_broad_start = Initial guess of the FWHM of S22 (n,m) species from
%                   previous fit of entire region.
%
% f1_change_down = Lower limit of the factor f1, needed for calculation of
%                  exction phonon sideband, based on previous fit of entire
%                  region.
%
% f1_change_up = Upper limit of the factor f1, needed for calculation of
%                exction phonon sideband, based on previous fit of entire
%                region.
%
% f1_change_start = Initial guess of the factor f1, needed for calculation
%                   of exction phonon sideband, based on previous fit of
%                   entire region.
%
% fwhm_gauss_b_up_22 = Upper limit of the broadening of the Gaussian FWHM
%                      for the exciton phonon sideband based on a previous
%                      fit of the entire region.
%
% fwhm_gauss_b_down_22 = Lower limit of the broadening of the Gaussian FWHM
%                      for the exciton phonon sideband based on a previous
%                      fit of the entire region.
%
% fwhm_gauss_b_start_22 = Initial guess of the broadening of the Gaussian
%                         FWHM for the exciton phonon sideband based on a
%                         previous fit of the entire region.
%
% ratio_up = Upper limit for the Gaussian / Lorentzian ratio for Voigtian
%            fit of S22 during assignment of metallic nanotubes.
%
% ratio_start = Initial guess for the Gaussian / Lorentzian ratio for
%               Voigtian fit of S22 during assignment of metallic
%               nanotubes.
%
% ratio_low = Lower limit for the Gaussian / Lorentzian ratio for Voigtian
%            fit of S22 during assignment of metallic nanotubes.
%
% doping = Consider doping effects: doping reduces S11 intensity. Therefore
%          S22 might become more intense than S11.
%
% GENERAL CONSTANTS
%
% h = Planck's constant in eV/s.
%
% c = Speed of light in m/s.
%
% c_len = Carbon-carbon bond length in nm.
%--------------------------------------------------------------------------


function Solution_Fitting
global Absorption lambda center h c fwhm_nm_v fwhm_nm_v22 cc Order diam_n fwhm_gauss center_22 cc2 s_ud


%--------------------------------------------------------------------------
% Definition of initial parameters required for the fit.
%--------------------------------------------------------------------------

% Peak fit

upper_lim=20;
upper_lim_peak=5;
lower_lim=0;
lower_lim_peak=-5;
broad_down = 0.8;
broad_up = 1.3;
broad_up_turnover=1.6;
turnover=1400;
height_guess=1;
height_init_low=0.8;
height_init=0.95;
height_init_up=1;
height_low_hg2 = 0.1;
height_start_hg2 = 1;
height_up_hg2 = 1.2;
height_low_hg1 = 0.1;
height_start_hg1 = .9;
height_up_hg1 = 1;
lorentz_up = 40;
ratio=1;
ratio_low=0.1;
ratio_up=2;
ratio_change=0.2;
max_shift=40;
path_length=0.2;

data_peakfit=[broad_down 1 broad_up;lower_lim_peak 0 upper_lim_peak; lower_lim 0 upper_lim; height_init_low height_init height_init_up; height_low_hg1 height_start_hg1 height_up_hg1; height_low_hg2 height_start_hg2 height_up_hg2];
Turnover=[turnover; broad_up_turnover; broad_down; height_guess];
data_voigt=[lorentz_up; ratio; ratio_low; ratio_up; ratio_change];

% Exciton phonon sideband

fwhm_gauss=40;
area_part=0;
area_part_low=-0.07;
area_part_up=0.07;
fwhm_gauss_b_up=2;
fwhm_gauss_b_down=0.5;
eps_shift_down=-0.005;
eps_shift_start=0;
eps_shift_up=0.005;
ph2_broad_d = 0.9;
ph2_broad_s = 1;
ph2_broad_up = 1.1;

data_EPS_peakfit=[fwhm_gauss_b_down fwhm_gauss fwhm_gauss_b_up;area_part_low area_part area_part_up;eps_shift_down eps_shift_start eps_shift_up; ph2_broad_d ph2_broad_s ph2_broad_up];

% Fit of entire region

height_ratio=4;
height_ratio_low=1;
height_ratio_up=5;
height_change_low=0.8;
height_change_up=1.2;
height_change=1;

data_heightratio=[height_ratio_low height_ratio height_ratio_up; height_change_low height_change height_change_up];

doping=0; % If doping==0, no doping is assumed and intensity and area of S22<S11. If doping==1, doping is assumed and intensity and area of S22 can become larger than their S11 contribution.

% Fit S22 and M11 based on previous fit of entire region

height_low=0.85;
height_up=1.05;
height_start=.9;
up_shift=5;
low_shift=-5;
fwhm_broad_down=0.75;
fwhm_broad_start=.8;
fwhm_broad_up=1.05;
f1_change_down=0.95;
f1_change_start=1;
f1_change_up=1.05;
fwhm_gauss_b_up_22=1.05;
fwhm_gauss_b_start_22=1;
fwhm_gauss_b_down_22=0.95;
ratio_up=1.05;
ratio_start=1;
ratio_low=0.95;

data_S22 = [fwhm_broad_down fwhm_broad_start fwhm_broad_up;low_shift 0 up_shift;height_low height_start height_up];
data_voigt_s22 = [ratio_low ratio_start ratio_up];
data_EPS = [fwhm_gauss_b_down_22 fwhm_gauss_b_start_22 fwhm_gauss_b_up_22; eps_shift_down eps_shift_start eps_shift_up;f1_change_down f1_change_start f1_change_up];
data_M11 = [0.6 1 broad_up; low_shift 0 upper_lim; 0.1 1.2 1.6];

% General constants

h=4.135667662e-15;
c=299792458;
c_len=0.142;


%--------------------------------------------------------------------------
% Sub-function "get_data" located in Section 1 performs:
% - user selection of absorption data
% - user defined background subtraction
% - user selection of reference data table
% - user selection of wavelength range to be fitted (S11 or S22)
% - initial estimation of the FWHM for the nanotubes provided in the data
%   table
% - initial assignment of peak intensities
%--------------------------------------------------------------------------


[lambda,Absorption,fwhm_nm,fwhm_nm_22,S11_sort,S22_sort,Names,n,m,diam,pos,height_init,colors,range,fit_region,start,End,protocol,time_stamp,fid,parallel_fit]=get_data(c_len,h,c);


%--------------------------------------------------------------------------
% Sub-function "select_nanotubes" located in Section 2 performs:
% - user selection of (n,m) species within the wavelength range of interest
%   to be fitted
% - modification of vectors defined in "get_data"
%
%
% Sub-function "check_S22" located in Section 3 performs:
% - user selection of a previous fit of the entire region
%
%
% Sub-function "select_nanotubes_metallic" located in Section 4 performs:
% - user selection of metallic (n,m) species within the wavelength range of
%   interest to be fitted
%--------------------------------------------------------------------------


switch fit_region
    case 1 % S11 region
        [vec,diam_n,colors_n,S11n,S22n,CNT_n,fwhm_nm_n,fwhm_nm_n22]=select_nanotubes(lambda(start:End),pos,height_init,S11_sort,S22_sort,colors,Names,diam,area_part,fwhm_gauss,Absorption(start:End),h,c,fwhm_nm,fwhm_nm_22,n,m,fit_region,protocol,fid);
        opt='n';
    case 2 % S22 region
        opt=input('Do you want to fit S22 based on a previous fit of the entire region and include metallic SWCNTs (y/n)? ','s');
        
        if opt=='n'
            [vec,diam_n,colors_n,S11n,S22n,CNT_n,fwhm_nm_n,fwhm_nm_n22]=select_nanotubes(lambda(start:End),pos,height_init,S11_sort,S22_sort,colors,Names,diam,area_part,fwhm_gauss,Absorption(start:End),h,c,fwhm_nm,fwhm_nm_22,n,m,fit_region,protocol,fid);
        else
            [height_22,center_22,FWHM_22,legstr_22,FWHM_g_22,f1_22,colors_22,w,y_sum,diam_22,phonon_pos2,fwhm_11,height_11,area_11,folder_name,legstr_help]=check_S22(lambda(start:End),h,c);
            [diam_m,colors_m,M11n,fwhm_nm_m11_n,CNT_n,height_m]=select_nanotubes_metallic(c_len,Absorption(start:End),lambda(start:End),h,c,y_sum,fid,protocol);
        end
end


%--------------------------------------------------------------------------
% Sub-function "check_peaks" located in Section 5 performs:
% - based on the total number N of (n,m) species defined in
%   "select_nanotubes" up to N number of peaks are found in the measurement
%   data
% - peaks not assigned to nanotube S11 or S22 transitions can be removed
% - additional peaks not automatically detected can be assigned
% - assignment of peaks in the user provided measurement data to SWCNTs in
%   the reference data table
% - addition of phonon sidebands where appropriate
% - modification of vectors defined in "select_nanotubes" based on peak
%   intensities
%--------------------------------------------------------------------------


switch fit_region
    case 1 % S11 region
        [Order,l_g,center,center_22,height2,up_bound,low_bound]=check_peaks(max_shift,Absorption(start:End),lambda(start:End),Names,vec,CNT_n,range(1),S11n,S22n,fit_region,upper_lim,upper_lim_peak,lower_lim,lower_lim_peak,fid,protocol);
    case 2 % S22 region
        if opt=='n'
            [Order,l_g,center,center_22,height2,up_bound,low_bound]=check_peaks(max_shift,Absorption(start:End),lambda(start:End),Names,vec,CNT_n,range(1),S11n,S22n,fit_region,upper_lim,upper_lim_peak,lower_lim,lower_lim_peak,fid,protocol);
        end
end


%--------------------------------------------------------------------------
% Sub-function "fit_data" located in Section 6 performs:
% - fit of absorption data with user specified (n,m) species and
%   line-shape
% - calculation of chi-square value to evaluate the quality of the fit
%
%
% Sub-function "fit_data_m" located in Section 7 performs:
% - fit of absorption data with user specified metallic (n,m) species and
%   previously defined line-shape
% - calculation of chi-square value to evaluate the quality of the fit
%
%
% Sub-function "plot_data" located in Section 8 performs:
% - plot the results from "fit_data"
%
%
% Sub-function "plot_data_m" located in Section 9 performs:
% - plot the results from "fit_data_m"
%
%
% Sub-function "check-fit" located in Section 10 performs:
% - user is asked to check whether they are satisfied with the quality of
%   the fit
%
%
% Sub-function "check-fit_m" located in Section 11 performs:
% - user is asked to check whether they are satisfied with the quality of
%   the fit of the S22 region with included metallic (n,m) species
%--------------------------------------------------------------------------


nn=1;
complete='n';
s_ud=0; % Variable that can be defined in "check_fit" to adjust the height of the absorption spectrum
PeakAssignment=[0; 0]; % If PeakAssignment(1 or 2) == 1, then redo "check_peaks" for S11 or S22 or both
c_er=0; % Count how often the internal loop for the fit of the entire region was entered. If it is the first time, do "check_peaks" for the S22 region, if not, evaluate "PeakAssignment".

% Pre-definition of variables for the fit of the entire region

if opt=='n'
    Order2=[];
    l_g2=[];
    up_bound2=[];
    low_bound2=[];
    w=[];
    height2_22=[];
    
    if fit_region~=2
        center_22=[];
    end
else
    metallic_peakassignment=0;
end

while 1>0
    
    if opt=='y'
        [x_fit,sse,colors_m,CNT_n,diam_m,M11n,fwhm_nm_m11_n]=fit_data_m(height_22+s_ud,center_22,FWHM_22,FWHM_g_22,f1_22,w,diam_22,M11n,fwhm_nm_m11_n,height_m+s_ud,lambda(start:End),Absorption(start:End)+s_ud,data_S22,phonon_pos2,h,c,data_EPS,data_M11,data_voigt_s22,data_voigt,metallic_peakassignment,colors_m,CNT_n,c_len,diam_m,y_sum,fid,protocol,fwhm_11,height_11,doping,parallel_fit);
        
        [plot_legstr,y]=plot_data_m(lambda(start:End),Absorption(start:End)+s_ud,colors_22,colors_m,center_22,M11n,phonon_pos2,x_fit,FWHM_22,fwhm_nm_m11_n,legstr_22,w,CNT_n,diam_22,h,c,FWHM_g_22,sse);
        
        [tt,data_S22,data_EPS,data_M11,data_voigt_s22,data_voigt,metallic_peakassignment,doping]=check_fit_m(data_S22,data_EPS,data_M11,Absorption(start:End),data_voigt_s22,data_voigt,w,fid,protocol,doping);
    else
        [x_fit,sse,cc,cc2,w,s1,e1,PeakAssignment,c_er,Order2,l_g2,height2_22,up_bound2,low_bound2,Order,l_g,center,height2,up_bound,low_bound]=fit_data(max_shift,data_voigt,data_peakfit,data_EPS_peakfit,data_heightratio,Turnover,PeakAssignment,fit_region,lambda,Absorption,height2+s_ud,Order,low_bound,up_bound,Names,vec,center,fwhm_nm_n,diam_n,l_g,fwhm_nm_n22,complete,start,End,h,c,S22n,S11n,CNT_n,upper_lim,upper_lim_peak,lower_lim,lower_lim_peak,c_er,Order2,l_g2,height2_22,up_bound2,low_bound2,w,fid,protocol,doping,parallel_fit);
        
        [legstr,plot_legstr,y]=plot_data(x_fit,colors_n,lambda,Absorption,Order,fwhm_nm_n,fwhm_gauss,diam_n,Names,vec,center,sse,cc,fwhm_nm_n22,cc2,fit_region,w,fwhm_nm_v22,fwhm_nm_v,complete,start,End);
        
        [nn,complete,tt,data_peakfit,data_EPS_peakfit,data_heightratio,Turnover,PeakAssignment,data_voigt,w,doping]=check_fit(data_peakfit,data_EPS_peakfit,data_heightratio,Turnover,nn,fit_region,complete,PeakAssignment,data_voigt,w,fid,protocol,doping);
    end
    
    if tt=='y'
        break
    end
end


%--------------------------------------------------------------------------
% Sub-function "save_variables_film" located in Section 12 performs:
% - save variables needed for a film fit in a new folder
%--------------------------------------------------------------------------


if opt=='n'
    [fwhm_nm_n,fwhm_nm_n22,diam_n,time_stamp]=save_variables_film(x_fit,center,sse,fwhm_nm_n,fwhm_nm_n22,diam_n,Order,plot_legstr,colors_n(Order,:),w,center_22,cc,cc2,fit_region,time_stamp,protocol);
end


%--------------------------------------------------------------------------
% Sub-function "save_data" located in Section 13 performs:
% - save data to a cell "Results" and save it as a MATLAB variable
%
%
% Sub-function "save_data_m" located in Section 14 performs:
% - save data to a cell "Results" and save it as a MATLAB variable
%--------------------------------------------------------------------------


if opt=='n'
    save_data(x_fit,w,fit_region,complete,center,cc,fwhm_nm_n,cc2,center_22,fwhm_nm_n22,h,c,colors_n,Order,diam_n,legstr,plot_legstr,nn,fwhm_gauss,time_stamp,path_length);
else
    save_data_m(area_11,M11n,center_22,diam_m,x_fit,fwhm_nm_m11_n,h,c,diam_22,FWHM_22,plot_legstr,phonon_pos2,FWHM_g_22,w,folder_name,legstr_help,protocol,time_stamp)
end


%--------------------------------------------------------------------------
% Sub-function "export_data" located in Section 15 performs:
% - save wavelength, absorption spectrum and the fit of each (n,m) species
%   in a .txt file.
%--------------------------------------------------------------------------


if opt=='n'
    export_data(plot_legstr,lambda(s1:e1),Absorption(s1:e1),y,time_stamp,[])
else
    export_data(plot_legstr,lambda(start:End),Absorption(start:End),y,'S22andM11',folder_name)
end


%--------------------------------------------------------------------------
%
%............................Sub-Functions.................................
%
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Section 1 - "get_data"
%--------------------------------------------------------------------------


function [lambda,Absorption,fwhm_nm,fwhm_nm_22,S11_sort,S22_sort,Names_n,n,m,diam,pos,height_init,colors,range,fit_region,start,End,protocol,time_stamp,fid,parallel_fit] = get_data(c_len,h,c)

% Check whether the Matlab version used is capable of performing parallel
% calculations. If yes, open a parallel pool that consists of of the max.
% number of cores minus one.
% Ask user whether a protocol file of the inputs shall be created or not.
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

protocol=input('Do you want to create a protocol file of your fit (y/n)?: ','s');
protocol=sscanf(protocol,'%s');

if exist('test.txt','file')~=0
    delete('test.txt');
end

disp('Select the file containing the absorption data');
filename = uigetfile('*.txt','Select the file containing the absorption data');

string='%f32 %f32';

fid=fopen(filename);
D = textscan(fid,string,'Delimiter',',','HeaderLines',2);
fclose(fid);

Data(:,1)=D{1};
Data(:,2)=D{2};

if Data(1,1)>Data(2,1)
    Data=flipud(Data);
end

figure('units','normalized','outerposition',[0 0 1 1]), plot(Data(:,1),Data(:,2),'k','linewidth',2)
xlabel('Wavelength (nm)','FontSize',33)
ylabel('Absorption (a.u.)','FontSize',33)
title('Measurement Data','FontSize',36)
set(gca,'LineWidth',2)
set(gca,'FontSize',30)

% In all cases the user must define a wavelength range of interest.
% Attention must be paid in the case that first the S11 and then the
% S22 region is fitted to ensure that all S22 transitions for (n,m) species
% selected for S11 actually are within the range defined for the background
% correction.
% Prompt the user to select a method of background subtraction. Method (1) 
% is adopted from Tian et al. (RSC Adv. 2015, 5, 102974). Method (2)
% is adopted from Nair et al. (Anal. Chem., 2006, 78, 7689-7696). Method (3)
% is based on the work of Naumov et al. (ACS Nano, 2011, 5, 3, 1639–1648).
% Method (4) performs no background subtraction. For method 1 to 3, the
% background subtraction is performed and constrained by functions defined
% in Section 1.1.

while 1>0
    q=input('Enter the starting and ending wavelengths for the region of interest, e.g. 400 1400. \nBe aware that the S11 and S22 transitions to be fitted should be in the region defined: ','s');
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
    if q(1)>q(2)
        warning('Please enter starting wavelength that is smaller than the ending wavelength.')
        ck=1;
    end
    if ck==0
        break
    end
end

rr=1;
bkg=1;

% Fano + Lorentzian Background
advanced=0;

while 1>0
    
    if bkg==1
        w=input('Please enter number of background subtraction:\n (1) - Fano + Lorentzian\n (2) - k/lambda^b\n (3) - A*exp(-b*lambda)\n (4) - No background\n (5) - Load Background\n','s');
        w=sscanf(w,'%i');
        bkg=2;
    end
    
    switch w
        case 1
            if rr==2
                while 1>0
                    q=input('Enter the starting and ending wavelengths for the region of interest, e.g. 400 1400. \nBe aware that the S11 and S22 transitions to be fitted should be in the region defined: ','s');
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
                    if q(1)>q(2)
                        warning('Please enter starting wavelength that is smaller than the ending wavelength.')
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
                    lb=[-4,4.2,0.4,max(yy(1:200))/16,4.8,1.5,0];
                    ub=[-2,4.7,0.8,max(yy(1:200))/9,6,3.5,max(yy(1:200))]; % define upper and lower boundary conditions
                    x_start=[-3,4.4,0.6,0.1,4.8,1,0.02];
                    
                    data_FL(:,1)=lb;
                    data_FL(:,2)=x_start;
                    data_FL(:,3)=ub;
                end
            end
            
            E=h*c./(xx(start:End,1)*1e-9); % transform nm in eV
            
            dlmwrite('test.txt',[E yy(start:End,1)]); % save data in test.txt so that the fit functions can access it easily
            
            if parallel_fit=='y'
                if rel>=2014
                    options = optimoptions('fmincon','Algorithm','interior-point','Display','off','UseParallel',true);
                else
                    options = optimoptions('fmincon','Algorithm','interior-point','Display','off','UseParallel','always');
                end
            else
                options = optimoptions('fmincon','Algorithm','interior-point','Display','off');
            end
            
            if parallel_fit=='y'
                try
                    options_l = optimoptions(@lsqnonlin,'Display','off','UseParallel',true);
                catch
                    options_l = optimoptions(@lsqnonlin,'Display','off');
                end
            else
                options_l = optimoptions(@lsqnonlin,'Display','off');
            end
            
            F_fano=@(x) double(yy(start:End,1)-x(4)*(x(1)+(E-x(2))/(x(3)/2)).^2./(1+((E-x(2))/(x(3)/2)).^2)-x(7)./(1+((E-x(5))/(0.5*x(6))).^2)); % Residual function
            x=lsqnonlin(F_fano,double(data_FL(:,2)'),[],[],options_l); % unrestricted fit with initial guess of starting values
            
            [x_bkg,~] = fmincon(@Fano,x,[],[],[],[],double(data_FL(:,1)'),double(data_FL(:,3)'),@conf_fano,options); % constrained background subtraction with lower and upper boundary conditions
            
            BKG=Data(start:End,2)-F_fano(x_bkg);
            
            delete('test.txt');
            
            figure('units','normalized','outerposition',[0 0 1 1]), hold on
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
                close all
                A=[xx(start:End,1),smooth(F_fano(x_bkg))]; % MATLAB function "smooth" is used to avoid the incorrect assignment of peaks to nanotubes
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
                    q=input('Enter the starting and ending wavelengths for the region of interest, e.g. 400 1400. \nBe aware that the S11 and S22 transitions to be fitted should be in the region defined: ','s');
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
                    if q(1)>q(2)
                        warning('Please enter starting wavelength that is smaller than the ending wavelength.')
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
            
            if parallel_fit=='y'
                if rel>=2014
                    options = optimoptions('fmincon','Algorithm','interior-point','Display','off','UseParallel',true);
                else
                    options = optimoptions('fmincon','Algorithm','interior-point','Display','off','UseParallel','always');
                end
            else
                options = optimoptions('fmincon','Algorithm','interior-point','Display','off');
            end
            
            xn = fmincon(@Nair,[.2 0.02],[],[],[],[],[],[],@conf_nair,options); % Approach of Nair et al. requires two starting values for k (0.2) and b (0.002).
            
            F_nair=@(x) double(yy(start:End,1)-x(1)./xx(start:End,1).^x(2)); % Background subtracted data!
            
            BKG=yy(start:End,1)-F_nair(xn);
            
            delete('test.txt');
            
            figure('units','normalized','outerposition',[0 0 1 1]), hold on
            
            plot(xx(start:End,1),yy(start:End,1),'k','linewidth',2)
            plot(xx(start:End,1),BKG,'r','linewidth',2)
            xlabel('Wavelength (nm)','FontSize',33)
            ylabel('Absorption (a.u.)','FontSize',33)
            legend('Measured Data','Background')
            set(gca,'LineWidth',2)
            set(gca,'FontSize',30)
            xlim(q)
            ylim([min(BKG)*0.9 max(BKG)*1.1])
            
            opt=input('Are you satisfied with the background subtraction (y/n)? ','s');
            if opt=='y'
                close all
                A=[xx(start:End,1),smooth(F_nair(xn))]; % MATLAB function "smooth" is used to avoid the incorrect assignment of peaks to nanotubes
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
            end
            
        case 3
            
            if rr==2
                while 1>0
                    q=input('Enter the starting and ending wavelengths for the region of interest, e.g. 400 1400. \nBe aware that the S11 and S22 transitions to be fitted should be in the region defined: ','s');
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
            
            if parallel_fit=='y'
                if rel>=2014
                    options = optimoptions('fmincon','Algorithm','interior-point','Display','off','UseParallel',true);
                else
                    options = optimoptions('fmincon','Algorithm','interior-point','Display','off','UseParallel','always');
                end
            else
                options = optimoptions('fmincon','Algorithm','interior-point','Display','off');
            end
            
            xn = fmincon(@Naumov,[.2 0.00002],[],[],[],[],[],[],@conf_naumov,options); % Approach of Naumov et al. requires two starting values for A (0.2) and b (0.002).
            
            F_naumov=@(x) double(yy(start:End,1)-x(1)*exp(-x(2).*xx(start:End,1))); % Background subtracted data!
            
            BKG=yy(start:End,1)-F_naumov(xn);
            
            delete('test.txt');
            
            figure('units','normalized','outerposition',[0 0 1 1]), hold on
            
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
                close all
                A=[xx(start:End,1), smooth(F_naumov(xn))]; % MATLAB function "smooth" is used to avoid the incorrect assignment of peaks to nanotubes
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
            end
            
        case 4
            
            xx=round(Data(1,1)):1:round(Data(end,1));
            yy=interp1(Data(:,1),Data(:,2),xx);
            xx=xx'; yy=yy';
            
            start=find(xx==q(1));
            End=find(xx==q(2));
            
            A=[xx(start:End,1),smooth(yy(start:End,1))]; % MATLAB function "smooth" is used to avoid the incorrect assignment of peaks to nanotubes
            
            close all
            
            break
            
        case 5
            
            disp('Select the file containing the background data');
            filename = uigetfile('*.txt','Select the file containing the background data');
            
            string='%f32 %f32';
            
            fid=fopen(filename);
            backG = textscan(fid,string,'Delimiter',',','HeaderLines',2);
            fclose(fid);
            
            BG(:,1)=backG{1};
            BG(:,2)=backG{2};
            
            if BG(1,1)>BG(2,1)
                BG=flipud(BG);
            end
            
            xx=round(Data(1,1)):1:round(Data(end,1));
            yy=interp1(Data(:,1),Data(:,2),xx);
            xx=xx'; 
            
            if size(yy,1)<size(yy,2)
                yy=yy';
            end
            
            start=find(xx==q(1));
            End=find(xx==q(2));
            
            xx_bg=round(BG(1,1)):1:round(BG(end,1));
            yy_bg=interp1(BG(:,1),BG(:,2),xx);
            xx_bg=xx_bg'; 
            
            if size(yy_bg,1)<size(yy_bg,2)
                yy_bg=yy_bg';
            end
            
            start_bg=find(xx_bg==q(1));
            End_bg=find(xx_bg==q(2));
            
            BKG=yy(start:End,1)-yy_bg(start_bg:End_bg);
            
            figure('units','normalized','outerposition',[0 0 1 1]), hold on
            
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
                close all
                A=[xx(start:End,1), smooth(yy(start:End,1)-yy_bg(start_bg:End_bg))]; % MATLAB function "smooth" is used to avoid the incorrect assignment of peaks to nanotubes
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
            end            
    end
    
end

% Prompt user to load the reference data table. This data table must
% consist of 3 columns with one header line. The column order is: 1 -
% string that contains the (n,m) species, 2 - S11 position in nm and 3 -
% S22 position in nm. (n,m) species must be sorted in ascending S11
% position. An example data set is provided in the file
% "Reference_Data.txt" for the HiPco semiconducting material based on the
% work of Bachilo et al. (Science, 2002, 298, 2361).

disp('Select the file containing the reference data for SWCNTs [(n,m), S11 and S22 position]');

ref = uigetfile('*.txt','Select the file containing the reference data for SWCNTs [(n,m), S11 and S22 position]');

fid=fopen(ref);
B = textscan(fid,'%s %f32 %f32 %f32 %f32','HeaderLines',1);
fclose(fid);

Names=B{1};

% Based on the data set provided, the "n" and "m" indicies of each nanotube
% are extracted and the diameter is calculated according to Pipes et al.
% (Compos. Sci. Technol., 2003, 63, 1349).

n=zeros(length(Names),1);
m=zeros(length(Names),1);
diam=zeros(length(Names),1);

for i=1:length(Names)
    a=Names{i,1};
    b=char(a);
    if strcmp(',',b(3))==1
        n(i)=str2double(b(2));
        m(i)=str2double(b(4));
    else
        n(i)=str2double(b(2:3));
        if strcmp(')',b(6))==1
            m(i)=str2double(b(5));
        else
            m(i)=str2double(b(5:6));
        end
    end
    clear a
    clear b
    diam(i)=sqrt(n(i)^2+m(i)^2+n(i)*m(i))*c_len/(pi)*sqrt(3);
end

% Prompt the user to decide whether they want to fit S11 or S22. In the
% event that S11 is chosen, the user will later on be given the option to
% fit the entire region of interest defined during background subtraction.
% Where the fit in the S22 region is constrained by the (n,m) species
% defined for the S11 region.

fit_region=input('Do you want to fit S11 (1) or S22 (2) ?: ','s');
fit_region=sscanf(fit_region,'%i');

switch fit_region
    case 1
        range=input('Please enter the wavelength range (in nm) in which the SWCNTs S11 shall be fitted, e.g. 820 1400. \nEnsure that this range falls within the region of interest defined during background subtraction: ','s');
        range=sscanf(range,'%i');
        [S11_sort,so]=sort(B{2});
        S22_sort=B{3};
        diam=diam(so);
        n=n(so);
        m=m(so);
        Names=Names(so);
        i_s=find(S11_sort>range(1));
        i_start=i_s(1);
    case 2
        range=input('Please enter the wavelength range (in nm) in which the SWCNTs S22 shall be fitted, e.g. 450 800. \nEnsure that this range falls within the region of interest defined during background subtraction: ','s');
        range=sscanf(range,'%i');
        S11_sort=B{2};
        [S22_sort,so]=sort(B{3});
        diam=diam(so);
        n=n(so);
        m=m(so);
        Names=Names(so);
        i_s=find(S22_sort>range(1));
        i_start=i_s(1);
end

lambda=A(:,1);

start=find(lambda==range(1));
End=find(lambda==range(2));

% A color scheme is defined for different (n,m) species; from blue (small
% diameter) to red (large diameter).

colors=zeros(1,3);

switch fit_region
    case 1
        for i=i_start:length(S11_sort)
            if S11_sort(i)<range(2)
                pos(i-i_start+1)=find(A(:,1)==S11_sort(i));
                colors(i-i_start+1,:)=[0+(S11_sort(i)-range(1))/(range(2)-range(1)),0,1-(S11_sort(i)-range(1))/(range(2)-range(1))];
            end
        end
    case 2
        for i=i_start:length(S22_sort)
            if S22_sort(i)<range(2)
                pos(i-i_start+1)=find(A(:,1)==S22_sort(i));
                colors(i-i_start+1,:)=[0+(S22_sort(i)-range(1))/(range(2)-range(1)),0,1-(S22_sort(i)-range(1))/(range(2)-range(1))];
            end
        end
end

Names_n=Names(i_start:length(pos)+i_start-1);
S11_sort=S11_sort(i_start:length(pos)+i_start-1);
S22_sort=S22_sort(i_start:length(pos)+i_start-1);
n=n(i_start:length(pos)+i_start-1);
m=m(i_start:length(pos)+i_start-1);
diam=diam(i_start:length(pos)+i_start-1);

% Calculation of full width at half maximum (FWHM) based on the work of
% Tune et al. (Energy Environ. Sci., 2013, 6, 2572-2577) with a correction
% factor of 0.5. This calculation takes the values of S11 and S22 from the
% reference data table and converts them into eV and applies the formula of
% Tune et al. The FWHM in eV is then reverted to nm by calculating the
% following difference:
% lambda1 = h*c / (E_S11 + fwhm_ev/2) * 1e9
% lambda2 = h*c / (E_S11 - fwhm_ev/2) * 1e9
% fwhm_nm = lambda2 - lambda1

fwhm_inp=input('Choice of FWHM:\n (1) - Default (29.86 meV for S11 and 57.96 meV for S22)\n (2) - Constant in eV\n (3) - Constant in nm\n (4) - User defined function\n  ','s');
fwhm_inp=sscanf(fwhm_inp,'%i');

switch fwhm_inp
    case 1
        E_S11=h*c./(1e-9*S11_sort);
        E_S22=h*c./(1e-9*S22_sort);
        
        fwhm_ev=29.86e-3;
        fwhm_ev_22=57.96e-3;
        
        fwhm_nm=fwhm_ev*h*c./(E_S11.^2-(fwhm_ev/2).^2)*1e9;
        fwhm_nm_22=fwhm_ev_22*h*c./(E_S22.^2-(fwhm_ev_22/2).^2)*1e9;
    case 2
        fwhm_user=input('Please enter two FWHM values in meV for S11 and S22, e.g. 20 30: ','s');
        fwhm_user=sscanf(fwhm_user,'%i');
        
        E_S11=h*c./(1e-9*S11_sort);
        E_S22=h*c./(1e-9*S22_sort);
        
        fwhm_ev=fwhm_user(1);
        fwhm_ev_22=fwhm_user(2);
        
        fwhm_nm=fwhm_ev*h*c./(E_S11.^2-(fwhm_ev/2).^2)*1e9;
        fwhm_nm_22=fwhm_ev_22*h*c./(E_S22.^2-(fwhm_ev_22/2).^2)*1e9;
    case 3
        fwhm_user=input('Please enter two FWHM values in nm for S11 and S22, e.g. 20 30: ','s');
        fwhm_user=sscanf(fwhm_user,'%i');
        
        for i=1:length(S11_sort)
            fwhm_nm(i)=fwhm_user(1);
            fwhm_nm_22(i)=fwhm_user(2);
        end
    case 4
        fct=input('Please specify whether you have a diameter dependent function that calculates the FWHM in\n (1) - nm\n (2) - eV\n (3) - cm^-1\n Or an absorption position dependent function that calculates the FWHM in\n (4) - nm\n (5) - eV\n (6) - cm^-1\n ','s');
        fct=sscanf(fct,'%i');
        
        switch fct
            case 1
                str = input('Give an equation in d (in nm) for the S11-region, e.g. 2*d^2-3*d+4\n If you want a constant function, enter d*0 + constant: ','s');
                f = inline(str,'d');
                fwhm_nm = feval(f,diam);
                
                str2 = input('Give an equation in d (in nm) for the S22-region, e.g. 2*d^2-3*d+4\n If you want a constant function, enter d*0 + constant: ','s');
                f2 = inline(str2,'d') ;
                fwhm_nm_22 = feval(f2,diam) ;
            case 2
                str = input('Give an equation in d (in nm) for the E11-region, e.g. 2*d^2-3*d+4\n If you want a constant function, enter d*0 + constant: ','s');
                f = inline(str,'d');
                fwhm_ev = feval(f,diam);
                
                str2 = input('Give an equation in d (in nm) for the E22-region, e.g. 2*d^2-3*d+4\n If you want a constant function, enter d*0 + constant: ','s');
                f2 = inline(str2,'d') ;
                fwhm_ev_22 = feval(f2,diam) ;
                
                E_S11=h*c./(1e-9*S11_sort);
                E_S22=h*c./(1e-9*S22_sort);
                
                fwhm_nm=fwhm_ev*h*c./(E_S11.^2-(fwhm_ev/2).^2)*1e9;
                fwhm_nm_22=fwhm_ev_22*h*c./(E_S22.^2-(fwhm_ev_22/2).^2)*1e9;
            case 3
                str = input('Give an equation in d (in nm) for the S11-region, e.g. 2*d^2-3*d+4\n If you want a constant function, enter d*0 + constant: ','s');
                f = inline(str,'d');
                fwhm_cm = feval(f,diam);
                
                str2 = input('Give an equation in d (in nm) for the S22-region, e.g. 2*d^2-3*d+4\n If you want a constant function, enter d*0 + constant: ','s');
                f2 = inline(str2,'d') ;
                fwhm_cm_22 = feval(f2,diam) ;
                
                fwhm_nm=1e7*(-(1e7./S11_sort+fwhm_cm./2).^(-1)+(1e7./S11_sort-fwhm_cm./2).^(-1));
                fwhm_nm_22=1e7*(-(1e7./S22_sort+fwhm_cm_22./2).^(-1)+(1e7./S22_sort-fwhm_cm_22./2).^(-1));
            case 4
                str = input('Give an equation in s1 for the S11-region, e.g. 2*s1^2-3*s1+4\n If you want a constant function, enter s1*0 + constant: ','s');
                f = inline(str,'s1');
                fwhm_nm = feval(f,S11_sort);
                
                str2 = input('Give an equation in s2 for the S22-region, e.g. 2*s2^2-3*s2+4\n If you want a constant function, enter s2*0 + constant: ','s');
                f2 = inline(str2,'s2') ;
                fwhm_nm_22 = feval(f2,S22_sort);
            case 5
                E_S11=h*c./(1e-9*S11_sort);
                E_S22=h*c./(1e-9*S22_sort);
                
                str = input('Give an equation in e1 for the E11-region, e.g. 2*e1^2-3*e1+4\n If you want a constant function, enter e1*0 + constant: ','s');
                f = inline(str,'e1');
                fwhm_ev = feval(f,E_S11);
                
                str2 = input('Give an equation in e2 for the E22-region, e.g. 2*e2^2-3*e2+4\n If you want a constant function, enter e2*0 + constant: ','s');
                f2 = inline(str2,'e2') ;
                fwhm_ev_22 = feval(f2,E_S22);
                
                fwhm_nm=fwhm_ev*h*c./(E_S11.^2-(fwhm_ev/2).^2)*1e9;
                fwhm_nm_22=fwhm_ev_22*h*c./(E_S22.^2-(fwhm_ev_22/2).^2)*1e9;
            case 6
                Cm_S11=1e7./S11_sort;
                Cm_S22=1e7./S22_sort;
                
                str = input('Give an equation in s1 for the S11-region, e.g. 2*s1^2-3*s1+4\n If you want a constant function, enter s1*0 + constant: ','s');
                f = inline(str,'s1');
                fwhm_cm = feval(f,Cm_S11);
                
                str2 = input('Give an equation in s2 for the S22-region, e.g. 2*s2^2-3*s2+4\n If you want a constant function, enter s2*0 + constant: ','s');
                f2 = inline(str2,'s2') ;
                fwhm_cm_22 = feval(f2,Cm_S22);
                
                fwhm_nm=1e7*(-(Cm_S11+fwhm_cm./2).^(-1)+(Cm_S11-fwhm_cm./2).^(-1));
                fwhm_nm_22=1e7*(-(Cm_S22+fwhm_cm_22./2).^(-1)+(Cm_S22-fwhm_cm_22./2).^(-1));
        end
end


% Initial peak intensity of the individual SWCNTs is determined based on
% the solution absorption measurement and the upper limit of the selected
% wavelength range.

height_init=A(pos(1:end),2);

Absorption=double(A(:,2));

% Save user input in a protocol file.

if protocol=='y'
    bla=now;
    time_stamp=datestr(bla);
    time_stamp=strrep(time_stamp,' ','_');
    time_stamp=strrep(time_stamp,':','-');
    
    fid=fopen(['Protocol_',time_stamp,'.txt'],'a');
    array{1}=['Filename: ',filename];
    array{2}=['Wavelength regime for background subtraction (nm): ',num2str(q')];
    array{3}=['Background subtraction method: ',num2str(w)];
    array{4}=['Reference data: ',ref];
    
    if fit_region==1
        array{5}=['Wavelength regime for S11 region (in nm): ',num2str(range')];
    else
        array{5}=['Wavelength regime for S22 region (in nm): ',num2str(range')];
    end
    
    fprintf(fid,'%s \n',array{:});
else
    time_stamp=[];
    fid=[];
end


%--------------------------------------------------------------------------
% Section 1.1 - Background Subtraction
%--------------------------------------------------------------------------


% Method (1) based on Tian et al.

function err = Fano(x)

A=dlmread('test.txt');
c = A(:,2)-x(4)*(x(1)+(A(:,1)-x(2))/(x(3)/2)).^2./(1+((A(:,1)-x(2))/(x(3)/2)).^2)-x(7)./(1+((A(:,1)-x(5))/(0.5*x(6))).^2);
err = double(c'*c);

function [c,ceq] = conf_fano(x) % Constraint function, that forces the background to be smaller than the absorption data for every single wavelength

A=dlmread('test.txt');
% Nonlinear inequality constraints
c = double(x(4)*(x(1)+(A(:,1)-x(2))/(x(3)/2)).^2./(1+((A(:,1)-x(2))/(x(3)/2)).^2)+x(7)./(1+((A(:,1)-x(5))/(0.5*x(6))).^2) - A(:,2));
% Nonlinear equality constraints
ceq = [];

% Method (2) based on Nair et al.

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

% Method (3) based on Naumov et al.

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


%..........................................................................


%--------------------------------------------------------------------------
% Section 2 - "select_nanotubes"
%--------------------------------------------------------------------------


function [vec,diam_n,colors_n,S11n,S22n,CNT_n,fwhm_nm_n,fwhm_nm_n22]=select_nanotubes(lambda,pos,height_init,S11_sort,S22_sort,colors,Names,diam,area_part,fwhm_gauss,Absorption,h,c,fwhm_nm,fwhm_nm_22,n,m,fit_region,protocol,fid)

% Using the center position, initial heights and FWHM (obtained in
% "get_data") of the (n,m) species within the user defined wavelength
% regime, Lorentzian profiles and Gaussian exciton phonon sidebands (EPS)
% are displayed to assist the user in choosing the appropriate (n,m)
% distribution.
%
% Lorentzian line shapes are defined as:
% y_L = height / (1 + ((lambda-S11)/(0.5*FWHM_L))^2)
%
% Gaussian line shapes are defined as:
% y_G = area_L * spectral_weight * exp(-2 * (lambda - EPS position)^2 / (FWHM_G / sqrt(ln(4)) )^2) / (FWHM_G * sqrt(pi / (2 * ln(4))))
% with area_L = height * pi * FWHM_L / 2
% and spectral weight transfer = 0.017 + 0.1 / diam_S11 + correction_factor
%
% After displaying all detected (n,m) species within the wavelength range,
% the user is given the option to remove certain nanotubes. In this regard,
% a photoluminescence contour map can be helpful.

switch fit_region
    case 1 % S11
        S_sort=S11_sort;
        FWHM=fwhm_nm;
    case 2 % S22
        S_sort=S22_sort;
        FWHM=fwhm_nm_22;
end

lamb_phonon = h*c./(h*c./S_sort*1e9+0.2)*1e9; % Peak position of EPS

figure('units','normalized','outerposition',[0 0 1 1]),hold on

y=zeros(length(lambda),length(pos)); % Lorentzian peaks
y_g=zeros(length(lambda),length(pos)); % Gaussian EPS peaks
legstr=cell(1,length(pos));

for j=1:length(pos)
    y(:,j)=height_init(j) ./ (1 + ((lambda-S_sort(j))/(0.5*FWHM(j))).^2);
    plot(lambda,y(:,j),'color',colors(j,:),'linewidth',1);
    text(double(S_sort(j)),double(height_init(j)),num2str(j),'color',[0 .4 0],'fontsize',20)
    legstr{j}=[num2str(j),' = ',Names{j,1}]; % save SWCNT chirality in a string and display it in the legend of the plot
end
for j=1:length(pos)
    y_g(:,j)=height_init(j)*pi*FWHM(j)/2*(0.017+0.1/diam(j)+area_part)*exp(-2*(lambda-lamb_phonon(j)).^2/(fwhm_gauss/sqrt(log(4)))^2)/(fwhm_gauss*sqrt(pi/(2*log(4))));
end

plot(lambda,Absorption,'k','linewidth',2)

for j=1:length(pos)
    plot(lambda,y_g(:,j));
end

title('(n,m)-Species Available in Selected Region','FontSize',28)
legend([legstr,{'Measured Spectrum'}])

xlabel('Wavelength (nm)','FontSize',24)
ylabel('Absorption (a.u.)','FontSize',24)
set(gca,'LineWidth',2)
set(gca,'FontSize',22)

clear legstr

w=input('Based on the data table provided the following SWCNTs could be found within this wavelength range. \nPlease enter number of peak(s) corresponding to SWCNTs you want to remove, e.g. 1 4 2. \nPress "Enter" if you want to keep them all: ','s');
w=sscanf(w,'%i');

close all

% After selecting the (n,m) species to be removed from the spectrum, the
% vector "vec" contains all SWCNTs to be considered. This is then used to
% adjust the vectors originally defined in "get_data".

vec=1:length(pos);
vec(w)=[];

posn=pos(vec);
diam_n=diam(vec);
colors_n=colors(vec,:);
S11n=S11_sort(vec);
S22n=S22_sort(vec);
fwhm_nm_n=fwhm_nm(vec);
fwhm_nm_n22=fwhm_nm_22(vec);
CNT_n(:,1)=n(vec);
CNT_n(:,2)=m(vec);

if protocol=='y'
    array{1}=['Number of SWCNTs that were not considered during fit: ',num2str(w')];
    fprintf(fid,'%s \n',array{:});
end


%..........................................................................


%--------------------------------------------------------------------------
% Section 3 - "check_S22"
%--------------------------------------------------------------------------


function [height_22,center_22,FWHM_22,legstr_22,FWHM_g_22,f1_22,colors,w,y,diam_22,phonon_pos2,fwhm_11,height_11,area_11,folder_name,legstr_help]=check_S22(lambda,h,c)

% After prompting the user to select the folder containing the results of
% a previous fit of the entire region, all variables required for the fit
% of the S22 region, including metallic SWCNTs, are loaded.

disp('Please select the folder containing the variables for the film fit, that was automatically created during the solution fit.')
folder_name=uigetdir(pwd,'Please select the folder containing the variables for the film fit, that was automatically created during the solution fit.');

oldD=cd(folder_name);

D=dir;

for i=1:size(D,1)
    a=strfind(D(i,1).name,'Results');
    if isempty(a)==0
        file=D(i,1).name;
        r=load(file);
        Results=r.Results; % All relevant results from the solution fit
        break
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

cd(oldD);

% Save the results of thesolution absorption fit.

switch w
    case {1,2} % Lorentzian or Gaussian Fit
        
        numb_CNT=size(Results,1)-length(phonon_pos)-length(phonon_pos2); % number of CNTs that were used to fit the solution data
        
        fwhm_11=[];
        height_11=[];
        area_11=[];
        
        center_22=[];
        FWHM_22=[];
        area_22=[];
        legstr_22=[];
        height_22=[];
        diam_22=[];
        legstr_help=[];
        y=zeros(length(lambda),numb_CNT-1);
        
        if size(Results,2)==8
            
            for i=2:numb_CNT
                center_22(i-1)=Results{i,3};
                FWHM_22(i-1)=Results{i,4};
                area_22(i-1)=Results{i,7};
                height_22(i-1)=Results{i,6};
                diam_22(i-1)=Results{i,2};
                legstr_22{i-1}=legstr_sol{i-1};
                if w==1
                    y(:,i-1)=height_22(i-1) ./ (1 + ((lambda-center_22(i-1))/(0.5*FWHM_22(i-1))).^2);
                else
                    y(:,i-1)=height_22(i-1)*exp(-log(2)*((lambda-center_22(i-1))/(FWHM_22(i-1)/2)).^2);
                end
            end
            
            FWHM_g_22=[]; f1_22=[];
            
            phonon_pos2=phonon_pos;
            
            for j=1:length(phonon_pos2)
                FWHM_g_22(j)=Results{numb_CNT+j,4};
                f1_22=Results{numb_CNT+j,6};
                
                y(:,j+length(center_22)) = area_22(phonon_pos2(j))*(0.017+0.1/diam_22(phonon_pos2(j)) + f1_22)*exp(-2*(lambda(1:end)-h*c/(1e-9*(h*c/((center_22(phonon_pos2(j)))*1e-9)+0.2))).^2/(FWHM_g_22(j)/sqrt(log(4)))^2)/((FWHM_g_22(j))*sqrt(pi/(2*log(4))));
            end
            
        else
            
            n=0;
            for i=2:numb_CNT
                if isempty(Results{i,8})==0
                    center_22(i-1)=Results{i,8};
                    FWHM_22(i-1)=Results{i,9};
                    area_22(i-1)=Results{i,12};
                    height_22(i-1)=Results{i,11};
                    diam_22(i-1)=Results{i,2};
                    legstr_22{i-1-n}=legstr_sol{numb_CNT+i-2-n};
                    if w==1
                        y(:,i-1)=height_22(i-1) ./ (1 + ((lambda-center_22(i-1))/(0.5*FWHM_22(i-1))).^2);
                    else
                        y(:,i-1)=height_22(i-1)*exp(-log(2)*((lambda-center_22(i-1))/(FWHM_22(i-1)/2)).^2);
                    end
                else
                    n=n+1;
                    legstr_help{n}=strrep(legstr_sol{i-1},'S_1_1-','');
                end
                fwhm_11(i-1)=Results{i,4};
                height_11(i-1)=Results{i,6};
                area_11(i-1)=Results{i,7};
            end
            
            FWHM_g_22=[]; f1_22=[];
            
            for j=1:length(phonon_pos2)
                FWHM_g_22(j)=Results{numb_CNT+length(phonon_pos)+j,9};
                f1_22=Results{numb_CNT+length(phonon_pos)+j,11};
                
                y(:,j+length(center_22)) = area_22(phonon_pos2(j))*(0.017+0.1/diam_22(phonon_pos2(j)) + f1_22)*exp(-2*(lambda(1:end)-h*c/(1e-9*(h*c/((center_22(phonon_pos2(j)))*1e-9)+0.2))).^2/(FWHM_g_22(j)/sqrt(log(4)))^2)/((FWHM_g_22(j))*sqrt(pi/(2*log(4))));
            end
        end
        
    case 3 % Voigtian
        
        numb_CNT=size(Results,1)-length(phonon_pos)-length(phonon_pos2); % number of CNTs that were used to fit the solution data
        
        fwhm_11=[];
        height_11=[];
        area_11=[];
        
        center_22=[];
        FWHM_22=[];
        area_22=[];
        legstr_22=[];
        height_22=[];
        diam_22=[];
        legstr_help=[];
        y=zeros(length(lambda),numb_CNT-1);
        
        if size(Results,2)==11 % S22 was fitted
            
            for i=2:numb_CNT
                center_22(i-1)=Results{i,3};
                FWHM_22(i-1,1)=Results{i,4};
                FWHM_22(i-1,2)=Results{i,5};
                FWHM_22(i-1,3)=Results{i,6};
                area_22(i-1)=Results{i,10};
                height_22(i-1)=Results{i,9};
                diam_22(i-1)=Results{i,2};
                legstr_22{i-1}=legstr_sol{i-1};
                wg=FWHM_22(i-1,2);
                XX=sqrt(log(2))*(lambda-center_22(i-1))/(wg/2); % wg/2 is used because the FWHM is saved and not the HWHM
                YY=sqrt(log(2))*FWHM_22(i-1,1)/wg;
                y(:,i-1)=height_22(i-1).*real(complexErrorFunction(XX,YY))./real(complexErrorFunction(0,YY));
            end
            
            FWHM_g_22=[]; f1_22=[];
            
            for i=1:length(phonon_pos)
                FWHM_g_22(i)=Results{numb_CNT+i,5};
                f1_22=Results{numb_CNT+i,9};
                
                y(:,i+length(center_22)) = area_22(phonon_pos(i))*(0.017+0.1/diam_22(phonon_pos(i)) + f1_22)*exp(-2*(lambda(1:end)-h*c/(1e-9*(h*c/((center_22(phonon_pos(i)))*1e-9)+0.2))).^2/(FWHM_g_22(i)/sqrt(log(4)))^2)/((FWHM_g_22(i))*sqrt(pi/(2*log(4))));
            end
            
            phonon_pos2=phonon_pos;
            
        else % Entire region was fitted
            
            n=0;
            for i=2:numb_CNT
                if isempty(Results{i,11})==0
                    center_22(i-1)=Results{i,11};
                    FWHM_22(i-1,1)=Results{i,12};
                    FWHM_22(i-1,2)=Results{i,13};
                    FWHM_22(i-1,3)=Results{i,14};
                    area_22(i-1)=Results{i,18};
                    height_22(i-1)=Results{i,17};
                    diam_22(i-1)=Results{i,2};
                    legstr_22{i-1-n}=legstr_sol{numb_CNT+i-2-n};
                    wg=FWHM_22(i-1,2);
                    XX=sqrt(log(2))*(lambda-center_22(i-1))/(wg/2); % wg/2 is used because the FWHM is saved and not the HWHM
                    YY=sqrt(log(2))*FWHM_22(i-1,1)/wg;
                    y(:,i-1)=height_22(i-1).*real(complexErrorFunction(XX,YY))./real(complexErrorFunction(0,YY));
                else
                    n=n+1;
                    legstr_help{n}=strrep(legstr_sol{i-1},'S_1_1-','');
                end
                fwhm_11(i-1,1)=Results{i,4};
                fwhm_11(i-1,2)=Results{i,5};
                fwhm_11(i-1,3)=Results{i,6};
                height_11(i-1)=Results{i,8};
                area_11(i-1)=Results{i,10};
            end
            
            FWHM_g_22=[]; f1_22=[];
            
            for i=1:length(phonon_pos2)
                FWHM_g_22(i)=Results{numb_CNT+length(phonon_pos)+i,13};
                f1_22=Results{numb_CNT+length(phonon_pos)+i,17};
                
                y(:,i+length(center_22)) = area_22(phonon_pos2(i))*(0.017+0.1/diam_22(phonon_pos2(i)) + f1_22)*exp(-2*(lambda(1:end)-h*c/(1e-9*(h*c/((center_22(phonon_pos2(i)))*1e-9)+0.2))).^2/(FWHM_g_22(i)/sqrt(log(4)))^2)/((FWHM_g_22(i))*sqrt(pi/(2*log(4))));
            end
        end
end

cd(oldD);


%..........................................................................


%--------------------------------------------------------------------------
% Section 4 - "select_nanotubes_metallic"
%--------------------------------------------------------------------------


function [diam_m,colors_mn,M11n,fwhm_nm_m11_n,CNT_n,height_m]=select_nanotubes_metallic(c_len,Absorp,lambda,h,c,y2,fid_p,protocol)

% This section repeats the procedure outlined in Section 2 to assign
% metallic carbon nanotubes.

disp('Select the file containing the reference data for metallic SWCNTs [(n,m), M11 position]');

ref = uigetfile('*.txt','Select the file containing the reference data for metallic SWCNTs [(n,m), M11 position]');

fid=fopen(ref);
B = textscan(fid,'%s %f32 %f32','HeaderLines',1);
fclose(fid);

Names=B{1};

n=zeros(length(Names),1);
m=zeros(length(Names),1);
diam=zeros(length(Names),1);

for i=1:length(Names)
    a=Names{i,1};
    b=char(a);
    if strcmp(',',b(3))==1
        n(i)=str2double(b(2));
        m(i)=str2double(b(4));
    else
        n(i)=str2double(b(2:3));
        if strcmp(')',b(6))==1
            m(i)=str2double(b(5));
        else
            m(i)=str2double(b(5:6));
        end
    end
    clear a
    clear b
    diam(i)=sqrt(n(i)^2+m(i)^2+n(i)*m(i))*c_len/(pi)*sqrt(3);
end

M11_sort_nm=B{2};
M11_sort_nm=round(M11_sort_nm);
M11_sort_ev=B{3};

fwhm_inp=input('Choice of FWHM:\n (1) - Default (93.42 meV for M11)\n (2) - Constant in eV\n (3) - Constant in nm\n (4) - User defined function\n  ','s');
fwhm_inp=sscanf(fwhm_inp,'%i');

switch fwhm_inp
    case 1
        fwhm_ev_m11=93.42e-3;
        fwhm_nm_m11=fwhm_ev_m11*h*c./(M11_sort_ev.^2-(fwhm_ev_m11/2).^2)*1e9;        
    case 2
        fwhm_user=input('Please enter the FWHM in meV for M11, e.g. 80: ','s');
        fwhm_user=sscanf(fwhm_user,'%i');
                
        fwhm_ev=fwhm_user(1);
        
        fwhm_nm_m11=fwhm_ev*h*c./(M11_sort_ev.^2-(fwhm_ev/2).^2)*1e9;
    case 3
        fwhm_user=input('Please enter the FWHM in nm for M11, e.g. 20: ','s');
        fwhm_user=sscanf(fwhm_user,'%i');
        
        for i=1:length(M11_sort_nm)
            fwhm_nm_m11(i)=fwhm_user(1);
        end
    case 4
        fct=input('Please specify whether you have diameter (1) or absorption position (in eV) (2) dependent function: ','s');
        fct=sscanf(fct,'%i');
        
        switch fct
            case 1
                str = input('Give an equation in d for the M11-region, e.g. 2*d^2-3*d+4\n If you want a constant function, enter e1*0 + constant: ','s');
                f = inline(str,'d');
                fwhm_nm_m11 = feval(f,diam);
            case 2                
                str = input('Give an equation in m1 for the M11-region, e.g. 2*m1^2-3*m1+4\n If you want a constant function, enter e1*0 + constant: ','s');
                f = inline(str,'e1');
                fwhm_ev = feval(f,M11_sort_ev);
                                
                fwhm_nm_m11=fwhm_ev*h*c./(M11_sort_ev.^2-(fwhm_ev/2).^2)*1e9;
        end
end

i_s=find(M11_sort_nm>lambda(1));
i_start=i_s(1);

colors_m=zeros(1,3);

for i=1:length(M11_sort_nm)
    if M11_sort_nm(i)<lambda(end)
        pos(i)=find(lambda==M11_sort_nm(i));
        colors_m(i,:)=[1,(M11_sort_nm(i)-lambda(1))/(lambda(end)-lambda(1)),1-(M11_sort_nm(i)-lambda(1))/(lambda(end)-lambda(1))];
    end
end

figure('units','normalized','outerposition',[0 0 1 1]),hold on

for j=i_start:length(pos)
    y(:,j-i_start+1)=.1*max(Absorp) ./ (1 + ((lambda-M11_sort_nm(j))/(0.5*fwhm_nm_m11(j))).^2);
    plot(lambda,y(:,j-i_start+1),'color',colors_m(j,:),'linewidth',1);
    text(double(M11_sort_nm(j)),double(.1*max(Absorp)),num2str(j-i_start+1),'color',[0 .4 0],'fontsize',20)
    legstr{j-i_start+1}=[num2str(j-i_start+1),' = ',Names{j,1}]; % save SWCNT chirality in a string and display it in the legend of the plot
end

plot(lambda,Absorp,'k','linewidth',2)
plot(lambda,sum(y2,2),'color',[0 0.4 0],'linewidth',2)

for j=1:size(y2,2)
    plot(lambda,y2(:,j),'color',[0.6 0.6 0.6],'linewidth',1);
end

residual=abs(Absorp-sum(y2,2));

title('(n,m)-Metallic Species Available in Selected Region','FontSize',26)
legend([legstr,{'Measured Spectrum'},{'Previously Calculated Spectrum'},{'Previously Fitted (n,m) Species'}])

xlabel('Wavelength (nm)','FontSize',24)
ylabel('Absorption (a.u.)','FontSize',24)
set(gca,'LineWidth',2)
set(gca,'FontSize',22)

ww=input('Based on the data table provided the following SWCNTs could be found within this wavelength range. \nPlease enter number of peak(s) corresponding to SWCNTs you want to remove, e.g. 1 4 2. \nPress "Enter" if you want to keep them all: ','s');
ww=sscanf(ww,'%i');

vec=i_start:length(pos);
vec(ww)=[];

pos_m=pos(vec);
height_m=residual(pos_m);
diam_m=diam(vec);
colors_mn=colors_m(vec,:);
M11n=M11_sort_nm(vec);
fwhm_nm_m11_n=fwhm_nm_m11(vec);
CNT_n(:,1)=n(vec);
CNT_n(:,2)=m(vec);

if protocol == 'y'
    array{1}=['Reference file for metallic SWCNTs: ',ref];
    array{2}=['Number of metallic SWCNTs that were not considered during the fit: ',num2str(ww')];
    
    fprintf(fid_p,'%s \n',array{:});
end


%..........................................................................


%--------------------------------------------------------------------------
% Section 5 - "check_peaks"
%--------------------------------------------------------------------------


function [Order,l_g,center,center_22,Height,up_bound,low_bound]=check_peaks(max_shift,Absorption,lambda,Names,vec,CNT_n,start_l,S11n,S22n,fit_region,upper_lim,upper_lim_peak,lower_lim,lower_lim_peak,fid,protocol)

% Check whether S11 or S22 region shall be fitted

switch fit_region
    case 1 % S11
        S_n=S11n;
    case 2 % S22
        S_n=S22n;
end

% Find all peaks in the user provided wavelength range and sort them based
% on their absorption intensity (high-low). Based on the number N of (n,m)
% species defined by the user, plot only the N (or less) most intense
% peaks. Each peak found is assigned a numerical value indicating its
% intensity (1 - most intense peak, 2 - second most intense peak, etc.).
% The user is given the option to remove peaks. This is for example
% important when a known exciton phonon sideband has been assigned to an
% S11 or S22 transition of a nanotube. The code checks whether one nanotube
% got assigned to multiple peaks. If so, the user is warned and asked to
% remove the corresponding peak. Position of peaks is stored in vector
% "locs_n" and the peak height in vector "pks_n".

[pks,locs]=findpeaks(Absorption);

[~,b]=sort(pks,'descend');

locs_n=locs(b);
pks_n=pks(b);

locs_n(length(vec)+1:end)=[];
pks_n(length(vec)+1:end)=[];
b(length(vec)+1:end)=[];

rid_cnts=[];

while 1>0
    
    tt=figure('units','normalized','outerposition',[0 0 1 1]); hold on
    plot(lambda,Absorption,'k','linewidth',2)
    plot(locs_n+start_l-1,pks_n,'o','linewidth',2,'MarkerSize',10)
    for i=1:length(locs_n)
        text(double(locs_n(i)+start_l+2),double(pks_n(i)),num2str(i),'color',[0 .4 0],'fontsize',20)
    end
    xlabel('Wavelength (nm)','FontSize',24)
    ylabel('Absorption (a.u.)','FontSize',24)
    legend('Measured Data')
    title('Automatic Assignment of Peaks to SWCNTs','FontSize',26)
    set(gca,'LineWidth',2)
    set(gca,'FontSize',22)
    
    opt=input('Do you want to remove one or more peaks that were assigned to a SWCNT transition (y/n)? ','s');
    
    if opt=='y'
        rid=input('Please enter number of peak(s) not corresponding to SWCNT transitions that you want to remove, e.g. 1 4 2: ','s');
        rid_vec=sscanf(rid,'%i');
        b(rid_vec)=[];
        locs_n=locs(b);
        pks_n=pks(b);
        rid_cnts=[rid_cnts,rid_vec'];
    else
        for j=1:length(locs_n)
            peak=repmat(locs_n(j)+start_l-1,length(vec),1);
            min_dist=peak-S_n;
            [~,bb]=min(abs(min_dist));
            if min_dist(bb)<-10 && bb>1 % Prevent blue-shift of more than 10 nm
                bb=bb-1;
            end
            dist(j)=min_dist(bb);
            lok(j)=bb;
        end
        if isempty(locs_n)==0
            [n, bin] = histc(lok, unique(lok));
            multiple = find(n > 1);
            if isempty(find(ismember(bin, multiple), 1))==0
                warning(['This assignment would result in one SWCNT transition being assigned to multiple peaks: ',num2str(find(ismember(bin, multiple))),'. Maybe an exciton phonon sideband was mistaken for a SWCNT. In order to avoid an error you have to remove one or more peaks!'])
            else
                pks_shift=[];
                for j=1:length(dist)
                    if dist(j)>max_shift
                        pks_shift=[pks_shift,j];
                    end
                end
                if isempty(pks_shift)==1
                    close(tt)
                    break
                else
                    warning(['This assignment would result in SWCNT transitions: ',num2str(pks_shift),' being shifted more than ',num2str(max_shift),' nm. This migth be unphysical!'])
                    w_q=input('Do you want to remove (1) these peaks or continue (2) with the manual addition of peaks?\n','s');
                    w_q=sscanf(w_q,'%i');
                    if w_q==2
                        close(tt)
                        break
                    end
                end
                
            end
        else
            close(tt)
            break
        end
        clear lok
        clear n
        clear bin
        clear multiple
    end
    
    close(tt)
    
end

% Following the automatic assignment and user defined removal of peaks, the
% user is now given the option to add additional peaks by specifying a
% wavelength value. Additional peaks and their associated height are saved
% in "locs_n" and "pks_n". If no additional peaks are assigned by the user
% the code checks for nanotubes being assigned to multiple peaks by
% calculating the distance of each (n,m) transition to the peaks assigned
% by the user and returns an error if true.

added_cnts=[];

while 1>0
    
    tt=figure('units','normalized','outerposition',[0 0 1 1]); hold on
    plot(lambda,Absorption,'k','linewidth',2)
    plot(locs_n+start_l-1,pks_n,'o','linewidth',2,'MarkerSize',10)
    for i=1:length(locs_n)
        text(double(locs_n(i)+start_l+2),double(pks_n(i)),num2str(i),'color',[0 .4 0],'fontsize',20)
    end
    
    xlabel('Wavelength (nm)','FontSize',24)
    ylabel('Absorption (a.u.)','FontSize',24)
    legend('Measured Data')
    title('Manual Assignment of Peaks to SWCNTs','FontSize',26)
    set(gca,'LineWidth',2)
    set(gca,'FontSize',22)
    
    opt=input('Do you want to add more peaks that correspond to SWCNT transitions (y/n)? ','s');
    
    if opt=='y'
        add=input('Please enter as many positions in "nm" that correspond to nanotube transitions that are missing, e.g. 450 1030: ','s');
        add_vec=sscanf(add,'%i');
        if size(locs_n,1)>size(locs_n,2)
            locs_n=locs_n';
            pks_n=pks_n';
        end
        locs_n_old=locs_n;
        pks_n_old=pks_n;
        for i=1:length(add_vec)
            locs_n=[locs_n,add_vec(i)-start_l+1];
            pks_n=[pks_n,Absorption(locs_n(end))];
        end
        [~,b]=sort(pks_n,'descend');
        locs_n=locs_n(b);
        pks_n=pks_n(b);
        
        % Check if assignment would result in any error or inconsistency.
        
        lok=zeros(length(locs_n),1);
        pik=zeros(length(locs_n),1);
        
        for j=1:length(locs_n)
            peak=repmat(locs_n(j)+start_l-1,length(vec),1);
            min_dist=peak-S_n;
            [~,bb]=min(abs(min_dist));
            if min_dist(bb)<-10 && bb>1
                bb=bb-1;
            end
            pik(j)=min_dist(bb);
            lok(j)=bb;
        end
        [n, bin] = histc(lok, unique(lok));
        multiple = find(n > 1);
        
        if isempty(find(ismember(bin, multiple), 1))==0
            warning(['This assignment would result in one SWCNT being assigned to multiple peaks: ',num2str(find(ismember(bin, multiple))'),'. The original vector is restored and you can continue.'])
            locs_n=locs_n_old;
            pks_n=pks_n_old;
        else
            pks_shift=[];
            for j=1:length(pik)
                if pik(j)>max_shift
                    pks_shift=[pks_shift,j];
                end
            end
            if isempty(pks_shift)==0
                warning(['This assignment would result in SWCNT transitions: ',num2str(pks_shift),' being shifted more than ',num2str(max_shift),' nm. This migth be unphysical!'])
                w_q=input('Do you want to restore (1) the old peak assignment or continue (2) with the check of (n,m) assignment?\n','s');
                w_q=sscanf(w_q,'%i');
                if w_q==1
                    locs_n=locs_n_old;
                    pks_n=pks_n_old;
                end
            end
            added_cnts=[added_cnts,add_vec'];
        end
        
    else
        lok=zeros(length(locs_n),1);
        pik=zeros(length(locs_n),1);
        
        for j=1:length(locs_n)
            peak=repmat(locs_n(j)+start_l-1,length(vec),1);
            min_dist=peak-S_n;
            [~,bb]=min(abs(min_dist));
            if min_dist(bb)<-10 && bb>1
                bb=bb-1;
            end
            pik(j)=min_dist(bb);
            lok(j)=bb;
        end
        [n, bin] = histc(lok, unique(lok));
        multiple = find(n > 1);
        
        if isempty(find(ismember(bin, multiple), 1))==0
            warning(['This assignment would result in one SWCNT being assigned to multiple peaks: ',num2str(find(ismember(bin, multiple))'),'. The original vector is restored and you can continue.'])
            locs_n=locs_n_old;
            pks_n=pks_n_old;
        else
            pks_shift=[];
            for j=1:length(pik)
                if pik(j)>max_shift
                    pks_shift=[pks_shift,j];
                end
            end
            if isempty(pks_shift)==1
                close(tt)
                break
            else
                warning(['This assignment would result in SWCNT transitions: ',num2str(pks_shift),' being shifted more than ',num2str(max_shift),' nm. This migth be unphysical!'])
                w_q=input('Do you want to restore (1) the old peak assignment or continue (2) with the check of (n,m) assignment?\n','s');
                w_q=sscanf(w_q,'%i');
                if w_q==2
                    close(tt)
                    break
                else
                    locs_n=locs_n_old;
                    pks_n=pks_n_old;
                end
            end
        end
    end
    
    close(tt)
end

% Peaks that were assigned to SWCNTs are plotted and the user is given the
% option to change this assignment. In the event that the user specifies
% an (n,m) species outside the wavelength range the code returns a warning.

change_cnts=[];

while 1>0
    
    tt=figure('units','normalized','outerposition',[0 0 1 1]); hold on
    plot(lambda,Absorption,'k','linewidth',2)
    plot(locs_n+start_l-1,pks_n,'o','linewidth',2,'MarkerSize',10)
    for i=1:length(locs_n)
        text(double(locs_n(i)+start_l-2),double(pks_n(i)*0.8),num2str(i),'color',[0 .4 0],'fontsize',20)
        text(double(locs_n(i)+start_l+2),double(pks_n(i)),Names{vec(lok(i))},'color','k','fontsize',20)
    end
    
    xlabel('Wavelength (nm)','FontSize',24)
    ylabel('Absorption (a.u.)','FontSize',24)
    legend('Measured Data')
    title('Peak Assignment to SWCNTs','FontSize',26)
    set(gca,'LineWidth',2)
    set(gca,'FontSize',22)
    
    opt=input('Do you want to change the (n,m) assignment of one or more SWCNTs (y/n)? ','s');
    
    if opt=='n'
        close(tt)
        break
    else
        
        opt=input('Enter number of peak(s) you would like to change, followed\nby the correct "n" and "m" of the SWCNT: e.g. peak 1 was incorrectly assigned as (8,6) \nand should be (6,5). Then you have to enter: 1 6 5 \nYour entry: ','s');
        ww=sscanf(opt,'%i');
        ww=ww';
        close(tt)
        
        for i=1:length(ww)/3
            [li,lp(i)]=ismember([ww(3*i-1),ww(3*i)],CNT_n, 'rows');
            if li==0
                warning('CNT not found!')
            else
                lok(ww(1+3*(i-1)))=lp(i);
                change_cnts=[change_cnts,ww(1+3*(i-1):3*i)];
            end
        end
    end
end

% Plot the final peak assignment and give the user the option to assign an
% exciton phonon sideband to the peaks specified by "locs_n" and "pks_n".
% The peaks assigned with an EPS are stored in the new vector "l_g".

opt=input('Do you want to consider exciton phonon sidebands (EPS) in your fit (y/n)? ','s');

if opt == 'y'
    
    tt=figure('units','normalized','outerposition',[0 0 1 1]); hold on
    plot(lambda,Absorption,'k','linewidth',2)
    plot(locs_n+start_l-1,pks_n,'o','linewidth',2,'MarkerSize',10)
    for i=1:length(locs_n)
        text(double(locs_n(i)+start_l-2),double(pks_n(i)*0.8),num2str(i),'color',[0 .4 0],'fontsize',20)
        text(double(locs_n(i)+start_l+2),double(pks_n(i)),Names{vec(lok(i))},'color','k','fontsize',20)
    end
    
    xlabel('Wavelength (nm)','FontSize',24)
    ylabel('Absorption (a.u.)','FontSize',24)
    legend('Measured Data')
    title('Inclusion of Exciton Phonon Sideband(s)','FontSize',26)
    set(gca,'LineWidth',2)
    set(gca,'FontSize',22)
    
    opt=input('Enter the peak number(s) of the SWCNT(s), which should have an exciton phonon sideband (EPS),\n e.g. 1 2 3. Enter 0 if all SWCNT(s) under consideration should have an EPS: ','s');
    
    w2=sscanf(opt,'%i');
    
    close(tt)
    
    if w2==0
        
        l_g=zeros(1,length(S_n));
        
        for i=1:length(S_n)
            l_g(i)=i;
        end
        
    else
        
        l_g=zeros(1,length(w2));
        
        for i=1:length(w2)
            l_g(i)=lok(w2(i));
        end
    end
    
else
    l_g=[];
    w2=[];
end

% Using the vector "Shift" the peak positions stored in the reference data
% table are arbitrarily shifted to match the absorption spectrum.

Shift=zeros(1,length(vec));

for i=1:length(locs_n)
    Shift(i)=locs_n(i)+start_l-1-S_n(lok(i));
end

% Order CNT peaks based on the absorption intensity with unassigned peaks
% listed last and store this in vector "Order".

Order(1,1:length(lok))=lok;

for i=1:length(vec)
    if isempty(find(i==Order, 1))==1
        Order=[Order,i];
    end
end

if w2==0
    l_g=Order;
end

% Modification of vectors based on the peak intensities stored in "Order".

n=1;
switch fit_region
    case 1
        center=S_n(Order')+Shift';
        center_22=S22n(Order');
        Height=zeros(length(center),1);
        for i=1:length(center)
            Height(i)=Absorption(lambda==center(i));
        end
        
        if protocol=='y'
            array{1}=['Removed peaks from S11 region: ',num2str(rid_cnts)];
            array{2}=['Added peak(s) to S11 region at (in nm): ',num2str(added_cnts)];
            array{3}=['Changed peak(s) and correct nanotube(s) for S11 region: ',num2str(change_cnts)];
            array{4}=['Peak(s) with EPS in S11 region: ',num2str(w2')];
            
            fprintf(fid,'%s \n',array{:});
        end
    case 2
        if isempty(S11n)==0
            center=S11n(Order');
        else
            center=[];
        end
        center_22=S_n(Order')+Shift';
        Height=zeros(length(center_22),1);
        for i=1:length(center_22)
            pos22{i}=find(lambda==center_22(i));
            if isempty(pos22{i})==0
                Height(i)=Absorption(pos22{i});
            else
                vec2(n)=i;
                n=n+1;
                warning(['S22 transition of SWCNT ',Names{vec(Order(i))},' is outside the user specified range and will not be considered during the fit.'])
            end
        end
        if n>1
            center_22(vec2)=0;
        end
        
        if protocol=='y'
            array{1}=['Removed peaks from S22 region: ',num2str(rid_cnts)];
            array{2}=['Added peak(s) to S22 region at (in nm): ',num2str(added_cnts)];
            array{3}=['Changed peak(s) and correct nanotube(s) for S22 region: ',num2str(change_cnts)];
            array{4}=['Peak(s) with EPS in S22 region: ',num2str(w2')];
            
            fprintf(fid,'%s \n',array{:});
        end
end

% Definition of upper and lower boundary limits for the shift of the peak
% position during the fit based on the user defined variables.

up_bound=ones(1,length(Order))*upper_lim;
up_bound(1:length(lok))=upper_lim_peak; % The peaks that were found should be maintained

low_bound=zeros(1,length(Order))+lower_lim;
low_bound(1:length(lok))=lower_lim_peak; % The peaks that were found should be maintained


%..........................................................................


%--------------------------------------------------------------------------
% Section 6 - "fit_data"
%--------------------------------------------------------------------------


function [x_fit,sse,cc,cc2,w,s1,e1,PeakAssignment,c_er,Order2,l_g2,height2_22,up_bound2,low_bound2,Order,l_g,center,height2,up_bound,low_bound]=fit_data(max_shift,data_voigt,data_peakfit,data_EPS_peakfit,data_heightratio,Turnover,PeakAssignment,fit_region,lambda,Absorption,height2,Order,low_bound,up_bound,Names,vec,center,fwhm_nm_n,diam_n,l_g,fwhm_nm_n22,complete,start,End,h,c,S22n,S11n,CNT_n,upper_lim,upper_lim_peak,lower_lim,lower_lim_peak,c_er,Order2,l_g2,height2_22,up_bound2,low_bound2,w,fid,protocol,doping,parallel_fit)
global center_22 fwhm_nm_v fwhm_nm_v22

% The user is asked to specify the line-shape that shall be fitted to the
% (n,m) species. If the user specified a Lorentzian or Gaussian line-shape
% ("1" or "2") the subfunction "generate_string", located in Section 6.1 is
% used to create a function handle for the calculation of the residuals of
% the fit. If the user specified a Voigtian line-shape ("3") the
% subfunction "generate_startvalues_voigt", located in Section 6.2 is used
% to create the initial starting values for the fit of the absorption
% spectrum with Voigtians. The fit of Voigtian line-shapes itself is
% performed in Section 6.3 with the functions "Voigt" and
% "complexErrorFunction".

if isempty(w) == 1 % Choose line profile only in the very first round of the fit
    w=input('How do you want to fit your nanotube transitions (exciton phonon sidebands are always fitted with a Gaussian line-shape)?:\n(1) - Lorentzian\n(2) - Gaussian\n(3) - Voigt\n ','s');
    w=sscanf(w,'%i');
end

if PeakAssignment(1)==1 % Redo peak assignment for S11
    [Order,l_g,center,center_22,height2,up_bound,low_bound]=check_peaks(max_shift,Absorption(start:End),lambda(start:End),Names,vec,CNT_n,lambda(start),S11n,S22n,fit_region,upper_lim,upper_lim_peak,lower_lim,lower_lim_peak,fid,protocol);
    PeakAssignment(1)=0;
end

bll=find(low_bound==(low_bound(1)));
low_bound(bll(1):bll(end))=data_peakfit(2,1);
low_bound(bll(end)+1:end)=data_peakfit(3,1);
up_bound(bll(1):bll(end))=data_peakfit(2,3);
up_bound(bll(end)+1:end)=data_peakfit(3,3);

switch w
    case {1,2}
        if parallel_fit=='y'
            try
                options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','Display','off','UseParallel',true);
            catch
                options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','Display','off');
            end
        else
            options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','Display','off');
        end
        
        switch fit_region
            case 1 % Construct function handle for S11 region
                if isempty(complete)==1 || complete == 'n'
                    tic
                    [plotString,lb,ub,x0,cc,cc2]=generate_string(data_peakfit,data_EPS_peakfit,data_heightratio,Turnover,lambda(start:End),Absorption(start:End),height2,Order,low_bound,up_bound,center,fwhm_nm_n,diam_n,l_g,w,[],[],[],[],[],[],h,c,doping);
                    s1=start; e1=End;
                else % Construct function handle for entire region
                    if c_er==0 || PeakAssignment(2)==1
                        start2=find(lambda==round(1e9/(1/min(S22n)*1e9+.2/(h*c)))-10);
                        if isempty(start2)==1 % make sure, that starting wavelength is inside the user defined wavelength range
                            start2=1;
                            start_l=lambda(1);
                        else
                            start_l=round(1e9/(1/min(S22n)*1e9+.2/(h*c)))-10;
                        end
                        End2=find(lambda==max(S22n)+10);
                        [Order2,l_g2,~,center_22,height2_22,up_bound2,low_bound2]=check_peaks(max_shift,Absorption(start2:End2),lambda(start2:End2),Names,vec,CNT_n,start_l,[],S22n,2,upper_lim,upper_lim_peak,lower_lim,lower_lim_peak,fid,protocol);
                        for i=1:length(Order)
                            Order3(i)=find(Order2==Order(i));
                        end
                        center_22=center_22(Order3);
                        
                        PeakAssignment(2)=0;
                        c_er=1;
                    end
                    
                    tic
                    
                    [plotString,lb,ub,x0,cc,cc2]=generate_string(data_peakfit,data_EPS_peakfit,data_heightratio,Turnover,lambda,Absorption,height2,Order,low_bound,up_bound,center,fwhm_nm_n,diam_n,l_g,w,fwhm_nm_n22,height2_22,center_22,low_bound2,up_bound2,l_g2,h,c,doping);
                    s1=1; e1=length(lambda);
                end
            case 2 % Construct function handle for S22 region
                if PeakAssignment(2)==1 % Redo peak assignment for S11
                    [Order,l_g,center,center_22,height2,up_bound,low_bound]=check_peaks(max_shift,Absorption(start:End),lambda(start:End),Names,vec,CNT_n,lambda(start),S11n,S22n,fit_region,upper_lim,upper_lim_peak,lower_lim,lower_lim_peak,fid,protocol);
                    PeakAssignment(2)=0;
                end
                bll=find(low_bound==(low_bound(1)));
                low_bound(bll(1):bll(end))=data_peakfit(2,1);
                low_bound(bll(end)+1:end)=data_peakfit(3,1);
                up_bound(bll(1):bll(end))=data_peakfit(2,3);
                up_bound(bll(end)+1:end)=data_peakfit(3,3);
                
                tic
                
                [plotString,lb,ub,x0,cc,cc2]=generate_string(data_peakfit,data_EPS_peakfit,data_heightratio,Turnover,lambda(start:End),Absorption(start:End),[],Order,low_bound,up_bound,center,[],diam_n,l_g,w,fwhm_nm_n22,height2,center_22,[],[],[],h,c,doping);
                s1=start; e1=End;
        end
        
        % Convert the string "plotString" into an actual function handle
        % "Fit_sum" that calculates the residuals of the fit. The fit is
        % then performed with "lsqnonlin" and constrained with "lb" and
        % "ub". The initial guesses are stored in "x0". Finally, the
        % quality of the fit is calculated by the normalized sum of squared
        % errors (sse).
        
        Fit_sum=str2func(plotString);
        x_fit=lsqnonlin(Fit_sum,abs(x0),double(lb),double(ub),options);
        sse=sum(Fit_sum(x_fit).^2)/sum((Absorption(s1:e1)-mean(Absorption(s1:e1))).^2);
        
        toc
    case 3
        if parallel_fit=='y'
            try
                options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','Display','off','FinDiffRelStep',1e-3,'UseParallel',true);
            catch
                options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','Display','off','FinDiffRelStep',1e-3);
            end
        else
            options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','Display','off','FinDiffRelStep',1e-3);
        end
        
        fwhm_nm_v=fwhm_nm_n(Order)/2; % convert the FWHM of S11 to a HWHM and make it a global variable to parse it easily to fit and plot functions
        fwhm_nm_v22=fwhm_nm_n22(Order)/2; % the same for FWHM of S22
        
        switch fit_region
            case 1 % Fit S11 region
                if isempty(complete)==1 || complete == 'n'
                    tic
                    [lb,ub,x0,cc,cc2]=generate_startvalues_voigt(data_voigt,height2,Order,low_bound,up_bound,center,l_g,data_EPS_peakfit,[],[],[],[],[],[],lambda(start:End),Absorption(start:End),Turnover,data_peakfit,diam_n(Order),h,c,doping);
                    x_fit=lsqnonlin(@(x) Voigt(x,fwhm_nm_v,[],lambda(start:End),Absorption(start:End),[],cc,doping,diam_n,Order,data_EPS_peakfit(1,2),center,[]),double(x0),double(lb),double(ub),options);                    
                    sse=sum(Voigt(x_fit,fwhm_nm_v,[],lambda(start:End),Absorption(start:End),[],cc,doping,diam_n,Order,data_EPS_peakfit(1,2),center,[]).^2)/sum((Absorption(start:End)-mean(Absorption(start:End))).^2);                    
                    toc
                    s1=start; e1=End;
                else % Fit entire region
                    if c_er==0 || PeakAssignment(2)==1
                        start2=find(lambda==round(1e9/(1/min(S22n)*1e9+.2/(h*c)))-40);
                        if isempty(start2)==1 % make sure, that starting wavelength is inside the user defined wavelength range
                            start2=1;
                            start_l=lambda(1);
                        else
                            start_l=round(1e9/(1/min(S22n)*1e9+.2/(h*c)))-40;
                        end
                        End2=find(lambda==max(S22n)+40);
                        
                        [Order2,l_g2,~,center_22,height2_22,up_bound2,low_bound2]=check_peaks(max_shift,Absorption(start2:End2),lambda(start2:End2),Names,vec,CNT_n,start_l,[],S22n,2,upper_lim,upper_lim_peak,lower_lim,lower_lim_peak,fid,protocol);
                        for i=1:length(Order)
                            Order3(i)=find(Order2==Order(i));
                        end
                        center_22=center_22(Order3);
                        
                        PeakAssignment(2)=0;
                        c_er=1;
                    end
                    
                    tic
                    [lb,ub,x0,cc,cc2]=generate_startvalues_voigt(data_voigt,height2,Order,low_bound,up_bound,center,l_g,data_EPS_peakfit,low_bound2,up_bound2,data_heightratio,l_g2,height2_22,center_22,lambda,Absorption,Turnover,data_peakfit,diam_n(Order),h,c,doping);
                    x_fit=lsqnonlin(@(x) Voigt(x,fwhm_nm_v,fwhm_nm_v22,lambda,Absorption,cc2,cc,doping,diam_n,Order,data_EPS_peakfit(1,2),center,center_22),double(x0),double(lb),double(ub),options);
                    sse=sum(Voigt(x_fit,fwhm_nm_v,fwhm_nm_v22,lambda,Absorption,cc2,cc,doping,diam_n,Order,data_EPS_peakfit(1,2),center,center_22).^2)/sum((Absorption-mean(Absorption)).^2);
                    toc
                    s1=1; e1=length(lambda);
                end
            case 2 % Fit S22 region
                if PeakAssignment(2)==1 % Redo peak assignment for S11
                    [Order,l_g,center,center_22,height2,up_bound,low_bound]=check_peaks(max_shift,Absorption(start:End),lambda(start:End),Names,vec,CNT_n,lambda(start),S11n,S22n,fit_region,upper_lim,upper_lim_peak,lower_lim,lower_lim_peak,fid,protocol);
                    PeakAssignment(2)=0;
                end
                bll=find(low_bound==(low_bound(1)));
                low_bound(bll(1):bll(end))=data_peakfit(2,1);
                low_bound(bll(end)+1:end)=data_peakfit(3,1);
                up_bound(bll(1):bll(end))=data_peakfit(2,3);
                up_bound(bll(end)+1:end)=data_peakfit(3,3);
                
                tic
                [lb,ub,x0,cc,cc2]=generate_startvalues_voigt(data_voigt,[],Order,low_bound,up_bound,center,l_g,data_EPS_peakfit,[],[],[],[],height2,center_22,lambda(start:End),Absorption(start:End),Turnover,data_peakfit,diam_n(Order),h,c,doping);
                x_fit=lsqnonlin(@(x) Voigt(x,[],fwhm_nm_v22,lambda(start:End),Absorption(start:End),[],cc,doping,diam_n,Order,data_EPS_peakfit(1,2),center,center_22),double(x0),double(lb),double(ub),options);
                sse=sum(Voigt(x_fit,[],fwhm_nm_v22,lambda(start:End),Absorption(start:End),[],cc,doping,diam_n,Order,data_EPS_peakfit(1,2),center,center_22).^2)/sum((Absorption(start:End)-mean(Absorption(start:End))).^2);
                toc
                s1=start; e1=End;
        end
end


%--------------------------------------------------------------------------
% Section 6.1 - "generate_string"
%--------------------------------------------------------------------------


function [plotString,lb,ub,x0,cc,cc2]=generate_string(data_peakfit,data_EPS_peakfit,data_heightratio,Turnover,lambda,Absorption,height_initn,Order,low_bound,up_bound,center,fwhm_nm_n,diam_n,l_g,w,fwhm_nm_n22,height2_22,center_22,low_bound2,up_bound2,l_g2,h,c,doping)

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


% Determine the region to be fitted based on the user choice stored in
% "fit_region".

complete=0; % Fit entire spectrum with S11 and S22

if isempty(height_initn)==1 % Fit only S22
    Cent=center_22;
    Height=height2_22;
    FWHM=fwhm_nm_n22;
    complete=1;
end

if isempty(height2_22)==1 % Fit only S11
    Cent=center;
    Height=height_initn;
    FWHM=fwhm_nm_n;
    complete=1;
end

if complete==0 % Fit entire spectrum with S11 and S22
    Cent=center;
    Height=height_initn;
    FWHM=fwhm_nm_n;
end

% Store constraints of the fit in the vectors "lb" and "ub" for the lower,
% upper boundary conditions, respectively. The initial starting values for
% the non-linear fit are stored in "x0". The initial Lorentzian or Gaussian
% line-shape is strored in "y2" and is needed if "height_guess" = 2. The
% initial EPS are stored in "y_gn" and are needed if any phonon sidebands
% were assigned by the user and "height_guess" = 2.

lb=[];
ub=[];
x0=[];
plotString='';
y2=zeros(length(lambda),length(Cent)); % Store all Lorentzian/Gaussian initial fits in this vector
y_gn=zeros(length(lambda),length(Cent));

% Creation of function handle for the fit of S11 or S22 (as defined by the
% user). Each peak is defined by three components that can vary during
% the fit: x(1)=height, x(2)=Shift in x-direction, x(3)=broadening of FWHM.

for j=1:length(Cent)
    
    % Define the line-shape of the most intense peak and store it in
    % "plotString". Fill the vectors "lb", "ub" and "x0" with lower/upper
    % boundary limit constraints for this first peak and an initial guess
    % of its height (from "check_peaks"), shift of center position and FWHM
    % broadening. If method (2) of guessing the height for an optimal fit
    % was chosen, calculate the line profile based on the initial values
    % for height and FWHM and subtract it from the absorption spectrum. If
    % any exciton phonon was specified by the user, calculate its
    % line-shape and also subtract it from the absorption spectrum.
    
    if j==1
        if w==1 % Lorentzian
            plotString = strcat(plotString, ['@(x) double([',num2str(Absorption.',' %f32'),'] -x(',num2str(j*3-2),')./(1+(([',num2str(lambda(1:end).',' %f32'),']-',num2str(Cent(j)),'-x(',num2str(j*3-1),'))/(0.5*',num2str(FWHM(Order(j))),'*x(',num2str(j*3),'))).^2)']);
            y2(:,1)=Height(j) ./ (1 + ((lambda-Cent(j))/(0.5*FWHM(Order(j)))).^2);
            y_gn(:,1)=Height(j)*pi*FWHM(Order(j))/2*(0.017+0.1/diam_n(Order(j))+data_EPS_peakfit(2,2))*exp(-2*(lambda-(h*c/(h*c/Cent(j)*1e9+0.2))*1e9).^2/(data_EPS_peakfit(1,2)/sqrt(log(4)))^2)/(data_EPS_peakfit(1,2)*sqrt(pi/(2*log(4))));
        else    % Gaussian
            plotString = strcat(plotString, ['@(x) double([',num2str(Absorption.',' %f32'),'] -x(',num2str(j*3-2),')*exp(-log(2)*(([',num2str(lambda(1:end).',' %f32'),']-',num2str(Cent(j)),'-x(',num2str(j*3-1),'))/(',num2str(FWHM(Order(j))),'*x(',num2str(j*3),')/2)).^2)']);
            y2(:,1)=Height(j)*exp(-log(2)*((lambda-Cent(j))/(FWHM(Order(j))/2)).^2);
            y_gn(:,1)=Height(j)* FWHM(Order(j)) * sqrt(pi/(2*log(4)))*(0.017+0.1/diam_n(Order(j))+data_EPS_peakfit(2,2))*exp(-2*(lambda-(h*c/(h*c/Cent(j)*1e9+0.2))*1e9).^2/(data_EPS_peakfit(1,2)/sqrt(log(4)))^2)/(data_EPS_peakfit(1,2)*sqrt(pi/(2*log(4))));
        end
        
        lb=[Height(j)*data_peakfit(4,1),low_bound(j),data_peakfit(1,1)];
        ub=[Height(j)*data_peakfit(4,3),up_bound(j),data_peakfit(1,3)];
        x0=[Height(j)*data_peakfit(4,2),0,data_peakfit(1,2)];
        
        if Turnover(4,1) == 2
            if isempty(l_g)==1 || complete==0
                Absorption=Absorption-y2(:,1);
            else
                Absorption=Absorption-y2(:,1)-y_gn(:,1);
            end
        end
        
    else
        
        % Repeat the process outlined for the most intense peak for any
        % other peak in set of (n,m) species defined by the user.
        
        if w==1 % Lorentzian
            plotString = strcat(plotString, ['-x(',num2str(j*3-2),')./(1+(([',num2str(lambda(1:end).',' %f32'),']-',num2str(Cent(j)),'-x(',num2str(j*3-1),'))/(0.5*',num2str(FWHM(Order(j))),'*x(',num2str(j*3),'))).^2)']);
            
            y_gn(:,j)=Height(j)*pi*FWHM(Order(j))/2*(0.017+0.1/diam_n(Order(j))+data_EPS_peakfit(2,2))*exp(-2*(lambda-(h*c/(h*c/Cent(j)*1e9+0.2))*1e9).^2/(data_EPS_peakfit(1,2)/sqrt(log(4)))^2)/(data_EPS_peakfit(1,2)*sqrt(pi/(2*log(4))));
            
            if Turnover(4,1) == 2
                if isempty(l_g)==1 || complete==0
                    y2(:,j)=Absorption(lambda==Cent(j)) ./ (1 + ((lambda-Cent(j))/(0.5*FWHM(Order(j)))).^2);
                else
                    y2(:,j)=Absorption(lambda==Cent(j)) ./ (1 + ((lambda-Cent(j))/(0.5*FWHM(Order(j)))).^2)+y_gn(:,j);
                end
            end
        else % Gaussian
            plotString = strcat(plotString, ['-x(',num2str(j*3-2),')*exp(-log(2)*(([',num2str(lambda(1:end).',' %f32'),']-',num2str(Cent(j)),'-x(',num2str(j*3-1),'))/(',num2str(FWHM(Order(j))),'*x(',num2str(j*3),')/2)).^2)']);
            
            y_gn(:,j)=Height(j)* FWHM(Order(j)) * sqrt(pi/(2*log(4)))*(0.017+0.1/diam_n(Order(j))+data_EPS_peakfit(2,2))*exp(-2*(lambda-(h*c/(h*c/Cent(j)*1e9+0.2))*1e9).^2/(data_EPS_peakfit(1,2)/sqrt(log(4)))^2)/(data_EPS_peakfit(1,2)*sqrt(pi/(2*log(4))));
            
            if Turnover(4,1)==2
                if isempty(l_g)==1 || complete==0
                    y2(:,j)=Absorption(lambda==Cent(j))*exp(-log(2)*((lambda-Cent(j))/(FWHM(Order(j))/2)).^2);
                else
                    y2(:,j)=Absorption(lambda==Cent(j))*exp(-log(2)*((lambda-Cent(j))/(FWHM(Order(j))/2)).^2)+y_gn(:,Order(j));
                end
            end
        end
        
        % The height guess and broadening of the FWHM is also dependent on
        % the center position of the nanotube peak. If the center position
        % is larger than a user specified "turnover" value, the broadening
        % of FWHM is increased to 1.6. Also the estimation of the initial
        % height is performed based on the inital heights stored in
        % "Height" that were determined in "check_peaks", which represents
        % method (1) of "height_guess".
        
        if Cent(j) > Turnover(1,1)
            
            lb=[lb,double(Height(Order(j)))*data_peakfit(5,1),low_bound(j),Turnover(3,1)];
            
            if double(abs(Absorption(lambda==Cent(j))))==0 % Make sure that upper bound is larger than zero
                ub=[ub,0.01,up_bound(j),Turnover(2,1)];
            else
                ub=[ub,double(Height(j))*data_peakfit(5,3),up_bound(j),Turnover(2,1)];
            end
            
            x0=[x0,double(Height(j))*data_peakfit(5,2),0,1];
            
        else
            
            lb=[lb,double(abs(Absorption(lambda==Cent(j))))*data_peakfit(6,1),low_bound(j),data_peakfit(1,1)];
            
            if double(abs(Absorption(lambda==Cent(j))))==0
                ub=[ub,0.01,up_bound(j),data_peakfit(1,3)];
            else
                ub=[ub,double(abs(Absorption(lambda==Cent(j))))*data_peakfit(6,3),up_bound(j),data_peakfit(1,3)];
            end
            
            x0=[x0,double(abs(Absorption(lambda==Cent(j))))*data_peakfit(6,2),0,data_peakfit(1,2)];
            
        end
        
        Absorption=Absorption-y2(:,j);
        
    end
end

% Creation of a function handle to fit the S22 region after the function
% handle for the fit of the S11 region was created.

if complete==0
    n=0;
    for j=1:length(center_22)
        if center_22(j) > 0
            if j==1
                if w==1 % Lorentzian
                    plotString = strcat(plotString, ['-x(',num2str(j*3-2-3*n),')/(x(',num2str(length(center)*3+j*3-2-3*n),'))./(1+(([',num2str(lambda(1:end).',' %f32'),']-',num2str(center_22(j)),'-x(',num2str(length(center)*3+j*3-1-3*n),'))/(0.5*',num2str(fwhm_nm_n22(Order(j))),'*x(',num2str(length(center)*3+j*3-3*n),'))).^2)']);
                else    % Gaussian
                    plotString = strcat(plotString, ['-x(',num2str(j*3-2-3*n),')/(x(',num2str(length(center)*3+j*3-2-3*n),'))*exp(-log(2)*(([',num2str(lambda(1:end).',' %f32'),']-',num2str(center_22(j)),'-x(',num2str(length(center)*3+j*3-1-3*n),'))/(',num2str(fwhm_nm_n22(Order(j))),'*x(',num2str(length(center)*3+j*3-3*n),')/2)).^2)']);
                end
                
                if doping==0 % no doping effects
                    if data_peakfit(1,1)<data_peakfit(1,1)*FWHM(Order(j))/fwhm_nm_n22(Order(j))
                        lb=[lb,data_heightratio(1,1),low_bound2(j),data_peakfit(1,1)];
                        ub=[ub,data_heightratio(1,3),up_bound2(j),data_peakfit(1,1)*FWHM(Order(j))/fwhm_nm_n22(Order(j))];
                    else
                        lb=[lb,data_heightratio(1,1),low_bound2(j),data_peakfit(1,1)*FWHM(Order(j))/fwhm_nm_n22(Order(j))];
                        ub=[ub,data_heightratio(1,3),up_bound2(j),data_peakfit(1,1)];
                    end
                    x0=[x0,data_heightratio(1,2),0,data_peakfit(1,2)];
                else
                    lb=[lb,1e-5,low_bound2(j),data_peakfit(1,1)];
                    ub=[ub,data_heightratio(1,3),up_bound2(j),data_peakfit(1,3)];
                    x0=[x0,data_heightratio(1,2),0,data_peakfit(1,2)];
                end
            else
                if w==1 % Lorentzian
                    plotString = strcat(plotString, ['-x(',num2str(j*3-2-3*n),')/(x(',num2str(length(center)*3+1-3*n),')*x(',num2str(length(center)*3+3*j-2-3*n),'))./(1+(([',num2str(lambda(1:end).',' %f32'),']-',num2str(center_22(j)),'-x(',num2str(length(center)*3+j*3-1-3*n),'))/(0.5*',num2str(fwhm_nm_n22(Order(j))),'*x(',num2str(length(center)*3+j*3-3*n),'))).^2)']);
                else    % Gaussian
                    plotString = strcat(plotString, ['-x(',num2str(j*3-2-3*n),')/(x(',num2str(length(center)*3+1-3*n),')*x(',num2str(length(center)*3+3*j-2-3*n),'))*exp(-log(2)*(([',num2str(lambda(1:end).',' %f32'),']-',num2str(center_22(j)),'-x(',num2str(length(center)*3+j*3-1-3*n),'))/(',num2str(fwhm_nm_n22(Order(j))),'*x(',num2str(length(center)*3+j*3-3*n),')/2)).^2)']);
                end
                
                if doping==0 % no doping effects
                    if data_peakfit(1,1)<data_peakfit(1,1)*FWHM(Order(j))/fwhm_nm_n22(Order(j))
                        lb=[lb,data_heightratio(2,1),low_bound2(j),data_peakfit(1,1)];
                        ub=[ub,data_heightratio(2,3),up_bound2(j),data_peakfit(1,1)*FWHM(Order(j))/fwhm_nm_n22(Order(j))];
                    else
                        lb=[lb,data_heightratio(2,1),low_bound2(j),data_peakfit(1,1)*FWHM(Order(j))/fwhm_nm_n22(Order(j))];
                        ub=[ub,data_heightratio(2,3),up_bound2(j),data_peakfit(1,1)];
                    end
                    x0=[x0,data_heightratio(2,2),0,data_peakfit(1,2)];
                else
                    lb=[lb,data_heightratio(2,1),low_bound2(j),data_peakfit(1,1)];
                    ub=[ub,data_heightratio(2,3),up_bound2(j),data_peakfit(1,3)];
                    x0=[x0,data_heightratio(2,2),0,data_peakfit(1,2)];
                end
            end
        else
            n=n+1;
        end
    end
end

% Creation of function handle for the fit of an exciton phonon sideband for
% the S11 or S22 region.

x0_len=length(x0);
cc=[]; cc2=[];

for j=1:length(l_g)
    cc(j)=find(Order==l_g(j));
    
    % Fit parameters:
    % x(1)=shift of diam vs. area ratio curve from Perebeinos et al.
    % PRL 94, 027402 (2005)
    % x(2)=shift of center around 0.2 eV above associated S11 or S22 peak
    % x(3)=Broadening of FWHM
    
    if j==1
        
        if w==1 % Lorentzian
            plotString = strcat(plotString,['-x(',num2str(cc(j)*3-2),')*pi*',num2str(FWHM(Order(cc(j)))),'*x(',num2str(cc(j)*3),')/2*(0.017+0.1/',num2str(diam_n(Order(cc(j)))),' + x(',num2str(x0_len+j),'))*exp(-2*([',num2str(lambda(1:end).',' %f32'),']-',num2str(h*c),'/(1e-9*(',num2str(h*c),'/((',num2str(Cent(cc(j))),'+x(',num2str(cc(j)*3-1),'))*1e-9)+0.2+x(',num2str(x0_len+2),')))).^2/(',num2str(data_EPS_peakfit(1,2)),'*x(',num2str(x0_len+3),')/sqrt(log(4)))^2)/((',num2str(data_EPS_peakfit(1,2)),'*x(',num2str(x0_len+3),'))*sqrt(pi/(2*log(4))))']);
        else % Gaussian
            plotString = strcat(plotString,['-x(',num2str(cc(j)*3-2),')*',num2str(FWHM(Order(cc(j)))),'*x(',num2str(cc(j)*3),')*sqrt(pi/(2*log(4)))*(0.017+0.1/',num2str(diam_n(Order(cc(j)))),' + x(',num2str(x0_len+j),'))*exp(-2*([',num2str(lambda(1:end).',' %f32'),']-',num2str(h*c),'/(1e-9*(',num2str(h*c),'/((',num2str(Cent(cc(j))),'+x(',num2str(cc(j)*3-1),'))*1e-9)+0.2+x(',num2str(x0_len+2),')))).^2/(',num2str(data_EPS_peakfit(1,2)),'*x(',num2str(x0_len+3),')/sqrt(log(4)))^2)/((',num2str(data_EPS_peakfit(1,2)),'*x(',num2str(x0_len+3),'))*sqrt(pi/(2*log(4))))']);
        end
        
        lb=[lb,data_EPS_peakfit(2,1),data_EPS_peakfit(3,1),data_EPS_peakfit(1,1)];
        ub=[ub,data_EPS_peakfit(2,3),data_EPS_peakfit(3,3),data_EPS_peakfit(1,3)];
        x0=[x0,data_EPS_peakfit(2,2),data_EPS_peakfit(3,2),1];
        
    else
        
        % Fit parameters:
        % x(2+2*(j-1)) = shift of center around 0.2 eV above
        % associated S11
        % x(3+2*(j-1)) = Percentage of FWHM broadening of first
        % peak - can vary between 90 and 110 %.
        
        if w==1
            plotString = strcat(plotString,['-x(',num2str(cc(j)*3-2),')*pi*',num2str(FWHM(Order(cc(j)))),'*x(',num2str(cc(j)*3),')/2*(0.017+0.1/',num2str(diam_n(Order(cc(j)))),'+x(',num2str(x0_len+1),'))*exp(-2*([',num2str(lambda(1:end).',' %f32'),']-',num2str(h*c),'/(1e-9*(',num2str(h*c),'/((',num2str(Cent(cc(j))),'+x(',num2str(cc(j)*3-1),'))*1e-9)+0.2+x(',num2str(x0_len+4+2*(j-2)),')))).^2/(',num2str(data_EPS_peakfit(1,2)),'*x(',num2str(x0_len+3),')*x(',num2str(x0_len+5+2*(j-2)),')/sqrt(log(4)))^2)/((',num2str(data_EPS_peakfit(1,2)),'*x(',num2str(x0_len+3),')*x(',num2str(x0_len+5+2*(j-2)),'))*sqrt(pi/(2*log(4))))']);
        else
            plotString = strcat(plotString,['-x(',num2str(cc(j)*3-2),')*',num2str(FWHM(Order(cc(j)))),'*x(',num2str(cc(j)*3),')*sqrt(pi/(2*log(4)))*(0.017+0.1/',num2str(diam_n(Order(cc(j)))),'+x(',num2str(x0_len+1),'))*exp(-2*([',num2str(lambda(1:end).',' %f32'),']-',num2str(h*c),'/(1e-9*(',num2str(h*c),'/((',num2str(Cent(cc(j))),'+x(',num2str(cc(j)*3-1),'))*1e-9)+0.2+x(',num2str(x0_len+4+2*(j-2)),')))).^2/(',num2str(data_EPS_peakfit(1,2)),'*x(',num2str(x0_len+3),')*x(',num2str(x0_len+5+2*(j-2)),')/sqrt(log(4)))^2)/((',num2str(data_EPS_peakfit(1,2)),'*x(',num2str(x0_len+3),')*x(',num2str(x0_len+5+2*(j-2)),'))*sqrt(pi/(2*log(4))))']);
        end
        
        lb=[lb,data_EPS_peakfit(3,1),data_EPS_peakfit(4,1)];
        ub=[ub,data_EPS_peakfit(3,3),data_EPS_peakfit(4,3)];
        x0=[x0,data_EPS_peakfit(3,2),data_EPS_peakfit(4,2)];
        
    end
end

% Creation of function handle for the fit of an exciton phonon sideband
% for the S22 region after assigning possible phonon sidebands to the S11
% region.

if complete == 0 && isempty(l_g2)==0
    
    for j=1:length(l_g2)
        if isempty(cc)
            cc_len=0;
        else
            cc_len=3+2*(length(cc)-1);
        end
        
        cc2(j)=find(Order==l_g2(j));
        if j==1
            
            % Fit parameters: as described above
            
            if w==1 % Lorentzian
                plotString = strcat(plotString,['-x(',num2str(cc2(j)*3-2),')/(x(',num2str(3*(length(center))+cc2(j)*3-2),'))*pi*',num2str(fwhm_nm_n22(Order(cc2(j)))),'*x(',num2str((length(center))*3+cc2(j)*3),')/2*(0.017+0.1/',num2str(diam_n(Order(cc2(j)))),' + x(',num2str(6*length(center)-3*n+cc_len+j),'))*exp(-2*([',num2str(lambda(1:end).',' %f32'),']-',num2str(h*c),'/(1e-9*(',num2str(h*c),'/((',num2str(center_22(cc2(j))),'+x(',num2str(3*(length(center))+cc2(j)*3-1),'))*1e-9)+0.2+x(',num2str(6*length(center)-3*n+cc_len+1+j),')))).^2/(',num2str(data_EPS_peakfit(1,2)),'*x(',num2str(6*length(center)-3*n+cc_len+2+j),')/sqrt(log(4)))^2)/((',num2str(data_EPS_peakfit(1,2)),'*x(',num2str(6*length(center)-3*n+cc_len+2+j),'))*sqrt(pi/(2*log(4))))']);
            else % Gaussian
                plotString = strcat(plotString,['-x(',num2str(cc2(j)*3-2),')/(x(',num2str(3*(length(center))+cc2(j)*3-2),'))*',num2str(fwhm_nm_n22(Order(cc2(j)))),'*x(',num2str((length(center))*3+cc2(j)*3),')*sqrt(pi/(2*log(4)))*(0.017+0.1/',num2str(diam_n(Order(cc2(j)))),' + x(',num2str(6*length(center)-3*n+cc_len+j),'))*exp(-2*([',num2str(lambda(1:end).',' %f32'),']-',num2str(h*c),'/(1e-9*(',num2str(h*c),'/((',num2str(center_22(cc2(j))),'+x(',num2str(3*(length(center))+cc2(j)*3-1),'))*1e-9)+0.2+x(',num2str(6*length(center)-3*n+cc_len+1+j),')))).^2/(',num2str(data_EPS_peakfit(1,2)),'*x(',num2str(6*length(center)-3*n+cc_len+2+j),')/sqrt(log(4)))^2)/((',num2str(data_EPS_peakfit(1,2)),'*x(',num2str(6*length(center)-3*n+cc_len+2+j),'))*sqrt(pi/(2*log(4))))']);
            end
            
            lb=[lb,data_EPS_peakfit(2,1),data_EPS_peakfit(3,1),data_EPS_peakfit(1,1)];
            ub=[ub,data_EPS_peakfit(2,3),data_EPS_peakfit(3,3),data_EPS_peakfit(1,3)];
            x0=[x0,data_EPS_peakfit(2,2),data_EPS_peakfit(3,2),1];
            
        else
            
            % Fit parameters: as described above
            
            if w==1
                plotString = strcat(plotString,['-x(',num2str(cc2(j)*3-2),')/(x(',num2str(3*(length(center))+cc2(1)*3-2),')*x(',num2str(3*(length(center))+cc2(j)*3-2),'))*pi*',num2str(fwhm_nm_n22(Order(cc2(j)))),'*x(',num2str((length(center))*3+cc2(j)*3),')/2*(0.017+0.1/',num2str(diam_n(Order(cc2(j)))),'+x(',num2str(6*length(center)-3*n+cc_len+1),'))*exp(-2*([',num2str(lambda(1:end).',' %f32'),']-',num2str(h*c),'/(1e-9*(',num2str(h*c),'/((',num2str(center_22(cc2(j))),'+x(',num2str(3*(length(center))+cc2(j)*3-1),'))*1e-9)+0.2+x(',num2str(6*length(center)-3*n+cc_len+4+2*(j-2)),')))).^2/(',num2str(data_EPS_peakfit(1,2)),'*x(',num2str(6*length(center)-3*n+3+cc_len),')*x(',num2str(6*length(center)-3*n+cc_len+5+2*(j-2)),')/sqrt(log(4)))^2)/((',num2str(data_EPS_peakfit(1,2)),'*x(',num2str(6*length(center)-3*n+3+cc_len),')*x(',num2str(6*length(center)-3*n+cc_len+5+2*(j-2)),'))*sqrt(pi/(2*log(4))))']);
            else
                plotString = strcat(plotString,['-x(',num2str(cc2(j)*3-2),')/(x(',num2str(3*(length(center))+cc2(1)*3-2),')*x(',num2str(3*(length(center))+cc2(j)*3-2),'))*',num2str(fwhm_nm_n22(Order(cc2(j)))),'*x(',num2str((length(center))*3+cc2(j)*3),')*sqrt(pi/(2*log(4)))*(0.017+0.1/',num2str(diam_n(Order(cc2(j)))),'+x(',num2str(6*length(center)-3*n+cc_len+1),'))*exp(-2*([',num2str(lambda(1:end).',' %f32'),']-',num2str(h*c),'/(1e-9*(',num2str(h*c),'/((',num2str(center_22(cc2(j))),'+x(',num2str(3*(length(center))+cc2(j)*3-1),'))*1e-9)+0.2+x(',num2str(6*length(center)-3*n+cc_len+4+2*(j-2)),')))).^2/(',num2str(data_EPS_peakfit(1,2)),'*x(',num2str(6*length(center)-3*n+3+cc_len),')*x(',num2str(6*length(center)-3*n+cc_len+5+2*(j-2)),')/sqrt(log(4)))^2)/((',num2str(data_EPS_peakfit(1,2)),'*x(',num2str(6*length(center)-3*n+3+cc_len),')*x(',num2str(6*length(center)-3*n+cc_len+5+2*(j-2)),'))*sqrt(pi/(2*log(4))))']);
            end
            
            lb=[lb,data_EPS_peakfit(3,1),data_EPS_peakfit(4,1)];
            ub=[ub,data_EPS_peakfit(3,3),data_EPS_peakfit(4,3)];
            x0=[x0,data_EPS_peakfit(3,2),data_EPS_peakfit(4,2)];
            
        end
    end
end

plotString=strcat(plotString,[')']);


%--------------------------------------------------------------------------
% Section 6.2 - "generate_startvalues_voigt"
%--------------------------------------------------------------------------


function [lb,ub,x0,cc,cc2]=generate_startvalues_voigt(data_voigt,height_initn,Order,low_bound,up_bound,center,l_g,data_EPS_peakfit,low_bound2,up_bound2,data_heightratio,l_g2,height2_22,center_22,lambda,Absorption,Turnover,data_peakfit,diam_n,h,c,doping)
global fwhm_nm_v fwhm_nm_v22

% Define the region to be fitted based on the users choice. Calculate the
% initial Lorentzian HWHM and save it as variable "FWHM".

complete = 0;

if isempty(height_initn)==1 % Fit S22
    complete=1;
    Height=height2_22;
    FWHM=fwhm_nm_v22/(0.5436+sqrt(data_voigt(2)^2+.2166));
    Cent=center_22;
end

if isempty(height2_22)==1 % Fit S11
    complete=1;
    Height=height_initn;
    FWHM=fwhm_nm_v/(0.5436+sqrt(data_voigt(2)^2+.2166));
    Cent=center;
end

if complete==0 % Fit entire region
    Height=height_initn;
    FWHM=fwhm_nm_v/(0.5436+sqrt(data_voigt(2)^2+.2166));
    Cent=center;
end

% Start values for S11 or S22. The Voigt line-shape for one (n,m) species
% is defined by the following parameters:
% x(1)=Height
% x(2)=shift of center
% x(3)=Broadening of HWHM
% x(4)=Ratio of Gaussian / Lorentzian HWHM
% The initial guess of these parameters is stored in "x0" and constrained
% during the fit by the upper and lower boundary values stored in "ub" and
% "lb". The initial calculation of the height follows the routine outlined
% in Section 6.1 "generate_string".

lb=[];
ub=[];
x0=[];

for j=1:length(center)
    
    if j==1
        lb=[lb,Height(j)*data_peakfit(4,1),low_bound(j),data_peakfit(1,1),data_voigt(3)];
        ub=[ub,Height(j)*data_peakfit(4,3),up_bound(j),data_peakfit(1,3),data_voigt(4)];
        x0=[x0,Height(j)*data_peakfit(4,2),0,data_peakfit(1,2),data_voigt(2)];
        
        wL(j)=FWHM(j);
        wg(j)=wL(j)*data_voigt(2);
        XX=sqrt(log(2))*(lambda-Cent(j))/wg(j);
        YY=sqrt(log(2))*wL(j)/wg(j);
        y(:,j)=Height(j).*real(complexErrorFunction(XX,YY))./real(complexErrorFunction(0,YY));
        
        if Turnover(4,1) == 2
            if isempty(l_g)==1 || complete==0
                Absorption=Absorption-y(:,1);
            else
                y_g(:,j)=trapz(lambda,y(:,j))*(0.017+0.1/diam_n(j)+data_EPS_peakfit(2,2))*exp(-2*(lambda-h*c/(1e-9*(h*c/(Cent(j)*1e-9)+0.2))).^2/(data_EPS_peakfit(1,2)/sqrt(log(4)))^2)/(data_EPS_peakfit(1,2)*sqrt(pi/(2*log(4))));
                Absorption=Absorption-y(:,j)-y_g(:,j);
            end
        end
    else
        
        if Cent(j) > Turnover(1,1)
            
            lb=[lb,Height(j)*data_peakfit(5,1),low_bound(j),data_peakfit(1,1),1-data_voigt(5)];
            
            if double(abs(Absorption(lambda==Cent(j))))==0 % Make sure that upper bound is larger than zero
                ub=[ub,0.01,up_bound(j),Turnover(2,1),1+data_voigt(5)];
            else
                ub=[ub,Height(j)*data_peakfit(5,3),up_bound(j),Turnover(2,1),1+data_voigt(5)];
            end
            
            x0=[x0,Height(j)*data_peakfit(5,2),0,data_peakfit(1,2),1];
            
        else
            
            lb=[lb,double(abs(Absorption(lambda==Cent(j))))*data_peakfit(6,1),low_bound(j),data_peakfit(1,1),1-data_voigt(5)];
            
            if double(abs(Absorption(lambda==Cent(j))))==0
                ub=[ub,0.01,up_bound(j),data_peakfit(1,3),1+data_voigt(5)];
            else
                ub=[ub,double(abs(Absorption(lambda==Cent(j))))*data_peakfit(6,3),up_bound(j),data_peakfit(1,3),1+data_voigt(5)];
            end
            
            x0=[x0,double(abs(Absorption(lambda==Cent(j))))*data_peakfit(6,2),0,data_peakfit(1,2),1];
            
        end
        
        wL(j)=FWHM(j);
        wg(j)=wL(j)*data_voigt(2);
        XX=sqrt(log(2))*(lambda-Cent(j))/wg(j);
        YY=sqrt(log(2))*wL(j)/wg(j);
        y(:,j)=Absorption(lambda==Cent(j)).*real(complexErrorFunction(XX,YY))./real(complexErrorFunction(0,YY));
        
        if Turnover(4,1) == 2
            if isempty(l_g)==1 || complete==0
                Absorption=Absorption-y(:,j);
            else
                y_g(:,j)=trapz(lambda,y(:,j))*(0.017+0.1/diam_n(j)+data_EPS_peakfit(2,2))*exp(-2*(lambda-h*c/(1e-9*(h*c/(Cent(j)*1e-9)+0.2))).^2/(data_EPS_peakfit(1,2)/sqrt(log(4)))^2)/(data_EPS_peakfit(1,2)*sqrt(pi/(2*log(4))));
                Absorption=Absorption-y(:,j)-y_g(:,j);
            end
        end
        
    end
end

% Calculate the start values and upper and lower boundary conditions for
% the fit of the S22 region, after the values for the S11 region were
% determined.

if complete == 0
    n=0;
    for j=1:length(center_22)
        if center_22(j) > 0
            if j==1
                if doping ==0
                    lb=[lb,data_heightratio(1,1),low_bound2(j),data_peakfit(1,1),data_voigt(3)];
                else
                    lb=[lb,1e-5,low_bound2(j),data_peakfit(1,1),data_voigt(3)];
                end
                ub=[ub,data_heightratio(1,3),up_bound2(j),data_peakfit(1,3),data_voigt(4)];
                x0=[x0,data_heightratio(1,2),0,data_peakfit(1,2),data_voigt(2)];
            else
                if doping ==0
                    lb=[lb,data_heightratio(2,1),low_bound2(j),data_peakfit(1,1),1-data_voigt(5)];
                else
                    lb=[lb,1e-5,low_bound2(j),data_peakfit(1,1),1-data_voigt(5)];
                end
                ub=[ub,data_heightratio(2,3),up_bound2(j),data_peakfit(1,3),1+data_voigt(5)];
                x0=[x0,data_heightratio(2,2),0,data_peakfit(1,2),1];
            end
        else
            n=n+1;
        end
    end
end

% Calculate the starting values for the exciton phonon sideband for S11 or
% S22 region.

cc=[]; cc2=[];

for j=1:length(l_g)
    cc(j)=find(Order==l_g(j));
    if j==1
        lb=[lb,data_EPS_peakfit(2,1),data_EPS_peakfit(3,1),data_EPS_peakfit(1,1)];
        ub=[ub,data_EPS_peakfit(2,3),data_EPS_peakfit(3,3),data_EPS_peakfit(1,3)];
        x0=[x0,data_EPS_peakfit(2,2),data_EPS_peakfit(3,2),1];
    else
        lb=[lb,data_EPS_peakfit(3,1),data_EPS_peakfit(4,1)];
        ub=[ub,data_EPS_peakfit(3,3),data_EPS_peakfit(4,3)];
        x0=[x0,data_EPS_peakfit(3,2),data_EPS_peakfit(4,2)];
    end
end

if complete==0 && isempty(l_g2)==0
    
    for j=1:length(l_g2)
        cc2(j)=find(Order==l_g2(j));
        if j==1
            lb=[lb,data_EPS_peakfit(2,1),data_EPS_peakfit(3,1),data_EPS_peakfit(1,1)];
            ub=[ub,data_EPS_peakfit(2,3),data_EPS_peakfit(3,3),data_EPS_peakfit(1,3)];
            x0=[x0,data_EPS_peakfit(2,2),data_EPS_peakfit(3,2),1];
        else
            lb=[lb,data_EPS_peakfit(3,1),data_EPS_peakfit(4,1)];
            ub=[ub,data_EPS_peakfit(3,3),data_EPS_peakfit(4,3)];
            x0=[x0,data_EPS_peakfit(3,2),data_EPS_peakfit(4,2)];
        end
    end
    
end


%--------------------------------------------------------------------------
% Section 6.3 - "Voigt" and "complexErrorFunction"
%--------------------------------------------------------------------------


function diff=Voigt(x0,fwhm_11,fwhm_22,lambda,Absorption,cc2,cc,doping,diam_n,Order,fwhm_gauss,center,center_22)

h=4.135667662e-15;
c=299792458;

% Unlike Lorentzian and Gaussian fits, no string is created to fit
% Voigtians, but a direct function handle that returns the
% difference between the absorption spectrum and the fit.
%
% Fit parameters of the Voigtian line-shape are:
% x(1)=height
% x(2)=shift of center
% x(3)=Broadening of HWHM
% x(4)=Ratio of Gaussian / Lorentzian HWHM

complete=2;

if isempty(fwhm_11)==1 % Fit S22
    complete=1;
    FWHM=fwhm_22;
    Cent=center_22;
    if isempty(cc)==1 % No EPS
        y=zeros(length(lambda),length(x0)/4);
        wg=zeros(length(x0)/4);
        wL=zeros(length(x0)/4);
        ll=length(x0);
    else
        ll=length(x0)-(3+2*(length(cc)-1));
        y=zeros(length(lambda),ll/4+length(cc));
        wg=zeros(ll/4);
        wL=zeros(ll/4);
    end
end

if isempty(fwhm_22)==1 % Fit S11
    complete=1;
    FWHM=fwhm_11;
    Cent=center;
    if isempty(cc)==1 % No EPS
        y=zeros(length(lambda),length(x0)/4);
        wg=zeros(length(x0)/4);
        wL=zeros(length(x0)/4);
        ll=length(x0);
    else
        ll=length(x0)-(3+2*(length(cc)-1));
        y=zeros(length(lambda),ll/4+length(cc));
        wg=zeros(ll/4);
        wL=zeros(ll/4);
    end
end

if complete==2 % Fit S11 and S22
    FWHM=fwhm_11;
    Cent=center;
    if isempty(cc)==1 && isempty(cc2)==1 % No EPS
        ll=length(x0)/2;
        y=zeros(length(lambda),length(x0)/4);
        wg=zeros(length(x0)/4);
        wL=zeros(length(x0)/4);
    else
        if isempty(cc)==1
            cc_len=0;
            cc_j=0;
        else
            cc_len=3+2*(length(cc)-1);
            cc_j=length(cc);
        end
        if isempty(cc2)==1
            cc_len2=0;
            cc_j2=0;
        else
            cc_len2=3+2*(length(cc2)-1);
            cc_j2=length(cc2);
        end
        ll=length(x0)-cc_len-cc_len2;
        y=zeros(length(lambda),ll/4+cc_j+cc_j2);
        wg=zeros(ll/4);
        wL=zeros(ll/4);
    end
end

Lamb=lambda;
FWHM_g=fwhm_gauss;

% Calculate Voigt line-shape for S11 or S22 region

for j=1:length(Cent)
    if j==1
        wL(j)=FWHM(j)*x0(4*j-1)/(0.5436+sqrt(x0(4*j)^2+.2166));
        wg(j)=wL(j)*x0(4*j);
    else
        wL(j)=FWHM(j)*x0(4*j-1)/(0.5436+sqrt((x0(4)*x0(4*j))^2+.2166));
        wg(j)=wL(j)*x0(4*j)*x0(4);
    end
    wv(j)=0.5436*wL(j)+sqrt(0.2166*wL(j)^2+wg(j)^2);
    XX=sqrt(log(2))*(Lamb-Cent(j)-x0(4*j-2))/wg(j);
    YY=sqrt(log(2))*wL(j)/wg(j);
    y(:,j)=x0(4*j-3).*real(complexErrorFunction(XX,YY))./real(complexErrorFunction(0,YY));
end

% Calculate Voigt line-shape for S22 region, after the ones for the S11
% region were calculated.

if complete==2
    n=0;
    for j=1:length(center_22)
        if center_22(j)>0
            if j==1
                wL(j+length(Cent)-n)=fwhm_22(j)*x0(4*(length(Cent)-n+j)-1)/(0.5436+sqrt(x0(4*(length(Cent)-n+j))^2+.2166));
                wg(j+length(Cent)-n)=wL(j+length(Cent)-n)*x0(4*(length(Cent)-n+j));
            else
                wL(j+length(Cent)-n)=fwhm_22(j)*x0(4*(length(Cent)-n+j)-1)/(0.5436+sqrt((x0(4*(length(Cent)-n+j))*x0(4*(length(Cent)-n+1)))^2+.2166));
                wg(j+length(Cent)-n)=wL(j+length(Cent)-n)*x0(4*(length(Cent)-n+j))*x0(4*(length(Cent)-n+1));
            end
            
            wv(j+length(Cent)-n)=0.5436*wL(j+length(Cent)-n)+sqrt(0.2166*wL(j+length(Cent)-n)^2+wg(j+length(Cent)-n)^2);
            if doping == 0
                if wv(j+length(Cent)-n)>wv(j-n)
                    wg(j+length(Cent)-n)=wg(j-n);
                    wL(j+length(Cent)-n)=wL(j-n);
                end
            end
            XX=sqrt(log(2))*(Lamb-center_22(j)-x0(4*(length(Cent)-n+j)-2))/wg(j+length(Cent)-n);
            YY=sqrt(log(2))*wL(j+length(Cent)-n)/wg(j+length(Cent)-n);
            if j==1
                y(:,j+length(Cent)-n)=x0(4*j-3)/x0(4*(j+length(Cent)-n)-3).*real(complexErrorFunction(XX,YY))./real(complexErrorFunction(0,YY));
            else
                y(:,j+length(Cent)-n)=x0(4*j-3)/(x0(4*(1+length(Cent)-n)-3)*x0(4*(j+length(Cent)-n)-3)).*real(complexErrorFunction(XX,YY))./real(complexErrorFunction(0,YY));
            end
        else
            n=n+1;
        end
    end
end

% Calculate the line-shape of the Gaussian exciton phonon sideband for the
% S11 or S22 region based on the area of the corresponding nanotube peak.

area=zeros(2*length(cc),1);

for j=1:length(cc)
    area(j)=x0(4*cc(j)-3)*wg(cc(j))*sqrt(pi/log(2))/real(complexErrorFunction(0,sqrt(log(2))*wL(cc(j))/wg(cc(j))));
    if j==1
        y(:,ll/4+j)=area(j)*(0.017+0.1/diam_n(Order(cc(j)))+x0(ll+1))*exp(-2*(Lamb-h*c/(1e-9*(h*c/((Cent(cc(j))+x0((cc(j)*4-2)))*1e-9)+0.2+x0(ll+2)))).^2/(FWHM_g*x0(ll+3)/sqrt(log(4)))^2)/(FWHM_g*x0(ll+3)*sqrt(pi/(2*log(4))));
    else
        y(:,ll/4+j)=area(j)*(0.017+0.1/diam_n(Order(cc(j)))+x0(ll+1))*exp(-2*(Lamb-h*c/(1e-9*(h*c/((Cent(cc(j))+x0((cc(j)*4-2)))*1e-9)+0.2+x0(ll+2+2*(j-1))))).^2/(FWHM_g*x0(ll+3)*x0(ll+3+2*(j-1))/sqrt(log(4)))^2)/(FWHM_g*x0(ll+3)*x0(ll+3+2*(j-1))*sqrt(pi/(2*log(4))));
    end
end

% Calculate the line-shape of the Gaussian exciton phonon sideband for the
% S22 region based on the area of the corresponding nanotube peak.

if complete==2 && isempty(cc2)==0
    for j=1:length(cc2)        
        if j==1
            area(j)=x0(4*cc2(j)-3)/x0(4*(cc2(j)+length(Cent)-n)-3)*wg(cc2(j)+length(Cent)-n)*sqrt(pi/log(2))/real(complexErrorFunction(0,sqrt(log(2))*wL(cc2(j)+length(Cent)-n)/wg(cc2(j)+length(Cent)-n)));
            y(:,2*length(Cent)-n+cc_j+j)=area(j)*(0.017+0.1/diam_n(Order(cc2(j)))+x0(ll+1+cc_len))*exp(-2*(Lamb-h*c/(1e-9*(h*c/((center_22(cc2(j))+x0((cc2(j)*4-2+length(Cent)*4)))*1e-9)+0.2+x0(ll+2+cc_len)))).^2/(FWHM_g*x0(ll+3+cc_len)/sqrt(log(4)))^2)/(FWHM_g*x0(ll+3+cc_len)*sqrt(pi/(2*log(4))));
        else
            area(j)=x0(4*cc2(j)-3)/(x0(4*(cc2(1)+length(Cent)-n)-3)*x0(4*(cc2(j)+length(Cent)-n)-3))*wg(cc2(j)+length(Cent)-n)*sqrt(pi/log(2))/real(complexErrorFunction(0,sqrt(log(2))*wL(cc2(j)+length(Cent)-n)/wg(cc2(j)+length(Cent)-n)));
            y(:,2*length(Cent)-n+cc_j+j)=area(j)*(0.017+0.1/diam_n(Order(cc2(j)))+x0(ll+1+cc_len))*exp(-2*(Lamb-h*c/(1e-9*(h*c/((center_22(cc2(j))+x0((cc2(j)*4-2+length(Cent)*4)))*1e-9)+0.2+x0(ll+2+2*(j-1)+cc_len)))).^2/(FWHM_g*x0(ll+3+cc_len)*x0(ll+3+2*(j-1)+cc_len)/sqrt(log(4)))^2)/(FWHM_g*x0(ll+3+cc_len)*x0(ll+3+2*(j-1)+cc_len)*sqrt(pi/(2*log(4))));
        end
    end
end

% Calculate the residuals by subtracting the sum of the Voigtian
% line-shapes and possible phonon sidebands from the absorption spectrum.

diff=double(Absorption-sum(y,2));


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
% Section 7 - "fit_data_m"
%--------------------------------------------------------------------------


function [x_fit,sse,cm,CNTn,dm,M11_n,FWHM_M11]=fit_data_m(height_22,center_22,FWHM_22,FWHM_g_22,f1_22,w,diam_22,M11n,fwhm_nm_m11_n,height_m,lambda,Absorption,data_S22,phonon_pos2,h,c,data_EPS,data_M11,data_voigt_s22,data_voigt,metallic_peakassignment,colors_m,CNT_n,c_len,diam_m,y,fid,protocol,fwhm_11,height_11,doping,parallel_fit)

% If in a previous fit of the entire region a Lorentzian or Gaussian
% line-shape was defined ("1" or "2") the subfunction "generate_string_m",
% located in Section 7.1 is used to create a function handle for the
% calculation of the residuals of the fit. If the user specified a Voigtian
% line-shape ("3") the subfunction "generate_startvalues_voigt_m", located
% in Section 7.2 is used to create the initial starting values for the fit
% of the absorption spectrum with Voigtians. The fit of Voigtian
% line-shapes itself is performed in Section 7.3 with the function
% "Voigt_m" and Section 6.3 with the function "complexErrorFunction".

if metallic_peakassignment==1
    [dm,cm,M11_n,FWHM_M11,CNTn,height_m]=select_nanotubes_metallic(c_len,Absorption,lambda,h,c,y,fid,protocol);
else
    cm=colors_m;
    CNTn=CNT_n;
    dm=diam_m;
    M11_n=M11n;
    FWHM_M11=fwhm_nm_m11_n;
end

switch w
    case{1,2}
        if parallel_fit=='y'
            try
                options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','Display','off','UseParallel',true);
            catch
                options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','Display','off');
            end
        else
            options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','Display','off');
        end
        
        [plotString,x0,lb,ub]=generate_string_m(height_22,center_22,FWHM_22,FWHM_g_22,f1_22,w,diam_22,lambda,Absorption,data_S22,M11_n,FWHM_M11,height_m,phonon_pos2,h,c,data_EPS,data_M11,height_11,doping);
        
        Fit_sum=str2func(plotString);
        x_fit=lsqnonlin(Fit_sum,x0,double(lb),double(ub),options);
        sse=sum(Fit_sum(x_fit).^2)/sum((Absorption-mean(Absorption)).^2);
    case 3
        if parallel_fit=='y'
            try
                options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','Display','off','FinDiffRelStep',1e-3,'UseParallel',true);
            catch
                options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','Display','off','FinDiffRelStep',1e-3);
            end
        else
            options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','Display','off','FinDiffRelStep',1e-3);
        end
        
        data_S22=[data_S22; data_voigt_s22];
        [lb,ub,x0]=generate_startvalues_voigt_m(center_22,height_22,data_S22,FWHM_22,M11_n,height_m,data_M11,data_voigt,f1_22,data_EPS,phonon_pos2,height_11,doping);
        
        x_fit=lsqnonlin(@(x) Voigt_m(x,FWHM_M11,M11_n,lambda,Absorption,phonon_pos2,FWHM_g_22,diam_22,fwhm_11,doping,center_22),double(x0),double(lb),double(ub),options);
        sse=sum(Voigt_m(x_fit,FWHM_M11,M11_n,lambda,Absorption,phonon_pos2,FWHM_g_22,diam_22,fwhm_11,doping,center_22).^2)/sum((Absorption-mean(Absorption)).^2);
end


%--------------------------------------------------------------------------
% Section 7.1 - "generate_string_m"
%--------------------------------------------------------------------------


function [plotString,x0,lb,ub]=generate_string_m(Height,center_22,FWHM_22,FWHM_g_22,f1_22,w,diam_22,lambda,Absorption,data_S22,M11n,fwhm_nm_m11_n,height_m,phonon_pos2,h,c,data_EPS,data_M11,height_11,doping)

lb=[];
ub=[];
x0=[];
plotString='';

% Creation of a function handle for the fit of S22 and M11 (as defined by the
% user). Each peak is defined by three components that can vary during
% the fit: x(1)=height, x(2)=Shift in x-direction, x(3)=broadening of FWHM.

n=0;

for j=1:length(center_22)
    if center_22(j)>0
        if j==1
            if w==1 % Lorentzian
                plotString = strcat(plotString, ['@(x) double([',num2str(Absorption.',' %f32'),'] -x(',num2str((j-n)*3-2),')./(1+(([',num2str(lambda(1:end).',' %f32'),']-',num2str(center_22(j)),'-x(',num2str((j-n)*3-1),'))/(0.5*',num2str(FWHM_22(j)),'*x(',num2str((j-n)*3),'))).^2)']);
            else    % Gaussian
                plotString = strcat(plotString, ['@(x) double([',num2str(Absorption.',' %f32'),'] -x(',num2str((j-n)*3-2),')*exp(-log(2)*(([',num2str(lambda(1:end).',' %f32'),']-',num2str(center_22(j)),'-x(',num2str((j-n)*3-1),'))/(',num2str(FWHM_22(j)),'*x(',num2str((j-n)*3),')/2)).^2)']);
            end
        else
            
            % Repeat the process outlined for the most intense peak for any
            % other peak in set of (n,m) species defined by the user.
            
            if w==1 % Lorentzian
                plotString = strcat(plotString, ['-x(',num2str((j-n)*3-2),')./(1+(([',num2str(lambda(1:end).',' %f32'),']-',num2str(center_22(j)),'-x(',num2str((j-n)*3-1),'))/(0.5*',num2str(FWHM_22(j)),'*x(',num2str((j-n)*3),'))).^2)']);
            else % Gaussian
                plotString = strcat(plotString, ['-x(',num2str((j-n)*3-2),')*exp(-log(2)*(([',num2str(lambda(1:end).',' %f32'),']-',num2str(center_22(j)),'-x(',num2str((j-n)*3-1),'))/(',num2str(FWHM_22(j)),'*x(',num2str((j-n)*3),')/2)).^2)']);
            end
        end
        
        if isempty(height_11)==1
            lb=[lb,Height(j)*data_S22(3,1),data_S22(2,1),data_S22(1,1)];
            ub=[ub,Height(j)*data_S22(3,3),data_S22(2,3),data_S22(1,3)];
            x0=[x0,Height(j)*data_S22(3,2),data_S22(2,2),data_S22(1,2)];
        else
            if doping ==0
                if height_11(j)/Height(j)>data_S22(3,3)
                    lb=[lb,Height(j)*data_S22(3,1),data_S22(2,1),data_S22(1,1)];
                    ub=[ub,Height(j)*data_S22(3,3),data_S22(2,3),data_S22(1,3)];
                    x0=[x0,Height(j)*data_S22(3,2),data_S22(2,2),data_S22(1,2)];
                else
                    lb=[lb,Height(j)*data_S22(3,1),data_S22(2,1),data_S22(1,1)];
                    ub=[ub,Height(j),data_S22(2,3),1];
                    x0=[x0,Height(j)*data_S22(3,2),data_S22(2,2),data_S22(1,2)];
                end
            else
                lb=[lb,Height(j)*data_S22(3,1),data_S22(2,1),data_S22(1,1)];
                ub=[ub,Height(j)*data_S22(3,3),data_S22(2,3),data_S22(1,3)];
                x0=[x0,Height(j)*data_S22(3,2),data_S22(2,2),data_S22(1,2)];
            end
        end
    else
        n=n+1;
    end
end

for j=1+length(center_22):length(M11n)+length(center_22)
    
    if w==1 % Lorentzian
        plotString = strcat(plotString, ['-x(',num2str((j-n)*3-2),')./(1+(([',num2str(lambda(1:end).',' %f32'),']-',num2str(M11n(j-length(center_22))),'-x(',num2str((j-n)*3-1),'))/(0.5*',num2str(fwhm_nm_m11_n(j-length(center_22))),'*x(',num2str((j-n)*3),'))).^2)']);
    else % Gaussian
        plotString = strcat(plotString, ['-x(',num2str((j-n)*3-2),')*exp(-log(2)*(([',num2str(lambda(1:end).',' %f32'),']-',num2str(M11n(j-length(center_22))),'-x(',num2str((j-n)*3-1),'))/(',num2str(fwhm_nm_m11_n(j-length(center_22))),'*x(',num2str((j-n)*3),')/2)).^2)']);
    end
    
    lb=[lb,height_m(j-length(center_22))*data_M11(3,1),data_M11(2,1),data_M11(1,1)];
    ub=[ub,height_m(j-length(center_22))*data_M11(3,3),data_M11(2,3),data_M11(1,3)];
    x0=[x0,height_m(j-length(center_22))*data_M11(3,2),data_M11(2,2),data_M11(1,2)];
end

x0_len=length(x0);

for j=1:length(phonon_pos2)
    
    % Fit parameters:
    % x(1)=shift of diam vs. area ratio curve from Perebeinos et al.
    % PRL 94, 027402 (2005)
    % x(2)=shift of center around 0.2 eV above associated S11 or S22 peak
    % x(3)=Broadening of FWHM
    
    if j==1
        
        if w==1 % Lorentzian
            plotString = strcat(plotString,['-x(',num2str(phonon_pos2(j)*3-2),')*pi*',num2str(FWHM_22(phonon_pos2(j))),'*x(',num2str(phonon_pos2(j)*3),')/2*(0.017+0.1/',num2str(diam_22(phonon_pos2(j))),' + x(',num2str(x0_len+j),'))*exp(-2*([',num2str(lambda(1:end).',' %f32'),']-',num2str(h*c),'/(1e-9*(',num2str(h*c),'/((',num2str(center_22(phonon_pos2(j))),'+x(',num2str(phonon_pos2(j)*3-1),'))*1e-9)+0.2+x(',num2str(x0_len+2),')))).^2/(',num2str(FWHM_g_22(j)),'*x(',num2str(x0_len+3),')/sqrt(log(4)))^2)/((',num2str(FWHM_g_22(j)),'*x(',num2str(x0_len+3),'))*sqrt(pi/(2*log(4))))']);
        else % Gaussian
            plotString = strcat(plotString,['-x(',num2str(phonon_pos2(j)*3-2),')*',num2str(FWHM_22(phonon_pos2(j))),'*x(',num2str(phonon_pos2(j)*3),')*sqrt(pi/(2*log(4)))*(0.017+0.1/',num2str(diam_22(phonon_pos2(j))),' + x(',num2str(x0_len+j),'))*exp(-2*([',num2str(lambda(1:end).',' %f32'),']-',num2str(h*c),'/(1e-9*(',num2str(h*c),'/((',num2str(center_22(phonon_pos2(j))),'+x(',num2str(phonon_pos2(j)*3-1),'))*1e-9)+0.2+x(',num2str(x0_len+2),')))).^2/(',num2str(FWHM_g_22(j)),'*x(',num2str(x0_len+3),')/sqrt(log(4)))^2)/((',num2str(FWHM_g_22(j)),'*x(',num2str(x0_len+3),'))*sqrt(pi/(2*log(4))))']);
        end
        
        if f1_22(1)<0
            ub_f1=f1_22(1)*data_EPS(3,1);
            lb_f1=f1_22(1)*data_EPS(3,3);
        else
            lb_f1=f1_22(1)*data_EPS(3,1);
            ub_f1=f1_22(1)*data_EPS(3,3);
        end
        
        lb=[lb,lb_f1,data_EPS(2,1),data_EPS(1,1)];
        ub=[ub,ub_f1,data_EPS(2,3),data_EPS(1,3)];
        x0=[x0,f1_22(1)*data_EPS(3,2),data_EPS(2,2),data_EPS(1,2)];
        
    else
        
        % Fit parameters:
        % x(2+2*(j-1)) = shift of center around 0.2 eV above
        % associated S11
        % x(3+2*(j-1)) = Percentage of FWHM broadening of first
        % peak - can vary between 90 and 110 %.
        
        if w==1 % Lorentzian
            plotString = strcat(plotString,['-x(',num2str(phonon_pos2(j)*3-2),')*pi*',num2str(FWHM_22(phonon_pos2(j))),'*x(',num2str(phonon_pos2(j)*3),')/2*(0.017+0.1/',num2str(diam_22(phonon_pos2(j))),' + x(',num2str(x0_len+1),'))*exp(-2*([',num2str(lambda(1:end).',' %f32'),']-',num2str(h*c),'/(1e-9*(',num2str(h*c),'/((',num2str(center_22(phonon_pos2(j))),'+x(',num2str(phonon_pos2(j)*3-1),'))*1e-9)+0.2+x(',num2str(x0_len+4+2*(j-2)),')))).^2/(',num2str(FWHM_g_22(j)),'*x(',num2str(x0_len+3),')*x(',num2str(x0_len+5+2*(j-2)),')/sqrt(log(4)))^2)/((',num2str(FWHM_g_22(j)),'*x(',num2str(x0_len+3),')*x(',num2str(x0_len+5+2*(j-2)),'))*sqrt(pi/(2*log(4))))']);
        else % Gaussian
            plotString = strcat(plotString,['-x(',num2str(phonon_pos2(j)*3-2),')*',num2str(FWHM_22(phonon_pos2(j))),'*x(',num2str(phonon_pos2(j)*3),')*sqrt(pi/(2*log(4)))*(0.017+0.1/',num2str(diam_22(phonon_pos2(j))),' + x(',num2str(x0_len+1),'))*exp(-2*([',num2str(lambda(1:end).',' %f32'),']-',num2str(h*c),'/(1e-9*(',num2str(h*c),'/((',num2str(center_22(phonon_pos2(j))),'+x(',num2str(phonon_pos2(j)*3-1),'))*1e-9)+0.2+x(',num2str(x0_len+4+2*(j-2)),')))).^2/(',num2str(FWHM_g_22(j)),'*x(',num2str(x0_len+3),')*x(',num2str(x0_len+5+2*(j-2)),')/sqrt(log(4)))^2)/((',num2str(FWHM_g_22(j)),'*x(',num2str(x0_len+3),')*x(',num2str(x0_len+5+2*(j-2)),'))*sqrt(pi/(2*log(4))))']);
        end
        
        lb=[lb,data_EPS(2,1),data_EPS(1,1)];
        ub=[ub,data_EPS(2,3),data_EPS(1,3)];
        x0=[x0,data_EPS(2,2),data_EPS(1,2)];
        
    end
end

plotString=strcat(plotString,[')']);


%--------------------------------------------------------------------------
% Section 7.2 - "generate_startvalues_voigt_m"
%--------------------------------------------------------------------------


function [lb,ub,x0]=generate_startvalues_voigt_m(center_22,Height,data_S22,FWHM_22,M11n,height_m,data_M11,data_voigt,f1_22,data_EPS,phonon_pos2,height_11,doping)

lb=[];
ub=[];
x0=[];

% Creation of a function handle for the fit of S22 and M11 (as defined by
% the user). Each peak is defined by three components that can vary during
% the fit: x(1)=height, x(2)=Shift in x-direction, x(3)=broadening of FWHM
% and x(4)=Ratio of Gaussian / Lorentzian HWHM.

for j=1:length(center_22)
    if center_22(j)>0
        if isempty(height_11)==1
            lb=[lb,Height(j)*data_S22(3,1),center_22(j)+data_S22(2,1),0.5*FWHM_22(j,3)*data_S22(1,1),FWHM_22(j,2)/FWHM_22(j,1)*data_S22(4,1)];
            ub=[ub,Height(j)*data_S22(3,3),center_22(j)+data_S22(2,3),0.5*FWHM_22(j,3)*data_S22(1,3),FWHM_22(j,2)/FWHM_22(j,1)*data_S22(4,3)];
            x0=[x0,Height(j)*data_S22(3,2),center_22(j)+data_S22(2,2),0.5*FWHM_22(j,3)*data_S22(1,2),FWHM_22(j,2)/FWHM_22(j,1)*data_S22(4,2)];
        else
            if doping==0
                if height_11(j)/Height(j)>data_S22(3,3)
                    lb=[lb,Height(j)*data_S22(3,1),center_22(j)+data_S22(2,1),0.5*FWHM_22(j,3)*data_S22(1,1),FWHM_22(j,2)/FWHM_22(j,1)*data_S22(4,1)];
                    ub=[ub,Height(j)*data_S22(3,3),center_22(j)+data_S22(2,3),0.5*FWHM_22(j,3)*data_S22(1,3),FWHM_22(j,2)/FWHM_22(j,1)*data_S22(4,3)];
                    x0=[x0,Height(j)*data_S22(3,2),center_22(j)+data_S22(2,2),0.5*FWHM_22(j,3)*data_S22(1,2),FWHM_22(j,2)/FWHM_22(j,1)*data_S22(4,2)];
                else
                    lb=[lb,Height(j)*data_S22(3,1),center_22(j)+data_S22(2,1),0.5*FWHM_22(j,3)*data_S22(1,1),FWHM_22(j,2)/FWHM_22(j,1)*data_S22(4,1)];
                    ub=[ub,Height(j),center_22(j)+data_S22(2,3),0.5*FWHM_22(j,3)*data_S22(1,3),FWHM_22(j,2)/FWHM_22(j,1)*data_S22(4,3)];
                    x0=[x0,Height(j)*data_S22(3,2),center_22(j)+data_S22(2,2),0.5*FWHM_22(j,3)*data_S22(1,2),FWHM_22(j,2)/FWHM_22(j,1)*data_S22(4,2)];
                end
            else
                lb=[lb,Height(j)*data_S22(3,1),center_22(j)+data_S22(2,1),0.5*FWHM_22(j,3)*data_S22(1,1),FWHM_22(j,2)/FWHM_22(j,1)*data_S22(4,1)];
                ub=[ub,Height(j)*data_S22(3,3),center_22(j)+data_S22(2,3),0.5*FWHM_22(j,3)*data_S22(1,3),FWHM_22(j,2)/FWHM_22(j,1)*data_S22(4,3)];
                x0=[x0,Height(j)*data_S22(3,2),center_22(j)+data_S22(2,2),0.5*FWHM_22(j,3)*data_S22(1,2),FWHM_22(j,2)/FWHM_22(j,1)*data_S22(4,2)];
            end
        end
    end
end

count=1;
for j=1+length(center_22):length(M11n)+length(center_22)
    if count==1
        lb=[lb,height_m(j-length(center_22))*data_M11(3,1),data_M11(2,1),data_M11(1,1),data_voigt(3)];
        ub=[ub,height_m(j-length(center_22))*data_M11(3,3),data_M11(2,3),data_M11(1,3),data_voigt(4)];
        x0=[x0,height_m(j-length(center_22))*data_M11(3,2),data_M11(2,2),data_M11(1,2),data_voigt(2)];
        count=2;
    else
        lb=[lb,height_m(j-length(center_22))*data_M11(3,1),data_M11(2,1),data_M11(1,1),1-data_voigt(5)];
        ub=[ub,height_m(j-length(center_22))*data_M11(3,3),data_M11(2,3),data_M11(1,3),1+data_voigt(5)];
        x0=[x0,height_m(j-length(center_22))*data_M11(3,2),data_M11(2,2),data_M11(1,2),1];
    end
end

x0_len=length(x0);

for j=1:length(phonon_pos2)
    
    % Fit parameters:
    % x(1)=shift of diam vs. area ratio curve from Perebeinos et al.
    % PRL 94, 027402 (2005)
    % x(2)=shift of center around 0.2 eV above associated S11 or S22 peak
    % x(3)=Broadening of FWHM
    
    if j==1
        
        if f1_22(1)<0
            ub_f1=f1_22(1)*data_EPS(3,1);
            lb_f1=f1_22(1)*data_EPS(3,3);
        else
            lb_f1=f1_22(1)*data_EPS(3,1);
            ub_f1=f1_22(1)*data_EPS(3,3);
        end
        
        lb=[lb,lb_f1,data_EPS(2,1),data_EPS(1,1)];
        ub=[ub,ub_f1,data_EPS(2,3),data_EPS(1,3)];
        x0=[x0,f1_22(1)*data_EPS(3,2),data_EPS(2,2),data_EPS(1,2)];
        
    else
        
        % Fit parameters:
        % x(2+2*(j-1)) = shift of center around 0.2 eV above
        % associated S11
        % x(3+2*(j-1)) = Percentage of FWHM broadening of first
        % peak - can vary between 90 and 110 %.
        
        lb=[lb,data_EPS(2,1),data_EPS(1,1)];
        ub=[ub,data_EPS(2,3),data_EPS(1,3)];
        x0=[x0,data_EPS(2,2),data_EPS(1,2)];
        
    end
end


%--------------------------------------------------------------------------
% Section 7.3 - "Voigt_m"
%--------------------------------------------------------------------------


function diff=Voigt_m(x0,fwhm_nm_m11_n,M11n,lambda,Absorption,phonon_pos2,FWHM_g_22,diam_22,fwhm_11,doping,center_22)

h=4.135667662e-15;
c=299792458;

% Unlike Lorentzian and Gaussian fits, no string is created for the fit of
% Voigtians, but a direct function handle that returns the
% difference between the absorption spectrum and the fit.
%
% Fit parameters of the Voigtian line-shape are:
% x(1)=height
% x(2)=shift of center
% x(3)=Lorentzian HWHM
% x(4)=Ratio of Gaussian / Lorentzian HWHM

Lamb=lambda;

% Calculate Voigt line-shape for S22 region

n=0;
for j=1:length(center_22)
    if center_22(j)>0
        
        wv(j-n)=x0(4*(j-n)-1);
        wL(j-n)=wv(j-n)/(0.5436+sqrt(x0(4*(j-n))^2+.2166));
        wg(j-n)=wL(j-n)*x0(4*(j-n));
        
        if doping == 0
            if isempty(fwhm_11)==0
                if wv(j-n)>fwhm_11(j-n,3)
                    wL(j-n)=fwhm_11(j-n,1);
                    wg(j-n)=fwhm_11(j-n,2);
                end
            end
        end
        XX=sqrt(log(2))*(Lamb-x0(4*(j-n)-2))/wg(j-n);
        YY=sqrt(log(2))*wL(j-n)/wg(j-n);
        y(:,j)=x0(4*(j-n)-3).*real(complexErrorFunction(XX,YY))./real(complexErrorFunction(0,YY));
    else
        n=n+1;
    end
end

% Calculate Voigt line-shape for M11 nanotubes

for j=1:length(M11n)
    if j==1
        wL(j+length(center_22)-n)=0.5*fwhm_nm_m11_n(j)*x0(4*(length(center_22)+j-n)-1)/(0.5436+sqrt(x0(4*(length(center_22)+j-n))^2+.2166));
        wg(j+length(center_22)-n)=wL(j+length(center_22)-n)*x0(4*(length(center_22)+j-n));
    else
        wL(j+length(center_22)-n)=0.5*fwhm_nm_m11_n(j)*x0(4*(length(center_22)+j-n)-1)/(0.5436+sqrt((x0(4*(length(center_22)+j-n))*x0(4*(length(center_22)+1-n)))^2+.2166));
        wg(j+length(center_22)-n)=wL(j+length(center_22)-n)*x0(4*(length(center_22)+j-n))*x0(4*(length(center_22)+1-n));
    end
    
    XX=sqrt(log(2))*(Lamb-M11n(j)-x0(4*(length(center_22)+j-n)-2))/wg(j+length(center_22)-n);
    YY=sqrt(log(2))*wL(j+length(center_22)-n)/wg(j+length(center_22)-n);
    y(:,j+length(center_22)-n)=x0(4*(j+length(center_22)-n)-3).*real(complexErrorFunction(XX,YY))./real(complexErrorFunction(0,YY));
end

% Calculate the line-shape of the Gaussian exciton phonon sideband for the
% S11 or S22 region based on the area of the corresponding nanotube peak.

area=zeros(2*length(phonon_pos2),1);

ll=(length(center_22)+length(M11n)-n)*4;

for j=1:length(phonon_pos2)
    area(j)=trapz(lambda,y(:,phonon_pos2(j)));
    if j==1
        y(:,ll/4+j)=area(j)*(0.017+0.1/diam_22(phonon_pos2(j))+x0(ll+1))*exp(-2*(Lamb-h*c/(1e-9*(h*c/((x0((phonon_pos2(j)*4-2)))*1e-9)+0.2+x0(ll+2)))).^2/(FWHM_g_22(j)*x0(ll+3)/sqrt(log(4)))^2)/(FWHM_g_22(j)*x0(ll+3)*sqrt(pi/(2*log(4))));
    else
        y(:,ll/4+j)=area(j)*(0.017+0.1/diam_22(phonon_pos2(j))+x0(ll+1))*exp(-2*(Lamb-h*c/(1e-9*(h*c/((x0((phonon_pos2(j)*4-2)))*1e-9)+0.2+x0(ll+2+2*(j-1))))).^2/(FWHM_g_22(j)*x0(ll+3)*x0(ll+3+2*(j-1))/sqrt(log(4)))^2)/(FWHM_g_22(j)*x0(ll+3)*x0(ll+3+2*(j-1))*sqrt(pi/(2*log(4))));
    end
end

% Calculate the residuals by subtracting the sum of the Voigtian
% line-shapes and possible phonon sidebands from the absorption spectrum.

diff=double(Absorption-sum(y,2));


%..........................................................................


%--------------------------------------------------------------------------
% Section 8 - "plot_data"
%--------------------------------------------------------------------------


function [legstr,plot_legstr,y]=plot_data(x_fit,colors_n,lambda,Absorption,Order,fwhm_nm_n,fwhm_gauss,diam_n,Names,vec,center,sse,cc,fwhm_nm_n22,cc2,fit_region,w,fwhm_nm_v22,fwhm_nm_v,complete,start,End)
global center_22

% Based on the line-shape chosen by the user and the results of the fit
% performed in Section 6 (stored in "x_fit"), the individual (n,m) species
% are plotted along with the calculated spectrum (sum of individual (n,m)
% fits) and the absorption spectrum initially loaded by the user. If
% Lorentzian or Gaussian line-shapes were defined, they are plotted with
% "plot_fit" defined in Section 8.1. If Voigtian line-shapes were chosen,
% they are plotted with "plot_fit_voigt" defined in Section 8.2.

switch w
    case {1,2}
        switch fit_region
            case 1 % Fit S11 region
                if isempty(complete)==1 || complete == 'n'
                    [legstr,plot_legstr,y]=plot_fit(x_fit,colors_n,lambda(start:End),Absorption(start:End),Order,fwhm_nm_n,fwhm_gauss,diam_n,Names,vec,center,w,sse,cc,[],[],[]);
                else % Fit entire region
                    [legstr,plot_legstr,y]=plot_fit(x_fit,colors_n,lambda,Absorption,Order,fwhm_nm_n,fwhm_gauss,diam_n,Names,vec,center,w,sse,cc,center_22,fwhm_nm_n22,cc2);
                end
            case 2 % Fit S22 region
                [legstr,plot_legstr,y]=plot_fit(x_fit,colors_n,lambda(start:End),Absorption(start:End),Order,[],fwhm_gauss,diam_n,Names,vec,[],w,sse,cc,center_22,fwhm_nm_n22,[]);
        end
    case 3
        switch fit_region
            case 1 % Fit S11 region
                if isempty(complete)==1 || complete == 'n'
                    [legstr,plot_legstr,y]=plot_fit_voigt(x_fit,colors_n,lambda(start:End),Absorption(start:End),Order,fwhm_gauss,diam_n,Names,vec,center,sse,cc,[],[],[],fwhm_nm_v);
                else % Fit entire region
                    [legstr,plot_legstr,y]=plot_fit_voigt(x_fit,colors_n,lambda,Absorption,Order,fwhm_gauss,diam_n,Names,vec,center,sse,cc,cc2,center_22,fwhm_nm_v22,fwhm_nm_v);
                end
            case 2 % Fit S22 region
                [legstr,plot_legstr,y]=plot_fit_voigt(x_fit,colors_n,lambda(start:End),Absorption(start:End),Order,fwhm_gauss,diam_n,Names,vec,[],sse,cc,cc2,center_22,fwhm_nm_v22,[]);
        end
end


%--------------------------------------------------------------------------
% Section 8.1 - "plot_fit"
%--------------------------------------------------------------------------


function [legstr,plot_legstr,y]=plot_fit(x_fit,colors_n,lambda,Absorption,Order,fwhm_nm_n,fwhm_gauss,diam_n,Names,vec,center,w,sse,cc,center_22,fwhm_nm_n22,cc2)
global h c

% This subfunction is used to plot the result of a Gaussian or Lorentzian
% fit. The structure of the code follows the rountine outlined in Section
% 6.1 "generate_string".

figure('units','normalized','outerposition',[0 0 1 1]), hold on

complete=0; % Plot S11 and S22 region

if isempty(center_22)==1 % Plot only S11
    complete=1;
    Cent=center;
    FWHM=fwhm_nm_n;
    if isempty(cc)==1
        y=zeros(length(lambda),length(center));
    else
        y=zeros(length(lambda),length(Cent)+length(cc));
    end
    x0_len=3*length(Cent);
    region='S_1_1-';
end

if isempty(center)==1 % Plot only S22
    complete=1;
    Cent=center_22;
    FWHM=fwhm_nm_n22;
    if isempty(cc)==1
        y=zeros(length(lambda),length(center));
    else
        y=zeros(length(lambda),length(Cent)+length(cc));
    end
    x0_len=3*length(Cent);
    region='S_2_2-';
end

if complete==0 % Plot S11 and S22
    Cent=center;
    FWHM=fwhm_nm_n;
    if isempty(cc)==1
        y=zeros(length(lambda),2*length(center)-sum(center_22==0));
    else
        y=zeros(length(lambda),2*length(center)-sum(center_22==0)+length(cc)+length(cc2));
    end
    x0_len=6*length(Cent)-3*sum(center_22==0);
    region='S_1_1-';
end

% An individual color is assigned to each peak in S11 or S22.

for j=1:length(colors_n)
    
    if w==1 % Plot Lorentzian
        y(:,j)=x_fit(j*3-2) ./ (1 + ((lambda-Cent(j)-x_fit(j*3-1))/(0.5*FWHM(Order(j))*x_fit(j*3))).^2);
        title_name='Lorentzian Fit with ';
    else    % Plot Gaussian
        y(:,j)=x_fit(j*3-2)*exp(-log(2)*((lambda-Cent(j)-x_fit(j*3-1))/(FWHM(Order(j))/2*x_fit(j*3))).^2);
        title_name='Gaussian Fit with ';
    end
    
    plot(lambda,y(:,j),'color',colors_n(Order(j),:),'linewidth',1);
    
    legstr{1,j}=[Names{vec(Order(j)),1}];
    plot_legstr{1,j}=[region,Names{vec(Order(j)),1}];
end

% The same color used for the (n,m) species in region S11 is used in the
% S22 region.

if complete == 0
    n=0;
    for j=1:length(center_22)
        if center_22(j) > 0
            if j==1
                if w==1 % Plot Lorentzian
                    y(:,length(center)+j-n)=x_fit(j*3-2-3*n)/(x_fit(length(center)*3+j*3-2-3*n)) ./ (1 + ((lambda-center_22(j)-x_fit(length(center)*3+j*3-1-3*n))/(0.5*fwhm_nm_n22(Order(j))*x_fit(length(center)*3+j*3-3*n))).^2);
                else    % Plot Gaussian
                    y(:,length(center)+j-n)=x_fit(j*3-2-3*n)/(x_fit(length(center)*3+j*3-2-3*n))*exp(-log(2)*((lambda-center_22(j)-x_fit(length(center)*3+j*3-1-3*n))/(fwhm_nm_n22(Order(j))/2*x_fit(length(center)*3+j*3-3*n))).^2);
                end
            else
                if w==1 % Plot Lorentzian
                    y(:,length(center)+j-n)=x_fit(j*3-2-3*n)/(x_fit(length(center)*3+1-3*n)*x_fit(length(center)*3+j*3-2-3*n)) ./ (1 + ((lambda-center_22(j)-x_fit(length(center)*3+j*3-1-3*n))/(0.5*fwhm_nm_n22(Order(j))*x_fit(length(center)*3+j*3-3*n))).^2);
                else    % Plot Gaussian
                    y(:,length(center)+j-n)=x_fit(j*3-2-3*n)/(x_fit(length(center)*3+1-3*n)*x_fit(length(center)*3+j*3-2-3*n))*exp(-log(2)*((lambda-center_22(j)-x_fit(length(center)*3+j*3-1-3*n))/(fwhm_nm_n22(Order(j))/2*x_fit(length(center)*3+j*3-3*n))).^2);
                end
            end
            plot(lambda,y(:,length(center)+j-n),'color',colors_n(Order(j),:),'linewidth',1);
            
            legstr{1,length(legstr)+1}=[Names{vec(Order(j)),1}];
            plot_legstr{1,length(plot_legstr)+1}=['S_2_2-',Names{vec(Order(j)),1}];
        else
            n=n+1;
        end
    end
end

% A different color is assigned to exciton phonon sidebands (EPS) compared
% to the color assigned to an (n,m) species.

for j=1:length(cc)
    if j==1
        if w==1 % Lorentzian
            y(:,x0_len/3+j)=x_fit(cc(j)*3-2)*pi*FWHM(Order(cc(j)))*x_fit(cc(j)*3)/2*(0.017+0.1/diam_n(Order(cc(j)))+x_fit(x0_len+1))*exp(-2*(lambda-h*c/(1e-9*(h*c/((Cent(cc(j))+x_fit((cc(j)*3-1)))*1e-9)+0.2+x_fit(x0_len+2)))).^2/(fwhm_gauss*x_fit(x0_len+3)/sqrt(log(4)))^2)/(fwhm_gauss*x_fit(x0_len+3)*sqrt(pi/(2*log(4))));
        else % Gaussian
            y(:,x0_len/3+j)=x_fit(cc(j)*3-2)*FWHM(Order(cc(j)))*x_fit(cc(j)*3)*sqrt(pi/(2*log(4)))*(0.017+0.1/diam_n(Order(cc(j)))+x_fit(x0_len+1))*exp(-2*(lambda-h*c/(1e-9*(h*c/((Cent(cc(j))+x_fit((cc(j)*3-1)))*1e-9)+0.2+x_fit(x0_len+2)))).^2/(fwhm_gauss*x_fit(x0_len+3)/sqrt(log(4)))^2)/(fwhm_gauss*x_fit(x0_len+3)*sqrt(pi/(2*log(4))));
        end
        plot(lambda,y(:,x0_len/3+j),'color',[colors_n(Order(cc(j)),1),.8,colors_n(Order(cc(j)),3)],'linewidth',1);
        
        legstr{1,length(legstr)+1}=[Names{vec(Order(j)),1}];
        plot_legstr{1,length(plot_legstr)+1}=['EPS-',region,Names{vec(Order(cc(j))),1}];
    else
        if w==1
            y(:,x0_len/3+j)=x_fit(cc(j)*3-2)*pi*FWHM(Order(cc(j)))*x_fit(cc(j)*3)/2*(0.017+0.1/diam_n(Order(cc(j)))+x_fit(x0_len+1))*exp(-2*(lambda-h*c/(1e-9*(h*c/((Cent(cc(j))+x_fit((cc(j)*3-1)))*1e-9)+0.2+x_fit(x0_len+4+2*(j-2))))).^2/(fwhm_gauss*x_fit(x0_len+3)*x_fit(x0_len+5+2*(j-2))/sqrt(log(4)))^2)/(fwhm_gauss*x_fit(x0_len+3)*x_fit(x0_len+5+2*(j-2))*sqrt(pi/(2*log(4))));
        else
            y(:,x0_len/3+j)=x_fit(cc(j)*3-2)*FWHM(Order(cc(j)))*x_fit(cc(j)*3)*sqrt(pi/(2*log(4)))*(0.017+0.1/diam_n(Order(cc(j)))+x_fit(x0_len+1))*exp(-2*(lambda-h*c/(1e-9*(h*c/((Cent(cc(j))+x_fit((cc(j)*3-1)))*1e-9)+0.2+x_fit(x0_len+4+2*(j-2))))).^2/(fwhm_gauss*x_fit(x0_len+3)*x_fit(x0_len+5+2*(j-2))/sqrt(log(4)))^2)/(fwhm_gauss*x_fit(x0_len+3)*x_fit(x0_len+5+2*(j-2))*sqrt(pi/(2*log(4))));
        end
        plot(lambda,y(:,x0_len/3+j),'color',[colors_n(Order(cc(j)),1),.8,colors_n(Order(cc(j)),3)],'linewidth',1);
        
        legstr{1,length(legstr)+1}=[Names{vec(Order(cc(j))),1}];
        plot_legstr{1,length(plot_legstr)+1}=['EPS-',region,Names{vec(Order(cc(j))),1}];
    end
end

% The same color used for the EPS in region S11 is taken for the EPS in
% S22.

if complete==0 && isempty(cc2)==0
    
    for j=1:length(cc2)
        if isempty(cc)
            cc_len=0;
            cc_j=0;
        else
            cc_len=3+2*(length(cc)-1);
            cc_j=length(cc);
        end
        
        if j==1
            
            if w==1 % Lorentzian
                y(:,x0_len/3+j+cc_j)=x_fit(cc2(j)*3-2)/(x_fit(length(center)*3+cc2(j)*3-2))*pi*fwhm_nm_n22(Order(cc2(j)))*x_fit(length(center)*3+cc2(j)*3)/2*(0.017+0.1/diam_n(Order(cc2(j)))+x_fit(x0_len+cc_len+j))*exp(-2*(lambda-h*c/(1e-9*(h*c/((center_22(cc2(j))+x_fit((length(center)*3+cc2(j)*3-1)))*1e-9)+0.2+x_fit(x0_len+cc_len+1+j)))).^2/(fwhm_gauss*x_fit(x0_len+2+cc_len+j)/sqrt(log(4)))^2)/(fwhm_gauss*x_fit(x0_len+2+cc_len+j)*sqrt(pi/(2*log(4))));
            else % Gaussian
                y(:,x0_len/3+j+cc_j)=x_fit(cc2(j)*3-2)/(x_fit(length(center)*3+cc2(j)*3-2))*fwhm_nm_n22(Order(cc2(j)))*x_fit(length(center)*3+cc2(j)*3)*sqrt(pi/(2*log(4)))*(0.017+0.1/diam_n(Order(cc2(j)))+x_fit(x0_len+cc_len+j))*exp(-2*(lambda-h*c/(1e-9*(h*c/((center_22(cc2(j))+x_fit((length(center)*3+cc2(j)*3-1)))*1e-9)+0.2+x_fit(x0_len+cc_len+1+j)))).^2/(fwhm_gauss*x_fit(x0_len+2+cc_len+j)/sqrt(log(4)))^2)/(fwhm_gauss*x_fit(x0_len+2+cc_len+j)*sqrt(pi/(2*log(4))));
            end
            plot(lambda,y(:,x0_len/3+j+cc_j),'color',[colors_n(Order(cc2(j)),1),.8,colors_n(Order(cc2(j)),3)],'linewidth',1);
            
            legstr{1,length(legstr)+1}=[Names{vec(Order(cc2(j))),1}];
            plot_legstr{1,length(plot_legstr)+1}=['EPS-S_2_2-',Names{vec(Order(cc2(j))),1}];
        else
            if w==1
                y(:,x0_len/3+j+cc_j)=x_fit(cc2(j)*3-2)/(x_fit(length(center)*3+cc2(1)*3-2)*x_fit(length(center)*3+cc2(j)*3-2))*pi*fwhm_nm_n22(Order(cc2(j)))*x_fit(length(center)*3+cc2(j)*3)/2*(0.017+0.1/diam_n(Order(cc2(j)))+x_fit(x0_len+1+cc_len))*exp(-2*(lambda-h*c/(1e-9*(h*c/((center_22(cc2(j))+x_fit(length(center)*3+cc2(j)*3-1))*1e-9)+0.2+x_fit(x0_len+4+cc_len+2*(j-2))))).^2/(fwhm_gauss*x_fit(x0_len+3+cc_len)*x_fit(x0_len+5+cc_len+2*(j-2))/sqrt(log(4)))^2)/(fwhm_gauss*x_fit(x0_len+3+cc_len)*x_fit(x0_len+5+cc_len+2*(j-2))*sqrt(pi/(2*log(4))));
            else
                y(:,x0_len/3+j+cc_j)=x_fit(cc2(j)*3-2)/(x_fit(length(center)*3+cc2(1)*3-2)*x_fit(length(center)*3+cc2(j)*3-2))*fwhm_nm_n22(Order(cc2(j)))*x_fit(length(center)*3+cc2(j)*3)*sqrt(pi/(2*log(4)))*(0.017+0.1/diam_n(Order(cc2(j)))+x_fit(x0_len+1+cc_len))*exp(-2*(lambda-h*c/(1e-9*(h*c/((center_22(cc2(j))+x_fit(length(center)*3+cc2(j)*3-1))*1e-9)+0.2+x_fit(x0_len+4+cc_len+2*(j-2))))).^2/(fwhm_gauss*x_fit(x0_len+3+cc_len)*x_fit(x0_len+5+cc_len+2*(j-2))/sqrt(log(4)))^2)/(fwhm_gauss*x_fit(x0_len+3+cc_len)*x_fit(x0_len+5+cc_len+2*(j-2))*sqrt(pi/(2*log(4))));
            end
            plot(lambda,y(:,x0_len/3+j+cc_j),'color',[colors_n(Order(cc2(j)),1),.8,colors_n(Order(cc2(j)),3)],'linewidth',1);
            
            legstr{1,length(legstr)+1}=[Names{vec(Order(cc2(j))),1}];
            plot_legstr{1,length(plot_legstr)+1}=['EPS-S_2_2-',Names{vec(Order(cc2(j))),1}];
        end
    end
    
end

% Plot the absorption spectrum supplied by the user, the calculated
% absorption spectrum based on the sum of all fits and the chi-square
% value.

plot(lambda,Absorption,'k','linewidth',2)

plot(lambda,sum(y,2),'r','linewidth',2)

legend([plot_legstr,{'Measured Spectrum'},{'Calculated Spectrum'}])

title([title_name,': nSSE = ',num2str(sse)],'FontSize',26)

xlabel('Wavelength (nm)','FontSize',24)
ylabel('Absorption (a.u.)','FontSize',24)
set(gca,'LineWidth',2)
set(gca,'FontSize',22)


%--------------------------------------------------------------------------
% Section 8.2 - "plot_fit_voigt"
%--------------------------------------------------------------------------


function [legstr,plot_legstr,y]=plot_fit_voigt(x0,colors_n,lambda,Absorption,Order,fwhm_gauss,diam_n,Names,vec,center,sse,cc,cc2,center_22,fwhm_22,fwhm_11)
global h c

% This subfunction is used to plot the result of a Voigtian fit. The
% structure of the code follows the rountine outlined in Section
% 6.3 "Voigt".

figure('units','normalized','outerposition',[0 0 1 1]), hold on

complete=2;

if isempty(fwhm_11)==1 % Plot S22
    complete=1;
    FWHM=fwhm_22;
    Cent=center_22;
    if isempty(cc)==1 % No EPS
        y=zeros(length(lambda),length(x0)/4);
        wg=zeros(length(x0)/4);
        ll=length(x0);
    else
        ll=length(x0)-(3+2*(length(cc)-1));
        y=zeros(length(lambda),ll/4+length(cc));
        wg=zeros(ll/4);
    end
    region='S_2_2-';
end

if isempty(fwhm_22)==1 % Plot S11
    complete=1;
    FWHM=fwhm_11;
    Cent=center;
    if isempty(cc)==1 % No EPS
        y=zeros(length(lambda),length(x0)/4);
        wg=zeros(length(x0)/4);
        ll=length(x0);
    else
        ll=length(x0)-(3+2*(length(cc)-1));
        y=zeros(length(lambda),ll/4+length(cc));
        wg=zeros(ll/4);
    end
    region='S_1_1-';
end

if complete==2 % Plot S11 and S22
    FWHM=fwhm_11;
    Cent=center;
    if isempty(cc)==1 && isempty(cc2)==1 % No EPS
        ll=length(x0)/2;
        y=zeros(length(lambda),length(x0)/4);
        wg=zeros(length(x0)/4);
    else
        if isempty(cc)==1
            cc_len=0;
            cc_j=0;
        else
            cc_len=3+2*(length(cc)-1);
            cc_j=length(cc);
        end
        if isempty(cc2)==1
            cc_len2=0;
            cc_j2=0;
        else
            cc_len2=3+2*(length(cc2)-1);
            cc_j2=length(cc2);
        end
        ll=length(x0)-cc_len-cc_len2;
        y=zeros(length(lambda),ll/4+cc_j+cc_j2);
        wg=zeros(ll/4);
    end
    region='S_1_1-';
end

Lamb=lambda;
FWHM_g=fwhm_gauss;

% Plot S11 or S22 region and assign a different color to each (n,m) species

for j=1:length(Cent)
    if j==1
        wL(j)=FWHM(j)*x0(4*j-1)/(0.5436+sqrt(x0(4*j)^2+.2166));
        wg(j)=wL(j)*x0(4*j);
    else
        wL(j)=FWHM(j)*x0(4*j-1)/(0.5436+sqrt((x0(4*j)*x0(4*1))^2+.2166));
        wg(j)=wL(j)*x0(4*j)*x0(4);
    end
    
    XX=sqrt(log(2))*(Lamb-Cent(j)-x0(4*j-2))/wg(j);
    YY=sqrt(log(2))*wL(j)/wg(j);
    y(:,j)=x0(4*j-3).*real(complexErrorFunction(XX,YY))./real(complexErrorFunction(0,YY));
    
    title_name='Voigtian Fit with ';
    
    plot(lambda,y(:,j),'color',colors_n(Order(j),:),'linewidth',1);
    
    legstr{1,j}=[Names{vec(Order(j)),1}];
    plot_legstr{1,j}=[region,Names{vec(Order(j)),1}];
end

% The same color used for the (n,m) species in region S11 is used in the
% S22 region.

if complete==2
    n=0;
    for j=1:length(center_22)
        if center_22(j)>0
            if j==1
                wL(j+length(Cent)-n)=fwhm_22(j)*x0(4*(length(Cent)-n+j)-1)/(0.5436+sqrt(x0(4*(length(Cent)-n+j))^2+.2166));
                wg(j+length(Cent)-n)=wL(j+length(Cent)-n)*x0(4*(length(Cent)-n+j));
            else
                wL(j+length(Cent)-n)=fwhm_22(j)*x0(4*(length(Cent)-n+j)-1)/(0.5436+sqrt((x0(4*(length(Cent)-n+j))*x0(4*(length(Cent)-n+1)))^2+.2166));
                wg(j+length(Cent)-n)=wL(j+length(Cent)-n)*x0(4*(length(Cent)-n+j))*x0(4*(length(Cent)-n+1));
            end
            
            XX=sqrt(log(2))*(Lamb-center_22(j)-x0(4*(length(Cent)-n+j)-2))/wg(j+length(Cent)-n);
            YY=sqrt(log(2))*wL(j+length(Cent)-n)/wg(j+length(Cent)-n);
            if j==1
                y(:,j+length(Cent)-n)=x0(4*j-3)/x0(4*(j+length(Cent)-n)-3).*real(complexErrorFunction(XX,YY))./real(complexErrorFunction(0,YY));
            else
                y(:,j+length(Cent)-n)=x0(4*j-3)/(x0(4*(1+length(Cent)-n)-3)*x0(4*(j+length(Cent)-n)-3)).*real(complexErrorFunction(XX,YY))./real(complexErrorFunction(0,YY));
            end
            
            plot(lambda,y(:,length(Cent)-n+j),'color',colors_n(Order(j),:),'linewidth',1);
            
            legstr{1,length(legstr)+1}=[Names{vec(Order(j)),1}];
            plot_legstr{1,length(plot_legstr)+1}=['S_2_2-',Names{vec(Order(j)),1}];
        else
            n=n+1;
        end
    end
end

area=zeros(2*length(cc),1);

% A different color is assigned to exciton phonon sidebands (EPS) compared
% to the color assigned to an (n,m) species.

for j=1:length(cc)
    area(j)=x0(4*cc(j)-3)*wg(cc(j))*sqrt(pi/log(2))/real(complexErrorFunction(0,sqrt(log(2))*wL(cc(j))/wg(cc(j))));
    if j==1
        y(:,ll/4+j)=area(j)*(0.017+0.1/diam_n(Order(cc(j)))+x0(ll+1))*exp(-2*(Lamb-h*c/(1e-9*(h*c/((Cent(cc(j))+x0((cc(j)*4-2)))*1e-9)+0.2+x0(ll+2)))).^2/(FWHM_g*x0(ll+3)/sqrt(log(4)))^2)/(FWHM_g*x0(ll+3)*sqrt(pi/(2*log(4))));
    else
        y(:,ll/4+j)=area(j)*(0.017+0.1/diam_n(Order(cc(j)))+x0(ll+1))*exp(-2*(Lamb-h*c/(1e-9*(h*c/((Cent(cc(j))+x0((cc(j)*4-2)))*1e-9)+0.2+x0(ll+2+2*(j-1))))).^2/(FWHM_g*x0(ll+3)*x0(ll+3+2*(j-1))/sqrt(log(4)))^2)/(FWHM_g*x0(ll+3)*x0(ll+3+2*(j-1))*sqrt(pi/(2*log(4))));
    end
    
    plot(lambda,y(:,ll/4+j),'color',[colors_n(Order(cc(j)),1),.8,colors_n(Order(cc(j)),3)],'linewidth',1);
    
    legstr{1,length(legstr)+1}=[Names{vec(Order(cc(j))),1}];
    plot_legstr{1,length(plot_legstr)+1}=['EPS-',region,Names{vec(Order(cc(j))),1}];
end

% The same color used for the EPS in region S11 is taken for the EPS in
% S22.

if complete==2 && isempty(cc2)==0
    for j=1:length(cc2)        
        if j==1
            area(j)=x0(4*cc2(j)-3)/x0(4*(cc2(j)+length(Cent)-n)-3)*wg(cc2(j)+length(Cent)-n)*sqrt(pi/log(2))/real(complexErrorFunction(0,sqrt(log(2))*wL(cc2(j)+length(Cent)-n)/wg(cc2(j)+length(Cent)-n)));
            y(:,2*length(Cent)-n+cc_j+j)=area(j)*(0.017+0.1/diam_n(Order(cc2(j)))+x0(ll+1+cc_len))*exp(-2*(Lamb-h*c/(1e-9*(h*c/((center_22(cc2(j))+x0((cc2(j)*4-2+length(Cent)*4)))*1e-9)+0.2+x0(ll+2+cc_len)))).^2/(FWHM_g*x0(ll+3+cc_len)/sqrt(log(4)))^2)/(FWHM_g*x0(ll+3+cc_len)*sqrt(pi/(2*log(4))));
        else
            area(j)=x0(4*cc2(j)-3)/(x0(4*(cc2(1)+length(Cent)-n)-3)*x0(4*(cc2(j)+length(Cent)-n)-3))*wg(cc2(j)+length(Cent)-n)*sqrt(pi/log(2))/real(complexErrorFunction(0,sqrt(log(2))*wL(cc2(j)+length(Cent)-n)/wg(cc2(j)+length(Cent)-n)));
            y(:,2*length(Cent)-n+cc_j+j)=area(j)*(0.017+0.1/diam_n(Order(cc2(j)))+x0(ll+1+cc_len))*exp(-2*(Lamb-h*c/(1e-9*(h*c/((center_22(cc2(j))+x0((cc2(j)*4-2+length(Cent)*4)))*1e-9)+0.2+x0(ll+2+2*(j-1)+cc_len)))).^2/(FWHM_g*x0(ll+3+cc_len)*x0(ll+3+2*(j-1)+cc_len)/sqrt(log(4)))^2)/(FWHM_g*x0(ll+3+cc_len)*x0(ll+3+2*(j-1)+cc_len)*sqrt(pi/(2*log(4))));
        end
        
        plot(lambda,y(:,2*length(Cent)-n+cc_j+j),'color',[colors_n(Order(cc2(j)),1),.8,colors_n(Order(cc2(j)),3)],'linewidth',1);
        
        legstr{1,length(legstr)+1}=[Names{vec(Order(cc2(j))),1}];
        plot_legstr{1,length(plot_legstr)+1}=['EPS-S_2_2-',Names{vec(Order(cc2(j))),1}];
    end
end

% Plot the absorption spectrum supplied by the user, the calculated
% absorption spectrum based on the sum of all fits and the chi-square
% value.

plot(lambda,Absorption,'k','linewidth',2)
plot(lambda,sum(y,2),'r','linewidth',2)

legend([plot_legstr,{'Measured Spectrum'},{'Calculated Spectrum'}])

title([title_name,': nSSE = ',num2str(sse)],'FontSize',26)

xlabel('Wavelength (nm)','FontSize',24)
ylabel('Absorption (a.u.)','FontSize',24)
set(gca,'LineWidth',2)
set(gca,'FontSize',22)


%..........................................................................


%--------------------------------------------------------------------------
% Section 9 - "plot_data_m"
%--------------------------------------------------------------------------


function [plot_legstr,y]=plot_data_m(lambda,Absorption,colors_22,colors_m,center_22,M11n,phonon_pos2,x_fit,FWHM_22,fwhm_nm_m11_n,legstr_22,w,CNTs_n,diam_22,h,c,FWHM_g_22,sse)

% Based on the line-shape chosen by the user previously and the results of
% the fit performed in Section 7 (stored in "x_fit"), the individual (n,m)
% species are plotted along with the calculated spectrum (sum of individual
% (n,m) fits) and the absorption spectrum initially loaded by the user. If
% Lorentzian or Gaussian line-shapes were defined, they are plotted with
% "plot_fit_m" defined in Section 9.1. If Voigtian line-shapes were chosen,
% they are plotted with "plot_fit_voigt_m" defined in Section 9.2.

switch w
    case {1,2}
        [plot_legstr,y]=plot_fit_m(lambda,Absorption,colors_22,colors_m,center_22,M11n,phonon_pos2,x_fit,FWHM_22,fwhm_nm_m11_n,legstr_22,w,CNTs_n,diam_22,h,c,FWHM_g_22,sse);
    case 3
        [plot_legstr,y]=plot_fit_voigt_m(lambda,center_22,x_fit,colors_22,legstr_22,M11n,fwhm_nm_m11_n,phonon_pos2,diam_22,h,c,FWHM_g_22,Absorption,sse,CNTs_n,colors_m);
        
end


%--------------------------------------------------------------------------
% Section 9.1 - "plot_fit_m"
%--------------------------------------------------------------------------


function [plot_legstr,y]=plot_fit_m(lambda,Absorption,colors_22,colors_m,center_22,M11n,phonon_pos2,x,FWHM_22,fwhm_nm_m11_n,legstr_22,w,CNTs_n,diam_22,h,c,FWHM_g_22,sse)

% This subfunction is used to plot the result of a Gaussian or Lorentzian
% fit. The structure of the code follows the rountine outlined in Section
% 7.1 "generate_string_m".

y=zeros(length(lambda),length(center_22)-sum(center_22==0)+length(M11n)+length(phonon_pos2));

figure('units','normalized','outerposition',[0 0 1 1]), hold on

n=0;

for j=1:length(center_22)
    if center_22(j)>0
        if w==1 % Lorentzian
            y(:,j) = x((j-n)*3-2)./(1+((lambda(1:end)-center_22(j)-x((j-n)*3-1))/(0.5*FWHM_22(j)*x((j-n)*3))).^2);
            title_name='Lorentzian Fit with ';
        else % Gaussian
            y(:,j) = x((j-n)*3-2)*exp(-log(2)*((lambda(1:end)-center_22(j)-x((j-n)*3-1))/(0.5*FWHM_22(j)*x((j-n)*3))).^2);
            title_name='Gaussian Fit with ';
        end
        
        plot(lambda,y(:,(j-n)),'color',colors_22(j,:),'linewidth',1)
        plot_legstr{(j-n)}=legstr_22{j-n};
    else
        n=n+1;
    end
end

for j=1+length(center_22):length(M11n)+length(center_22)
    
    if w==1 % Lorentzian
        y(:,(j-n)) = x((j-n)*3-2)./(1+((lambda(1:end)-M11n(j-length(center_22))-x((j-n)*3-1))/(0.5*fwhm_nm_m11_n(j-length(center_22))*x((j-n)*3))).^2);
    else % Gaussian
        y(:,(j-n)) = x((j-n)*3-2)*exp(-log(2)*((lambda(1:end)-M11n(j-length(center_22))-x((j-n)*3-1))/(0.5*fwhm_nm_m11_n(j-length(center_22))*x((j-n)*3))).^2);
    end
    
    plot(lambda,y(:,(j-n)),'color',colors_m(j-length(center_22),:),'linewidth',1)
    plot_legstr{(j-n)}=['M_1_1 - (',num2str(CNTs_n(j-length(center_22),1)),',',num2str(CNTs_n(j-length(center_22),2)),')'];
end

if isempty(M11n)==1
    x0_len=(length(center_22)-n)*3;
else
    x0_len=(j-n)*3;
end

for j=1:length(phonon_pos2)
    
    % Fit parameters:
    % x(1)=shift of diam vs. area ratio curve from Perebeinos et al.
    % PRL 94, 027402 (2005)
    % x(2)=shift of center around 0.2 eV above associated S11 or S22 peak
    % x(3)=Broadening of FWHM
    
    if j==1
        
        if w==1 % Lorentzian
            y(:,j+length(center_22)-n+length(M11n)) = x(phonon_pos2(j)*3-2)*pi*FWHM_22(phonon_pos2(j))*x(phonon_pos2(j)*3)/2*(0.017+0.1/diam_22(phonon_pos2(j)) + x(x0_len+1))*exp(-2*(lambda(1:end)-h*c/(1e-9*(h*c/((center_22(phonon_pos2(j))+x(phonon_pos2(j)*3-1))*1e-9)+0.2+x(x0_len+2)))).^2/(FWHM_g_22(j)*x(x0_len+3)/sqrt(log(4)))^2)/((FWHM_g_22(j)*x(x0_len+3))*sqrt(pi/(2*log(4))));
        else % Gaussian
            y(:,j+length(center_22)-n+length(M11n)) = x(phonon_pos2(j)*3-2)*FWHM_22(phonon_pos2(j))*x(phonon_pos2(j)*3)*sqrt(pi/(2*log(4)))*(0.017+0.1/diam_22(phonon_pos2(j)) + x(x0_len+1))*exp(-2*(lambda(1:end)-h*c/(1e-9*(h*c/((center_22(phonon_pos2(j))+x(phonon_pos2(j)*3-1))*1e-9)+0.2+x(x0_len+2)))).^2/(FWHM_g_22(j)*x(x0_len+3)/sqrt(log(4)))^2)/((FWHM_g_22(j)*x(x0_len+3))*sqrt(pi/(2*log(4))));
        end
        
    else
        
        % Fit parameters:
        % x(2+2*(j-1)) = shift of center around 0.2 eV above
        % associated S11
        % x(3+2*(j-1)) = Percentage of FWHM broadening of first
        % peak - can vary between 90 and 110 %.
        
        if w==1 % Lorentzian
            y(:,j+length(center_22)-n+length(M11n)) = x(phonon_pos2(j)*3-2)*pi*FWHM_22(phonon_pos2(j))*x(phonon_pos2(j)*3)/2*(0.017+0.1/diam_22(phonon_pos2(j)) + x(x0_len+1))*exp(-2*(lambda(1:end)-h*c/(1e-9*(h*c/((center_22(phonon_pos2(j))+x(phonon_pos2(j)*3-1))*1e-9)+0.2+x(x0_len+4+2*(j-2))))).^2/(FWHM_g_22(j)*x(x0_len+3)*x(x0_len+5+2*(j-2))/sqrt(log(4)))^2)/((FWHM_g_22(j)*x(x0_len+3)*x(x0_len+5+2*(j-2)))*sqrt(pi/(2*log(4))));
        else % Gaussian
            y(:,j+length(center_22)-n+length(M11n)) = x(phonon_pos2(j)*3-2)*FWHM_22(phonon_pos2(j))*x(phonon_pos2(j)*3)*sqrt(pi/(2*log(4)))*(0.017+0.1/diam_22(phonon_pos2(j)) + x(x0_len+1))*exp(-2*(lambda(1:end)-h*c/(1e-9*(h*c/((center_22(phonon_pos2(j))+x(phonon_pos2(j)*3-1))*1e-9)+0.2+x(x0_len+4+2*(j-2))))).^2/(FWHM_g_22(j)*x(x0_len+3)*x(x0_len+5+2*(j-2))/sqrt(log(4)))^2)/((FWHM_g_22(j)*x(x0_len+3)*x(x0_len+5+2*(j-2)))*sqrt(pi/(2*log(4))));
        end
    end
    
    plot(lambda,y(:,j+length(center_22)-n+length(M11n)),'color',[colors_22(phonon_pos2(j),1),.8,colors_22(phonon_pos2(j),3)],'linewidth',1)
    plot_legstr{j+length(center_22)-n+length(M11n)}=['EPS-',legstr_22{phonon_pos2(j)}];
    
end

plot(lambda,Absorption,'k','linewidth',2)

plot(lambda,sum(y,2),'r','linewidth',2)

legend([plot_legstr,{'Measured Spectrum'},{'Calculated Spectrum'}])

title([title_name,': nSSE = ',num2str(sse)],'FontSize',26)

xlabel('Wavelength (nm)','FontSize',24)
ylabel('Absorption (a.u.)','FontSize',24)
set(gca,'LineWidth',2)
set(gca,'FontSize',22)


%--------------------------------------------------------------------------
% Section 9.2 - "plot_fit_voigt_m"
%--------------------------------------------------------------------------


function [plot_legstr,y]=plot_fit_voigt_m(lambda,center_22,x0,colors_22,legstr_22,M11n,fwhm_nm_m11_n,phonon_pos2,diam_22,h,c,FWHM_g_22,Absorption,sse,CNTs_n,colors_m)

Lamb=lambda;

figure('units','normalized','outerposition',[0 0 1 1]), hold on

% Calculate Voigtian line-shape for S22 region

n=0;

for j=1:length(center_22)
    if center_22(j)>0
        wv(j-n)=x0(4*(j-n)-1);
        wL(j-n)=wv(j-n)/(0.5436+sqrt(x0(4*(j-n))^2+.2166));
        wg(j-n)=wL(j-n)*x0(4*(j-n));
        XX=sqrt(log(2))*(Lamb-x0(4*(j-n)-2))/wg(j-n);
        YY=sqrt(log(2))*wL(j-n)/wg(j-n);
        y(:,j-n)=x0(4*(j-n)-3).*real(complexErrorFunction(XX,YY))./real(complexErrorFunction(0,YY));
        
        title_name='Voigtian Fit with';
        
        plot(lambda,y(:,j-n),'color',colors_22(j,:),'linewidth',1)
        plot_legstr{j-n}=legstr_22{j-n};
    else
        n=n+1;
    end
end

% Calculate Voigtian line-shape for M11 nanotubes

for j=1:length(M11n)
    if j==1
        wL(j+length(center_22)-n)=0.5*fwhm_nm_m11_n(j)*x0(4*(length(center_22)+j-n)-1)/(0.5436+sqrt(x0(4*(length(center_22)+j-n))^2+.2166));
        wg(j+length(center_22)-n)=wL(j+length(center_22)-n)*x0(4*(length(center_22)+j-n));
    else
        wL(j+length(center_22)-n)=0.5*fwhm_nm_m11_n(j)*x0(4*(length(center_22)+j-n)-1)/(0.5436+sqrt((x0(4*(length(center_22)+j-n))*x0(4*(length(center_22)+1-n)))^2+.2166));
        wg(j+length(center_22)-n)=wL(j+length(center_22)-n)*x0(4*(length(center_22)+j-n))*x0(4*(length(center_22)+1-n));
    end
    
    XX=sqrt(log(2))*(Lamb-M11n(j)-x0(4*(length(center_22)+j-n)-2))/wg(j+length(center_22)-n);
    YY=sqrt(log(2))*wL(j+length(center_22)-n)/wg(j+length(center_22)-n);
    y(:,j+length(center_22)-n)=x0(4*(j+length(center_22)-n)-3).*real(complexErrorFunction(XX,YY))./real(complexErrorFunction(0,YY));
    
    plot(lambda,y(:,j+length(center_22)-n),'color',colors_m(j,:),'linewidth',1)
    plot_legstr{j+length(center_22)-n}=['M_1_1 - (',num2str(CNTs_n(j,1)),',',num2str(CNTs_n(j,2)),')'];
end

% Calculate the line-shape of the Gaussian exciton phonon sideband for the
% S11 or S22 region based on the area of the corresponding nanotube peak.

area=zeros(2*length(phonon_pos2),1);

ll=(length(center_22)+length(M11n)-n)*4;

for j=1:length(phonon_pos2)
    area(j)=trapz(lambda,y(:,phonon_pos2(j)));
    if j==1
        y(:,ll/4+j)=area(j)*(0.017+0.1/diam_22(phonon_pos2(j))+x0(ll+1))*exp(-2*(Lamb-h*c/(1e-9*(h*c/((x0((phonon_pos2(j)*4-2)))*1e-9)+0.2+x0(ll+2)))).^2/(FWHM_g_22(j)*x0(ll+3)/sqrt(log(4)))^2)/(FWHM_g_22(j)*x0(ll+3)*sqrt(pi/(2*log(4))));
    else
        y(:,ll/4+j)=area(j)*(0.017+0.1/diam_22(phonon_pos2(j))+x0(ll+1))*exp(-2*(Lamb-h*c/(1e-9*(h*c/((x0((phonon_pos2(j)*4-2)))*1e-9)+0.2+x0(ll+2+2*(j-1))))).^2/(FWHM_g_22(j)*x0(ll+3)*x0(ll+3+2*(j-1))/sqrt(log(4)))^2)/(FWHM_g_22(j)*x0(ll+3)*x0(ll+3+2*(j-1))*sqrt(pi/(2*log(4))));
    end
    
    plot(lambda,y(:,j+length(center_22)+length(M11n)-n),'color',[colors_22(phonon_pos2(j),1),.8,colors_22(phonon_pos2(j),3)],'linewidth',1)
    plot_legstr{j+length(center_22)+length(M11n)-n}=['EPS-',legstr_22{phonon_pos2(j)}];
end

plot(lambda,Absorption,'k','linewidth',2)

plot(lambda,sum(y,2),'r','linewidth',2)

legend([plot_legstr,{'Measured Spectrum'},{'Calculated Spectrum'}])

title([title_name,': nSSE = ',num2str(sse)],'FontSize',26)

xlabel('Wavelength (nm)','FontSize',24)
ylabel('Absorption (a.u.)','FontSize',24)
set(gca,'LineWidth',2)
set(gca,'FontSize',22)


%..........................................................................


%--------------------------------------------------------------------------
% Section 10 - "check_fit"
%--------------------------------------------------------------------------


function [nn,complete,tt,d_pf,d_eps,d_hr,TO,pa,d_v,w_n,dp]=check_fit(data_peakfit,data_EPS_peakfit,data_heightratio,Turnover,nn,fit_region,complete,PeakAssignment,data_voigt,w,fid,protocol,doping)
global Absorption s_ud fwhm_gauss

% In the event that the user chose to fit the S11 region, they are
% presented with the option to fit the entire region (previously defined
% during background subtraction in "get_data") based on the (n,m)
% distribution of S11.
%
% NOTE: if S22 of an (n,m) species is outside the range defined during
% background subtraction it will be ignored during the fit!
%
% The user is presented with the fit and asked whether they are satisfied
% with its quality. If not, the user is given various options to improve
% the quality of the fit.

tt='n';
d_pf=data_peakfit;
d_eps=data_EPS_peakfit;
d_hr=data_heightratio;
TO=Turnover;
pa=PeakAssignment;
d_v=data_voigt;
w_n=w;
dp=doping;

if nn==1 && fit_region==1 % After first fit and only if S11 was fitted
    complete=input('Do you also want fit the S22 region? \nIf "yes", the quality of the S11 fit might change \nand the computation time might increase drastically (y/n)? ','s');
    complete=sscanf(complete,'%s');
end

if complete == 'n' || fit_region==2 || nn==2
    
    tt=input('Are you satisfied with the fit (y/n)? ','s');
    
    if tt=='n'
        line_prof='\n(1) - Change the line profile.';
        shift='\n(2) - Shift measured absorption spectrum up or down.';
        change_bc='\n(3) - Change boundary conditions for the fit.';
        change_eps='\n(4) - Change boundary conditions for exciton phonon sideband, if there are any.';
        change_height='\n(5) - Change the boundary conditions for the height assignment.';
        
        ASK=['Please choose what you want to adjust:',line_prof,shift,change_bc,change_eps,change_height];
        
        height_ratio='Change the peak intensity ratio to constrain the S22 fit based on the S11 fit.';
        voigt_fit='Change parameters of Voigt fit.';
        redo_s11='Redo Peak Assignment for S11.';
        redo_s22='Redo Peak Assignment for S22.';
        dop='Consider doping effects (S22 might become more intense than S11).';
        
        if fit_region==1
            if complete=='n' % Fit only the S11 region
                switch w
                    case {1,2}
                        ASK=[ASK,'\n(6) - ',redo_s11,'\n'];
                        oo=input(ASK,'s');
                        oo=sscanf(oo,'%i');
                        if oo==6
                            oo=8;
                        end
                    case 3
                        ASK=[ASK,'\n(6) - ',voigt_fit,'\n(7) - ',redo_s11,'\n'];
                        oo=input(ASK,'s');
                        oo=sscanf(oo,'%i');
                        if oo==6 || oo==7
                            oo=oo+1;
                        end
                end
            else
                switch w
                    case {1,2}
                        ASK=[ASK,'\n(6) - ',height_ratio,'\n(7) - ',redo_s11,'\n(8) - ',redo_s22,'\n'];
                        oo=input(ASK,'s');
                        oo=sscanf(oo,'%i');
                        if oo==7 || oo==8
                            oo=oo+1;
                        end
                    case 3
                        ASK=[ASK,'\n(6) - ',height_ratio,'\n(7) - ',voigt_fit,'\n(8) - ',redo_s11,'\n(9) - ',redo_s22,'\n'];
                        oo=input(ASK,'s');
                        oo=sscanf(oo,'%i');
                end
            end
        end
        
        if fit_region==2
            switch w
                case {1,2}
                    ASK=[ASK,'\n(6) - ',redo_s22,'\n(7) - ',dop,'\n'];
                    oo=input(ASK,'s');
                    oo=sscanf(oo,'%i');
                    if oo==6
                        oo=9;
                    else
                        oo=10;
                    end
                case 3
                    ASK=[ASK,'\n(6) - ',voigt_fit,'\n(7) - ',redo_s22,'\n(8) - ',dop,'\n'];
                    oo=input(ASK,'s');
                    oo=sscanf(oo,'%i');
                    if oo==6
                        oo=oo+1;
                    else
                        oo=oo+2;
                    end
            end
        end
        
        switch oo
            case 1
                w_n=input('Please choose a line-shape (exciton phonon sidebands are always fitted with a Gaussian line-shape):\n(1) - Lorentzian\n(2) - Gaussian\n(3) - Voigt\n ','s');
                w_n=sscanf(w_n,'%i');
            case 2
                s_ud=input('Please enter shift in percentage of the maximum value,\n e.g. a down shift of -1% or an up-shift of 3%: ','s');
                s_ud=sscanf(s_ud,'%i');
                s_ud=s_ud/100;
                s_ud=max(Absorption)*s_ud;
                Absorption=Absorption+s_ud;
            case 3
                cnames = {'Lower Limit','Initial Guess','Upper Limit'};
                rnames = {'FWHM','Shift of Center (constrained peak)','Shift of Center (un-constrained peak)', 'Height of First Peak', 'Height Guess 1','Height Guess 2'};
                vals=[true true true true true true true true true true true true true true true true true true];
                
                while 1>0
                    f=figure('Position',[100 100 800 150]); t = uitable('Position',[0 0 800 150],'Data', data_peakfit,'ColumnName',cnames,'RowName',rnames,'ColumnEditable', vals);
                    vv=version('-release');
                    year=str2double(vv(1:4));
                    if year>=2014
                        t.Position(3) = t.Extent(3);
                        t.Position(4) = t.Extent(4);
                    end
                    
                    opt=input('Please change the variables in the table and press d if done: ','s');
                    if opt=='d'
                        d_pf=t.Data;
                        close(f)
                        break
                    end
                end
                
            case 4
                cnames = {'Lower Limit','Initial Guess','Upper Limit'};
                rnames = {'FWHM','Correction Factor f1','Shift of Center (eV)','Change of FWHM of Second EPS'};
                vals=[true true true true true true true true true true true true];
                
                while 1>0
                    f=figure('Position',[100 100 800 150]); t = uitable('Position',[0 0 800 150],'Data', data_EPS_peakfit,'ColumnName',cnames,'RowName',rnames,'ColumnEditable', vals);
                    vv=version('-release');
                    year=str2double(vv(1:4));
                    if year>=2014
                        t.Position(3) = t.Extent(3);
                        t.Position(4) = t.Extent(4);
                    end
                    
                    opt=input('Please change the variables in the table and press d if done: ','s');
                    if opt=='d'
                        d_eps=t.Data;
                        close(f)
                        fwhm_gauss=d_eps(1,2);
                        break
                    end
                end
                
            case 5
                
                cnames = {};
                rnames = {'Turnover (nm)','FWHM Broadening Upper Limit','FWHM Broadening Lower Limit','Height Guess (1 or 2)'};
                vals=[true true true true];
                
                while 1>0
                    f=figure('Position',[100 100 800 150]); t = uitable('Position',[0 0 800 150],'Data', Turnover,'ColumnName',cnames,'RowName',rnames,'ColumnEditable', vals);
                    vv=version('-release');
                    year=str2double(vv(1:4));
                    if year>=2014
                        t.Position(3) = t.Extent(3);
                        t.Position(4) = t.Extent(4);
                    end
                    
                    opt=input('Please change the variables in the table and press d if done: ','s');
                    if opt=='d'
                        TO=t.Data;
                        close(f)
                        break
                    end
                end
                
            case 6
                
                cnames = {'Lower Limit','Initial Guess','Upper Limit'};
                rnames = {'Height Ratio','Change of Initial Height Ratio'};
                vals=[true true true true true true];
                
                while 1>0
                    f=figure('Position',[100 100 800 150]); t = uitable('Position',[0 0 800 150],'Data', data_heightratio,'ColumnName',cnames,'RowName',rnames,'ColumnEditable', vals);
                    vv=version('-release');
                    year=str2double(vv(1:4));
                    if year>=2014
                        t.Position(3) = t.Extent(3);
                        t.Position(4) = t.Extent(4);
                    end
                    
                    opt=input('Please change the variables in the table and press d if done: ','s');
                    if opt=='d'
                        d_hr=t.Data;
                        close(f)
                        break
                    end
                end
                
            case 7
                
                cnames = {};
                rnames = {'Maximum Lorentzian HWHM (nm)'; 'Ratio of Gaussian / Lorentzian HWHM'; 'Ratio of HWHM_G / HWHM_L Lower Limit'; 'Ratio of HWHM_G / HWHM_L Upper Limit'; 'Change of HWHM_G / HWHM_L'};
                vals=[true true true true true];
                
                while 1>0
                    f=figure('Position',[100 100 800 150]); t = uitable('Position',[0 0 800 150],'Data', data_voigt,'ColumnName',cnames,'RowName',rnames,'ColumnEditable', vals);
                    vv=version('-release');
                    year=str2double(vv(1:4));
                    if year>=2014
                        t.Position(3) = t.Extent(3);
                        t.Position(4) = t.Extent(4);
                    end
                    
                    opt=input('Please change the variables in the table and press d if done: ','s');
                    if opt=='d'
                        d_v=t.Data;
                        close(f)
                        break
                    end
                end
                
            case 8
                
                pa(1)=1;
                
            case 9
                
                pa(2)=1;
                
            case 10
                opt=input('Do you want to consider doping effects (y/n)?: ','s');
                if opt=='y'
                    dp=1;
                else
                    dp=0;
                end
        end
    else
        if protocol=='y'
            array{1}=['Chosen line profile: ',num2str(w)];
            array{2}=['Shift of absorption spectrum (in %): ',num2str(s_ud)];
            array{3}=['Boundary condition (BC) of FWHM: ',num2str(data_peakfit(1,:))];
            array{4}=['BC of Shift of Center (constrained peak): ',num2str(data_peakfit(2,:))];
            array{5}=['BC of Shift of Center (un-constrained peak): ',num2str(data_peakfit(3,:))];
            array{6}=['BC of Height of First Peak: ',num2str(data_peakfit(4,:))];
            array{7}=['BC of Height Guess 1: ',num2str(data_peakfit(5,:))];
            array{8}=['BC of Height Guess 2: ',num2str(data_peakfit(6,:))];
            array{9}=['BC of EPS FWHM: ',num2str(data_EPS_peakfit(1,:))];
            array{10}=['BC of EPS Correction Factor f1: ',num2str(data_EPS_peakfit(2,:))];
            array{11}=['BC of EPS Shift of Center (eV): ',num2str(data_EPS_peakfit(3,:))];
            array{12}=['BC of Change of FWHM of Second EPS: ',num2str(data_EPS_peakfit(4,:))];
            array{13}=['BC of Turnover value (nm): ',num2str(Turnover(1))];
            array{14}=['BC of FWHM Broadening Upper Limit: ',num2str(Turnover(2))];
            array{15}=['BC of FWHM Broadening Lower Limit: ',num2str(Turnover(3))];
            array{16}=['Height Guess: ',num2str(Turnover(4))];
            
            if fit_region==1
                if complete=='y'
                    array{17}='Fitted region: S11 and S22';
                    array{18}=['BC of Allowed Change of Height Ratio: ',num2str(data_heightratio(2,:))];
                    array{19}=['Doping effects considered (y=1/n=0)? ',num2str(dp)];
                else
                    array{17}='Fitted region: S11';
                end
            else
                array{17}='Fitted region: S22';
            end
            
            if w==3
                array{20}=['Upper limit of Lorentzian HWHM: ',num2str(data_voigt(1))];
                array{21}=['BC of ratio of Gaussian/Lorentzian HWHM: ',num2str([data_voigt(3) data_voigt(2) data_voigt(4) data_voigt(5)])];
            end
            
            fprintf(fid,'%s \n',array{:});
            
            fclose('all');
        end
    end
    
else
    nn=2;
end


%..........................................................................


%--------------------------------------------------------------------------
% Section 11 - "check_fit_m"
%--------------------------------------------------------------------------


function [tt,d_S22,d_EPS,d_M11,d_vs22,d_vm,mp,dp]=check_fit_m(data_S22,data_EPS,data_M11,Absorption,data_voigt_s22,data_voigt,w,fid,protocol,doping)
global s_ud

% The user is presented with the fit and asked whether they are satisfied
% with its quality. If not, the user is given various options to improve
% the quality of the fit.

tt=input('Are you satisfied with the fit (y/n)? ','s');
d_S22=data_S22;
d_EPS=data_EPS;
d_M11=data_M11;
d_vs22=data_voigt_s22;
d_vm=data_voigt;
mp=0;
dp=doping;

if tt=='n'
    oo=input('Please choose what you want to adjust:\n(1) - Shift measured absorption spectrum up or down.\n(2) - Change boundary conditions for the fit of S22.\n(3) - Change boundary conditions for exciton phonon sideband, if there are any.\n(4) - Change the boundary conditions for the fit of M11.\n(5) - Change assignment of metallic nanotubes.\n(6) - Consider doping effects (S22 might become more intense than S11)?\n','s');
    oo=sscanf(oo,'%i');
    
    switch oo
        case 1
            s_ud=input('Please enter shift in percentage of the maximum value,\n e.g. a down shift of -1% or an up-shift of 3%: ','s');
            s_ud=sscanf(s_ud,'%i');
            s_ud=s_ud/100;
            s_ud=max(Absorption)*s_ud;
        case 2
            if w==3
                cnames = {'Lower Limit','Initial Guess','Upper Limit'};
                rnames = {'FWHM','Shift of Center','Height', 'Gaussian / Lorentzian HWHM'};
                vals=[true true true true true true true true true true true true];
                
                while 1>0
                    f=figure('Position',[100 100 800 150]); t = uitable('Position',[0 0 800 150],'Data', [data_S22; data_voigt_s22],'ColumnName',cnames,'RowName',rnames,'ColumnEditable', vals);
                    vv=version('-release');
                    year=str2double(vv(1:4));
                    if year>=2014
                        t.Position(3) = t.Extent(3);
                        t.Position(4) = t.Extent(4);
                    end
                    
                    opt=input('Please change the variables in the table and press d if done: ','s');
                    if opt=='d'
                        d_S22=t.Data(1:3,:);
                        d_vs22=t.Data(4,:);
                        close(f)
                        break
                    end
                end
            else
                cnames = {'Lower Limit','Initial Guess','Upper Limit'};
                rnames = {'FWHM','Shift of Center','Height'};
                vals=[true true true true true true true true true];
                
                while 1>0
                    f=figure('Position',[100 100 800 150]); t = uitable('Position',[0 0 800 150],'Data', data_S22,'ColumnName',cnames,'RowName',rnames,'ColumnEditable', vals);
                    vv=version('-release');
                    year=str2double(vv(1:4));
                    if year>=2014
                        t.Position(3) = t.Extent(3);
                        t.Position(4) = t.Extent(4);
                    end
                    
                    opt=input('Please change the variables in the table and press d if done: ','s');
                    if opt=='d'
                        d_S22=t.Data;
                        close(f)
                        break
                    end
                end
            end
            
        case 3
            cnames = {'Lower Limit','Initial Guess','Upper Limit'};
            rnames = {'FWHM','Shift of Center','Factor f1'};
            vals=[true true true true true true true true true];
            
            while 1>0
                f=figure('Position',[100 100 800 150]); t = uitable('Position',[0 0 800 150],'Data', data_EPS,'ColumnName',cnames,'RowName',rnames,'ColumnEditable', vals);
                vv=version('-release');
                year=str2double(vv(1:4));
                if year>=2014
                    t.Position(3) = t.Extent(3);
                    t.Position(4) = t.Extent(4);
                end
                
                opt=input('Please change the variables in the table and press d if done: ','s');
                if opt=='d'
                    d_EPS=t.Data;
                    close(f)
                    break
                end
                
            end
            
        case 4
            
            if w==3
                cnames = {'Lower Limit','Initial Guess','Upper Limit'};
                rnames = {'FWHM','Shift of Center','Height', 'Gaussian / Lorentzian HWHM'};
                vals=[true true true true true true true true true true true true true true true true];
                
                while 1>0
                    f=figure('Position',[100 100 800 150]); t = uitable('Position',[0 0 800 150],'Data', [data_M11(1,1:3) NaN;data_M11(2,1:3) NaN;data_M11(3,1:3) NaN; data_voigt(2:5)],'ColumnName',cnames,'RowName',rnames,'ColumnEditable', vals);
                    vv=version('-release');
                    year=str2double(vv(1:4));
                    if year>=2014
                        t.Position(3) = t.Extent(3);
                        t.Position(4) = t.Extent(4);
                    end
                    
                    opt=input('Please change the variables in the table and press d if done: ','s');
                    if opt=='d'
                        d_M11=t.Data(1:3,1:3);
                        d_vm=[data_voigt(1), t.Data(4,:)];
                        close(f)
                        break
                    end
                end
                
            else
                cnames = {'Lower Limit','Initial Guess','Upper Limit'};
                rnames = {'FWHM','Shift of Center','Height'};
                vals=[true true true true true true true true true];
                
                while 1>0
                    f=figure('Position',[100 100 800 150]); t = uitable('Position',[0 0 800 150],'Data', data_M11,'ColumnName',cnames,'RowName',rnames,'ColumnEditable', vals);
                    vv=version('-release');
                    year=str2double(vv(1:4));
                    if year>=2014
                        t.Position(3) = t.Extent(3);
                        t.Position(4) = t.Extent(4);
                    end
                    
                    opt=input('Please change the variables in the table and press d if done: ','s');
                    if opt=='d'
                        d_M11=t.Data;
                        close(f)
                        break
                    end
                end
            end
        case 5
            mp=1;
        case 6
            opt=input('Do you want to consider doping effects (y/n)?: ','s');
            if opt=='y'
                dp=1;
            else
                dp=0;
            end
    end
else
    if protocol=='y'
        array{1}=['Shift of S22 region (in %): ',num2str(s_ud)];
        array{2}=['Boundary Condition (BC) of EPS FWHM for S22 SWCNT transitions: ',num2str(data_EPS(1,:))];
        array{3}=['BC of EPS Shift of Center for S22 SWCNT transitions: ',num2str(data_EPS(2,:))];
        array{4}=['BC of EPS Factor f1 for S22 SWCNT transitions: ',num2str(data_EPS(3,:))];
        
        k=5;
        if w==3
            array{k}=['BC of Ratio of Gaussian/Lorentzian HWHM for S22 SWCNT transitions: ',num2str(data_voigt_s22)];
            array{k+1}=['BC of Ratio of Gaussian/Lorentzian HWHM for M11 SWCNT transitions: ',num2str([data_voigt(3) data_voigt(2) data_voigt(4)])];
            k=k+2;
        end
        array{k}=['BC of FWHM for S22 SWCNT transitions: ',num2str(data_S22(1,:))];
        array{k+1}=['BC of Shift of Center for S22 SWCNT transitions: ',num2str(data_S22(2,:))];
        array{k+2}=['BC of Height for S22 SWCNT transitions: ',num2str(data_S22(3,:))];
        array{k+3}=['BC of FWHM for M11 SWCNT transitions: ',num2str(data_M11(1,:))];
        array{k+4}=['BC of Shift of Center for M11 SWCNT transitions: ',num2str(data_M11(2,:))];
        array{k+5}=['BC of Height for M11 SWCNT transitions: ',num2str(data_M11(3,:))];
        array{k+6}=['Doping effects considered (y=1/n=0)? ',num2str(dp)];
        
        fprintf(fid,'%s \n',array{:});
        fclose('all');
    end
end


%..........................................................................


%--------------------------------------------------------------------------
% Section 12 - "save_variables_film"
%--------------------------------------------------------------------------


function [fwhm_nm_n,fwhm_nm_n22,diam_n,time_stamp]=save_variables_film(x_fit,center,sse,fwhm_nm_n,fwhm_nm_n22,diam_n,Order,legstr,colors_n,w,center_22,cc,cc2,fit_region,ts,protocol)

if isempty(ts)==1
    bla=now;
    time_stamp=datestr(bla);
    time_stamp=strrep(time_stamp,' ','_');
    time_stamp=strrep(time_stamp,':','-');
else
    time_stamp=ts;
end

str=['VariablesFilmFitting',time_stamp];
mkdir(str)
current=pwd;

if protocol=='y'
    
    string='%s ';
    for i=1:25
        string=[string, '%s '];
    end
    
    fid=fopen(['Protocol_',time_stamp,'.txt']);
    F = textscan(fid,string);
    fclose('all');
    
    delete(['Protocol_',time_stamp,'.txt']);
    
    bla=[F{:}];
    for i=1:size(bla,1)
        r_peaks1(i)=strcmp(strjoin(bla(i,1:5)), 'Removed peaks from S11 region:');
        r_peaks2(i)=strcmp(strjoin(bla(i,1:5)), 'Removed peaks from S22 region:');
    end
    
    idx_r_peaks1=find(r_peaks1==1);
    idx_r_peaks2=find(r_peaks2==1);
    
    bll=bla(1:6,:);
    
    if isempty(idx_r_peaks1)==0
        bll=[bll; bla(idx_r_peaks1(end):idx_r_peaks1(end)+3,:)];
    end
    if isempty(idx_r_peaks2)==0
        bll=[bll; bla(idx_r_peaks2(end):idx_r_peaks2(end)+3,:)];
    end
    if isempty(idx_r_peaks2)==1
        start=idx_r_peaks1(end)+4;
    else
        start=idx_r_peaks2(end)+4;
    end
    
    bll=[bll; bla(start:end,:)];
    
    for i=1:size(bll,1)
        array{i}=strjoin(bll(i,:));
    end
    
    fid=fopen(['Protocol_',time_stamp,'.txt'],'w');
    if fid<0
        fclose('all');
        clear fid
        fid=fopen(['Protocol_',time_stamp,'.txt'],'w');
    end
    fprintf(fid,'%s \n',array{:});
    fclose('all');
    
    tt=ispc;
    if tt==1
        movefile(['Protocol_',time_stamp,'.txt'],[current,['\',str]])
    else
        movefile(['Protocol_',time_stamp,'.txt'],[current,['/',str]])
    end
end

tt=ispc;
if tt==1
    oldD=cd([current,['\',str]]);
else
    oldD=cd([current,['/',str]]);
end

save('x_fit_solution','x_fit');
save('center','center')
save('center_22','center_22')
save('sse','sse')

fwhm_nm_n=fwhm_nm_n(Order);
fwhm_nm_n22=fwhm_nm_n22(Order);
save('fwhm_nm_n','fwhm_nm_n')
save('fwhm_nm_n22','fwhm_nm_n22')

diam_n=diam_n(Order);
save('diam_n','diam_n')

save('legstr','legstr')

save('colors_n','colors_n')

save('w','w')

save('fit_region','fit_region')

save('cc','cc')
save('cc2','cc2')

cd(oldD);


%..........................................................................


%--------------------------------------------------------------------------
% Section 13 - "save_data"
%--------------------------------------------------------------------------


function save_data(x_fit,w,fit_region,complete,center,cc,fwhm_nm_n,cc2,center_22,fwhm_nm_n22,h,c,colors_n,Order,diam_n,legstr,plot_legstr,nn,fwhm_gauss,time_stamp,path_length)

s11=0;

switch w
    case {1,2}
        switch fit_region
            case 1 % Store values for S11 region
                if isempty(complete)==1 || complete == 'n'
                    Results{1,1}='CNTs';
                    Results{1,2}='Diameter (nm)';
                    Results{1,3}='Center S11 (nm)';
                    Results{1,4}='FWHM S11 (nm)';
                    Results{1,5}='FWHM Broadening S11';
                    Results{1,6}='Height / Area-Frac S11';
                    Results{1,7}='Area S11';
                    Results{1,8}='Rel. Concentration (%) S11';
                    Results{1,9}='Concentration (g/l)';
                    s11=1;
                    x_len=0;
                else % Store values for entire region
                    Results{1,1}='CNTs';
                    Results{1,2}='Diameter (nm)';
                    Results{1,3}='Center S11 (nm)';
                    Results{1,4}='FWHM S11 (nm)';
                    Results{1,5}='FWHM Broadening S11';
                    Results{1,6}='Height / Area-Frac S11';
                    Results{1,7}='Area S11';
                    Results{1,8}='Center S22 (nm)';
                    Results{1,9}='FWHM S22 (nm)';
                    Results{1,10}='FWHM Broadening S22';
                    Results{1,11}='Height / Area-Frac S22';
                    Results{1,12}='Area S22';
                    Results{1,13}='Ratio of Height S11/S22';
                    Results{1,14}='Rel. Concentration (%) S11';
                    Results{1,15}='Rel. Concentration (%) S22';
                    Results{1,16}='Concentration (g/l)';
                    x_len=length(center_22);
                end
                Cent=center;
                FWHM=fwhm_nm_n;
                CC=cc;
            case 2 % Store values for S22 region
                Results{1,1}='CNTs';
                Results{1,2}='Diameter (nm)';
                Results{1,3}='Center S22 (nm)';
                Results{1,4}='FWHM S22 (nm)';
                Results{1,5}='FWHM Broadening S22';
                Results{1,6}='Height / Area-Frac S22';
                Results{1,7}='Area S22';
                Results{1,8}='Rel. Concentration (%) S22';
                x_len=0;
                Cent=center_22;
                FWHM=fwhm_nm_n22;
                CC=cc;
        end
        n=0;
        for i=1:length(Cent)+length(cc)+length(cc2)
            if i<=length(Order)
                Results{i+1,2}=diam_n(i);
                Results{i+1,3}=Cent(i)+x_fit(i*3-1);
                Results{i+1,4}=FWHM(i)*x_fit(i*3);
                Results{i+1,5}=x_fit(i*3);
                Results{i+1,6}=x_fit(i*3-2);
                if w==1 % Lorentzian
                    Results{i+1,7}=x_fit(i*3-2)*FWHM(i)*x_fit(i*3)*pi/2;
                    area(i)=x_fit(i*3-2)*FWHM(i)*x_fit(i*3)*pi/2;
                else % Gaussian
                    Results{i+1,7}=x_fit(i*3-2)*FWHM(i)*x_fit(i*3)*sqrt(pi/(2*log(4)));
                    area(i)=x_fit(i*3-2)*FWHM(i)*x_fit(i*3)*sqrt(pi/(2*log(4)));
                end
                
                if nn==2
                    if center_22(i)>0
                        Results{i+1,8}=center_22(i)+x_fit(length(center)*3+(i-n)*3-1);
                        Results{i+1,9}=fwhm_nm_n22(i)*x_fit(length(center)*3+(i-n)*3);
                        Results{i+1,10}=x_fit(length(center)*3+(i-n)*3);
                        if i==1
                            Results{i+1,11}=x_fit(i*3-2)/(x_fit(length(center)*3+(i-n)*3-2));
                        else
                            Results{i+1,11}=x_fit(i*3-2)/(x_fit(length(center)*3+(1-n)*3-2)*x_fit(length(center)*3+(i-n)*3-2));
                        end
                        if w==1 % Lorentzian
                            Results{i+1,12}=Results{i+1,11}*fwhm_nm_n22(i)*x_fit(length(center)*3+(i-n)*3)*pi/2;
                            area22(i)=Results{i+1,11}*fwhm_nm_n22(i)*x_fit(length(center)*3+(i-n)*3)*pi/2;
                        else % Gaussian
                            Results{i+1,12}=Results{i+1,11}*fwhm_nm_n22(i)*x_fit(length(center)*3+(i-n)*3)*sqrt(pi/(2*log(4)));
                            area22(i)=Results{i+1,11}*fwhm_nm_n22(i)*x_fit(length(center)*3+(i-n)*3)*sqrt(pi/(2*log(4)));
                        end
                        Results{i+1,13}=Results{i+1,6}/Results{i+1,11};
                    else
                        n=n+1;
                    end
                end
                Results{i+1,1}=legstr{i};
            else
                Results{i+1,1}=plot_legstr{i+x_len-n};
            end
        end
        
        for i=1:length(center)
            if nn==2
                Results{i+1,14}=area(i)/sum(area)*100;
                if Results{i+1,2}<1
                    Results{i+1,16}=Results{i+1,6}/(path_length*(2210+1950*exp((0.85-Results{i+1,2})/0.158)))*12.0107;
                else
                    Results{i+1,16}=Results{i+1,6}/(path_length*(0.183+0.163*exp((0.85-Results{i+1,2})/0.158)))*1e-3;
                end
                if area22(i)>0
                    Results{i+1,15}=area22(i)/sum(area22)*100;
                end
            else
                Results{i+1,8}=area(i)/sum(area)*100;
                if s11==1 % Only absorption cross-sections of S11 are given in Sanchez et al., "(n,m)-Specific Absorption Cross Sections of Single-Walled Carbon Nanotubes Measured by Variance Spectroscopy"
                    if Results{i+1,2}<1
                        Results{i+1,9}=Results{i+1,6}/(path_length*(2210+1950*exp((0.85-Results{i+1,2})/0.158)))*12.0107;
                    else
                        Results{i+1,9}=Results{i+1,6}/(path_length*(0.183+0.163*exp((0.85-Results{i+1,2})/0.158)))*1e-3;
                    end
                end
            end
        end
        
        if nn==2
            ll=6*length(center)-3*n;
        else
            ll=3*length(center);
        end
        
        kk=0; % if there is no EPS assigned to the S11 transition
        
        for j=1:length(CC) % Store only if there was a phonon sideband
            Results{length(colors_n)+1+j,3}=h*c/(1e-9*(h*c/((Cent(CC(j))+x_fit((CC(j)*3-1)))*1e-9)+0.2+x_fit(ll+2*j)));
            if j==1
                Results{length(colors_n)+1+j,4}=fwhm_gauss*x_fit(ll+3);
                Results{length(colors_n)+1+j,5}=x_fit(ll+3);
            else
                Results{length(colors_n)+1+j,4}=fwhm_gauss*x_fit(ll+3)*x_fit(ll+3+2*(j-1));
                Results{length(colors_n)+1+j,5}=x_fit(ll+3)*x_fit(ll+3+2*(j-1));
            end
            Results{length(colors_n)+1+j,6}=x_fit(ll+1);
            Results{length(colors_n)+1+j,7}=area(CC(j))*(0.017+0.1/diam_n(CC(j))+Results{length(colors_n)+1+j,6});
            kk=3+2*(j-1);
        end
        
        if nn==2
            for j=1:length(cc2) % Store only if there was a phonon sideband
                Results{length(colors_n)+1+j+length(CC),8}=h*c/(1e-9*(h*c/((center_22(cc2(j))+x_fit(length(center)*3+cc2(j)*3-1))*1e-9)+0.2+x_fit(6*length(center)-3*n+kk+2*j)));
                if j==1
                    Results{length(colors_n)+1+j+length(CC),9}=fwhm_gauss*x_fit(6*length(center)+3-3*n+kk);
                    Results{length(colors_n)+1+j+length(CC),10}=x_fit(6*length(center)+3-3*n+kk);
                else
                    Results{length(colors_n)+1+j+length(CC),9}=fwhm_gauss*x_fit(6*length(center)+3-3*n+kk)*x_fit(6*length(center)+3-3*n+kk+2*(j-1));
                    Results{length(colors_n)+1+j+length(CC),10}=x_fit(6*length(center)+3-3*n+kk)*x_fit(6*length(center)+3-3*n+kk+2*(j-1));
                end
                Results{length(colors_n)+1+j+length(CC),11}=x_fit(6*length(colors_n)-3*n+1+kk);
                Results{length(colors_n)+1+j+length(CC),12}=area22(cc2(j))*(0.017+0.1/diam_n(cc2(j))+Results{length(colors_n)+1+j+length(CC),11});
            end
        end
    case 3                
        switch fit_region
            case 1 % Store values for S11 region
                if isempty(complete)==1 || complete == 'n'
                    Results{1,1}='CNTs S11';
                    Results{1,2}='Diameter (nm)';
                    Results{1,3}='Center (nm) S11';
                    Results{1,4}='FWHM_l (nm) S11';
                    Results{1,5}='FWHM_g (nm) S11';
                    Results{1,6}='FWHM_v (nm) S11';
                    Results{1,7}='Broadening of FWHM_v / EPS';
                    Results{1,8}='FWHM_g / FWHM_l S11';
                    Results{1,9}='Height / Area-Frac S11';
                    Results{1,10}='Area S11';
                    Results{1,11}='Rel. Concentration (%) S11';
                    Results{1,12}='Concentration (g/l)';
                    x_len=0;
                    s11=1;
                else % Store values for entire region
                    Results{1,1}='CNTs S11';
                    Results{1,2}='Diameter (nm)';
                    Results{1,3}='Center (nm) S11';
                    Results{1,4}='FWHM_l (nm) S11';
                    Results{1,5}='FWHM_g (nm) S11';
                    Results{1,6}='FWHM_v (nm) S11';
                    Results{1,7}='Broadening of FWHM_v / EPS';
                    Results{1,8}='FWHM_g / FWHM_l S11';
                    Results{1,9}='Height / Area-Frac S11';
                    Results{1,10}='Area S11';
                    Results{1,11}='Center (nm) S22';
                    Results{1,12}='FWHM_l (nm) S22';
                    Results{1,13}='FWHM_g (nm) S22';
                    Results{1,14}='FWHM_v (nm) S22';
                    Results{1,15}='Broadening of FWHM_v / EPS';
                    Results{1,16}='FWHM_g / FWHM_l S22';
                    Results{1,17}='Height / Area-Frac S22';
                    Results{1,18}='Area S22';
                    Results{1,19}='Ratio of Height S11/S22';
                    Results{1,20}='Rel. Concentration (%) S11';
                    Results{1,21}='Rel. Concentration (%) S22';
                    Results{1,22}='Concentration (g/l)';
                    x_len=length(center_22);
                end
                Cent=center;
                FWHM=fwhm_nm_n;
                CC=cc;
                
            case 2 % Store values for S22 region
                Results{1,1}='CNTs';
                Results{1,2}='Diameter (nm)';
                Results{1,3}='Center S22 (nm)';
                Results{1,4}='FWHM_l S22 (nm)';
                Results{1,5}='FWHM_g / FWHM EPS S22 (nm)';
                Results{1,6}='FWHM_v S22 (nm)';
                Results{1,7}='Broadening of FWHM_v / EPS';
                Results{1,8}='FWHM_g / FWHM_l S22';
                Results{1,9}='Height / Area-Frac S22';
                Results{1,10}='Area S22';
                Results{1,11}='Rel. Concentration (%) S22';
                x_len=0;
                Cent=center_22;
                FWHM=fwhm_nm_n22;
                CC=cc;
        end
        
        n=0;
        for i=1:length(Cent)+length(cc)+length(cc2)
            if i<=length(Order)
                Results{i+1,2}=diam_n(i);
                Results{i+1,3}=Cent(i)+x_fit(4*i-2);
                if i==1
                    Results{i+1,4}=FWHM(i)*x_fit(4*i-1)/(0.5436+sqrt(x_fit(4*i)^2+.2166));
                    Results{i+1,5}=Results{i+1,4}*x_fit(4*i);
                    Results{i+1,8}=x_fit(4*i);
                else
                    Results{i+1,4}=FWHM(i)*x_fit(4*i-1)/(0.5436+sqrt((x_fit(4*i)*x_fit(4))^2+.2166));
                    Results{i+1,5}=Results{i+1,4}*x_fit(4*i)*x_fit(4);
                    Results{i+1,8}=x_fit(4*i)*x_fit(4);
                end                
                Results{i+1,6}=0.5436*Results{i+1,4}+sqrt(0.2166*Results{i+1,4}^2+Results{i+1,5}^2);
                Results{i+1,7}=x_fit(4*i-1);                
                Results{i+1,9}=x_fit(4*i-3);
                Results{i+1,10}=Results{i+1,9}*Results{i+1,5}/2*sqrt(pi/log(2))/real(complexErrorFunction(0,sqrt(log(2))*Results{i+1,4}/Results{i+1,5}));
                area(i)=Results{i+1,10};
                                
                if nn==2
                    if center_22(i)>0
                        Results{i+1,11}=center_22(i)+x_fit(4*(i-n)-2+4*length(Cent));
                        if i==1
                            Results{i+1,12}=fwhm_nm_n22(i)*x_fit(4*(i-n)-1+4*length(Cent))/(0.5436+sqrt(x_fit(4*(i-n)+4*length(Cent))^2+.2166));
                            Results{i+1,13}=Results{i+1,12}*x_fit(4*(length(Cent)-n+i));
                            Results{i+1,16}=x_fit(4*(length(Cent)-n+i));
                        else
                            Results{i+1,12}=fwhm_nm_n22(i)*x_fit(4*(i-n)-1+4*length(Cent))/(0.5436+sqrt((x_fit(4*(i-n)+4*length(Cent))*x_fit(4*(1-n)+4*length(Cent)))^2+.2166));
                            Results{i+1,13}=Results{i+1,12}*x_fit(4*(length(Cent)-n+i))*x_fit(4*(length(Cent)-n+1));
                            Results{i+1,16}=x_fit(4*(length(Cent)-n+i))*x_fit(4*(length(Cent)-n+1));
                        end                        
                        Results{i+1,14}=0.5436*Results{i+1,12}+sqrt(0.2166*Results{i+1,12}^2+Results{i+1,13}^2);
                        Results{i+1,15}=x_fit(4*(i-n)-1+4*length(Cent));
                        if i==1
                            Results{i+1,17}=x_fit(4*i-3)/x_fit(4*(i+length(Cent)-n)-3);
                        else
                            Results{i+1,17}=x_fit(4*i-3)/(x_fit(4*(1+length(Cent)-n)-3)*x_fit(4*(i+length(Cent)-n)-3));
                        end
                        Results{i+1,18}=Results{i+1,17}*Results{i+1,13}/2*sqrt(pi/log(2))/real(complexErrorFunction(0,sqrt(log(2))*Results{i+1,12}/Results{i+1,13}));  
                        area22(i)=Results{i+1,18};
                        Results{i+1,19}=Results{i+1,9}/Results{i+1,17};
                    else
                        n=n+1;
                    end
                end
                Results{i+1,1}=legstr{i};
            else
                Results{i+1,1}=plot_legstr{i+x_len-n};
            end
        end
        for i=1:length(Cent)
            if nn==2
                Results{i+1,20}=area(i)/sum(area)*100;
                if Results{i+1,2}<1
                    Results{i+1,22}=Results{i+1,9}/(path_length*(2210+1950*exp((0.85-Results{i+1,2})/0.158)))*12.0107;
                else
                    Results{i+1,22}=Results{i+1,9}/(path_length*(0.183+0.163*exp((0.85-Results{i+1,2})/0.158)))*1e-3;
                end
                if area22(i)>0
                    Results{i+1,21}=area22(i)/sum(area22)*100;
                end
            else
                Results{i+1,11}=area(i)/sum(area)*100;
                if s11==1
                    if Results{i+1,2}<1
                        Results{i+1,12}=Results{i+1,9}/(path_length*(2210+1950*exp((0.85-Results{i+1,2})/0.158)))*12.0107;
                    else
                        Results{i+1,12}=Results{i+1,9}/(path_length*(0.183+0.163*exp((0.85-Results{i+1,2})/0.158)))*1e-3;
                    end
                end
            end
        end
        
        kk=0;
        
        if nn==2
            for j=1:length(cc2) % Store only if there was a phonon sideband
                ll=length(x_fit)-(3+2*(length(cc2)-1));
                Results{length(colors_n)+1+j+length(CC),11}=h*c/(1e-9*(h*c/((center_22(cc2(j))+x_fit((cc2(j)*4-2)))*1e-9)+0.2+x_fit(ll+2)));
                if j==1
                    Results{length(colors_n)+1+j+length(CC),13}=fwhm_gauss*x_fit(ll+3);
                    Results{length(colors_n)+1+j+length(CC),15}=x_fit(ll+3);
                else
                    Results{length(colors_n)+1+j+length(CC),13}=fwhm_gauss*x_fit(ll+3)*x_fit(ll+3+2*(j-1));
                    Results{length(colors_n)+1+j+length(CC),15}=x_fit(ll+3)*x_fit(ll+3+2*(j-1));
                end
                
                Results{length(colors_n)+1+j+length(CC),17}=x_fit(ll+1);
                Results{length(colors_n)+1+j+length(CC),18}=Results{cc2(j)+1,9}*(0.017+0.1/diam_n(cc2(j))+x_fit(ll+1));
                
                kk=3+2*(j-1);
            end
        end
        
        for j=1:length(CC) % Store only if there was a phonon sideband
            ll=length(x_fit)-(3+2*(length(CC)-1))-kk;
            Results{length(colors_n)+1+j,3}=h*c/(1e-9*(h*c/((Cent(CC(j))+x_fit((CC(j)*4-2)))*1e-9)+0.2+x_fit(ll+2)));
            if j==1
                Results{length(colors_n)+1+j,5}=fwhm_gauss*x_fit(ll+3);
                Results{length(colors_n)+1+j,7}=x_fit(ll+3);
            else
                Results{length(colors_n)+1+j,5}=fwhm_gauss*x_fit(ll+3)*x_fit(ll+3+2*(j-1));
                Results{length(colors_n)+1+j,7}=x_fit(ll+3)*x_fit(ll+3+2*(j-1));
            end
            
            Results{length(colors_n)+1+j,9}=x_fit(ll+1);
            Results{length(colors_n)+1+j,10}=Results{CC(j)+1,9}*(0.017+0.1/diam_n(CC(j))+x_fit(ll+1));
        end
end

str=['VariablesFilmFitting',time_stamp];
current=pwd;

tt=ispc;
if tt==1
    oldD=cd([current,['\',str]]);
else
    oldD=cd([current,['/',str]]);
end

save(['Results_',time_stamp],'Results')

cd(oldD)


%..........................................................................


%--------------------------------------------------------------------------
% Section 14 - "save_data_m"
%--------------------------------------------------------------------------


function save_data_m(area_11,M11n,center_22,diam_m,x_fit,fwhm_nm_m11_n,h,c,diam_n,FWHM_22,plot_legstr,phonon_pos2,FWHM_g_22,w,folder_name,legstr_help,protocol,time_stamp)

switch w
    case {1,2}
        Results{1,1}='CNTs';
        Results{1,2}='Diameter (nm)';
        Results{1,3}='Center (nm)';
        Results{1,4}='FWHM (nm)';
        Results{1,5}='FWHM Broadening';
        Results{1,6}='Height / Area-Frac';
        Results{1,7}='Area';
        Results{1,8}='Rel. Concentration S22/M11 (%)';
        if isempty(area_11)==0
            Results{1,9}='Rel. Concentration S11 (%)';
        end
        
        for i=1:length(M11n)
            Results{i+length(center_22)+1,2}=diam_m(i);
            Results{i+length(center_22)+1,3}=M11n(i)+x_fit((i+length(center_22)-sum(center_22==0))*3-1);
            Results{i+length(center_22)+1,4}=fwhm_nm_m11_n(i)*x_fit((i+length(center_22)-sum(center_22==0))*3);
            Results{i+length(center_22)+1,5}=x_fit((i+length(center_22)-sum(center_22==0))*3);
            Results{i+length(center_22)+1,6}=x_fit((i+length(center_22)-sum(center_22==0))*3-2);
            if w==1 % Lorentzian
                Results{i+length(center_22)+1,7}=x_fit((i+length(center_22)-sum(center_22==0))*3-2)*fwhm_nm_m11_n(i)*x_fit((i+length(center_22)-sum(center_22==0))*3)*pi/2;
                area_m(i)=x_fit((i+length(center_22)-sum(center_22==0))*3-2)*fwhm_nm_m11_n(i)*x_fit((i+length(center_22)-sum(center_22==0))*3)*pi/2;
            else % Gaussian
                Results{i+length(center_22)+1,7}=x_fit((i+length(center_22)-sum(center_22==0))*3-2)*fwhm_nm_m11_n(i)*x_fit((i+length(center_22)-sum(center_22==0))*3)*sqrt(pi/(2*log(4)));
                area_m(i)=x_fit((i+length(center_22)-sum(center_22==0))*3-2)*fwhm_nm_m11_n(i)*x_fit((i+length(center_22)-sum(center_22==0))*3)*sqrt(pi/(2*log(4)));
            end
        end
        
        n=0;
        for i=1:length(center_22)+length(M11n)+length(phonon_pos2)
            if i<=length(center_22)
                if center_22(i)>0
                    Results{i+1,2}=diam_n(i);
                    Results{i+1,3}=center_22(i)+x_fit((i-n)*3-1);
                    Results{i+1,4}=FWHM_22(i)*x_fit((i-n)*3);
                    Results{i+1,5}=x_fit((i-n)*3);
                    Results{i+1,6}=x_fit((i-n)*3-2);
                    if w==1 % Lorentzian
                        Results{i+1,7}=x_fit((i-n)*3-2)*FWHM_22(i)*x_fit((i-n)*3)*pi/2;
                        area(i)=x_fit((i-n)*3-2)*FWHM_22(i)*x_fit((i-n)*3)*pi/2;
                    else % Gaussian
                        Results{i+1,7}=x_fit((i-n)*3-2)*FWHM_22(i)*x_fit((i-n)*3)*sqrt(pi/(2*log(4)));
                        area(i)=x_fit((i-n)*3-2)*FWHM_22(i)*x_fit((i-n)*3)*sqrt(pi/(2*log(4)));
                    end
                    
                    Results{i+1,1}=plot_legstr{i-n};
                else
                    n=n+1;
                    Results{i+1,1}=legstr_help{n};
                end
            else
                Results{i+1,1}=plot_legstr{i-n};
            end
        end
        
        if isempty(area_11)==1 % Just S22 was fitted
            for i=1:length(M11n)
                Results{i+length(center_22)+1,8}=area_m(i)/(sum(area_m)+sum(area))*100;
            end
        else % Entire region was fitted
            for i=1:length(M11n)
                Results{i+length(center_22)+1,8}=area_m(i)/(sum(area_m)+sum(area_11))*100;
            end
        end
        
        for i=1:length(center_22)
            if isempty(area_11)==0
                if center_22(i)>0
                    Results{i+1,8}=area(i)/sum(area)*100;
                end
                Results{i+1,9}=area_11(i)/(sum(area_m)+sum(area_11))*100;
            else
                if center_22(i)>0
                    Results{i+1,8}=area(i)/(sum(area_m)+sum(area))*100;
                end
            end
        end
        
        ll=3*(length(center_22)+length(M11n)-n);
        
        for j=1:length(phonon_pos2) % Store only if there was a phonon sideband
            Results{ll/3+n+1+j,3}=h*c/(1e-9*(h*c/((center_22(phonon_pos2(j))+x_fit((phonon_pos2(j)*3-1)))*1e-9)+0.2+x_fit(ll+2*j)));
            if j==1
                Results{ll/3+n+1+j,4}=FWHM_g_22(j)*x_fit(ll+3);
                Results{ll/3+n+1+j,5}=x_fit(ll+3);
            else
                Results{ll/3+n+1+j,4}=FWHM_g_22(j)*x_fit(ll+3)*x_fit(ll+3+2*(j-1));
                Results{ll/3+n+1+j,5}=x_fit(ll+3)*x_fit(ll+3+2*(j-1));
            end
            Results{ll/3+n+1+j,6}=x_fit(ll+1);
            Results{ll/3+n+1+j,7}=area(phonon_pos2(j))*(0.017+0.1/diam_n(phonon_pos2(j))+Results{ll/3+n+1+j,6});
        end
        
    case 3        
        Results{1,1}='CNTs';
        Results{1,2}='Diameter (nm)';
        Results{1,3}='Center (nm)';
        Results{1,4}='FWHM_l (nm)';
        Results{1,5}='FWHM_g (nm)';
        Results{1,6}='FWHM_v (nm)';
        Results{1,7}='FWHM Broadening';
        Results{1,8}='FWHM_g / FWHM_l S22';
        Results{1,9}='Height / Area-Frac';
        Results{1,10}='Area';
        Results{1,11}='Rel. Concentration S22/M11 (%)';
        if isempty(area_11)==0
            Results{1,12}='Rel. Concentration S11 (%)';
        end
        
        for i=1:length(M11n)
            Results{i+length(center_22)+1,2}=diam_m(i);
            Results{i+length(center_22)+1,3}=M11n(i)+x_fit((i+length(center_22)-sum(center_22==0))*4-2);
            if i==1
                Results{i+length(center_22)+1,4}=fwhm_nm_m11_n(i)*x_fit(4*(length(center_22)-sum(center_22==0)+i)-1)/(0.5436+sqrt(x_fit(4*(length(center_22)-sum(center_22==0)+i))^2+.2166));
                Results{i+length(center_22)+1,5}=Results{i+length(center_22)+1,4}*x_fit(4*(length(center_22)-sum(center_22==0)+i));
                Results{i+length(center_22)+1,8}=x_fit(4*(length(center_22)-sum(center_22==0)+i));
            else
                Results{i+length(center_22)+1,4}=fwhm_nm_m11_n(i)*x_fit(4*(length(center_22)-sum(center_22==0)+i)-1)/(0.5436+sqrt((x_fit(4*(length(center_22)-sum(center_22==0)+i))*x_fit(4*(length(center_22)-sum(center_22==0)+1)))^2+.2166));
                Results{i+length(center_22)+1,5}=Results{i+length(center_22)+1,4}*x_fit(4*(length(center_22)-sum(center_22==0)+i))*x_fit(4*(length(center_22)-sum(center_22==0)+1));
                Results{i+length(center_22)+1,8}=x_fit(4*(length(center_22)-sum(center_22==0)+i))*x_fit(4*(length(center_22)-sum(center_22==0)+1));
            end
            Results{i+length(center_22)+1,6}=0.5436*Results{i+length(center_22)+1,4}+sqrt(0.2166*Results{i+length(center_22)+1,4}^2+Results{i+length(center_22)+1,5}^2);
            Results{i+length(center_22)+1,7}=x_fit(4*(length(center_22)-sum(center_22==0)+i)-1);            
            Results{i+length(center_22)+1,9}=x_fit(4*(length(center_22)-sum(center_22==0)+i)-3);
            Results{i+length(center_22)+1,10}=Results{i+length(center_22)+1,9}*Results{i+length(center_22)+1,5}/2*sqrt(pi/log(2))/real(complexErrorFunction(0,sqrt(log(2))*Results{i+length(center_22)+1,4}/Results{i+length(center_22)+1,5}));
            area_m(i)=Results{i+length(center_22)+1,10};
        end
        
        n=0;
        for i=1:length(center_22)+length(M11n)+length(phonon_pos2)
            if i<=length(center_22)
                if center_22(i)>0
                    Results{i+1,2}=diam_n(i);
                    Results{i+1,3}=x_fit(4*(i-n)-2);                    
                    Results{i+1,6}=2*x_fit(4*(i-n)-1);
                    Results{i+1,4}=Results{i+1,6}/(0.5436+sqrt(x_fit(4*(i-n))^2+.2166));
                    Results{i+1,5}=Results{i+1,4}*x_fit(4*(i-n));     
                    Results{i+1,7}=Results{i+1,6}/FWHM_22(i,3);
                    Results{i+1,8}=x_fit(4*(i-n));
                    Results{i+1,9}=x_fit(4*(i-n)-3);
                    Results{i+1,10}=Results{i+1,9}*Results{i+1,5}/2*sqrt(pi/log(2))/real(complexErrorFunction(0,sqrt(log(2))*Results{i+1,4}/Results{i+1,5}));
                    area(i)=Results{i+1,10};
                    
                    Results{i+1,1}=plot_legstr{i-n};
                else
                    n=n+1;
                    Results{i+1,1}=legstr_help{n};
                end
            else
                Results{i+1,1}=plot_legstr{i-n};
            end
        end
        
        if isempty(area_11)==0
            for i=1:length(M11n)
                Results{i+length(center_22)+1,11}=area_m(i)/(sum(area_m)+sum(area_11))*100;
            end
        else
            for i=1:length(M11n)
                Results{i+length(center_22)+1,11}=area_m(i)/(sum(area_m)+sum(area))*100;
            end
        end
        
        for i=1:length(center_22)
            if isempty(area_11)==0
                if center_22(i)>0
                    Results{i+1,11}=area(i)/sum(area)*100;
                end
                Results{i+1,12}=area_11(i)/(sum(area_m)+sum(area_11))*100;
            else
                if center_22(i)>0
                    Results{i+1,11}=area(i)/(sum(area_m)+sum(area))*100;
                end
            end
        end
        
        ll=4*(length(center_22)+length(M11n)-n);
        
        for j=1:length(phonon_pos2) % Store only if there was a phonon sideband
            Results{ll/4+n+1+j,3}=h*c/(1e-9*(h*c/((x_fit(4*phonon_pos2(j)-2))*1e-9)+0.2+x_fit(ll+2*j)));
            if j==1
                Results{ll/4+n+1+j,5}=FWHM_g_22(j)*x_fit(ll+3);
                Results{ll/4+n+1+j,7}=x_fit(ll+3);
            else
                Results{ll/4+n+1+j,5}=FWHM_g_22(j)*x_fit(ll+3)*x_fit(ll+3+2*(j-1));
                Results{ll/4+n+1+j,7}=x_fit(ll+3)*x_fit(ll+3+2*(j-1));
            end
            Results{ll/4+n+1+j,9}=x_fit(ll+1);
            Results{ll/4+n+1+j,10}=area(phonon_pos2(j))*(0.017+0.1/diam_n(phonon_pos2(j))+Results{ll/4+n+1+j,9});
        end
end

oldD=cd(folder_name);

save('Result_S22andM11','Results')

cd(oldD);

if protocol=='y'
    string='%s ';
    for i=1:25
        string=[string, '%s '];
    end
    
    fid=fopen(['Protocol_',time_stamp,'.txt']);
    F = textscan(fid,string);
    fclose('all');
    
    delete(['Protocol_',time_stamp,'.txt']);
    
    bla=[F{:}];
    for i=1:size(bla,1)
        r_peaks1(i)=strcmp(strjoin(bla(i,1:5)), 'Reference file for metallic SWCNTs:');
    end
    
    idx_r_peaks1=find(r_peaks1==1);
    
    bll=bla(1:5,:);
    
    bll=[bll; bla(idx_r_peaks1(end):end,:)];
    
    for i=1:size(bll,1)
        array{i}=strjoin(bll(i,:));
    end
    
    fid=fopen('Protocol_S22andM11.txt','a');
    fprintf(fid,'%s \n',array{:});
    fclose('all');
    
    movefile('Protocol_S22andM11.txt',folder_name)
end


%..........................................................................


%--------------------------------------------------------------------------
% Section 15 - "export_data"
%--------------------------------------------------------------------------


function export_data(legstr,lambda,Absorption,y,time_stamp,folder_name)

opt=input('Do you want to export the fitted data as a .txt file (y/n)? ','s');

if opt=='y'
    
    strr=['VariablesFilmFitting',time_stamp];
    current=pwd;
    
    tt=ispc;
    
    if isempty(folder_name)==1
        if tt==1
            oldD=cd([current,['\',strr]]);
        else
            oldD=cd([current,['/',strr]]);
        end
    else
        oldD=cd(folder_name);
    end 
    
    legstr=[{'Wavelength' 'Measurement-Data' 'Calculated-Data'},legstr];
    
    if size(lambda,1)<size(lambda,2)
        lambda=lambda';
        Absorption=Absorption';
    end
    
    y=[lambda,Absorption,sum(y,2),y];
    
    str=[];
    fileID = fopen(['FittedSpectrum_',time_stamp,'.txt'],'w');
    for i=1:length(legstr)
        fprintf(fileID,['%',num2str(6*i),'s '],legstr{1,i});
        str=[str,'%',num2str(6*i),'f '];
    end
    fprintf(fileID,'\n');
    fprintf(fileID,[str,'\n'],y');
    fclose(fileID);
    
    cd(oldD)
end


