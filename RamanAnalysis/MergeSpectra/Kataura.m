%% Initialize KATAURA structure
clc;
clear;


KATAURA.RBM = [];
KATAURA.WL1 = [];
KATAURA.WL2 = [];
KATAURA.WL3 = [];
KATAURA.WL4 = [];
KATAURA.D = [];
KATAURA.InverseD = [];
KATAURA.Theta = [];
KATAURA.M = [];
KATAURA.N = [];
KATAURA.Chirality = [];
KATAURA.Type = [];

% Kataura plot calculation
for m = 3:20
    for n = 0:m
        [rbm, wl1, wl2, wl3, wl4, diam, theta, type] = CalculateKataura([n, m]);
        KATAURA.RBM = [KATAURA.RBM, rbm];
        KATAURA.WL1 = [KATAURA.WL1, wl1];
        KATAURA.WL2 = [KATAURA.WL2, wl2];
        KATAURA.WL3 = [KATAURA.WL3, wl3];
        KATAURA.WL4 = [KATAURA.WL4, wl4];
        KATAURA.D = [KATAURA.D, diam];
        KATAURA.InverseD = [KATAURA.InverseD, 1/diam];
        KATAURA.Theta = [KATAURA.Theta, theta];
        KATAURA.M = [KATAURA.M, m];
        KATAURA.N = [KATAURA.N, n];
        KATAURA.Type = [KATAURA.Type, {sprintf('%s', type)}];
        KATAURA.Chirality = [KATAURA.Chirality, {sprintf('(%d,%d)', m, n)}]; % Store formatted string in cell array
    end
end


%% Plotting
colors = {[0 0.4470 0.7410], [0.4660 0.6740 0.1880], [0.6350 0.0780 0.1840], [0.9290 0.6940 0.1250]}; % Elegant colors
markers = {'v', 'diamond', 'o', 's'}; % Markers for each energy level

figure;
hold on;
title('Excitation/Emission Wavelengths vs. RBM');
xlabel('RBM frequency (cm^{-1})');
ylabel('Wavelength (nm)');
for i = 1:length(KATAURA.RBM)
    if strcmp(KATAURA.Type(i), 'M')
        color = colors{1}; % Metallic color
    else
        color = colors{3}; % Semiconducting color
    end
    
    % Loop through each energy level and plot with different marker
    for energy = 1:4
        marker = markers{energy};
        scatter(KATAURA.RBM(i), KATAURA.(['WL' num2str(energy)])(i), 50, color, marker, 'filled');
        text(KATAURA.RBM(i), KATAURA.(['WL' num2str(energy)])(i), KATAURA.Chirality{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    end
end


%% Kataura Calculation
function [nuRBM,wl1,wl2,wl3,wl4,diam,theta, type]=CalculateKataura(P)

n=P(1);
m=P(2);

diam=0.144*sqrt(3)*(sqrt(n.^2+m.^2+n.*m))/pi;
theta=atan(sqrt(3)*m./(m+2*n));
nuRBM=(223.5./diam)+12.5;

a=1.049;%eV nm
b=0.456;
c=0.812;%nm-1
    
           if mod(n-m,3)==0 % METALLIC TUBES               
%                 freq1=(1*10^7)./(150+(370*diam))+(-2000*cos(3*theta))./diam.^2;
%                 wl1=(1*10^7)./freq1; 
%                 freq2=10000000./(150+370*diam)+3000*cos(3*theta)./diam.^2;
%                 wl2=(1*10^7)./freq2;
                energy1=((a*3./diam).*(1+(b*log10(c./(3./diam)))))-0.18*cos(3*theta)/diam^2;%(in eV)
                wl1=1240/energy1;
                energy2=((a*3./diam).*(1+(b*log10(c./(3./diam)))))+0.29*cos(3*theta)./diam^2;%(in eV)
                wl2=1240/energy2;
                energy3=((a*6./diam).*(1+(b*log10(c./(6./diam)))))-0.6*cos(3*theta)/diam^2;%(in eV)
                wl3=1240/energy3;
                energy4=((a*6./diam).*(1+(b*log10(c./(6./diam)))))+0.87*cos(3*theta)./diam^2;%(in eV)
                wl4=1240/energy4;
                type = 'M';
           end
           if mod(n-m,3)==1 % semiconducting tube 
                freq1=10000000./(157.5+1066.9*diam)-710*cos(3*theta)./diam.^2;
                wl1=(1*10^7)./freq1; 
                freq2=10000000./(145.6+575.7*diam)+1375*cos(3*theta)./diam.^2;
                wl2=(1*10^7)./freq2;
                energy3=((a*4/diam).*(1+(b*log10(c./(4/diam)))))-0.42*cos(3*theta)./diam^2+(0.0596*4/diam);%(in eV)
                wl3=1240/energy3;
                energy4=((a*5/diam).*(1+(b*log10(c./(5/diam)))))+0.4*cos(3*theta)./diam^2+(0.0596*5/diam);%(in eV)
                wl4=1240/energy4;
                type = 'S';

           end
           if mod(n-m,3)==2 % semiconducting tube 
                freq1=10000000./(157.5+1066.9*diam)+369*cos(3*theta)./diam.^2;
                wl1=(1*10^7)./freq1; 
                freq2=10000000./(145.6+575.7*diam)-1475*cos(3*theta)./diam.^2;
                wl2=(1*10^7)./freq2;
                energy3=((a*4/diam).*(1+(b*log10(c./(4/diam)))))+0.42*cos(3*theta)/diam^2+(0.0596*4/diam);%(in eV)
                wl3=1240/energy3;
                energy4=((a*5/diam).*(1+(b*log10(c./(5/diam)))))-0.4*cos(3*theta)/diam^2+(0.0596*5/diam);%(in eV)
                wl4=1240/energy4;
                type = 'S';
           end
end

