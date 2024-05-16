% Kataura plot

[nuRBM,wl1,wl2,wl3,wl4,diam,theta]=KatauraF([7,5])

function [nuRBM,wl1,wl2,wl3,wl4,diam,theta]=KatauraF(P)

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
           end
end
