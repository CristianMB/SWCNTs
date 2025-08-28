clc;
clear;
addpath('X:\SWCNTs');

addpath('X:\Measurements Data\Raman');
%addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\');
import UsefulFunctions.*;

%rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements Data\Raman\';
rootpath = 'X:\Measurements Data\Raman';

%All paths as default
path_20240320 = [rootpath,'\20240320\'];
path_20240321 = [rootpath,'\20240321\'];
path_20240514 = [rootpath,'\20240514\'];
path_20240515 = [rootpath,'\20240515\'];
path_20240517 = [rootpath,'\20240517\'];

%Select the paths of interest

paths = {
    path_20240514
    path_20240515
    };

%All paths as default
path_20240426 = [rootpath,'\20240426\'];
path_20240514 = [rootpath,'\20240514\'];
path_20240515 = [rootpath,'\20240515\'];

%Select the paths of interest

paths = {
    path_20240426
    path_20240514
    path_20240515
    };


ReadRamanFromPaths(paths);

%%% --------------------------------------------------------
%%%           Correcting for sensitivity detector and inclination and
%%%           cutting edges away
%%% --------------------------------------------------------

% DATA_20240514.S2L650D = clip_spectrum(DATA_20240514.S2L650D,10,10);
% DATA_20240514.S2L650D = remove_inclination(DATA_20240514.S2L650D,650);
% DATA_20240514.S2L650D = correct_instrument_response(DATA_20240514.S2L650D, 650);
% DATA_20240514.S2L650D = remove_bg_poly(DATA_20240514.S2L650D);

MyData={
            DATA_20240515.BAL570C
            DATA_20240515.BAL570D
            DATA_20240515.BAL570G
            DATA_20240515.BAL570RA
            DATA_20240515.BAL570RB
            DATA_20240515.BBL570C
            DATA_20240515.BBL570D
            DATA_20240515.BBL570G
            DATA_20240515.BBL570RA
            DATA_20240515.BBL570RB

            DATA_20240426.BAL650D,
            DATA_20240426.BAL650G,
            DATA_20240426.BAL650R,
            DATA_20240426.BAL650RB,
            DATA_20240426.BBL650G,
            DATA_20240426.BBL650D,
            DATA_20240426.BBL650RA,
            DATA_20240426.BBL650RB
            DATA_20240514.BAL650C1
            DATA_20240514.BAL650C2
            DATA_20240514.BBL650C1
            DATA_20240514.BBL650C2  
        };
    
WL = 650;

for i=1:length(MyData)
    current = MyData{i};  % Access the cell array element once
    current = clip_spectrum(current, 40, 40);
    current = remove_inclination(current, WL);
    current = correct_instrument_response(current, WL);
    current = remove_bg_poly(current);
    MyData{i} = current;  % Save the result back to the cell array
end        

% plotRaman(MyData, 0.0)


%% FIGURESSSS GDBAND
%TODO - Consider a new feature called DS.R which specifies the ROI
%DS.R can be R, G, D, C
%Define also a Zoom variable for rescaling in different areas.
% G-band

% yex_min = 550;
% yex_max = 750;
% yra_min = -1.5;
% yra_max =  1.1;
% x_min = 130;
% x_max = 210;
% wl = 650;

GBand = [3 8 12 15];
GFactor = 1.0;
DBand = [2 7 11 16];
DFactor = 1.0;
CBand = [1 6 19 21];
CFactor = 1.0;

Peaks = [155.0, 160.0, 165.0, 175.0, 182.0];

figure; clf; set(gcf,'color','w','position',[10 10 800 500])
s=0;

for i=GBand
   current = MyData{i};
   hold on;
   indexes=find(or(current.X<1482,current.X>1633));
   p=polyfit(current.X(indexes),current.Y(indexes),1);
   current.Y=current.Y-polyval(p,current.X);
   NORM(s+1)=max(current.Y);
   plot(current.X,GFactor*(current.Y/ NORM(s+1))-s)
   s=s+1;
end

% D-band
s=0;
for i=DBand
   current = MyData{i};
   hold on;
   indexes=find(or(current.X<1200,current.X>1400));
   p=polyfit(current.X(indexes),current.Y(indexes),1);
   current.Y=current.Y-polyval(p,current.X);
   plot(current.X,DFactor*(current.Y/ NORM(s+1))-s);
   s=s+1;
end

% C-band
s=0;
for i=CBand
   current = MyData{i};
   hold on;
   indexes=find(or(current.X<1633,current.X>5000));
   p=polyfit(current.X(indexes),current.Y(indexes),1);
   current.Y=current.Y-polyval(p,current.X);
   plot(current.X,CFactor*(current.Y/ NORM(s+1))-s);
   s=s+1;

end

set(gca,'FontSize',16);
set(gca,'TickDir','in');
set(gca,'Box','on');
set(gca,'TickLength',[0.01 0.01]);
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
axis([1300 1900 -3.5 1.5])

%text(1300.5, 0.7,'D-band region (X5)')
%text(1600.2, 0.7,'G-band region')
% xline(1593,'r--','LineWidth',0.5);
% xline(1400,'k','LineWidth',1.5);
% xline(1700,'k','LineWidth',1.5);

xlabel('Raman Shift (cm^{-1})')
ylabel('Norm. Raman')

%% FIGURESSSS RBMs
% 
% yex_min = 550;
% yex_max = 750;
% yra_min = -1.5;
% yra_max =  1.1;
% x_min = 130;
% x_max = 210;
% wl = 650;
% RamanSpectra = [7 8];
% Peaks = [155.0, 160.0, 165.0, 175.0, 182.0];
% 
% figure; clf; set(gcf,'color','w','position',[10 10 800 500])
% subplot('Position',[0.1 0.2 0.35 0.7])  
% s=0;
% for i= RamanSpectra
%    current = MyData{i};
%    hold on;
%    indexes=find(or(current.X<230,current.X>439));
%    p=polyfit(current.X(indexes),current.Y(indexes),1);
%    current.Y=current.Y-polyval(p,current.X);
%    plot(current.X,(current.Y/max(current.Y))+s);
%     s=s-1;
% end
% set(gca,'FontSize',16);
% set(gca,'TickDir','in');set(gca,'Box','on');set(gca,'TickLength',[0.02 0.01]);
% set(gca,'XMinorTick','on');set(gca,'YMinorTick','on');
% xlabel('Raman Shift (cm^{-1})')
% ylabel('Norm. Raman')
% axis([x_min x_max yra_min yra_max])
% for i=Peaks
%     xline(i,'r--','LineWidth',0.5);
% end
% subplot('Position',[0.6 0.2 0.35 0.7])
% hold on;
% n=0:1:20;
% m=0:1:20;
% for j=1:19
%     for i=1:19
%         s=i+(j-1)*19;
%         if ge(n(i),m(j))
%         P(:,s)=[n(i) m(j)];
%         else 
%         P(:,s)=[0 0];
%         end
%     end
% end
% k=find(ne(P(1,:),0));
% PP=P(:,k);
% N=size(k);
% for i=1:N(2);
%     [nuRBM(i),wl1(i),wl2(i),wl3(i),wl4(i),diam(i),theta(i),type(i)] = CalculateKataura(PP(:,i));
% end
% for i=1:N(2)
%         if (strcmp(type(i), 'S') == 1)
%            marker = 'o';
%         else
%            marker = 'v';
%         end
%         
%         if (nuRBM(i) < 200 && nuRBM(i) > 140) && ...
%            ((wl1(i) >= yex_min && wl1(i) <= yex_max) || ...
%             (wl2(i) >= yex_min && wl2(i) <= yex_max) || ...
%             (wl3(i) >= yex_min && wl3(i) <= yex_max) || ...
%             (wl4(i) >= yex_min && wl4(i) <= yex_max))
%     
%            scatter(nuRBM(i),wl1(i),marker,'r');
%            text(nuRBM(i),wl1(i),strcat('(',num2str(PP(1,i)),',',num2str(PP(2,i)),')'),'color','r');
%            
%            scatter(nuRBM(i),wl2(i),marker,'m');
%            text(nuRBM(i),wl2(i),strcat('(',num2str(PP(1,i)),',',num2str(PP(2,i)),')'),'color','m');
%            
%            scatter(nuRBM(i),wl3(i),marker,'g');
%            text(nuRBM(i),wl3(i),strcat('(',num2str(PP(1,i)),',',num2str(PP(2,i)),')'),'color','g');
%            
%            scatter(nuRBM(i),wl4(i),marker,'b');
%            text(nuRBM(i),wl4(i),strcat('(',num2str(PP(1,i)),',',num2str(PP(2,i)),')'),'color','b')
%         end
% end
% axis([x_min x_max yex_min yex_max])
% set(gca,'FontSize',16);
% set(gca,'TickDir','in');set(gca,'Box','on');set(gca,'TickLength',[0.02 0.01]);
% set(gca,'XMinorTick','on');set(gca,'YMinorTick','on');
% xlabel('Raman Shift (cm^{-1})','FontSize',16);
% ylabel('Excitation Wavelength (nm)','FontSize',16);
% yline(wl,'k--','LineWidth',2.0);
% for i=Peaks
%     xline(i,'r--','LineWidth',0.5);
% end
