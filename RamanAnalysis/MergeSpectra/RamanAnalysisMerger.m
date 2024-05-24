clc;
clear;
addpath('C:\Users\cborja\OneDrive - Universiteit Antwerpen\SWCNTs\');
%addpath('C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\SWCNTs\');
import UsefulFunctions.*;

%rootpath = 'C:\Users\Cristian Borja\OneDrive - Universiteit Antwerpen\Measurements Data\Raman\';
rootpath = 'C:\Users\cborja\OneDrive - Universiteit Antwerpen\Measurements Data\Raman\';

%All paths as default
path_20240320 = [rootpath,'20240320\'];
path_20240321 = [rootpath,'20240321\'];
path_20240514 = [rootpath,'20240514\'];
path_20240515 = [rootpath,'20240515\'];
path_20240517 = [rootpath,'20240517\'];

%Select the paths of interest

paths = {
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
        DATA_20240515.BBL570C
        DATA_20240515.BBL570D
        DATA_20240515.BBL570G

        DATA_20240514.EAL650R
        DATA_20240514.WAL650R
        }
    

for i=1:length(MyData)
    current = MyData{i};  % Access the cell array element once
    current = clip_spectrum(current, 40, 40);
    current = remove_inclination(current, 650);
    current = correct_instrument_response(current, 650);
    current = remove_bg_poly(current);
    MyData{i} = current;  % Save the result back to the cell array
end        

% plotRaman(MyData, 0.0)


%% FIGURESSSS GDBAND
% figure; clf; set(gcf,'color','w','position',[10 10 800 500])
% %TODO - Consider a new feature called DS.R which specifies the ROI
% %DS.R can be R, G, D, C
% %Define also a Zoom variable for rescaling in different areas.
% % G-band
% s=0;
% for i=[3 6]
%    current = MyData{i}
%    hold on;
%    indexes=find(or(current.X<1482,current.X>1633));
%    p=polyfit(current.X(indexes),current.Y(indexes),1);
%    current.Y=current.Y-polyval(p,current.X);
%    NORM(s+1)=max(current.Y);
%    plot(current.X,(current.Y/ NORM(s+1))-s);
%    s=s+1;
% xlabel('Raman Shift (cm^{-1})')
% ylabel('Norm. Raman')
% end
% % D-band
% s=0;
% for i=[2 5]
%    current = MyData{i}
%    hold on;
%    indexes=find(or(current.X<1200,current.X>1400));
%    p=polyfit(current.X(indexes),current.Y(indexes),1);
%    current.Y=current.Y-polyval(p,current.X);
%    plot(current.X,5*(current.Y/ NORM(s+1))-s);
%    s=s+1;
% xlabel('Raman Shift (cm^{-1})')
% ylabel('Norm. Raman')
% end
% % C-band
% s=0;
% for i=[1 4]
%    current = MyData{i}
%    hold on;
%    indexes=find(or(current.X<1633,current.X>5000));
%    p=polyfit(current.X(indexes),current.Y(indexes),1);
%    current.Y=current.Y-polyval(p,current.X);
%    plot(current.X,1*(current.Y/ NORM(s+1))-s);
%    s=s+1;
% xlabel('Raman Shift (cm^{-1})')
% ylabel('Norm. Raman')
% end
% set(gca,'FontSize',16);
% set(gca,'TickDir','in');
% set(gca,'Box','on');
% set(gca,'TickLength',[0.01 0.01]);
% set(gca,'XMinorTick','on');
% set(gca,'YMinorTick','on');
% axis([1200 2000 -2.0 1.5])
% %text(1300.5, 0.7,'D-band region (X5)')
% %text(1600.2, 0.7,'G-band region')
% xline(1593,'r--','LineWidth',0.5);
% xline(1400,'k','LineWidth',1.5);
% xline(1700,'k','LineWidth',1.5);


%% FIGURESSSS RBMs

figure; clf; set(gcf,'color','w','position',[10 10 800 500])

subplot('Position',[0.1 0.2 0.35 0.7])  
s=0;
for i=[7 8]
   current = MyData{i}
   hold on;
   indexes=find(or(current.X<230,current.X>439));
   p=polyfit(current.X(indexes),current.Y(indexes),1);
   current.Y=current.Y-polyval(p,current.X);
   plot(current.X,(current.Y/max(current.Y))+s);
    s=s-1;
xlabel('Raman Shift (cm^{-1})')
ylabel('Norm. Raman')
end

set(gca,'FontSize',16);
set(gca,'TickDir','in');set(gca,'Box','on');set(gca,'TickLength',[0.02 0.01]);
set(gca,'XMinorTick','on');set(gca,'YMinorTick','on');
axis([130 210 -2.0 1.3])

xline(175.0,'r--','LineWidth',0.5);
xline(165.0,'r--','LineWidth',0.5)
xline(160.0,'r--','LineWidth',0.5);
xline(155.0,'r--','LineWidth',0.5);
xline(182.0,'r--','LineWidth',0.5);



subplot('Position',[0.6 0.2 0.35 0.7])
hold on;
n=0:1:20;
m=0:1:20;

for j=1:19
    for i=1:19
        s=i+(j-1)*19;
        if ge(n(i),m(j))
        P(:,s)=[n(i) m(j)];
        else 
        P(:,s)=[0 0];
        end
    end
end
%enkel NTs met (n,m) n>=m erin houden
k=find(ne(P(1,:),0));
PP=P(:,k);
N=size(k);
% alles berekenen voor deze NTs
for i=1:N(2);
    [nuRBM(i),wl1(i),wl2(i),wl3(i),wl4(i),diam(i),theta(i),type(i)] = CalculateKataura(PP(:,i))
end

for i=1:N(2)
        if (and(nuRBM(i) < 200,nuRBM(i)>100))
           plot(nuRBM(i),wl1(i),'o','color','r')
           text(nuRBM(i),wl1(i),strcat('(',num2str(PP(1,i)),',',num2str(PP(2,i)),')'),'color','r')
           
           plot(nuRBM(i),wl2(i),'o','color','g')
           text(nuRBM(i),wl2(i),strcat('(',num2str(PP(1,i)),',',num2str(PP(2,i)),')'),'color','g')
           
           plot(nuRBM(i),wl3(i),'o','color','b')
           text(nuRBM(i),wl3(i),strcat('(',num2str(PP(1,i)),',',num2str(PP(2,i)),')'),'color','b')
           
           plot(nuRBM(i),wl4(i),'o','color',[0 0.5 1])
           text(nuRBM(i),wl4(i),strcat('(',num2str(PP(1,i)),',',num2str(PP(2,i)),')'),'color',[0 0.5 1])
        end
end


axis([100 220 550 750])
set(gca,'FontSize',16);
set(gca,'TickDir','in');set(gca,'Box','on');set(gca,'TickLength',[0.02 0.01]);
set(gca,'XMinorTick','on');set(gca,'YMinorTick','on');
xlabel('Raman Shift (cm^{-1})','FontSize',16);
ylabel('Excitation Wavelength (nm)','FontSize',16);

yline(650,'m--','LineWidth',1.5);
xline(175.0,'r--','LineWidth',0.5);
xline(165.0,'r--','LineWidth',0.5)
xline(160.0,'r--','LineWidth',0.5);
xline(155.0,'r--','LineWidth',0.5);
xline(182.0,'r--','LineWidth',0.5);
