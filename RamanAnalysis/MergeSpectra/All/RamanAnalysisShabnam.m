clear all;
close all;
NAMES={'a','b','c','d','e','f','g','h','i','m','n','o','p','q','rb','s'};

for i=1:16
test=RdExp(strcat('S25424',NAMES{i},'.m3d'));
RamanShift(:,i)=test(:,1);
    for j=1:1024
    Q=sort(test(j,3:12));
    DATA(j,i)=mean(Q(2:9));
    end
end



    %%% --------------------------------------------------------
    %%%           Correcting for sensitivity detector and inclination and
    %%%           cutting edges away
    %%% --------------------------------------------------------
    
    % clip edges of spectra - 
    ClipLeft = 75; %how many pixels to clip from the edges
    ClipRight = 25;

    RamanShift(1:ClipLeft,:)=[];RamanShift(1024-ClipLeft-ClipRight:1024-ClipLeft,:)=[];
    DATA(1:ClipLeft,:)=[];DATA(1024-ClipLeft-ClipRight:1024-ClipLeft,:)=[];

    % remove spectrum inclination
    for k = 1:16
        DATAb(:,k) = remove_inclination(RamanShift(:,k),DATA(:,k),[1:length(DATA)],650);
    end
   
    % correct for the Dilor XY instrument response
    for k = 1:16
        xl = 10^7./(10^7./650 - RamanShift(:,k)); % x values in nm
        DATAb(:,k)= correct_instrument_response(xl,DATAb(:,k));    
    end

%%
figure; clf; set(gcf,'color','w','position',[10 10 800 500])

subplot('Position',[0.1 0.2 0.35 0.7])  
s=0;
for i=[1 2 10 11 16]
   hold on;
   QQ=find(or(RamanShift(:,i)<230,RamanShift(:,i)>439));
   p=polyfit(RamanShift(QQ,i),DATAb(QQ,i),1);
   DATAb(:,i)=DATAb(:,i)-polyval(p,RamanShift(:,i));
    plot(RamanShift(:,i),(DATAb(:,i)/max(DATAb(:,i)))+s);
    s=s-1;
xlabel('Raman Shift (cm^{-1})')
ylabel('Norm. Raman')
end
set(gca,'FontSize',16);
set(gca,'TickDir','in');set(gca,'Box','on');set(gca,'TickLength',[0.02 0.01]);
set(gca,'XMinorTick','on');set(gca,'YMinorTick','on');
axis([210 350 -4.2 1.1])

xline(283.125,'r--','LineWidth',0.5);xline(285.375,'r--','LineWidth',0.5)
xline(298.125,'r--','LineWidth',0.5);xline(300.375,'r--','LineWidth',0.5);
xline(264.375,'r--','LineWidth',0.5);

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
    [nuRBM(i),wl1(i),wl2(i),wl3(i),wl4(i),diam(i),theta(i)]=Kataura(PP(:,i));
end


for i=1:N(2)
    if mod(PP(1,i)-PP(2,i),3)==0
%    plot(wl1(i),diam(i),'o','color',[1 0 0])
%    plot(wl2(i),diam(i),'o','color',[1 0 0])
    else
   if and(and(nuRBM(i) < 350,nuRBM(i)>20),and(wl2(i) < 750,wl2(i)>500))
   plot(nuRBM(i),wl2(i),'o','color',[0 0.5 0])
   text(nuRBM(i),wl2(i),strcat('(',num2str(PP(1,i)),',',num2str(PP(2,i)),')'),'color','b')
   end
   %    plot(wl3(i),diam(i),'o','color',[0 0.5 0])
%    plot(wl4(i),diam(i),'o','color',[1 0 1])
    end
end
axis([200 350 500 750])
set(gca,'FontSize',16);
set(gca,'TickDir','in');set(gca,'Box','on');set(gca,'TickLength',[0.02 0.01]);
set(gca,'XMinorTick','on');set(gca,'YMinorTick','on');
xlabel('Raman Shift (cm^{-1})','FontSize',16);
ylabel('Excitation Wavelength (nm)','FontSize',16);

yline(650,'m--','LineWidth',1.5);

xline(283.125,'r--','LineWidth',0.5);xline(285.375,'r--','LineWidth',0.5)
xline(298.125,'r--','LineWidth',0.5);xline(300.375,'r--','LineWidth',0.5);
xline(264.375,'r--','LineWidth',0.5);