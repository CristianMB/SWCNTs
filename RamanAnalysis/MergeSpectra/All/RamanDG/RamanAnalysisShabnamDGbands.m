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
[3 4 5 6 7 8 9 12 13 14 15]
% G-band
s=0;
for i=[3 7 8 12 14]
   hold on;
   QQ=find(or(RamanShift(:,i)<1482,RamanShift(:,i)>1633));
   p=polyfit(RamanShift(QQ,i),DATAb(QQ,i),1);
   DATAb(:,i)=DATAb(:,i)-polyval(p,RamanShift(:,i));
   NORM(s+1)=max(DATAb(:,i));
   plot(RamanShift(:,i),(DATAb(:,i)/ NORM(s+1))-s);
   s=s+1;
xlabel('Raman Shift (cm^{-1})')
ylabel('Norm. Raman')
end
% D-band
s=0;
for i=[4 6 9 13 15]
   hold on;
   QQ=find(or(RamanShift(:,i)<1250,RamanShift(:,i)>1404));
   p=polyfit(RamanShift(QQ,i),DATAb(QQ,i),1);
   DATAb(:,i)=DATAb(:,i)-polyval(p,RamanShift(:,i));
   plot(RamanShift(:,i),5*(DATAb(:,i)/ NORM(s+1))-s);
   s=s+1;
xlabel('Raman Shift (cm^{-1})')
ylabel('Norm. Raman')
end

set(gca,'FontSize',16);
set(gca,'TickDir','in');set(gca,'Box','on');set(gca,'TickLength',[0.02 0.01]);
set(gca,'XMinorTick','on');set(gca,'YMinorTick','on');
axis([1220 1680 -4.2 1.1])
text(1300.5, 0.7,'D-band region (X5)')
text(1600.2, 0.7,'G-band region')
xline(1589.5,'r--','LineWidth',0.5);
xline(1458,'k','LineWidth',5);