clear all; close all;
% read in data and plot
load dataTGA.mat
DTGA=TGA;
DTGA(:,2)=gradient(TGA(:,2),TGA(:,1));
DTGA(:,4)=gradient(TGA(:,4),TGA(:,3));
DTGA(:,6)=gradient(TGA(:,6),TGA(:,5));

X1=DTGA(600:end,1);
Y1=-datasmooth(DTGA(600:end,2),5,'binom')
N1=4;
Center=[203 596 628 729];
fwhm=[62 84 46 118];
x0=[Center fwhm];
UB=[Center+10 fwhm*3];
LB=[Center-10 fwhm*0.5];
[result1]=lsqcurvefit(@(x0,X1) FitGaussians(x0,X1,N1,Y1),x0,X1,Y1,LB,UB);
[FIT1,L1,A1]=FitGaussians(result1,X1,N1,Y1);

%
X2=DTGA(600:end,5);
Y2=-datasmooth(DTGA(600:end,6),3,'binom')
N2=5;
Center=[259 440 594 620 721];
fwhm=[95 181 101 52 122];
x02=[Center fwhm];
UB2=[Center+10 fwhm*3];
LB2=[Center-10 fwhm*0.5];
[result2]=lsqcurvefit(@(x02,X2) FitGaussians(x02,X2,N2,Y2),x02,X2,Y2,LB2,UB2);
[FIT2,L2,A2]=FitGaussians(result2,X2,N2,Y2);

X3=DTGA(600:end,3);
Y3=-datasmooth(DTGA(600:end,4),3,'binom')
N3=5;
Center=[246 260 342 545 685];
fwhm=[16 15 80  99 88];
x03=[Center fwhm];
UB3=[Center+10 fwhm*2];
LB3=[Center-10 fwhm*0.5];
[result3]=lsqcurvefit(@(x03,X3) FitGaussians(x03,X3,N3,Y3),x03,X3,Y3,LB3,UB3);
[FIT3,L3,A3]=FitGaussians(result3,X3,N3,Y3);

%%

% plot
figure(1);clf;set(gcf,'color','w','Position',[50 50 800 700])
subplot(3,2,1)
plot(TGA(:,1),TGA(:,2),'k','LineWidth',2); hold on;
xlabel('Temperature (°C)')
ylabel('Weight Loss (%)')
% legend('EA-SWCNTs','Box','off','Location','southwest','IconColumnWidth',10)
xlim([30 805])
ylim([0 110])
    set(gca,'TickDir','in','XMinortick','on','YMinorTick','on','TickLength',[0.03 0.05],'LineWidth',1.5,'FontSize',12)


subplot(3,2,2)
plot(X1,Y1,'k','LineWidth',2);  hold on;
for i=1:N1
plot(X1,L1(:,i).*A1(i)+L1(:,N1+1).*A1(N1+1),'Color',[0.5 0.5 0.5],'LineWidth',1.5)
end
plot(X1,FIT1,'g','LineWidth',2);
xlabel('Temperature (°C)')
ylabel('DTGA (a.u.)')
xlim([120 805])
ylim([0 0.6])
xline(result1(N1),'--')
    set(gca,'TickDir','in','XMinortick','on','YMinorTick','on','TickLength',[0.03 0.05],'LineWidth',1.5,'FontSize',12)


subplot(3,2,3)
plot(TGA(:,5),TGA(:,6),'k','LineWidth',2); hold on;
xlabel('Temperature (°C)')
ylabel('Weight Loss (%)')
% legend('TTF@SWCNTs','Box','off','Location','southwest','IconColumnWidth',10)
xlim([30 805])
ylim([0 110])
    set(gca,'TickDir','in','XMinortick','on','YMinorTick','on','TickLength',[0.03 0.05],'LineWidth',1.5,'FontSize',12)


subplot(3,2,4)
plot(X2,Y2,'k','LineWidth',2);  hold on;
for i=1:N2
    if i==1
        plot(X2,L2(:,i).*A2(i)+L2(:,N2+1).*A2(N2+1),'r','LineWidth',3)
    elseif i==2
        plot(X2,L2(:,i).*A2(i)+L2(:,N2+1).*A2(N2+1),'b','LineWidth',3)
    else
    plot(X2,L2(:,i).*A2(i)+L2(:,N2+1).*A2(N2+1),'Color',[0.5 0.5 0.5],'LineWidth',1.5)
    end
    end
plot(X2,FIT2,'g','LineWidth',2);
xlabel('Temperature (°C)')
ylabel('DTGA (a.u.)')
xlim([120 805])
ylim([0 0.4])
TOTALWEIGHT = sum(A2(1:N2))+TGA(end,6)+(100-TGA(600,6))
TTFoutside=100.*A2(1)/TOTALWEIGHT
TTFinside=100.*A2(2)/TOTALWEIGHT
text(230,0.13,strcat(num2str(round(TTFoutside,2)),'%'),'color','r')
text(390,0.13,strcat(num2str(round(TTFinside,2)),'%'),'color','b')
xline(result1(N1),'--')
    set(gca,'TickDir','in','XMinortick','on','YMinorTick','on','TickLength',[0.03 0.05],'LineWidth',1.5,'FontSize',12)





subplot(3,2,5)
plot(TGA(:,3),TGA(:,4),'k','LineWidth',2); hold on;
xlabel('Temperature (°C)')
ylabel('Weight Loss (%)')
% legend('TCNQ@SWCNTs','Box','off','Location','southwest','IconColumnWidth',10)
xlim([30 805])
ylim([0 110])
    set(gca,'TickDir','in','XMinortick','on','YMinorTick','on','TickLength',[0.03 0.05],'LineWidth',1.5,'FontSize',12)


subplot(3,2,6)
[FIT3,L3,A3]=FitGaussians(result3,X3,N3,Y3);
plot(X3,Y3,'k','LineWidth',2);  hold on;
for i=1:N3
    if ismember(i,[1 2])
        plot(X3,L3(:,i).*A3(i)+L3(:,N3+1).*A3(N3+1),'r','LineWidth',2)
    elseif i==3
        plot(X3,L3(:,i).*A3(i)+L3(:,N3+1).*A3(N3+1),'b','LineWidth',2)
    else
    plot(X3,L3(:,i).*A3(i)+L3(:,N3+1).*A3(N3+1),'Color',[0.5 0.5 0.5],'LineWidth',1.5)
    end
end
plot(X3,FIT3,'g','LineWidth',2)
xlabel('Temperature (°C)')
ylabel('DTGA (a.u.)')
xlim([120 805])
ylim([0 0.9])
xline(result1(N1),'--')


TOTALWEIGHT = sum(A3(1:N3))+TGA(end,4)+(100-TGA(600,4))

TCNQoutside=100.*(A3(1)+A3(2))/TOTALWEIGHT;
TCNQinside=100.*A3(3)/TOTALWEIGHT;
text(200,0.8,strcat(num2str(round(TCNQoutside,2)),'%'),'color','r')
text(300,0.13,strcat(num2str(round(TCNQinside,2)),'%'),'color','b')
set(gca,'TickDir','in','XMinortick','on','YMinorTick','on','TickLength',[0.03 0.05],'LineWidth',1.5,'FontSize',12)
