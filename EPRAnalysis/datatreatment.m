addpath('X:\SWCNTs\SpecialMatlabFunctions\easyspin-6.0.10\');

[X1,Y1,par1]=eprload('X:\Measurements Data\EPR\20250806\s060825c.DSC')
[X2,Y2,par2]=eprload('X:\Measurements Data\EPR\20250806\s060825l.DSC')
[X0,Y0,par0]=eprload('X:\Measurements Data\EPR\20250806\s060825a.DSC')

freq=9.44*10^9;
X0=X0.*freq/par0.MWFQ;
X1=X1.*freq/par1.MWFQ;
X2=X2.*freq/par2.MWFQ;

figure;
A=0.4;
B=0.5;
%plot(X0,Y0,'k');hold on
plot(X1,(Y1-A*Y0)*0.1,'r');hold on
plot(X2,Y2-B*Y0,'b');