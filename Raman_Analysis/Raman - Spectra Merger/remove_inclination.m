function Y = remove_inclination(X0,Y0,P0,laser)
%Function removes inclination effect for the spectrum for Dilor XY spectrometer

load InclinationCoeff.mat IncPoly3;

XL = 10^7./(10^7./laser - X0); % change x from cm-1 to nm

center = XL(P0==512);

p(1) = interp1(IncPoly3(:,8), IncPoly3(:,2),center,'spline');
p(2) = interp1(IncPoly3(:,8), IncPoly3(:,3),center,'spline');
p(3) = interp1(IncPoly3(:,8), IncPoly3(:,4),center,'spline');
p(4) = interp1(IncPoly3(:,8), IncPoly3(:,5),center,'spline');

y= @(p,x) p(1)*x.^3+p(2)*x.^2+p(3)*x+p(4);

Y = Y0./y(p,XL);




