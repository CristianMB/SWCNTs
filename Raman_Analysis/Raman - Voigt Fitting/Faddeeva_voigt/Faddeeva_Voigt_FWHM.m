function F=Faddeeva_Voigt_FWHM(X,x0,FWHM,A)
% FWHM=0.5346FL+sqrt(0.216598FL^2+FG^2) Voigt Profile
% A=FL/FWHM

if A == 0 % then FL = 0 and thus FWHM is in fact FG (gauss)
    sigma = FWHM/(2*sqrt(2*log(2)));
    F=(1/(sigma*sqrt(2*pi)))*exp(-(X-x0).^2./(2*sigma^2));
elseif A == 1 % then FG=0 and thus FWHM is FL
    gamma=FWHM/2;
    F=(gamma/pi)./(((X-x0).^2)+gamma^2);
elseif and(gt(A,0),lt(A,1))
    FL=A*FWHM;
    gamma=FL/2;
    FG=sqrt((FWHM-0.5346*FL)^2-0.216598*FL^2);
    sigma=FG/(2*sqrt(2*log(2)));
    F=Faddeeva_Voigt(X,x0,sigma,gamma);
else % not valid A must be between 0 and 1
    error('A value not valid, should be between 0 and 1');
    F=zeros(size(X));
end
 
end
