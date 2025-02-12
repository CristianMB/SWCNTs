function [X,Y,P] = clip_spectrum(X0,Y0,P0,ClipLeft0,ClipRight0)
%function clips spectrum at the left and right sides by calculating its derivative; 
%ClipLeft and ClipRight determine number of pixels to remove from the center of left and right derivative peaks,respectively; 
%By default, ClipLeft = 8, ClipRight = 9;

load InclinationCoeff.mat IncPoly3;

LeftLim = 130; %all peaks of the spectrum's derivative curve should be to the left of LeftLim - determined visually from graph with all derivative curves
RightLim = 900; %all peaks are to the right of RightLim

ClipRangeL = ClipLeft0; %number of pixels to add/substract to the deriv. peak
ClipRangeR = ClipRight0;

arb = 3; %tolerance = PeakHeight/AverageBaselineValue; determines how high should be the derivative peak, in respect to average value, to be treated as real; 
ClipDef = 150; %default number of pixels to clip if no clear peak present in the derivative

f.x = X0;
f.y = Y0;
f.p = P0;
     
deriv = diff(f.y); %finding derivative for current spectrum
deriv = [deriv;nan];
    
cond = (f.p >= LeftLim)&(f.p<=RightLim);
mdn = median(deriv(cond), 'omitnan'); %calculating median of derivative
    
condL = (f.p <= LeftLim);
dleft = max(deriv(condL)) - mdn; %calculating difference between peak and median
    
condR = (f.p>=RightLim);
dright = mdn - min(deriv(condR));
dev = mad(deriv,1); %calculating average deviation from the median
    
rl = dleft/dev; %calculating how different (peak-median) from average deviation from median 
rr = dright/dev;
    
if rl>=arb
   pleft=f.p(deriv == max(deriv(condL)));
   if length(pleft)>=2 %to make sure there is no  such value in the right part condR (happened once)
       pleft = min(pleft);
   end
   ClipLeft = pleft + ClipRangeL;        
else
   ClipLeft = 1+ClipDef;
end    
    
if rr>=arb
   pright=f.p(deriv == min(deriv(condR)));
   if length(pright)>=2
       pright = max(pright);
   end
   ClipRight = pright - ClipRangeR;
else
   ClipRight = 1024-ClipDef;
end
    
%here we determine the convergence range (in terms of pixels) of the cubic polynomial IncPoly3 interpolated for the x_mean of the spectrum;
PixFirst = interp1(IncPoly3(:,9),IncPoly3(:,7),f.x(f.p==512));
PixLast = interp1(IncPoly3(:,9),IncPoly3(:,6),f.x(f.p==512));
    
if ClipLeft<PixFirst %choose the highest starting pixel for clipping
   ClipLeft = PixFirst;
end 
if ClipRight >PixLast %choose the smallest ending pixel for clipping
   ClipRight = PixLast;
 end
    
condition = (f.p >= ClipLeft)&(f.p <= ClipRight);
X = f.x(condition);
Y = f.y(condition);
P = f.p(condition);
    
