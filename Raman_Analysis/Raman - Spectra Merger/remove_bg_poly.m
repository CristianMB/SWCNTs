function [X,Y] = remove_bg_poly(X0,Y0)
%removes background according to the polynomial method (Zhao et al. Applied Spectroscopy 61 11 2007 - 10.1366/000370207782597003)

spectrum0 = [X0,Y0];
EndCycle = false; %declaring logical variable (= criterion for ending the cycle) to start the first iteration of the while loop
i=1; %iteration number

fitfunc = 'poly1'; %order of the polynomial function
options = fitoptions; % initializing "fitoptions" object for removing points during the fit after the first iteration
% options = fitoptions(fitfunc, 'Normalize','On','Robust', 'LAR'); %to activate if normalization is needed; however, it works worse than without normalization

while not(EndCycle)
   
   %fitting spectrum with high-order polynomial function
   if i == 1 %choose which spectrum to fit for current iteration
       current_spectrum = spectrum0(:,2); %take intensity array of input spectrum
   else
       current_spectrum = spectrum{i-1};
   end
   
   [cpoly{i},gofs{i}] = fit(spectrum0(:,1), current_spectrum, fitfunc, options);
   
%    if options.Normalize == 'on' %code to use when normalization is on for fit function
%        meanx = mean(spectrum0(:,1));
%        stdx = std(spectrum0(:,1));
%        tempx = (spectrum0(:,1)-meanx)./stdx;
%    else
%        tempx = spectrum0(:,1);
%    end
   polyfit{i} = polyval(coeffvalues(cpoly{i}),spectrum0(:,1)); %estimating polynomial function at every x value
   
   %calculating residual
   res{i} = current_spectrum - polyfit{i}; 
   
   %calculating standard deviation of residual
   dev{i} = std(res{i},1);
   sum = polyfit{i}+dev{i};
   
   %reconstructing model input for the next iteration
   spectrum{i} = sum;
   if i == 1    %remove peaks from fitting for first iteration
      PointsToExclude = current_spectrum > sum;
      options.Exclude = excludedata(spectrum0(:,1),current_spectrum,'indices',PointsToExclude);
   end
   
   cond_mod = current_spectrum < sum; %determine at which frequencies spectrum is less then polynomial function
   spectrum{i}(cond_mod) = current_spectrum(cond_mod);
   
   
   if i~= 1
       EndCycle = abs((dev{i}-dev{i-1})/dev{i})<0.05; %calculate condition to end cycle
   end
   i = i+1; %change to next iteration
end

X = spectrum0(:,1);
Y = spectrum0(:,2) - polyfit{i-1}; % spectral intensity without background contribution

