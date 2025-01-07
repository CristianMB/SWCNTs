function [FIT,L,A] = RamanFit(xdata,ydata,x0,NumSpec,NumLor)

NumPoints=size(xdata,1);
L=zeros(NumPoints,NumLor+2);

%Lorentzian Matrix
    for k=1:NumLor
        Gamma=x0(2*k);
        Pos=x0(2*k-1);
        L(:,k)=(1/pi).*(Gamma/2)./((xdata-Pos).^2+(Gamma/2).^2); 
    end

%We create two extra columns, one for a constant background and another one
%for a slope
        L(:,NumLor+1)=ones(NumPoints,1);  
        L(:,NumLor+2)=xdata;

FIT = zeros(size(ydata));
A = zeros(NumLor+2,NumSpec);

% We calculate analitically the amplitudes at each loop of the optimisation
% simultaneous fit means each time we use the same basis functions but we
% calculate for each spectrum different amplitudes

for i=1:NumSpec
        A(:,i) = mldivide(L,ydata(:,i));
        % check which amplitudes of the Lorentzians are negative 
            negPeaks = A(1:NumLor,i) < 0;
            while any(negPeaks) % repeat until no more negative amplitudes
                    % Remove the corresponding
                    % basis functions from L and recalculate.
                    posPeaks = A(1:NumLor,i) > 0;
                    % only take the posPeaks and add the linear baseline
                    % ones
                    A([posPeaks; logical(1); logical(1)],i) = mldivide(L(:,[posPeaks; logical(1); logical(1)]),ydata(:,i));
                    A(negPeaks,i) = 0;
                    % Check again for negative amplitudes.
                    negPeaks = A(1:NumLor,i)<0;
            end
        
%We create the fit function
    FIT(:,i)=L*A(:,i);
end
    
end
