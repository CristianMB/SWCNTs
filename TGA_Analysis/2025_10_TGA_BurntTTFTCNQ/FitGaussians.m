function [FIT,L,A]=FitGaussians(x0,X,N,ExpY)


L=zeros(length(X),N+1);
A = zeros(2*N+2,1);
FIT=zeros(size(ExpY));
numcols=size(ExpY,2);
for k=1:numcols
%We put the lorentzians in a matrix
    for i=1:N
        L(:,i)=gaussian(X,x0(i),x0(N+i)); % the empty one
    end
        L(:,i+1)=ones(length(X),1);

                 A = mldivide(L,ExpY);
                 negPeaks = A(1:N)<0;
                  while any(negPeaks) % repeat until no more negative amplitudes
                    % Negative amplitudes were found. Remove the corresponding
                    % basis functions from L and recalculate.
                    posPeaks = A(1:N) > 0; 
                    posPeaks = [posPeaks; true];
                    A(posPeaks) = mldivide(L(:,posPeaks),ExpY);
                    A(negPeaks) = 0;
                    % Check again for negative amplitudes.
                    negPeaks = A(1:N)<0;
                  end

    FIT=L*A;

end
end
    
    

