%%My first script :)

global a;
a=0.144; %Interatomic spacing


% Initialize empty arrays
l11_array = [];
l22_array = [];
m_array = [];
n_array = [];
D_array = [];
InvreseD_array = [];
RBM_array = [];


%Calculate the parameters required to plot
for n = 4:15 
    for m = 1:n
        result = mod(n-m,3);
        if (result==1 || result==2)
            [l11, l22] = calculate_lambdas(n, m);
            m_array = [m_array, m];
            n_array = [n_array, n];
            l11_array = [l11_array, l11]; 
            l22_array = [l22_array, l22];
            D_array = [D_array,calculateD(n,m)];
            InvreseD_array = [InvreseD_array,1/calculateD(n,m)];
            RBM_array = [RBM_array,calculateRBM(n,m)];
            tupleArray = calculateTouple(m_array,n_array);
            stokes_shift = l11_array-l22_array;
        end
    end
end
    

figure;
scatter(InvreseD_array, RBM_array, 'MarkerFaceColor', 'g')
title('RBM vs. 1/d');
for i = 1:length(InvreseD_array)
    text(InvreseD_array(i), RBM_array(i), tupleArray(i), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end
xlabel('Inverse of d (nm^{-1})')
ylabel('RMB frequency (cm^{-1})')


figure;
scatter(RBM_array, l11_array, 'MarkerFaceColor', 'b'); hold on
scatter(RBM_array, l22_array, 'MarkerFaceColor', 'r')
title('Excitation/Emission Wavelengths vs. RBM');
for i = 1:length(RBM_array)
    text(RBM_array(i), l11_array(i), tupleArray(i), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end
for i = 1:length(RBM_array)
    text(RBM_array(i), l22_array(i), tupleArray(i), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end
xlabel('RMB frequency (cm^{-1})')
ylabel('Wavelength (nm)')
    

figure;
scatter(D_array, l11_array, 'MarkerFaceColor', 'b'); hold on
scatter(D_array, l22_array, 'MarkerFaceColor', 'r')
title('Excitation/Emission Wavelengths vs. RBM');
for i = 1:length(D_array)
    text(D_array(i), l11_array(i), tupleArray(i), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end
for i = 1:length(D_array)
    text(D_array(i), l22_array(i), tupleArray(i), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end
xlabel('Diameter (nm)')
ylabel('Wavelength (nm)')


figure;
scatter(l11_array, l22_array, 'MarkerFaceColor', 'r')
title('Excitation vs. Emission Wavelength');
for i = 1:length(l11_array)
    text(l11_array(i), l22_array(i), tupleArray(i), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end
xlabel('Emission Wavelength (nm)')
ylabel('Excitation Wavelength (nm)')


figure;
scatter(D_array, stokes_shift, 'MarkerFaceColor', 'b')
title('Stokes shift vs. Diameter');
for i = 1:length(D_array)
    text(D_array(i), stokes_shift(i), tupleArray(i), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end
xlabel('Diameter (nm)')
ylabel('Stokes shift (nm)')


function tupleArray = calculateTouple(arr1, arr2)
    
    % Preallocate cell array for tuple
    tupleArray = cell(size(arr1));
    % Convert and concatenate element-wise
    for i = 1:numel(arr1)
        tupleArray{i} = ['(',num2str(arr1(i)), ',', num2str(arr2(i)),')'];
    end
end
    
function RBM = calculateRBM(n,m)
    Acons = 223.5;
    Bcons = 12.5;
    RBM = (Acons/calculateD(n,m)) + Bcons;
end

function D = calculateD(n,m)
    global a;
    factor = n*n + n*m + m*m;
    D = (a/pi)*sqrt(3*factor);
end

function A1 = calculateA1(n,m)
    % Calculate (n-m) mod 3
    result = mod(n - m, 3);
    % Check the result
    if result == 1
        A1 = -710.0;
    elseif result == 2
        A1 = 369.0;
    end
end

function A2 = calculateA2(n,m)
    % Calculate (n-m) mod 3
    result = mod(n - m, 3);
    % Check the result
    if result == 1
        A2 = 1375.0;
    elseif result == 2
        A2 = -1475.0;
    end
end

function [lambda11, lambda22] = calculate_lambdas(n, m)
    factor = n*n + n*m + m*m;
    cos_alpha = (n+ 0.5*m)/sqrt(factor);
    D = calculateD(n,m);
    
    b=10000000.0; %
    c=157.5; %
    e=1066.9; %
    f=145.6; %
    g=575.7; %
     
    A1 = calculateA1(n,m);
    A2 = calculateA2(n,m);

    cos3factor = 4*(cos_alpha)^3-3*cos_alpha;
    
    lambda11 = b/((b/(c+e*D)) + (A1*cos3factor/(D*D)));
    lambda22 = b/((b/(f+g*D)) + (A2*cos3factor/(D*D)));
end


