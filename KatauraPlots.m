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
mod_array = [];


%Calculate the parameters required to plot
for m = 1:16
    for n = 0:m
        result = mod(n-m,3);
        if (result==1 || result==2)
            [l11, l22] = calculate_lambdas(n, m);
            mod_array = [mod_array, result];
            m_array = [m_array, m];
            n_array = [n_array, n];
            l11_array = [l11_array, l11]; 
            l22_array = [l22_array, l22];
            D_array = [D_array,calculateD(n,m)];
            InvreseD_array = [InvreseD_array,1/calculateD(n,m)];
            RBM_array = [RBM_array,calculateRBM(n,m)];
            tupleArray = calculateTouple(m_array,n_array);
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
scatter(D_array, RBM_array, 'MarkerFaceColor', 'g')
title('RBM Frequency vs. Diameter');
xlabel('Carbon Nanotube Diameter (nm)')
ylabel('RMB Frequency (cm^{-1})')
for i = 1:length(D_array)
    text(D_array(i), RBM_array(i), tupleArray(i));
end



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
scatter(l11_array, l22_array, [], mod_array, 'filled')
title('Excitation vs. Emission Wavelength');
for i = 1:length(l11_array)
    text(l11_array(i), l22_array(i), tupleArray(i), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end
colorbar
xlabel('Emission Wavelength (nm)')
ylabel('Excitation Wavelength (nm)')





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




clc;
clear;

% Initialize KATAURA structure
KATAURA.RBM = [];
KATAURA.WL1 = [];
KATAURA.WL2 = [];
KATAURA.WL3 = [];
KATAURA.WL4 = [];
KATAURA.D = [];
KATAURA.InverseD = [];
KATAURA.Theta = [];
KATAURA.M = [];
KATAURA.N = [];
KATAURA.Chirality = [];
KATAURA.Type = [];

% Kataura plot calculation
for m = 1:16
    for n = 0:m
        [rbm, wl1, wl2, wl3, wl4, diam, theta, type] = CalculateKataura([n, m]);
        KATAURA.RBM = [KATAURA.RBM, rbm];
        KATAURA.WL1 = [KATAURA.WL1, wl1];
        KATAURA.WL2 = [KATAURA.WL2, wl2];
        KATAURA.WL3 = [KATAURA.WL3, wl3];
        KATAURA.WL4 = [KATAURA.WL4, wl4];
        KATAURA.D = [KATAURA.D, diam];
        KATAURA.InverseD = [KATAURA.InverseD, 1/diam];
        KATAURA.Theta = [KATAURA.Theta, theta];
        KATAURA.M = [KATAURA.M, m];
        KATAURA.N = [KATAURA.N, n];
        KATAURA.Chirality = [KATAURA.Chirality, {sprintf('(%d,%d)', m, n)}];
        KATAURA.Type = [KATAURA.Type, type];
    end
end

% Plotting
figure;
hold on;
for i = 1:length(KATAURA.RBM)
    if strcmp(KATAURA.Type{i}, 'Metallic')
        marker = 'o'; % Circle for metallic
    else
        marker = 's'; % Square for semiconducting
    end
    scatter(KATAURA.RBM(i), KATAURA.WL1(i), marker, 'MarkerFaceColor', 'b');
    scatter(KATAURA.RBM(i), KATAURA.WL2(i), marker, 'MarkerFaceColor', 'g');
    scatter(KATAURA.RBM(i), KATAURA.WL3(i), marker, 'MarkerFaceColor', 'r');
    scatter(KATAURA.RBM(i), KATAURA.WL4(i), marker, 'MarkerFaceColor', 'y');
    text(KATAURA.RBM(i), KATAURA.WL1(i), KATAURA.Chirality{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    text(KATAURA.RBM(i), KATAURA.WL2(i), KATAURA.Chirality{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    text(KATAURA.RBM(i), KATAURA.WL3(i), KATAURA.Chirality{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end
title('Excitation/Emission Wavelengths vs. RBM');
xlabel('RBM frequency (cm^{-1})');
ylabel('Wavelength (nm)');
hold off;

function [nuRBM, wl1, wl2, wl3, wl4, diam, theta, type] = CalculateKataura(P)
    n = P(1);
    m = P(2);

    diam = 0.144 * sqrt(3) * (sqrt(n^2 + m^2 + n * m)) / pi;
    theta = atan(sqrt(3) * m / (m + 2 * n));
    nuRBM = (223.5 / diam) + 12.5;

    a = 1.049; % e


    
    
    clc;
clear;

% Initialize KATAURA structure
KATAURA.RBM = [];
KATAURA.WL1 = [];
KATAURA.WL2 = [];
KATAURA.WL3 = [];
KATAURA.WL4 = [];
KATAURA.D = [];
KATAURA.InverseD = [];
KATAURA.Theta = [];
KATAURA.M = [];
KATAURA.N = [];
KATAURA.Chirality = [];
KATAURA.Type = []; % New field to store the type of nanotube

% Kataura plot calculation
for m = 1:16
    for n = 0:m
        [rbm, wl1, wl2, wl3, wl4, diam, theta, type] = CalculateKataura([n, m]);
        KATAURA.RBM = [KATAURA.RBM, rbm];
        KATAURA.WL1 = [KATAURA.WL1, wl1];
        KATAURA.WL2 = [KATAURA.WL2, wl2];
        KATAURA.WL3 = [KATAURA.WL3, wl3];
        KATAURA.WL4 = [KATAURA.WL4, wl4];
        KATAURA.D = [KATAURA.D, diam];
        KATAURA.InverseD = [KATAURA.InverseD, 1/diam];
        KATAURA.Theta = [KATAURA.Theta, theta];
        KATAURA.M = [KATAURA.M, m];
        KATAURA.N = [KATAURA.N, n];
        KATAURA.Chirality = [KATAURA.Chirality, {sprintf('(%d,%d)', m, n)}];
        KATAURA.Type = [KATAURA.Type, type];
    end
end

% Plotting with different symbols for metallic and semiconducting tubes
figure;
hold on;
title('Excitation/Emission Wavelengths vs. RBM');
xlabel('RBM frequency (cm^{-1})');
ylabel('Wavelength (nm)');

colors = {[0 0.4470 0.7410], [0.4660 0.6740 0.1880], [0.6350 0.0780 0.1840], [0.9290 0.6940 0.1250]}; % Elegant colors

for i = 1:length(KATAURA.RBM)
    if KATAURA.Type(i) == 0
        marker = 'o'; % Circle for metallic
    else
        marker = 's'; % Square for semiconducting
    end
    scatter(KATAURA.RBM(i), KATAURA.WL1(i), marker, '  text(KATAURA.RBM(i), KATAURA.WL1(i), KATAURA.Chirality{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    text(KATAURA.RBM(i), KATAURA.WL2(i), KATAURA.Chirality{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    text(KATAURA.RBM(i), KATAURA.WL3(i), KATAURA.Chirality{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');', colors{1});
    scatter(KATAURA.RBM(i), KATAURA.WL2(i), marker, 'MarkerEdgeColor', colors{2});
    scatter(KATAURA.RBM(i), KATAURA.WL3(i), marker, 'MarkerEdgeColor', colors{3});
    scatter(KATAURA.RBM(i), KATAURA.WL4(i), marker, 'MarkerEdgeColor', colors{4});
    text(KATAURA.RBM(i), KATAURA.WL1(i), KATAURA.Chirality{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 8);
end

hold off;

% Function definition
function [nuRBM, wl1, wl2, wl3, wl4, diam, theta, type] = CalculateKataura(P)
    n = P(1);
    m = P(2);

    diam = 0.144 * sqrt(3) * (sqrt(n^2 + m^2 + n * m)) / pi;
    theta = atan(sqrt(3) * m / (m + 2 * n));
    nuRBM = (223.5 / diam) + 12.5;

    a = 1.049; % eV nm
    b = 0.456;
    c = 0.812; % nm^-1

    if mod(n - m, 3) == 0 % METALLIC TUBES
        energy1 = ((a * 3 / diam) * (1 + (b * log10(c / (3 / diam))))) - 0.18 * cos(3 * theta) / diam^2;
        wl1 = 1240 / energy1;
        energy2 = ((a * 3 / diam) * (1 + (
