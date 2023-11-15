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


if 1==1
    for n = 4:20 % Adjust the range as needed
        for m = 1:n % Adjust the range as needed
            result = mod(n-m,3);
            if (result==1 || result==2)
                [l11, l22] = calculate_lambdas(n, m);
                m_array = [m_array, m];
                n_array = [n_array, n];
                l11_array = [l11_array, l11]; % Save value in the array
                l22_array = [l22_array, l22]; % Save value in the array
                D_array = [D_array,calculateD(n,m)];
                InvreseD_array = [InvreseD_array,1/calculateD(n,m)];
                RBM_array = [RBM_array,calculateRBM(n,m)];
            end
        end
    end
    
    figure;
    subplot(1, 1, 1);
    scatter(InvreseD_array, RBM_array)
    xlabel('Inverse of "d" (nm-1)')
    ylabel('RMB frequency (cm-1) ')

   
    figure;
    subplot(1, 1, 1);
    scatter(RBM_array, l11_array); hold on
    scatter(RBM_array, l22_array)
    xlabel('RMB (cm-1)')
    ylabel('WL (nm)')
    
    figure;
    subplot(1, 1, 1);
    scatter(l11_array, l22_array)
    xlabel('Lambda_{11} (nm)')
    ylabel('Lambda_{22} (nm)')

end

keySet = {  'Water1',
            'Water2', 
            'Water3', 
            'Water4'};
        
valueSet = {[327.2, 327.2],
            [327.2, 327.2], 
            [327.2, 327.2], 
            [327.2, 327.2]};
        
M = containers.Map(keySet, valueSet);


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


