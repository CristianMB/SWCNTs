% Step 1: Load the absorption data from a .txt file
disp('Select the file containing the absorption data');
filename = uigetfile('*.txt', 'Select the file containing the absorption data');

% Format for reading the data (assuming two columns: wavelength, absorption)
string = '%f32 %f32'; 

% Open the file and read the data
fid = fopen(filename);
D = textscan(fid, string, 'Delimiter', ',', 'HeaderLines', 2);
fclose(fid);

% Assign the loaded data to the 'Data' variable
Data(:,1) = D{1}; % Wavelength
Data(:,2) = D{2}; % Absorption

% Ensure the data is in ascending order (if needed, flip the data)
if Data(2,1) < Data(1,1)
    Data = flipud(Data);
end

% Step 2: Interpolate the data
xx=round(Data(1,1)):1:round(Data(end,1));
yy=interp1(Data(:,1),Data(:,2),xx);
xx=xx'; yy=yy';


% Step 10: Plot the results
figure; hold on;
plot(xx, yy, 'k', 'LineWidth', 1); % Original data
xlabel('Wavelength (nm)', 'FontSize', 15);
ylabel('Absorption (a.u.)', 'FontSize', 15);
legend('Measured Data', 'Background');
set(gca, 'LineWidth', 1);
set(gca, 'FontSize', 15);
xlim(q)

% Step 3: Define the range for background subtraction
% q is assumed to be a two-element vector specifying the start and end wavelengths
q = [264, 2500];  % You should define these values

start=find(xx==q(1));
End=find(xx==q(2));

% Step 4: Write temporary data to a text file (for Naumov method)
dlmwrite('test.txt',[xx(start:End,1) yy(start:End,1)]);

% Step 5: Optimization options for fmincon
options = optimoptions('fmincon','Algorithm','interior-point','Display','off');

% Step 6: Perform the optimization to find A and b (initial guess: [0.2, 0.00002])
xn = fmincon(@Naumov,[.2 0.00002],[],[],[],[],[],[],@conf_naumov,options); % % Approach of Naumov et al. requires two starting values for A (0.2) and b (0.002).

% Step 7: Define the background subtraction model
F_naumov=@(x) double(yy(start:End,1)-x(1)*exp(-x(2).*xx(start:End,1))); % Background subtracted data!

% Step 8: Subtract the background
BKG = yy(start:End,1) - F_naumov(xn);

% Step 9: Delete temporary file after use
delete('test.txt');

% Step 10: Plot the results
figure; hold on;
plot(xx(start:End,1), yy(start:End,1), 'k', 'LineWidth', 1); % Original data
plot(xx(start:End,1), BKG, 'r', 'LineWidth', 1); % Background subtracted data
xlabel('Wavelength (nm)', 'FontSize', 15);
ylabel('Absorption (a.u.)', 'FontSize', 15);
legend('Measured Data', 'Background');
set(gca, 'LineWidth', 1);
set(gca, 'FontSize', 15);
xlim(q)
ylim([min(BKG)*0.9 max(BKG)*1.1])


% Step 11: Define Naumov error function for fmincon
function err = Naumov(x) % Calculates the difference between the absorption data and the background. The MATLAB function "fmincon" tries to minimize this difference by fitting x(1)=A and x(2)=b

A=dlmread('test.txt');
c = A(:,2)-x(1)*exp(-x(2).*A(:,1));
err = double(sum(c));

end

function [c,ceq] = conf_naumov(x) % Constraint function, that forces the background to be smaller than the absorption data for every single wavelength

A=dlmread('test.txt');
% Nonlinear inequality constraints
c = double(x(1)*exp(-x(2).*A(:,1))-A(:,2));
% Nonlinear equality constraints
ceq = [];
end

