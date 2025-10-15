%% Step 1: Specify folder and get list of files
folder_path = 'X:\Measurements Data\TGA\';  % folder containing your txt files
addpath('X:\SWCNTs\TGA_Analysis\');  % folder containing your txt files

file_list = dir(fullfile(folder_path, '*.txt'));

file_names = {
    '20250825-01 TCNQ filled.txt'
    '20250826-01 P2 unfilled.txt'
    '20250826-02 TTF filled.txt'
    '20251002-01 CB-HS1.txt'
    '20251014-02 TTF_B.txt'
    '20251014-03 TCNQ_B.txt'
};

aliases = {'TCNQ_GasPhase','P2_Reference','TTF_GasPhase','TTF_VacInfil','TTF_Burnt','TCNQ_Burnt'};


num_files = length(file_names);

if num_files == 0
    error('No .txt files found in the selected folder.');
end

%% Step 2: Read each file and extract Temperature & Weight
all_data = cell(num_files,1);  % each cell = [Temperature, Weight]

for k = 1:num_files
    file_path = fullfile(folder_path, file_names{k});
    
    % Read file as text
    fid = fopen(file_path, 'r');
    file_text = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);
    lines = file_text{1};
    
    % Skip headers and last non-numeric line
    data_lines = lines(3:end-1);
    num_rows = length(data_lines);
    
    temp_weight = zeros(num_rows,2);  % preallocate
    
    % Loop through each line
    for i = 1:num_rows
        line = strtrim(data_lines{i});             % remove leading/trailing spaces
        line = strrep(line, ',', '.');             % replace comma with dot
        tokens = regexp(line, '\s+', 'split');     % split by spaces
        
        % Extract Temperature (column 2) and Weight (column 5)
        temp_weight(i,1) = str2double(tokens{2});
        temp_weight(i,2) = str2double(tokens{5});
    end
    
    all_data{k} = temp_weight;
end

%% Step 3: Combine all samples into a single matrix (2 columns per file)
% First, find the maximum number of rows
max_rows = max(cellfun(@(x) size(x,1), all_data));

% Initialize combined matrix with NaNs
combined_matrix = NaN(max_rows, 2*num_files);


for k = 1:num_files
    n_rows = size(all_data{k},1);
    combined_matrix(1:n_rows, 2*k-1:2*k) = all_data{k};  % Temperature & Weight
end


%% Display first few rows of combined matrix

TGA = combined_matrix;
TGA(isnan(TGA)) = 0;

%%FIX FOR LAST MEASUREMENT

mass = TGA(:,8);           % current mass in mg
initial_mass = mass(1);     % take the first row as initial mass
TGA(:,8) = 100*(1-((initial_mass-mass)/initial_mass));

mass = TGA(:,10);           % current mass in mg
initial_mass = mass(1);     % take the first row as initial mass
TGA(:,10) = 100*(1-((initial_mass-mass)/initial_mass));

mass = TGA(:,12);           % current mass in mg
initial_mass = mass(1);     % take the first row as initial mass
TGA(:,12) = 100*(1-((initial_mass-mass)/initial_mass));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read in data and plot


% DTGA=TGA;
% DTGA(:,2)=gradient(TGA(:,2),TGA(:,1));
% DTGA(:,4)=gradient(TGA(:,4),TGA(:,3));
% DTGA(:,6)=gradient(TGA(:,6),TGA(:,5));
% DTGA(:,8)=gradient(TGA(:,8),TGA(:,7));

num_samples = size(TGA,2)/2;

%Individual Plots
for k = 1:num_samples
    % Extract Temperature & Weight
    temp = TGA(600:end, 2*k-1);
    weight = TGA(600:end, 2*k);
    
    % Only keep non-zero rows for plotting/fitting
    valid_idx = (temp ~= 0 & weight ~= 0);
    temp_plot = temp(valid_idx);
    weight_plot = weight(valid_idx);
    
    % Compute DTGA
    dWdT  = gradient(weight_plot, temp_plot);
    
    alias_name = aliases{k};
    DATA.(alias_name).T = temp_plot;
    DATA.(alias_name).W = weight_plot;
    DATA.(alias_name).dWdT = dWdT;
    DATA.(alias_name).N = alias_name;  % store alias as string
    
    % Smooth DTGA (optional)
    Y = -datasmooth(dWdT, 3, 'binom');
    DATA.(alias_name).dWdT = dWdT;
        
    % Plot TGA and DTGA separately
%     figure(k); clf; set(gcf,'color','w','Position',[100 100 700 500])
% 
%     subplot(2,1,1)
%     plot(temp_plot, weight_plot, 'k', 'LineWidth', 2);
%     xlabel('Temperature (°C)')
%     ylabel('Weight Loss (%)')
%     title(sprintf('Sample %d TGA', k))
%     set(gca,'TickDir','in','XMinorTick','on','YMinorTick','on','LineWidth',1.5,'FontSize',12)
%     xlim([100, 1000])
%     ylim([0, 100])
%     
%     subplot(2,1,2)
%     plot(temp_plot, -dWdT, 'b', 'LineWidth', 2);
%     xlabel('Temperature (°C)')
%     ylabel('DTGA (a.u.)')
%     title(sprintf('Sample %d DTGA', k))
%     set(gca,'TickDir','in','XMinorTick','on','YMinorTick','on','LineWidth',1.5,'FontSize',12)
%     xlim([100, 1000])
end
% close all;

% figure()
% plot(DATA.TTF_GasPhase.T, -DATA.TTF_GasPhase.dWdT)
% hold on
% plot(DATA.TTF_VacInfil.T, -DATA.TTF_VacInfil.dWdT)
% 
% 
% figure()
% plot(DATA.TTF_GasPhase.T, DATA.TTF_GasPhase.W)
% hold on
% plot(DATA.TTF_VacInfil.T, DATA.TTF_VacInfil.W)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


X1=DATA.P2_Reference.T;
Y1=-datasmooth(DATA.P2_Reference.dWdT,5,'binom');
N1=4;
Center=[203 596 628 729];
fwhm=[62 84 46 118];
x0=[Center fwhm];
UB=[Center+10 fwhm*3];
LB=[Center-10 fwhm*0.5];
[result1]=lsqcurvefit(@(x0,X1) FitGaussians(x0,X1,N1,Y1),x0,X1,Y1,LB,UB);
[FIT1,L1,A1]=FitGaussians(result1,X1,N1,Y1);

%
X2=DATA.TTF_GasPhase.T;
Y2=-datasmooth(DATA.TTF_GasPhase.dWdT,3,'binom');
N2=5;
Center=[259 440 594 620 721];
fwhm=[95 181 101 52 122];
x02=[Center fwhm];
UB2=[Center+10 fwhm*3];
LB2=[Center-10 fwhm*0.5];
[result2]=lsqcurvefit(@(x02,X2) FitGaussians(x02,X2,N2,Y2),x02,X2,Y2,LB2,UB2);
[FIT2,L2,A2]=FitGaussians(result2,X2,N2,Y2);

X3=DATA.TCNQ_GasPhase.T;
Y3=-datasmooth(DATA.TCNQ_GasPhase.dWdT,3,'binom');
N3=5;
Center=[246 260 342 545 685];
fwhm=[16 15 80  99 88];
x03=[Center fwhm];
UB3=[Center+10 fwhm*2];
LB3=[Center-10 fwhm*0.5];
[result3]=lsqcurvefit(@(x03,X3) FitGaussians(x03,X3,N3,Y3),x03,X3,Y3,LB3,UB3);
[FIT3,L3,A3]=FitGaussians(result3,X3,N3,Y3);


X4=DATA.TTF_VacInfil.T;
Y4=-datasmooth(DATA.TTF_VacInfil.dWdT,3,'binom');
N4=5;
Center=[276 440 594 620 721];
fwhm=[200 181 101 52 122];
x04=[Center fwhm];
UB4=[Center+10 fwhm*3];
LB4=[Center-10 fwhm*0.5];
[result4]=lsqcurvefit(@(x04,X4) FitGaussians(x04,X4,N4,Y4),x04,X4,Y4,LB4,UB4);
[FIT4,L4,A4]=FitGaussians(result4,X4,N4,Y4);


X5=DATA.TTF_Burnt.T;
Y5=-datasmooth(DATA.TTF_Burnt.dWdT,3,'binom');
N5=5;
Center=[276 440 594 620 721];
fwhm=[200 181 101 52 122];
x05=[Center fwhm];
UB5=[Center+10 fwhm*3];
LB5=[Center-10 fwhm*0.5];
[result5]=lsqcurvefit(@(x05,X5) FitGaussians(x05,X5,N5,Y5),x05,X5,Y5,LB5,UB5);
[FIT5,L5,A5]=FitGaussians(result5,X5,N5,Y5);


X6=DATA.TCNQ_Burnt.T;
Y6=-datasmooth(DATA.TCNQ_Burnt.dWdT,3,'binom');
N6=5;
Center=[276 440 594 620 721];
fwhm=[200 181 101 52 122];
x06=[Center fwhm];
UB6=[Center+10 fwhm*3];
LB6=[Center-10 fwhm*0.5];
[result6]=lsqcurvefit(@(x06,X6) FitGaussians(x06,X6,N6,Y6),x06,X6,Y6,LB6,UB6);
[FIT6,L6,A6]=FitGaussians(result6,X6,N6,Y6);

%%


figure(1);clf;set(gcf,'color','w','Position',[50 50 800 700])

subplot(num_samples,2,1)
plot(DATA.P2_Reference.T,DATA.P2_Reference.W,'k','LineWidth',2); hold on;
xlabel('Temperature (°C)')
ylabel('Weight Loss (%)')
% legend('EA-SWCNTs','Box','off','Location','southwest','IconColumnWidth',10)
xlim([30 805])
ylim([0 110])
    set(gca,'TickDir','in','XMinortick','on','YMinorTick','on','TickLength',[0.03 0.05],'LineWidth',1.5,'FontSize',12)


subplot(num_samples,2,2)
plot(X1,Y1,'k','LineWidth',2);  hold on;
for i=1:N1
plot(X1,L1(:,i).*A1(i)+L1(:,N1+1).*A1(N1+1),'Color',[0.5 0.5 0.5],'LineWidth',1.5)
end
plot(X1,FIT1,'g','LineWidth',2);
xlabel('Temperature (°C)')
ylabel('DTGA (a.u.)')
xlim([120 805])
ylim([0 0.6])
xline(result1(N1),'--')
    set(gca,'TickDir','in','XMinortick','on','YMinorTick','on','TickLength',[0.03 0.05],'LineWidth',1.5,'FontSize',12)


subplot(num_samples,2,3)
plot(DATA.TTF_GasPhase.T,DATA.TTF_GasPhase.W,'k','LineWidth',2); hold on;
xlabel('Temperature (°C)')
ylabel('Weight Loss (%)')
% legend('TTF@SWCNTs','Box','off','Location','southwest','IconColumnWidth',10)
xlim([30 805])
ylim([0 110])
    set(gca,'TickDir','in','XMinortick','on','YMinorTick','on','TickLength',[0.03 0.05],'LineWidth',1.5,'FontSize',12)


subplot(num_samples,2,4)
plot(X2,Y2,'k','LineWidth',2);  hold on;
for i=1:N2
    if i==1
        plot(X2,L2(:,i).*A2(i)+L2(:,N2+1).*A2(N2+1),'r','LineWidth',3)
    elseif i==2
        plot(X2,L2(:,i).*A2(i)+L2(:,N2+1).*A2(N2+1),'b','LineWidth',3)
    else
    plot(X2,L2(:,i).*A2(i)+L2(:,N2+1).*A2(N2+1),'Color',[0.5 0.5 0.5],'LineWidth',1.5)
    end
end
plot(X2,FIT2,'g','LineWidth',2);
xlabel('Temperature (°C)')
ylabel('DTGA (a.u.)')
xlim([120 805])
ylim([0 0.4])

TOTALWEIGHT = sum(A2(1:N2))+TGA(end,6)+(100-TGA(600,6))
TTFoutside=100.*A2(1)/TOTALWEIGHT
TTFinside=100.*A2(2)/TOTALWEIGHT
text(230,0.13,strcat(num2str(round(TTFoutside,2)),'%'),'color','r')
text(390,0.13,strcat(num2str(round(TTFinside,2)),'%'),'color','b')
xline(result1(N1),'--')
    set(gca,'TickDir','in','XMinortick','on','YMinorTick','on','TickLength',[0.03 0.05],'LineWidth',1.5,'FontSize',12)


subplot(num_samples,2,5)
plot(DATA.TTF_VacInfil.T,DATA.TTF_VacInfil.W,'k','LineWidth',2); hold on;
xlabel('Temperature (°C)')
ylabel('Weight Loss (%)')
% legend('TCNQ@SWCNTs','Box','off','Location','southwest','IconColumnWidth',10)
xlim([30 805])
ylim([0 110])
    set(gca,'TickDir','in','XMinortick','on','YMinorTick','on','TickLength',[0.03 0.05],'LineWidth',1.5,'FontSize',12)


subplot(num_samples,2,6)
plot(X4,Y4,'k','LineWidth',2);  hold on;
for i=1:N4
    if i==1
        plot(X4,L4(:,i).*A4(i)+L4(:,N4+1).*A4(N4+1),'r','LineWidth',3)
    elseif i==2
        plot(X4,L4(:,i).*A4(i)+L4(:,N4+1).*A4(N4+1),'b','LineWidth',3)
    else
    plot(X4,L4(:,i).*A4(i)+L4(:,N4+1).*A4(N4+1),'Color',[0.5 0.5 0.5],'LineWidth',1.5)
    end
end
plot(X4,FIT4,'g','LineWidth',2);
xlabel('Temperature (°C)')
ylabel('DTGA (a.u.)')
xlim([120 805])
ylim([0 0.45])
xline(result1(N1),'--')

TOTALWEIGHT = sum(A4(1:N4))+TGA(end,8)+(100-TGA(600,8))
TTFoutside_VI=100.*A4(1)/TOTALWEIGHT
TTFinside_VI=100.*A4(2)/TOTALWEIGHT
text(230,0.13,strcat(num2str(round(TTFoutside_VI,2)),'%'),'color','r')
text(390,0.13,strcat(num2str(round(TTFinside_VI,2)),'%'),'color','b')
xline(result1(N1),'--')
set(gca,'TickDir','in','XMinortick','on','YMinorTick','on','TickLength',[0.03 0.05],'LineWidth',1.5,'FontSize',12)


subplot(num_samples,2,7)
plot(DATA.TTF_Burnt.T,DATA.TTF_Burnt.W,'k','LineWidth',2); hold on;
xlabel('Temperature (°C)')
ylabel('Weight Loss (%)')
% legend('TTF@SWCNTs','Box','off','Location','southwest','IconColumnWidth',10)
xlim([30 805])
ylim([0 110])
    set(gca,'TickDir','in','XMinortick','on','YMinorTick','on','TickLength',[0.03 0.05],'LineWidth',1.5,'FontSize',12)


subplot(num_samples,2,8)
plot(X5,Y5,'k','LineWidth',2);  hold on;
for i=1:N5
    if ismember(i,[1 2])
        plot(X5,L5(:,i).*A5(i)+L5(:,N5+1).*A5(N5+1),'r','LineWidth',2)
    elseif i==3
        plot(X5,L5(:,i).*A5(i)+L5(:,N5+1).*A5(N5+1),'b','LineWidth',2)
    else
    plot(X5,L5(:,i).*A5(i)+L5(:,N5+1).*A5(N5+1),'Color',[0.5 0.5 0.5],'LineWidth',1.5)
    end
end
plot(X5,FIT5,'g','LineWidth',2)
xlabel('Temperature (°C)')
ylabel('DTGA (a.u.)')
xlim([120 805])
ylim([0 0.9])
xline(result1(N1),'--')

TOTALWEIGHT = sum(A5(1:N5))+TGA(end,5)+(100-TGA(600,5))
TTFOutside=100.*(A5(1)+A5(2))/TOTALWEIGHT;
TTFInside=100.*A5(3)/TOTALWEIGHT;
text(200,0.8,strcat(num2str(round(TTFOutside,2)),'%'),'color','r')
text(300,0.13,strcat(num2str(round(TTFInside,2)),'%'),'color','b')
set(gca,'TickDir','in','XMinortick','on','YMinorTick','on','TickLength',[0.03 0.05],'LineWidth',1.5,'FontSize',12)
    

subplot(num_samples,2,9)
plot(DATA.TCNQ_GasPhase.T,DATA.TCNQ_GasPhase.W,'k','LineWidth',2); hold on;
xlabel('Temperature (°C)')
ylabel('Weight Loss (%)')
% legend('TCNQ@SWCNTs','Box','off','Location','southwest','IconColumnWidth',10)
xlim([30 805])
ylim([0 110])
    set(gca,'TickDir','in','XMinortick','on','YMinorTick','on','TickLength',[0.03 0.05],'LineWidth',1.5,'FontSize',12)


subplot(num_samples,2,10)
plot(X3,Y3,'k','LineWidth',2);  hold on;
for i=1:N3
    if ismember(i,[1 2])
        plot(X3,L3(:,i).*A3(i)+L3(:,N3+1).*A3(N3+1),'r','LineWidth',2)
    elseif i==3
        plot(X3,L3(:,i).*A3(i)+L3(:,N3+1).*A3(N3+1),'b','LineWidth',2)
    else
    plot(X3,L3(:,i).*A3(i)+L3(:,N3+1).*A3(N3+1),'Color',[0.5 0.5 0.5],'LineWidth',1.5)
    end
end
plot(X3,FIT3,'g','LineWidth',2)
xlabel('Temperature (°C)')
ylabel('DTGA (a.u.)')
xlim([120 805])
ylim([0 0.9])
xline(result1(N1),'--')

TOTALWEIGHT = sum(A3(1:N3))+TGA(end,4)+(100-TGA(600,4))
TCNQoutside=100.*(A3(1)+A3(2))/TOTALWEIGHT;
TCNQinside=100.*A3(3)/TOTALWEIGHT;
text(200,0.8,strcat(num2str(round(TCNQoutside,2)),'%'),'color','r')
text(300,0.13,strcat(num2str(round(TCNQinside,2)),'%'),'color','b')
set(gca,'TickDir','in','XMinortick','on','YMinorTick','on','TickLength',[0.03 0.05],'LineWidth',1.5,'FontSize',12)


subplot(num_samples,2,11)
plot(DATA.TCNQ_GasPhase.T,DATA.TCNQ_GasPhase.W,'k','LineWidth',2); hold on;
xlabel('Temperature (°C)')
ylabel('Weight Loss (%)')
% legend('TCNQ@SWCNTs','Box','off','Location','southwest','IconColumnWidth',10)
xlim([30 805])
ylim([0 110])
    set(gca,'TickDir','in','XMinortick','on','YMinorTick','on','TickLength',[0.03 0.05],'LineWidth',1.5,'FontSize',12)


subplot(num_samples,2,12)
plot(X3,Y3,'k','LineWidth',2);  hold on;
for i=1:N3
    if ismember(i,[1 2])
        plot(X3,L3(:,i).*A3(i)+L3(:,N3+1).*A3(N3+1),'r','LineWidth',2)
    elseif i==3
        plot(X3,L3(:,i).*A3(i)+L3(:,N3+1).*A3(N3+1),'b','LineWidth',2)
    else
    plot(X3,L3(:,i).*A3(i)+L3(:,N3+1).*A3(N3+1),'Color',[0.5 0.5 0.5],'LineWidth',1.5)
    end
end
plot(X3,FIT3,'g','LineWidth',2)
xlabel('Temperature (°C)')
ylabel('DTGA (a.u.)')
xlim([120 805])
ylim([0 0.9])
xline(result1(N1),'--')

TOTALWEIGHT = sum(A3(1:N3))+TGA(end,4)+(100-TGA(600,4))
TCNQoutside=100.*(A3(1)+A3(2))/TOTALWEIGHT;
TCNQinside=100.*A3(3)/TOTALWEIGHT;
text(200,0.8,strcat(num2str(round(TCNQoutside,2)),'%'),'color','r')
text(300,0.13,strcat(num2str(round(TCNQinside,2)),'%'),'color','b')
set(gca,'TickDir','in','XMinortick','on','YMinorTick','on','TickLength',[0.03 0.05],'LineWidth',1.5,'FontSize',12)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2);clf;set(gcf,'color','w','Position',[50 50 800 700])

subplot(2,2,1)
plot(DATA.TTF_GasPhase.T,DATA.TTF_GasPhase.W,'k','LineWidth',2); hold on;
xlabel('Temperature (°C)')
ylabel('Weight Loss (%)')
% legend('TTF@SWCNTs','Box','off','Location','southwest','IconColumnWidth',10)
xlim([30 805])
ylim([0 110])
    set(gca,'TickDir','in','XMinortick','on','YMinorTick','on','TickLength',[0.03 0.05],'LineWidth',1.5,'FontSize',12)


subplot(2,2,2)
plot(X2,Y2,'k','LineWidth',2);  hold on;
for i=1:N2
    if i==1
        plot(X2,L2(:,i).*A2(i)+L2(:,N2+1).*A2(N2+1),'r','LineWidth',3)
    elseif i==2
        plot(X2,L2(:,i).*A2(i)+L2(:,N2+1).*A2(N2+1),'b','LineWidth',3)
    else
    plot(X2,L2(:,i).*A2(i)+L2(:,N2+1).*A2(N2+1),'Color',[0.5 0.5 0.5],'LineWidth',1.5)
    end
end
plot(X2,FIT2,'g','LineWidth',2);
xlabel('Temperature (°C)')
ylabel('DTGA (a.u.)')
xlim([120 805])
ylim([0 0.4])
TOTALWEIGHT = sum(A2(1:N2))+TGA(end,6)+(100-TGA(600,6))
TTFoutside=100.*A2(1)/TOTALWEIGHT
TTFinside=100.*A2(2)/TOTALWEIGHT
text(230,0.13,strcat(num2str(round(TTFoutside,2)),'%'),'color','r')
text(390,0.13,strcat(num2str(round(TTFinside,2)),'%'),'color','b')
xline(result1(N1),'--')
    set(gca,'TickDir','in','XMinortick','on','YMinorTick','on','TickLength',[0.03 0.05],'LineWidth',1.5,'FontSize',12)



subplot(2,2,3)
plot(DATA.TTF_VacInfil.T,DATA.TTF_VacInfil.W,'k','LineWidth',2); hold on;
xlabel('Temperature (°C)')
ylabel('Weight Loss (%)')
% legend('TCNQ@SWCNTs','Box','off','Location','southwest','IconColumnWidth',10)
xlim([30 805])
ylim([0 110])
    set(gca,'TickDir','in','XMinortick','on','YMinorTick','on','TickLength',[0.03 0.05],'LineWidth',1.5,'FontSize',12)


subplot(2,2,4)
plot(X4,Y4,'k','LineWidth',2);  hold on;
for i=1:N4
    if i==1
        plot(X4,L4(:,i).*A4(i)+L4(:,N4+1).*A4(N4+1),'r','LineWidth',3)
    elseif i==2
        plot(X4,L4(:,i).*A4(i)+L4(:,N4+1).*A4(N4+1),'b','LineWidth',3)
    else
    plot(X4,L4(:,i).*A4(i)+L4(:,N4+1).*A4(N4+1),'Color',[0.5 0.5 0.5],'LineWidth',1.5)
    end
end
plot(X4,FIT4,'g','LineWidth',2);
xlabel('Temperature (°C)')
ylabel('DTGA (a.u.)')
xlim([120 805])
ylim([0 0.45])
xline(result1(N1),'--')

TOTALWEIGHT = sum(A4(1:N4))+TGA(end,8)+(100-TGA(600,8))
TTFoutside_VI=100.*A4(1)/TOTALWEIGHT
TTFinside_VI=100.*A4(2)/TOTALWEIGHT
text(230,0.13,strcat(num2str(round(TTFoutside_VI,2)),'%'),'color','r')
text(390,0.13,strcat(num2str(round(TTFinside_VI,2)),'%'),'color','b')
xline(result1(N1),'--')
set(gca,'TickDir','in','XMinortick','on','YMinorTick','on','TickLength',[0.03 0.05],'LineWidth',1.5,'FontSize',12)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3);clf;set(gcf,'color','w','Position',[50 50 800 700])

subplot(2,2,1)
plot(DATA.TTF_GasPhase.T,DATA.TTF_GasPhase.W,'k','LineWidth',2); hold on;
xlabel('Temperature (°C)')
ylabel('Weight Loss (%)')
% legend('TTF@SWCNTs','Box','off','Location','southwest','IconColumnWidth',10)
xlim([30 805])
ylim([0 110])
    set(gca,'TickDir','in','XMinortick','on','YMinorTick','on','TickLength',[0.03 0.05],'LineWidth',1.5,'FontSize',12)


subplot(2,2,2)
plot(X2,Y2,'k','LineWidth',2);  hold on;
for i=1:N2
    if i==1
        plot(X2,L2(:,i).*A2(i)+L2(:,N2+1).*A2(N2+1),'r','LineWidth',3)
    elseif i==2
        plot(X2,L2(:,i).*A2(i)+L2(:,N2+1).*A2(N2+1),'b','LineWidth',3)
    else
    plot(X2,L2(:,i).*A2(i)+L2(:,N2+1).*A2(N2+1),'Color',[0.5 0.5 0.5],'LineWidth',1.5)
    end
end
plot(X2,FIT2,'g','LineWidth',2);
xlabel('Temperature (°C)')
ylabel('DTGA (a.u.)')
xlim([120 805])
ylim([0 0.4])
TOTALWEIGHT = sum(A2(1:N2))+TGA(end,6)+(100-TGA(600,6))
TTFoutside=100.*A2(1)/TOTALWEIGHT
TTFinside=100.*A2(2)/TOTALWEIGHT
text(230,0.13,strcat(num2str(round(TTFoutside,2)),'%'),'color','r')
text(390,0.13,strcat(num2str(round(TTFinside,2)),'%'),'color','b')
xline(result1(N1),'--')
    set(gca,'TickDir','in','XMinortick','on','YMinorTick','on','TickLength',[0.03 0.05],'LineWidth',1.5,'FontSize',12)



subplot(2,2,3)
plot(DATA.TTF_Burnt.T,DATA.TTF_Burnt.W,'k','LineWidth',2); hold on;
xlabel('Temperature (°C)')
ylabel('Weight Loss (%)')
% legend('TCNQ@SWCNTs','Box','off','Location','southwest','IconColumnWidth',10)
xlim([30 805])
ylim([0 110])
    set(gca,'TickDir','in','XMinortick','on','YMinorTick','on','TickLength',[0.03 0.05],'LineWidth',1.5,'FontSize',12)


subplot(2,2,4)
plot(X5,Y5,'k','LineWidth',2);  hold on;
for i=1:N5
    if i==1
        plot(X4,L4(:,i).*A5(i)+L5(:,N5+1).*A5(N5+1),'r','LineWidth',3)
    elseif i==2
        plot(X5,L5(:,i).*A5(i)+L5(:,N5+1).*A5(N5+1),'b','LineWidth',3)
    else
    plot(X5,L5(:,i).*A5(i)+L5(:,N5+1).*A5(N5+1),'Color',[0.5 0.5 0.5],'LineWidth',1.5)
    end
end
plot(X5,FIT5,'g','LineWidth',2);
xlabel('Temperature (°C)')
ylabel('DTGA (a.u.)')
xlim([120 805])
ylim([0 0.45])
xline(result1(N1),'--')

TOTALWEIGHT = sum(A5(1:N5))+TGA(end,10)+(100-TGA(600,10))
TTFoutside_VI=100.*A5(1)/TOTALWEIGHT
TTFinside_VI=100.*A5(2)/TOTALWEIGHT
text(230,0.13,strcat(num2str(round(TTFoutside_VI,2)),'%'),'color','r')
text(390,0.13,strcat(num2str(round(TTFinside_VI,2)),'%'),'color','b')
xline(result1(N1),'--')
set(gca,'TickDir','in','XMinortick','on','YMinorTick','on','TickLength',[0.03 0.05],'LineWidth',1.5,'FontSize',12)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4);clf;set(gcf,'color','w','Position',[50 50 800 700])

subplot(2,2,1)
plot(DATA.TCNQ_GasPhase.T,DATA.TCNQ_GasPhase.W,'k','LineWidth',2); hold on;
xlabel('Temperature (°C)')
ylabel('Weight Loss (%)')
% legend('TCNQ@SWCNTs','Box','off','Location','southwest','IconColumnWidth',10)
xlim([30 805])
ylim([0 110])
    set(gca,'TickDir','in','XMinortick','on','YMinorTick','on','TickLength',[0.03 0.05],'LineWidth',1.5,'FontSize',12)


subplot(2,2,2)
plot(X3,Y3,'k','LineWidth',2);  hold on;
for i=1:N3
    if i==1
        plot(X3,L3(:,i).*A3(i)+L3(:,N3+1).*A3(N3+1),'r','LineWidth',3)
    elseif i==2
        plot(X3,L3(:,i).*A3(i)+L2(:,N3+1).*A3(N3+1),'b','LineWidth',3)
    else
    plot(X3,L3(:,i).*A3(i)+L3(:,N3+1).*A3(N3+1),'Color',[0.5 0.5 0.5],'LineWidth',1.5)
    end
end
plot(X3,FIT3,'g','LineWidth',2);
xlabel('Temperature (°C)')
ylabel('DTGA (a.u.)')
xlim([120 805])
ylim([0 0.4])
TOTALWEIGHT = sum(A2(1:N2))+TGA(end,4)+(100-TGA(600,4))
TTFoutside=100.*A2(1)/TOTALWEIGHT
TTFinside=100.*A2(2)/TOTALWEIGHT
text(230,0.13,strcat(num2str(round(TTFoutside,2)),'%'),'color','r')
text(390,0.13,strcat(num2str(round(TTFinside,2)),'%'),'color','b')
xline(result1(N1),'--')
    set(gca,'TickDir','in','XMinortick','on','YMinorTick','on','TickLength',[0.03 0.05],'LineWidth',1.5,'FontSize',12)



subplot(2,2,3)
plot(DATA.TCNQ_Burnt.T,DATA.TCNQ_Burnt.W,'k','LineWidth',2); hold on;
xlabel('Temperature (°C)')
ylabel('Weight Loss (%)')
% legend('TCNQ@SWCNTs','Box','off','Location','southwest','IconColumnWidth',10)
xlim([30 805])
ylim([0 110])
    set(gca,'TickDir','in','XMinortick','on','YMinorTick','on','TickLength',[0.03 0.05],'LineWidth',1.5,'FontSize',12)


subplot(2,2,4)
plot(X6,Y6,'k','LineWidth',2);  hold on;
for i=1:N6
    if i==1
        plot(X6,L6(:,i).*A6(i)+L6(:,N6+1).*A6(N6+1),'r','LineWidth',3)
    elseif i==2
        plot(X6,L6(:,i).*A6(i)+L6(:,N6+1).*A6(N6+1),'b','LineWidth',3)
    else
    plot(X6,L6(:,i).*A6(i)+L6(:,N6+1).*A6(N6+1),'Color',[0.5 0.5 0.5],'LineWidth',1.5)
    end
end
plot(X6,FIT6,'g','LineWidth',2);
xlabel('Temperature (°C)')
ylabel('DTGA (a.u.)')
xlim([120 805])
ylim([0 0.45])
xline(result1(N1),'--')

TOTALWEIGHT = sum(A6(1:N6))+TGA(end,12)+(100-TGA(600,12))
TTFoutside_VI=100.*A6(1)/TOTALWEIGHT
TTFinside_VI=100.*A6(2)/TOTALWEIGHT
text(230,0.13,strcat(num2str(round(TTFoutside_VI,2)),'%'),'color','r')
text(390,0.13,strcat(num2str(round(TTFinside_VI,2)),'%'),'color','b')
xline(result1(N1),'--')
set(gca,'TickDir','in','XMinortick','on','YMinorTick','on','TickLength',[0.03 0.05],'LineWidth',1.5,'FontSize',12)

