%% Initialize
clc
close all

fs = 16000;                       % sampling freq
N = 4*fs;                         % number of samples
f = 121;                           % Signal frequency
A = [0.5 1 2 4 8 10];             % Amplitude
lhh = 1;                          % Low High Hot

%% FE method
for ii = 1:6
    [XFE(:,:,ii),eigvalFE(ii,:)] = part2a(fs,N,f,A(ii),"FE",lhh);
end
%% MP method
for ii = 1:6
    [XMP(:,:,ii),eigvalMP(ii,:)] = part2a(fs,N,f,A(ii),"MP",lhh);
end
%% calculate THD for FEM
NDFT = length(XFE(4,2*fs+1:end,1));
w = hann(NDFT)';
wsum = sum(w);
fv = (0:NDFT/2-1)*fs/NDFT;
for ii = 1:6
    XF = fft(XFE(4,2*fs+1:end,ii))/wsum;
    [FErms_121(ii,:),FEdB_121(ii,:)] = vel2pres(XF(1:NDFT/2),fv);
    [THDFE121(ii),~] = Task_3(FErms_121(ii,:),f,fs,N,10,1);
    dispFE121(ii,:) = [min(XFE(3,2*fs+1:end,ii)) max(XFE(3,2*fs+1:end,ii))];
end

%% calculate THD for MP
NDFT = length(XMP(4,2*fs+1:end,1));
w = hann(NDFT)';
wsum = sum(w);
fv = (0:NDFT/2-1)*fs/NDFT;
for ii = 1:6
    XF = fft(XMP(4,2*fs+1:end,ii))/wsum;
    [MPrms_121(ii,:),MPdB_121(ii,:)] = vel2pres(XF(1:NDFT/2),fv);
    [THDMP121(ii),~] = Task_3(MPrms_121(ii,:),f,fs,N,10,1);
    dispMP121(ii,:) = [min(XMP(3,2*fs+1:end,ii)) max(XMP(3,2*fs+1:end,ii))];
end

%% Compare simulated THD (FE) Vs Measured THD
voltvec = [0.5 1 2 4 8 10];
figure
load('30Hz_sim_FE_MP.mat')
plot(voltvec,THDFE30,'r','LineStyle','-','Marker','o','LineWidth',1.2)
hold on
load('60Hz_sim_FE_MP.mat')
plot(voltvec,THDFE60,'b','LineStyle','-','Marker','^','LineWidth',1.2)
hold on
load('121Hz_sim_FE_MP.mat')
plot(voltvec,THDFE121,'k','LineStyle','-','Marker','d','LineWidth',1.2)
hold on
load('THDMeas.mat')
plot(voltvec,THDr30,'r','LineStyle','--','Marker','o','LineWidth',1.2)
hold on
plot(voltvec,THDr60,'b','LineStyle','--','Marker','^','LineWidth',1.2)
hold on
plot(voltvec,THDr121,'k','LineStyle','--','Marker','d','LineWidth',1.2)
legend('Sim 30Hz','Sim 60Hz','Sim 121Hz','Meas 30Hz','Meas 60Hz','Meas 121Hz')
hold off
xlabel('Voltage (Volts)')
ylabel('THD_{R}(%)')
title('THD_{R} - Simulated (Forward Euler) Vs Measured')
xticks(voltvec)
grid minor
set(gca,'FontSize',16)
set(gcf,'position',[50 50 800 600])
clear

%% Compare simulated THD (MP) Vs Measured THD
voltvec = [0.5 1 2 4 8 10];
figure
load('30Hz_sim_FE_MP.mat')
plot(voltvec,THDMP30,'r','LineStyle','-','Marker','o','LineWidth',1.2)
hold on
load('60Hz_sim_FE_MP.mat')
plot(voltvec,THDMP60,'b','LineStyle','-','Marker','^','LineWidth',1.2)
hold on
load('121Hz_sim_FE_MP.mat')
plot(voltvec,THDMP121,'k','LineStyle','-','Marker','d','LineWidth',1.2)
hold on
load('THDMeas.mat')
plot(voltvec,THDr30,'r','LineStyle','--','Marker','o','LineWidth',1.2)
hold on
plot(voltvec,THDr60,'b','LineStyle','--','Marker','^','LineWidth',1.2)
hold on
plot(voltvec,THDr121,'k','LineStyle','--','Marker','d','LineWidth',1.2)
legend('Sim 30Hz','Sim 60Hz','Sim 121Hz','Meas 30Hz','Meas 60Hz','Meas 121Hz')
hold off
xlabel('Voltage (Volts)')
ylabel('THD_{R}(%)')
title('THD_{R} - Simulated (Mid-point) Vs Measured')
xticks(voltvec)
grid minor
set(gca,'FontSize',16)
set(gcf,'position',[50 50 800 600])
clear

%% Plot the maximum excursions FE vs MP
voltvec = [0.5 1 2 4 8 10];
figure
load('30Hz_sim_FE_MP.mat')
plot(voltvec,dispFE30(:,2)*1e+3,'r','LineStyle','-','Marker','o','LineWidth',1.2)
hold on
load('60Hz_sim_FE_MP.mat')
plot(voltvec,dispFE60(:,2)*1e+3,'b','LineStyle','-','Marker','^','LineWidth',1.2)
hold on
load('121Hz_sim_FE_MP.mat')
plot(voltvec,dispFE121(:,2)*1e+3,'k','LineStyle','-','Marker','d','LineWidth',1.2)
hold on
plot(voltvec,dispMP30(:,1)*1e+3,'r','LineStyle','--','Marker','o','LineWidth',1.2)
hold on
plot(voltvec,dispMP60(:,1)*1e+3,'b','LineStyle','--','Marker','^','LineWidth',1.2)
hold on
plot(voltvec,dispMP121(:,1)*1e+3,'k','LineStyle','--','Marker','d','LineWidth',1.2)
legend('FE 30Hz','FE 60Hz','FE 121Hz','MP 30Hz','MP 60Hz','MP 121Hz')
hold off
xlabel('Voltage (Volts)')
ylabel('Diaphragm Excursion (mm)')
title('Diaphragm excursion Simulated - Forward Euler Vs Mid point')
xticks(voltvec)
grid minor
set(gca,'FontSize',16)
set(gcf,'position',[50 50 800 600])
clear

%% Plot the levels of fundamentals
voltvec = [0.5 1 2 4 8 10];
fs = 96000;
fvmeas = 0:fs/2-1;
ind = find(fvmeas==30);
load('PressMeas.mat')
for ii = 1:length(voltvec)
    pmeas(1,ii)= pS11(ii).prmsdB(ind);
end
figure
subplot(2,2,1)
load('30Hz_sim_FE_MP.mat')
ind = find(fv==30);
plot(voltvec,FEdB_30(:,ind),'r','LineStyle','-','Marker','o','LineWidth',1.2)
hold on
plot(voltvec,pmeas(1,:),'k','LineStyle','--','Marker','^','LineWidth',1.2)
xlabel('Voltage (Volts)')
ylabel('SPL (dB re 1Pa)')
title('FE, f1 = 30 Hz')
xticks(voltvec)
grid minor
set(gca,'FontSize',15)

ind = find(fvmeas==60);
for ii = 1:length(voltvec)
    pmeas(1,ii)= pS12(ii).prmsdB(ind);
end
subplot(2,2,2)
load('60Hz_sim_FE_MP.mat')
ind = find(fv==60);
plot(voltvec,FEdB_60(:,ind),'r','LineStyle','-','Marker','o','LineWidth',1.2)
hold on
plot(voltvec,pmeas(1,:),'k','LineStyle','--','Marker','^','LineWidth',1.2)
xlabel('Voltage (Volts)')
ylabel('SPL (dB re 1Pa)')
title('FE, f1 = 30 Hz')
xticks(voltvec)
grid minor
set(gca,'FontSize',15)


ind = find(fvmeas==121);
for ii = 1:length(voltvec)
    pmeas(1,ii)= pS13(ii).prmsdB(ind);
end
subplot(2,2,3)
load('121Hz_sim_FE_MP.mat')
ind = find(fv==121);
plot(voltvec,FEdB_121(:,ind),'r','LineStyle','-','Marker','o','LineWidth',1.2)
hold on
plot(voltvec,pmeas(1,:),'k','LineStyle','--','Marker','^','LineWidth',1.2)
xlabel('Voltage (Volts)')
ylabel('SPL (dB re 1Pa)')
title('FE, f1 = 121 Hz')
xticks(voltvec)
legend('Simulated','Measured')
grid minor
set(gca,'FontSize',15)
set(gcf,'position',[50 50 800 600])
clear
