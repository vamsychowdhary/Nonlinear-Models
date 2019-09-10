%% Initialize
clc
close all

fs = 96000;                       % sampling freq
N = 4*fs;                         % number of samples
f = 30;                           % Signal frequency
A = [0.5 1 2 4 8 10];             % Amplitude
lhh = 1;                          % Low High Hot

%% ODE45 Solver
t = 0:1/fs:(N-1)/fs;            % time vector
x0 = [0;0;0;0];
loadode2cdata;

[~,XODE30] = ode45(@ode2c,t,x0);
XODE30 = XODE30';
%% calculate THD for ODE
NDFT = length(XODE30(4,2*fs+1:end));
w = hann(NDFT)';
wsum = sum(w);
fv = (0:NDFT/2-1)*fs/NDFT;
XF = fft(XODE30(4,2*fs+1:end))/wsum;
[ODErms_30,ODEdB_30] = vel2pres(XF(1:NDFT/2),fv);
[THDODE_30,~] = Task_3(ODErms_30,f,fs,NDFT,10,1);
dispODE_30 = [min(XODE30(3,2*fs+1:end)) max(XODE30(3,2*fs+1:end))];

%% Plot
t = 0:1/fs:(N-1)/fs;            % time vector
tend = 5/f + 2;
tsind = find(t>=2,1);
teind = find(t>=tend,1);
fhar = f*(1:8);
for ii = 1:length(fhar)
    ind(ii) = find(fv==fhar(ii));
end
load '30Hz_2V.mat'
figure
subplot(2,2,1)
plot(t(tsind:teind)*1e+3,XODE30(3,tsind:teind)*1e+3,'LineWidth',1.4)
xlabel('Time (ms)')
ylabel('Displacement (mm)')
title(['f = ',int2str(f),' Hz, Input voltage = 2V'])
grid minor
set(gca,'FontSize',16)
subplot(2,2,3)
stem(fv(ind),ODEdB_30(ind),'r','LineWidth',1.4)
xlabel('Frequency (Hz)')
ylabel('SPL (dB re 20uPa)')
title('Pressure spectrum')
grid minor
set(gca,'FontSize',16)
legend(['THD(%) =',num2str(THDODE_30)])

load '30Hz_8V.mat'
subplot(2,2,2)
plot(t(tsind:teind)*1e+3,XODE30(3,tsind:teind)*1e+3,'LineWidth',1.4)
xlabel('Time (ms)')
ylabel('Displacement (mm)')
title(['f = ',int2str(f),' Hz, Input voltage = 8V'])
grid minor
set(gca,'FontSize',16)
subplot(2,2,4)
stem(fv(ind),ODEdB_30(ind),'r','LineWidth',1.4)
xlabel('Frequency (Hz)')
ylabel('SPL (dB re 20uPa)')
title('Pressure spectrum')
grid minor
set(gca,'FontSize',16)
legend(['THD(%) =',num2str(THDODE_30)])
clear