%% Initialize
clc
close all

fs = 96000;                       % sampling freq
N = 4*fs;                         % number of samples
f = 121;                           % Signal frequency
A = [2 8];             % Amplitude
lhh = 1;                          % Low High Hot

% Load the generated coeffficients
load ('Cms_Bl_right_coeffs_8th.mat')

%% FE method
% for ii = 1:2
%     [XFE(:,:,ii),eigvalFE(ii,:)] = part2b(fs,N,f,A(ii),"FE",lhh,Cms_coeff,Bl_coeff,Le_coeff);
% end
%% MP method
for ii = 1:2
    [XMP(:,:,ii),eigvalMP(ii,:)] = part2b(fs,N,f,A(ii),"MP",lhh,Cms_coeff,Bl_coeff,Le_coeff);
end

%% calculate THD for FE
% NDFT = length(XFE(4,2*fs+1:end,1));
% w = hann(NDFT)';
% wsum = sum(w);
% fv = (0:NDFT/2-1)*fs/NDFT;
% for ii = 1:6
%     XF = fft(XFE(4,2*fs+1:end,ii))/wsum;
%     [FErms_30(ii,:),FEdB_30(ii,:)] = vel2pres(XF(1:NDFT/2),fv);
%     [THDFE_30(ii),~] = Task_3(FErms_30(ii,:),f,fs,NDFT,10,1);
%     dispFE_30(ii,:) = [min(XFE(3,2*fs+1:end,ii)) max(XFE(3,2*fs+1:end,ii))];
% end

%% calculate THD for MP
NDFT = length(XMP(4,2*fs+1:end,1));
w = hann(NDFT)';
wsum = sum(w);
fv = (0:NDFT/2-1)*fs/NDFT;
for ii = 1:2
    XF = fft(XMP(4,2*fs+1:end,ii))/wsum;
    [MPrms_30(ii,:),MPdB_30(ii,:)] = vel2pres(XF(1:NDFT/2),fv);
    [THDMP_30(ii),~] = Task_3(MPrms_30(ii,:),f,fs,NDFT,10,1);
    dispMP_30(ii,:) = [min(XMP(3,2*fs+1:end,ii)) max(XMP(3,2*fs+1:end,ii))];
end

%%
t = 0:1/fs:(N-1)/fs;            % time vector
tend = 5/f + 2;
tsind = find(t>=2,1);
teind = find(t>=tend,1);
fhar = f*(1:8);
for ii = 1:length(fhar)
    ind(ii) = find(fv==fhar(ii));
end
figure
subplot(2,2,1)
plot(t(tsind:teind)*1e+3,XMP(3,tsind:teind,1)*1e+3,'LineWidth',1.4)
xlabel('Time (ms)')
ylabel('Displacement (mm)')
title(['f = ',int2str(f),' Hz, Input voltage = 2V'])
grid minor
set(gca,'FontSize',16)
subplot(2,2,3)
stem(fv(ind),MPdB_30(1,ind),'r','LineWidth',1.4)
xlabel('Frequency (Hz)')
ylabel('SPL (dB re 20uPa)')
title('Pressure spectrum')
grid minor
set(gca,'FontSize',16)
legend(['THD(%) =',num2str(THDMP_30(1))])

subplot(2,2,2)
plot(t(tsind:teind)*1e+3,XMP(3,tsind:teind,2)*1e+3,'LineWidth',1.4)
xlabel('Time (ms)')
ylabel('Displacement (mm)')
title(['f = ',int2str(f),' Hz, Input voltage = 8V'])
grid minor
set(gca,'FontSize',16)
subplot(2,2,4)
stem(fv(ind),MPdB_30(2,ind),'r','LineWidth',1.4)
xlabel('Frequency (Hz)')
ylabel('SPL (dB re 20uPa)')
title('Pressure spectrum')
grid minor
set(gca,'FontSize',16)
legend(['THD(%) =',num2str(THDMP_30(2))])