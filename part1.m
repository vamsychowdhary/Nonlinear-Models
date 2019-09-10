close all
clear
clc

importData;

Re = [Low.El.Re High.El.Re Hot.El.Re];
Le = [Low.El.Le High.El.Le Hot.El.Le];
Bl = [Low.Me.Bl High.Me.Bl Hot.Me.Bl];
Mms = [Low.Me.Mms High.Me.Mms Hot.Me.Mms];
Cms = [Low.Me.Cms High.Me.Cms Hot.Me.Cms];
Rms = [Low.Me.Rms High.Me.Rms Hot.Me.Rms];
clear High Hot Low Ni4 Ni8

% Create matrices
for i=1:3
    F(:,:,i) = [-Re(i)/Le(i) 0 -Bl(i)/Le(i);
                0 0 1;
                Bl(i)/Mms(i) -1/(Cms(i)*Mms(i)) -Rms(i)/Mms(i)];
    G(:,:,i) = [1/Le(i);0;0];
end
%% Calculate the state vector

fs = 96000;                      % Sampling frequency
s = 4;                           % length of signal in seconds
N = s*fs;                        % number of samples
t = 0:1/fs:(N-1)/fs;             % time vector
A = 0.15;                        % Amplitude .15 1.5 .15
nn = 3;                          % 1  2   3
Fv = unique(ceil(Fv));

XFE = zeros(3,N);            % State vector
XMP = zeros(3,N);            % State vector

eg = sum(A*sin(2*pi*Fv*t)); % Input voltage
 
% Forward Euler
for ii = 1:N-1 % Calculate for each sample
    XFE(:,ii+1) = XFE(:,ii)+(1/fs)*(F(:,:,nn)*XFE(:,ii)+G(:,:,nn)*eg(ii));
end
    
% Midpoint Method
for ii = 1:N-1 % Calculate for each sample
    Xtemp = XMP(:,ii)+(1/(2*fs))*(F(:,:,nn)*XMP(:,ii)+G(:,:,nn)*eg(ii));
    XMP(:,ii+1) = XMP(:,ii)+(1/fs)*(F(:,:,nn)*Xtemp+G(:,:,nn)*eg(ii));
end

%% Randomly chosen stable period after 2 secs
clear XFEi XFEx XFEu XMPi XMPx XMPu

XFEi = XFE(1,2*fs+1:end);
XFEx = XFE(2,2*fs+1:end);
XFEu = XFE(3,2*fs+1:end);
 
XMPi = XMP(1,2*fs+1:end);
XMPx = XMP(2,2*fs+1:end);
XMPu = XMP(3,2*fs+1:end);
%% Calculate the spectrum

et = eg(2*fs+1:end);
    
NN = length(XFEi);
w = hann(NN)';
wsum = sum(w);

% Calculate the spectra for FE state vector
Vf = 2*abs(fft(et.*w)/wsum);
IFEf = 2*abs(fft(XFEi.*w)/wsum);
XFEf = 2*abs(fft(XFEx.*w)/wsum);
UFEf = 2*abs(fft(XFEu.*w)/wsum);
HFEf = (XFEf./Vf)*1e+3;

% Calculate the spectra for FE state vector
IMPf = 2*abs(fft(XMPi.*w)/wsum);
XMPf = 2*abs(fft(XMPx.*w)/wsum);
UMPf = 2*abs(fft(XMPu.*w)/wsum);
HMPf = (XMPf./Vf)*1e+3;

fx = (0:NN/2-1)*fs/NN;
Nfreq = find(fx>=5000,1);
fx = fx(1:Nfreq);

% Find the indices of only the signal frequencies
for ii = 1:length(Fv)
    ind(ii) = find(fx==Fv(ii));
end


% Load the measured TF
load('LPM_Hot.mat')
TFmeas = LPMhot_Hx(:,2);
TFfit = LPMhot_Hx(:,4);
fmeas = LPMhot_Hx(:,1);

% Compare the Transfer function from Forward Euler and Mid point
figure
semilogx(fmeas,TFfit,'b','LineWidth',1.1,'LineStyle','-')
hold on
semilogx(fx(ind),HFEf(ind),'r','LineWidth',1.1,'LineStyle','--')
hold on
semilogx(fx(ind),HMPf(ind),'k','LineWidth',1.1,'LineStyle','-.')
title('Transfer function H(f) = X(f)/V(f)')
xlabel('Frequency (Hz)')
ylabel('Magnitude (mm/V)')
xlim([5 3000])
grid minor
set(gca,'FontSize',16)
set(gcf,'position',[50 50 800 600]);
legend('Klippel Fit','Forward Euler','Mid Point')