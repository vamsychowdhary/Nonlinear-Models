clc
close all
clear

%% Compliance polynomial functions
xunits = 1e-3; % to account for mm
xdim = linspace(-6,6,200)*xunits;
sig1 = 2.5*xunits; % changes the slope of left side - increases in value decreases slope
sig2 = 2.5*xunits;  % changes the slope of right side - increases in value decreases slope
c1 = 2*xunits;   % c1 > c2 shifts the center to the right
c2 = 1*xunits;
Cms = xunits*gauss2mf(xdim,[sig1 c1 sig2 c2]);

Cms_coeff = fliplr(polyfit(xdim,Cms,8));

%% Force factor polynomial functions
sig1 = 2.5*xunits; % changes the slope of left side - increases in value decreases slope
sig2 = 2.5*xunits;  % changes the slope of right side - increases in value decreases slope
c1 = 2*xunits;   % c1 > c2 shifts the center to the right
c2 = 1*xunits;
Bl = 7*(gauss2mf(xdim,[sig1 c1 sig2 c2]));

Bl_coeff = fliplr(polyfit(xdim,Bl,8));
%% Inductance polynomial functions
a = -1*1e+3; % 
c = 0*1e-3; %
Le = 1e-3*(0.4*sigmf(xdim,[a c])+0.4);

Le_coeff = fliplr(polyfit(xdim,Le,8));
%% Plot

figure
subplot(3,1,1)
plot(xdim*1e+3,Cms*1e+3,'color','r','LineWidth',1.2)
ylabel('Compliance (mm/N)')
set(gca,'FontSize',14)
subplot(3,1,2)
plot(xdim*1e+3,Bl,'color','b','LineWidth',1.2)
ylabel('Bl (Tm)')
set(gca,'FontSize',14)
subplot(3,1,3)
plot(xdim*1e+3,Le*1e+3,'color','k','LineWidth',1.2)
ylabel('L_{e} (mH)')
xlabel('Excursion (mm)')
set(gca,'FontSize',14)
set(gcf,'position',[50 50 800 600])