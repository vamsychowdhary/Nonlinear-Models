addpath('MeasData')

%load 'Fv.mat'

[~,~,ExLinearLow] = xlsread('LinearParametersLow.xlsx');
[~,~,ExLinearHigh] = xlsread('LinearParametersHigh.xlsx');
[~,~,ExLinearHot] = xlsread('LinearParametersHot.xlsx');
[~,~,ExNonlinear4] = xlsread('NonlinearParameters4.xlsx');
[~,~,ExNonlinear8] = xlsread('NonlinearParameters8.xlsx');

% Structure linear low data
Low.El.Re = cell2mat(ExLinearLow(2,2));
Low.El.Le = cell2mat(ExLinearLow(3,2))*1e-3;
Low.El.L2 = cell2mat(ExLinearLow(4,2))*1e-3;
Low.El.R2 = cell2mat(ExLinearLow(5,2));
Low.El.Cmes = cell2mat(ExLinearLow(6,2))*1e-6;
Low.El.Lces = cell2mat(ExLinearLow(7,2))*1e-3;
Low.El.Res = cell2mat(ExLinearLow(8,2));
Low.El.fs = cell2mat(ExLinearLow(9,2));

Low.Me.Mms = cell2mat(ExLinearLow(13,2))*1e-3;
Low.Me.MmsSd = cell2mat(ExLinearLow(14,2))*1e-3;
Low.Me.Rms = cell2mat(ExLinearLow(15,2));
Low.Me.Cms = cell2mat(ExLinearLow(16,2))*1e-3;
Low.Me.Kms = cell2mat(ExLinearLow(17,2))*1e+3;
Low.Me.Bl = cell2mat(ExLinearLow(18,2));
Low.Me.Lambda = cell2mat(ExLinearLow(19,2));

Low.Loss.Qtp = cell2mat(ExLinearLow(22,2));
Low.Loss.Qms = cell2mat(ExLinearLow(23,2));
Low.Loss.Qes = cell2mat(ExLinearLow(24,2));
Low.Loss.Qts = cell2mat(ExLinearLow(25,2));

Low.Other.Vas = cell2mat(ExLinearLow(28,2));
Low.Other.n0 = cell2mat(ExLinearLow(29,2));
Low.Other.Lm = cell2mat(ExLinearLow(30,2));
Low.Other.Lnom = cell2mat(ExLinearLow(31,2));
Low.Other.Zrms = cell2mat(ExLinearLow(33,2));
Low.Other.Hxrms = cell2mat(ExLinearLow(34,2));
Low.Other.Rseries = cell2mat(ExLinearLow(36,2));
Low.Other.Sd = cell2mat(ExLinearLow(37,2))*1e-4;

% Structure linear high data
High.El.Re = cell2mat(ExLinearHigh(2,2));
High.El.Le = cell2mat(ExLinearHigh(3,2))*10^-3;
High.El.L2 = cell2mat(ExLinearHigh(4,2))*10^-3;
High.El.R2 = cell2mat(ExLinearHigh(5,2));
High.El.Cmes = cell2mat(ExLinearHigh(6,2))*10^-6;
High.El.Lces = cell2mat(ExLinearHigh(7,2))*10^-3;
High.El.Res = cell2mat(ExLinearHigh(8,2));
High.El.fs = cell2mat(ExLinearHigh(9,2));

High.Me.Mms = cell2mat(ExLinearHigh(13,2))*10^-3;
High.Me.MmsSd = cell2mat(ExLinearHigh(14,2))*10^-3;
High.Me.Rms = cell2mat(ExLinearHigh(15,2));
High.Me.Cms = cell2mat(ExLinearHigh(16,2))*10^-3;
High.Me.Kms = cell2mat(ExLinearHigh(17,2))*10^3;
High.Me.Bl = cell2mat(ExLinearHigh(18,2));
High.Me.Lambda = cell2mat(ExLinearHigh(19,2));

High.Loss.Qtp = cell2mat(ExLinearHigh(22,2));
High.Loss.Qms = cell2mat(ExLinearHigh(23,2));
High.Loss.Qes = cell2mat(ExLinearHigh(24,2));
High.Loss.Qts = cell2mat(ExLinearHigh(25,2));

High.Other.Vas = cell2mat(ExLinearHigh(28,2));
High.Other.n0 = cell2mat(ExLinearHigh(29,2));
High.Other.Lm = cell2mat(ExLinearHigh(30,2));
High.Other.Lnom = cell2mat(ExLinearHigh(31,2));
High.Other.Zrms = cell2mat(ExLinearHigh(33,2));
High.Other.Hxrms = cell2mat(ExLinearHigh(34,2));
High.Other.Rseries = cell2mat(ExLinearHigh(36,2));
High.Other.Sd = cell2mat(ExLinearHigh(37,2))*10^-2;

% Structure linear hot data
Hot.El.Re = cell2mat(ExLinearHot(2,2));
Hot.El.Le = cell2mat(ExLinearHot(3,2))*10^-3;
Hot.El.L2 = cell2mat(ExLinearHot(4,2))*10^-3;
Hot.El.R2 = cell2mat(ExLinearHot(5,2));
Hot.El.Cmes = cell2mat(ExLinearHot(6,2))*10^-6;
Hot.El.Lces = cell2mat(ExLinearHot(7,2))*10^-3;
Hot.El.Res = cell2mat(ExLinearHot(8,2));
Hot.El.fs = cell2mat(ExLinearHot(9,2));

Hot.Me.Mms = cell2mat(ExLinearHot(13,2))*10^-3;
Hot.Me.MmsSd = cell2mat(ExLinearHot(14,2))*10^-3;
Hot.Me.Rms = cell2mat(ExLinearHot(15,2));
Hot.Me.Cms = cell2mat(ExLinearHot(16,2))*10^-3;
Hot.Me.Kms = cell2mat(ExLinearHot(17,2))*10^3;
Hot.Me.Bl = cell2mat(ExLinearHot(18,2));
Hot.Me.Lambda = cell2mat(ExLinearHot(19,2));

Hot.Loss.Qtp = cell2mat(ExLinearHot(22,2));
Hot.Loss.Qms = cell2mat(ExLinearHot(23,2));
Hot.Loss.Qes = cell2mat(ExLinearHot(24,2));
Hot.Loss.Qts = cell2mat(ExLinearHot(25,2));

Hot.Other.Vas = cell2mat(ExLinearHot(28,2));
Hot.Other.n0 = cell2mat(ExLinearHot(29,2));
Hot.Other.Lm = cell2mat(ExLinearHot(30,2));
Hot.Other.Lnom = cell2mat(ExLinearHot(31,2));
Hot.Other.Zrms = cell2mat(ExLinearHot(33,2));
Hot.Other.Hxrms = cell2mat(ExLinearHot(34,2));
Hot.Other.Rseries = cell2mat(ExLinearHot(36,2));
Hot.Other.Sd = cell2mat(ExLinearHot(37,2))*10^-2;

Bly = [1 1e3 1e3^2 1e3^3 1e3^4 1e3^5 1e3^6 1e3^7 1e3^8]';
Ly = [1e-3 1 1e3 1e3^2 1e3^3 1e3^4 1e3^5 1e3^6 1e3^7]';
Cy = [1e-3 1 1e3 1e3^2 1e3^3 1e3^4 1e3^5 1e3^6 1e3^7]';

% Structure nonlinear 4 data
Nl4.Bl = cell2mat(ExNonlinear4(24:28,2)).*Bly(1:5);
Nl4.L = cell2mat(ExNonlinear4(30:34,2)).*Ly(1:5);
Nl4.C = cell2mat(ExNonlinear4(36:40,2)).*Cy(1:5);
Nl4.K = cell2mat(ExNonlinear4(42:45,2));
Nl4.f = cell2mat(ExNonlinear4(47:48,2));

% Structure nonlinear 8 data
Nl8.Bl = cell2mat(ExNonlinear8(24:32,2)).*Bly;
Nl8.L = cell2mat(ExNonlinear8(34:42,2)).*Ly;
Nl8.C = cell2mat(ExNonlinear8(44:52,2)).*Cy;
Nl8.f = cell2mat(ExNonlinear8(59:60,2));



clear ExLinearHigh ExLinearLow ExLinearHot ExNonlinear4 ExNonlinear8 Bly Ly Cy