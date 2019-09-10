function [X,eigval] = part2b(fs,N,f,A,meth,lhh,Cms_coeff,Bl_coeff,Le_coeff)

% Load data
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

Re = [Low.El.Re High.El.Re Hot.El.Re];
R2 = [Low.El.R2 High.El.R2 Hot.El.R2];
L2 = [Low.El.L2 High.El.L2 Hot.El.L2];

Mms = [Low.Me.Mms High.Me.Mms Hot.Me.Mms];
Rms = [Low.Me.Rms High.Me.Rms Hot.Me.Rms];

Bl8 = Bl_coeff';
Le8 = Le_coeff';
Cms8 = Cms_coeff';
ordCms = 8;        % Order of the non-linearity
ordBl = 8;        % Order of the non-linearity
ordLe = 8;        % Order of the non-linearity

CmsLin = [Low.Me.Cms High.Me.Cms Hot.Me.Cms];
BlLin = [Low.Me.Bl High.Me.Bl Hot.Me.Bl];
LeLin = [Low.El.Le High.El.Le Hot.El.Le];

t = 0:1/fs:(N-1)/fs;            % time vector
X = zeros(4,N);                 % State vector
eg = A*sin(2*pi*f*t);           % Input voltage
F = zeros(4,4);                 % Initialize F

% Initialize the vector to store eigen vlaues
eigval = zeros(1,N-1);

% Calculate the state vector

% Forward Euler
if (meth == "FE")
    
    for ii = 1:N-1 % Calculate for each sample

        xn = X(3,ii); % Dispalcement at the present time sample
        in = X(1,ii); % current at the present time sample
        iL2 = X(2,ii); % current at L2 in the present time sample

        Le0 = Le(0,ordLe,LeLin,lhh,Le8);
        Lexn = Le(xn,ordLe,LeLin,lhh,Le8);
        dLexn = dLe(xn,ordLe,Le8);
        Blxn = Bl(xn,ordBl,BlLin,lhh,Bl8);
        Cmxn = Cm(xn,ordCms,CmsLin,lhh,Cms8);

        R2xn = (R2(lhh)/Le0)*Lexn;
        L2xn = (L2(lhh)/Le0)*Lexn;

        F(1,1) = -(Re(lhh)+R2xn)/Lexn;
        F(1,2) = R2xn/Lexn;
        F(1,4) = (-1/Lexn)*(in*dLexn+Blxn);
        F(2,1) = R2xn/L2xn;
        F(2,2) = -F(2,1);
        F(2,4) = -(iL2/L2xn)*((L2(lhh)/Le0)*dLexn);
        F(3,4) = 1;
        F(4,1) = Blxn/Mms(lhh);
        F(4,3) = -1/(Cmxn*Mms(lhh));
        F(4,4) = -Rms(lhh)/Mms(lhh);

        G = [1/Lexn;0;0;0];
        
        % calculate the eigen values of F to check they are not > 1
        eigval(ii) = max(abs(eig(eye(4)+F/fs)));
        
        X(:,ii+1) = X(:,ii)+1/fs*(F*X(:,ii)+G*eg(ii));
    end
    
else
    
    % Midpoint Method
    for ii = 1:N-1 % Calculate for each sample

        xn = X(3,ii); % Dispalcement at the present time sample
        in = X(1,ii); % current at the present time sample
        iL2 = X(2,ii); % current at L2 in the present time sample

        Le0 = Le(0,ordLe,LeLin,lhh,Le8);
        Lexn = Le(xn,ordLe,LeLin,lhh,Le8);
        dLexn = dLe(xn,ordLe,Le8);
        Blxn = Bl(xn,ordBl,BlLin,lhh,Bl8);
        Cmxn = Cm(xn,ordCms,CmsLin,lhh,Cms8);

        R2xn = (R2(lhh)/Le0)*Lexn;
        L2xn = (L2(lhh)/Le0)*Lexn;

        F(1,1) = -(Re(lhh)+R2xn)/Lexn;
        F(1,2) = R2xn/Lexn;
        F(1,4) = (-1/Lexn)*(in*dLexn+Blxn);
        F(2,1) = R2xn/L2xn;
        F(2,2) = -F(2,1);
        F(2,4) = -(iL2/L2xn)*((L2(lhh)/Le0)*dLexn);
        F(3,4) = 1;
        F(4,1) = Blxn/Mms(lhh);
        F(4,3) = -1/(Cmxn*Mms(lhh));
        F(4,4) = -Rms(lhh)/Mms(lhh);

        G = [1/Lexn;0;0;0];
        
        % calculate the eigen values of F to check they are not > 1
        eigval(ii) = max(abs(eig(eye(4)+F/fs)));

        Xtemp = X(:,ii)+(1/(2*fs))*(F*X(:,ii)+G*eg(ii));
        X(:,ii+1) = X(:,ii)+(1/fs)*(F*Xtemp+G*eg(ii));
    end

end

% Nonlinear Functions

% Calculate Le for different xn for the same order
function Leout = Le(xn,ordLe,LeLin,lhh,Le8)
    if ordLe == 0
        Leout = LeLin(lhh);
        return
    else
        Lec = Le8(1:ordLe+1);
    end
    Lex = sum(Lec.*(xn.^(0:ordLe))',1);
    Leout = Lex;
end

% Calculate dLe for different xn for the same order
function dLeout = dLe(xn,ordLe,Le8)
    if ordLe == 0
        dLeout = 0;
        return
    else
        Lec = Le8(2:1+ordLe);
    end
    Lex = sum(Lec.*(1:ordLe)'.*(xn.^(0:ordLe-1))',1);
    dLeout = Lex;
end

% Calculate Bl for different xn for the same order
function Blout = Bl(xn,ordBl,BlLin,lhh,Bl8)
    if ordBl == 0
        Blout = BlLin(lhh);
        return
    else
        Blc = Bl8(1:ordBl+1);
    end
    Blout = sum(Blc.*(xn.^(0:ordBl))',1);
end

% Calculate Cm for different xn for the same order
function Cmout = Cm(xn,ordCms,CmsLin,lhh,Cms8)
    if ordCms <=0
        Cmout = CmsLin(lhh);
        return
    else
        Cmsc = Cms8(1:ordCms+1);
    end
    Cmout = sum(Cmsc.*(xn.^(0:ordCms))',1);
end
end