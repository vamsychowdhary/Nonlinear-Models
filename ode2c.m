function dydt = ode2c(t,X)

global Re R2 L2 Mms Rms Bl8 Le8 Cms8 ordCms ordBl ordLe CmsLin BlLin LeLin

f = 30;
A = 8;
lhh = 1;
eg = A*sin(2*pi*f*t);           % Input voltage
F = zeros(4,4);                 % Initialize F


% Update F and G matrices

    xn = X(3); % Dispalcement at the present time sample
    in = X(1); % current at the present time sample
    iL2 = X(2); % current at L2 in the present time sample

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
        
    FX1 = F(1,1)*X(1)+F(1,2)*X(2)+F(1,4)*X(4);
    FX2 = F(2,1)*X(1)+F(2,2)*X(2)+F(2,4)*X(4);                                        
    FX3 = X(4);
    FX4 = F(4,1)*X(1)+F(4,3)*X(3)+F(4,4)*X(4);
    FX = [FX1; FX2; FX3; FX4];
    
    GE = G.*eg;
    
    dydt = FX + GE;

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