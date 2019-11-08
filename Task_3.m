function [THDr,THDf] = thd(Y,f0,fs,N,nmax,Nsl)

% Frequency axis
f = (0:N/2-1)*fs/N;

% compute the numerator for THDr and THDf
num = 0;
for i = 2:nmax
    fh = i*f0;
    % Break the loop if order exceeds maximum frequency available
    if (fh > fs/2)
        break;
    end
    [~,ind] = sort(abs(f-fh));
    num = num + sum(Y(ind(1:Nsl)).^2);
end

% compute the denominator for THDr
denr = 0;
for i = 1:nmax
    fh = i*f0;
    % Break the loop if order exceeds maximum frequency available    
    if (fh > fs/2)
        break;
    end
    [~,ind] = sort(abs(f-fh));
    denr = denr + sum(Y(ind(1:Nsl)).^2);
end

% compute the denominator for THDf
[~,ind] = sort(abs(f-f0));
denf = sum(Y(ind(1:Nsl)).^2);

% Compute THDr
THDr = 100*sqrt(num/denr);

% Compute THDf
THDf = 100*sqrt(num/denf);
end
