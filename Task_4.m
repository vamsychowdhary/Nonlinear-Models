function [TIMDr,TIMDf] = Task_4(Y,f1,f2,fs,N,nmax,Nsl)

% Frequency axis
f = (0:N/2-1)*fs/N;

% compute the harmonics below f2
negf2 = 0;
for i = 1:nmax
    fh = f2-i*f1;
    % Break the loop if the harmonic frequency is less than 0
    if (fh < 0)
        break;
    end
    [~,ind] = sort(abs(f-fh));
    negf2 = negf2 + sum(Y(ind(1:Nsl)).^2);
end

% compute the harmonics above f2
posf2 = 0;
for i = 1:nmax
    fh = f2+i*f1;
    % Break the loop if the harmonic frequency is greater than fs/2
    if (fh > fs/2)
        break;
    end
    [~,ind] = sort(abs(f-fh));
    posf2 = posf2 + sum(Y(ind(1:Nsl)).^2);
end

% compute the value at f2
[~,ind] = sort(abs(f-f2));
atf2 = sum(Y(ind(1:Nsl)).^2);

% Now build the numerator for both THDr and THDf
num = negf2+posf2;

% Denominator for THDr
denr = negf2+atf2+posf2;

% Denominator for THDf
denf = atf2;

% Compute THDr
TIMDr = 100*sqrt(num/denr);

% Compute THDf
TIMDf = 100*sqrt(num/denf);
end