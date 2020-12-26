function X = dtft(x, w)
%dtft  Calculates X(e^jw) using inputs impulse response vector and
%frequency vector
X = zeros(1,length(w));
for k = 1:length(w)
    for r = 1:length(x)
        X(k)= X(k)+(x(r)*exp(-1j*w(k)*(r-1)));
    end
end
end

