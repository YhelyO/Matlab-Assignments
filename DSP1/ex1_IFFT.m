function [x] = ex1_IFFT(X)
%we used our fft function and the fact that it is a near identical
%algorithem (just need to divide by length(X))
    x = conj(ex1_FFT(conj(X)))/length(X);
end