function [X] = ex1_FFT(x)
%ex1_FFT algorithem 
    g = x(1:2:end);
    h = x(2:2:end);
    N = length(x);
    k = 0:N-1;
    L = N/2;
    X = 0:1:N-1;
    l = 0:1:L-1;
    W_lk_L = exp(-1j*2*pi* l'*k/L);
    W_k_N = exp(-1j*2*pi*k/N);
    %   we used matrix multiplication to compute the algorithem
    for k = 1:1:N
        X(:,k) = g*W_lk_L(:,k) + W_k_N(:,k) * h * W_lk_L(:,k) ;
    end
end