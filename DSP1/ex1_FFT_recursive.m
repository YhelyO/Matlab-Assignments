function F = ex1_FFT_recursive(f)
%just to be in sync with the ex1_IFFT_recursive function
F = fft_rec(f) ;
end

function F = fft_rec(f)
n = length(f);
if (n == 1) %if the length is 1 return the function
  F = f;
else
  %divide into two vectors, evens and odds
  f_even = f(1:2:n);
  f_odd = f(2:2:n);
  %recursive calling, in a tree style 
  X1 = fft_rec(f_even);
  X2 = fft_rec(f_odd).*Wn(n);  
  F1 = X1 + X2;
  F2 = X1 - X2;
  F = [F1 F2];
end
end

function W = Wn(n)
%just to be more neat
m = n/2;
W = exp(-2*pi*1i.*(0:1:m-1)/n);
end
