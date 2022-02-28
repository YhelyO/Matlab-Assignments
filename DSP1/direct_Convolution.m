function y = direct_Convolution(x, h)
 %just the straightforward convolution algorithm
    N = length(x);
    p = length(h);
    Nk = N + p -1;
    y = zeros(1,Nk);
    for i = 1:N
          for k = 1:p
           y(i+k-1) = y(i+k-1) + h(k)*x(i);
          end
    end
  
end
