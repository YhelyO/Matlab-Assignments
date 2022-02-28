function [F] = ex1(f)
%EX1 Summary of this function goes here
%   Detailed explanation goes here
F_rec = ex1_FFT_recursive(f);
F = ex1_FFT(f);
IF_rec = ex1_IFFT_recursive(F_rec);
IF = ex1_IFFT(F);
x = 1:1:length(F);
%FFT plot
figure
subplot(4,1,1)
plot(x,real(fft(f)));
title('fft')
axis tight;

subplot(4,1,2)
plot(x,F)
title('my fft')
axis tight;

subplot(4,1,3)
plot(x,real(F_rec))
title('my rec fft')
axis tight;

subplot(4,1,4)
plot(x,real(fft(f) - F))
legend({'fft- my rec fft '},'Location','southwest')
title('\delta plot')
axis tight;

x0=300;
y0=-100;
width=1300;
height=800;
set(gcf,'position',[x0,y0,width,height])

%IFFT plot
figure
subplot(4,1,1)
plot(x,ifft(F))
title('ifft')
axis tight;

subplot(4,1,2)
plot(x,IF)
title('my ifft')
axis tight;

subplot(4,1,3)
plot(x,IF_rec)
title('my rec ifft')
axis tight;

subplot(4,1,4)
plot(x,ifft(F) - IF)
legend({'fft- my rec fft '},'Location','southwest')
title('\delta plot')
axis tight;

x0=300;
y0=-100;
width=1300;
height=800;
set(gcf,'position',[x0,y0,width,height])
end

