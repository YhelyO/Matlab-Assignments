function [] = ex3()
%EX3 Summary of this function goes here
%   Detailed explanation goes here
load('sig_x.mat');
figure
plot(-30:61/274499:30,fftshift(fft(x)))
load('filter_1.mat');
F1 = xx;
load('filter_2.mat');
F2 = xx;
figure
subplot(2,1,1)
plot(-30:1:30,abs(fftshift(fft(F1))));
title('fft of filter 1');
xlabel('freq');
ylabel('amplitude');
axis tight;
subplot(2,1,2)
plot(-30:1:30,abs(fftshift(fft(F2))));
title('fft of filter 2');
xlabel('freq');
ylabel('amplitude');
axis tight;

j=1;
for K = 62:1000:270000
    tic
    ova_conv1 = OVA(x, F1, K); 
    time(j) = toc;
    j=j+1;
end
 
[M,I] = min(time);
plot(time)
 
K = 62 + (I-1)*1000;

d_conv1 = direct_Convolution(x,F1);
d_conv1 = d_conv1(1:1000);

mtlb_conv1 = conv(x,F1);
mtlb_conv1 = mtlb_conv1(1:1000);

ova_conv1 = OVA(x, F1, K);
ova_conv1 = ova_conv1(1:1000);

ovs_conv1 = OVS(x, F1, K);
ovs_conv1 = ovs_conv1(1:1000);

d_conv2 = direct_Convolution(x,F2);
d_conv2 = d_conv2(1:1000);

mtlb_conv2 = conv(x,F2);
mtlb_conv2 = mtlb_conv2(1:1000);

ova_conv2 = OVA(x, F2, K);
ova_conv2 = ova_conv2(1:1000);

ovs_conv2 = OVS(x, F2, K);
ovs_conv2 = ovs_conv2(1:1000);


figure
plot(1:100, real(d_conv1(1:100)),1:100, real(mtlb_conv1(1:100)), 1:100, real(ova_conv1(1:100)), 1:100, real(ovs_conv1(1:100)));
legend({'direct','matlab','ova','ovs'},'Location','southwest')


figure
subplot(4,1,1)
plot(real(d_conv1));
title('direct Convolution of x and filter 1');
xlabel('C1[n]');
ylabel('n');
axis tight;
 
subplot(4,1,2)
plot(real(mtlb_conv1));
title('conv of x and filter 1');
xlabel('C1[n]');
ylabel('n');
axis tight;


subplot(4,1,3)
plot(real(ova_conv1));
title('ova conv of x and filter 1');
xlabel('C1[n]');
ylabel('n');
axis tight;

subplot(4,1,4)
plot(real(ovs_conv1));
title('ovs conv of x and filter 1');
xlabel('C1[n]');
ylabel('n');
axis tight;

x0=300;
y0=-100;
width=1300;
height=800;
set(gcf,'position',[x0,y0,width,height])

figure
subplot(4,1,1)
plot(real(d_conv2));
title('direct Convolution of x and filter 2');
xlabel('C2[n]');
ylabel('n');
axis tight;
 
subplot(4,1,2)
plot(real(mtlb_conv2));
title('conv of x and filter 2');
xlabel('C2[n]');
ylabel('n');
axis tight;


subplot(4,1,3)
plot(real(ova_conv2));
title('ova conv of x and filter 2');
xlabel('C2[n]');
ylabel('n');
axis tight;

subplot(4,1,4)
plot(real(ovs_conv2));
title('ovs conv of x and filter 2');
xlabel('C2[n]');
ylabel('n');
axis tight;

x0=300;
y0=-100;
width=1300;
height=800;
set(gcf,'position',[x0,y0,width,height])
end





