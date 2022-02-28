function [] = ex3_runtime()
load('sig_x.mat');
load('filter_1.mat');
F1 = xx;
load('filter_2.mat');
F2 = xx;
j=1;
for K = 100:1000:270000
    %ova runtime calculation
    tic
    ova_conv1 = OVA(x, F1, K); 
    T = toc;
    time_ova_f1(j) = T ;
    
    tic
    ova_conv1 = OVA(x, F2, K); 
    T = toc;
    time_ova_f2(j) = T ;
    %ovs runtime calculation
    tic
    ovs_conv1 = OVS(x, F1, K); 
    T = toc;
    time_ovs_f1(j) = T ;
    
    tic
    ovs_conv1 = OVS(x, F2, K); 
    T = toc;
    time_ovs_f2(j) = T ;
    
    %direct runtime calculation
    tic 
    direct_conv2 = direct_Convolution(x, F2);
    T = toc;
    time_direct_f2(j) = T;
    
    tic 
    direct_conv1 = direct_Convolution(x, F1);
    T = toc;
    time_direct_f1(j) = T;
    j=j+1;
end


fprintf("min runtime direct f1: %d\n",min(time_direct_f1) )

fprintf("min runtime direct f2: %d\n",min(time_direct_f2) )


K = 100:1000:270000;
%ova plots
[opt_runtime_ova_f1,index_ova_f1] = min(time_ova_f1);
fprintf('optimal runtime ova f1: %d .  K = %d\n',opt_runtime_ova_f1, K(index_ova_f1))

figure
plot(K, time_ova_f1)
title('ova runtime vs K (conv of x and filter 1)');
xlabel('K');
ylabel('runtime');
axis tight;

[opt_runtime_ova_f2,index_ova_f2] = min(time_ova_f2);
fprintf('optimal runtime ova f2: %d .  K = %d\n',opt_runtime_ova_f2, K(index_ova_f2))

figure
plot(K, time_ova_f2)
title('ova runtime vs K (conv of x and filter 2)');
xlabel('K');
ylabel('runtime');
axis tight;


%ovs plots
[opt_runtime_ovs_f1,index_ovs_f1] = min(time_ovs_f1);
fprintf('optimal runtime ovs f1: %d .  K = %d\n',opt_runtime_ovs_f1, K(index_ovs_f1))

figure
plot(K, time_ovs_f1)
title('ovs runtime vs K (conv of x and filter 1)');
xlabel('K');
ylabel('runtime');
axis tight;

[opt_runtime_ovs_f2,index_ovs_f2] = min(time_ovs_f2);
fprintf('optimal runtime ovs f2: %d .  K = %d\n',opt_runtime_ovs_f2, K(index_ovs_f2))

figure
plot(K, time_ovs_f2)
title('ovs runtime vs K (conv of x and filter 2)');
xlabel('K');
ylabel('runtime');
axis tight;

%all together plot
figure
subplot(2,1,1)
plot(K, time_ova_f1, K, time_ovs_f1)
title('ova and ovs runtime vs K (conv of x and filter 1)');
xlabel('K');
ylabel('runtime');
legend({'ova','ovs'},'Location','southwest')
axis tight;

subplot(2,1,2)
plot(K, time_ova_f2, K, time_ovs_f2)
title('ova and ovs runtime vs K (conv of x and filter 2)');
xlabel('K');
ylabel('runtime');
legend({'ova','ovs'},'Location','southwest')
axis tight;

end

