clear all
close all
rng(315468363)
beta=0.5;
n=200; 
span=6;
sps=8;
M=4;
Ts=1/10^6;
t= 0:Ts/8:n*Ts-Ts/8 ;
Rs=1*10^6;

%% Create 200 Symbols
data= randi([0 M-1],n,1);
b=sqrt(2).*pskmod(data,M,pi/4);
b_padded=upsample(b,sps);

figure
plot(t,real(b_padded))
ylim([-1.5,1.5])
title('Re\{b[n]\}')
figure
plot(t,imag(b_padded))
ylim([-1.5,1.5])
title('Im\{b[n]\}')


%% Create Transmission Filter: g_TX=RRC
g_TX=rcosdesign(beta,span,sps);
t_TX=((-length(g_TX)-1)/2+1)*(Ts/sps):Ts/sps:((length(g_TX)-1)/2)*(Ts/sps);

figure
plot(t_TX,g_TX)
title('gTX')
xlabel('t')


%% Convolotion Of Symbols (b) and g_TX
s=conv(b_padded,g_TX,'same');
fs=1/(Ts/8);%sample time 8 mhz
w=linspace(-fs/2,fs/2,length(s));%frequency axis
t_conv= 0:Ts/8:(length(s)-1)*Ts/8 ;
fmax=2*10^6; %f max for bw
sf=fft(s);

figure 
subplot(2,1,1)
plot(t_conv, real(s))
axis('tight')
title('Re\{s[n]\}')
xlabel('t [sec]')
subplot(2,1,2)
plot(t_conv, imag(s))
axis('tight')
title('Im\{s[n]\}')
xlabel('t [sec]')

figure
plot(w,fftshift(sf))
title('s_b_a_s_e_b_a_n_d(f)')
xlabel('frequency [Hz]')


%% Up-Conversion
fc=2.125*10^6;
sp=real(sqrt(2).*s.'.*exp(j*2*pi*fc.*t_conv));
f_pb=linspace(-fc-fmax,fc+fmax,length(sp));

figure
plot(t_conv,sp)
title('s_p_a_s_s_b_a_n_d[n]')
xlabel('t [sec]')
axis('tight')
figure
sf_up=fft(sp);
plot(f_pb,fftshift(sf_up))
title('s_p_a_s_s_b_a_n_d(f)')
xlabel('frequency [Hz]')
axis('tight')


%% Down-Conversion
sc_before_lp=sqrt(2).*cos(2*pi*fc.*t_conv).*sp;
ss_before_lp=-sqrt(2).*sin(2*pi*fc.*t_conv).*sp;
lp=fir1(100*sps, Rs/(fs/2),'low');
sc_after_lp=conv(sc_before_lp,lp,'same');
ss_after_lp=conv(ss_before_lp,lp,'same');
r_n=sc_after_lp+j*ss_after_lp;
w=linspace(-fs/2,fs/2,length(r_n));%frequency axis

figure
subplot(2,1,1)
plot(real(r_n))% two parts of signal after down conversion
title('Re\{r[n]\}');
xlabel('n')
xlim([0,length(r_n)])
subplot(2,1,2)
plot(imag(r_n));
title('Im\{r[n]\}');
xlabel('n')
xlim([0,length(r_n)])
r_f=fft(r_n);
figure
plot(w,fftshift(r_f));
title('Down-Conversion of r(f)');
xlabel('Frequency [Hz]')


%compare with Raised Cosine
Filter=rcosdesign(beta,span,sps,'normal');
r_hat=conv(b_padded,Filter,'same');



index=(length(r_n)-length(r_hat))/2;
figure
hold on
plot(real(r_hat));
plot(real(r_n(1:length(r_n))));
title('Re\{r[n]\}')
axis('tight')
legend('RC filter','down coversion of r[n]')
xlabel('n')
figure
hold on
plot(imag(r_hat));
plot(imag(r_n(1:length(r_n))));
title('Im\{r[n]\}')
axis('tight')
legend('RC filter','down coversion of r[n]')
xlabel('n')
figure 
hold on
title('RC filter, sps=3');
stem(-6:1:6,rcosdesign(0,4,3));


%% Create Matching Filter Without Channel Effect (alpha=0)
p_match_filter=g_TX;%(RRC is symetry and real number), sqrt 2 is not needed
z_n=conv(r_n,p_match_filter,'same');
z_n=z_n(1:sps:end);
index=length(z_n)/2;
z1=z_n(index-n/2+1:end-index+n/2);

figure
hold on
title('Real part of r(n) and z(n)');
plot(real(z1));
plot(real(b));
legend('Re\{z[n]\}','Re\{b[n]\}')
xlabel('n')
ylim([-1.25,1.25])

figure
hold on
title('Image part of r(n) and z(n)');
plot(imag(z1));
plot(imag(b));
legend('Im\{z[n]\}','Im\{b[n]\}')
xlabel('n')
ylim([-1.25,1.25])
     


%% Create Noise for SNR=1
N0=sum(g_TX.^2).*Ts/sps;%integral on p^2 goes to numeric integral with step size TS
sigma=N0*fs;
w_n=sigma/2*(randn(size(sp))+j.*randn(size(sp)));

figure
subplot(2,1,1)
plot(real(w_n));
title('Real\{z[n]\}')
xlabel('n')
axis('tight')
subplot(2,1,2)
plot(imag(w_n));
title('Image\{z[n]\}')
xlabel('n')
axis('tight')


%% Add Noise To sp And Do Down-Conversion
y_n_only_noise=sp+w_n;
yc_only_noise_before_lp=sqrt(2).*cos(2*pi*fc.*t_conv).*y_n_only_noise;
ys_only_noise_before_lp=-sqrt(2).*sin(2*pi*fc.*t_conv).*y_n_only_noise;
lp=fir1(100*sps, Rs/(fs/2),'low');
yc_only_noise_after_lp=conv(yc_only_noise_before_lp,lp,'same');
ys_only_noise_after_lp=conv(ys_only_noise_before_lp,lp,'same');
r_n_only_noise=yc_only_noise_after_lp+j*ys_only_noise_after_lp;
r_f_only_noise=fft(r_n_only_noise);

figure
subplot(2,1,1)
plot(real(y_n_only_noise))
title('Re\{y[n]\}+ Gaussian Noise')
axis('tight')
subplot(2,1,2)
plot(imag(y_n_only_noise))
title('Ima\{y[n]\}+ Gaussian Noise')
axis('tight')

figure
plot(w, fftshift(r_f_only_noise));
title('r(f) + Gaussian Noise and Down-Coversion')
xlabel('Frequency [Hz]')

%% Pass r[n] With Noise in Reciever(Matching Filter)
z_n_only_noise=conv(r_n_only_noise,p_match_filter,'same');
z_n_only_noise=z_n_only_noise(1:sps:end);

index=length(z_n_only_noise)/2;
z_only_noise=z_n_only_noise(index-n/2+1:end-index+n/2);
figure
subplot(2,1,1)
plot(real(z_only_noise));
title('Re\{z[n]\} + Gaussian Noise');
xlabel('n')

subplot(2,1,2)
plot(imag(z_only_noise));
title('Im\{z[n]\} + Gaussian Noise');
xlabel('n')
    


%% Add Channel Effect, alpha=-0.5
alpha=-0.5;
sp_with_channel=[sp,zeros(1,sps)]+alpha*[zeros(1,sps),sp];
sp_with_channel=sp_with_channel(sps/2:end-sps/2-1);
y_n=sp_with_channel+w_n;

% Up-Conversion of y to r
yc_before_lp=sqrt(2).*cos(2*pi*fc.*t_conv).*y_n;
ys_before_lp=-sqrt(2).*sin(2*pi*fc.*t_conv).*y_n;
lp=fir1(100*sps, Rs/(fs/2),'low');
yc_after_lp=conv(yc_before_lp,lp,'same');
ys_after_lp=conv(ys_before_lp,lp,'same');
r_n_channel=yc_after_lp+j*ys_after_lp;
r_f_channel=fft(r_n_channel);

%p_match_filter_channel is p*(-t)
p_match_filter_channel=sqrt(2).*([zeros(1,sps),flip(g_TX)]+alpha.*exp(j*2*pi*fc*Ts).*[flip(g_TX),zeros(1,sps)]);
p_match_filter_channel=p_match_filter_channel(sps+1:end);

z_n_channel=conv(r_n_channel,p_match_filter_channel);
z_n_channel=z_n_channel(1:8:end);

index=length(z_n_channel)/2;
z=z_n_channel(index-n/2+1:end-index+n/2);


figure
plot(t_TX,real(p_match_filter_channel))
title("g_T_X(t) , \alpha=0.5")
xlabel('t [second]')  
figure
subplot(2,1,1)
plot(real(r_n_channel))% two parts of signal after down conversion
title('Re\{r[n]\}  , \alpha=0.5');
xlabel('n')
subplot(2,1,2)
plot(imag(r_n_channel));
title('Image\{r[n]\}  , \alpha=0.5');
xlabel('n')

figure
w=linspace(-fs/2,fs/2,length(r_f_channel));%frequency axis
plot(w,fftshift(real(r_f_channel)));
title('r(f)  , \alpha=0.5, Down-Conversion');

figure
index=length(z_n_channel)/2;
z=z_n_channel(index-n/2:end-index+n/2);
subplot(2,1,1)
plot(real(b));
title('Re\{b[n]\}')
ylim([-1.5,1.5])
xlabel('n')
subplot(2,1,2)
plot(real(z));
title('Re\{z[n]\}')
axis('tight');

figure
subplot(2,1,1)
plot(imag(b));
title('Im\{b[n]\}')
ylim([-1.5,1.5])
xlabel('n')
subplot(2,1,2)
plot(imag(z));
title('Im\{z[n]\}')
axis('tight');


%% Calculation Of h[m]
g_channel_pb=zeros(1,length(t_TX));
g_channel_pb(ceil(length(t_TX)/2))=1;
g_channel_pb(ceil(length(t_TX)/2)+sps)=alpha;
gc_before_lp=sqrt(2).*cos(2*pi*fc.*t_TX).*g_channel_pb;
gs_before_lp=-sqrt(2).*sin(2*pi*fc.*t_TX).*g_channel_pb;
lp=fir1(100*sps, Rs/(fs/2),'low');
gc_after_lp=conv(gc_before_lp,lp,'same');
gs_after_lp=conv(gs_before_lp,lp,'same');
g_bb=gc_after_lp+j*gs_after_lp;
p_t=conv(g_bb,g_TX);
h_m=conv(p_match_filter_channel,p_t);
h_m=h_m(1:sps:end);
index=length(h_m)/2;
h_m=h_m(index-3:index+3);

figure
subplot(2,1,1)
plot(real(g_bb))% two parts of signal after down conversion
title('Re\{g_c[n]\}, Baseband');
xlabel('n')
axis('tight');
subplot(2,1,2);
plot(imag(g_bb));% two parts of signal after down conversion
title('Im\{g_c[n]\}, Baseband');
xlabel('n')
axis('tight');

figure
stem([-3:3],h_m)
title('h[m]');
axis('tight');
xlabel('n')


%% MLSE
h1=h_m(ceil(length(h_m)/2)-1);
lambda=zeros(M,1);
weights=zeros(M,n);
index=zeros(M,n);
symbols=[1+1i,-1+1i,-1-1i,1-1i];
symbol0=1+1i;
estimation=zeros(1,n+1);

for j=1:M
   lambda(j)=real(conj(symbols(j))*b(1))-real(conj(symbols(j))*symbol0*h1); %calculate first 4 forward arrows (matrics)
end
weights(:,1)=lambda;
index(:,1)=1;       
for m=2:n
    for j=1:M     %to which symbol we going to (in this matric calculation,@m) ====b[n] in the formula
       for i=1:M     %from which symbol we came from (in this matric calculation,@m) ===b[n-1] in the formula
           lambda(i)=real(conj(symbols(j))*b(m))-real(conj(symbols(j))*symbols(i)*h1); %all together this 2 for loops 
		    %                    b*[n]      /\               b*[n]       b[n-1]          calculate the 16 forward arrows (matrics)
       end  %                               ||                                          in the m time point
	   %                                    ||
	   %                         the real symbol that we recive at this time point m ===z[m] in the formula
       final=weights(:,m-1)+lambda; %sum up all the options we currently left with 
       weights(j,m)=max(final);  %decide which i (symbol we came from) is most likely 
	                             %(given we going to symbol j, from which symbol i is the most likely that we came from)
       [~,max_index]=max(final); %saving which symbol it was (his index i in the symbol(i) vector)
       index(j,m)=max_index;     %saving it here to all iterations 
	   %(this is the most important result those are our OUR DECISION/ESTIMATE SYMBOL)
	   % in the index matrix we save: given the arrow is to symbol j (the next symbol), in time point m, 
	   %so in this case the most likely symbol  that we came from (OUR DECISION/ESTIMATE SYMBOL) is i (j,m)=i
	   end 
end

for i=1:M  %M=4 . n=200
   if (round(b(n))==round(symbols(i))) %checking what is the last symbol (n=200) that was transmitted
      last_max_index=i;  % this is the i of this symbol
   end
end
for m=n+1:-1:1 % going down from 200
   estimation(m)=symbols(last_max_index); %here we set the estimate symbol,
                                  % first one (n=200) is easy cuse we already found it in the for loop
   if(m~=1) %just before the loop on m is end
       last_max_index=index(last_max_index,m-1);
	    %here we set the i (to choose from the symbol(i) 4X1 vector of the QPSK symbols) 
	    %and we do it by the index matrix (that we already said it is the most important result)
	    %and we giving her the j of forward symbole that she need 
		%(that we already know because we going *down* from 200 that we found already
		%and the time point m that she also need
		%and she (the index matrix) giving us back the most likely symbol 
		%that was transmited one step before (===i), 
		%by given the symbol transmitted (or estimated) now (===j)
		%because this is just how we defined her
		%note: we just called the matrix: 'she' because we are very respectfal :)
   end
end

BER =0;   
b=[1+1i,b.'];
for i=1:n+1
   if(round(real(estimation(i)))==round(real(b(i))))
       BER=BER+1;
   end
   if(round(imag(estimation(i)))==round(imag(b(i))))
       BER=BER+1;
   end   
end
correct=BER;
BER=BER/(2*(n+1))*100;
BER

figure
subplot(2,1,1)
plot(real(estimation))
axis('tight')
ylim([-1.5,1.5])
title('Real Part Estimation')
subplot(2,1,2)
plot(real(b))
axis('tight')
ylim([-1.5,1.5])
title('Real\{[b]\}')

figure
subplot(2,1,1)
plot(imag(estimation))
axis('tight')
ylim([-1.5,1.5])
title('Imaginry Part Estimation')
subplot(2,1,2)
plot(imag(b))
axis('tight')
ylim([-1.5,1.5])
title('Image\{[b]\}')

