clear all
close all
rng(315468363)
beta=0.6;
n=2000; %n= number of symbols, should be 200.
T=1;
Fs=100; % num of samples per T
Ts=T/100;
span=4;%20*sps;%6;
M=4; %because we use QPSK modulation
m=2;
K=8;
Eb=1.73;
sigma_b=sqrt(2);  
sig_b_sq=2;
%% Create Symbols
data= randi([0 M-1],n,1);
b=sqrt(2).*pskmod(data,M,pi/4);
b_t=b.';
if 0
    figure
    stem(1:10, b_t(1:10));
    axis ('tight');
    ylim([-0.5,0.5])
    title('b[n]');
end
%% Create filters and plot b in different parts of the process
% gTX
g_TX=rcosdesign(beta,span,Fs/span,'sqrt');
b_upsample=upsample(b,length(g_TX));%8 sample per symbol (with zeros)
s=conv(b_upsample,g_TX);
if 0 % 1 for plots
    figure 
    subplot(2,1,1)
    plot(real(s(1:10*length(g_TX))))
    axis('tight')
    title('Re\{s[n]\}')
    subplot(2,1,2)
    plot( imag(s(1:10*length(g_TX))))
    axis('tight')
    title('Im\{s[n]\}')
end

% Channel
g_channel = zeros(1,floor(1.10*length(g_TX)));
g_channel(1)=-0.8;
g_channel(length(g_TX))=1;
g_channel(floor(1.10*length(g_TX)))=-0.3;
y=conv(s,g_channel);
if 0
    figure
    plot(real(y(1:10*length(g_TX))));
    axis ('tight');
    ylim([-0.5,0.5])
    title('y[n] = s[n]*g_c[n]');
end

% gRX = gTX
g_RX=g_TX;
r_n=conv(y,g_RX,'same');
if 0
    figure
    plot(real(r_n(1:10*length(g_TX))));
    axis('tight')
    title('Re\{r[n]\}');
end
%% b Upsample m=2
b_upsample = upsample(b,2);
b_upsample = [zeros(2,1); b_upsample];

%% f[k]
f=@(n) double(-0.8*kroneckerDelta(sym(n))+kroneckerDelta(sym(n-m))-0.3*kroneckerDelta(sym(n-1.5*m)));
n=[0:3];
r_n=conv(b_upsample,f(n),'same');

if 0
    figure
    plot(real(r_n(1:20)));
    axis('tight')
    ylim([-1.8,1.8])
    title('b[n] after passing through f[k]');
    figure
    stem(f(n));
    axis('tight')
    title('f[k]');
end

%% Noise
SNR=[-5:5:40,inf];
N0=Eb./(10.^(SNR./10));
r_with_noise=zeros(length(SNR),length(r_n));
for i=1:length(N0)
   %%normrnd(mu,sigma) generates a random number 
   %%from the normal distribution with mean parameter mu 
   %%and standard deviation parameter sigma.
   r_with_noise(i,:)=r_n.' + (normrnd(0,sqrt(N0(i)/2),1,length(r_n))+1i*normrnd(0,sqrt(N0(i)/2),1,length(r_n)));%awgn(r_n,SNR(i)); 

end

%% Perror function 
Perror = @(SIR) 1-(1-qfunc(sqrt(SIR))).^2; % just An anonymous function to calc perror
U=[1,-0.8,0,0,0;
   -0.3,0,0,0,0;
   0,1,-0.8,0,0;
   0,-0.3,0,0,0;
   0,0,1,-0.8,0;
   0,0,-0.3,0,0;
   0,0,0,1,-0.8;
   0,0,0,-0.3,0]; 
Cw=eye(8); %% eye is identity matrix 
u0=U(:,3);

%% Direct Decision
DIRECT_z=zeros(length(SNR),length(b));
DIRECT_Perror=zeros(length(SNR),1);
% Y = sign(x) returns an array Y the same size as x, where each element of Y is:
% 1 if the corresponding element of x is greater than 0.
% 0 if the corresponding element of x equals 0.
% -1 if the corresponding element of x is less than 0.
% x./abs(x) if x is complex.

%decide based on sign only
DIRECT_z=sign(real(r_with_noise(:,3:2:end)))+1i*sign(imag(r_with_noise(:,3:2:end)));
tmp_error=0;
for i=1:length(SNR)%%for each SNR 
    for j=3:2:length(b_upsample) %%loop 
        j_ = ceil(j/2); % to get 1 when j = 1 
        %real part error
        tmp_error=tmp_error+...
            sum(sign(real(b(j_-1)))~=sign(real(DIRECT_z(i,j_-1))));
        %image part error
        tmp_error=tmp_error+...
            sum(sign(imag(b(j_-1)))~=sign(imag(DIRECT_z(i,j_-1)))); 
    end 
    %%devide by 2*length(b) because we have re and im
    DIRECT_Perror(i)=tmp_error/(2*length(b));
    tmp_error=0;%reset before next SNR 
end

%% zero forcing
e=[0;0;1;0;0];
czf=U*inv(U'*U)*e;
r_n_sample=[zeros(length(SNR),2), r_with_noise , zeros(length(SNR),3)];
ZF_z=zeros(length(SNR),length(b));
PerrorZF=zeros(length(SNR),1);
tmp_error=0;
for i=1:length(SNR) %for each SNR
    for j=1:2:length(r_n)-(K-2) %loop thrugh rn
        j_ = ceil(j/2); % to get 1 when j = 1 
        %vector multiplication with czf and 8 samples of rn each time and
        %stride 2 
        ZF_z(i,j_)=(czf')*(r_n_sample(i,j:j+7).');
        %comparing z and b to calc the perror (devide by length at the end)
        tmp_error=tmp_error+sum(sign(real(b(j_)))~= sign(real(ZF_z(i,j_))) ||...
                                sign(imag(b(j_)))~= sign(imag(ZF_z(i,j_))));
    end
    PerrorZF(i)=tmp_error/(length(b));
    tmp_error=0;
end


%% ZF bound
ZF_SIR=zeros(1,length(N0));
ZF_bound=Perror((sig_b_sq*czf'*u0)./((czf'*Cw*czf).*(N0./2)));

%% MMSE-LE
p=sig_b_sq.*u0;
MMSE_z = zeros(length(SNR),length(b));%each row is different SNR so 11 rows 
MMSE_Perror = zeros(length(SNR),1); %for each SNR
MMSE_SIR = zeros(1,length(N0));
for i=1:length(SNR)
    R=sigma_b^2*U*U'+Cw*N0(i)/2;
    MMSE_c=inv(R)*p;
    for j=1:2:length(r_n)-(K-2)%loop thrugh rn
        j_ = ceil(j/2); % to get 1 when j = 1 
        %vector multiplication with c_mmse and 8 samples of rn each time and
        %stride 2 
        MMSE_z(i,j_) = (MMSE_c')*(r_n_sample(i,j:j+7).');
        %comparing z and b to calc the perror (devide by length at the end)
        tmp_error=tmp_error+ sum(sign(real(b(j_)))~=sign(real(MMSE_z(i,j_))) ||...
                                 sign(imag(b(j_)))~=sign(imag(MMSE_z(i,j_))));
            
    end
    MMSE_SIR(i)=sig_b_sq*(MMSE_c'*u0)^2/...
        (sig_b_sq*(sum(abs(MMSE_c'*U(:,[1:2,4:5])).^2)) + (MMSE_c'*Cw*(N0(i)/2)*MMSE_c));
    MMSE_Perror(i)=tmp_error/(length(b));
    tmp_error=0;
end

%% MMSE bound
MMSE_bound=Perror(MMSE_SIR);

%% MMSE-DFE
Uf=U(:,3:end);
Up=fliplr(U(:,1:2));%flip array
DFE_z=zeros(length(SNR),length(b));
DFE_Perror=zeros(length(SNR),1);
DFE_SIR=zeros(1,length(N0));
for i=1:length(SNR) %for each SNR
    if(i~=length(SNR)) %calc c feedfoward
        cff = inv(sigma_b^2*Uf*Uf'+Cw*N0(i)/2)*sigma_b^2*u0;
    else %for SNR=inf
        cff = inv(sigma_b^2*Uf*Uf'+Cw*N0(i)/2+0.00001*eye(8))*sigma_b^2*u0;
    end
    % c feedback is cff multiplid by Up 
    cfb=-cff'*Up;
    for j = 1:2:length(r_n)-(K-2)%loop thrugh rn
        j_ = ceil(j/2); % to get 1 when j = 1  
        if (j_ == 1) % we dont have feedback so just cff 
            DFE_z(i,j_)=(cff')*(r_n_sample(i,j:j+7).');
        elseif j_ == 2
            % calc a row vec with sign of last sym in first column (and 0 for the next senerio)
            sign_of_last_sym = [sign(real(DFE_z(i,j_-1))) + 1i*sign(imag(DFE_z(i,j_-1))),0];
            %%%      cff part \/
            DFE_z(i,j_)=(cff')*(r_n_sample(i,j:j+7).')+...
                sum(cfb.*sign_of_last_sym);
            %%% now we sum cfb * the sign of the last symbol 
        else
            % calc a row vec with sign of last 2 sym in each column
            sign_of_last_two_sym =[sign(real(DFE_z(i,j_-1)))+  1i*sign(imag(DFE_z(i,j_-1))),...
                                   sign(real(DFE_z(i,j_-2))) + 1i*sign(imag(DFE_z(i,j_-2)))];
            %%%      cff part \/
            DFE_z(i,j_)=(cff')*(r_n_sample(i,j:j+7).')+...
                               sum(cfb.*sign_of_last_two_sym);  
            %%% now we sum cfb * the sign of the last 2  symbols 
        end
        %to calculate error we compare decided sign and original sign
        tmp_error=tmp_error + sum(sign(real(b(j_)))~=sign(real(DFE_z(i,j_))) ||...
                                  sign(imag(b(j_)))~=sign(imag(DFE_z(i,j_))));
    end
    DFE_SIR(i) = sigma_b^2*(cff'*u0)^2 /...
                (sigma_b^2*(sum(abs(cff'*U(:,4:5)).^2)) + cff'*Cw*(N0(i)/2)*cff );
    DFE_Perror(i)=tmp_error/(length(b));
    tmp_error=0;
end

%% DFE bound
DFE_bound=Perror(DFE_SIR);

%% MF bound
MF_bound(:) = Perror((sig_b_sq*norm(u0)^2)./(N0./2));

%% plotting perror and bounds
if 1
    figure
    hold on
    plot(SNR, DIRECT_Perror+10^-6,'--r')
    plot(SNR, PerrorZF+10^-6,'--g')
    plot(SNR, [ZF_bound]+10^-6,'-g')
    plot(SNR, MMSE_Perror+10^-6,'--b')
    plot(SNR, MMSE_bound+10^-6,'-b')
    plot(SNR, DFE_Perror+10^-6,'--m')
    plot(SNR, DFE_bound+10^-6,'-m')
    plot(SNR, MF_bound+10^-6,'-')
    set(gca, 'YScale', 'log')
    legend('Direct: Perror',...
           'ZF-LE: Perror','ZF-LE: bound',...
           'MMSE-LE: Perror ','MMSE-LE: bound',...
           'MMSE-DFE: Perror','MMSE-DFE: bound',...
           'MF: bound')
    title('Probability of Error')
    xlabel('E_b/N_0')
    ylabel('P_{error}')
    %grid on
    axis tight
end

