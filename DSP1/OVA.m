function y = OVA(x, h, K)
     N = length(x);
     P = length(h);
     L = K - P + 1; 
     h = [h zeros(1,(K - P))]; %we need to add K - P zeros after the filter
     H = fft(h); 
     t = 1;
     y = zeros (1,(P+N-1)); %y has P+N samples (let matlab know length, not necessary)
     while t < N - L %while loop to calculate each frame
        X = [x(t:(t+L-1)) zeros(1,P-1)];
        wt = ifft((fft(X)).* H);
        yt = y(t:t+K-1) + wt;
        y (t:t+K-1) = yt; %add the frame to the finel y
        t = t + L; %move to next frame
     end 
   
end