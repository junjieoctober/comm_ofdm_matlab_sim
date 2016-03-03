function ofdm_example_802_11a()

close all; 
fs = 1e6; 
Ts = 1/fs; 
N  = 64; 
t = 0 : Ts : (N-1)*Ts; 

f1 = fs/10; 

nsc = 53; 

np = .01; 
% Long Preamble
LongPreamble = complex(...
    [ 1,  1, -1, -1,  1,  1, -1,  1, -1,  1,  1,  1,...
    1,  1,  1, -1, -1,  1,  1, -1,  1, -1,  1,  1,  1,  1, 0,...
    1, -1, -1,  1,  1, -1,  1, -1,  1, -1, -1, -1, -1, -1,...
    1,  1, -1, -1,  1, -1,  1, -1,  1,  1,  1,  1].', 0);

LongPreamble(25) = 1 + 1i; 
LongPreamble(26) = 1 - 1i; 
LongPreamble(28) = 1 - 1i; 
LongPreamble(29) = 1 - 1i; 


df = fs/N; 
subc = size(LongPreamble,1)-1; 
K = (-subc/2:subc/2)';
y = exp(-j*2*pi*df*K*t);  %53*N matrix

B=repmat(LongPreamble,1,53);
y1 = zeros(size(y)); 

for i = 1 : nsc
    y1(i,:) = LongPreamble(i)*y(i,:);
end


%y1 = (y'*B)';


y2 = sum(y1);


plot(t,real(y2),'*-r');
hold on; 
plot(t,imag(y2),'+-b');
legend('real','imag'); 
title('raw data'); 


r = y2' + np*complex(randn(N,1),randn(N,1)); 
hold on;
plot(t,r,'+-b');
legend('tx','rx'); 

dec = zeros(size(LongPreamble)); 

for i = 1 : nsc
    dec(i) = y(i,:)*r;
end

dec = dec/N; 

figure;
plot(real(LongPreamble),'*-');
hold on; 
plot((real(dec)),'+-r'); 
legend('truth preamble real','decoded real');

figure;
plot(imag(LongPreamble),'*-');
hold on; 
plot((imag(dec)),'+-r'); 
legend('truth preamble imag','decoded imag');

Fr = fftshift(ifft(r))/N; 
figure; 
plot(real(Fr),'+b-'); 
hold on; 
plot(imag(Fr),'dr-'); 

return; 

%y = LongPreamble*y(:,
zeropad = (N-subc)/2; 
y = [zeros(zeropad,1);y;zeros(zeropad-1,1)]; 
%y = fftshift(y); 
Y = ifft(y); 

Y = Y + np*complex(randn(N,1),randn(N,1)); 

rx = fft(Y); 

F = (-N/2: (N-1)/2)/N*fs; 
plot(F,fftshift(abs(Y)),'*-');
hold on;
plot(F,real(Y),'+r-'); 
grid on; 
plot(F,real(rx),'dk-'); 
plot(F,imag(rx),'oc-'); 

