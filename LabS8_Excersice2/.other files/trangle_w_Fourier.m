
k=[1, 3, 5, 15];
a_k = -2./(pi^2.*k.^2);

% Compute the ratios a3/a1, a5/a1 and a15/a1.

a3_a1=a_k(1,2)/a_k(1,1);
a5_a1=a_k(1,3)/a_k(1,1);
a15_a1 =a_k(1,4)/a_k(1,1);

% fs=1000; 
% tStop=1;
% T=0.2; 
% Amp=9;
% 
% tt=0:(1/fs):tStop;
% qq=rem(tt,T);
% xx=Amp*(abs(qq-(0.5*T))-0.25*T);
% 
% figure('Name','Triangle Wave.');
% clf;
% plot(tt,xx);title('Triangle Wave.');
% xlabel('time, [s]')
% ylabel('Amplitude')