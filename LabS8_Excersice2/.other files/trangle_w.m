
fs=1000; 
tStop=1;
T=0.2; 
Amp=9;

tt=0:(1/fs):tStop;
qq=rem(tt,T);
xx=Amp*(abs(qq-(0.5*T))-0.25*T);

figure('Name','Triangle Wave.');
clf;
plot(tt,xx);title('Triangle Wave.');
xlabel('time, [s]')
ylabel('Amplitude')