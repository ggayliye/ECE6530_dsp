% Skylar S and Kyle G
%DSP Fall 2024
% Question 1 complete
% Question 2
%This code segment extracts our data from the text file.
filename = 'C:\Users\computer\Downloads\dataset+for+adl+recognition+with+wrist+worn+accelerometer\HMP_Dataset\Climb_stairs\Accelerometer-2011-03-24-10-24-39-climb_stairs-f1.txt';
delimiterIn = ' ';
headerlinesIn = 1;
A = importdata(filename,delimiterIn,headerlinesIn);

%Code to display data
% for k = [1,2,3]
%    disp(A.colheaders{1, k})
%    disp(A.data(:, k))
%    disp(' ')
% end

% Question 3
%This code segment plots our data pre conversion

x = A.data(:,1);

y = A.data(:,2);

z = A.data(:,3);

figure(1)

plot3(x,y,z)


%This code segments plots our normalized data.

g=9.8;
C= exp(A.data);
B= -1.5*g + (3*g)*(C/63);

x2 = B(:,1);

y2 = B(:,2);

z2 = B(:,3);

figure(2)
plot3(x2,y2,z2)



% problem 4

%30 Hz sampled cos signal
% We can see for this plot that we have a completely incorrect frequency,
% which makes sense because we are sampling at the same frequency as our
% signal.

f=30; %30 Hz

%t= [-.5,.5];

t = linspace(-.5,.5,30) %30 Hz sampled cos signal

%t = 0:0.01:(2*pi);
%c = cos(2*pi*f*t);
c = cos(2*pi*f*t);

figure(3)

plot(t,c)







 %60 Hz sampled cos signal
% We can see for this plot which is at nyquist frequency we finally get the correct frequency, but still
% have some major aliasing going on.

f=30; %30 Hz

%t= [-.5,.5];

t2 = linspace(-.5,.5,60) %60 Hz sampled cos signal

%t = 0:0.01:(2*pi);
%c = cos(2*pi*f*t);
c2 = cos(2*pi*f*t2);  

figure(4)

plot(t2,c2)






 %80 Hz sampled cos signal
% We see that at 80 we are getting rid of the majority of aliasing but
% still have a pretty rought signal.

f=30; %30 Hz

%t= [-.5,.5];

t3 = linspace(-.5,.5,80) %80 Hz sampled cos signal

%t = 0:0.01:(2*pi);
%c = cos(2*pi*f*t);
c3 = cos(2*pi*f*t3);

figure(5)

plot(t3,c3)



% Question 5

f=30; %30 Hz

%t= [-.5,.5];

t4 = linspace(-.5,.5,32) %32 Hz sampled cos signal

%t = 0:0.01:(2*pi);
%c = cos(2*pi*f*t);

c4 = cos(2*pi*f*t4);

% Because the arrays of our 32 hz sampled cosine signal and our x vector data values do
% not match, I actually had to first normalize and then interpolate them to 
% add them and put them on the same graph like asked for
% in the homework. This is interesting because the acceleration data from a
% stair climb should be occilatory, perhaps we were told to pick 32 hz
% because the sampling of the stair climb was at roughly the frequency of
% the acceleration change?


R1=rescale(c4)
R2=rescale(x2')
R1(numel(R2)) = 0;

%R3 = resize(R1', 270)

%yi = interp1(t4, R1, R2, 'linear','extrap');

output_final= R1 + R2

figure(6)
% 
% plot(t4, R1, '-b')
% hold on
% plot(R2, yi, 'pg')
% hold off
plot(linspace(-.5,.5,270),output_final)


% In the end, I ended up plotting the normalized values of both against
% each other in extended arrays. You can see the undersampled cosine wave
% at the beginning of the graph. 













