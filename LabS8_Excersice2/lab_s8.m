% ========================================================================== 
% ECE 5530/ ECE 6530 Digital Signal Processing, ECE Dep., The University of Utah 
% Group Project
% Team: Kyle G. Gayliyev, Skylar Stockham, Eddie Franco
% ========================================================================== 
close all;
clear;

%% Sampling, Aliasing and Reconstruction : Lab S8 Exercise.
%%  There are 2 Lab Exercises : 2.1 Spectrogram for a Chirp that Aliases 
%%  and 2.2 Spectrogram of Periodic Signal
%% Overvew
% The objective of this lab is to study further the spectral content of 
% signals analyzed via the spectrogram. There are several specific steps 
% that will be considered in this lab: 
%
% 1. Synthesize a linear-FM chirp with a MATLAB M-file, and display its 
% spectrogram. Choose the chrip parameters so that aliasing will happen. 
% 
% 2. Synthesize a periodic triangle wave with a MATLAB M-file, and display 
% its spectrogram. Relate the harmonic line spectrum to the fundamental 
% period of the triangle wave. 3. Compare spectrograms using different 
% scales for amplitude: decibels (dB) for amplitude versus linear amplitude. 
% 
% 4. Examine details of the harmonic lines in the dB spectrogram of the 
% triangle wave. 
% 
% 5. Spectrogram: make a spectrogram of your voice signal, and relate the 
% harmonic line spectrum to your previous measurement of pitch period.

% We'll learn more about the connection between the time-domain definition 
% of the signal and its frequency-domain content.

 %% 2.1 Spectrogram for a Chirp that Aliases
% TASK: Use the code provided in the pre-Lab section as a starting point 
% in order to write a MATLAB script or function that will synthesize 
% a “chirp” signal. Then use that M-file in this section.

% a) What happens when we make a signal that “chirps” up to a very high 
% frequency, and the instantaneous frequency goes past half the sampling 
% rate? 
% Generate a chirp signal that starts at 1000 Hz when t=0 s, and 
% chirps up to 11,000 Hz, at t= 4 s. Use fs= 4000 Hz. 
% Determine the parameters needed in (1).

fs = 4000; % -Number of time samples per second
dt = 1/fs; %time domain conversion
tStart = 0;  % in sec
tStop = 4; % in sec
tt = tStart:dt:tStop; % total duration with sampling intervals

f1 = 1000; % [Hz], also fzero =f1 since t=0sec.
f2= 11000; % [Hz]

fzero=f1;

% Calculate slope parameter (mu)
mu = (f2-f1)/(2*(tStop-tStart));

phi = 2*pi*rand; %-- random phase

psi=2*pi*mu.*(tt.^2)+2*pi*fzero.*tt+phi; % chirp signal equation
chirp_signal = real( 7.7*exp(1i*psi) );
% soundsc( chirp_signal, fs ); %-- uncomment to hear the sound

figure('Name','Chirp Signal in Time Domain.');
clf;
plot(tt,chirp_signal); title('Chirp Signal in Time Domain.');
xlabel('t, [s]');
ylabel('Amplitude');

% Plot the spectrogram of the chirp signal
figure('Name','Spectrogram Of The Chirp Signal.');
clf;
spectrogram(chirp_signal, 256, 250, 256, fs, 'yaxis');
title('Spectrogram of the Chirp Signal');
xlabel('Time (s)');
ylabel('Frequency (kHz)');
colorbar;

% When a signal "chirps" up to a very high frequency, and the instantaneous 
% frequency goes past half the sampling rate it means the signal's 
% frequency gradually increases over time, reaching a significantly higher 
% frequency at the end of the chirp; if the final frequency exceeds half 
% the sampling rate of the system recording it, the signal will experience 
% aliasing, meaning the high frequencies will be misinterpreted as lower 
% frequencies due to the limitations of digital sampling, distorting the 
% signal's representation in the recorded data.

% Aliasing is a phenomenon where high-frequency components appear as lower 
% frequencies. In the spectrogram,  you will notice that once the chirp 
% frequency exceeds 2000Hz, it folds back into the frequency range below  
% 2000Hz, creating mirrored frequency components.This results in distortion 
% because the true frequency of the signal is not correctly represented in 
% the sampled signal.


% (b) Generate the chirp signal in MATLAB and make a spectrogram with a 
% short section length, LSECT, to verify that you have the correct starting 
% and ending frequencies.2 For your chosen LSECT, determine the section 
% duration TSECT in secs

% Spectrogram parameters
L_SECT = 128;  % Section length (number of samples per segment)
overlap = L_SECT - 10;  % Overlap of sections, chosen to provide smoother transitions
nfft = 256;  % Number of FFT points

% Plot the spectrogram of the chirp signal
figure('Name','Spectrogram of the Chirp Signal with Short Section Length.');
clf;
spectrogram(chirp_signal, L_SECT, overlap, nfft, fs, 'yaxis');
title('Spectrogram of the Chirp Signal with Short Section Length');
xlabel('Time (s)');
ylabel('Frequency (kHz)');
colorbar;

% Calculate section duration (T_SECT) in seconds
T_SECT = L_SECT / fs;
disp(['Section duration T_SECT: ', num2str(T_SECT), ' seconds']);

% (c) Explain why the instantaneous frequency seen in the spectrogram is 
% goes up and down between zero and fs=2, i.e., it does not chirp up to 
% 11,000Hz. 
% There are two effects that should be accounted for in your explanation. 
% Note: If possible listen to the signal to verify that the spectrogram is 
% faithfully representing the audio signal that you hear.

% Effect #1: Aliasing Effect
% The sampling frequency is 4000 Hz, which means the Nyquist frequency 
% is: 2000 Hz.
% When the instantaneous frequency of the chirp exceeds the Nyquist 
% frequency, aliasing occurs.
% Aliasing is the phenomenon where high frequencies are misrepresented as 
% lower frequencies 
% when sampled at a rate that is insufficient to capture them accurately.
% In the case of this chirp signal: As the frequency of the chirp exceeds 
% 2000Hz, 
% it is folded back into the lower frequency range (below the Nyquist 
% frequency).
% This folding effect causes the frequency to appear to go down once it 
% exceeds 2000Hz, 
% resulting in the up and down pattern seen in the spectrogram.
% This behavior is consistent with the aliasing theorem, which states 
% that frequencies above 
% the Nyquist limit appear at incorrect, lower frequencies in the sampled 
% signal.

% Effect #2:. Periodic Nature of Frequency Representation in the Spectrogram
% The Discrete Fourier Transform (DFT), which is used in calculating the 
% spectrogram, treats frequency content as periodic.
% The frequency axis is considered to be periodic, with the range extending 
% from 0 to f_s
% but practically folding at f_s/2.
% Frequencies between 0 and f s/2 are represented correctly.
% Frequencies beyond f_s/2 (up to f_s) are aliased and appear as if they 
% are in the range 0 to f s/2.
% This periodic folding creates a mirror effect in the spectrogram.
% Therefore, the instantaneous frequency appears to go back down after 
% reaching f s/2. 
% This is a direct result of how the DFT operates and how the spectrogram 
% interprets frequencies beyond the Nyquist frequenc

%% 2.2 Spectrogram of Periodic Signal
% A periodic signal is known to have a Fourier Series, which is usually 
% described as a harmonic line spectrum because the only frequencies 
% present in the spectrum are integer multiples of the fundamental 
% frequency. With the spectrogram, it is easy to exhibit this harmonic 
% line characteristic.
    %% 2.2.1 Spectrogram of Periodic Triangle Wave
% (a) Write a simple MATLAB script that will generate a periodic triangle 
% wave once the period is given. The DClevel of the triangle wave should 
% be zero, and the peak amplitude should be equal to 0.5. Here is a MATLAB
% one-liner that can form the basis of this script: tt_=0:(1/f_s):t_Stop;
% qq=rem(tt_,T);xx=Amp*(abs(qq-(0.5*T))-0.25*T); The values of f_s, t_Stop, 
% T, Amp will have to be determined.

% Parameters for the triangle wave
f_s = 1000;        % Sampling frequency in Hz
t_Stop = 1;        % Duration of the signal in seconds
T = 0.2;          % Period of the triangle wave in seconds
Amp = 2;          % Amplitude scaling factor to achieve peak amplitude of 0.5

% Generate the time vector
tt_ = 0:(1/f_s):t_Stop;  % Time vector from 0 to t_Stop with step size 1/f_s

% Generate the triangle wave using the provided one-liner
qq = rem(tt_, T);  % Folding time into the period T (creates sawtooth shape)
xx = Amp * (abs(qq - (0.5 * T)) - 0.25 * T);  % Create the triangle wave

% Plot the triangle wave in the time domain
figure('Name','Triangle Wave in Time Domain.');
clf;
plot(tt_, xx, 'b', 'LineWidth', 1.5);
title('Triangle Wave in Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% (b) Generate a triangle wave with a period of 10msec, using a sampling 
% rate of f_s D 10000Hz. The duration should be 3secs. Then make a plot 
% of a short section of the signal consisting of 3–5 periods to verify that 
% you have the correct time waveform.

% Parameters for the triangle wave
f_s1 = 10000;       % Sampling frequency in Hz
t_Stop1 = 3;        % Duration of the signal in seconds
T_1 = 0.01;         % Period of the triangle wave in seconds (10 ms)
Amp_1 = 1;          % Amplitude scaling factor for peak amplitude of 1

% Generate the time vector
tt_1 = 0:(1/f_s1):t_Stop1;  
% Time vector from 0 to t_Stop1 with step size 1/f_s1

% Generate the triangle wave using the given one-liner
qq_ = rem(tt_1, T_1);  % Folding time into the period T_1 (creates sawtooth shape)
xx_ = Amp_1 * (abs(qq_ - (0.5 * T_1)) - 0.25 * T_1);  % Create the triangle wave

% Plot the entire triangle wave in the time domain
figure('Name','Triangle Wave in Time Domain.');
clf;
plot(tt_1, xx_, 'b', 'LineWidth', 1);
title('Triangle Wave in Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
xlim([0 t_Stop1]);
ylim([-0.003 0.003]);

% Plot a short section of the triangle wave consisting of 3-5 periods
period_count = 5;  % Number of periods to plot (3-5 periods)
time_window = period_count * T_1;  % Time window to cover 5 periods

% Find the indices to plot
indices_to_plot = find(tt_1 <= time_window);

% Plot the short section
figure('Name','Short Section of the Triangle Wave (3-5 Periods).');
clf;
plot(tt_1(indices_to_plot), xx_(indices_to_plot), 'r', 'LineWidth', 1.5);
title('Short Section of the Triangle Wave (3-5 Periods)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
xlim([0 time_window]);
ylim([-0.003 0.003]);


% (c) Make a spectrogram with a long section duration.
% It is important to pick a section duration that is equal to an integer 
% number of periods of the periodic triangular waveform created in the 
% previous part. Define TSECT to get exactly 5 periods, and then determine 
% the section length LSECT (an integer) to be used in plotspec.

% Define section duration for spectrogram
T_SECT_ = 5 * T_1;  % Section duration for exactly 5 periods
L_SECT_ = T_SECT_ * f_s1;  % Section length in terms of samples

% Plot the spectrogram of the triangle wave with long section duration
figure('Name','Spectrogram of the Periodic Triangle Wave (Long Section Duration).');
clf;
spectrogram(xx, L_SECT_, L_SECT_/2, L_SECT_, f_s1, 'yaxis');
title('Spectrogram of the Periodic Triangle Wave (Long Section Duration)');
xlabel('Time (s)');
ylabel('Frequency (kHz)');
colorbar;

% (d) You should expect to see a “harmonic line spectrum” in the 
% spectrogram. Since frequency is along the vertical axis,the harmonic 
% lines will appear as horizontal lines in the spectrogram. Make a list of 
% all the harmonic frequencies that you can see in the spectrogram.

harmonic_frequencies = [100, 300, 500, 700, 900, 1100, 1300, 1500, 1700, ... 
    1900, 2100, 2300, 2500, 2700, 2900, 3100, 3300, 3500, 3700, 3900, ...
                        4100, 4300, 4500, 4700, 4900];


% (e) Determine the fundamental frequency for the harmonic lines.

% Given period of the triangle wave (in seconds)
T_2 = 0.01;  % Period in seconds (10 ms)

% Calculate the fundamental frequency
f0 = 1 / T_2;

% Display the fundamental frequency
disp(['The fundamental frequency is: ', num2str(f0), ' Hz']);

% (f) Measure the amplitudes of the first and third harmonic lines by using
% MATLAB’s Data Cursor after zooming in on those parts of the spectrogram image. 
% Record the values for the amplitudes and compute the ratio.
% Example power values from Data Cursor (in dB)
P1_dB = -71.23;  % Power of the first harmonic in dB
P3_dB = -90.29;  % Power of the third harmonic in dB

% Convert power values to linear amplitudes
A1_linear = 10^(P1_dB / 20);
A3_linear = 10^(P3_dB / 20);

% Calculate the ratio of the third harmonic amplitude to the first harmonic amplitude
ratio = A3_linear / A1_linear;

% Display the ratio
disp(['The ratio of the amplitude of the third harmonic to the first harmonic is: ', num2str(ratio)]);
    
    %% 2.2.2 Decibels (dB): Seeing Small Values in the Spectrogram
% Where did all the harmonics go?
% The answer is that the higher harmonics have amplitudes that are too 
% small to be seen in a spectrogram that displays values with a linear 
% amplitude. Instead, a logarithmic amplitude scale is needed.

% TASK : Answer the following questions about decibels:

% a) In the language of dB, a factor of two is “6 dB.” In other words, 
% if B2 is 6 dB bigger than B1, then it is twice as big (approximately). 
% Explain why this statement is true.

% Answer: In the language of dB, a factor of two corresponds to an 
% increase of 6 dB. In other words, if B2 is 6 dB larger than B1, then  B2
% is approximately twice as large as B1.
% The relationship between a value in linear scale and decibels is 
% given by: 
% Power in dB = 20log10(B2/B1).
% To understand why 6 dB corresponds to a factor of two, let’s substitute 
% and evaluate this expression.
% Linear to dB Conversion: Given that the relationship between two values 
% B1 and B2 in terms of dB is:
% 20log10(B2/B1) =6dB. Solving for the Ratio, we get:
% B2/B1=~2.
% This means that an increase of 6 dB corresponds to a value that is twice 
% as large in linear terms.

% (b) Determine the dB difference between a1 and a3. In other words, a3 is 
% how many dB below a1. Furthermore, explain why the dB difference depends 
% only on the k indices.

% Given Fourier coefficients formula for a triangle wave
k1 = 1;  % First harmonic index
k3 = 3;  % Third harmonic index

% Calculate a1 and a3 using the given formula
a1 = -2 / (pi^2 * k1^2);
a3 = -2 / (pi^2 * k3^2);

% Calculate the ratio in linear terms
ratio_linear = abs(a1) / abs(a3);

% Calculate the dB difference
dB_difference = 20 * log10(ratio_linear);

% Display the result
disp(['The dB difference between a1 and a3 is: ', num2str(dB_difference), ' dB']);

% The dB difference depends only on the k indices because the amplitude of 
% each harmonic a_k is inversely proportional to the square of k
% This means that the amplitude decreases as  k increases, and this 
% relationship is governed by k2. Therefore, when taking the ratio of 
% two coefficients a_1 and a_3
% Since this ratio depends purely on the squares of the indices, 
% the resulting dB difference also depends only on these indices. 
% This is why the dB difference between a_1 and a_3 depends only on the 
% values of k and not on the specific waveform amplitude values.

% (c) Determine (in dB) how far a15 is below a1 for the periodic triangular 
% wave.

% Given Fourier coefficients formula for a triangle wave
k1 = 1;  % First harmonic index
k15 = 15;  % Fifteenth harmonic index

% Calculate a1 and a15 using the given formula
a1 = -2 / (pi^2 * k1^2);
a15 = -2 / (pi^2 * k15^2);

% Calculate the ratio in linear terms
ratio_linear = abs(a1) / abs(a15);

% Calculate the dB difference
dB_difference = 20 * log10(ratio_linear);
% This confirms that a_15 is approximately 43.52 dB below a_1 for the 
% periodic triangular wave.

% Display the result
disp(['The dB difference between a1 and a15 is: ', num2str(dB_difference), ' dB']);

%% 2.2.3 Spectrogram in dB

% A variation of the SP-First function plotspec has been written to 
% incorporate the dB amplitude scale. This new function is called 
% plotspecDB, and its calling template is shown below.


% a) Create a “dB-Spectrogram” for the 10-msec periodic triangular wave 
% generated in Sect. 2.2.1. Use a dBrange equal to 80 dB. Notice that 
% many more spectrum lines are now visible. List the frequencies of all 
% the harmonic spectrum lines, or give a general formula.

% b) Generate a second triangle wave by changing the period to 20 msec. 
% Then make the dB-Spectrogram of this 20-msec triangle wave, being careful 
% to select the section duration as an integer number of periods. 
% From the spectrogram, determine the fundamental frequency and also the 
% frequency of the highest harmonic line. Also, determine the harmonic 
% number for the highest frequency, e.g., the 17th or 31st, etc.

% c) For the 20-msec triangle wave, measure the amplitudes (in dB) of the 
% first and third harmonic lines by using MATLAB’s Data Cursor after 
% zooming in on those parts of the spectrogram image. Compare the dB 
% difference to the ratio obtained in Section 2.2.1, part (f).

% d) Change the period to 4 msec and make another dB-Spectrogram. 
% Be careful to select the section duration as an integer number of periods. 
% This period is shorter but the frequency separation of the harmonic lines 
% is greater. Notice that this inverse relationship was also true when 
% comparing the 20 msec case to the 10 msec case.


function him = plotspecDB(xx,fsamp,Lsect,DBrange)
%PLOTSPECDB plot a Spectrogram as an image
% (display magnitude in decibels)
% usage: him = plotspec(xx,fsamp,Lsect,DBrange)
% him = handle to the image object
% xx = input signal
% fsamp = sampling rate
% Lsect = section length (integer, power of 2 is a good choice)
% amount of data to Fourier analyze at one time
% DBrange = defines the minimum dB value; max is always 0 dB
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %% Questions. 
%  % Answer the following questions for the lighthouse picture,
%  % downsampled by two.
%  % a) Describe how the aliasing appears visually. Compare the original to 
%  % the downsampled image. Which parts of the image show the aliasing 
%  % effects most dramatically?
%  % b) This part is challenging: explain why the aliasing happens in the 
%  % lighthouse image by using a “frequency domain” explanation. 
%  % In other words, estimate the frequency of the features that are 
%  % being aliased. Give this frequency as a number in cycles per pixel. 
%  % (Note that the fence provides a sort of “spatial chirp” where the 
%  % spatial frequency increases from left to right.) Can you relate your 
%  % frequency estimate to the Sampling Theorem?
%  % You might try zooming in on a very small region of both the original 
%  % and downsampled images.
% %% Answers. 
% % First, let's downsample the lighthouse picture by two. To do so,
% % we upload the "lighthouse.mat" given file. Then, we apply the following
% % formula to download it and save the downloaded file into a new
% % variable.
% 
% % Upload the "lighthouse.mat" given file
% % lighthouse = importdata('lighthouse.mat');
% load lighthouse.mat;
% 
% % Downsampling the lighthouse image by factor of 2
% %Formula: wp = ww(1:p:end,1:p:end); %downsampling the "ww" by a factor of p.
% downsampled_lh_2 = xx(1:2:end,1:2:end); % size : 163x213
% whos;
% 
% %% 3.1a) 
% % Let's draw the image so that we can see and compare the results.
% % Plot the original lighthouse image 
% 
% figure('Name','Lighthouse Original Image.');
% clf;
% imshow(xx); title('Lighthouse Original Image.');
% 
% % Plot the downsampled lighthouse image 
% figure('Name','Lighthouse Downlampled Image by Factor of 2.');
% clf;
% imshow(downsampled_lh_2);
% title('Lighthouse Downlampled Image by Factor of 2.');
% 
% % Comparing both on the same figure
% figure('Name','Comparing Lighthouse Original VS Downlampled Image by Factor of 2.');
% clf;
% 
% subplot(1, 2, 1); imshow(xx); 
% title('Original Image');
% subplot(1, 2, 2); imshow(downsampled_lh_2); 
% title('Downsampled Image by Factor of 2.');
% 
% % subplot(2,1,1);
% % imshow(xx);
% % subplot(2,1,2);
% % imshow(downsampled_lh);
% 
% % Aliasing appears visually as jagged edges, stair-step patterns, or 
% % a "pixelated" look, particularly along sharp transitions between 
% % colors or high-contrast areas in the image.
% % In this example, the downsampled image appears less smooth and more 
% % distorted compared to the original lighthouse picture. Due to aliasing, 
% % some part of the fence, and the lighthouse itself appear blurred or 
% % distorted in the downsampled version. The aliasing effect is most 
% % noticeable in areas with fine details, such as the roof shingles, 
% % the fence, and the edges of the windows and doors.
% 
% % In summary, where high spatial frequencies are present in the downsampled 
% % image, aliasing is visually represented as a distorted, jagged pattern. 
% % The aliasing effect is especially apparent in regions with fine details, 
% % such the fence, the windows and door margins, and the roof shingles, 
% % when comparing the downsampled and original lighthouse image.
% 
% %% 3.1b) 
% 
% % Frequency Domain Analysis
% reg_original = xx(100:150, 100:150); %Extract the regions of interest
% reg_downsampled = downsampled_lh_2(50:75, 50:75);
% 
% % Perform 2D Fourier Transform to analyze the frequency
% fft_original = fft2(reg_original); 
% fft_downsampled = fft2(reg_downsampled);
% 
% % Frequency Estimation and Comparation
% figure('Name','Frequency Estimation and Comparation.');
% subplot(1, 2, 1); imagesc(abs(fftshift(fft_original))); 
% title('Original Image Frequency');
% subplot(1, 2, 2); imagesc(abs(fftshift(fft_downsampled))); 
% title('Downsampled Image by Factor of 2 Frequency.');
% 
% % The lighthouse image's aliasing is caused by the downsampled image's 
% % inability to completely capture the high frequency information in the 
% % original image. This can be explained in the frequency domain by the 
% % fact that the frequency of the aliased features is more than half of the 
% % sampling rate. This is called the Nyquist rate, and it says that the 
% % sampling rate must to be more than twice the maximum frequency in the 
% % original image to avoid aliasing.
% 
% % The lighthouse image's peak frequency is around one cycle per pixel 
% % because the spatial frequency rises from left to right (the fence 
% % creates a kind of "spatial chirp"). The maximum frequency that the 
% % downsampled image can record is 0.5 cycles per pixel because it only 
% % has half as many pixels as the original. Aliasing would occur at 
% % frequencies higher than this.
% 
% % Because the downsampled image has a lower sampling rate (fewer pixels), 
% % it is unable to capture the high frequency information seen in the 
% % original image, which is why aliasing occurs. The Sampling Theorem, which 
% % asserts that the sample rate must be more than twice the maximum 
% % frequency in the signal in order to adequately capture it, is broken 
% % by this.
% 
% % The reason for the aliasing is that the high frequency information 
% % in the original image is not captured by the downsampled version. 
% % About 150 cycles per pixel is the highest frequency found in the original 
% % image. The sample Theorem states that the sample rate must be at least 
% % twice the highest frequency in the image, or 2*150 = 300 cycles per 
% % pixel, in order to prevent aliasing. Aliasing comes from the downsampled 
% % image's inability to capture high frequency information due to its 
% % sampling rate of only 150 cycles per pixel.
% 
% %% 3.2 Reconstruction of Images
% 
% % When an image has been sampled, we can fill in the missing samples by 
% % doing interpolation. For images, this would be analogous which is part 
% % of the reconstruction process in a D-to-A converter. We could use a 
% % “square pulse” or a “triangular pulse” or other pulse shapes for the 
% % reconstruction.
% 
% % For the reconstruction experiment, we'll use the lighthouse image, 
% % down-sampled by a factor of 3.
% 
% % Downsampling the lighthouse image by factor of 3
% downsampled_lh_3 = xx(1:3:end,1:3:end); % size : 109x142
% whos;
% 
% % Plot the downsampled lighthouse image 
% figure('Name','Lighthouse Downlampled Image by Factor of 3.');
% clf;
% imshow(downsampled_lh_3);
% title('Lighthouse Downlampled Image by Factor of 3.');
% % The goal is to reconstruct an approximation to the original lighthouse 
% % image, which is the size of 256x256, from the smaller down-sampled image.
% 
% %% 3.2a) 
% % The simplest interpolation would be reconstruction with a square pulse 
% % which produces a “zero-order % hold.” Here is a method that works for 
% % a one-dimensional signal (i.e., one row or one column of the image),
% % assuming that we start with a row vector xr1, and the result is the
% % row vector xr1hold.
% 
% xr1 = (-2).^(0:6);
% L = length(xr1); % Length of the row vector xr1.
% nn = ceil((0.999:1:4*L)/4); %<-- Rounds up to the integer part
% xr1hold = xr1(nn); % The result row vector
% 
% % TASK1: Plot the vector xr1hold to verify that it is a zero-order hold 
% % version derived from xr1.
% 
% figure('Name','Zero-order Hold Version.');
% clf;
% subplot(2, 1, 1); stem(xr1hold); hold on; plot(xr1hold); hold off;
% title('Plot of The Zero-order Hold Version.');
% subplot(2, 1, 2); imshow(xr1hold);
% title('Image of The Zero-order Hold Version.');
% 
% 
% % TASK2: Explain what values are contained in the indexing vector nn.
% % Answer: The indexing vector nn contains the indices used for zero-order 
% % hold interpolation. The values that are contained in the indexing vector 
% % nn are the integer values. 
% 
% % TASK3: If xr1hold is treated as an interpolated version of xr1, then 
% % what is the interpolation factor?
% % Answer: The interpolation factor is determined by the length of nn 
% % relative to the original row length. Interpolation_factor = length(nn)/L;
% % Using this formula, 
% 
% Interpolation_factor = length(nn)/L;
% 
% fprintf('The Interpolation Factor =');
% disp(Interpolation_factor);
% 
% % The Interpolation Factor =     4.
% 
% 
% 
% %% 3.2b) 
% % Process all the rows of downsampled_lh_3 to fill in the missing points. 
% % Use the zero-order hold idea from part (a) with an interpolation factor 
% % of 3. Call the result xholdrows. Display xholdrows as an image, and 
% % compare it to the downsampled image downsampled_lh_3; compare the size 
% % of the images as well as their content.
% 
% xholdrows = zeros(size(downsampled_lh_3));
% for i = 1:size(downsampled_lh_3, 1)
%     xr1 = downsampled_lh_3(i, :);  % Take one row at a time
%     nn = ceil((0.999/1/4 * length(xr1)) / 4);
%     xholdrows(i, :) = xr1(nn);
% end
% 
% % Display xholdrows as an image
% figure('Name','Zero-Order Hold Interpolation for Rows.');
% clf;
% imshow(xholdrows);
% title('Zero-Order Hold Interpolation for Rows.');
% 
% % Compare xholdrows to downsampled_lh_3 in terms of size and content.
% 
% % It saves the interpolated rows in the matrix xholdrows and uses the same 
% % zero-order hold interpolation method as in part (a) for every row.
% 
% 
% % The final xholdrows matrix is shown as a picture, representing the full 
% % image with interpolated rows. This illustrates how zero-order hold 
% % interpolation is used to fill in the missing spots.
% 
% size_downsampled_lh_3 = length(downsampled_lh_3);
% size_xholdrows= length(xholdrows);
% 
% fprintf(' The Size of downsampled_lh_3 =');
% disp(size_downsampled_lh_3);
% fprintf(' The Size of xholdrows =');
% disp(size_xholdrows);
% 
% % xholdrows and downsampled_lh_3 both have some amount of data.
% 
% % xholdrows and downsampled_lh_3 picture comparation:
% 
% figure('Name','xholdrows and downsampled_lh_3 Image Comparation.');
% clf;
% subplot(1, 2, 1); imshow(downsampled_lh_3); 
% title('Downsampled Lighthouse Image By Factor of 3.');
% subplot(1, 2, 2); imshow(xholdrows);  
% title('xholdrows Image.');
% 
% %% 3.2c) 
% % Process all the columns of xholdrows to fill in the missing points 
% % in each column and and call the result xhold. 
% % Compare the result (xhold) to the original image lighthouse. 
% 
% % Fill in the missing points for all columns in xholdrows using zero-order
% % hold interpolation:
% 
% xhold = zeros(size(xholdrows));
% for i = 1:size(xholdrows, 2)
%     xr1 = xholdrows(:, i);  % Take one column at a time
%     nn = ceil((0.999/1/4 * length(xr1)) / 4);
%     xhold(:, i) = xr1(nn);
% end
% 
% % Looping through columns,the code iterates through all the columns of the
% % xholdrows matrix obtained in part (b).
% % Like section (a), the zero-order hold interpolation approach is applied 
% % for each column, and the interpolated columns are stored in the matrix 
% % xhold. The interpolation for the full image is finished in this stage.
% 
% figure('Name','xhold VS Original Lighthouse Image Comparation.');
% clf;
% subplot(1, 2, 1); imshow(xx); 
% title('Original Lighthouse Image.');
% subplot(1, 2, 2); imshow(xhold);  
% title('xhold Image.');
% 
% %% 3.2d) 
% % Linear interpolation can be done in MATLAB using the interp1 function 
% % (that’s “interp-one”). Its default mode is linear interpolation, which 
% % is equivalent to using the ’*linear’ option, but interp1 can also do 
% % other types of polynomial interpolation.
% 
% % TASK: For the example of a 1-D signal below, what is the interpolation 
% % factor when converting xr1_1 to xr1linear?
% 
% n1 = 0:6;
% xr1_1 = (-2).^n1;
% tti = 0:0.1:6; %-- locations between the n1 indices
% xr1linear = interp1(n1,xr1_1,tti); %-- function is INTERP-ONE
% 
% figure('Name','Signal Example For Section 3.2d.');
% clf;
% stem(tti,xr1linear); title('a 1-D Signal Example For Section 3.2d.');
% 
% % Interpolation Factor Calculation
% 
% L_xr1_1 = length(xr1_1); % Length of the row vector xr1_1.
% Interpolation_factor_xr1_xr1linear = length(tti)/L_xr1_1;
% 
% fprintf('The Interpolation Factor When Converting xr1_1 to xr1linear =');
% disp(Interpolation_factor_xr1_xr1linear);
% 
% % Answer: The Interpolation Factor Converting xr1_1 to xr1linear = 8.7143.
% 
% %% 3.2e) 
% % In the case of the lighthouse image, we need to carry out a linear 
% % interpolation operation on both the rows and columns of the down-sampled 
% % image downsampled_lh_3. This requires two calls to the interp1 function,
% % because one call will only process all the columns of a matrix.
% 
% % Name the interpolated output image xxlinear.
% 
% 
% % Define the interpolation factor 
% interpolation_f = 3;
% 
% % Initialize a matrix for an interpolated image
% xxlinear = [];
% 
% % Perform linear interpolation for each row
% % Iterate through rows of the downsampled image, perform linear 
% % interpolation on each row separately.
% for i = 1:size(downsampled_lh_3, 1)
%     row = downsampled_lh_3(i, :);
%     xxlinear_row = interp1(1:length(row), row, 1:(1/interpolation_f):length(row));
%     xxlinear = [xxlinear; xxlinear_row];
% end
% 
% %% 3.2f) 
% % Task: Compare xxlinear to the original image lighthouse. Comment on the 
% % visual appearance of the “reconstructed” image versus the original; 
% % point out differences and similarities. 
% % Can the reconstruction (i.e., zooming) process remove the aliasing 
% % effects from the down-sampled lighthouse image?
% 
% % xxlinear and original lighthouse image comparation:
% 
% figure('Name','xxlinear and Original Lighthouse Image Comparation.');
% clf;
% subplot(2, 1, 1); imshow(xx); 
% title('Original Lighthouse Image.');
% subplot(2, 1, 2); imshow(xxlinear);  
% title('xxlinear (Linear Interpolation) Image.');
% 
% 
% % Compared the “reconstructed” image versus the original, we observe that
% % the “reconstructed” image is low quality. It has aliasing. It looks wider
% % and bigger, there are lots of quality issues on the fence, lighthouse
% % structure, house structure, and telescope. The number of items are
% % similar.
% 
% %  The linear interpolation reconstruction (i.e., zooming) process didn't 
% % remove the aliasing effects from the down-sampled lighthouse image.
% 
% % No, a simple reconstruction process like zooming (interpolation) cannot 
% % completely remove aliasing effects from a down-sampled image like a 
% % lighthouse picture; once information is lost due to downsampling, it 
% % cannot be fully recovered by simply zooming in, as the missing details 
% % are not present in the original sampled data.
% 
% % When an image is downsampled, high-frequency details (like fine lines 
% % or textures) can be misinterpreted as lower frequencies due to the 
% % reduced sampling rate, creating the "jagged" aliasing effect.
% 
% % While zooming (which uses interpolation algorithms to fill in missing 
% % pixels) can make an image appear larger, it cannot generate new 
% % information that was lost during downsampling.
% 
% %% 3.2g) 
% 
% % TASK: Compare the quality of the linear interpolation result to the 
% % zero-order hold result. 
% % Point out regions where they differ and try to justify this difference 
% % by estimating the local frequency content. In other words, look for 
% % regions of “low-frequency” content and “high-frequency” content 
% % and see how the interpolation quality is dependent on this factor.
% 
% % Differences between zero-order hold and linear interpolation
% figure('Name','Zero-order Hold and Linear Interpolation Difference.');
% clf;
% subplot(2, 1, 1); imshow(xhold); 
% title('Zero-order Hold Image.');
% subplot(2, 1, 2); imshow(xxlinear);  
% title('Linear Interpolation Image.');
% 
% % Low-frequency content refers to areas of an image with gradual changes 
% % in color or intensity, like a flat background or large smooth surfaces.
% 
% % High-frequency content refers to areas with sharp changes in color or 
% % intensity, like edges between objects, fine details, or textures.
% 
% % Compared to xhold, the linear interpolation image (xxlinear) is more 
% % closely matches the original lighthouse image. Although the 
% % reconstruction procedure can eliminate some aliasing issues from the 
% % down sampled lighthouse image, we can still distinguish the two images. 
% % This is because it is carried out linearly, for every column and row, 
% % this method performs far better. Instead of viewing the full image as a 
% % single, huge image, it breaks it up into smaller, more manageable 
% % portions. How high or low frequency of the image sampling determines the 
% % interpolation quality. 
% % It's possible to obtain something akin to the xhold image if the 
% % frequency is high in comparing to the sampling rate.
% 
% % Are edges low frequency or high frequency features?
% % Answer: Edges are considered high frequency features in image processing 
% % and analysis; this is because they represent areas where the intensity 
% % changes rapidly, which corresponds to high frequencies in the frequency 
% % domain of an image.
% % High frequency means rapid change, low frequency means smooth change.
% % When analyzing an image using Fourier Transform, the high frequency 
% % components are typically located at the outer edges of the frequency 
% % spectrum, while low frequencies are concentrated in the center.
% 
% % Are the fence posts low frequency or high frequency features? 
% % Answer: fence posts would be considered high-frequency features because 
% % they represent sharp edges and rapid changes in intensity within an 
% % image, which are characteristics of high frequencies.
% % Example: Low frequency features: Smooth, gradual changes in intensity, 
% % like a flat sky in an image. 
% % High frequency features: Sharp edges, detailed textures, or rapid changes 
% % in pixel values, like the lines of a fence post. 
% 
% % Is the background a low frequency or high frequency feature?
% % Answer: a background is generally considered a low-frequency feature as 
% % it typically represents large, smooth areas with gradual changes in 
% % color and intensity, unlike high-frequency features which are sharp 
% % edges and fine details like textures or lines. 
% % Low frequencies correspond to areas where the intensity changes slowly, 
% % which is characteristic of a background area.
% %  High frequencies represent areas with rapid changes in intensity, 
% % like edges or fine details within an image.
