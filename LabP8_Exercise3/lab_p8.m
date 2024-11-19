% ========================================================================== 
% ECE 5530/ ECE 6530 Digital Signal Processing, ECE Dep., The University of Utah 
% Group Project
% Team: Kyle G. Gayliyev, Skylar Stockham, Eddie Franco
% ========================================================================== 
close all;
clear;

%% Sampling, Aliasing and Reconstruction : Lab P8 Exercise: 3 Lab Exercises 
 %% 3.1 Down-Sampling 
 %% Questions. 
 % Answer the following questions for the lighthouse picture,
 % downsampled by two.
 % a) Describe how the aliasing appears visually. Compare the original to 
 % the downsampled image. Which parts of the image show the aliasing 
 % effects most dramatically?
 % b) This part is challenging: explain why the aliasing happens in the 
 % lighthouse image by using a “frequency domain” explanation. 
 % In other words, estimate the frequency of the features that are 
 % being aliased. Give this frequency as a number in cycles per pixel. 
 % (Note that the fence provides a sort of “spatial chirp” where the 
 % spatial frequency increases from left to right.) Can you relate your 
 % frequency estimate to the Sampling Theorem?
 % You might try zooming in on a very small region of both the original 
 % and downsampled images.
%% Answers. 
% First, let's downsampler the lighthouse picture by two. To do so,
% we upload the "lighthouse.mat" given file. Then, we apply the following
% formula to download it and save the downloaded file into a new
% variable.

% Upload the "lighthouse.mat" given file
% lighthouse = importdata('lighthouse.mat');
load lighthouse.mat;

% Downsampling the lighthouse image by factor of 2
downsampled_lh = xx(1:2:end,1:2:end); % size : 163x213
whos;

%% 3.1a) 
% Let's draw the image so that we can see and compare the results.
% Plot the original lighthouse image 

figure('Name','Lighthouse Original Image');
clf;
imshow(xx);

% Plot the downsampled lighthouse image 
figure('Name','Lighthouse Downlampled Image');
clf;
imshow(downsampled_lh);

% Comparing both on the same figure
figure('Name','Comparing Lighthouse Original VS Downlampled Image');
clf;

subplot(1, 2, 1); imshow(xx); title('Original Image');
subplot(1, 2, 2); imshow(downsampled_lh); title('Downsampled Image');

% subplot(2,1,1);
% imshow(xx);
% subplot(2,1,2);
% imshow(downsampled_lh);

% Aliasing appears visually as jagged edges, stair-step patterns, or 
% a "pixelated" look, particularly along sharp transitions between 
% colors or high-contrast areas in the image.
% In this example, the downsampled image appears less smooth and more 
% distorted compared to the original lighthouse picture. Due to aliasing, 
% some part of the fence, and the lighthouse itself appear blurred or 
% distorted in the downsampled version. The aliasing effect is most 
% noticeable in areas with fine details, such as the roof shingles, 
% the fence, and the edges of the windows and doors.

%% 3.1b) 

% Frequency Domain Analysis
reg_original = xx(100:150, 100:150); %Extract the regions of interest
reg_downsampled = downsampled_lh(50:75, 50:75);

% Perform 2D Fourier Transform to analyze the frequency
fft_original = fft2(reg_original); 
fft_downsampled = fft2(reg_downsampled);

% Frequency Estimation and Comparation
figure('Name','Frequency Estimation and Comparation');
subplot(1, 2, 1); imagesc(abs(fftshift(fft_original))); title('Original Image Frequency');
subplot(1, 2, 2); imagesc(abs(fftshift(fft_downsampled))); title('Downsampled Image Frequency');

% The lighthouse image's aliasing is caused by the downsampled image's 
% inability to completely capture the high frequency information in the 
% original image. This can be explained in the frequency domain by the 
% fact that the frequency of the aliased features is more than half of the 
% sampling rate. This is called the Nyquist rate, and it says that the 
% sampling rate must to be more than twice the maximum frequency in the 
% original image to avoid aliasing.

% The lighthouse image's peak frequency is around one cycle per pixel 
% because the spatial frequency rises from left to right (the fence 
% creates a kind of "spatial chirp"). The maximum frequency that the 
% downsampled image can record is 0.5 cycles per pixel because it only 
% has half as many pixels as the original. Aliasing would occur at 
% frequencies higher than this.

% Because the downsampled image has a lower sampling rate (fewer pixels), 
% it is unable to capture the high frequency information seen in the 
% original image, which is why aliasing occurs. The sample Theorem, which 
% asserts that the sample rate must be more than twice the maximum 
% frequency in the signal in order to adequately capture it, is broken 
% by this.

% The reason for the aliasing is that the high frequency information 
% in the original image is not captured by the downsampled version. 
% About 150 cycles per pixel is the highest frequency found in the original 
% image. The sample Theorem states that the sample rate must be at least 
% twice the highest frequency in the image, or 2*150 = 300 cycles per 
% pixel, in order to prevent aliasing. Aliasing comes from the downsampled 
% image's inability to capture high frequency information due to its 
% sampling rate of only 150 cycles per pixel.

%% 3.2 Reconstruction of Images



% fprintf(' The minimum value of N =');
% disp(N);
