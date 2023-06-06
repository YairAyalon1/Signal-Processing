clear all;
clc;
 
% Final work - Matlab code
% Made by: Yair Ayalon & Avraham Remer
 
%% Part 1
Fs=16000; % Sample frequency 
%% Audio recording adding (Recording_1)
[y,Fs] = audioread('Recording_1_original_for_matlab_code.wav'); % Getting the audio recording
 
%% Audio recording adding (Recording_2)
[z,Fs] = audioread('Recording_2_original_for_matlab_code.wav'); % Getting the audio recording
 
%% Adding Gaussian noise to first recording (y)
mean=600; % Mandatory noise requierements
sigma=200; % Mandatory noise requierements
N_1=length(y);
f_1=0:N_1-1;
gauss_distribution_1=normpdf(f_1,mean,sigma)+normpdf(f_1,N_1-mean,sigma); % Creating gaussian distribution in frequency axis 
gauss_noise_freq_1=rand(1,N_1).*gauss_distribution_1; % Creating gaussian noise running on top of gaussian distribution in frequency axis 
gauss_noise_time_1=ifft(gauss_noise_freq_1); % Transforming the noise to time axis
gauss_noise_time_normalized_1=gauss_noise_time_1/(var(real(gauss_noise_time_1))^0.5);
SNR=10; %SNR ratio
y_noise=zeros(1,length(y')); % Initial matrix represents the signal after adding gaussian noise
p_1=y';
for i=1:length(y')
    y_noise(i)=(2*real(gauss_noise_time_normalized_1(i))+p_1(i)); % Sum of the noise twice and the original recording (twice in order to be able to listen to the noise)    
end
 
%% Ploting recording signal and noised recording (y)
% subplot(2,1,1);
% plot(abs(fft(y)),'b');
% axis([0 10000 -0.04 800]);
% title('First recording - original');
% ylabel('Magnitude');
% subplot(2,1,2);
% plot(abs(fft(y_noise)),'r');
% axis([0 10000 -0.04 800]);
% title('First recording - with noise');
% ylabel('Magnitude');
% xlabel('f [Hz]');
 
%% Adding Gaussian noise to second recording (z)
mean=600; % Mandatory noise requierements
sigma=200; % Mandatory noise requierements
N_2=length(z);
f_2=0:N_2-1;
gauss_distribution_2=normpdf(f_2,mean,sigma)+normpdf(f_2,N_2-mean,sigma); % Creating gaussian distribution in frequency axis 
gauss_noise_freq_2=rand(1,N_2).*gauss_distribution_2; % Creating gaussian noise running on top of gaussian distribution in frequency axis 
gauss_noise_time_2=ifft(gauss_noise_freq_2); % Transforming the noise to time axis
gauss_noise_time_normalized_2=gauss_noise_time_2/(var(real(gauss_noise_time_2))^0.5);
SNR=10; %SNR ratio
z_noise=zeros(1,length(z')); % Initial matrix represents the signal after adding gaussian noise
p_2=z';
for i=1:length(z')
    z_noise(i)=(2*real(gauss_noise_time_normalized_2(i))+p_2(i)); % Sum of the noise twice and the original recording (twice in order to be able to listen to the noise)   
end
 
%% Ploting recording signal and noised recording (z)
% subplot(2,1,1);
% plot(abs(fft(z)),'b');
% axis([0 10000 -0.04 800]);
% title('Second recording - original');
% ylabel('Magnitude');
% subplot(2,1,2);
% plot(abs(fft(z_noise)),'r');
% axis([0 10000 -0.04 800]);
% title('Second recording - with noise');
% ylabel('Magnitude');
% xlabel('f [Hz]');
 
%% Echo creating 
echo_every_sec=1; % Echo dealay in seconds - every 1 sec in order to be able to hear the echo problem
echo_every_sec_in_samples=echo_every_sec*Fs; % Every _ samples the echo will apear
A= [1 ...
    zeros(1,echo_every_sec_in_samples-1)...
    -0.5...
    zeros(1,echo_every_sec_in_samples-1)...
    0.25...
    zeros(1,echo_every_sec_in_samples-1)...
    -0.125]; % Echo creating - every 1 seconds and 180 degrees phase
 
%% Recordings with gaussian noise && echo
y_echo=conv(A,y_noise); % First recording with echo & gaussian noise
z_echo=conv(A,z_noise); % Second recording with echo & gaussian noise
 
%% LMS filter design (first recording)
N=length(y');
order=20; %Filter oreder
mu=1;
h=zeros(order,1);
h_mat=zeros(order,N-order);
err=zeros(1,N-order);
for i=order:(N-order) 
    x=(y_noise(i:-1:i-order+1))';
    k_1(i)=h'*x;
    d=y(i);
    err(i)=d-k_1(i);
    h=h+(mu*conj(err(i))*x)/(x'*x);
    h_mat(:,i)=h;
end 
% k_1 is the filtered signal
 
%% LMS filter design (second recording)
N=length(z');
order=20; %Filter oreder
mu=1;
h=zeros(order,1);
h_mat=zeros(order,N-order);
err=zeros(1,N-order);
for i=order:(N-order) 
    x=(z_noise(i:-1:i-order+1))';
    k_2(i)=h'*x;
    d=z(i);
    err(i)=d-k_2(i);
    h=h+(mu*conj(err(i))*x)/(x'*x);
    h_mat(:,i)=h;
end
% k_2 is the filtered signal
 
%% Recordings playing (first - y)
% sound(y,Fs); % Playing first signal - original
% sound(y_noise,Fs); % Playing first signal with gaussian noise
% sound(y_echo,Fs); % Playig first signal after adding gaussian noise && echo
% sound(k_1,Fs); % Playing first signal after LMS filter
 
%% Recordings playing (second - z)
% sound(z,Fs); % Playing second signal - original
% sound(z_noise,Fs); % Playing second signal with gaussian noise
% sound(z_echo,Fs); % Playig second signal after adding gaussian noise && echo
% sound(k_2,Fs); % Playing second signal after LMS filter
 
%% Part 2
%% First signal modulation by cosinos wave
t=0:length(y')-1;
y_m=y'.*cos(2*pi*150*t); % First recording after modulation
% sound(y_m,Fs); % Playing first recording after modulation by cosinos wave
 
%% Highpass filter 
y_m_f=highpass(y_m,200,Fs); % Highpass filtering to first recording after modulation
% sound(y_m_f,Fs); % Playing first modulated recording after highpass filter
 
%% Adding gauss noise to y_m_f
mean=600; % Mandatory noise requierements
sigma=200; % Mandatory noise requierements
N_3=length(y_m_f);
f_3=0:N_3-1;
gauss_distribution_3=normpdf(f_3,mean,sigma)+normpdf(f_3,N_3-mean,sigma); % Creating gaussian distribution in frequency axis 
gauss_noise_freq_3=rand(1,N_3).*gauss_distribution_3; % Creating gaussian noise running on top of gaussian distribution in frequency axis 
gauss_noise_time_3=ifft(gauss_noise_freq_3); % Transforming the noise to time axis
gauss_noise_time_normalized_3=gauss_noise_time_3/(var(real(gauss_noise_time_3))^0.5);
SNR=10; %SNR ratio
y_m_f_n=zeros(1,length(y_m_f')); % Initial matrix represents the signal after adding gaussian noise
p_3=y_m_f';
for i=1:length(y_m_f')
    y_m_f_n(i)=(2*real(gauss_noise_time_normalized_3(i))+p_3(i)); % Sum of the noise twice and the original recording (twice in order to be able to listen to the noise)   
end
% sound(y_m_f_n,Fs); % Playing y_m_f signal with noise
 
%% Adding echo to y_m_f_n
y_m_f_n_e=conv(A,y_m_f_n);
% sound(y_m_f_n_e,Fs); % Plating y_m_f_n signal with echo problem
 
%% RLS filter design
[y,Fs] = audioread('Recording_1_original_for_matlab_code.wav'); % Getting the original audio recording 
Order=30; % Filter order
Delta=0.1; % Initial input covariance estimate
P0=(1/Delta)*eye(Order,Order); % Initial setting for the P matrix
rlsfilt=dsp.RLSFilter(Order,'InitialInverseCovariance',P0);
for J=1:10
    b=y_m_f_n_e(); % Noise
    S=y();
    D=[[y_m_f,zeros(1,48000)]+[S',zeros(1,48000)]]; % Fixing array dimentions by adding zeros vector
    [a,e]=rlsfilt(b,D);
end
% sound(e,Fs); % Playing first recording after RLS filter
 
%% Ploting - part 2
% subplot(3,1,1);
% plot(abs(fft(y_m)),'b');
% axis([0 10000 -0.04 800]);
% title('First recording - modulated');
% ylabel('Magnitude');
% subplot(3,1,2);
% plot(abs(fft(y_m_f)),'r');
% axis([0 10000 -0.04 800]);
% title('First recording - after HPF');
% ylabel('Magnitude');
% subplot(3,1,3);
% plot(abs(fft(y_m_f_n)),'g');
% title('First recording - after HPF with noise');
% ylabel('Magnitude');
% axis([0 10000 -0.04 800]);
% xlabel('f [Hz]');
 
%% Conclusion of the filters functionality (time axis)
%% Ploting signal before and after LMS filter (part 1 - first recording)
% subplot(3,1,1);
% plot((y),'b');
% axis([0 123904 -1 1]);
% title('First recording - original');
% ylabel('Magnitude');
% subplot(3,1,2);
% plot(y_noise,'r');
% axis([0 123904 -1 1]);
% title('First recording - with noise');
% ylabel('Magnitude');
% subplot(3,1,3);
% plot(k_1,'g');
% title('First recording - after LMS filter');
% ylabel('Magnitude');
% axis([0 123904 -1 1]);
% xlabel('N (samples)');

%% Ploting signal before and after LMS filter (part 1 - Second recording)
% subplot(3,1,1);
% plot((z),'b');
% axis([0 128000 -1 1]);
% title('Second recording - original');
% ylabel('Magnitude');
% subplot(3,1,2);
% plot(z_noise,'r');
% axis([0 128000 -1 1]);
% title('Second recording - with noise');
% ylabel('Magnitude');
% subplot(3,1,3);
% plot(k_2,'g');
% title('Second recording - after LMS filter');
% ylabel('Magnitude');
% axis([0 128000 -1 1]);
% xlabel('N (samples)');

%% Ploting signal before and after RLS filter (part 2)
% subplot(3,1,1);
% plot((y),'b');
% axis([0 123904 -1 1]);
% title('First recording - original');
% ylabel('Magnitude');
% subplot(3,1,2);
% plot(y_m_f_n,'r');
% axis([0 123904 -1 1]);
% title('First recording - after HPF & with noise');
% ylabel('Magnitude');
% subplot(3,1,3);4
% plot(e,'g');
% title('First recording - after RLS filter');
% ylabel('Magnitude');
% axis([0 123904 -1 1]);
% xlabel('N (samples)');