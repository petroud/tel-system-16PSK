clc;
close all;
clear all;
%----------------------------%
% Author: Petrou Dimitrios
% Year: 2023  
% TU of Crete
% Telecommunication Systems I
%----------------------------%


%% Constant Declaration
N = 100;
Nf = 2048;
T = 10^(-2);
over = 10;
Ts = T/over;
Fs = 1/Ts;
A = 4;
a = 0.5;
F0 = 200;

%% A1
bit_seq = (sign(randn(N,4))+1)/2;
%% A2
Xn = bits_to_PSK_16(bit_seq);

XI=Xn(:,1);
XQ=Xn(:,2);

%% A3
%Generate SRRC pulse
[ph,t] = srrc_pulse(T,Ts,A,a);
F = -Fs/2:Fs/Nf:Fs/2-Fs/Nf;

%Create the X_delta signals for XI and XQ, to be used for convolution.
XI_d = Fs * upsample(XI, over);
t_XI_d = (0:Ts:N/(1/T)-Ts);

XQ_d = Fs * upsample(XQ, over);
t_XQ_d = (0:Ts:N/(1/T)-Ts);

%Convolve the delta signals with the pulse and create the time axis.
t_XIt = (t(1)+t_XI_d(1):Ts:t(end)+t_XI_d(end));
t_XQt = (t(1)+t_XQ_d(1):Ts:t(end)+t_XQ_d(end));

XI_t = conv(ph, XI_d)*Ts;
XQ_t = conv(ph, XQ_d)*Ts; 
t_total_I = length(t_XIt)*Ts;
t_total_Q = length(t_XQt)*Ts;

%Perform Fourier Transform and calculate periodgram
FX_i = fftshift(fft(XI_t,Nf)*Ts);
PX_i = (abs(FX_i).^2)/t_total_I;

FX_q = fftshift(fft(XQ_t,Nf)*Ts);
PX_q = (abs(FX_q).^2)/t_total_Q;


%Plot the output and the periodgrams 
figure('Name', 'A1 A2 A3 (a)');
subplot(3,1,1);
plot(t_XI_d, XI_d);
title('X_{i,n}');

subplot(3,1,2);
plot(t_XIt, XI_t);
title('X_{i,n} filterd using SRRC pulse');
xlabel('t')

subplot(3,1,3);
plot(F, PX_i);
title('Periodgram of X_{i,n} filtered using SRRC pulse');
ylabel('P_x(F)');
xlabel('Frequency(Hz)');


figure('Name', 'A1 A2 A3 (b)');
subplot(3,1,1);
plot(t_XQ_d, XQ_d);
title('X_{Q,n}');

subplot(3,1,2);
plot(t_XQt, XQ_t);
title('X_{Q,n} filterd using SRRC pulse');
xlabel('t')

subplot(3,1,3);
plot(F, PX_q);
title('Periodgram of X_{Q,n} filtered using SRRC pulse');
ylabel('P_x(F)');
xlabel('Frequency(Hz)');

%% A4

%Create XI(t) and XQ(t) by multiplying the filter outputs with appropriate
%phasors
XI = 2*(XI_t).*(cos(2*pi*F0*transpose(t_XIt)));
XI_total = length(t_XIt)*Ts;

XQ = -2*(XQ_t).*(sin(2*pi*F0*transpose(t_XQt)));
XQ_total = length(t_XQt)*Ts;

%Perform Fourier Transform and produce periodgrams
FXI = fftshift(fft(XI,Nf)*Ts);
PXI = (abs(FXI).^2)/XI_total;

FXQ = fftshift(fft(XQ,Nf)*Ts);
PXQ = (abs(FXQ).^2)/XQ_total;

%Plot the the XI,XQ and their periodgrams
figure('Name','A4')
subplot(4,1,1);
plot(t_XIt, XI);
title('X_I(t)');
xlabel('t');
subplot(4,1,2);
plot(F, PXI);
title('Periodgram of X_I(t)');
ylabel('P_X(F)');
xlabel('Frequency (Hz)');
subplot(4,1,3);
plot(t_XQt, XQ);
title('X_Q(t)');
xlabel('t');
subplot(4,1,4);
plot(F, PXQ);
title('Periodgram of X_Q(t)');
ylabel('P_X(F)');
xlabel('Frequency (Hz)');

%% A5
%Come up with the channel input X(t)
Xt = XI + XQ;
T_total = length(t_XQt)*Ts
%Perform Fourier Transform and compute Periodgram
Fxt = fftshift(fft(Xt,Nf)*Ts);
Pxt = (abs(Fxt).^2)/T_total;

figure('Name','A5');
subplot(2,1,1);
plot(t_XQt,Xt);
title('X = XI + XQ');
xlabel('t');
subplot(2,1,2);
plot(F,Pxt);
title('Periodgram for X');
ylabel('Px(F)');
xlabel('Frequency(Hz)');

%% Á7

%Create the Gaussian noise W(t)
SNR = 20;
sw2 = 1/(Ts*(10^(SNR/10)));
W = sqrt(sw2).*randn(length(Xt),1);

Y = Xt + W;

%Plot the signal containing the noise
figure('Name','A7');
plot(t_XQt, Xt);
hold on;
plot(t_XQt, Y);
hold off;
title('X(t) before and after Gaussian noise');
xlabel('Time');
ylabel('Magnitude');
legend('X(t)','Y(t) = X(t) + W(t)');

%% A8 

%Those are the YI,YQ signals coming from the multiplication of Yt with
%appropriate phasors
YI = Y.*(cos(2*pi*F0*transpose(t_XQt)));
YQ = Y.*(-sin(2*pi*F0*transpose(t_XQt)));

t_YI_total = length(t_XQt)*Ts;
t_YQ_total = length(t_XQt)*Ts;

%Perform Fourier Transform on them and produce periodgrams
FYI = fftshift(fft(YI,Nf)*Ts);
FYQ = fftshift(fft(YQ,Nf)*Ts);

PYI = (abs(FYI).^2)/t_YI_total;
PYQ = (abs(FYQ).^2)/t_YQ_total;

figure('Name','A8');
subplot(2,2,1);
plot(t_XQt,YI);
title('Y_I(t)');
subplot(2,2,3);
plot(F, PYI);
title('Periodgram of Y_I(t)');
ylabel('P_Y(F)');
xlabel('Frequency(Hz)');
subplot(2,2,2);
plot(t_XQt,YQ);
title('Y_Q(t)');
subplot(2,2,4);
plot(F, PYQ);
title('Periodgram of Y_Q(t)');
ylabel('P_Y(F)');
xlabel('Frequency(Hz)');


%% A9
%Using the SRRC pulse from A3 we convolve the the YI,YQ signals with it.
YI_2 = conv(ph, YI)*Ts;
YQ_2 = conv(ph, YQ)*Ts;
%Specify the correct time axis
t_YI_2 = (t(1)+t_XQt(1):Ts:t(end)+t_XQt(end));
t_YQ_2 = (t(1)+t_XQt(1):Ts:t(end)+t_XQt(end));

%Tail cut the filtered signal by 80 indexes
YI_cut = YI_2(80+1:(length(t_YI_2)-80));
t_YI_cut = t_YI_2(80+1:(length(t_YI_2)-80));

YQ_cut = YQ_2(80+1:(length(t_YQ_2)-80));
t_YQ_cut = t_YQ_2(80+1:(length(t_YQ_2)-80));

%Perform Fourier Transform on the tailcut signals
FYI_2 = fftshift(fft(YI_cut,Nf)*Ts);
FYQ_2 = fftshift(fft(YQ_cut,Nf)*Ts);

tyi_2_total = length(t_YI_cut)*Ts;
tyq_2_total = length(t_YQ_cut)*Ts;

PYI_2 = (abs(FYI_2).^2)/tyi_2_total;
PYQ_2 = (abs(FYQ_2).^2)/tyq_2_total;

%Plot the signals and their periodgrams
figure('Name','A9');
subplot(2,2,1);
plot(t_YI_cut,YI_cut);
title('Y_I(t) filtered with SRRC pulse');
subplot(2,2,3);
plot(F, PYI_2);
title('Periodgram of Y_I(t) filtered with SRRC pulse');
ylabel('P_Y(F)');
xlabel('Frequency(Hz)');
subplot(2,2,2);
plot(t_YQ_cut,YQ_cut);
title('Y_Q(t) filtered with SRRC pulse');
subplot(2,2,4);
plot(F, PYQ_2);
title('Periodgram of Y_Q(t) filtered with SRRC pulse');
ylabel('P_Y(F)');
xlabel('Frequency(Hz)');

%% A10
%Downsample and create the default Yi,n and Yq,n
YI_down = downsample(YI_cut, over);
YQ_down = downsample(YQ_cut, over);

%Create Yn
Yn=[YI_down, YQ_down];
%Scatterplot 
scatterplot(Yn);

%% A11
[est_X, est_bit_seq] = detect_PSK_16(Yn);
%% A12
symbolErrors = symbol_errors(est_X,Xn);
%% A13
bitErrors = bit_errors(est_bit_seq, bit_seq);

