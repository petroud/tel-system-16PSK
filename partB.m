clc;
close all;
clear all;
%----------------------------%
% Author: Petrou Dimitrios
% Year: 2023  
% TU of Crete
% Telecommunication Systems I
%----------------------------%


%% B1

N=100;
SNR=[-2:2:24];
T=10^(-2);
over=10;
Ts=T/over;
Fs=1/Ts;
A=4;
a=0.5;
F0=200;
K=1000;

%Generate SRRC pulse, all cycles should be 
%performed using the same SRRC filters.
[ph,t] = srrc_pulse(T,Ts,A,a);

disp('Please be patient...this might take a while');
%For all values of SNR set in the vector above, simulate
%the telecom system over and over again to calculate
%the Symbol Error Probability and Bit Error Probability
for cycle=1:length(SNR)
    
    symbolErrorsOfSamples = 0;
    bitErrorsOfSamples = 0;
    fprintf('Running %d experiments for cycle %d for SNR value=%ddB\n',K,cycle,SNR(cycle));
    
    for z=1:K
        bit_seq = (sign(randn(N,4))+1)/2;
        Xn = bits_to_PSK_16(bit_seq);
        Xi = Xn(:,1);
        Xq = Xn(:,2);
        
        %Calculate the deltas
        Xi_d = Fs*upsample(Xi,over);
        Xq_d = Fs*upsample(Xq,over);
        
        t_Xi_d = (0:Ts:N/(1/T)-Ts);
        t_Xq_d = (0:Ts:N/(1/T)-Ts);
        
        %Convolve with the SRRC filter
        Xi_t = conv(ph,Xi_d);
        Xq_t = conv(ph,Xq_d);
        
        t_Xi_t = (t(1)+t_Xi_d(1):Ts:t(end)+t_Xi_d(end));
        t_Xq_t = (t(1)+t_Xq_d(1):Ts:t(end)+t_Xq_d(end));
        
        %Attach the signals to phasors
        XI = 2 * (Xi_t) .* (cos(2*pi*F0*transpose(t_Xi_t)));
        XQ = -2 * (Xq_t) .* (sin(2*pi*F0*transpose(t_Xq_t)));
        
        %Create the channel input signal
        X = XI + XQ;
        
        %Generate Gaussian White Noise and attach to signal;
        sw2 = 1/(Ts*(10^(SNR(cycle)/10)));
        W = sqrt(sw2).*randn(length(X),1);
        %Calculate variance
        sigma2 = Ts*sw2/2;
        
        Y = X + W;
        
        %Generate Yi,Yq
        Yi_t = Y.*(cos(2*pi*F0*transpose(t_Xq_t)));
        Yq_t = Y.*(-sin(2*pi*F0*transpose(t_Xq_t)));

        %Convole the 2 signals with the SRRC filter
        YI_f = conv(ph, Yi_t);
        YQ_f = conv(ph, Yq_t);
        
        t_YI_f = (t(1)+t_Xq_t(1):Ts:t(end)+t_Xq_t(end));
        t_YQ_f = (t(1)+t_Xq_t(1):Ts:t(end)+t_Xq_t(end));
        
        %Tailcut the two convolutions
        Yi_cut = YI_f(80+1:(length(t_YI_f)-80));
        Yq_cut = YQ_f(80+1:(length(t_YQ_f)-80));
        
        %Downsample and produce Yi,n and Yq,n
        YI = downsample(Yi_cut,over);
        YQ = downsample(Yq_cut,over);
        
        %Create the derived Yn sequence of symbols
        Yn = [YI YQ];
        
        %Get the estimated bit sequence
        [est_X, est_bit_seq] = detect_PSK_16(Yn);
        
        symbolErrorsOfSamples = symbolErrorsOfSamples + symbol_errors(est_X, Xn);
        bitErrorsOfSamples = bitErrorsOfSamples + bit_errors(est_bit_seq, bit_seq);
    end
    
    %Update statistics 
    %Calculate probabilities of passed cycle
    SEP(1,cycle) = symbolErrorsOfSamples/(N*K);
    BEP(1,cycle) = bitErrorsOfSamples/(N*K*4);
    
    %Calculate smart upper bound for SEP and bit lower bound for BEP
    SUB(cycle) = 2*Q(1./sqrt(sigma2)*sin(pi/16));
    BLB(cycle) = SUB(cycle)/4;    
end


%Plot the outcomes of the experiment
figure('Name','B2');
semilogy(SNR, SUB);
hold on;
semilogy(SNR, SEP);
hold off;
xlabel('SNR(dB)');
ylabel('Symbol Error Probability (SEP)');
legend('Smart Upper Bound','Exp Est SEP');
title('Monte Carlo method for SEP estimation');

figure('Name','B3');
semilogy(SNR, BLB);
hold on;
semilogy(SNR, BEP);
hold off;
xlabel('SNR(dB)');
ylabel('Bit Error Probability (BEP)');
legend('Lower Bound','Exp Est BEP');
title('Monte Carlo method for BEP estimation');

