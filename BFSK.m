clear all; close all; clc;

% BFSk

Fc = 10000; %carier frequency = 10kHz

Fs = 16 * Fc; % 16 times oversampled -> sample/sec

Ts=1/Fs; % period of sampling

data_rate = 1000; % data rate = 1kbps

N = 1024; % number of data bits

sample_per_bit = Fs / data_rate; % 160 sample/bit

amplitude = 12;

T=N/data_rate; % total time for the whole data

t = 0: Ts : T-Ts; % time interval for sampling signal

[b_low,a_low] = butter(6, 0.2); % 6th order LP butterworth filter with 0.2 normalized cutoff frequency

[b_high,a_high] = butter(6, 0.2, 'high'); % 6th order HP butterworth filter with 0.2 normalized cutoff frequency


Carrier = amplitude .* cos(2*pi*Fc*t); % generate carrier frequency signal
Carrier_2 = amplitude .* cos(2*pi*2*Fc*t); % generate a higher carrier frequecy signal for FSK


sample_length = T/Ts; % total number of samples
            
SNR_dB = 0:5:50; % Different SNR values from 0 dB to 50 dB (in multiples of 5 dB)
% SNR_dB = 10 log (Signal_Power/Noise_Power)        
% SNR = Signal_Power/Noise_Power = 10^(SNR_dB/10)
SNR = (10.^(SNR_dB/10));

% MODIFY THE VARIABLE BELOW TO CHOOSE AT WHICH SNR VALUE 
plot_SNR_dB = 5; % choose the SNR in dB to plot the result


% create arrays to store the error result
% Error_RateOOK = zeros(length(SNR));
% Error_RateBPSK = zeros(length(SNR));
Error_RateFSK = zeros(length(SNR)); 

data = round(rand(1,N)); % 1024-bit data generation

sample_data = zeros(1, sample_length); % create the sample data stream

% fill the data stream
for i = 1: sample_length
    sample_data(i) = data(ceil(i/sample_per_bit));
end

%for each SNR value
for i = 1 : length(SNR)

        %--FSK--%
        % --Modulation--
        inv_sample_data = ~sample_data; % bit flip 
        % bit-1-->High Carrier Frequency Signal, bit-0-->Low Carrier Frequency Signal
        signal_FSK = Carrier_2 .* sample_data + Carrier.*inv_sample_data; 
        
        % generate noise
        signal_power_FSK = (norm(signal_FSK)^2)/sample_length;
		noise_power_FSK = signal_power_FSK ./SNR(i);
		noiseFSK = sqrt(noise_power_FSK/2) .*randn(1,sample_length);
        
        %transmission
		receiveFSK = signal_FSK+noiseFSK;
        
        % --Demodulation--
        %Non-coherent Detection: bandpass filters
        receiveFSK_HIGH = filtfilt(b_high,a_high,receiveFSK); % high pass filter
        receiveFSK_LOW = filtfilt(b_low,a_low,receiveFSK); % low pass filter
        %Square Law
        squaredFSK_LOW = receiveFSK_LOW.*receiveFSK_LOW;
        squaredFSK_HIGH = receiveFSK_HIGH.*receiveFSK_HIGH;
        squaredFSK = squaredFSK_HIGH - squaredFSK_LOW; % Merge both high and low signal
        %sample and decision device
        % sampledFSK = sample(squaredFSK, sample_per_bit, N);% bit-1 is positive and bit-0 is negative
        demodFSK = decision_device(squaredFSK,sample_length, 0);  % thresold = 0. 
        resultFSK=desample(demodFSK,sample_per_bit,N);
        
        

        %--Calculate Error--%
        %ErrorOOK = 0;
        %ErrorBPSK = 0;
        ErrorFSK = 0;
        
        for k = 1: N
            if(resultFSK(k) ~= data(k))
                ErrorFSK = ErrorFSK + 1;
            end
        end

    
    %Plotting variables
    if (plot_SNR_dB == SNR_dB(i))
            plot_signal = data;
            plot_mod_FSK = signal_FSK;
            plot_receive_FSK = receiveFSK;
            plot_demod_FSK = demodFSK;
            plot_decoded_FSK = resultFSK;
    end
    
%	Error_RateOOK(i) = (ErrorOOK)/nBits;
  %  Error_RateBPSK(i) = (ErrorBPSK)/nBits;
    Error_RateFSK(i) = (ErrorFSK)/N;
end


%Error plot
figure(1);
semilogy(SNR_dB, Error_RateFSK, 'r-*');
title('Error rate of BFSK for different SNR');
%legend('FSK');
ylabel('Pe');
xlabel('Eb/No')

%FSK
figure(2)
subplot(511);plot(plot_signal);title('Data Waveform');
subplot(512);plot( plot_mod_FSK,'k');title('Modulated Signal FSK');
subplot(513);plot(plot_receive_FSK, 'k');title('Received Signal FSK');
subplot(514);plot(plot_demod_FSK, 'k');title('Demodulated Signal FSK');
subplot(515);plot(plot_decoded_FSK);title('Decoded Signal');



%%--HELPER FUNCTION--%%
function desampled = desample(x,num_sample_per_bit,num_bit)
    desampled = zeros(1, num_bit);
    length(x);
    for n = 1: num_bit
        start=(1+(n-1)*num_sample_per_bit);
        end_=(start-1+num_sample_per_bit);
        a=x(start:end_);
        desampled(n) = median(a);
    end
end


%This function simulates the decision device
function binary_out = decision_device(sampled,num_bit,threshold)
    binary_out = zeros(1,num_bit);
    for n = 1:num_bit
        if(sampled(n) > threshold)
            binary_out(n) = 1;
        else 
            binary_out(n) = 0;
        end
    end
end