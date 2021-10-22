clear all; close all; clc;

%define carrier frequency and sampling frequency
Fc = 10000;                                             %carrier frequency
Fs = 16*Fc;                                             %16 times oversampled(sampling frequency) -> sample/sec

%define signal sender requirement
N = 1024;
dataRate = 1000;
samplePerBit = Fs/dataRate;                             %(16*10000)/1000, 160 sample/bit
amplitude = 12;                                         %assume the amplitude=12

%time
Ts = 1/Fs;                                              %period of sampling(time between two samples)
T = N/dataRate;                                         %total time needed for the whole data
t = 0: Ts: T-Ts;                                        %time interval for sampling signal

%define filter
[b_low, a_low] = butter(6,0.2);                         %6th order LP butterworth filter with 0.2 normalized cutoff frequency
[b_high, a_high] = butter(6,0.2, 'high');               %6th order HP butterworth filter with 0.2 normalized cutoff frequency

%carrier signal
carrier = amplitude .* cos(2*pi*Fc*t);                  %generate carrier frequency signal
carrier_2 = amplitude .* cos(2*pi*2*Fc*t);              %generate a higher carrier frequency signal for FSK

sampleLength = T/Ts;                                    %total number of samples

%SNR
SNR_db = 0:5:50;                                        %Different SNR values from 0 to 50 db
SNR = (10.^(SNR_db/10));                                %SNR_dB = 10 log (Signal_Power/Noise_Power)        
                                                        %SNR = Signal_Power/Noise_Power = 10^(SNR_dB/10)

%modify the variable below to choose at which SNR value to plot
plot_SNR_db = 5;                                        %choose SNR in dB to plot the result

%array for error result
errorRateOOK = zeros(length(SNR));
errorRateFSK = zeros(length(SNR));
errorRateBPSK = zeros(length(SNR));

%generate 1024 bits of data
data = round(rand(1,N));                                %1024-bits data generation

%sample data stream
sampleDataStream = zeros(1, sampleLength);              %generate the sample data stream
for i = 1:sampleLength
    sampleDataStream(i) = data(ceil(i/samplePerBit));
end

for i = 1:length(SNR)
    %------------OOK-------------%

    %OOK signal generation
    signalOOK = carrier .* sampleDataStream;            %generate OOK signal with carrier
    
    %generate noise for OOK signal
    signalPowerOOK = (norm(signalOOK)^2)/sampleLength;  %sum of squared signal amp over signal length
        noisePowerOOK = signalPowerOOK ./SNR(i);
        noiseOOK = sqrt(noisePowerOOK/2) .* randn(1,sampleLength);
    
    %OOK and noise transmission
        receiveOOK = signalOOK + noiseOOK;
    
    %detection of signal using square law device
    squaredOOK = receiveOOK .* receiveOOK;

    %low pass filter
    filteredOOK = filtfilt(b_low, a_low, squaredOOK);

    %sample and decision device, OOK threshold is 0.5*(A+0)
    demodulatedOOK = decisionDevice(filteredOOK, sampleLength, amplitude/2); 
    resultOOK = sample(demodulatedOOK, samplePerBit, N);

    %------------BPSK-------------%
    dataStreamBPSK = sampleDataStream .*2;
    signalBPSK = carrier .*dataStreamBPSK;

    %generate noise
    signalPowerBPSK = (norm(signalBPSK)^2)/sampleLength;
        noisePowerBPSK = signalPowerBPSK ./SNR(i);
        noiseBPSK = sqrt(noisePowerBPSK/2) .*randn(1,sampleLength);

    %transmission
        receiveBPSK = signalBPSK + noiseBPSK;

    %non-coherent detection --square law detection
    squaredBPSK = receiveBPSK .* carrier;  
    
    %output filter
    filteredBPSK = filtfilt(b_low, a_low, squaredBPSK);

    demodulatedBPSK = sample(filteredBPSK, samplePerBit, N);
    resultBPSK = decisionDevice(demodulatedBPSK, N, 0);


    %------------FSK-------------%
    %--------------Modulation--------------%
    inverseDataFSK = ~sampleDataStream;                 %bit flip
                                                        %bit-1-->High Carrier Frequency Signal, bit-0-->Low Carrier Frequency Signal
    signalFSK = carrier_2 .* sampleDataStream + carrier .* inverseDataFSK;

    %generate noise
    signalPowerFSK = (norm(signalFSK)^2)/sampleLength;
        noisePowerFSK = signalPowerFSK ./SNR(i);
        noiseFSK = sqrt(noisePowerFSK/2) .*randn(1,sampleLength);

    %transmission
        receiveFSK = signalFSK + noiseFSK;
    
    %--------------Demodulation--------------%
    %Non-coherent Detection: bandpass filters
    receiveFSKHigh = filtfilt(b_high,a_high,receiveFSK); %high pass filter
    receiveFSKLow = filtfilt(b_low,a_low,receiveFSK);    %low pass filter
    
    %Square Law
    squaredFSKLow = receiveFSKLow.*receiveFSKLow;
    squaredFSKHigh = receiveFSKHigh.*receiveFSKHigh;
    squaredFSK = squaredFSKHigh - squaredFSKLow;        %Merge both high and low signal
    
    %sample and decision device
    demodulatedFSK = decisionDevice(squaredFSK,sampleLength, 0);    %thresold = 0. 
    resultFSK = sample(demodulatedFSK,samplePerBit,N);

    %--------------calculation of error--------------%
    errorOOK = 0;
    errorBPSK = 0;
    errorFSK = 0;

    for k = 1: N
        if(resultOOK(k) ~= data(k))
            errorOOK = errorOOK + 1;
        end
        if(resultFSK(k) ~= data(k))
            errorFSK = errorFSK + 1;
        end
        if(resultBPSK ~= data(k))
            errorBPSK = errorBPSK + 1;
        end
    end

    %----------plotting variables---------%
    if(plot_SNR_db == SNR_db(i))
        plot_signal = data;

        plot_mod_OOK = signalOOK;
        plot_receive_OOK = receiveOOK;
        plot_demod_OOK = demodulatedOOK;
        plot_decoded_OOK = resultOOK;

        plot_mod_BPSK = signalBPSK;
        plot_receive_BPSK = receiveBPSK;
        plot_demod_BPSK = demodulatedBPSK;
        plot_decoded_BPSK = resultBPSK;

        plot_mod_FSK = signalFSK;
        plot_receive_FSK = receiveFSK;
        plot_demod_FSK = demodulatedFSK;
        plot_decoded_FSK = resultFSK;
    end

    errorRateOOK(i) = (errorOOK)/N;
    errorRateBPSK(i) = (errorBPSK)/N;
    errorRateFSK(i) = (errorFSK)/N;
end

%Error plot
figure(1);
semilogy(SNR_db, errorRateOOK, 'k-*');
hold on
semilogy(SNR_db, errorRateBPSK, 'c-*');
hold off
hold on
semilogy(SNR_db, errorRateFSK, 'r-*');
hold off
title('Error rate of OOK and BFSK for different SNR');
legend('OOK', 'BPSK','FSK');
ylabel('Pe');
xlabel('Eb/No')

%OOK
figure(2)
subplot(511);plot(plot_signal);title('Data Waveform');
subplot(512);plot( plot_mod_OOK,'k');title('Modulated Signal OOK');
subplot(513);plot(plot_receive_OOK, 'k');title('Received Signal OOK');
subplot(514);plot(plot_demod_OOK, 'k');title('Demodulated Signal OOK');
subplot(515);plot(plot_decoded_OOK);title('Decoded Signal');

%BPSK
figure(3)
subplot(511);plot(plot_signal);title('Data Waveform');
subplot(512);plot( plot_mod_BPSK,'k');title('Modulated Signal BPSK');
subplot(513);plot(plot_receive_BPSK, 'k');title('Received Signal BPSK');
subplot(514);plot(plot_demod_BPSK, 'k');title('Demodulated Signal BPSK');
subplot(515);plot(plot_decoded_BPSK);title('Decoded Signal');

%FSK
figure(4)
subplot(511);plot(plot_signal);title('Data Waveform');
subplot(512);plot( plot_mod_FSK,'k');title('Modulated Signal FSK');
subplot(513);plot(plot_receive_FSK, 'k');title('Received Signal FSK');
subplot(514);plot(plot_demod_FSK, 'k');title('Demodulated Signal FSK');
subplot(515);plot(plot_decoded_FSK);title('Decoded Signal');


%%--HELPER FUNCTION--%%
function desampled = sample(x,num_sample_per_bit,num_bit)
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
function binary_out = decisionDevice(sampled,num_bit,threshold)
    binary_out = zeros(1,num_bit);
    for n = 1:num_bit
        if(sampled(n) > threshold)
            binary_out(n) = 1;
        else 
            binary_out(n) = 0;
        end
    end
end











