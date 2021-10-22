clc, clear all, close all;

%OOK modulation
N = 1024; %number of bits
data_rate = 1000; %bits per second
Tb = 1/data_rate; %bit period
Fc = 10000; %carrier frequency
Fs = 16 * Fc; %sampling frequency
signal_length = Fs * N/data_rate + 1; %length of carrier signal
samples_per_bit = Fs/data_rate; %160 samples per bit
SNR_dB = 0:5:50; %range of SNR values to take in db
SNR = (10.^(SNR_dB/10));
amplitude = 12;     %Carrier amplitude for binary input '1'
t = 0:1/Fs:N/data_rate;   %time interval for carrier signal
error_rate_OOK = zeros(length(SNR)); 

binary_data = round(rand(1,N)); 
data_stream = zeros(1, signal_length); 
index = 0;
for i = 1 : signal_length-1
    if mod(i-1, samples_per_bit) == 0
        index = index + 1;
    end
    data_stream(i) = binary_data(index);
end
data_stream(signal_length) = data_stream(signal_length - 1);

Carrier = amplitude .* cos(2*pi*Fc*t); 
signal_OOK = Carrier .* data_stream;
signal_power_OOK = (norm(signal_OOK)^2) / signal_length;  %Sum of squared signal amp over signal length

reruns = 10;

for j = 1 : length(SNR) 
    average_errorOOK = 0;
    
    for k= 1 : reruns
        noise_power = signal_power_OOK / SNR(j); %noise variance, calculated from S/N = SNR. noise power = noise variance = std^2
        noise = sqrt(noise_power) .* randn(1,signal_length) + 0; %Formula for Noise = std * randn() + mean
        noisy_OOK = signal_OOK + noise;
        squared_OOK = noisy_OOK .* noisy_OOK;
        [LPF_b ,LPF_a] = butter(6, 0.2); %LPF 
        filtered_OOK = filtfilt(LPF_b, LPF_a, squared_OOK);
        sampled_OOK = sample(filtered_OOK, samples_per_bit, N);
        decoded_OOK = decision_device(sampled_OOK, N, amplitude/2);  %OOK threshold is 0.5*(A+0)
    
        error_OOK = 0;
        
        for l = 1: N - 1
            if(decoded_OOK(l) ~= binary_data(l))
                error_OOK = error_OOK + 1;
            end
        end
        average_errorOOK = error_OOK + average_errorOOK;
    end
    error_rate_OOK(j) = (average_errorOOK / reruns)/N;
    
    if SNR_dB(j) == 15 %customizable SNR
        figure('Name', 'OOK Error');
        semilogy (SNR_dB, error_rate_OOK,'k-*');

        figure('Name','ASK Modulation and Demodulation','NumberTitle','off');
        grid on;
        subplot(5,1,1);
        plot(binary_data);
        xlabel('Time(Sec)');
        ylabel('Amplitude(Volts)');
        title('BINARY DATA WAVEFORM');
        
        subplot(5,1,2);
        plot(signal_OOK);
        xlabel('Time(Sec)');
        ylabel('Amplitude(Volts)');
        title('ASK MODULATED SIGNAL');

        subplot(5,1,3);
        plot(noisy_OOK);
        xlabel('Time(Sec)');
        ylabel('Amplitude(Volts)');
        title('RECEIVED SIGNAL');

        subplot(5,1,4);
        stem(sampled_OOK);
        xlabel('Time(Sec)');
        ylabel('Amplitude(Volts)');
        title('DEMODULATED SIGNAL');

        subplot(5,1,5);
        plot(decoded_OOK);
        xlabel('Time(Sec)');
        ylabel('Amplitude(Volts)');
        title('DECODED SIGNAL');
    end
end

function sampled = sample(x,sampling_period,num_bit)
    sampled = zeros(1, num_bit);
    for n = 1: num_bit
        sampled(n) = x((2 * n - 1) * sampling_period / 2);
    end
end

function binary_out = decision_device(sampled,num_bit,threshold)
    binary_out = zeros(1,num_bit);
    for n = 1:num_bit
        if(sampled(n) >= threshold)
            binary_out(n) = 1;
        else 
            binary_out(n) = 0;
        end
    end
end