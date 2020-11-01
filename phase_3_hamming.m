%
% Data Assumptions
% 
% Amplitude of OOK is set to 5
% Hamming Code: Codewode 7 bits, Data 4 bits
% Sampling Frequency = 16 x Carrier Frequency
% 
% Generate data -> Hamming Code Encoding -> OOK/BPSK Modulation -> Noise -> Demodulation -> Comparison
% 

clear all; close all; clc;
carrier_freq = 10000; %10kHz
sample_freq = 16 * carrier_freq;
data_rate = 1000; %1kbps
data_length = 1024;
encoded_signal_length = 1792;
amp = 1;

% Low Pass 6th order Butterworth filter with 0.2 normalised cutoff freq
[b, a] = butter(6, 0.2);

% Time simulation
t = 0: 1/sample_freq : encoded_signal_length/data_rate;

% Carrier Signal Generation
carrier_signal = amp .* cos(2*pi*carrier_freq*t);

% Length of transmitted signal
signal_length = sample_freq*encoded_signal_length/data_rate + 1;

% SNR values to test
SNR_dB = 0:1:5;
SNR = (10.^(SNR_dB/10));

% Number of tests per SNR
test_samples = 20;

OOK_error_rate = zeros([length(SNR) 1]);
BPSK_error_rate = zeros([length(SNR) 1]);
        
% Generate Hamming encoded signals
data = round(rand(1,data_length));
hamming_signal= encode(data,7,4,'hamming/fmt');

signal = zeros(1, signal_length);
for k = 1: signal_length - 1
    signal(k) = hamming_signal(ceil(k*data_rate/sample_freq));
end
signal(signal_length) = signal(signal_length - 1);

% OOK Modulation
OOK_signal = carrier_signal .* signal;

% BPSK Modulation
BPSK_source_signal = signal .* 2 - 1;
BPSK_signal = carrier_signal .* BPSK_source_signal;

OOK_signal_power = (norm(OOK_signal)^2)/signal_length;
BPSK_signal_power = (norm(BPSK_signal)^2)/signal_length;

% For different SNR values, test over 20 samples
for i = 1 : length(SNR)
	OOK_average_error = 0;
    BPSK_average_error = 0;
    
	for j = 1 : test_samples
        
        % Generate noise
		noise_power = OOK_signal_power ./SNR(i);
		noise = sqrt(noise_power/2) .*randn(1,signal_length);
        
		% OOK Signal on Receiver's end
		OOK_received = OOK_signal+noise;
        
        % Start of OOK Detection
        OOK_squared = OOK_received .* 2.* carrier_signal;
        
        % Low Pass Filter
        OOK_filtered = filtfilt(b, a, OOK_squared);
        
		% BPSK Signal on Receiver's end
		BPSK_received = BPSK_signal+noise;
        
        % Non-coherent detection
        BPSK_squared = BPSK_received .* 2.* carrier_signal;
        % Low Pass Filter
        BPSK_output = filtfilt(b, a, BPSK_squared);
        
        % Demodulation by sample & threshold
        sample_period = sample_freq / data_rate;
        [OOK_sample, OOK_result] = sample_and_threshold(OOK_filtered, sample_period, (amp^2)/2, encoded_signal_length);
        [BPSK_sample, BPSK_result] = sample_and_threshold(BPSK_output, sample_period, 0, encoded_signal_length);
        

        OOK_decoded = decode(OOK_result,7,4,'hamming/fmt');
        BPSK_decoded = decode(BPSK_result,7,4,'hamming/fmt');
        
        OOK_error =  biterr(OOK_decoded, data) ./encoded_signal_length;
        BPSK_error = biterr(BPSK_decoded, data) ./encoded_signal_length;
        OOK_average_error = OOK_error + OOK_average_error;
        BPSK_average_error = BPSK_error + BPSK_average_error;
    end
    
	OOK_error_rate(i) = OOK_average_error / test_samples;
    BPSK_error_rate(i) = BPSK_average_error / test_samples;
    
%     Plot the 5db SNR signals
    if (SNR_dB(i) == 5)
        figure(2)
        subplot(2, 1, 1);
        plot(data, 'b');
        title("Original Data")
        xlim([0 2000])

        subplot(2, 1, 2);
        plot(hamming_signal, 'b');
        title("Hamming Encoded Data")
        xlim([0 2000])

        figure(4)
        subplot(4, 1, 1);
        spectrogram(OOK_signal,'yaxis')
        title("Transmitted OOK Modulated Signal")

        subplot(4, 1, 2);
        spectrogram(OOK_received,'yaxis')
        title("Received OOK Modulated Signal")

        subplot(4, 1, 3);
        plot(OOK_sample)
        title("OOK Demodulated Signal")

        subplot(4, 1, 4);
        plot(OOK_decoded)
        title("OOK Decoded Signal");

        figure(5)
        subplot(4, 1, 1);
        spectrogram(BPSK_signal,'yaxis')
        title("Transmitted BPSK Modulated Signal")

        subplot(4, 1, 2);
        spectrogram(BPSK_received,'yaxis')
        title("Received BPSK Modulated Signal")

        subplot(4, 1, 3);
        plot(BPSK_sample)
        title("BPSK Demodulated Signal")

        subplot(4, 1, 4);
        plot(BPSK_decoded)
        title("BPSK Decoded Signal");

        figure(6);
        subplot(3, 1, 1);
        plot(data);
        title("Original Data");
        xlim([0 1024]);
        ylim([0 1]);

        subplot(3, 1, 2);
        plot(OOK_decoded);
        title("OOK Decoded Data");
        xlim([0 1024]);

        subplot(3, 1, 3);
        plot(BPSK_decoded);
        title("BPSK Decoded Data");
        xlim([0 1024]);
    end
end

% Plot OOK vs DBSK bit error rate
figure(1)
p1 = semilogy(SNR_dB, OOK_error_rate,'r-*');
hold on
p2 = semilogy(SNR_dB, BPSK_error_rate, 'b-*');
hold off
ylabel('Bit Error Rate (BER)');
xlabel('SNR (dB)');
legend([p1(1) p2(1)],{'OOK','DPSK'})
xlim([0 50]);

