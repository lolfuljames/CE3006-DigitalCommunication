clear all; close all; clc;

% Codeword length (n) & Message length (k)
codeword_length = 7;
message_length = 4;

% Data length (1024) & Encoded data length (1024/4 * 7 = 1792)
signal_length = 1024;
encoded_signal_length = signal_length/message_length * codeword_length;

% Carrier parameters and signal (amplitude = 5) as required
carrier_freq = 10000;
sampling_freq = 16 * carrier_freq;
data_rate = 1000;
t = 0: 1/sampling_freq : encoded_signal_length/data_rate;
carrier_signal = 5 .* cos(2*pi*carrier_freq*t);

% Sampled signal length
sampled_signal_length = sampling_freq*encoded_signal_length/data_rate + 1;

% Number of samples
nb_samples = 20;
% Low-pass filter - 6th order, 0.2 cutoff frequency
[b, a] = butter(6, 0.2);
% Signal power  
signal_power = 1; 
% maximum SNR
MAX_dB = 50;

% Generator and Syndrome Table for Cyclic Code (Channel Encoding-Decoding)
% Generator polynomial & Parity-check matrix for cyclic encoding
genpoly = cyclpoly(codeword_length, message_length);
parmat = cyclgen(codeword_length, genpoly);
% Syndrome table for cyclic decoding
trt = syndtable(parmat);


% Generate SNR from - 0 - 50 dB with conversion
SNR = generate_SNR(MAX_dB, 5);

OOK_error_rate = zeros([length(SNR) 1]);
BPSK_error_rate = zeros([length(SNR) 1]);

% Original signal 
original_signal = round(rand(1,signal_length));
% Cyclic encoding
encoded_signal = encode(original_signal, codeword_length, message_length, 'cyclic/binary', genpoly);

% Sampling
sampled_signal = zeros(1, sampled_signal_length);
for k = 1: sampled_signal_length - 1
    sampled_signal(k) = encoded_signal(ceil(k*data_rate/sampling_freq));
end
sampled_signal(sampled_signal_length) = sampled_signal(sampled_signal_length - 1);

% Modulation: On-Off Keying
OOK_signal = carrier_signal .* sampled_signal;

% Modulation: Binary Phase Shift Keying 
BPSK_source_signal = sampled_signal .* 2 - 1; % put to -1 +1
BPSK_signal = carrier_signal .* BPSK_source_signal;

OOK_signal_power = (norm(OOK_signal)^2)/sampled_signal_length;
BPSK_signal_power = (norm(BPSK_signal)^2)/sampled_signal_length;

% For each value of SNR, test of 20 samples.
for i = 1 : length(SNR) 
    
	OOK_average_error = 0;
	BPSK_average_error = 0;
    result = zeros(1, nb_samples);
    
    for j = 1 : nb_samples
        
        % Generate White Gaussian Channel Noise for OOK and BPSK
        noise_power = OOK_signal_power ./ SNR(i);
        noise = sqrt(noise_power/2) .*randn(1,sampled_signal_length);
        
        % Received Signal with added Channel Noise for OOK and BPSK
        OOK_received = OOK_signal + noise;
        BPSK_received = BPSK_signal + noise;
        
        % Non-Coherent Detection: OOK (Lecture notes 04 - pg 18)
        % Squared
        OOK_squared = OOK_received .* 2 .* carrier_signal;
        % Low-Pass Filter
        OOK_filtered = filtfilt(b, a, OOK_squared);
  
        % Non-Coherent Detection: BPSK (Lecture notes 04 - pg 50)
        % Squared
        BPSK_squared = BPSK_received .* 2 .* carrier_signal;
        % Low-Pass Filter
        BPSK_output = filtfilt(b, a, BPSK_multiplied);
        
        % Demodulation: Sampling and Threshold
        sampling_period = sampling_freq/data_rate;
        [OOK_sample, OOK_result] = sample_and_threshold(OOK_filtered, sampling_period, 25/2, encoded_signal_length);
        [BPSK_sample, BPSK_result] = sample_and_threshold(BPSK_output, sampling_period, 0, encoded_signal_length);
        
        % Cyclic Decoding
        decoded_signal_OOK = decode(OOK_result, codeword_length, message_length, 'cyclic/binary', genpoly, trt);
        decoded_signal_BPSK = decode(BPSK_result, codeword_length, message_length, 'cyclic/binary', genpoly, trt);
        
        % Get Bit Errors
        OOK_error = biterr(decoded_signal_OOK, original_signal) ./ signal_length;
        BPSK_error = biterr(decoded_signal_BPSK, original_signal) ./ signal_length;
                
        OOK_average_error = OOK_error + OOK_average_error;
        BPSK_average_error = BPSK_error + BPSK_average_error;
    end
    
    % Plots for SNR @ 5dB
    if (i == 5)
        figure(2)
        subplot(2, 1, 1);
        plot(original_signal, 'b');
        title("Original Signal")
        xlim([0 2000])

        subplot(2, 1, 2);
        plot(encoded_signal, 'b');
        title("Cyclic Encoded Data")
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
        plot(decoded_signal_OOK)
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
        plot(decoded_signal_BPSK)
        title("BPSK Decoded Signal");

        figure(6);
        subplot(3, 1, 1);
        plot(original_signal);
        title("Original Data");
        xlim([0 1024]);
        ylim([0 1]);

        subplot(3, 1, 2);
        plot(decoded_signal_OOK);
        title("OOK Decoded Data");
        xlim([0 1024]);

        subplot(3, 1, 3);
        plot(decoded_signal_BPSK);
        title("BPSK Decoded Data");
        xlim([0 1024]);
    end
      
    OOK_error_rate(i) = OOK_average_error / nb_samples;
    BPSK_error_rate(i) = BPSK_average_error / nb_samples;
end

% Plot OOK vs BPSK bit error rate
SNR_dB = 0:5:50;
figure(1)
p1 = semilogy(SNR_dB, OOK_error_rate,'r-*');
hold on
p2 = semilogy(SNR_dB, BPSK_error_rate, 'b-*');
hold off
ylabel('Bit Error Rate (BER)');
xlabel('SNR (dB)');
legend([p1(1) p2(1)],{'OOK','DPSK'})
xlim([0 50]);