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
amp = 5;

%For FSK modulation
fsk_freq_1 = 30000;
fsk_freq_2 = 10000;

%Carrier signal generation
carrier_signal = amp .* cos(2*pi*carrier_freq*t);
fsk_carrier_signal_1 = amp .* cos(2*pi*fsk_freq_1*t);
fsk_carrier_signal_2 = amp .* cos(2*pi*fsk_freq_2*t);

% Sampled signal length
sampled_signal_length = sampling_freq*encoded_signal_length/data_rate + 1;
sampled_unencoded_signal_length = sampling_freq*signal_length/data_rate + 1;

% Number of samples
nb_samples = 20;

% Low-pass filter - 6th order, 0.2 cutoff frequency
[b, a] = butter(6, 0.2);

% maximum SNR
MAX_dB = 20;

% Generator and Syndrome Table for Cyclic Code (Channel Encoding-Decoding)
% Generator polynomial & Parity-check matrix for cyclic encoding
genpoly = cyclpoly(codeword_length, message_length);
parmat = cyclgen(codeword_length, genpoly);
% Syndrome table for cyclic decoding
trt = syndtable(parmat);


% SNR values to test
SNR_dB = 0:1:10;
SNR = (10.^(SNR_dB/10));

OOK_error_rate = zeros([length(SNR) 1]);
unencoded_OOK_error_rate = zeros([length(SNR) 1]);
BPSK_error_rate = zeros([length(SNR) 1]);
BFSK_error_rate = zeros([length(SNR) 1]);

% Original signal 
original_signal = round(rand(1,signal_length));
% Cyclic encoding
encoded_signal = encode(original_signal, codeword_length, message_length, 'cyclic/binary', genpoly);

% Sampling
sampled_signal = zeros(1, sampled_signal_length);
sampled_unencoded_signal = zeros(1, sampled_unencoded_signal_length);
for k = 1: sampled_signal_length - 1
    sampled_signal(k) = encoded_signal(ceil(k*data_rate/sampling_freq));
end
for k = 1: sampled_unencoded_signal_length - 1
    sampled_unencoded_signal(k) = original_signal(ceil(k*data_rate/sampling_freq));
end
sampled_signal(sampled_signal_length) = sampled_signal(sampled_signal_length - 1);
sampled_unencoded_signal(sampled_unencoded_signal_length) = sampled_unencoded_signal(sampled_unencoded_signal_length - 1);

% Modulation: On-Off Keying
OOK_signal = carrier_signal .* sampled_signal;
unencoded_OOK_signal = carrier_signal(1:sampled_unencoded_signal_length) .* sampled_unencoded_signal;

% Modulation: Binary Phase Shift Keying 
BPSK_source_signal = sampled_signal .* 2 - 1; % put to -1 +1
BPSK_signal = carrier_signal .* BPSK_source_signal;

%Modulation: Binary FSK
BFSK_source_signal_1 = fsk_carrier_signal_1 .* (sampled_signal == 1);
BFSK_source_signal_0 = fsk_carrier_signal_2 .* (sampled_signal == 0);
BFSK_signal = BFSK_source_signal_1 + BFSK_source_signal_0;

OOK_signal_power = (norm(OOK_signal)^2)/sampled_signal_length;
unencoded_OOK_signal_power = (norm(unencoded_OOK_signal)^2)/sampled_unencoded_signal_length;
BPSK_signal_power = (norm(BPSK_signal)^2)/sampled_signal_length;
BFSK_signal_power = (norm(BFSK_signal)^2)/sampled_signal_length;

% For each value of SNR, test of 20 samples.
for i = 1 : length(SNR) 
    
	OOK_average_error = 0;
    unencoded_OOK_average_error = 0;
	BPSK_average_error = 0;
    BFSK_average_error = 0;
    result = zeros(1, nb_samples);
    
    for j = 1 : nb_samples
        
        % Generate White Gaussian Channel Noise for OOK and BPSK
        noise_power_OOK = OOK_signal_power ./ SNR(i);
        noise_OOK = sqrt(noise_power_OOK) .*randn(1,sampled_signal_length);
        
        noise_power_unencoded_OOK = unencoded_OOK_signal_power ./ SNR(i);
        noise_unencoded_OOK =  sqrt(noise_power_unencoded_OOK) .*randn(1,sampled_unencoded_signal_length);

        noise_power_BPSK = BPSK_signal_power ./ SNR(i);
        noise_BPSK = sqrt(noise_power_BPSK) .*randn(1,sampled_signal_length);
        
        noise_power_BFSK = BFSK_signal_power ./ SNR(i);
        noise_BFSK = sqrt(noise_power_BFSK) .*randn(1,sampled_signal_length);
        
        % Received Signal with added Channel Noise for OOK and BPSK
        OOK_received = OOK_signal + noise_OOK;
        unencoded_OOK_received = unencoded_OOK_signal + noise_unencoded_OOK;
        BPSK_received = BPSK_signal + noise_BPSK;
        BFSK_received = BFSK_signal + noise_BFSK;
        
        % Non-Coherent Detection: OOK (Lecture notes 04 - pg 18)
        % Squared
        OOK_squared = OOK_received .* 2 .* carrier_signal;
        unencoded_OOK_squared = unencoded_OOK_received .* 2 .* carrier_signal(1:sampled_unencoded_signal_length);

        % Low-Pass Filter
        OOK_filtered = filtfilt(b, a, OOK_squared);
        unencoded_OOK_filtered = filtfilt(b, a, unencoded_OOK_squared);
  
        % Non-Coherent Detection: BPSK (Lecture notes 04 - pg 50)
        % Squared
        BPSK_squared = BPSK_received .* 2 .* carrier_signal;
        % Low-Pass Filter
        BPSK_output = filtfilt(b, a, BPSK_squared);
        
        %BFSK detection
        BFSK_carrier_1_corr = BFSK_received .* (2 .* fsk_carrier_signal_1);
        BFSK_branch_1_filtered = filtfilt(b, a, BFSK_carrier_1_corr);
        BFSK_carrier_2_corr = BFSK_received .* (2 .* fsk_carrier_signal_2);
        BFSK_branch_2_filtered = filtfilt(b, a, BFSK_carrier_2_corr);
        BFSK_differenced = BFSK_branch_1_filtered - BFSK_branch_2_filtered;
        
        % Demodulation: Sampling and Threshold
        sampling_period = sampling_freq/data_rate;
        [OOK_sample, OOK_result] = sample_and_threshold(OOK_filtered, sampling_period, 25/2, encoded_signal_length);
        [unencoded_OOK_sample, unencoded_OOK_result] = sample_and_threshold(unencoded_OOK_filtered, sampling_period, 25/2, signal_length);
        [BPSK_sample, BPSK_result] = sample_and_threshold(BPSK_output, sampling_period, 0, encoded_signal_length);
        [BFSK_sample, BFSK_result] = sample_and_threshold(BFSK_differenced, sampling_period, 0, encoded_signal_length);

        % Cyclic Decoding
        decoded_signal_OOK = decode(OOK_result, codeword_length, message_length, 'cyclic/binary', genpoly, trt);
        decoded_signal_BPSK = decode(BPSK_result, codeword_length, message_length, 'cyclic/binary', genpoly, trt);
        decoded_signal_BFSK = decode(BFSK_result, codeword_length, message_length, 'cyclic/binary', genpoly, trt);
        
        % Get Bit Errors
        OOK_error = biterr(decoded_signal_OOK, original_signal) ./ signal_length;
        unencoded_OOK_error = biterr(unencoded_OOK_result, original_signal) ./ signal_length;
        BPSK_error = biterr(decoded_signal_BPSK, original_signal) ./ signal_length;
        BFSK_error = biterr(decoded_signal_BFSK, original_signal) ./signal_length;
                
        OOK_average_error = OOK_error + OOK_average_error;
        unencoded_OOK_average_error = unencoded_OOK_error + unencoded_OOK_average_error;
        BPSK_average_error = BPSK_error + BPSK_average_error;
        BFSK_average_error = BFSK_error + BFSK_average_error;
    end
         
    OOK_error_rate(i) = OOK_average_error / nb_samples;
    BPSK_error_rate(i) = BPSK_average_error / nb_samples;
    BFSK_error_rate(i) = BFSK_average_error / nb_samples;
    unencoded_OOK_error_rate(i) = unencoded_OOK_average_error / nb_samples; 
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

        figure(3)
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

        figure(3)
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
        
        figure(5)
        subplot(4, 1, 1);
        spectrogram(BFSK_signal,'yaxis')
        title("Transmitted BFSK Modulated Signal")

        subplot(4, 1, 2);
        spectrogram(BFSK_received,'yaxis')
        title("Received BFSK Modulated Signal")

        subplot(4, 1, 3);
        plot(BFSK_sample)
        title("BFSK Demodulated Signal")

        subplot(4, 1, 4);
        plot(decoded_signal_BFSK)
        title("BFSK Decoded Signal");

        figure(6);
        subplot(4, 1, 1);
        plot(original_signal);
        title("Original Data");
        xlim([0 1024]);
        ylim([0 1]);

        subplot(4, 1, 2);
        plot(decoded_signal_OOK);
        title("OOK Decoded Data");
        xlim([0 1024]);

        subplot(4, 1, 3);
        plot(decoded_signal_BPSK);
        title("BPSK Decoded Data");
        xlim([0 1024]);
        
        subplot(4, 1, 4);
        plot(decoded_signal_BFSK);
        title("BFSK Decoded Data");
        xlim([0 1024]);
    end

end

% Plot OOK vs BPSK bit error rate
figure(1)
p1 = semilogy(SNR_dB, OOK_error_rate,'r-*');
hold on
p2 = semilogy(SNR_dB, BPSK_error_rate, 'b-*');
p3 = semilogy(SNR_dB, BFSK_error_rate, 'g-*');
p4 = semilogy(SNR_dB, unencoded_OOK_error_rate, 'k-*');
hold off
ylabel('Bit Error Rate (BER)');
xlabel('SNR (dB)');
legend([p1(1) p2(1) p3(1) p4(1)],{'Cyclic/OOK','Cyclic/BPSK',' Cyclic/BFSK','Unencoded/OOK'})
xlim([0 50]);