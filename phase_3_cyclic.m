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
% Low-pass filter & High-pass filter - 6th order, 0.2 cutoff frequency
[b, a] = butter(6, 0.2);
[d, c] = butter(6, 0.2, 'high');
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
SNR = generate_SNR(MAX_dB, 5); % SNR length = 20

OOK_error_rate = zeros(length(SNR));
BPSK_error_rate = zeros(length(SNR));

% For each value of SNR
for i = 1 : length(SNR) 
    
	OOK_average_error = 0;
	BPSK_average_error = 0;
    result = zeros(1, nb_samples);
    
    for j = 1 : nb_samples
        % Original signal 
        original_signal = round(rand(1,signal_length));
        % Cyclic encoding
        encoded_signal = encode(original_signal, codeword_length, message_length, 'cyclic/binary', genpoly);
        
        % Sampling
        sampled_signal = zeros(1, sampled_signal_length);
        for k = 1: signal_length - 1
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
        
        % Generate White Gaussian Channel Noise for OOK and BPSK
        OOK_noise_power = OOK_signal_power ./ SNR(i);
        OOK_noise = sqrt(OOK_noise_power/2) .*randn(1,signal_length);
        
        BPSK_noise_power = BPSK_signal_power ./ SNR(i);
        BPSK_noise = sqrt(BPSK_noise_power/2) .* randn(1, signal_length);
        
        % Received Signal with added Channel Noise for OOK and BPSK
        OOK_received = OOK_signal + OOK_noise;
        BPSK_received = BPSK_signal + BPSK_noise;
        
        % Non-Coherent Detection: OOK (Lecture notes 04 - pg 18)
        % Squared
        OOK_squared = OOK_received .* OOK_received;
        % Low-Pass Filter
        OOK_filtered = filtfilt(b, a, OOK_squared);
  
        % Non-Coherent Detection: BPSK (Lecture notes 04 - pg 50)
        % Squared
        BPSK_squared = BPSK_received .* BPSK_received;
        % Bandpass Filter
        BPSK_filtered = filtfilt(d, c, BPSK_squared);
        % Frequency Divider
        BPSK_divided = interp(BPSK_filtered, 2);
        BPSK_divided = BPSK_divided(1:length(BPSK_filtered));
        % Multiply with recovered carrier
        BPSK_multiplied = BPSK_divided .* BPSK_received;
        % Low-Pass Filter
        BPSK_output = filtfilt(b, a, BPSK_multiplied);
        
        % Demodulation: Sampling and Threshold
        sampling_period = sampling_freq/data_rate;
        [OOK_sample, OOK_result] = sample_and_threshold(OOK_filtered, sampling_period, amp/2, codeword_length);
        [BPSK_sample, BPSK_result] = sample_and_threshold(BPSK_output, sampling_period, 0, encoded_length);
        
        % Cyclic Decoding
        decoded_signal_OOK = decode(OOK_result, codeword_length, message_length, 'cyclic/binary', genpoly, trt);
        decoded_signal_BPSK = decode(BPSK_result, codeword_length, message_length, 'cyclic/binary', genpoly, trt);
        
        % Get Bit Errors
        OOK_error_ber = biterr(decoded_signal_OOK, orignal_signal);
        BPSK_error_ber = biterr(decoded_signal_BPSK, original_signal);
        
        OOK_error = 0;
        BPSK_error = 0;
        for k = 1: signal_length - 1
        	if(decoded_signal_OOK(k) ~= original_signal(k))
               OOK_error = OOK_error + 1;
         end
         if(decoded_signal_BPSK(k) ~= original_signal(k))
               BPSK_error = BPSK_error + 1;
         end
        end
		
        disp(OOK_error_ber)
        disp(BPSK_error_ber)
        disp(OOK_error)
        disp(BPSK_error)
		
        OOK_average_error = OOK_error + OOK_average_error;
        BPSK_average_error = BPSK_error + BPSK_average_error;
    end
    
    OOK_error_rate(i) = OOK_average_error / test_samples;
    BPSK_error_rate(i) = BPSK_average_error / test_samples;
end
