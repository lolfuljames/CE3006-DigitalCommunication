%
% Data Assumptions
% 
% Amplitude of OOK is set to 5
% Hamming Code: Codewode 7 bits, Data 4 bits
% Sampling Frequency = 16 x Carrier Frequency
%
clear all; close all; clc;
carrier_freq = 10000; %10kHz
sample_freq = 16 * carrier_freq;
data_rate = 1000; %1kbps
data_length = 1024;
codeword_length = 1792;
amp = 5;

%low pass butterworth filter
%6th order, 0.2 cutoff frequency
[b, a] = butter(6, 0.2);

%high pass butterworth filter
[d, c] = butter(6, 0.2, 'high');

%time
t = 0: 1/sample_freq : codeword_length/data_rate;

%Carrier
carrier_signal = amp .* cos(2*pi*carrier_freq*t);

%signal length
signal_length = sample_freq*codeword_length/data_rate + 1;

%SNR_dB = 10 log (Signal_Power/Noise_Power)                 
SNR_dB = 0:1:50;
%==> SNR = Signal_Power/Noise_Power = 10^(SNR_dB/10)
SNR = (10.^(SNR_dB/10));

% Set run times
test_samples = 20;

OOK_error_rate = zeros([length(SNR) 1]);
BPSK_error_rate = zeros([length(SNR) 1]);

%Different SNR value
for i = 1 : length(SNR)
	OOK_average_error = 0;
    BPSK_average_error = 0;
    
	for j = 1 : test_samples
        %Generate Data
        data = round(rand(1,data_length));
        %encode
        hamming_signal= encode(data,7,4,'hamming/fmt');

        signal = zeros(1, signal_length);
        for k = 1: signal_length - 1
            signal(k) = hamming_signal(ceil(k*data_rate/sample_freq));
        end
        signal(signal_length) = signal(signal_length - 1);

        %on-off keying
        OOK_signal = carrier_signal .* signal;

        %binary phase shift keying
        BPSK_source_signal = signal .* 2 - 1;
        BPSK_signal = carrier_signal .* BPSK_source_signal;

        OOK_signal_power = (norm(OOK_signal)^2)/signal_length;
        BPSK_signal_power = (norm(BPSK_signal)^2)/signal_length;
        
        %Generate Noise OOK
		OOK_noise_power = OOK_signal_power ./SNR(i);
		OOK_noise = sqrt(OOK_noise_power/2) .*randn(1,signal_length);
		%Received Signal OOK
		OOK_received = OOK_signal+OOK_noise;
        
        %OOK detection
        OOK_squared = OOK_received .* OOK_received;
        %low pass filter
        OOK_filtered = filtfilt(b, a, OOK_squared);
        
        %Generate Noise BPSK
		BPSK_noise_power = BPSK_signal_power ./SNR(i);
		BPSK_noise = sqrt(BPSK_noise_power/2) .*randn(1,signal_length);
		%Received Signal BPSK
		BPSK_received = BPSK_signal+BPSK_noise;
        
        %non-coherent detection
        BPSK_squared = BPSK_received .* BPSK_received;
        %high pass filter (supposingly band pass filter)
        BPSK_filtered = filtfilt(d, c, BPSK_squared);
        
        %frequency divider
        BPSK_divided = interp(BPSK_filtered, 2);
        BPSK_divided = BPSK_divided(1:length(BPSK_filtered));
        
        %Multiple and Low Pass Filter
        BPSK_multiplied = BPSK_divided .* BPSK_received;
        BPSK_output = filtfilt(b, a, BPSK_multiplied);
        
        %demodulate
        %sampling AND threshold
        sample_period = sample_freq / data_rate;
        [OOK_sample, OOK_result] = sample_and_threshold(OOK_filtered, sample_period, amp/2, codeword_length);
        [BPSK_sample, BPSK_result] = sample_and_threshold(BPSK_output, sample_period, 0, codeword_length);
        
		%Calculate the average error for every runtime
		%Avg_ErrorOOK = num_error(resultOOK, EncodeHamming, Num_Bit) + Avg_ErrorOOK;                   
        %Avg_ErrorBPSK = num_error(resultBPSK, EncodeHamming, Num_Bit) + Avg_ErrorBPSK;
        
        OOK_decoded = decode(OOK_result,7,4,'hamming/fmt');
        BPSK_decoded = decode(BPSK_result,7,4,'hamming/fmt');
        
        OOK_error = 0;
        BPSK_error = 0;
        for k = 1: data_length - 1
            if(OOK_decoded(k) ~= data(k))
                OOK_error = OOK_error + 1;
            end
            if(BPSK_decoded(k) ~= data(k))
                BPSK_error = BPSK_error + 1;
            end
        end
        OOK_error = OOK_error./codeword_length;
        BPSK_error = BPSK_error./codeword_length;
        OOK_average_error = OOK_error + OOK_average_error;
        BPSK_average_error = BPSK_error + BPSK_average_error;

    end
%     plot the 5db SNR signals
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
    
	OOK_error_rate(i) = OOK_average_error / test_samples;
    BPSK_error_rate(i) = BPSK_average_error / test_samples;
end

figure(1)
p1 = semilogy(SNR_dB, OOK_error_rate,'r-*');
hold on
p2 = semilogy(SNR_dB, BPSK_error_rate, 'b-*');
%axis([0 20 10^(-5) 1]);
hold off
ylabel('Bit Error Rate (BER)');
xlabel('SNR (dB)');
legend([p1(1) p2(1)],{'OOK','DPSK'})

