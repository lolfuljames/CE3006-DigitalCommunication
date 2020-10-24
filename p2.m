%
% Data Assumptions
% 
% Amplitude of OOK is set to 5
% Sampling Frequency = 16 x Carrier Frequency
%

clear all; close all; clc;
carrier_freq = 10000; %10kHz
sample_freq = 16 * carrier_freq;
data_rate = 1000; %1kbps
data_length = 1024;
amp = 5;

%low pass butterworth filter
%6th order, 0.2 cutoff frequency
[b, a] = butter(6, 0.2);

%high pass butterworth filter
[d, c] = butter(6, 0.2, 'high');

%time
t = 0: 1/sample_freq : data_length/data_rate;

%Carrier
carrier_signal = amp .* cos(2*pi*carrier_freq*t);

%signal length
signal_length = sample_freq*data_length/data_rate + 1;

%SNR_dB = 10 log (Signal_Power/Noise_Power)                 
SNR_dB = 0:1:20;
%==> SNR = Signal_Power/Noise_Power = 10^(SNR_dB/10)
SNR = (10.^(SNR_dB/10));

% Set run times
test_samples = 20;

OOK_error_rate = zeros(length(SNR));
BPSK_error_rate = zeros(length(SNR));

%Different SNR value
for i = 1 : length(SNR)
	OOK_average_error = 0;
    BPSK_average_error = 0;
    
	for j = 1 : test_samples
        %Generate Data
        data = round(rand(1,data_length));

        signal = zeros(1, signal_length);
        for k = 1: signal_length - 1
            signal(k) = data(ceil(k*data_rate/sample_freq));
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
        [OOK_sample, OOK_result] = sample_and_threshold(OOK_filtered, sample_period, amp/2, data_length);
        [BPSK_sample, BPSK_result] = sample_and_threshold(BPSK_output, sample_period, 0, data_length);
        
		%Calculate the average error for every runtime
		%Avg_ErrorOOK = num_error(resultOOK, EncodeHamming, Num_Bit) + Avg_ErrorOOK;                   
        %Avg_ErrorBPSK = num_error(resultBPSK, EncodeHamming, Num_Bit) + Avg_ErrorBPSK;
        
        OOK_error = 0;
        BPSK_error = 0;
        for k = 1: data_length - 1
            if(OOK_result(k) ~= data(k))
                OOK_error = OOK_error + 1;
            end
            if(BPSK_result(k) ~= data(k))
                BPSK_error = BPSK_error + 1;
            end
        end
        OOK_average_error = OOK_error + OOK_average_error;
        BPSK_average_error = BPSK_error + BPSK_average_error;

	end
	OOK_error_rate(i) = OOK_average_error / test_samples;
    BPSK_error_rate(i) = BPSK_average_error / test_samples;
end

figure(1)
plot1 = semilogy(SNR_dB, OOK_error_rate,'r-*');
hold on
plot2 = semilogy(SNR_dB, BPSK_error_rate, 'b-*');
%axis([0 20 10^(-5) 1]);
hold off
ylabel('Error Rate');
xlabel('SNR (dB)');
legend([plot1(1) plot2(1)],{'OOK','DPSK'})

figure(2)
plot(data, 'b');
title("Original Data")
xlim([0 1100])

figure(3)
subplot(4, 1, 1);
plot(OOK_signal, 'b');
title("Transmitted OOK Modulated Signal")

subplot(4, 1, 2);
plot(OOK_received, 'b');
title("Received OOK Modulated Signal")

subplot(4, 1, 3);
plot(OOK_sample)
title("OOK Demodulated Signal")

subplot(4, 1, 4);
plot(OOK_result);
title("OOK Decoded Signal");

figure(4)
subplot(4, 1, 1);
plot(BPSK_signal, 'b');
title("Transmitted BPSK Modulated Signal")

subplot(4, 1, 2);
plot(BPSK_received, 'b');
title("Received BPSK Modulated Signal")

subplot(4, 1, 3);
plot(BPSK_sample)
title("BPSK Demodulated Signal")

subplot(4, 1, 4);
plot(BPSK_result);
title("BPSK Decoded Signal");

figure(5);
subplot(3, 1, 1);
plot(data);
title("Original Data");
xlim([0 1024]);
ylim([0 1]);

subplot(3, 1, 2);
plot(OOK_result);
title("OOK Decoded Data");
xlim([0 1024]);

subplot(3, 1, 3);
plot(BPSK_result);
title("BPSK Decoded Data");
xlim([0 1024]);