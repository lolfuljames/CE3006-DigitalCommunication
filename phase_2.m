%
% Data Assumptions
% 
% Amplitude of OOK is set to 5
% Sampling Frequency = 16 x Carrier Frequency
% 
% Generate data -> OOK/BPSK Modulation -> Noise -> Demodulation -> Comparison
% 

clear all; close all; clc;
carrier_freq = 10000; %10kHz
sample_freq = 16 * carrier_freq;
data_rate = 1000; %1kbps
data_length = 1024;
amp = 5;

%For FSK modulation
fsk_freq_1 = 30000;
fsk_freq_2 = 10000;

% Low Pass 6th order Butterworth filter with 0.2 normalised cutoff freq
[b, a] = butter(6, 0.2);

% Time simulation
t = 0: 1/sample_freq : data_length/data_rate;

% Carrier Signal Generation
carrier_signal = amp .* cos(2*pi*carrier_freq*t);
fsk_carrier_signal_1 = amp .* cos(2*pi*fsk_freq_1*t);
fsk_carrier_signal_2 = amp .* cos(2*pi*fsk_freq_2*t);

% Length of transmitted signal
signal_length = sample_freq*data_length/data_rate + 1;

% SNR values to test             
SNR_dB = 0:1:20;
SNR = (10.^(SNR_dB/10));

% Number of tests per SNR
test_samples = 100;

OOK_error_rate = zeros([length(SNR) 1]);
BPSK_error_rate = zeros([length(SNR) 1]);
BFSK_error_rate = zeros([length(SNR) 1]);

%*********************Baseband***********************%
% Generate symbols(data) with NRZ-L
data = round(rand(1,data_length));
% Sampled signal generated from data
signal = zeros(1, signal_length);
for k = 1: signal_length - 1
    signal(k) = data(ceil(k*data_rate/sample_freq));
end
signal(signal_length) = signal(signal_length - 1);

%**************Baseband -> Bandpass ****************

% OOK Modulation
OOK_signal = carrier_signal .* signal;

% BPSK Modulation
BPSK_source_signal = signal .* 2 - 1;
BPSK_signal = carrier_signal .* BPSK_source_signal;

% BFSK Modulation
BFSK_source_signal_1 = fsk_carrier_signal_1 .* (signal == 1);
BFSK_source_signal_0 = fsk_carrier_signal_2 .* (signal == 0);
BFSK_signal = BFSK_source_signal_1 + BFSK_source_signal_0;


%***************Xmit through Channel ***************
%Simulates channel effect by adding in noise
OOK_signal_power = (norm(OOK_signal)^2)/signal_length;
BPSK_signal_power = (norm(BPSK_signal)^2)/signal_length;
BFSK_signal_power = (norm(BFSK_signal)^2)/signal_length;

noise_powers_OOK = OOK_signal_power ./ SNR;
noise_powers_BPSK = BPSK_signal_power ./ SNR;
noise_powers_BFSK = BFSK_signal_power ./ SNR;

% For different SNR values, test over 20 samples
for i = 1 : length(SNR)
	OOK_average_error = 0;
    BPSK_average_error = 0;
    BFSK_average_error = 0;
    
	for j = 1 : test_samples
        
        %Generate Noise
		noise_power_OOK = OOK_signal_power ./SNR(i);
		noise_OOK = sqrt(noise_power_OOK) .*randn(1,signal_length);
        
        noise_power_BPSK = BPSK_signal_power ./SNR(i);
        noise_BPSK = sqrt(noise_power_BPSK) .*randn(1,signal_length);
        
        noise_power_BFSK = BFSK_signal_power ./SNR(i);
        noise_BFSK = sqrt(noise_power_BFSK) .*randn(1,signal_length);
        
		%Received Signal OOK
		OOK_received = OOK_signal+noise_OOK;
        
		%Received Signal BPSK
		BPSK_received = BPSK_signal+noise_BPSK;
        
        %Received Signal BFSK
        BFSK_received = BFSK_signal+noise_BFSK;
        
        %*****************Receiver detection***************
        %OOK detection
        OOK_squared = OOK_received .* (2 .* carrier_signal);
        OOK_filtered = filtfilt(b, a, OOK_squared);
        
        %BPSK detection
        BPSK_squared = BPSK_received .* (2 .* carrier_signal);
        BPSK_filtered = filtfilt(b, a, BPSK_squared);
        
        %BFSK detection
        BFSK_carrier_1_corr = BFSK_received .* (2 .* fsk_carrier_signal_1);
        BFSK_branch_1_filtered = filtfilt(b, a, BFSK_carrier_1_corr);
        BFSK_carrier_2_corr = BFSK_received .* (2 .* fsk_carrier_signal_2);
        BFSK_branch_2_filtered = filtfilt(b, a, BFSK_carrier_2_corr);
        BFSK_differenced = BFSK_branch_1_filtered - BFSK_branch_2_filtered;
        
        %sampling AND threshold
        sample_period = sample_freq / data_rate;
        [OOK_sample, OOK_result] = sample_and_threshold(OOK_filtered, sample_period, amp^2/2, data_length);
        [BPSK_sample, BPSK_result] = sample_and_threshold(BPSK_filtered, sample_period, 0,data_length);
        [BFSK_sample, BFSK_result] = sample_and_threshold(BFSK_differenced, sample_period, 0, data_length);
        
		%Calculate the average error for every runtime
		%Avg_ErrorOOK = num_error(resultOOK, EncodeHamming, Num_Bit) + Avg_ErrorOOK;                   
        %Avg_ErrorBPSK = num_error(resultBPSK, EncodeHamming, Num_Bit) + Avg_ErrorBPSK;
        
        OOK_error = 0;
        BPSK_error = 0;
        BFSK_error = 0;
        for k = 1: data_length
            if(OOK_result(k) ~= data(k))
                OOK_error = OOK_error + 1;
            end
            if(BPSK_result(k) ~= data(k))
                BPSK_error = BPSK_error + 1;
            end
            if(BFSK_result(k) ~= data(k))
                BFSK_error = BFSK_error + 1;
            end
        end
        OOK_error = OOK_error./data_length;
        BPSK_error = BPSK_error./data_length;
        BFSK_error = BFSK_error./data_length;
        OOK_average_error = OOK_error + OOK_average_error;
        BPSK_average_error = BPSK_error + BPSK_average_error;
        BFSK_average_error = BFSK_error + BFSK_average_error;

    end
    
%    Plot the 5db SNR signals
    if (SNR_dB(i) == 5)
        %Plot of original data with respect to time
        figure(2)
        plot(data, 'b');
        title("Original Data")
        xlim([0 1100])
        
        % Plotting of spectrogram for modulated signal and corrupted
        % signal for OOK
        figure(3)
        subplot(2, 1, 1);
        spectrogram(OOK_signal, 'yaxis');
        title("Transmitted OOK Modulated Signal")

        subplot(2, 1, 2);
        spectrogram(OOK_received, 'yaxis');
        title("Received OOK Modulated Signal")
        
        % Plotting of spectrogram for modulated signal and corrupted
        % signal for BPSK
        figure(4)
        subplot(2, 1, 1);
        spectrogram(BPSK_signal, 'yaxis');
        title("Transmitted BPSK Modulated Signal")

        subplot(2, 1, 2);
        spectrogram(BPSK_received, 'yaxis');
        title("Received BPSK Modulated Signal")
        
        % Plotting of spectrogram for modulated signal and corrupted
        % signal for BFSK
        figure(5)
        subplot(2, 1, 1);
        spectrogram(BFSK_signal, 'yaxis');
        title("Transmitted BFSK Modulated Signal")

        subplot(2, 1, 2);
        spectrogram(BFSK_received, 'yaxis');
        title("Received BFSK Modulated Signal")
        
        %Plotting the baseband data and modulated data before transmission
        figure(7);
        subplot(4,1,1);
        plot(signal(1:1000));
        title("Baseband oversampled signal (snippet)");
        
        subplot(4, 1, 2);
        plot(OOK_signal(1:1000));
        title("OOK modulated signal (Snippet)");
        
        subplot(4, 1, 3);
        plot(BPSK_signal(1:1000));
        title("BPSK modulated signal (Snippet)");
        
        subplot(4, 1, 4);
        plot(BFSK_signal(1:1000));
        title("BFSK modulated signal (Snippet)");
        
        %Plotting of Received signal (corrupted with noise)
        figure(8);
        subplot(3, 1, 1);
        plot(OOK_received(1:1000));
        title("OOK received with noise (Snippet)");
        
        subplot(3,1,2);
        plot(BPSK_received(1:1000));
        title("BPSK received with noise (Snippet)");
        
        subplot(3,1,3);
        plot(BFSK_received(1:1000));
        title("BFSK received with noise (Snippet)")
        
        %Plotting of demodulated signal (mixed and passed through low pass
        %filter
        figure(9);
        subplot(3, 1, 1);
        plot(OOK_sample);
        title("OOK demodulated and sampled");
        
        subplot(3,1,2);
        plot(BPSK_sample);
        title("BPSK demodulated and sampled");
        
        subplot(3,1,3);
        plot(BFSK_sample);
        title("BFSK demodulated and sampled")
        
        %Plotting of results from threshold logic (decoded)
        figure(10);
        subplot(4, 1, 1);
        plot(data);
        title("Original Data");
        xlim([0 1024]);
        ylim([0 1]);

        subplot(4, 1, 2);
        plot(OOK_result);
        title("OOK Decoded Data");
        xlim([0 1024]);

        subplot(4, 1, 3);
        plot(BPSK_result);
        title("BPSK Decoded Data");
        xlim([0 1024]);
        
        subplot(4, 1, 4);
        plot(BFSK_result);
        title("BFSK Decoded Data");
        xlim([0 1024]);
        
        
    end
	OOK_error_rate(i) = OOK_average_error / test_samples;
    BPSK_error_rate(i) = BPSK_average_error / test_samples;
    BFSK_error_rate(i) = BFSK_average_error / test_samples;
end

% Calculate OOK theoretical BER
OOK_E1 = (1 / 2) * amp^2 / data_rate;
OOK_E0 = 0;
OOK_Eb = (1 / 2) * (OOK_E1 + OOK_E0);
OOK_No = noise_powers_OOK ./ data_rate ./ 2;
OOK_theory_rate = (1 / 2) .* erfc(sqrt(OOK_Eb ./ (2 .* OOK_No)));

% Calculate BFSK theoretical BER
BFSK_E1 = (1 / 2) * amp^2 / data_rate;
BFSK_E0 = (1 / 2) * amp^2 / data_rate;
BFSK_Eb = (1 / 2) * (BFSK_E1 + BFSK_E0);
BFSK_No = noise_powers_BFSK ./ data_rate ./ 2;
BFSK_theory_rate = (1 / 2) .* erfc(sqrt(BFSK_Eb ./ (2 .* BFSK_No)));

% Calculate BPSK theoretical BER
BPSK_E1 = (1 / 2) * amp^2 / data_rate;
BPSK_E0 = (1 / 2) * amp^2 / data_rate;
BPSK_Eb = (1 / 2) * (BPSK_E1 + BPSK_E0);
BPSK_No = noise_powers_BPSK ./ data_rate ./ 2;
BPSK_theory_rate = (1 / 2) .* erfc(sqrt(BPSK_Eb ./ BPSK_No));

% Plot OOK vs DBSK bit error rate
% Plot OOK vs BFSK and BPSK bit error rate
figure(1)
semilogy (SNR_dB, OOK_theory_rate,'r', 'linewidth', 1.5);
hold on
semilogy (SNR_dB, BPSK_theory_rate,'b', 'linewidth', 1.5);
hold on
semilogy (SNR_dB, BFSK_theory_rate,'g', 'linewidth', 1.5);
hold on
plot1 = semilogy(SNR_dB, OOK_error_rate,'r*');
hold on
plot2 = semilogy(SNR_dB, BPSK_error_rate, 'b*');
plot3 = semilogy(SNR_dB, BFSK_error_rate, 'g*');
hold off
ylabel('Bit Error Rate (BER)');
xlabel('SNR (dB)');
legend([plot1(1) plot2(1) plot3(1)],{'OOK','BPSK','BFSK'})
xlim([0 50]);