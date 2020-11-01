clear all; close all; clc;
%Assume the number of bits for transmission is 1024
data_length = 1024;
signal_power = 1; 
% Signal Power(dB) = 10 log (Signal_Power/Noise_Power)                 
% 0-50(dB) step size of 5
SNR_dB = 0:1:50;
% SNR = Signal_Power/Noise_Power = 10^(SNR_dB/10)

SNR = convert_dB_to_dec(SNR_dB, 'power');
% SNR = (10.^(SNR_dB/10));

% Noise Power = Signal_Power / SNR
noise_powers = signal_power ./ SNR;
threshold = 0;

% Generate data length of 1024 bits
data = generate_data(data_length); 

test_samples = 20;
%Different SNR value
% Average error for different SNR
bit_errors = [];
for i = 1 : length(SNR)
	bit_errors(i) = 0;
	for j = 1 : test_samples
    % Generate noise according to SNR
		noise = generate_noise(data_length, noise_powers(i));   
		% Received Signal is Data + Noise
		received_signal = data + noise;                      

		% Threshold Logic
		% data >= threshold, treats as 1
		% data < threshold, treats as 0 
		received_signal = 2*(received_signal >= threshold)-1;
		error_signal = received_signal ~= data;
			%Calculate bit error rate during transmission
		bit_errors(i) = bit_errors(i) + mean(error_signal);
	end
  bit_errors(i) = bit_errors(i)/test_samples;
end  

% Calculate theoretical BER
theory_rate = (1 / 2) * erfc(sqrt(SNR / 2));
    
%Graph and Plot the result           
figure(1)
semilogy (SNR_dB, theory_rate,'r', 'linewidth', 1.5);
ylabel('BER');
xlabel('SNR (dB)')
title('BER vs SNR (dB) - Step Size: 1');
hold on
semilogy (SNR_dB, bit_errors,'bx', 'linewidth', 2);
legend('Theoretical BER','Real BER');
axis([0 50 1/(test_samples*data_length) 1]);
hold off

%Graph and Plot the result           
figure(2)
ylabel('BER');
xlabel('SNR (dB)')
semilogy (SNR_dB(1:5:50), theory_rate(1:5:50),'r', 'linewidth', 1.5);
title('BER vs SNR (dB) - Step Size: 5');
hold on
semilogy (SNR_dB(1:5:50), bit_errors(1:5:50),'bx', 'linewidth', 2);
legend('Theoretical BER','Real BER');
axis([0 50 0 1]);
hold off

%data generation
figure(3) 
subplot(311)
plot(data);
xlim([0 1024]);
title('Generated Data (1024 Bits)')
%noise generation
subplot(312) 
plot(noise);
xlim([0 1024]);
title('Generated Noise (1024 Bits)')
%received data generation
subplot(313) 
plot(received_signal);
xlim([0 1024]);
title('Received Data (1024 Bits)')

