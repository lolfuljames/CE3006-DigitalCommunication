function [generated_noise] = generate_noise(data_length, signal_power, snr_i)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
noise_power = signal_power ./ snr_i;
generated_noise = sqrt(noise_power / 2) .* randn(1, data_length);
end

