function [generated_noise] = generate_noise(data_length, noise_power)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
generated_noise = sqrt(noise_power / 2) .* randn(1, data_length);
end

