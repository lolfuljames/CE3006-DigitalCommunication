function [generated_noise] = generate_noise(data_length, noise_power)
% This function generates an array of noise with variance equals to noise_power
%   params:
%   data_length - number of noise to be generated (length of the array)
%   noise_power - the required noise power (noise variance) for the generated noise
generated_noise = sqrt(noise_power) .* randn(1, data_length);
end

