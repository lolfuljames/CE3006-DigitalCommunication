function [generated_SNR] = generate_SNR(snr_max, steps)
%   This function generates the list of SNRs from 0:snr_max in step size of steps
generated_SNR = 10.^((0:steps:snr_max)/10);
end