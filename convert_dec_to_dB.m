function [output_dB] = convert_dec_to_dB(input_dec, type)
%convert_dB_to_dec Convert decimal to dB
%   input_dec - input to be converted
%   type - default to 'non_power'

if nargin == 1
  type = 'non_power'; % default type to be 'non_power'
end

if type == "power"
    output_dB = 10 * log10(input_dec);
else
    output_dB = 20 * log10(input_dec);
end

