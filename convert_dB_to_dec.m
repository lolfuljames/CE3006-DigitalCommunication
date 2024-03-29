function [output_dec] = convert_dB_to_dec(input_dB, type)
%convert_dB_to_dec Convert dB to decimal
%   input_dB - input to be converted
%   type - default to 'non_power'

if nargin == 1
  type = 'non_power'; % default type to be 'non_power'
end

if type == "power"
    output_dec = 10 .^ (input_dB ./ 10);
else
    output_dec = 20 .^ (input_dB ./ 10);
end

