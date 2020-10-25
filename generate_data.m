function [generated_data] = generate_data(data_length)
% This function generates an array of data
%   params:
%   data_length - number of data to be generated (length of the array)
generated_data = round(rand(1, data_length)) .* 2 - 1;
end