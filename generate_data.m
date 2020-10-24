function [generated_data] = generate_data(data_length)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
generated_data = round(rand(1, data_length)) .* 2 - 1;
end