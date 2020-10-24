function [sampled, result] = sample_and_threshold(x, sampling_period, threshold, num_bit)
    sampled = zeros(1, num_bit);
    result = sampled;
    for n = 1: num_bit
        sampled(n) = x((2 * n - 1) * sampling_period / 2);
        if(sampled(n) > threshold)
            result(n) = 1;
        else
            result(n) = 0;
        end
    end
end