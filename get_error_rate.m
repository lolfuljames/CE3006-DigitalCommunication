function error_rate = get_error_rate(actual, target)
    error = 0;
    for i=1:length(actual)
        if(actual(i) ~= target(i))
            error = error + 1;
        end
    end
    error_rate = error ./ length(actual);
end