function r = f_test(in, u)
    % 15 states -> 15 states
    r = in + ones(15,1) * u; 
end