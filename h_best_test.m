function r = h_best_test( in )
  % 15 states -> 3 measurements
  
  r =  [eye(6), zeros(6,9)] * in;

end

