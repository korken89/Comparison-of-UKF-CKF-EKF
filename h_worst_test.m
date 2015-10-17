function r = h_worst_test( in )
  % 15 states -> 6 measurements
  
  r =  [eye(6), zeros(6,9)] * in;

end

