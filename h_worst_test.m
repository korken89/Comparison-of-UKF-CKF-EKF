function r = h_worst_test( in )
  % 15 states -> 6 measurements
  
  r =  [eye(9), zeros(9,6)] * in;

end

