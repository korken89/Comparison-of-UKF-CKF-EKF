function r = h_best_test( in )
  % 15 states -> 3 measurements
  
  r =  [eye(3), zeros(3,12)] * in;

end

