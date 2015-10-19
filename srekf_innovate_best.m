function [ x_hat, Sp ] = srekf_innovate_best( x_in, u_in, z_in, Sp_in, Sq_in, Sr_in, f_function, h_function)  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % SR-EKF starts here!
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Prediction Update
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Model prediction (f-function)
  x_hat = f_function(x_in, u_in);

  % Calculate A matrix
  A = eye(length(x_in));
  
  % State covariance matrix update based on model
  [~, R] = qr([A * Sp_in, Sq_in]', 0);
  Sp = R';
  

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Measurement Update
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      
  % Measurement prediction function
  h = h_function(x_hat);  
  
  % Calculate error
  y = z_in - h;    
  
  % The H matrix maps the measurement to the states
  H = [eye(6), zeros(6,9)];

  % Measurement covariance update
  [~, R] = qr([H * Sp, Sr_in]', 0);
  Ss = R';
     
  % Calculate Kalman gain
  K = ( (Sp * Sp') * H' / Ss') / Ss;
  
  % Corrected model prediction
  x_hat = x_hat + K*y;      % Output state vector
  
  % Update state covariance with new knowledge
  U = (Ss \ H * Sp)';
  
  R = eye(size(Sp));
  for i=1:size(U, 2)
      R = cholupdate(R, U(:,i), '-');
  end
  Sp = R';
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % SR-EKF ends here!
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

