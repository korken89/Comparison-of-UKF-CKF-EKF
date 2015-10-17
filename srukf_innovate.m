function [ x_hat, Sp ] = srukf_innovate( x_in, u_in, z_in, Sp_in, Sq_in, Sr_in, f_function, h_function)  
  
  % UKF settings
  ukf_L = length(x_in);
  ukf_kappa = 3 - ukf_L;
  ukf_alpha = 0.9;
  ukf_beta = 2;
 
  ukf_lambda = ukf_alpha^2*(ukf_L + ukf_kappa) - ukf_L;
  ukf_gamma = sqrt(ukf_L + ukf_lambda);
  ukf_W0_c = ukf_lambda / (ukf_L + ukf_lambda) + (1 - ukf_alpha^2 + ukf_beta);
  ukf_W0_m = ukf_lambda / (ukf_L + ukf_lambda);
  ukf_Wi_m = 1 / (2*(ukf_L + ukf_lambda));
  ukf_Wi_c = ukf_Wi_m;
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % SR-UKF starts here!
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Prediction Update
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  % 1) Evaluate sigma points for i = 0 to m
  Xx = repmat(x_in, 1, length(x_in));
  X_k = [x_in, ( Xx + ukf_gamma * Sp_in ), ( Xx - ukf_gamma * Sp_in )];
  
  % 2) Propagate the sigma points through the nonlinear f transformation
  for i=1:size(X_k, 2)
    X_k(:,i) = f_function(X_k(:,i), u_in);
  end
 
  
  % 3) Estimate the predicted state
  x_hat = ukf_W0_m * X_k(:,1) + ukf_Wi_m * sum(X_k(:,2:end), 2);
  
  % 4) Estimate the square-root factor of the predicted error 
  %    covariance matrix
  Xx = repmat(x_hat, 1, length(x_in)*2);
  [~, R] = qr([sqrt(ukf_Wi_c) * ( X_k(:,2:end) - Xx ), Sq_in]', 0);
  Sp = cholupdate(R, sqrt(ukf_W0_c) * (X_k(:,1) - x_hat), '-')';
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Measurement Update
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % 1) Redraw sigma points for i = 0 to m
  Xx = repmat(x_hat, 1, length(x_in));
  X_k = [x_hat, ( Xx + ukf_gamma * Sp ), ( Xx - ukf_gamma * Sp )];
    
  % 2) Propagate the sigma points through the nonlinear h transformation
  for i=1:size(X_k, 2)
    Z_k(:,i) = h_function(X_k(:,i));
  end
    
  % 3) Estimate the predicted measurement
  z_hat = ukf_W0_m * Z_k(:,1) + ukf_Wi_m * sum(Z_k(:,2:end), 2);

  % 4) Estimate the square-root factor of the corrected innovation 
  %    covariance matrix
  Yy = repmat(z_hat, 1, length(x_in)*2);
  [~, R] = qr([sqrt(ukf_Wi_c) * (Z_k(:,2:end) - Yy), Sr_in]', 0);
  Szz = cholupdate(R, sqrt(ukf_W0_c) * (Z_k(:,1) - z_hat), '-')';
    
  % 5) Estimate the cross covariance matrix
  Xd = X_k - repmat(x_hat, 1, length(x_in)*2 + 1);
  Zd = Z_k - repmat(z_hat, 1, length(x_in)*2 + 1);
  
  Pxz = ( ukf_W0_c*  Xd(:,1) * Zd(:,1)' ) + ...
        ( ukf_Wi_c * Xd(:,2:end) * Zd(:,2:end)' );
    
  % 6) Estimate the Kalman gain
  Wk = (Pxz / Szz') / Szz;
    
  % 7) Estimate the updated state
  x_hat = x_hat + Wk * (z_in - z_hat);

  % 8) Estimate the square-root factor of the corresponding error 
  %    covariance matrix
  U = Wk * Szz;
  
  R = Sp';
  for i=1:size(U, 2)
      R = cholupdate(R, U(:,i), '-');
  end
  Sp = R';
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % SR-UKF ends here!
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

