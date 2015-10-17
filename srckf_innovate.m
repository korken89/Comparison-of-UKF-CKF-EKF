function [ x_hat, Sp ] = srckf_innovate( x_in, u_in, z_in, Sp_in, Sq_in, Sr_in, f_function, h_function)  
  
  % CKF settings
  ckf_m = 2 * length(x_in);
  ckf_xi = sqrt(ckf_m / 2);
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % SR-CKF starts here!
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Prediction Update
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  % 1) Evaluate cubature points for i = 0 to m
  Xx = repmat(x_in, 1, ckf_m/2);
  X_k = [( Xx + ckf_xi * Sp_in ), ( Xx - ckf_xi * Sp_in )];
  
  % 2) Propagate the cubature points through the nonlinear f transformation
  for i=1:size(X_k, 2)
    X_k(:,i) = f_function(X_k(:,i), u_in);
  end
  
  % 3) Estimate the predicted state
  x_hat = 1 / ckf_m * sum(X_k, 2);
  
  % 4) Estimate the square-root factor of the predicted error 
  %    covariance matrix
  Xx = repmat(x_hat, 1, ckf_m);
  [~, R] = qr([1 / sqrt(ckf_m) * ( X_k - Xx ), Sq_in]', 0);
  Sp = R';
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Measurement Update
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % 1) Redraw cubature points for i = 0 to m
  Xx = repmat(x_hat, 1, ckf_m/2);
  X_k = [( Xx + ckf_xi * Sp ), ( Xx - ckf_xi * Sp )];
    
  % 2) Propagate the cubature points through the nonlinear h transformation
  for i=1:size(X_k, 2)
    Z_k(:,i) = h_function(X_k(:,i));
  end
    
  % 3) Estimate the predicted measurement
  z_hat = 1 / ckf_m * sum(Z_k, 2);

  % 4) Estimate the square-root factor of the corrected innovation 
  %    covariance matrix
  Yy = repmat(z_hat, 1, ckf_m);
  [~, R] = qr([1 / sqrt(ckf_m) * (Z_k - Yy), Sr_in]', 0);    
  Szz = R';

  
  % 5) Estimate the cross covariance matrix
  Xd = 1 / sqrt(ckf_m) * (X_k - repmat(x_hat, 1, ckf_m));
  Zd = 1 / sqrt(ckf_m) * (Z_k - repmat(z_hat, 1, ckf_m));
  
  Pxz = Xd * Zd';
    
  % 6) Estimate the Kalman gain
  Wk = (Pxz / Szz') / Szz;
    
  % 7) Estimate the updated state
  x_hat = x_hat + Wk * (z_in - z_hat);

  % 8) Estimate the square-root factor of the corresponding error 
  %    covariance matrix
  [~, R] = qr([Xd - Wk * Zd, Wk * Sr_in]', 0);
  Sp = R';
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % SR-UKF ends here!
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

