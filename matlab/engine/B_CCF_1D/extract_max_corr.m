function [vel, max_corr , poly_coeffs] = extract_max_corr(v_vec,corr_vec,N)
% ======================================================================================================
% FUNCTION: [max_x, max_corr , poly_coeffs] = extract_max_corr(corr_vec,x_vec,N)
% ======================================================================================================
% This function fits a parabula to the peak of the CCF function.
% The function returns the velocity at the peak, it's CCF value, and the
% coefficients of the fitted parabula. 
%
% Written by Sahar Shahaf, April 2018, as a part of UNICOR.
% Last update: ***  By: *** 
%
% -------------------------------------------------------------------------------------------------------
%
% INPUT:
% v_vec       - a vector, the velocity vector.
% coor_vec    - a vector, the CCF value for each velocity.
% N           - Number of points to use in the parabolic fit  (default = 10).
%
%  
% OUTPUT:
% vel         - a number. The velocity at the peak of the fitted parabula.
% max_corr    - a number. The maximal CCF value (peak of the fitted parabula).
% poly_coeffs - a vector. The coefficients of the fitted parabula.  
%
% -------------------------------------------------------------------------------------------------------
%

% 1) initialize:
% --------------
if nargin == 2
    N = 3;
end

N_half = floor(N/2);

% 2) Fit a parabola to the peak:
% -----------------------------

% Finding the highest point
[~, max_ind] = max(corr_vec);

if max_ind-N_half<1
    max_ind  = N_half + 1;
elseif max_ind+N_half>length(corr_vec)
    max_ind  = max_ind - (N_half+1);
end
 
% Taking N points around the maximum
x                = v_vec(max_ind-N_half:max_ind+N_half);
y                = corr_vec(max_ind-N_half:max_ind+N_half);

% Fitting a polynom
XX               = [x(:) x(:).^2];
try
[poly_coeffs, ~] = robustfit(XX,y,'ols');
catch
   vel         = NaN; 
   max_corr    = NaN;
   poly_coeffs = NaN;
   return
end
% Finding the maximum
vel              = -poly_coeffs(2)/(2*poly_coeffs(3));
max_corr         = poly_coeffs(2)^2/(4*poly_coeffs(3)) - poly_coeffs(2)^2/(2*poly_coeffs(3)) + poly_coeffs(1);


end

