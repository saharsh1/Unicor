function [corr_vec, v_vec] = spectra_cross_corr(wv_tmp,sp_tmp,wv_obs,sp_obs,hcv_obs,v_max,dv,interp_method,hcv_tmp,vel_tmp)
% ======================================================================================================
% FUNCTION: 
% [corr_vec, v_vec] = spectra_cross_corr(wv_tmp,sp_tmp,wv_obs,sp_obs,hcv_obs,v_max,dv,interp_method,hcv_tmp,vel_tmp)
% ======================================================================================================
%
% This function cross-correlates observed spectra vs. template. See Tonry & Davis (1979).
% The function returns a vector of barycentric velocities, and the CCF score at each of these velocities.
% Written by Sahar Shahaf, April 2018, as a part of UNICOR.
%
% Last update: ***  By: *** 
%
% -------------------------------------------------------------------------------------------------------
%
% INPUT:
% wv_tmp        - a vector, wavelength vector of the template.
% sp_tmp        - a vector, de-blazed, continuum-normalized, intensity of each wavelength in "x_tmp"
% wv_obs        - a vector, wavelength vector of the observed spectra.
% sp_obs        - a vector, de-blazed, continuum-normalized, intensity of each wavelength in "x_obs"
% hcv_obs       - a number, heliocntric correction for the observed data [km/s].
% 
% OPTIONAL input:
% v_max         - a number, the maximum shift in velocity [km/sec]                 (default 100 km/s)
% dv            - a number, the resultion in velocity for cross correlation [km/s] (default 0.1 km/s)
% interp_method - a string. method of interpolation. 
%                 'trigo' for trigonometric interpolation (preliminary, in construction). 
%                 All Matlab's `interp1` methods are allowed as well (e.g. 'spline').
% hcv_tmp       - a number, the heliocntric correction for the template [km/s].
% vel_tmp       - a number, the velocity of the template [km/s].
%
% OUTPUT:
% corr_vec      - the correlation vector for each velocity in "v_vec"
% v_vec         - the vector of velocities (baricentric) in km/sec 
%
% ======================================================================================================



% 1) initialize.
%    -----------
if nargin == 5
    v_max             = 100;
    dv                = 0.1;
    interp_method     = 'spline';
    hcv_tmp           = 0;
    vel_tmp           = 0;
elseif  nargin == 7  
    interp_method     = 'spline';
    hcv_tmp           = 0;
    vel_tmp           = 0;      
end

% !!! - TODO: add a warning if median of spectra are not normalized (~=1) & divide by median. 


% speed of light in vacuum (km/sec):
c               = 299792.458; 
% Translating max shift to number of points
maxlags         = floor(v_max / dv);

% Make sure all vectors are column vectors.
wv_tmp           = wv_tmp(:);
sp_tmp           = sp_tmp(:);
wv_obs           = wv_obs(:);
sp_obs           = sp_obs(:);

% If the observed spectra covers a larger wavelength range than the model -
% cut it to fit the model.
if min(wv_obs) < min(wv_tmp) || max(wv_obs) > max(wv_tmp) ;
    Indx  =  find(wv_obs<min(wv_tmp),1,'last');
    Indx2 =  find(wv_obs>max(wv_tmp),1,'first');
    wv_obs(Indx2-40:end) = []; sp_obs(Indx2-40:end) = [];
    wv_obs(1:Indx+40)    = []; sp_obs(1:Indx+40)  = [];
   
end

% 2) Move to log scale on the wavelength & interpulate.
%    ------------------------------------------------
delta_log_lambda = log( 1 + dv / c);
wv_logscale      = log(wv_obs(1)):delta_log_lambda:log(wv_obs(end));
wv               = exp(wv_logscale);

% Moving the data to the log spacing ( this is not a model no special integration is required )
if ~strcmp(interp_method,'trigo') 
    sp_obs_LS       = interp1(wv_obs,sp_obs,wv(:),interp_method); 
    sp_obs_LS       = sp_obs_LS(:);
    sp_tmp_LS       = interp1(wv_tmp,sp_tmp,wv(:),interp_method); 
    sp_tmp_LS       = sp_tmp_LS(:);
else
    sp_obs_LS       = triginterp(wv(:),wv_obs,sp_obs); 
    sp_obs_LS       = sp_obs_LS(:);
    sp_tmp_LS       = triginterp(wv(:),wv_tmp,sp_tmp); 
    sp_tmp_LS       = sp_tmp_LS(:);
end

% Subtract median flux.
sp_obs_LS       = (sp_obs_LS - median(sp_obs_LS));  
sp_tmp_LS       = (sp_tmp_LS - median(sp_tmp_LS));

 
% 2) Run cross-correlation 
%    ---------------------- 
% Because we moved to log spacing when we did the interpolation we have
% basiclly evenly spaced in the log, so now we simply run cross corr.

% Running the cross-correlation
[corr_vec, lags] = xcorr(sp_obs_LS,sp_tmp_LS,maxlags,'coeff');
v_vec1           = lags*((delta_log_lambda)*c);

% 3) Heliocentric correction
%    -----------------------
v_vec            = v_vec1 + hcv_obs + vel_tmp + hcv_tmp;

% Transform to column vectors
corr_vec         = corr_vec(:);
v_vec            = v_vec(:);











