function [chi2red, v_vec] = spectra_cross_corr_chi2(wv_tmp,sp_tmp,wv_obs,sp_obs,v_bounds,dv,interp_method,plot_flag)
% ======================================================================================================
% FUNCTION: 
% [corr_vec, v_vec] = spectra_cross_corr(wv_tmp,sp_tmp,wv_obs,sp_obs,hcv_obs,v_max,dv,interp_method,hcv_tmp,vel_tmp)
% ======================================================================================================
%
% This function cross-correlates observed spectra vs. template. See Tonry & Davis (1979).
% The function returns a vector of non-barycentric velocities, and the chi^2 score at each of these velocities.
% Written by Sahar Shahaf, April 2018, as a part of UNICOR.
%
% Last update: ***  By: *** 
%
% -------------------------------------------------------------------------------------------------------
%
% INPUT:
% wv_tmp        - a vector, wavelength vector of the template.
% sp_tmp        - a vector, de-blazed, continuum-normalized, intensity of each wavelength in "wv_tmp" (typically ini physical unitd)
% wv_obs        - a vector, wavelength vector of the observed spectra.
% sp_obs        - a vector, intensity of each wavelength in "wv_obs" (in counts) 
% v_bounds      - a vector, the range in velocity to be fitted [v_min , v_max] - [km/sec]                 (default 100 km/s)
% dv            - a number, the resultion in velocity for cross correlation [km/s] (default 0.1 km/s)
%
% OPTIONAL input:
% interp_method - a string. method of interpolation. 
%                 All Matlab's `interp1` methods are allowed as well (e.g. 'spline' - the default).
% plot_flag     - if true, plot of the chi2 map is generated at the end of
%                 the process, and each step of the fit is plotted as well.
%
% OUTPUT:
% chi2red       - the chi2 reduced vector for each velocity in "v_vec"
% v_vec         - the vector of velocities (baricentric) in km/sec 
%
% ======================================================================================================



% 1) initialize.
%    -----------
if nargin == 6
    interp_method     = 'spline';
    plot_flag         = false;
elseif  nargin == 9 
    interp_method     = 'spline';
    hcv_tmp           = 0;
    vel_tmp           = 0;      
end
 
% speed of light in vacuum (km/sec):
c               = 299792.458; 
 
% Make sure all vectors are column vectors.
wv_tmp           = wv_tmp(:);
sp_tmp           = sp_tmp(:);
wv_obs           = wv_obs(:);
sp_obs           = sp_obs(:);

% 2) Move to log scale on the wavelength & interpulate.
%    ------------------------------------------------
delta_log_lambda = log( 1 + dv / c);
wv_logscale      = log(wv_obs(1)):delta_log_lambda:log(wv_obs(end));
wv               = exp(wv_logscale);

% Moving the data to the log spacing ( this is not a model no special integration is required ) 
sp_obs_LS       = interp1(wv_obs,sp_obs,wv(:),interp_method); 
medCounts       = nanmedian(sp_obs_LS);
NormFactObs     = sqrt( sum((sp_obs_LS - medCounts).^2 ));

sp_tmp_LS       = interp1(wv_tmp,sp_tmp,wv(:),interp_method); 
sp_tmp_LS       = sp_tmp_LS(:) - nanmedian(sp_tmp_LS);
NormFactTmp     = sqrt( sum( (sp_tmp_LS).^2 ) );

sp_tmp_LS       = sp_tmp_LS ./ NormFactTmp * NormFactObs + medCounts;


 
 
% 2) Run chi-squared fitting 
%    -----------------------
% Because we moved to log spacing when we did the interpolation we have
% basiclly evenly spaced in the log. To shift the velocity of the
% template, we simply need to move the template indices. 

shift_min        = round(v_bounds(1) / dv);
shift_max        = round(v_bounds(2) / dv) ;

% Calculate the chi2 - initialize the vectors
shifts           = shift_min : 1 : shift_max;
v_vec1           = shifts*(dv);
chi2red          = NaN(size(v_vec1));     %The vector of chi squared reduced.
template_shifted = NaN(size(sp_tmp_LS));

figure
for i  = 1:length(chi2red)
    
    shift            = shifts(i);
    % If the shift is towards red wavelengths (positive velocity) 
    RedShift         = max(1 + shift, 1);        
    
    % If the shift is towards blue wavelengths (positive velocity) 
    BlueShift        = abs(min(shift, 0)) + 1;
    
    % Shift the template
    template_shifted(RedShift:end-BlueShift+1) = sp_tmp_LS(BlueShift:end-RedShift+1);
    DataIdx                                    = ~isnan(template_shifted);
    template_shifted(~DataIdx)                 = [];
    template_shifted                           = template_shifted(:);
    
    % Cut the observed spectrum and the wavelength vecotor:
    Ydata     = sp_obs_LS(DataIdx);
    Xdata     = wv(DataIdx)';
    Err       = sqrt(abs(Ydata));   
    
    % Define the model:
    X0        = median(Xdata);
    ModelX    = [ones(size(template_shifted)), template_shifted,...
                    template_shifted.*(Xdata-X0), template_shifted.*(Xdata-X0).^2,...
                    template_shifted.*(Xdata-X0).^3, template_shifted.*(Xdata-X0).^4,...
                    template_shifted.*(Xdata-X0).^5, template_shifted.*(Xdata-X0).^6];
    
    % Create the design matrix:    
    DesignM   = ModelX./ repmat(Err,[1,size(ModelX,2)]);
    
    % define the normalized observations vector. Assuming poisson noise, this is the square root of the observation:
    Bvector   = Ydata./Err;
    
    % calculate the covariance matrix and solve the normal equation:
    errsM     = (DesignM'*DesignM)^(-1);
    solution  = errsM*DesignM'*Bvector;
    model     = ModelX*solution;
    

    chi2      = nansum( (model- Ydata).^2 ./ Err.^2 );
    DoF       = length(template_shifted)-size(ModelX,2);
    chi2red(i)= chi2/DoF;
    
    
    subplot(2,1,1)
    stairs(Xdata,Ydata,'k'); hold on; grid on;
 
    plot(Xdata,model,'r','linewidth',2);
    xlabel('\lambda','fontsize',16);
    title([num2str(dv*shift) ' ; ' num2str(solution') ]);
    subplot(2,1,2)
    plot(Xdata,Ydata-model,'.k'); hold on; grid on;
    %pause(0.15);
    pause;
    clf
%     
    template_shifted = NaN(size(sp_tmp_LS));

end


% 3) Heliocentric correction
%    -----------------------
v_vec            = v_vec1;% + hcv_obs + vel_tmp + hcv_tmp;

% Transform to column vectors
if min(chi2red) > 0
chi2red          = chi2red(:)/min(chi2red);
end
v_vec            = v_vec(:);

figure;
plot(v_vec,chi2red,'k','linewidth',2); hold on; grid on;
plot([min(v_vec),max(v_vec)],[1+(1/length(template_shifted)),1+(1/length(template_shifted))],'r','linewidth',2)










