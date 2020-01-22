function [template] = inst_broad(template,R) 
     
%    Generate instrumental broadening
%    ------------------------------------------------------------------
%     Inputs:
%     template - a structure that contains wv vector (template.wv) and spec
%     vector (template.sp)
%     inst_broaden - instrumental broadening [Angstrems] at the central wavelength.
%     Outputs:
%     template - broadened template

%     Interpolating the spectrum to an equally spaced spectrum in log(Lambda) axis.
      LL_diff  = min(diff(log(template.wv)))/10; % minimal step of Log(Lambda)
      X        = min(log(template.wv)) : LL_diff : max(log(template.wv)); % Equally spaced "log-lambda" vector.
      SP       = spline(log(template.wv),template.sp,X); % The interpolated template (which is Equally spaced in "log-lambda")

%     Create the Gaussian for the convolution:
%     Sigma is calculated as shown below. area under gaussian should be
%     normlized to 1 with miu = 0;
      
%     inst broaden is the FWHM. convert to sigma and to v/c.
      broaden_sigma        = (1/R)/2.35482;      
      broadening_gaussian  = GaussianCurve((-2*broaden_sigma:LL_diff:2*broaden_sigma),[broaden_sigma 0]); % Gaussian filter in log-lambda axis.
      normalization_factor = sum(broadening_gaussian);
      broadening_gaussian  = broadening_gaussian/normalization_factor;
      
%     Convolve the gaussian with the template 
      SP                   = conv(SP,broadening_gaussian,'same'); % Convolution in log-lambda axis LL_diff*
      template.sp          = spline(X,SP,log(template.wv)); % Interpolating the result to template.wv 
      template.name        = [template.name '-GB' num2str(R)];
           
end