function [output_wave,output_spect,sigma,excluded_wave]=exclude_outliers(wave,spectrum,Npix2filt,Nsig)
% function to find the outliers, exclude them and replace them
% by the median of closest neighbors.

% median filtering and calculate residus
fspect = medfilt1(spectrum,Npix2filt);
% for i=1:floor(Npix2filt/2)
%     fspect(i)       = nanmedian(spectrum(1:i+floor(Npix2filt/2)));
%     fspect(end-i+1) = nanmedian(spectrum(end-i+1-floor(Npix2filt/2):end));
% end
res = spectrum-fspect;

% sigma
sigma = 1.48*mad(res,1);

% find pixel above and below the threshold
pix = find(res<=Nsig*sigma);
exc = res > Nsig*sigma;


% replace the excluded values by the linear interpolation
output_spect  = interp1(wave(pix),spectrum(pix),wave,'linear');

% clean nan values
output_wave   = wave(~isnan(output_spect));
excluded_wave = wave(exc & ~isnan(output_spect));
output_spect  = output_spect(~isnan(output_spect));

end

