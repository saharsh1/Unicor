function [obs] = SpecPreProcess(obs,FiltPar,BellPar,ExcludeOL)

% This function prepares observed spectra for cross-correlation analysis
% using 'SpecPreProc_main' of Unicor2.
% Description of the procedure: 
% 
% 1)     The spectra is resampled and interpolated over a uniform wavelength grid.
% 2)     The uniformly sampled spectra is filtered using a Kaiser window. 
%   2.1) A high pass Kaiser window (producing the spHP output variable).
%   2.2) A low pass Kaiser window (producing the spLP output variable).
%        Note that spHP and spLP are both computed directly from the resampled spectra. 
% 3)     The low-pass response, spLP, approximates the instrumental response
%        function of the spectrograph (the Blaze). The de-Blazed spectra (spDBLZ)
%        is calculated by dividing the high pass response, spHP, by the Blaze, spLP.
% 4)     A Tukey (tapered cosine) window is applied on the high-pass and the deblazed spectra.
%
%
% INPUTS:
%       obs    A cell array. Each cell containing data structure for each measurement with the following fields:
%                   wv - wavelength matrix
%                   sp - spectra matrix (colomn 'i' of 'sp' corresponds the same colomn in 'wv')
%                   jd - JD of observation
%                   hcv - heliocentric correction velocity
%                   name - star name
%                   filename - the path of the measurement raw files
%                   bcv - baricentric correction velocity 
%        FiltPar  - structure. Parameters of the Kaiser window filter (see details below).
%        BellPar  - a number. if 0 < BellPar <=1  the Tukey (tapered cosine) window is applied. 
%                   A Tukey window is a rectangular window with the first and last BellPar/2 
%                   percent of the samples equal to parts of a cosine. Defalt: 0.2
%        ExcludeOL- structure of the exclude outliers. 
%                    ExcludeOn - Determines wheather to run the procedure. 
%                    Npix2filt - Number of pixels to be used in the median filter.
%                    Nsig      - Sigma threshold.
% 
% OUTPUT:
%       obs    A cell array. Each cell containing data structure for each measurement with the following fields:
%                   wv - wavelength matrix
%                   sp - Preprocessed spectra matrix (colomn 'i' of 'sp' corresponds the same colomn in 'wv')
%                   jd - JD of observation
%                   hcv - heliocentric correction velocity
%                   name - star name
%                   filename - the path of the measurement raw files
%                   bcv - heliocentric correction velocity 
% ---------------------------------------------------------------------------------------------
%
% Parameters:
%     Interpolation:
%     FiltPar.StepScale  - a number. The spacing of the wavelength vector, used for interpolation, in units of
%                          the median wavelengh step of the input vector wv_in. Default: 0.5
%     
%     Bandpass filtering:
%     The basic frequecies here are measured as the basic smallest frequency, in usints of the Nyquist frequency: 0.5 / (# of datapoints).
%     FiltPar.NfStop      - The stopping frequency, below which the high pass filter will remove all frequencies.  Default: 1                       
%     FiltPar.NfPass      - The passing frequency, above which the high pass filter will retain all frequencies.   Default: 75                           
%     FiltPar.Ripple      - The maximum ripple  (in dB) for in the bandpass.  Default: 0.5
%     FiltPar.Attenuation - The minimal attenuation (in dB) outside thebandpass.  Default: 75

obsA = obs;
N = 0;
Size = 1e4;

% This part implements SpecPreProc_main on the observed spec, taking into account the
% changing size of the output spec.
for i = 1:length(obs)
    
    WVA = zeros([size(obs{1,i}.wv,2),Size]); 
    SPA = WVA;
    
    WV = obs{1,i}.wv;
    SP = obs{1,i}.sp;
    
    for j = 1:size(obs{1,i}.wv,2)
        clc
        fprintf(['PreProc Obs#' num2str(i)  ' Ord#' int2str(j)]);
        
        if ExcludeOL.ExcludeOn
            [wvA,spA,~,~]=exclude_outliers(WV(:,j),SP(:,j),ExcludeOL.Npix2filt,ExcludeOL.Nsig);
        end
        
        [wvA,spA,~,~] = SpecPreProc_main(wvA,spA,false,FiltPar,BellPar);
        
      
        Delta = median(diff(wvA));
        WVA(j,:) = [wvA (max(wvA)+Delta):Delta:(max(wvA)+Delta*(Size-length(wvA)))];
        SPA(j,1:length(spA)) = spA;
        
        N = max(N,length(wvA));
    end
    
    obsA{1,i}.wv = WVA';
    obsA{1,i}.sp = SPA';
    
end

% Removing zero cells
for i = 1:length(obs)
      
    M1 = obsA{1,i}.wv;
    M2 = obsA{1,i}.sp;
    
    M1((N+1):end,:) = [];
    M2((N+1):end,:) = [];
    
    obsA{1,i}.wv = M1;
    obsA{1,i}.sp = M2;
    
end
 

obs = obsA;
end
