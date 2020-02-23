function [template] = TempPreProcess(template,OBSwv,DeltaV,FiltPar,BellPar)
% This function prepares template spectra for cross-correlation analysis
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
%       template  - A cell array. Each cell containing data structure for each measurement with the following fields:
%                   wv - wavelength matrix
%                   sp - spectra matrix (colomn 'i' of 'sp' corresponds the same colomn in 'wv')
%                   name - star name
%                   bcv - baricentric correction velocity (set to 0)
%                   vel - set to 0.
%        FiltPar  - structure. Parameters of the Kaiser window filter (see details below).
%        BellPar  - a number. if 0 < BellPar <=1  the Tukey (tapered cosine) window is applied. 
%                   A Tukey window is a rectangular window with the first and last BellPar/2 
%                   percent of the samples equal to parts of a cosine. Defalt: 0.2
%        Limits   - 2xn matrix (n=number of orders). Each colomn contains
%                   the minimal (first element) and maximal (second element)
%                   wavelength of a specific order.
% 
% OUTPUT:
%       template  - A cell array. Each cell containing data structure for each measurement with the following fields:
%                   wv - wavelength matrix
%                   sp - Preprocessed spectra matrix (colomn 'i' of 'sp' corresponds the same colomn in 'wv')
%                   name - star name
%                   bcv - barycentric correction velocity (set to 0)
%                   vel - set to 0.
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


% This part implements SpecPreProc_main on the template spec, taking into account the
% changing size of the output spec.

c = 299792.458; % Speed of light, km/s.

% Setting limits:
Limits      = [min(OBSwv) ; max(OBSwv)];
Limits(1,:) = Limits(1,:) - DeltaV/c*mean(Limits);
Limits(2,:) = Limits(2,:) + DeltaV/c*mean(Limits);

Size = 5e5; N=0;
WVA = zeros(Size,size(Limits,2)); 
SPA = WVA;


OrdN   = size(Limits,2);
Nvec   = zeros(1,OrdN);
wvCell = cell(OrdN,1);
spCell = cell(OrdN,1);

for i = 1:OrdN
        clc
        fprintf(['Temp PreProc Ord#' int2str(i)]);
        I             = (template.wv > Limits(1,i)) & (template.wv < Limits(2,i));
        wvTmp         = template.wv(I);
        spTmp         = template.sp(I);
        [wvTmp,spTmp] = SpecPreProc_main(wvTmp,spTmp,false,FiltPar,BellPar);
       
        wvCell{i}     = wvTmp;
        spCell{i}     = spTmp;
        Nvec(i)       = max(length(wvTmp) );
end

Nmax    = max(Nvec);
wv      = zeros(Nmax,OrdN);
sp      = zeros(Nmax,OrdN);

clc
for i = 1:OrdN
        
    
        wv(1:Nvec(i),i)     = wvCell{i};
        FillUpVec           = linspace(max(wvCell{i})+0.1,max(wvCell{i})+1.1,Nmax-Nvec(i));
        wv(Nvec(i)+1:end,i) = FillUpVec;
        sp(1:Nvec(i),i)     = spCell{i};
  
end

template.wv = wv;
template.sp = sp;



end
