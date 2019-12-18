function [wv,spHP,spLP,spDBLZ ] = SpecPreProc_main( wv_in,sp_in,PlotFlag,FiltPar,BellPar)
% FUNCTION: [wv,spHP,spLP,spDBLZ ] = SpecPreProc_main( wvin,spin,PlotFlag,FiltPar,BellPar)
%
% This function prepares observed spectra for cross-correlation analysis. 
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
% INPUT: wvin     - a vector. The sampled wavelength.
%        spin     - a vector. The measured spectra at wvin.
%        PlotFlag - logical. If true the function plots the resulting data. Default: false.
%        FiltPar  - structure. Parameters of the Kaiser window filter (see details below).
%        BellPar  - a number. if 0 < BellPar <=1  the Tukey (tapered cosine) window is applied. 
%                   A Tukey window is a rectangular window with the first and last BellPar/2 
%                   percent of the samples equal to parts of a cosine. Defalt: 0.2
% 
% OUTPUT: wv      - a vector. The evenly sampled wavelength vector.
%         spHP    - a vector. High-pass filtered spectra.
%         spLP    - a vector. Low-pass filtered spectra, used for deblazing.
%         spDBLZ  - a vector. Deblazed spectra.
%
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

%
% 1) Initialize
% -------------
fprintf('\n');
fprintf('SpecPreProc report: \n\n');
if nargin < 3
    PlotFlag   = false;
end

if nargin < 4
    FiltPar.StepScale   = 0.5;
    FiltPar.NfStop      = 1; 
    FiltPar.NfPass      = 75;   
    FiltPar.Ripple      = 0.25;       
    FiltPar.Attenuation = 75;    
    fprintf('No Bandpass filter input. Operating according to default settings. \n');
end

if nargin < 5
    BellPar = 0.2;
    fprintf('No cosine-bell input. Operating according to default settings. \n');
end

if isstruct(FiltPar)
if ~isfield(FiltPar,'Ripple')
    FiltPar.Ripple = 0.25;
    fprintf('Missing ''Ripple'' parameter. Operating according to default input (filter: kaiserwin).  \n');
end

if ~isfield(FiltPar,'Attenuation')
    FiltPar.Attenuation = 75;
    fprintf('Missing ''Attenuation'' parameter. Operating according to default input (filter: kaiserwin).  \n');
end
end
% 2) Interpolate the spectra to obtain uniform sampling.
% ------------------------------------------------------ 
span       = max(wv_in) - min(wv_in);
step       = FiltPar.StepScale*median(diff(wv_in)); 
Nsteps     = ceil( span/step );
wv         = linspace(wv_in(1), wv_in(end), Nsteps);
sp         = interp1(wv_in,sp_in,wv,'spline');
fprintf([' * Spectra interpolated (evenly spaced). Number of points:' num2str(Nsteps) '\n']);


% 3) Perform HP / LP filter.
% --------------------------
df         = 0.5/Nsteps;
fstop      = FiltPar.NfStop*df;
fpass      = FiltPar.NfPass*df;

% High pass filter:
hpFilt   = designfilt('highpassfir','StopbandFrequency',fstop, ...
           'PassbandFrequency',fpass,'PassbandRipple',FiltPar.Ripple, ...
           'StopbandAttenuation',FiltPar.Attenuation,'DesignMethod','kaiserwin');
spHP = filtfilt(hpFilt,sp);     
fprintf([' * High-pass kaiserwin filter applied. \n'...
         '   Stopband = %f, Passband = %f, Ripple = %f, Attenuation = %f.\n'],...
         fstop,fpass,FiltPar.Ripple,FiltPar.Attenuation);

% Low pass filter (notice that fstop & fpass switched places)     
lpFilt   = designfilt('lowpassfir','StopbandFrequency',fpass, ...
           'PassbandFrequency',fstop,'PassbandRipple',FiltPar.Ripple, ...
           'StopbandAttenuation',FiltPar.Attenuation,'DesignMethod','kaiserwin');
spLP = filtfilt(lpFilt,sp);
fprintf([' * Blaze calculated from low-pass kaiserwin filter. \n'...
         '   Stopband = %f, Passband = %f, Ripple = %f, Attenuation = %f.\n'],...
         fpass,fstop,FiltPar.Ripple,FiltPar.Attenuation);


% Calculate the deblazed spectra
spDBLZ              = spHP./spLP;
spDBLZ              = spDBLZ - median(spDBLZ);

% I1                  = find(sp_in > 75,1,'first');
% I2                  = find(sp_in > 75,1,'last');
% try
%     ind                 = wv < wv_in(I1) | wv > wv_in(I2);
%     spDBLZ(ind)         = 0;
% catch
%     spDBLZ              = zeros(size(spDBLZ));
% end
fprintf([' * De-Blzed spectrum calculated.\n'...
         '   WARNING: De-Blazed spectra are bad for your health! Use with caution!\n']);

% 4) Apply Cosine Bell
% --------------------
if BellPar > 0 && BellPar <= 1
    CosBell   = tukeywin(length(wv),BellPar);
    spDBLZ    = CosBell' .* spDBLZ;
    spHP      = CosBell' .* spHP;
    fprintf([' * Cosine Bell applied. Tukeywin ratio parameter = ' num2str(BellPar)  '.\n']);
else
    fprintf(' * Cosine Bell was not applied. \n');
end

if PlotFlag
 fvtool(lpFilt);
 fvtool(hpFilt);
end
 
%if PlotFlag
    %figure
%     subplot(3,1,1)
%     plot(wv,spLP,'r'); hold on; grid on;
%     plot(wv,sp,'k'); hold on; grid on;
%     title('Raw spectrum & Low-pass response (Blaze)');
%     xlabel('\lambda','fontsize',16);
%     ylabel('Counts','fontsize',14);
% 
%     subplot(3,1,2);
%     plot(wv,spHP,'k'); hold on; grid on;
%     title('High-pass response (Recommended)');
%     xlabel('\lambda','fontsize',16);
%     ylabel('Counts','fontsize',14);
% 
%     subplot(3,1,3);
%     plot(wv,spDBLZ,'k'); hold on; grid on;
%     title('De-Blazed spectrum');
%     xlabel('\lambda','fontsize',16);
%     ylabel('Normalized Flux','fontsize',14);
% else 
%     fprintf(' * Plot not generated by SpecPreProc_main. \n');
% end

end

