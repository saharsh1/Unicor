function [template] = vsini_matrot(template,vsini,epsilon)
% [nwl,rotated] = matrot1(wl,spec,vsini,epsilon)
% Creates roteted spectrum with Vsini and epsilon.
% 
% Input: wl      - the wavelength vector of the spectrum to rotate
%        spec    - the flux vector of the spectrum to rotate
%        vsini   - in km/s.
%        epsilon - linear limb darkenning (between 0 and 1).
% 
% use epsilon=0.6 for G stars
% Rotationally broaden spectra.
%
% Created: Shay Zucker
% Modified: 12.08.2010 LT
%
c = 299792.46;
oversamplingfactor = 4;
windowfactor = 1;
beta = vsini/c;

wl   = template.wv';
spec = template.sp';
%
% Interpolate to log-lambda space.
%
N = size(wl);
deltalog = log(wl(end,:)./wl(1,:))/(N(1)-1);

nwl     = zeros(N(1),N(2));
intspec = zeros(N(1),N(2));

for i = 1:length(deltalog)
   tmp = exp(log(wl(1,i)):deltalog(i):log(wl(end,i)));
   nwl(:,i) = tmp;
   nwl(1,i) = wl(1,i);
   nwl(end,i) = wl(end,i);
   
   intspec(:,i) = interp1(wl(:,i),spec(:,i),nwl(:,i),'linear'); % 22.4.10 LT: is it the best interpolation method?
    
end


fintspec = fft(intspec);

%
% We want the pixel size to be smaller by at least "oversamplingfactor" than beta.
%
resamplingfactor = ceil(deltalog/(beta/oversamplingfactor));
ndeltalog = deltalog/resamplingfactor;

if (resamplingfactor > 1)
    if (rem(N,2) == 0)
        nfintspec = [fintspec(1:(N/2)) fintspec(N/2+1)/2 zeros(1,N*(resamplingfactor-1)-1) ...
                     fintspec(N/2+1)/2 fintspec((N/2+2):N)];
    else
        nfintspec = [fintspec(1:((N+1)/2)) zeros(1,N*(resamplingfactor-1)) fintspec(((N+3)/2):N)];
    end 
else
    nfintspec = fintspec;
end
nintspec = real(ifft(nfintspec))*resamplingfactor;

%
% Prepare the broadening profile
%
edge = floor((beta/ndeltalog)*windowfactor); 
l = (-edge:edge)*ndeltalog; % the window
G = (2*(1-epsilon)*(1-(l/beta).^2).^0.5 + 0.5*pi*epsilon*(1-(l/beta).^2))/(pi*beta*(1-epsilon/3))*ndeltalog;
% locally normalized: 
G = G./sum(G);

%figure;
%plot(l,G,'.');

%
% Add values to the vector to take care of the limits.
%
nnintspec = [nintspec(1)*ones(1,edge) nintspec(:)' nintspec(end)*ones(1,edge)];

% 
% Create the rotated spectrum:
%
rotated = conv(G,nnintspec);
%
% return the wavelength to its original length:
%
rotated = rotated((2*edge+1):(length(rotated)-2*edge));
rotated = rotated(1:resamplingfactor:end);


template.sp = rotated;
template.wv = nwl;