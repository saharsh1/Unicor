function [v,CCF,corr_max,corr_RMS] = calc_multi_order_velocity(spect,template,par)
% ========================================================================== 
% FUNCTION: [v,CCF] = calc_multi_order_velocity(spect,template,par)
% ========================================================================== 
% 
% This function recieves multi-order spectra and a theoretical multi-order template,
% and returns velocity and cross-correlation function per order, for each observation.
% Written by Sahar Shahaf, April 2018, as a part of UNICOR.
% 
% Last update: 13/9/2018     By: Sahar
%
% --------------------------------------------------------------------------
%
% INPUT:
% spect    - data structre containing the observed spectrum. 
%            For multiple exposures, the input is a cell array of
%            structures in the following format: 
%            spect.wv      - wavelength matrix      [ lambda , order ] 
%            spect.sp      -  de-blazed, continuum-normalized,  
%                             intensity of each wavelength in spect.wv [ flux   , order ] 
%            spect.bcv     - barycentric correction in km/s.
%            spect.bjd     - time of observation.
%
%            for example, for the i-th exposure
%                         spect{i}.wv is the i-th wavelength matrix and 
%                         spect{i}.sp is the i-th spectrum.
%            if only one exposure is analyzed, then 
%                         spect.wv is the wavelength matrix and 
%                         spect.sp is the spectrum.
%
% template - contrary to spect, template is not a cell array, but a single
%            data structure, containing the theoreticl template.
%            template.wv   - wavelength matrix      [ lambda , order ] 
%            template.sp   -  de-blazed, continuum-normalized,  
%                             intensity of each wavelength in spect.wv [ flux   , order ] 
%            template.bcv  - barycentric correction in km/s.
%            template.vel  - velocity of the theoretical template in km/s
% par      - a structure of parameters  (see detailes below)
%
% OUTPUT:
% v        - A matrix of derived velocities, per order [ # obs , # order] 
% CCF      - a structure, contining four fields:
%            the CCF value("corr") at velocities("vels"), 
%            and the parameters of the parabola fitted
%            to the peak ("poly_coeffs").
% corr_max - A matrix containing the CCF peaks ([Nobs,Nord])
% corr_RMS - A matrix containing the CCF RMS (1,48*MAD), ([Nobs,Nord]).
%--------------------------------------------------------------------------
%
% Parameters: 
% -----------
% parameters for spectra_cross_corr.m:
% par.v_max and par.dv determine velocity vector for which correlation is calculated. 
% par.v_max         - maximum shift in velocity [km/sec]
% par.dv            - delta-v [km/s]
% par.interp_method - a string. method of interpolation. 
%                    'trigo' for trigonometric interpolation (preliminary, in construction). 
%                     All other interp1 methods are allowed as well (e.g. 'spline').
% parameters for extract_max_corr.m:
% par.Nparab        - number of points in the parabolic fit.
%
% --------------------------------------------------------------------------

    
% 1) Initialize.
% --------------

% Set default if no input was given.
if nargin == 2
    par.v_max         = 500;
    par.dv            = 0.1;
    par.Nparab        = 10;
    par.interp_method = 'spline';
elseif nargin ==3  
    if ~isfield(par,'Nparab')
        par.Nparab  = 10;
    end
    if ~isfield(par,'interp_method')
        par.interp_method = 'spline';
    end
end

% If input data is struct, convert it into a cell array.
if isstruct(spect)
    tmpVar   = spect;  clear('spect');  
    spect{1} = tmpVar; clear('tmpVar');
end

Nobs                 = length(spect);
Nord                 = size(spect{1}.wv,2);
v                    = NaN(Nobs,Nord);
CCF                  = cell(Nobs,1);
corr_max             = zeros(Nobs,Nord);
corr_RMS             = ones(Nobs,Nord);

% 2) Run multi-order cross correlation for each observation   
% ---------------------------------------------------------
for Iobs = 1 : Nobs
fprintf(['Obs #' int2str(Iobs) '. ']);
% initiailze variables.
    CCF_score            = zeros(2*floor(par.v_max/par.dv) + 1,Nord);
    CCF_vels             = zeros(2*floor(par.v_max/par.dv) + 1,Nord);

    vels                 = zeros(1,Nord);
    poly_coeffs          = zeros(3,Nord);


% 2.1) Calculate the cross-correlation per order.
% -----------------------------------------------
   
    fprintf('Calculating order:  ');
    for Iord = 1:Nord
        if mod(Iord,10)==0
        fprintf([int2str(Iord) ' ']);
        end
        wvOBS          = spect{Iobs}.wv(~isnan(spect{Iobs}.sp(:,Iord)),Iord);
        spOBS          = spect{Iobs}.sp(~isnan(spect{Iobs}.sp(:,Iord)),Iord);
  
        wvTMPL         = template.wv(~isnan(template.sp(:,Iord)),Iord);
        spTMPL         = template.sp(~isnan(template.sp(:,Iord)),Iord);
   
        [CCF_score(:,Iord), CCF_vels(:,Iord)]                = spectra_cross_corr(wvTMPL,spTMPL,wvOBS,spOBS,spect{Iobs}.bcv,par.v_max,par.dv,par.interp_method,template.bcv,template.vel);
        [vels(Iord),corr_max(Iobs,Iord),poly_coeffs(:,Iord)] = extract_max_corr(CCF_vels(:,Iord),CCF_score(:,Iord),par.Nparab);
        
        corr_RMS(Iobs,Iord) = 1.48*mad(CCF_score(:,Iord),1);
    end

% 2.2) Arrange output
% -------------------
CCF{Iobs}.corr        = CCF_score;
CCF{Iobs}.vels        = CCF_vels;
% CCF{Iobs}.peak        = corr_max;
% CCF{Iobs}.RMS         = corr_rms;
CCF{Iobs}.poly_coeffs = poly_coeffs;

v(Iobs,:)             = vels;
fprintf(' \n');



end

