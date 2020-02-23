function [template] = GetTemplate(Dir,TemplPars,SaveFlag)
% This function loads a PHOENIX synthetic template.
% If the template cannot be found on the supplied template directory (Dir),
% it is downloaded from the FTP of PHOENIX project. 
%
% INPUT:
% Dir       - string. The template-storing directory.
% TemplPars - structure. Continas the tempalte parameters (Teff,logg,z).
% SaveFlag  - logical. Determined whether to save the downloaded template.
%                      (defualt=false)
%
% OUTPUT:
% template - structure: 
%                   wv   - wavelength vector.
%                   sp   - synthetic PHOENIX spectrum.
%                   name - the name of the template.
%                   vel  - the velocity of the template (set to be 0).
%                   bcv  - barycentric correction (set to be 0).

% Initialize:
% -----------
if nargin < 3
    SaveFlag = false;
end

% Load PHOENIX template wavelength vector:
if exist(fullfile(Dir,'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits'),'file')
    wv = fitsread(fullfile(Dir,'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits'));
    DownloadWV = false;
else
    DownloadWV = true;
end

% Set Phoenix template name:
TemplName  = ['lte' num2str(TemplPars.Teff,'%05d') ...
              '-' num2str(TemplPars.logg,'%3.2f') ...
              num2str(TemplPars.z-1e-3,'%+2.1f')...
              '.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'];

% Load the PHOENIX template. 
% If the template is not available in the local machine (in the directory
% supplied to the code, Dir, it is possible to download the template from
% the PHOENIX ftp (only if DownloadFlag is true).
if exist(fullfile(Dir,TemplName),'file')
    sp = fitsread(fullfile(Dir,TemplName));
else

        
    % Setup an FTP connection.
    PHOENIXftp    = 'phoenix.astro.physik.uni-goettingen.de';
    PHOENIXfolder = ['/v2.0/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/'...
                     'Z' num2str(TemplPars.z-1e-3,'%+2.1f/')];
    ftpobj        = ftp(PHOENIXftp);
    pasv(ftpobj);
    hs            = struct(ftpobj);
    hs.jobject.setConnectTimeout(300);  
    
    % Download the spectra.
    cd(ftpobj,PHOENIXfolder);
    mget(ftpobj,TemplName,Dir);
    sp = fitsread(fullfile(Dir,TemplName));
    
    % If needed download the wavelength vector
    if DownloadWV
        cd(ftpobj,'/v2.0/HiResFITS/');
        mget(ftpobj,'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits',Dir);
        wv = fitsread(fullfile(Dir,'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits'));
    end
    
    % If saving is not required - remove the downloaded files.
    if ~SaveFlag
        delete(fullfile(Dir,TemplName));
        if DownloadWV
             delete(fullfile(Dir,'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits'));
        end
    end
      
    close(ftpobj);
end

    
% Cut bounds:
I             = wv > min(TemplPars.Bounds) & wv < max(TemplPars.Bounds);
wv            = wv(I);
sp            = sp(I);

% Saving into structure format:
template.wv   = wv;
template.sp   = sp;
template.name = ['PHOENIX-T' num2str(TemplPars.Teff,'%5d') '-logg'  num2str(TemplPars.logg,'%3.2f') '-Z' num2str(TemplPars.z-1e-3,'%+2.1f')];
template.vel  = 0;
template.bcv  = 0;
 

end