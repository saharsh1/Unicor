function [jd,vels,evels,CCF] = Unicor2execution(Dir,Fname)
%[jd,vels,evels,CCF] = Unicor2execution(Dir,Fname)
% Initialize:
% ==========
% First, load the configuration of Unicor2.0 from the ini file. A file
% named "UnicorConfig.ini" should be saved in the data directory. If no
% UnicorConfig fille is found, the defualt initialization is used.
%
% NOTE: See the details of the parameters listed in the initialization file, or
% at the relevant routine that uses them. For simplicity, the parameters in
% are loaded in this script without additional data supplied.
%

if nargin == 0
[Fname,Dir] = uigetfile({'*.mat', '*.mat'; '*.*','All Files (*.*)'},'Pick a file');
end

warning off

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

ini = IniConfig();
if ~isempty(dir(fullfile(Dir,'UnicorConfig.ini')));
    ini.ReadFile(fullfile(Dir,'UnicorConfig.ini'));
else
    ini.ReadFile(fullfile('UnicorConfigDefault.ini'));
end


[mFilePath,~,~] = fileparts(mfilename('fullpath'));
TemplateDir     = fullfile(mFilePath,'templates');


% 1) Pre Proccess the observed spectra:
% =====================================
%
% Load Unicor2.0 'obs' structure.
% ------------------------------ 
% Load the obs array.
% See the details of the "obs" cell array contents a "calc_multi_order_velocity.m" header.
if exist(fullfile(Dir,[Fname(1:end-4) '_PostProc.mat']),'file')
    PreProcCheck = input('Processed obs file found. Repeat PreProc? [y/n] >> ','s');
else
    PreProcCheck = 'y';
end

if PreProcCheck == 'y'
    load(fullfile(Dir,Fname)); 

% Set the Pre-proccessing parameters 
% -------------------------------------
    [keys, ~  ] = ini.GetKeys('FiltPar');
    values      = ini.GetValues('FiltPar', keys);
    FiltPar     = cell2struct(values,keys);
  
    [keys, ~  ] = ini.GetKeys('TukeyWin');
    values      = ini.GetValues('TukeyWin', keys);
    BellPar     = cell2mat(values);
    
    [keys, ~  ] = ini.GetKeys('ExcludeOutLiers');
    values      = ini.GetValues('ExcludeOutLiers', keys);
    ExcludeOL   = cell2struct(values,keys);
   
% Run Unicor2.0 Pre-processing routine, and save the data.
    obs         = SpecPreProcess(obs,FiltPar,BellPar,ExcludeOL);  
    save(fullfile(Dir,[Fname(1:end-4) '_PostProc.mat']),'obs','FiltPar','BellPar');
else
    load(fullfile(Dir,[Fname(1:end-4) '_PostProc.mat']));
end
% 2) Pre Proccess the theoretical template:
% ========================================
%

if exist(fullfile(Dir,[Fname(1:end-4) '_TempPostProc.mat']),'file')
    PreProcCheck = input('Processed template file found. Repeat PreProc? [y/n] >> ','s');
else
    PreProcCheck = 'y';
end

if PreProcCheck == 'y'
% Get the template parameters
% --------------------------- 
    [keys, ~  ] = ini.GetKeys('Template');
    values      = ini.GetValues('Template', keys);
    TemplPars   = cell2struct(values,keys);
  
% Set the Pre-proccessing parameters 
% ----------------------------------- 
    [keys, ~  ]  = ini.GetKeys('FiltPar Template');
    values       = ini.GetValues('FiltPar Template', keys);
    FiltPar      = cell2struct(values,keys);
  
    [keys, ~  ]  = ini.GetKeys('TukeyWin Template');
    values       = ini.GetValues('TukeyWin Template', keys);
    BellPar      = cell2mat(values);

    [keys, ~  ]  = ini.GetKeys('Bounds Template');
    values       = ini.GetValues('Bounds Template', keys);
    DeltaV       = cell2mat(values);

% Retrive theoretical template and preprocess
% -------------------------------------------
    template     = GetTemplate(TemplateDir,TemplPars,true);                         % Template devision into orders + Template binning (if needed).
% Convert to air wavelengths, if required:
    if isfield(TemplPars,'VaccumWV')
      if ~TemplPars.VaccumWV 
        % The following transformation is taken from Morton (2000, ApJ. Suppl., 130, 403)
        s2             = 10^8 .* template.wv.^(-2);
        f              = 1 + 0.0000834254 + 0.02406147 ./ (130 - s2) + 0.00015998 ./ (38.9 - s2);    
        template.wv     = template.wv./f;
      end  
    end
% Rotational broadening, if required:    
    if isfield(TemplPars,'vsini')
    template     = vsini_matrot(template,TemplPars.vsini,TemplPars.epsilon);
    end
% Gaussian broadening:
    template     = inst_broad(template,TemplPars.R);
    template     = TempPreProcess(template,obs{1}.wv,DeltaV,FiltPar,BellPar);       % Prepares PHOENIX temp for cross-correlation analysis using 'SpecPreProc_main' of Unicor2.

    save(fullfile(Dir,[Fname(1:end-4) '_TempPostProc.mat']),'template','FiltPar','BellPar');
else
     load(fullfile(Dir,[Fname(1:end-4) '_TempPostProc.mat']));
end

% 3) CCF calculation
% ==================
% This part recieves multi-order spectra and a theoretical multi-order template,
% and returns velocity and cross-correlation function per order, for each observation.
if exist(fullfile(Dir,[Fname(1:end-4) '_CCF.mat']),'file')
    PreProcCheck = input('CCF analysis found. Repeat PreProc? [y/n] >> ','s');
else
    PreProcCheck = 'y';
end

if PreProcCheck == 'y'
    [keys, ~  ]       = ini.GetKeys('Unicor CCF pars');
    values            = ini.GetValues('Unicor CCF pars', keys);
    CCFpars           = cell2struct(values,keys);
    
    [v,CCF,CCFpeaks,CCFRMS] = calc_multi_order_velocity(obs,template,CCFpars); % CCF and velocity calculation
    save(fullfile(Dir,[Fname(1:end-4) '_CCF.mat']),'v','CCF','CCFpars','CCFpeaks','CCFRMS');
else
    load(fullfile(Dir,[Fname(1:end-4) '_CCF.mat']));
end

%RV extraction
% ------------
% Sets the parameters for the RV extraction:
[keys, ~  ]  = ini.GetKeys('RV matrix pars');
values       = ini.GetValues('RV matrix pars', keys);
RVpars       = cell2struct(values,keys);

% Determine what is the required initial set of  orders to be used.
CCFpk_vec    = mean(CCFpeaks);
CCFRMS_vec   = mean(CCFRMS);
I0 = CCFpk_vec > RVpars.CCFminPeak & CCFpk_vec > RVpars.CCFpk2RMS*CCFRMS_vec;
[vels, evels] = MatrixRV(v,RVpars.IterN,RVpars.C,RVpars.SigThr,I0); 


for i = 1:length(obs)
    jd(i) = obs{i}.bjd;
end

warning on
end
