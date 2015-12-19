function [pp ss slargs] = ppdefs(varargin)
% PPDEFS make some default plot structs for spikeslpos
  
  a = struct(...
      'datadir', '/media/disk/data/',...
      'rat', [],...
      'day', [],...
      'presets', [],...
      'e_cl', [],...
      'cdat',[]);
  
  a = parseArgsLite(varargin,a);
  
  if isempty(a.rat) && isempty(a.day)
    warning('no rat/day given, using current dir for all data');
    ratdir = pwd;
  else
    ratdir = [a.datadir '/' a.rat '/' a.day '/'];
  end

  pp = {};
  ss = {};
  slargs = spikeslpos('getargs');
  
  if isempty(a.presets)
    warning('no presets requested, returning');
    return
  end
  
  if ~iscell(a.presets)
    a.presets = {a.presets};
  end
  
  % load in e_cl, get defaults
  load ([ratdir '/processed/latest/e_cl.mat']);
  pd = plotdefs('e', e_cl);
  assignin('base', 'e_cl', e_cl);
  slargs.ename = 'e_cl';
  
  % load in usecls, segs
  load ([ratdir '/processed/latest/usecls.mat']);
  assignin('base', 'usecls', usecls);

  try
    load ([ratdir '/analysis/latest/segs.mat']);
    assignin('base', 'segs', segs);
  end
  
  % loop over all requested presets, adding them to the pp/segs list
  for j = 1:numel(a.presets);
    
    ppnew = {};
    segsnew = {};
    
    switch a.presets{j}
     
     case [],
      % do nothing

     case 'default';
      ppnew = pd.default;

      %% plot each cdat_ca1
     case 'cdat_ca1'
      if ~exist('cdat_eeg_ca1', 'var')
        load ([ratdir '/processed/latest/cdat_eeg_ca1.mat']);
      end
      
      for k = 1:size(cdat_eeg_ca1.data,2);
        ppnew{k} = pd.lfp;
        ppnew{k}.dat.contvar = 'cdat_eeg_ca1';
        ppnew{k}.dat.chans = k;
        ppnew{k}.height = 0.5;
      end      
      
      assignin('base', 'cdat_eeg_ca1', cdat_eeg_ca1);
     
      %% plot cdat_ripptet
     case 'cdat_ripptet'
      if ~exist('cdat_eeg_ripptet', 'var')
        load ([ratdir '/processed/latest/cdat_eeg_ripptet.mat']);
      end
      
      for k = 1:size(cdat_eeg_ripptet.data,2);
        ppnew{k} = pd.lfp;
        ppnew{k}.dat.contvar = 'cdat_eeg_ripptet';
        ppnew{k}.dat.chans = k;
        ppnew{k}.height = 0.5;
      end      
      
      assignin('base', 'cdat_eeg_ripptet', cdat_eeg_ripptet);

      
      %% plot each cdat_all
     case 'cdat_all'
      if ~exist('cdat_eeg_all', 'var')
        load ([ratdir '/cont/cdat_eeg_all.mat']);
      end
      
      for k = 1:size(cdat_eeg_all.data,2);
        ppnew{k} = pd.lfp;
        ppnew{k}.dat.contvar = 'cdat_eeg_all';
        ppnew{k}.dat.chans = k;
        ppnew{k}.height = 0.5;
      end      
      
      assignin('base', 'cdat_eeg_all', cdat_eeg_all);

     case 'ripp_t'
      if ~exist('ripples', 'var')
        load ([ratdir '/analysis/latest/ripples.mat']);
      end
      segsnew = ripples.ripp_times(end);

     case 'estcls'
      slargs.clnos = usecls.est;

     case 'segscand'
      segsnew = {segs.cand_nosleep};
      
     case 'rippenv'
      if ~exist('ripples', 'var')
        load ([ratdir '/analysis/latest/ripples.mat']);
      end
      cdat_renv = contcombine(ripples.cdat_rippenv, ripples.cdat_rippenv_smooth);
      ppnew{1} = pd.lfp;
      ppnew{1}.dat.contdata = cdat_renv;
      ppnew{1}.contdrawopt.levels = ripples.thresh;
      
     otherwise
      error(['unrecognized ''preset'' : ' a.presets(j)]);
      
    end
    
    % assertion test
    if ~iscell(ppnew) || ~iscell(segsnew)
      error('each ppnew or segsnew must be a cell array of slplots');
    end
    
    pp = [pp ppnew];
    ss = [ss segsnew];
  end
  
  slargs.plots = pp;
  slargs.seglists = ss;