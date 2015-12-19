function [sla sld slcache arglist] = spikeslpos (varargin)
% SPIKESLPOS: plot place-centered rasters or position estimate overlaid
% on linearized animal trajectory, on a bed of steamed chard.
%
%  [sla sld arglist] = ... 
%         spikeslpos('ename','paul05run','runspeed', 15);
%
% After scrolling around, to, e.g. exclude clusters from the existing plot:
%
%  sla = spikeslpos('usefig', (gcf), 'excludecls', [20 29:33]); 
%
%  the sla struct can be saved and used to recreate the identical
%  plot, using the 'argstruct' option.
%
%  the 'arglist' output struct is useful for getting the possible args by
%  tab-completion at the Matlab command line--use 'getargs'
%
% ARGUMENTS:
%
% * = required
% [ = not implemented yet
% F = flag argument, true if provided
%
%
%   Only one of 'argstruct'/'usefig'/'usepanel' can be given, and must
%   appear first in arglist, subsequent args on command line will
%   override. Useful for changing existing view from command line.
%
%     'argstruct': a struct containing fields corresponding to the named
%                  arguments to spikeslpos. ([])
%        'usefig': use argstruct from specified figure, and plot into it ([])
%      'usepanel': use argstruct from specified uipanel, and plot into it ([])
%           'fig': draw plots into specified figure
%         'panel': draw plots into specified uipanel
%
%       'getargs': just return an empty argstruct. no value, must be
%                  first arg.
%
%         'plots': plots to display. See function slmkplots. Default to:
%                  sp.default as defined in plotdefs.m
%
%         'cache': pre-existing object cache to use ([] = use from figure
%                  appdata)
%    'flushcache': empty the figure's cache, recompute all objs from
%                  scratch once (false)
%    'cachelimit': max # of objects of each type to keep in the cache 
%                  (default 20)
%
% display args:
%       'whitebg': use white or black background for panel and all plots
% 'whitebg_plots': set background color for all plots
%   'whitebg_panel': set background color for panel
%
% general args:
% *       'ename': name of experimental struct in workspace as string
%
%         'clnos': clusters to plot - vector of indexes or logical array
%                  into e.cl (all cls)
%    'excludecls': clusters to exclude - vector/logical array into e.cl ([])
%      'tcthresh': further exclude clusters without at least this peak (or other
%                  function set by 'tcthreshfn') firing rate in 'place' 
%                  pfmap during 'run'.If [mx2] vector specified, only use 
%                  clusters with mean f.r. between the two values. 
%                  (0 = use all clusters)
%    'tcthreshfn': function handle of function to apply to place-rate
%                  maps for thresholding by 'tcthresh'.Try @mean or
%                  @spatialinfo. (@max = peak rate)
%      'randperm': randomly permute cl -> tc mapping on each redraw (false)
%        'clperm': use provided permutation of cl -> tc mapping (can't
%                  also use 'tcthresh' or 'excludecls') ([])
%
%       'timewin': [1x2] time to display (data start/end)
%      'viewlist': [mx2] array of start/end times. Only used when no
%                  'timewin' arg is provided ([])
%     'viewlisti': current index into array of viewlist (1)
%   'viewlistwin': width of window centered on viewlist ([] = width of
%                  viewlist entry)
% [  'viewlabels': cell array of names of views in viewlist ({})
%
%      'runspeed': runspeed filter for place field finding (from
%                  tcoptdefs: 15 cm/s)
%      'runtimes': [runstart runend] times to use for run analysis
%                  ([] = use e.runtimes)
%
% parameter reconstruction/estimation (parmest) args:
%           'pixwid': width of time bins, in pixels (4)
%      'timebinsize': size of timebin ([] = use pixwid)
%       'parmestopt': options for position reconstruction; _see mkparmestopt_
%   'parmestdistbin': size in cm of binsize to use for parameter estimate (10)
%
% position estimate display args
%    'parmestgwtime': filter to use in time (1 = no filtering)
%    'parmestgwdist': filter to use in distance (1 = no filtering)
%                     for, e.g. 60cm filter: gausswin(60/parmestdistbin,2.5)
%    'parmestsharpen': factor to raise parmest map to. ([] = no sharpening)
%      'parmestscale': factor to multiply parmest map by. ([] = no scaling)
%  'parmestautoscale': scale maps (same factor across fwd/rev/place) such
%                     that a fractional value (across all 3 maps) is
%                     max. Range 0-1, try 0.995 ([] = no scale)
%
% raster args:
%        'maxpks': plot raster at each of 'top' # of of peaks (1)
%      'rasterbg': plot light gray background in raster when no spikes (true)
% 'rastcolormode': how to choose colors for rasters: ('electrode') same
%                  color for all clusters on each tetrode; 'clno' cycle
%                  through colors by cluster number; 'firstpeak' cycle 
%                  through colors by order of first peak; 'parm',
%                  use one color for all rasts, determined from
%                  parmdrawopt ptcolor.
%   'peakfindbin': binsize to use for TC peak finding in cm (2cm)
%    'peakfindgw': size of gausswin to use for filtering TC peak finding (40cm)
%    'peakfindth': firing rate threshold to use for TC peak finding (2Hz)
%
% continuous data:
% 'contlimsequal': [boolean] use min/max of all contlims for each contplot
%
% overlays:
%      'seglists': {[mx2] or [mx1]} cell array of p segment lists, each containing
%                  m segments, specified as [start end] or [time].  ({})
%     'segcolors': [px3] array of colorspecs (prism(p)*0.4)
%      'plotsegs': [p] logical vector of whether to plot seglists (true(1,p))
%
%

% $Id$


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% constants

% approximate frequency of position record samples
posrecHz = 30;

%% plot layout stuff

% string displayed in title bar
titlebarhelp = ['SPIKESLPOS  { } [ ] = pan;  , . = zoom; b = rasterbg;' ...
  ' + - = # pks per cell; p = position; f = raster/parmest; h = crosshair'...
  ];

% positioning of plot within panel
lmarginpix = 40; % room for labels
rmarginpix = 3; % so it doesn't look cut off
botmarginpix = 28; % room for labels
topmarginpix = 20; % room for status text

% vertical gap (in pixels) between subplots
subaxgappix = 5;

%% Colormap stuff

%%% whitebg/blackbg

% color of panel bg
panelbg_white = [1 1 1];
panelbg_black = [0.2 0.2 0.2];

% colors of axes that overlap with the panel, and of plotboxes,
% regardless of each plot's whitebg status
axcol_white = [0.3 0.3 0.3];
axcol_black = [0.8 0.8 0.8];

%%%rasters
% number of cluster colors for rasters
nclustercols = 32;

% color of rasters when no spikes present
% would be nice to move this to mkraster, someday
rasterplotbg_black = [0 0 0]; 
rasterbgcolor_black = [0.15 0.15 0.15];
rasterplotbg_white = 1-rasterplotbg_black;
rasterbgcolor_white = 1-rasterbgcolor_black;


%% build colormap
cmap = [];

% grayscale
cmapgray = gray(64);

cmapgraystart = size(cmap,1)+1;
cmap = [cmap ; cmapgray];
cmapgrayend = size(cmap,1);


% inverse gray (black = high value)
cmapgrayinv = 1-gray(64);

cmapgrayinvstart = size(cmap,1)+1;
cmap = [cmap ; cmapgrayinv];
cmapgrayinvend = size(cmap,1);


% hot map for position estimate;
%cmaphot = brighten(jet(64),-0.3);
cmaphot = hot(64);

cmaphotstart = size(cmap,1)+1;
cmap = [cmap ; cmaphot];
cmaphotend = size(cmap,1);


% jet map for position estimate;
%cmaphot = brighten(jet(64),-0.3);
cmapjet = jet(256);

cmapjetstart = size(cmap,1)+1;
cmap = [cmap ; cmapjet];
cmapjetend = size(cmap,1);


% colored spike rasters
%
% raster is very finicky, since it uses indexes into the cmap to draw all
% those markers pixel by pixel. The first 4 entries in its map must be as
% below, followed by an hsv colormap.
% first four colors are reserved: black, white, 2 raster bgs (lighter and
% darker gray)
cmaprast = [rasterplotbg_black; rasterbgcolor_black; ...
            rasterplotbg_white; rasterbgcolor_white; ...
            hsv(nclustercols)];

cmapraststart = size(cmap,1)+1;
cmap = [cmap ; cmaprast];
cmaprastend = size(cmap,1);
nrastcolors = cmaprastend - cmapraststart + 1;


% get clims for each map, *assuming data over interval 0-1*
% newclim(BeginSlot,EndSlot,CDmin,CDmax,CmLength)

cmapl = size(cmap,1);

clim01gray = newclim(cmapgraystart,cmapgrayend,0,1,cmapl); %#ok
clim01grayinv = newclim(cmapgrayinvstart,cmapgrayinvend,0,1,cmapl); %#ok
clim01hot = newclim(cmaphotstart,cmaphotend,0,1,cmapl); %#ok
clim01jet = newclim(cmapjetstart,cmapjetend,0,1,cmapl); %#ok
climrast = newclim(cmapraststart, cmaprastend,1,nrastcolors,cmapl);

cmapwingray = [cmapgraystart cmapgrayend];
cmapwingrayinv = [cmapgrayinvstart cmapgrayinvend];
cmapwinhot = [cmaphotstart cmaphotend];
cmapwinjet = [cmapjetstart cmapjetend];
cmapwinrast = [cmapraststart cmaprastend];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% parse args:

% note defaults are set in code below, as many depend on other params, so
% we can't just use parseArgsLite

a = struct( ...
  'argstruct', [], ...
  'usefig', [], ...
  'fig', [],...
  'usepanel', [],...
  'panel', [],...
  'plots', [], ...
  'whitebg', [],...
  'whitebg_plots', [],...
  'whitebg_panel', [],...
  'cache', [],...
  'flushcache', false,...
  'cachelimit', 20,...
  'ename', [], ...
  'clnos', [], ...
  'excludecls', [], ...
  'randperm', false,...
  'clperm', [],...
  'timewin', [], ...
  'pixwid', [], ...
  'timebinsize', [], ...
  'viewlist', [], ...
  'viewlisti', [], ...
  'viewlistwin', [], ...
    ...
  'runspeed', [], ...
  'runtimes', [], ...
  'tcthresh', [], ...
  'tcthreshfn', [], ...
      ...
  'parmestdistbin', [], ...
  'parmestopt', [],...
      ...
  'parmestautoscale',[],...
  'parmestsharpen', [], ...
  'parmestscale', [], ...
  'parmestgwtime', [], ...
  'parmestgwdist', [], ...
      ...
  'maxpks', [], ...
  'rasterbg', [], ...
  'rastcolormode',[],...
  'peakfindbin', [], ...
  'peakfindgw', [], ...
  'peakfindth', [], ...
      ...
  'seglists', {{}}, ...
  'segcolors', [], ...
  'plotsegs', [], ...
      ...
  'contlimsequal', false...
      );
  
  % useful in base workspace for tab completion of argument names
  arglist = cell2struct(fieldnames(a),fieldnames(a));
  
  if isempty(varargin),
    error('no args provided');
  end
  
  % pass out an empty argstruct, if requested, on spikeslpos('getargs');
  if strcmp(varargin{1},'getargs')
    sla = a;
    sld = [];
    slcache = [];
    %arglist also passed out
    return;
  end


a = parseArgsLite(varargin, a);

% save a copy of args as passed in on command-line for debugging
a_in = a;

%% _argstruct_ / _usefig_ / _usepanel_

switch sum([~isempty(a.argstruct) ...
            ~isempty(a.usefig) ...
            ~isempty(a.usepanel)])
 case 0,
  % OK, just use args as provided
 case {2 3} 
  error([mfilename ':BadArgs'], ...
        'At most one of ''usefig'', ''usepanel'' and ''argstruct'' can be provided');
 case 1,
  
  for usestr = {'usefig' 'usepanel' 'argstruct'}
    usestr = usestr{1};
    matchi = strcmpi(varargin, usestr);
    if any(matchi)
      if matchi(1) ~= 1, 
        error([mfilename ':BadArgs'], ...
              ['''' usestr ''' must be first argument, if provided']);
      end
      if sum(matchi) > 1, 
        error([mfilename ':BadArgs'], ['Only one ''' usestr ''' may be provided']);
      end
    end
  end

  %% _usefig_ / _usepanel_
  if ~isempty(a.usefig)
    if ~ishandle(a.usefig) || ~strcmp(get(a.usefig,'type'), 'figure'),
      error([mfilename ':BadArgs'], '''usefig'' must be figure handle'); 
    end
    
    figp = a.usefig;
    % find uipanel, fullsize it
    panel = findobj(a.usefig,'type','uipanel');

    if length(panel) ~= 1;
      error('''usefig'' must specify a figure with exactly one panel, try ''usepanel''');
    end

    a.usefig = [];
    usepanel = panel;
  end
  
  if ~isempty(a.usepanel)
    if ~ishandle(a.usepanel) || ~strcmp(get(a.usepanel,'type'), 'uipanel'),
      error([mfilename ':BadArgs'], '''usepanel'' must be uipanel handle'); 
    end

    figp = ancestor(a.fig, 'figure');
    usepanel = a.usepanel;
    
  end

  if exist('usepanel', 'var') && ~isempty(usepanel)
    
    % use slargs from uipanel specified by 'usefig' or 'usepanel'
    a_tmp = getappdata(usepanel,'slargs');
    if isempty(a_tmp),
      error('no spikeslpos uipanel in figure specified by ''usepanel''');
    else      
      a = a_tmp;
    end
    
  end

  %% _argstruct_

  if ~isempty(a.argstruct),
    
    % don't allow for recursive argstructs
    if ~isempty(a.argstruct.argstruct);
      warning([mfilename ':BadArgs'],'Recursive argstructs not allowed--ignoring');
      a.argstruct = [];
    end
    
    % use argstruct as defaults, let additional params override.
    a = a.argstruct;
  end

  % use usefig, usepanel or argstruct to get new defaults (a) then parse
  % the rest of the args
  a = parseArgsLite(varargin(3:end),a);

end % end usefig/usepanel/argstruct, now we have new 'a'

%% _fig_ / _panel_ to draw into
if ~isempty(a.fig)
  figp = figure(a.fig);
  clf(figp);
  panel = [];
end

if ~isempty(a.panel)
  if ~ishandle(a.panel) || ~strcmp(get(a.panel, 'type'), 'uipanel')
    error('''panel'' must be a handle to a uipanel');
  end
  figp = ancestor(a.panel, 'figure');
  panel = a.panel;
end

% Make a new figure/panel, if none specified
if ~exist('panel', 'var') || isempty(panel)
  figp = gcf;
  panel = findobj(figp,'type', 'uipanel');
  if ~isempty(panel),
    delete(panel(2:end));
  else
    panel = uipanel;
  end
  set(panel(1), 'position', [0 0 1 1],...
                'parent', figp,...
                'bordertype', 'none');
end

% $$$   
% $$$   panel = uipanel('parent', figp,...
% $$$                   'position', [0 0 1 1],...
% $$$                   'bordertype', 'none');


% flag the panel as being worked on ('locked') so that we can ignore
% resizes. This is necessary b/c matlab ignores the 'Interruptible' flag
% for resizes, and triggers a resize when we change the units of the
% panel
setappdata(panel,'locked', true);

% fetch appdata right after identifying panel to use
sld = getappdata(panel,'sldata'); % getappdata returns [] if no appdata yet



% _cache_/_flushcache_
if a.flushcache
  slcache = [];
  a.flushcache = false; % only do it once
elseif ~isempty(a.cache),
  slcache = a.cache;
else
  slcache = getappdata(figp,'slcache');
end

% don't keep the cache in the argstruct
a.cache = [];


% OK, we have now parsed the various argstructs, and figured out which figure
% and panel to plot into

% Display args:

% _whitebg_/_whitebg_panel_/_white_bg_plots_

if ~isempty(a.whitebg)
  a.whitebg_panel = a.whitebg;
  a.whitebg_plots = a.whitebg;
end

if isempty(a.whitebg_panel)
  if ~isempty(a.whitebg_plots),
    % use the requested plots arg for panel, if no other arg provided
    a.whitebg_panel = a.whitebg_plots;
  else
    warning(['Panel bg color not specified, using whitebg: use either ' ...
             '''whitebg'' or ''whitebg_panel'' or ''whitebg_plots''']);
    a.whitebg_panel = 1;
  end
end

% Generic args:

% _ename_
% use just the name of e in appdata (below) to save memory when we store the params
% in the figure's appdata or for saved argstructs

if ischar(a.ename) && evalin('base', ['exist(''' a.ename ''', ''var'')']),
  e = evalin('base',a.ename);
end

if ~exist('e','var') || ~isstruct(e),
  eargnum = find(strcmpi('ename', varargin), 1, 'first')+1;
  if isempty(eargnum) || isempty(inputname(eargnum)),
    error([mfilename ':BadArgs'], '''ename'' must be name of valid expt (as a string)');
  else
    a.ename = inputname(eargnum);
    e = evalin('base',a.ename);
  end
end

% _clnos_ / _excludecls_
% clusters to use/exclude from all analyses

clnos = a.clnos;
if islogical(clnos),
  clnos = find(clnos);
end

excludecls = a.excludecls;
if islogical(excludecls),
  excludecls = find(excludecls);
end

% use all clusters by default
if isempty(clnos),
  clnos = (1:length(e.cl));
end

clnos = setdiff(clnos, excludecls);

% _randperm_ / _clperm_
% shuffle cl -> tc mapping

if a.randperm && ~isempty(a.clperm),
  error([mfilename ':BadArgs'], 'Only one of ''randperm'' and ''clperm'' can be provided.');
end

if ~isempty(a.clperm) && (~isempty(a.tcthresh) || ~isempty(a.excludecls))
  error([mfilename ':BadArgs'], '''clperm'' and ''tcthresh''/''excludecls'' not currently supported together');
end

% defaults to empty if no shuffle
% randperm will set clperm later if requested
if a.randperm,
  if islogical(clnos),
    clperm = randperm(sum(clnos));
  else
    clperm = randperm(length(clnos));
  end
else
  clperm = a.clperm;
end

% _timewin_ / _viewlist_ / _viewlisti_ /_viewlistwin_
% note these are requested times, may be cropped/bumped below

% timewin overrides viewlist (so that the user can move around after
% loading a view, but keep their place in the viewlist), so don't provide
% one if you want to use viewlist.

if ~isempty(a.viewlist)
  % start at beginning of viewlist
  if isempty(a.viewlisti) || a.viewlisti > size(a.viewlist,1),
    a.viewlisti = 1;
  end
end

if ~isempty(a.timewin),  % use timewin if provided
  tstart = a.timewin(1);
  tend = a.timewin(2);
else
  
  if ~isempty(a.viewlist) % use viewlist if no timewin
    switch size(a.viewlist,2),
     case 1,
      % single-column viewlist specifies centers
      if isempty(a.viewlistwin),
        warning(['for single-column viewlists, please provide a viewlistwin, ' ...
                 ' using 1 second window']);
        viewlistwin = 1;
      else
        viewlistwin = a.viewlistwin;
      end
     case 2,
      viewlistwin = a.viewlistwin;
     otherwise,
      error('''viewlist'' must be a 1- or 2-column array');
    end
    % handles 1 and 2 column case
    tstart = min(a.viewlist(a.viewlisti,:));
    tend = max(a.viewlist(a.viewlisti,:));
    if ~isempty(viewlistwin)
      % use fixed window size, if provided
      vlmid = (tstart + tend) / 2;
      tstart = vlmid - viewlistwin/2;
      tend = vlmid + viewlistwin/2;
    end
    
  else % use sldata if no timewin or viewlist
     if isstruct(sld),
       % changed from x to t, why was it x? this is requested time,
       % should stay the same. b/c we didn't used to have a dstart!
       tstart = sld.tstart;
       tend = sld.tend;
     
     else % user has given us no idea--plot around zero
      tstart = [];
      tend = [];
     end
  end
end


% _pixwid_ / _timebinsize_

% default minimum width of timebin for parmest in pixels
if isempty(a.pixwid);
  a.pixwid = 4;
end

% $$$ %_tcthresh_ / _tcthreshfn_
% $$$ 
% $$$ if ~isempty(a.tcthreshfn) 
% $$$   if ~isa(a.tcthreshfn, 'function_handle'),
% $$$     error('''tcthreshfn'' must be a function handle (use @-notation)');
% $$$   end
% $$$   tcthreshfn = a.tcthreshfn;
% $$$ else
% $$$   tcthreshfn = func2str(@max);
% $$$ end
  
% _plots_
if isempty(a.plots)
  sp = plotdefs('e', e);
  a.plots = sp.default;
end

if ~iscell(a.plots),
  a.plots = {a.plots};
end

% _runspeed_
% if provided, will override all runspeeds in 'plots' loop

% _runtimes_
% default is e.runtimes
if isempty(a.runtimes)
  runtimes = e.runtimes;
else
  runtimes = a.runtimes;
end

% Parameter estimation (parmest) prefs:

% _parmestdistbin_ 
if isempty(a.parmestdistbin),
  a.parmestdistbin = 10; % 10 cm bins
end

% _parmestgwtime_
if isempty(a.parmestgwtime),
  a.parmestgwtime = gausswin(1); % no filter 
end

% _parmestgwdist_
if isempty(a.parmestgwdist),
  a.parmestgwdist = gausswin(1); % no filter
end

% _parmestopt_
if isempty(a.parmestopt),
  a.parmestopt = mkparmestopt('preset','bayes_zhang');
end

% _parmestsharpen_
% defaults to [] 

% _parmestscale_
% defaults to [] 

% _parmestautoscale_
if ~isempty(a.parmestautoscale) && ...
      (a.parmestautoscale <= 0 || ...
       a.parmestautoscale > 1),
  error('''parmestautoscale'' must be >0 and <=1, use [] to turn off autoscaling');
end

% defaults to [] 


% Raster prefs:

% _maxpks_
if isempty(a.maxpks) || a.maxpks < 0;
  a.maxpks = 1;
end

% _rasterbg_
if isempty(a.rasterbg),
  a.rasterbg = false;
end

% _rastcolormode_
if isempty(a.rastcolormode),
  a.rastcolormode = 'electrode';
end

% _peakfindbin_ (cm)
if isempty(a.peakfindbin)
  a.peakfindbin = 2;
end

% _peakfindgw_ (cm)
if isempty(a.peakfindgw)
  a.peakfindgw = 40;
end

% _peakfindth_ (Hz)
if isempty(a.peakfindth)
  a.peakfindth = 2;
end

% Overlay args:

% Segments = lists of time windows, to plotted as translucent colored boxes
% across all plots

% _seglists_
if ~iscell(a.seglists)
    error([mfilename ':Badargs'], 'Argument ''seglists'' must be cell array');
end

nseglists = length(a.seglists);

for k = 1:nseglists,
  if ~isnumeric(a.seglists{k}) || size(a.seglists{k},2) > 2,
    error([mfilename ':Badargs'], 'Each seglist must be 1- or 2-column numeric array');
  end
end

%   _segcolors_
if ~isempty(a.segcolors),
  if size(a.segcolors,2) ~= 3 || ~(size(a.segcolors,1) >= nseglists),
    error([mfilename ':Badargs'], ['''segcolors'' must be 3-column array ' ...
                        'with at least as many rows as seglists']);
  end
  segcolors = a.segcolors;
else
  segcolors = hsv(nseglists) * 0.8;
% $$$   segcolors = hsv(nseglists) * 0.4;
end

%  _plotsegs_
if ~isempty(a.plotsegs),
  if ~islogical(a.plotsegs) || length(a.plotsegs(:)) ~= nseglists,
    error([mfilename ':Badargs'], '''plotsegs'' must be logical vector with as many elements as seglists');
  end
  plotsegs = a.plotsegs;
else
  plotsegs = true(1,nseglists);
end


% Continuous/specgram data plotting - now part of a.plots framework
% see matlab/cont/* , matlab/obj/cont/*


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Switch to 'wait' cursor

% Tell the user they got to wait while we do some stuff!

oldptr = (get(figp,'pointer'));
if strcmp(oldptr, 'watch')
  oldptr = 'arrow';
end
set(figp,'pointer','watch');
pause(0.01); % necessary for pointer update (why.)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set up figure/panel/axes
% (see end of m-file for post-plotting stuff)

%%% figure stuff

% zbuffer blinks, is slow, software OpenGL chokes on complex figures. Only
% 'painters' has fast drawmode (no z-depth calculations, so we just
% draw/redraw in right order)

% have to set these on the figure
set (figp, 'renderer', 'painters');
%set (figp, 'renderer', 'zbuffer');
%opengl software;
%set (figp, 'renderer', 'opengl');

set (figp, 'doublebuffer', 'on');
set (figp, 'backingstore', 'off');

%no menu bar
set (figp, 'menubar', 'none');
drawnow;

% helpful titlebar
%set(figp,'numbertitle','off');
set(figp,'name', [a.ename ' - ' titlebarhelp]);


% set up callback functions
set(figp,'keypressfcn', 'slkeypress(gcbf)');
set(figp,'windowbuttondownfcn', 'slmousedown(gcbf)');
set(figp,'windowbuttonmotionfcn', 'slmousemove(gcbf)');
set(figp,'windowbuttonupfcn', 'slmouseup(gcbf)');
set(figp,'interruptible','off');
set(figp,'resizefcn','slresize(gcbf)');

% assign colormap
set(figp,'Colormap', cmap);

%%% panel stuff


%black bg
if a.whitebg_panel,
  set (panel, 'BackgroundColor', panelbg_white);
else
  set (panel, 'BackgroundColor', panelbg_black);
end


%%% main/overlay axes 'ax'

if isstruct(sld)
  ax = sld.ax;
else
  ax = axes('parent', panel);
  set (ax, 'units', 'pixels');
  
  % axes on top of data
  set (ax, 'layer', 'top');
  
  % no depth ordering (draw in correct order)
  set (ax, 'drawmode', 'fast');

  % transparent, so we can see the subaxes
  set (ax, 'color', 'none');

  % this axis provides the time ticks
  set (ax, 'tickdir', 'out');
  set (ax, 'ticklength', [0.005 0.025]);

  % hide the y-axis (can't turn it off!)
  set (ax, 'ycolor', get(panel,'BackgroundColor'));
  set (ax, 'yaxislocation', 'right');
  set (ax, 'yticklabel', {});
  set (ax, 'ytick', []);
  
  % don't show a box around the plot
  set (ax, 'box', 'off');
end

% make time tick axis legible
if a.whitebg_panel
  set (ax, 'xcolor', axcol_white);
else
  set (ax, 'xcolor', axcol_black);
end

hold(ax, 'on');


% position 'ax' 
% setting panel units triggers a resize! Pain in the ass!
% $$$ if ~strcmp(get(panel,'units'), 'pixels');
% $$$   disp(['bad units = ' get(panel,'units')])
% $$$   error('not getting set correctly');
% $$$ end
  
set(panel, 'units', 'pixels');
panelpos = get(panel, 'position');

% have to set it back so that it gets sized with figure resizes
set(panel, 'units', 'normalized');

axlpix = lmarginpix;
axbotpix = botmarginpix;
axwidpix = panelpos(3)-lmarginpix-rmarginpix;
axhtpix = panelpos(4)-botmarginpix-topmarginpix;

% set position on init/panel resize
if ~isstruct(sld) || ~all(sld.panelpos == panelpos),
  set(ax, 'position', [axlpix axbotpix axwidpix axhtpix]);
end

%%% subplots for data

nplots = length(a.plots);

% set plots' whitebg, if requested
if ~isempty(a.whitebg_plots),
  for k = 1:nplots
    a.plots{k}.whitebg = a.whitebg_plots;
  end
end

% initialize subaxes:

if isstruct(sld) 
  if length(sld.subaxes) == nplots,
    % clear and reuse subaxes
    subaxes = sld.subaxes;
    for k = 1:nplots,
      cla(subaxes(k));
      set(subaxes(k),'ytickmode','auto');
      set(subaxes(k),'yticklabelmode','auto');
      legend(subaxes(k),'off');
    end
  else
    try
      delete (sld.subaxes);
    catch
      % no biggie if they're already gone
      warning('Couldn''t delete subaxes');
    end
    sld.subaxes = [];
  end
end

if ~isstruct(sld) || isempty(sld.subaxes),
  % make new subaxes
  for k = 1:nplots;
    subax = axes('parent', panel);
    subaxes(k) = subax;
    set (subax, 'drawmode', 'fast');
    set (subax, 'units', 'pixels');

    % set up axis tickmarks
    set (subax, 'tickdir', 'out');
    set (subax, 'ticklength', [0.005 0.025]);   
  end
end

% do this setup every time, in case plot or panel bg color has changed
for k = 1:nplots
  
  subax = subaxes(k);
  hold(subax, 'on');
  
  % hide x-axis (can't turn it off!)
  set (subax, 'xcolor', get(panel,'BackgroundColor'));
  set (subax, 'xticklabel', {});
  set (subax, 'xtick', []);
  

  if ~a.whitebg_panel, 
    % make text legible by setting color of plot axes relative to *panel* bg
    % color, irrespective of bg of each plot
    axcol = axcol_black;
  else
    axcol = axcol_white;
  end

  %show a box around the subax
  set (subax, 'box', 'on');
  
  % show axes
  set (subax, 'layer', 'top');

  % make text legible against panel bg
  set (subax, 'xcolor', axcol);
  set (subax, 'ycolor', axcol);

end


%%% size/resize subaxes

% start plotting from the top of ax
subaxtop = axbotpix + axhtpix;

% vertical pixels available
useableaxhtpix = axhtpix - ((nplots-1) * subaxgappix);

% sum of plots.height values
% can't use sum a.plots.height b/c plots is now a cell array
totalht = 0;
for k = 1:nplots,
  totalht = totalht + a.plots{k}.height;
end

for k = 1:nplots;

  subax = subaxes(k);

  % proportion of ax to take up
  subaxht = (a.plots{k}.height/totalht) * useableaxhtpix;
  subaxbot = subaxtop - subaxht;

  % position subax:
  set(subax, 'position', [axlpix subaxbot axwidpix subaxht]);

  % next subaxtop (add gap)
  subaxtop = subaxbot - subaxgappix;
end

%set(figp,'currentaxes',ax);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% start/end time handling

% tstart/tend (t*) are the requested start and end times (from timewin,
%                  viewlist, etc)
% dstart/dend (d*) are the limits of the data to be plotted
% xstart/xend (x*) are the x-limits of the plot
%
% (dstart_cl/dend_cl) are the limits of spike data from e.cl

% Logic flow: calculate d*, x* from requested t*
%
% Get data limits (min/max timestamp of spikes in all clusters)
%
% if t* empty, use data limits for d*, x*; if no data, x* = [0 1];
%
% if timewin doesn't overlap with data, revert to last t* from appdata, or
% error if not available
%
% if t* overlaps at all with data, plot requested window, but don't plot
% rasters/parmest in 'void' area (keep track of dstart/dend)


%%% Get data limits
  
maxts_cl = -Inf;
mints_cl = +Inf;

% If we're dealing with the same clusters from the same experiment,
% then reuse values in sldata
if isstruct(sld) && ...
      strcmp(sld.ename, a.ename) && ...
      all(size(sld.clnos) == size(clnos)) && ...
      all(sld.clnos == clnos),
  maxts_cl = sld.maxts_cl;
  mints_cl = sld.mints_cl;
  
else 
  for j = clnos
    cltsj = getepd(e.cl(j), 'time');
    if cltsj
      % timestamps are sorted, just use first/last one
      mints_cl = min([mints_cl cltsj(1)]);
      maxts_cl = max([maxts_cl cltsj(end)]);
    end
  end
end

maxts = maxts_cl;
mints = mints_cl;

% find extent of contdata, too (but may not have contdata yet)
for k = 1:length(a.plots);
  plotk = a.plots{k};
  if isfield(plotk.dat, 'contdata') && ~isempty(plotk.dat.contdata)
    mints = min([mints plotk.dat.contdata.tstart]);
    maxts = max([maxts plotk.dat.contdata.tend]);
  end
end

% find extent of position records
pos_timei = getepindex(e.pos, 'time');
if ~isempty(pos_timei)
    mints = min([mints e.pos.dat(1,pos_timei)]);
    maxts = max([maxts e.pos.dat(end,pos_timei)]);
end

% deal with case where clnos = [] or there are no spikes in any cluster
if mints_cl == +Inf, mints_cl = 0; end
if maxts_cl == -Inf, maxts_cl = 0; end

if (isempty(tstart) || isempty(tend)),
  
  %% no timewin requested, plot all data
  dstart = mints;
  dstart_cl = mints_cl;
  dend = maxts;
  dend_cl = maxts_cl;

  xstart = mints;
  xend = maxts;  
 
  if xstart == xend; 
    % 0 or 1 spikes, no requested timewin, plot can't have zero width
    xend = xend + 1;
  end
else
  
  %% requested tstart/tend out of data range, fix if possible
  if (tstart > maxts || tend < mints);
    % if we are scrolling from a valid view, use last requested tstart/tend
    if isstruct(sld) && isfield(sld,'tstart') && isfield(sld,'tend'),
      tstart  = sld.tstart;
      tend = sld.tend;
    else
      % otherwise complain
      error([mfilename ':TimeOutOfRange'], 'Start/End time out of range.');
    end
  end

  
  %% at least one of tstart/tend within data range

  % plot should be of requested width
  xstart = tstart;
  xend = tend;  
  
  % crop dstart to data range
  dstart = max([tstart mints]);
  dstart = min([dstart maxts]);
  
  dstart_cl = max([tstart mints_cl]);
  dstart_cl = min([dstart_cl maxts_cl]);

  % crop dend to data range
  dend = max([tend mints]);
  dend = min([dend maxts]);
  
  dend_cl = max([tend mints_cl]);
  dend_cl = min([dend_cl maxts_cl]);

end

if dend < dstart || dend_cl < dstart_cl,
  error([mfilename ':EndTimeBeforeStartTime'],...
        'Data End Time is before or equal to Data Start Time');
end

% get minimum timebin based on pixwid (had to defer until we know how big
% the axis will be

axwidt = xend - xstart;

% width in time of one pixel
pixwid1 = max(axwidt/axwidpix);
mintimebin = a.pixwid * pixwid1;

% $$$ % place firing threshold ('pfthresh')
% $$$ if ~isempty(a.pfthresh),
% $$$   cli = pfthreshfn(ppfdats') > a.pfthresh(1);
% $$$   if length(a.pfthresh(:)) > 1,
% $$$     cli = cli & (pfthreshfn(ppfdats') < a.pfthresh(2));
% $$$   end
% $$$ else
% $$$   cli = true(1,length(clnos));
% $$$ end
% $$$ 
% $$$ clnos = clnos(cli);


%%%%% distribute arguments from user/spikeslpos into plot structs

% set cache limit
slcache = mkcache('cache',slcache,'limit', a.cachelimit);

ptypenos = struct();

for k = 1:nplots;

  ptype = a.plots{k}.type;

  % tell each plot what its plotnumber is
  a.plots{k}.plotno = k;

  % keep a count of how many of this type of plot we've seen, and tell
  % the plot
  if ~isfield(ptypenos, ptype),
    ptypenos.(ptype) = 0;
  end
  ptypenos.(ptype) = ptypenos.(ptype) + 1;
  a.plots{k}.plottypeno = ptypenos.(ptype);
  
  % set timewin of all plots (data to use)
  if strcmp(ptype,{'raster' 'parmest'}),
    a.plots{k}.dat.timewin = [dstart_cl dend_cl];
  else
    a.plots{k}.dat.timewin = [dstart dend];
  end
  
  % set timewin (xlim of final plot)
  a.plots{k}.dat.timewin_plot = [xstart xend];


  % set exp struct name of all plots
  if isfield(a.plots{k}.dat, 'ename'),
    a.plots{k}.dat.e = [];
    if isempty(a.plots{k}.dat.ename),
      a.plots{k}.dat.ename = a.ename;
    end
  end
  
  % set clnos of all plots
  if isfield(a.plots{k}.dat, 'clnos') && ~isempty(a.clnos)
    a.plots{k}.dat.clnos = clnos;
    a.plots{k}.dat.clperm = clperm;
  end
  
  % set runtimes of all plots
  if isfield(a.plots{k}, 'tunecurveopt')
    a.plots{k}.tunecurveopt.runtimes = runtimes;
  end
  
% $$$   % set runspeed if provided
% $$$   if ~isempty(a.runspeed) && ...
% $$$         isfield(a.plots{k}, 'tunecurveopt') &&...
% $$$         ~isempty(a.plots{k}.tunecurveopt.filtparam),
% $$$     a.plots{k}.tunecurveopt.filtparamval = a.runspeed;
% $$$   end
  
  %% set colormaps/limits for all structs
  switch a.plots{k}.type
   case 'parmest'
    a.plots{k}.parmestdrawopt.cmap = cmap;
    a.plots{k}.parmestdrawopt.cmapwin = cmapwingray;
    a.plots{k}.parmestdrawopt.cmapwin_whitebg = cmapwingrayinv;
   case 'raster'
    a.plots{k}.rasterdrawopt.cmap = cmap;
    a.plots{k}.rasterdrawopt.clim = climrast;
    a.plots{k}.rasterdrawopt.cmapwin = cmapwinrast;
    a.plots{k}.rasteropt.nclcolors = nclustercols;
   case 'cont'
    a.plots{k}.contdrawopt.cmap = cmap;
    a.plots{k}.contdrawopt.cmapwin = cmapwinhot;
   case 'specgram'
    a.plots{k}.specgramdrawopt.cmap = cmap;
    a.plots{k}.specgramdrawopt.cmapwin = cmapwinjet;
  end
  
  %% pass in parmest params
  if strcmp(a.plots{k}.type, 'parmest'),
    %    a.plots{k}.dat.timebinsize = max([mintimebin a.timebinsize]);
    a.plots{k}.parmestdrawopt.autoscale = a.parmestautoscale;
    a.plots{k}.parmestdrawopt.sharpen = a.parmestsharpen;
    a.plots{k}.parmestdrawopt.scale = a.parmestscale;
    a.plots{k}.parmestdrawopt.gwtime = a.parmestgwtime;
    a.plots{k}.parmestdrawopt.gwdist = a.parmestgwdist;
  end
  
  %% pass in raster params
  if strcmp(a.plots{k}.type, 'raster'),
    a.plots{k}.rasteropt.maxpks = a.maxpks;
    %    a.plots{k}.rasteropt.rasterbg = a.rasterbg;
    if isempty(a.plots{k}.rasteropt.colormode),
      a.plots{k}.rasteropt.colormode = a.rastcolormode;
    end
  end
  
end
  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loop across subplots

for k = 1:nplots
  subax = subaxes(k);

  % shouldn't be necessary, as we pass in 'ax'
  %     set(figp,'currentaxes',subax);

  slcache = drawplot('plot', a.plots{k},...
                     'cache', slcache,...
                     'ax', subax,...
                     'mintimebin', mintimebin);
  

  if ~all([xstart xend] == xlim(subax)),
    error(['xlims not set correctly for plot ' num2str(k)]);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Do all 'overlay' plots, text

% special case of redraw with no clearing of ax--leave overlays alone
%set(figp,'CurrentAxes',ax);

cla(ax);
hold(ax,'on');

yl = ylim(ax);

%%% plot segments

for k = 1:nseglists
  if plotsegs(k), % do we want to plot these? 
    sl = a.seglists{k};
    segcol = segcolors(k,:);
    
    if size(sl,2) > 1,% get all bouts that are in range (no need to truncate)
      splot = sl((sl(:,1) < xend & sl(:,2) > xstart),:);
    else
      splot = sl(sl < xend & sl > xstart);
    end
    
    if ~isempty(splot),
      for seg = splot', % for iterates over columns
        if length(seg) == 1 || diff(seg) == 0,
          % draw lines for single-time events
          line ('parent', ax,...
                'XData', [seg(1) seg(1)],...
                'YData', yl,...
                'color', segcol,...
                'linewidth', 2,...
                'erasemode', 'xor');
        else
          switch(get(figp, 'renderer'))
           case 'opengl'
            rectangle ('parent', ax,...
                       'position', [seg(1) yl(1) seg(2)-seg(1) yl(2)-yl(1)], ...
                       'edgecolor', segcol,...
                       'linewidth', 2,...
                       'facealpha', 0.5,...
                       'facecolor', segcol);
           otherwise
            rectangle ('parent', ax,...
                       'position', [seg(1) yl(1) seg(2)-seg(1) yl(2)-yl(1)], ...
                       'edgecolor', segcol,...
                       'linewidth', 2,...
                       'erasemode', 'xor',...
                       'facecolor', 'none');
          end

        end
      end
    end
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% status text

statustext = [];
if a.randperm,
  statustext = ['RANDPERM  |  ' statustext];
end

text('parent', ax,...
     'string', statustext, ...
     'fontunits', 'pixels', ...
     'fontsize', 12, ...
     'units', 'normalized', ...
     'position', [1 1], ...
     'horizontalalignment', 'right', ...
     'verticalalignment', 'bottom',...
     'color', 'w');


% set xlim of overlay 
xlim(ax, [xstart xend]);



% put 'ax' on the bottom for ui purposes: in callbacks, the root
% 'currentobject' will now be in the subaxes, not in 'ax'. Since we use
% painters, all ax plots stay on top.
uistack(ax,'bottom');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% save appdata/cleanup

%%% set new appdata for use by callback fn's
% store argstruct 'a' in figure's appdata so that we can use it in 
% callback fn's.

setappdata(panel,'slargs',a);
sla = a;

% useful display params to store
sld.ax = ax;
sld.subaxes = subaxes;
sld.lastsubax = [];
sld.ename = a.ename;
sld.panelpos = panelpos;
sld.mints = mints;
sld.maxts = maxts;
sld.mints_cl = mints_cl;
sld.maxts_cl = maxts_cl;
sld.tstart = tstart;
sld.tend = tend;
sld.dstart = dstart;
sld.dend = dend;
sld.dstart_cl = dstart_cl;
sld.dend_cl = dend_cl;
sld.xstart = xstart;
sld.xend = xend;
sld.runtimes = runtimes;

% data analysis params to store
sld.clnos = clnos;

% prefs for callbackfns
sld.prefmousezoom = true;

setappdata(panel,'sldata',sld);
setappdata(figp,'slcache',slcache);

%%%% Tell the user their wait is over
set(figp,'pointer',oldptr);

