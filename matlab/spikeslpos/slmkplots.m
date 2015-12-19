function plots = slmkplots(varargin)
%SLMKPLOTS-make a plot definition struct for passing into spikeslposs
% plots = slmkplots(plotargs)
%
% arglist or several arglists, each a cell array containing param/value pairs
%
%     'height': arbitrary units, all plots will be scaled to fill fig (1)
%       'type': type of display 'raster'/'posest'/'cont'/('none')
%      'yflip': change direction of y-axis (false)/true
%    'whitebg': plot has white or black background? false/(true)
%     'dat.*' : OBJ-specific data
%     'OBJopt': OBJopt options structure ([] = default from mkOBJopt)
% 'OBJdrawopt': OBJdrawopt options structure ([] = default from mkOBJdrawopt)
%
%  e.g.:
%       plot = slmkplots('type','parmest','opt',peopt);
%
%
%       slmkplots(...
%          {'type','parmest'},...
%          {'type','parmest'}, ...
%          {},... % type='none', height=1
%          {'type','raster', 'height', 2});
%
% $Id$

plotdefaults = struct( ...
    'height', 1, ...
    'type', 'none',...
    'yflip', false,...
    'whitebg', false,...
    'plotno', [],... % number of plot in spikeslpos figure
    'plottypeno', [],... % number of this type of plot in spikeslpos figure
    'dat', struct()); % so we can test for isfield

% just easier variable name
plotargs = varargin;

% default is to create a single default plot with args from plotdefaults
if isempty(plotargs),
  plotargs = {{}};
end

if ~iscell(plotargs{1}),
  % if arglist, convert to cell array of arglists for processing below
  plotargs = {plotargs};
end

for k = 1:length(plotargs);
  plots{k} = parseArgsLite(plotargs{k},plotdefaults);

  plots{k}.dat.timewin = []; % provided by slp
  plots{k}.dat.timewin_plot = []; % provided by slp - xlim of plot

  switch plots{k}.type
   case 'parmest'
    plots{k}.dat.e = []; % provided by slp
    plots{k}.dat.ename = ''; % provided by slp
    plots{k}.dat.clnos = [];
    plots{k}.dat.clperm = [];
    plots{k}.dat.timebinsize = []; 
    plots{k}.dat.overlapsperbin = []; 
    plots{k}.parmestopt = mkparmestopt();
    plots{k}.parmestdrawopt = mkparmestdrawopt();
    plots{k}.tunecurveopt = mktunecurveopt();
    plots{k}.parmdrawopt = mkparmdrawopt();
    plots{k}.whitebg = true;
    
   case 'raster'
    plots{k}.dat.e = []; % provided by slp/drawplot
    plots{k}.dat.ename = ''; % provided by slp
    plots{k}.dat.clnos = [];
    plots{k}.dat.clperm = [];
    plots{k}.rasteropt = mkrasteropt();
    plots{k}.rasterdrawopt = mkrasterdrawopt();
    plots{k}.tunecurveopt = mktunecurveopt();
    plots{k}.tcpeakopt = mktcpeakopt();
    plots{k}.parmdrawopt = mkparmdrawopt();
    
   case 'parm'
    plots{k}.dat.e = []; % provided by slp
    plots{k}.dat.ename = ''; % provided by slp
    plots{k}.tunecurveopt = mktunecurveopt();
    plots{k}.parmdrawopt = mkparmdrawopt();
    
   case 'cont'
    plots{k}.dat.contdata = [];
    plots{k}.dat.contvar = [];
    plots{k}.dat.contname = [];
    plots{k}.dat.chans = [];
    plots{k}.dat.chanlabels = [];
    plots{k}.contopt = mkcontopt();
    plots{k}.contdrawopt = mkcontdrawopt();
   
   case 'specgram'
    plots{k}.dat.contdata = [];
    plots{k}.dat.contvar = [];
    plots{k}.dat.contname = [];
    plots{k}.dat.chans = [];
    plots{k}.dat.chanlabels = [];
    plots{k}.contopt = mkcontopt();
    plots{k}.specgramopt = mkspecgramopt();
    plots{k}.specgramdrawopt = mkspecgramdrawopt();
   
   case 'none'
    
   otherwise
    error([mfilename ':BadArgs'], ['Bad value for ''type'' arg: ' plots{k}.display]); 
  end
end

% if only making one, just return the plot struct
if length(plots) == 1,
  plots = plots{1};
end