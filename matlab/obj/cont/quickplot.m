function [p, ax, cache, hs] = quickplot(c, varargin)
% QUICKPLOT - throw up a continuous data plot
%
% [p ax cache hs] = cplot(cdat, <param/value pairs>)
%
% cdat: contstruct to plot (optional if 'plotstruct' provided)
% 
% optional parameters
% 'plotstruct': output from slmkplots, or this function; can be used to
%               update params easily. (Can provide cdat!)
%      'chans': select channels to plot by number
% 'chanlabels': select channels to plot by name
%    'timewin': [tsart tend], or 'valid' (default: [] = plot all data)
%         'ax': axes to draw into
%    'whitebg': white or black background
%       'ylim': obvs
%      'color': colors to use for channels
%        'dpi': dpi for figures to be printed (subsample appropriately)
%  'subsample': whether to subsample data for plotting (default = true)
%      'cache': object cache
%
% Examples:
%
% [p] = quickplot(cdat, 'chans', 1, 'timewin', [10 20], 'ylim', [-5 5], 'yscale', 'log')
% % advance timewin
% [p] = quickplot([], 'plotstruct', p, 'timewin', [20 30]);


% TODO:
% -Specgrams
% -Segs
% -Multiple Cdats: contcombine and show on same axes
% -Pass in filtopts/envopts?
% -Just make a wrapper around spikeslpos instead. Get a GUI for free
% MAYBE:
% -optionally relabel channels (or just do that to cdats before plotting)

a = struct(...
    'plotstruct', [],...
    'chans', [],...
    'chanlabels', [],...
    'timewin', [],...
    'ax', [],...
    'whitebg', true,...
    'ylim', [],...
    'color', [],...
    'yscale', [],...
    'dpi', [],...
    'subsample', true,...
    'cache', mkcache());

a = parseArgsLite(varargin,a);

if ~ishandle(a.ax)
    error('''a.ax'' is not a valid handle (figure deleted?)');
end

if isempty(a.ax),
    a.ax = gca;
end

% handle this special case later
if strcmp(a.timewin, 'valid')
    tw_valid = true;
    a.timewin = [];
else
    tw_valid = false;
end

% reuse plotstruct cont object if provided
if isempty(c),
    if isempty(a.plotstruct) || ~strcmp(a.plotstruct.type, 'cont') || ~strcmp(a.plotstruct.dat.type, 'cont');
        error('Must provide a cdat (or a contobject as part of ''plotstruct''');
    end
    cont_template = a.plotstruct.dat;
else
    cont_template = [];
end

if ~isempty(a.plotstruct),
    contdrawopt = a.plotstruct.contdrawopt;
else
    contdrawopt = mkcontdrawopt;
end
contdrawopt.datalim =  a.ylim;
contdrawopt.subsample = a.subsample;
contdrawopt.color = a.color;

cont = mkcont('template', cont_template,...
    'contdata', c, ...
    'timewin', a.timewin,...
    'chans', a.chans,...
    'chanlabels', a.chanlabels);

% special case of 'valid' timewin
if tw_valid,
    cdat = cont.contdata;
    cont.timewin = [cdat.tstart+(cdat.nbad_start/cdat.samplerate)...
                       cdat.tend-(cdat.nbad_end/cdat.samplerate)];
elseif isempty(cont.timewin)
        cont.timewin = [cont.contdata.tstart cont.contdata.tend];
end

% make plot struct and call drawplot
p = slmkplots('type', 'cont');
p.dat = cont;
p.contdrawopt = contdrawopt;
p.whitebg = a.whitebg;

cla(a.ax);
[a.cache, hs] = drawplot(...
    'plot', p,...
    'cache', a.cache,...
    'ax', a.ax,...
    'dpi', a.dpi);

xlabel('Time (s)');

if ~isempty(c.units);
    ylabel(['(' c.units ')']);
end

ax = a.ax;
cache = a.cache;

