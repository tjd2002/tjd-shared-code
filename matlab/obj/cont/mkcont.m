function [cont filt] = mkcont(varargin)
% MKCONT make a cont object to be plotted by drawcont
  
  cont = struct(...
      'type', 'cont',...
      ...
      'contvar', [],... % name of contdata variable in base workspace
      'chans', [],... % channels (columns) of contdata to plot
      'chanlabels',[],... % names of channels to plot
      'timewin', [],... % time window (use all)
      'contopt',[],... % see mkcontopt.m
      ...
      'contdata', [],... % contdata struct to use (input & output)
      'timewini', [],... % indexes into rows of contdata.data
      ...
      'template', [],...
      'cache', [],...
      'cache_hit', false);

  cont = obj_reparse(cont, varargin);

  if isempty(cont.contopt),
    cont.contopt = mkcontopt();
  end
  
  % update timewin on cached object
  cachecont = obj_cachesearch(cont);
  
  if ~cachecont.cache_hit,
    
    % get full contdata struct from base workspace or inputs
    cont = sub_getcontdata(cont);
    
    % select channels
    if ~isempty(cont.chans) || ~isempty(cont.chanlabels),
        cont.contdata = contchans(cont.contdata, ...
                                  'chans', cont.chans,...
                                  'chanlabels', cont.chanlabels);
        % we successfully selected them and replaced the contdata struct
        cont.chans = [];
        cont.chanlabels = [];
    else
        % keep all channels
    end
    
    % filter if requested
    if ~isempty(cont.contopt.filtopt),
      
      % use old name as newname to help caching
      [cont.contdata filt] = contfilt(cont.contdata,...
                                      'newname', cont.contdata.name,... 
                                      'filtopt', cont.contopt.filtopt,...
                                      'autoresample', cont.contopt.autoresample,...
                                      'cache', cont.cache);
    end

    % calculate the envelope if needed
    if ~isempty(cont.contopt.envopt),
      % use 'nosuffix' to help caching
      cont.contdata = contenv(cont.contdata, 'envopt', cont.contopt.envopt, ...
                              'nosuffix', true);
    end
    
  else
    % we got a cache hit, update and use it

    % use the newly-requested timewin
    cachecont.timewin = cont.timewin;
    
    cont = cachecont;
    
  end
  
  % we get a new timewini even if there is a cache_hit, since the cached
  % object is a whole contdata object
  cont = sub_gettimewini(cont);
  cont = obj_cleanup(cont);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunctions

% get contdata from base workspace, if necessary
function cont = sub_getcontdata(cont)

  cdat = cont.contdata;
  cvar = cont.contvar;
  
  switch sum([isempty(cdat) isempty(cvar)]),
   case 2,
    error(['One of ''contdata'' or ''contvar'' must be provided for each ' ...
           '''cont'' struct']);
   case 0,
    error(['Only one of ''contdata'' or ''contvar'' may be provided for each ' ...
           '''cont'' struct']);
   otherwise
    % OK!
  end
  
  % full contdata struct in cont
  if ~isempty(cdat),
    if ~isstruct(cdat),
      error('''contdata'' must be a struct');
    else
      cont.contdata = cdat; %#ok
    end
  else
    
    % name of contdata var in base workspace provided
    if ischar(cvar) && evalin('base', ['exist(''' cvar ''', ''var'')']),
      cont.contdata = evalin('base',cvar);
    else
      error(['''contvar'' must be name of contdata struct in base ' ...
             'workspace']);
    end
  end
  
  
function cont = sub_gettimewini(cont)
  
% timewini is in sample indexes; max/min keeps us within range of data;
% floor/ceil ensures we go just past start and end of range to be
% plotted, to cover whole plot
  
% this used to use xstart/xend, rather than timewin... why?

  if isempty(cont.timewin),
    cont.timewini = [1 size(cont.contdata.data,1)];
  else
    
    cont.timewini(1) = max([1, ...
                        floor((cont.timewin(1) - cont.contdata.tstart) ...
                              * cont.contdata.samplerate)]);
    
    cont.timewini(2) = min([size(cont.contdata.data,1), ...
                        ceil((cont.timewin(2) - cont.contdata.tstart) ...
                             * cont.contdata.samplerate)+1]);
  end
  
  if any(cont.timewini < 1) || (cont.timewini(1) >= cont.timewini(2)),
    error('bad ''timewini'' indexes for cont (requested time out of data range?)');
  end

  
