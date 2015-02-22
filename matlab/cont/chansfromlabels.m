function chans = chansfromlabels(clabels, alabels)
% CHANSFROMLABELS: get cont channel index from chanlabel/s
%
% chansfromlabels(c, alabels)
% chansfromlabels(clabels, alabels)
%  c: contstruct, with 'chanlabels' field
%  clabels: cell array of chanlabels
%  alabels: string or cell array of strings to match.
%
% Order of inputs is preserved, so can be used for reordering channels.
%
% Examples:
% idx = chansfromlabels(cdat, 'Ref1X');
% idx = chansfromlabels(cdat, {'Ref1X' 'Ref1Y' 'Dat1'});
% idx = chansfromlabels(cdat.chanlabels, {'Ref1X' 'Ref1Y' 'Dat1'});

  % make a 1x1 cell array out of a string input
  if ischar(alabels) && ~isempty(alabels),
    alabels = {alabels};
  end
  
  % allow passing in a cont struct
  if isstruct(clabels) && isfield(clabels, 'chanlabels');
      clabels = clabels.chanlabels;
  end
  
  chans = [];
  if ~isempty(alabels) && ~isempty(clabels)
      
    if length(unique(clabels)) ~= length(clabels),
      error('repeated chanlabels in cont struct');
    end
    
    if length(unique(alabels)) ~= length(alabels),
      error('repeated chanlabels requested');
    end

    for chanstr = alabels(:)'
        ch = find(strcmp(chanstr,clabels));
        if isempty(ch)
            error(['requested chanlabel not found: ' chanstr{:}]);
        end
        chans = [chans ch];
    end
    
  end
