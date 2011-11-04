function chans = chansfromlabels(clabels, alabels)

  % make a 1x1 cell array out of a string
  if ischar(alabels) && ~isempty(alabels),
    alabels = {alabels};
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
