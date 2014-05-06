function cout = contreref(c, varargin);
% CONTREFERENCE - rereference cdat from one channel to another
%
%  cout = contreref(c, [name/value pairs])
%
% Case 1) When all channels in a group are recorded against ground (I.e.
% recorded single-ended), we can reference by subtracting that channel from
% all others.
%
% Case 2) When all channels but one in a cdat are recorded against a
% reference, and *that reference* is recorded against ground, we can
% recover the single-ended recordings and re-reference as in case 1.
%
% Inputs: (* means required, -> indicates default value)
%   * c - a cont struct to rereference
%   'oldrefchan'/'oldrefchanlabel' - the channel in c currently used as a
%       reference (and recorded against ground). [] == all channels
%       recorded against ground (i.e. single-ended recording).
%   'newrefchan'/'newrefchanlabel' - the channel in c to use as a reference
%       for cout. [] == return 'dereferenced' cdat: all channels vs. gnd.
%
%   
% Outputs:
%   cout - the resulting cont struct
%
% Example:
%  Rereference LFP channels from one channel to another
%
%

%  Tom Davidson <tjd@alum.mit.edu> 2003-2014

% TODO

a = struct('name', [],...
  'oldrefchan', [],...
  'oldrefchanlabel', [],...
  'newrefchan', [],...
  'newrefchanlabel', []);

a = parseArgsLite(varargin,a);

% validate inputs

if ~isempty(a.oldrefchanlabel)
  if ~isempty(a.oldrefchan)
    error('Can only provide 1 of ''oldrefchan''/''oldrefchanlabel''');
  end
  oldrefchan = chansfromlabels(c, a.oldrefchanlabel);
else
  oldrefchan = a.oldrefchan;
end

if ~isempty(a.newrefchanlabel)
  if ~isempty(a.newrefchan)
    error('Can only provide 1 of ''newrefchan''/''newrefchanlabel''');
  end
  newrefchan = chansfromlabels(c, a.newrefchanlabel);
else
  newrefchan = a.newrefchan;
end

nchans = size(c.data,2);

%% First rerefernce to Ground
if isempty(oldrefchan)
  fprintf('Assuming all channels are recorded against ground.\n');
  c_gnd = c;
  
else
  fprintf('Re-referencing all channels to ground (assuming that reference channel #%d, ''%s'' is recorded against ground, and that all other chans were referenced to chan #%d.\n',...
    oldrefchan, c.chanlabels{oldrefchan}, oldrefchan);
  
  c_gnd = c;
  
  for k = 1:nchans
    if k == oldrefchan,
      c_gnd.data(:,k) = c.data(:,k);
    else
      c_gnd.data(:,k) = c.data(:,k)+c.data(:,oldrefchan);
    end
    c_gnd.chanlabels{k} = [c.chanlabels{k} '_vs_Gnd'];
  end
  c_gnd = contdatrange(c_gnd);
end

%% Then rereference to new ref, if requested

if isempty(newrefchan)
  % return c_gnd, all chans recorded against ground
  cout = c_gnd;
  return
  
else
  fprintf('Re-referencing all channels to new reference ''%s''. (Assuming that all channels now referenced to ground)\n', ...
    c.chanlabels{newrefchan});
  
  c_reref = c; 
     
  for k = 1:nchans
    if k == newrefchan,
      % do nothing, leave new ref as Gnd-referenced
      c_reref.chanlabels{k} = [c.chanlabels{k} '_vs_Gnd'];
    else      
      c_reref.data(:,k) = c.data(:,k)-c.data(:,newrefchan);
      c_reref.chanlabels{k} = [c.chanlabels{k} '_vs_' c.chanlabels{newrefchan}];
    end
  end
  c_reref = contdatarange(c_reref);

end

cout = c_reref;
