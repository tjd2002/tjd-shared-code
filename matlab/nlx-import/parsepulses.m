function [out] = parsepulses(varargin)
% PARSEPULSES Parse output of NlxLoadEV into bouts of laser stim
%
% Finds 'bouts' of repeating 'trains' of 'pulses'
% Bouts are separated by at least minboutgap
% Trains are identified by repetition
%
% 'ev': struct returned from NlxLoadEV
% 'chans': which channel to use
% 'minboutgap': gap to separate bouts
% 'tol': tolerance factor (def 10000)
% 'traintol': tolerance factor for stdev of train start times (def 1000)
%
%
% TODO
%  -group pulses that occur in sequence, not across a whole experiment?
%  -summary figure with overlapping pulse traces aligned by train
%   (interactive mode?)

p = inputParser();
p.addRequired('ev');
p.addRequired('chans');
p.addParamValue('minboutgap', 1, @isnumeric);
p.addParamValue('tol', 10000, @isscalar); % tolerance factor: pulse/gap
% durations are rounded to
% nearest 1/tol
p.addParamValue('traintol', 1000, @isscalar); % tolerance factor for stdev of
% train start times

p.addParamValue('maxtrainlength', 8); % max number of pulses in a train
% to look for

p.parse(varargin{:});
a = p.Results;

if numel(a.chans) ~= 1,
    error('only one pulse channel currently supported');
end

lp = a.ev.pulses{a.chans};

if isempty(lp),
    out.args = a;
    out.bouts = [];
    out.groups = [];
    return
end


% get pulse duration, time gap to next pulse
lp = [lp subf_getdurgap(lp)];

% separate into bouts (periods with pulse gap < 'minboutgap')
m = 1;
b(m).pulses = lp(1,:); % first pulse always begins a new bout;
for n = 2:size(lp,1)
    lastp = lp(n-1,:);
    thisp = lp(n,:);
    gap = lastp(4);
    if gap < a.minboutgap
        b(m).pulses = [b(m).pulses; thisp];
    else % start new bout with pulse after gap
        m = m + 1; % inc bout counter
        b(m).pulses = thisp;
    end
end


for m = 1:numel(b)
    bp = b(m).pulses;
        
    % Get overall duration of bout
    b(m).timewin = [bp(1,1) bp(end,2)];
    b(m).dur = diff(b(m).timewin);
    
    if size(bp,1) == 1,
        % only 1 pulse in bout = no train, or pattern
        b(m).trainfreq = [];
        b(m).trainpattern = round(bp(3)*a.tol)/a.tol;
        b(m).trainfirstpulsei = 1;

    else
        
        % Figure out number of pulses in a train
        % want at least 2 instances of train to call it a 'pattern' (also need this
        % for massaging of pattern below)
        maxtlen = min(ceil(size(bp,1)./2), a.maxtrainlength);
        for lag = 1:maxtlen
            % find the smallest chunklength that results in all roughly equal
            % (within traintol) gaps between start of trains.
            % conv(A, [1 zeros -1], 'valid') does sliding diff with varying lags.
            lagdiff = conv(bp(:,1), [1 zeros(1, lag-1) -1], 'valid');
            tstdev(lag) = std(lagdiff);
        end
        
        tlen = find(tstdev<1/a.traintol, 1, 'first');
        if isempty(tlen)
            warning('No pulse train identified within bout for pulses:');
            % disp(bp);
            %         warning('No pulse trains < ''traintol'', using minimum instead');
            %         [~, tlen] = min(tstdev);
        end
        
        
        % if we found a repeating pattern, write it out
        if ~isempty(tlen),
            
            %% get pattern of median durs/gaps (median helps ignore last long gap, noise
            %% at start)
            
            % select only durs/gaps for full repetitions
            patreps = floor(size(bp,1)/tlen);
            durgaps = bp(1:tlen*patreps,3:4);
            
            % reshape into a 3D array of pattern instances (tlen x patreps x 2)
            patarr = reshape(durgaps, tlen, patreps, 2);
            
            % collect median value for each element in the pattern
            pat = reshape(median(patarr,2),[],2,1);
            
            % get pattern repetition rate, in Hz, before rounding
            b(m).trainfreqhz = 1/sum(pat(:));
            
            % round durs/gaps to tolerance (usu 0.1ms)
            pat = round(pat*a.tol)/a.tol;
            
            % shift pattern so that longest gap comes last
            [~, maxgapi] = max(pat(:,2));
            b(m).trainpattern = circshift(pat, tlen-maxgapi);
            b(m).trainfirstpulsei = mod(maxgapi,tlen)+1; % index into first pulse in first train
            
            % write out 'quality' of train
            b(m).trainstdev = tstdev(tlen);
        else
            % none found
            b(m).trainfreq = [];
            b(m).trainpattern = [];
            b(m).trainfirstpulsei = 1;
        end
     
    end
    
    %% save out all pulse and gap durations seen during the bout (helpful
    %% for debugging/assessing quality)
    
    % round durs to nearest 1/tol
    [b(m).unique.pulsedurs , ~, undurj] = unique(round(bp(:,3).*a.tol)./a.tol);
    
    % count # of pulses ofeach duration
    for k = 1:numel(b(m).unique.pulsedurs)
        b(m).unique.durcounts(k) = sum(undurj==k);
    end
    
    % last gap is gap to next bout, don't include it
    [b(m).unique.pulsegaps , ~, ungapj] = unique(round(bp(1:end-1,4).*a.tol)./a.tol);
    
    for k = 1:numel(b(m).unique.pulsegaps)
        b(m).unique.gapcounts(k) = sum(ungapj==k);
    end
end

%% group bouts by stim train pattern

% build array of train pattern signatures, 1-per-row, padded with
% zeros, so we can use 'unique' to find groups
for m = 1:numel(b)
    tp = b(m).trainpattern(:)';
    if ~isempty(tp)
        tpr(m, 1:numel(tp)) = tp;
    else
        tpr(m,:) = 0;
    end
end
% group bouts by train pattern
[~, firsti stimi] = unique(tpr, 'rows', 'first');

% get original pattern arrays, not 1-per-row version
patterns = {b(firsti).trainpattern};

% renumber groups by occurrence of first bout ('unique' doesn't do this)
[~, si] = sort(firsti); % si gives reordering vector
patterns = patterns(si);
rnk(si) = 1:numel(si);
stimi = rnk(stimi);


%% write out groups
for k = 1:numel(patterns)
    g(k).pattern = patterns{k}; %#ok
    bouti = find(stimi==k);
    g(k).bouti = bouti; %#ok
    for m = bouti(:)'
        b(m).groupi = k;
    end
end

%% make pretty name for each group
for k = 1:numel(g),
    gpatms = 1000 .* g(k).pattern;
    gtraindurms = sum(gpatms(:));
    gtrainfreqhz = 1000/gtraindurms;
    
    % no pattern detected
    if isempty(g(k).pattern),
      g(k).name = 'No detected pattern';
    
    % single pulse
    elseif size(g(k).pattern,1)==1
        if numel(g(k).pattern) == 1,
            % non-repeating, only 1 pulse in bout
            g(k).name = sprintf('%gms, non-repeating',...
                gpatms(1,1)); % pulse dur
        else
            % repeating within bout        
            g(k).name = sprintf('%gms @ %gHz',...
                gpatms(1,1),... % pulse dur
                gtrainfreqhz); % train repetition rate
        end
        
    % repeating train of regular pulses
    elseif all(gpatms(:,1)==gpatms(1,1)) && ... % pulse durs equal
            all(gpatms(1:end-1,2)==gpatms(1,2)) % within-train gaps equal
        gsubfreqhz = 1000/sum(gpatms(1,:));
        g(k).name = sprintf('(%u x %gms @ %gHz) @ %gHz',...
            size(gpatms,1),... % within-train repetitions
            gpatms(1,1),... % pulse dur
            gsubfreqhz,... % within-train freq
            gtrainfreqhz); % train repetition rate
        
        % something more complicated: just write out pattern of pulse durs/gaps
    else
        g(k).name = ['[' ...
            sprintf('%gms (%gms) ', gpatms(:))... % repeats
            ']' ...
            sprintf('@ %gHz', gtrainfreqhz)];
    end
end


out.args = a;
out.bouts = b;
out.groups = g;

function durgap = subf_getdurgap(datxi)
% recompute durgap array
if isempty(datxi),
    durgap = [];
else
    % the 'diff' between crossings is alternately a gap or a duration
    % col 1 = samps btw. up xing & dn xing (= durations)
    % col 2 = samps btw. down xing & up xing (= gaps)
    datxi = shiftdim(datxi,1);
    durgap = [diff(datxi(:)); Inf];
    durgap = shiftdim(reshape(durgap,2,[]),1);
end


