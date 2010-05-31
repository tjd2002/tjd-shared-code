function c = imcontsegs(varargin)
% IMCONT make contdata struct from chunks of timeseries  data 
%
% possible inputs:
%  
% -name/labels/units
% -samplerate (cdat will have same samplerate as data)
% -seg start times
% -multi-chan (same segs)
% -align/interpolate
% -default value (NaN)
%
% TODO
% -existing cdat (select/crop segs) easy&inefficient with timestamps
% -compact repr'n?

  
  a = struct(...
      'timewin', [],... % timewin for cdat
      'starts_t',[],... % seg starts, (s)
      'data', {{}},... % seg data (cell array)
      'rowdata', false,... % data in rows
      'samplerate', [],... 
      'alignok', false,... % suppress align warning?
      'name', '',... % cdat name
      'units', '',...
      'chanvals', [],...
      'chanlabels', {{}},...
      'dataunits', '' ,...
      'datatype', 'single'); % warning: may be converted to double at
                             % various points in filtering, etc...
  
  a = parseArgsLite(varargin,a);
  
  c = mkcdat('name', a.name,...
             'chanvals', a.chanvals,...
             'chanlabels', a.chanlabels,...
             'samplerate', a.samplerate,...
             'units', a.units);

  
  nsegs = numel(a.starts_t);

  if ~a.rowdata
    nchans = size(a.data{1},2);
  else
    nchans = size(a.data{1},1);
  end
  
  if numel(a.data) ~= numel(a.starts_t)
    error('Must provide same # of seg starts and data chunks');
  end

  if ~issorted(a.starts_t)
    [a.starts_t sorti] = sort(a.starts_t);
    a.data = a.data(sorti);
  end

  % get timewin to import
  if isempty(a.timewin)
    a.timewin(1) = a.starts_t(1);
    a.timewin(2) = a.starts_t(end)+ (numel(a.data{end})-1)/a.samplerate;
    warning(sprintf('No ''timewin'' provided, using extent of data: [%d %d]', a.timewin)); %#ok
  end
    
  % align data to sampling rate (if starts_t are not aligned)
  startsi = (a.starts_t-a.timewin(1)).*a.samplerate + 1;
  aligns_t = (startsi-round(startsi))/a.samplerate;
  startsi = round(startsi);
  
  % warn if jitter of > 1% of samplerate was introduced
  if ~a.alignok && any(aligns_t>1/a.samplerate/100),
    warning(['start times rounded > 1% of sample period, largest align = ' ...
             num2str(max(abs(aligns_t(:)))*1000) ' ms']);
  end
  startsi = round(startsi);
  
  for k = 1:nsegs
    if a.rowdata
      a.data{k} = a.data{k}';
    end
    dat(startsi(k):(startsi(k)+size(a.data{k},1)-1),:) = a.data{k};
  end
  
  c.data = dat;
  c.tstart = a.timewin(1);
  c.tend = c.tstart+(size(dat,1)-1)/a.samplerate;
  
  c = contdatarange(c);