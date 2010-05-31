function c = contresamp(c,varargin)
% CONTRESAMP up- or downsample data in a cont struct, with antialiasing
%
%  [cout] = contresamp(c, [param/value pair args])
%
%  Usage note: to resample a cont struct to a specific samplerate (rather
%  than by a fraction), use the continterp function with the 'samplerate'
%  option. 
%  
% Inputs:
%  * c - a cont struct
%  *'resample'- fraction to up/downsample signal
%
%  Infrequently-used options:
%   'tol' - fractional tolerance for detecting integer resampling factors
%   'res_filtlen' - length of resampling filter to apply in non-decimation
%      case (defaults to 10)
%
% Example: To downsample a cont struct to a 5x lower sampling frequency
%
%  cout = contresamp(c, 'resample', 1/5);

% Tom Davidson <tjd@alum.mit.edu> 2003-2010 
  
  a = struct(...
      'resample',[],...
      'tol', 0.001,...
      'res_filtlen', 10);
  
  a = parseArgsLite(varargin,a);

  % keep it, for later
  datatype = class(c.data);
  [nrows ncols] = size(c.data);

  if strcmp(datatype, 'logical')
    error('Can''t filter data of type ''logical''');
  end
  
  if ~isempty(a.resample),
    
    if a.resample == 1,
      disp('Resample factor == 1, nothing to do...');
      return
    end
      
    
    if abs((fix(1/a.resample) / (1/a.resample)) - 1) < a.tol,
      % integer factor, we can use decimate
      dec_f = fix(1/a.resample);
      res_f = 1/dec_f;

      % use a 30-pt FIR filter to conserve memory ('decimate' usually uses
      % IIR/filtfilt)
      filtlen = 30; % the default for decimate, just being explicit
      
      % pre-allocate decimated array, preserving datatype of c.data
      data_dec = zeros(ceil(nrows/dec_f), ncols, datatype);
      disp('decimating...');
      for col = 1:ncols,
        % cast to double and back since 2007a's decimate doesn't like single
        % datatype (per Greg Hale)
        data_dec(:,col) = cast(decimate(double(c.data(:,col)),...
                                        dec_f, ...
                                        filtlen,...
                                        'fir'), ...
                               datatype);
      end

      assert(strcmp(class(data_dec), datatype), ...
             'wrong datatype for decimated data');

      c.data = data_dec;
      clear data_dec;
      
      c.samplerate = c.samplerate/dec_f;

      % decimate has 0 group delay (i.e. first point represents same
      %start time), but since it takes every 'rth' point after the first, the last
      %sample in the sequence will likely not be from the same time as the
      %last sample in the input. Recalculate its actual time from the samplerate
      c.tend = c.tstart + ((size(c.data,1)-1) ./c.samplerate);
      
    else
      
      % 'resample'
      % use a wider tolerance than the default (1e-6) to get smaller terms
      % for resampling (use continterp for very precise control over sampling
      % rates/times)
      [res_num res_den] = rat(a.resample, a.resample.*a.tol);
      res_f = res_num/res_den;
      
      filtlen = a.res_filtlen;
      
      % pre-allocate
      data_res = zeros(ceil(nrows*res_f), ncols);
      for col = 1:ncols,
        disp('resampling...');
        % upfirdn (called by resample) can't handle 'single' type data. bug
        % filed with mathworks 11/14/06, confirmed by Mathworks as fixed in R14SP2)
        data_res(:,col) = resample(double(c.data(:,col)),...
                                   res_num,res_den,...
                                   filtlen);
      end
      c.data = cast(data_res,datatype); % back to original data type
      clear data_res;
      
      % we could also recalc this by determining where new tend is (old tend -
      % some samples). But for, e.g. 1e7 samples, at around 1kHz, float
      % error is around 1e-6, i.e. 1 microsecond, and this is way simpler. I
      % think we're okay.
      c.samplerate = c.samplerate*res_f;
      
      % for resample, the filtlen param is proportional to the length filter
      % used. doc/help resample say filtlen samples of *input* are used, so
      % we're going to mark as bad that many at start and end, then
      % multiply by 

      
    end
      
    % # of unreliable samples in resampled signal = filter length/dec_f
    c.nbad_start = ceil(c.nbad_start * res_f) + ceil(filtlen * res_f);
    c.nbad_end = ceil(c.nbad_end * res_f) + ceil(filtlen * res_f);

    % calculate new data range
    c = contdatarange(c);

  else
    disp('no resample requested');
  end
  
  