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

  % data integrity check
  contcheck(c);
  
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
        % cast to double and back since decimate doesn't like single
        % datatype (since R2007a or earlier)
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

      % samplerate, tstart, and tend are updated below
      
    else
      
      % 'resample'
      % use a wider tolerance than the default (1e-6) to get smaller terms
      % for resampling (use continterp for very precise control over sampling
      % rates/times)
      [res_num res_den] = rat(a.resample, a.resample.*a.tol);
      res_f = res_num/res_den;
      filtlen = a.res_filtlen;
      
      % pre-allocate
      % (note must use res_num/res_den rather than res_f for this calculation,
      % since this is what 'resample' uses internally, and in some cases
      % (e.g. 6724096*47/46) float error would cause the 2 to give a different
      % answer).
      data_res = zeros(ceil(nrows*res_num/res_den), ncols, datatype);
      for col = 1:ncols,
        disp('resampling...');
        % resample can't handle 'single' type data. bug
        % filed with mathworks 11/14/06. R14SP2 and later correctly error
        % on 'single' inputs)
        data_res(:,col) = cast(resample(double(c.data(:,col)),...
                                        res_num,res_den,...
                                        filtlen),...
                               datatype);
      end
      c.data = data_res;
      clear data_res;
    
          
      % samplerate, tstart, tend updated below

    end

    % both 'resample' and 'decimate' have 0 group delay (so tstart
    % remains the same), but the time of the last sample (tend) changes since it
    % must occur at an integer number of sample intervals from tstart. We
    % calculate this time from the size of the data and the new samplerate.
    % (This method is subject to some float error, but for, e.g. 1E7 samples,
    % this error is only 1 microsecond, which is acceptable).
    
    c.samplerate = c.samplerate*res_f;
    c.tend = c.tstart + ((size(c.data,1)-1) ./c.samplerate);
    
    % the filtlen param is proportional to the length filter used. doc/help
    % resample say filtlen samples of *input* are used, so we're going to
    % mark as bad that many at start and end, then multiply by res_f
    c.nbad_start = ceil(c.nbad_start * res_f) + ceil(filtlen * res_f);
    c.nbad_end = ceil(c.nbad_end * res_f) + ceil(filtlen * res_f);

    % calculate new data range
    c = contdatarange(c);

  else
    disp('no resample requested');
  end
  
  % data integrity check
  contcheck(c);
