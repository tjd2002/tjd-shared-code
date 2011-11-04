function [sg cont] = mkspecgram(varargin)
% MKSPECGRAM make a specgram object 
% 
%  NOTE: unlike matlab's specgram fn, we report the times at the
%  *centers* of the bins in sg.b
%
%  NOTE: we also expand the requested timewin so that we have
%  data for all the requested times
  
  sg = struct(...
      'type', 'specgram',...
      ...
      'contdata', [],... % contdata struct to operate on
      'contvar', [],... % name of contdata variable in base workspace
      'chans', [],... % channels (columns) of contdata to analyze
      'chanlabels',[],... % names of channels to analyze
      'timewin', [],...
      'contopt', [],... % contdata filtering, etc
      'specgramopt',[],... % method, etc
      ...
      'b', [],...  % specgram
      'f', [],...  % frequency at *center* of each bin in sg.b
      't', [],...  % time at *center* of each bin in sg.b
      't_window', [],... % time width of fft window 
      't_overlap', [],... % time width of overlapping region
      'label', [],... % sg.contdata.chanlabels, saved for labeling plots
          ...
      'template', [],...
      'cache', [],...
      'cache_hit', false);

  sg = obj_reparse(sg, varargin);

  sgopt = sg.specgramopt;

  if ~xor(isempty(sgopt.window_time), isempty(sgopt.window_samp)),
    error(['exactly one of ''window_time'' or ''window_samp'' must be ' ...
           'provided in specgramopt']);
  end
  
  sg = obj_cachesearch(sg);
  
  if sg.cache_hit,
    sg = obj_cleanup(sg);
    cont = [];
    return;
  end
  
  %%% no cache hit, calculate fresh
  
  cont = mkcont('contdata', sg.contdata,...
                'contvar', sg.contvar,...
                'chans', sg.chans,...
                'chanlabels', sg.chanlabels,...
                'timewin', sg.timewin,...
                'contopt', sg.contopt,...
                'cache', sg.cache);

  if size(cont.contdata.data,2) ~= 1,
    error('mkspecgram can only operate on one channel at a time');
  end
  
  samplerate = cont.contdata.samplerate;

  % convert window and overlap from seconds/fractions to samples
  if ~isempty(sgopt.window_time),
    window_samp = round(sgopt.window_time * samplerate);
  else
    window_samp = sgopt.window_samp;
  end
  
  overlap_samp = round(sgopt.overlap_frac * window_samp);
  bin_samp = window_samp - overlap_samp;
  
  % expand analysis window so we get (centered) bins out past the edges
  % of the requested timewindow
  spectimewini = cont.timewini + [-window_samp window_samp];

  spectimewini(1) = max([1 spectimewini(1)]);
  spectimewini(2) = min([size(cont.contdata.data,1) spectimewini(2)]);
  
  
  if diff(spectimewini)/bin_samp > sgopt.maxcols,
    sg = obj_cleanup(sg);
    return;
  end
  
  specdata = cont.contdata.data(spectimewini(1):spectimewini(2));
  
  switch sgopt.method,
   case 'specgram',
    [sg.b sg.f sg.t] = specgram(...
        specdata, ...
        sgopt.nfft, ...
        samplerate,...
        window_samp,...
        overlap_samp);
    
   case 'multitaper',
    [sg.b sg.f sg.t] = mtpsg(...
        specdata, ...
        sgopt.nfft, ...
        samplerate,...
        window_samp,...
        overlap_samp,...
        sgopt.mt_NW,...
        sgopt.mt_detrend);
  end
  
  
  if sgopt.multiply_by_f,
      sg.b = bsxfun(@times, sg.b, sg.f);
  end

  % actual time widths of window/overlap (useful for image/surf plot offsets)
  sg.t_window = window_samp / samplerate;
  sg.t_overlap = overlap_samp / samplerate;
  
  % save chan label for later plotting
  sg.label = cont.contdata.chanlabels;
  
  %%% make sg.b and sg.f reflect *centers of bins*
  
  % sg.b is already 'centered' on sg.f values
  
  % sg.t currently reflects start times of bins--we want it to reflect
  % the *centers*, like sg.f.:
  %
  %    sg.t = sg.t + 
  %           actual time of first sample used +
  %           1/2-window offset
  offset = cont.contdata.tstart + (spectimewini(1)-1)/samplerate + ... 
         sg.t_window/2; 
  
  % have to use linspace (rather than adding offset to each value) so
  % that when we later try to detect uniformity (in draw2d), diff(sg.t,2)
  % <= eps(max(sg.t)).
  sg.t = linspace(sg.t(1) + offset, sg.t(end) + offset, length(sg.t));
  
  %%% ???
  % when MATLAB plots a specgram, it shifts it by 1/2 (the window width
  % less the overlap). I think this is incorrect, since we want to know
  % the center of the range of data used to generate the FFT.
  
  % strip out contdata, e.g.
  sg = obj_cleanup(sg);