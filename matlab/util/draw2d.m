function [hs, zlims] = draw2d(varargin)
% DRAW2D plot 2-D data, w/ shading, linear/log scale... (chooses
% image/surf as needed). See code for input options
%
% Author: Tom Davidson (tjd@alum.mit.edu) 2007
%
% $Id$
  
% todo: -subsampling/decimating of large data ?
%       -surf with rgb triplets (for log scale or non-uniform plotting)

  hs = [];
  
  def = struct(...
      'ax', [],... % axes to draw into (default gca)
      ...
      'data', [],... % the data to plot, can be 2d or 3d array
      ...
      'xhistctrs',[],... % location of centers of pixels
      'xhistedges',[],... % location of edges of pixels
      'xhistendedges', [],... % location of outer edges of 1st and last pixels
      'xhistendctrs', [],... % location of centers of first and last pixels
      ...
      'yhistctrs',[],... % ditto, but for y axes
      'yhistedges',[],...
      'yhistendedges', [],...
      'yhistendctrs', [],...
      ...
      'zlimscale', false,... % scale z data to max
      'zcolnorm', false,... % renormalize columns to sum to 1,..
      'whitebg', false,... % used only for 3d slices
      ...
      'contour', false,... % draw data as a contour plot
... %      'ncontours', 5,... % not supported yet
      'contourvals', [],... % default # of contour values
      'ylog', false,...  % draw data with a logarithmic scale in x- or y-axis
      'xlog', false,...  % 
      'ydir', 'normal', ... % y-axis ascending ('normal') or descending ('reverse')
      'shading', 'flat',...% cf. matlab's 'shading' command. 'flat', 'interp'
      'ylim', [],... % ylim of final plot
      ...
      'forcesurf', false,... % use surf even when image would be faster
      'forceRGB', false,... % convert indexed to RGB using cmap (requires cmap!)
      'RGBchan', [],... % force mapping of 2D data into R,G,B (1,2,3) 
          ...           % note: 3D data will always be mapped into RGB
      'cmap', [],... % colormap to use (will change figure colormap
          ...        % unless forceRGB)
      'cmapwin',[],... % range of cmap to use for plot (default: use whole cmap)
      'cmapwin_whitebg',[],... % range of cmap to use for plot if whitebg
      'datalims',[],... % data range to scale to cmapwin (default use data min/max)
      'noclip', false,... % don't clip data that may fall outside cmapwin
          ...
      'nodatastring', 'no data to plot');
  
  
  
  a = parseArgsLite(varargin,def);

  if ndims(a.data) > 3, 
    error('Data can''t be more than 3 dimensions');
  end
  
  if ~isempty(a.RGBchan) && (a.RGBchan > 3 || a.RGBchan < 1),
    error('RGBdata option must be between 1 and 3, if given');
  end
  
  % convert indexed values to their RGB equivalent from a cmap (helpful
  % for figures with multiple colormaps)b
  if a.forceRGB,
    if isempty(a.cmap),
      error('''cmap'' must be provided for forceRGB');
    end

    % scale the data to the length of the colormap
    if isempty(a.datalims),
      valid_data = ~isinf(a.data(:)) & ~isnan(a.data(:));
      dl = [min(a.data(valid_data)) max(a.data(valid_data))];
    else
      dl = a.datalims;
    end

    % clip data to new datalims
    a.data(a.data > dl(2)) = dl(2);
    a.data(a.data < dl(1)) = dl(1);

    % put data on 0-1
    a.data = (a.data-dl(1)) ./ diff(dl);
    
    % put data on 1-length(cmap)
    cmaplen = size(a.cmap,1);
    a.data = round((a.data.*(cmaplen-1))+1); % (0->cmaplen-1)+1 = 1->cmaplen
    
    a.data = ind2rgb(a.data,a.cmap);
    a.cmap = [];
  end
  
  
  [nrows ncols nslices] = size(a.data);
  
  zlims = a.datalims;
  
  if nslices == 1 && ~isempty(a.RGBchan),
    % if we are forcing RGBchan of 1D data, expand the array to m x n x 3
    a.data(:,:,[2 3]) = 0;
    nslices = 3;
    
    % allow user to put 2D data in red, green, or blue slice
    a.data = circshift(a.data, [0 0 -(a.RGBchan-1)]);
  end
  
  if ~isempty(a.data),
    % can only provide one each of histedges/histctrs/histendctrs
    switch sum([~isempty(a.xhistedges)...
                ~isempty(a.xhistctrs)...
                ~isempty(a.xhistendedges)...
                ~isempty(a.xhistendctrs)]),
     case 0,
      warning(['none of x{histedges,histctrs,histendedges,histendctrs} provided, using ' ...
               'matrix size']);
      a.xhistctrs = 1:ncols;
     case 1,
      % ok!
     otherwise
      error('only one of x{histedges,histctrs,histendedgeshistendctrs} may be provided');
    end

    switch sum([~isempty(a.yhistedges)...
                ~isempty(a.yhistctrs)...
                ~isempty(a.yhistendedges)...
                ~isempty(a.yhistendctrs)]),
     case 0,
      warning(['none of y{histedges,histctrs,histendedges,histendctrs} provided, using ' ...
               'matrix size']);
      a.yhistctrs = 1:nrows;
     case 1,
      % ok!
     otherwise
      error('only one of y{histedges,histctrs,histendedges,histendctrs} may be provided');
    end
  end
  
  % get axes
  if isempty(a.ax),
    a.ax = gca;
    
  end
  
  % get axis size
  oldunits = get(a.ax,'units');
  set(a.ax,'units','pixels');
  axpospix = get(a.ax,'position');
  set(a.ax,'units', oldunits);
  
  % get axis grid state (surf messes with it) so we can reset it later
  oldgrid = get(a.ax, {'XGrid' 'XMinorGrid' 'YGrid' 'YMinorGrid'});
   
  % set/get the colormap
  if ~isempty(a.cmap),
    cmap = colormap(a.ax,a.cmap);
  else
    cmap = colormap(a.ax);
  end

  % no data, write a message to the user
  if isempty(a.data) || isvector(a.data), 
    if ~isempty(a.nodatastring),
      % draw advisory text box
      txtcol = [0.7 0.7 0.7];
      h = text (axpospix(3)/2, axpospix(4)/2, ... % centered
                a.nodatastring,...
                'parent', a.ax,...
                'units','pixels',...
                'verticalalignment','middle',...
                'horizontalalignment','center',...
                'color', txtcol,...
                'edgecolor', txtcol);
      hs = [hs; h];
      return;
    end
  end
  
  %%% 'unify' histbinedges, histbinctrs, histendctrs
  %
  % We verify above that we only have one each for x and y axes to start
  % with, and that data is not empty. Here we 'cross-convert' whatever the
  % user has given us so that we have all edges/centers/ends, etc. for
  % later operations
  
  [a.xhistedges a.xhistctrs a.xhistendedges a.xhistendctrs xuniform] = ...
      subf_getallhist(a.xhistedges(:), a.xhistctrs(:), ...
                      a.xhistendedges(:), a.xhistendctrs(:), ...
                      ncols);
  
  [a.yhistedges a.yhistctrs a.yhistendedges a.yhistendctrs yuniform] = ...
      subf_getallhist(a.yhistedges(:), a.yhistctrs(:), ...
                      a.yhistendedges(:), a.yhistendctrs(:), ...
                      nrows);
  
  if size(a.data,1) ~= length(a.yhistctrs) ||...
        size(a.data,2) ~= length(a.xhistctrs),
    error(['size of data not consistent with size of histedges or ' ...
           'histctrs']);
  end
  
  if ~xuniform || ~yuniform,
    %    error('non-uniform image bins not yet supported');
    warning('non-uniform image bins plotted as if uniform !!');
  end

  % Deal with 3-D data case (plot up to 3 slices in R,G,B channels at
  % each pixel)
  if nslices > 1,
    if nslices > 3, 
      error('can only have 3 ''slices'' i.e.: size(a.data,3) must be < 3');
    end
  
    if a.contour || ~strcmpi(a.shading, 'flat') || a.ylog || a.xlog
      error(['for 3-D data, contour plots, non-flat shading and log scales ' ...
             'are not supported']);
    end

    if islogical(a.data)
      % kind of a hack to get color mixtures to work with 3-D logical arrays,
      % without writing special-case code. Works fine, though R/G/B by
      % themselves are undersaturated.
      a.data = double(a.data);
      zlims = repmat([0 3], nslices, 1);
      a.zlimscale = false;
    end
    
    if isempty(zlims),
      zlims = repmat([0 1],nslices,1);
    end
    
    % if only one set of zlims provided, use it for all slices
    if size(zlims,1) == 1,
      zlims = repmat(zlims,nslices,1);
    end
    
    if a.zlimscale,
      % scale data in each slice
      if a.whitebg,
        % sum of values must be < 1, so that we can attribute the
        % 'remaining' luminance to the other pixels 
        zlims = repmat([0 max(max(sum(a.data,3)))], nslices, 1);
      else
        zlims = repmat([0 max(a.data(:))], nslices,1);
      end
    end
    
    % scale z data from zlim -> [0,1] so that we can map it into RGB
    % pixel intensities...
    for z = 1:nslices,
      zlim = zlims(z,:);
      a.data(:,:,z) = (a.data(:,:,z) - zlim(1)) ./ zlim(2);
    end
    % ...clip to 0->1
    a.data(a.data > 1) = 1;
    a.data(a.data < 0) = 0;
    
    % if requested, renormalize it so each column sums to 1
    if a.zcolnorm,
      colnorm = sum(sum(a.data,1),3);
      for k = 1:nrows,
        for l = 1:nslices;
          a.data(k,:,l) = a.data(k,:,l) ./ colnorm;
        end
      end
    end
    % otherwise NaNs get drawn as black, which is ugly with whitebg
    a.data(isnan(a.data)) = 0;

    if nslices == 3,
      im = a.data;
    else
      if a.whitebg,
        im = ones(nrows, ncols, 3); % truecolor array
        im(:,:,[2 3]) = im(:,:,[2 3]) - repmat(a.data(:,:,1),[1 1 2]); % red
        im(:,:,[1 2]) = im(:,:,[1 2]) - repmat(a.data(:,:,2),[1 1 2]); % blue
      else
        %blackbg
        im = zeros(nrows,ncols, 3);
        im(:,:,1) = a.data(:,:,1); % red
        im(:,:,3) = a.data(:,:,2); % blue
      end
    end
    
    % occasionally we get small outliers, due to rounding error. Then
    % 'image' barfs. Ruegggh! ... kah kah keh.
    im(im>1) = 1; 
    im(im<0) = 0;
    
    h = image(a.xhistendctrs, a.yhistendctrs, im,...
              'parent',a.ax);
    hs = [hs; h];
    
  else
    if isempty(a.datalims),
      valid_data = ~isinf(a.data(:)) & ~isnan(a.data(:));
      dl = [min(a.data(valid_data)) max(a.data(valid_data))];
    else
      dl = a.datalims;
    end

    
    % set/get the colormap range to use
    if ~isempty(a.cmapwin),
      if a.whitebg && ~isempty(a.cmapwin_whitebg)
        cmapwin = a.cmapwin_whitebg;
      else
        cmapwin = a.cmapwin;
      end
    else
      cmapwin = [1 size(cmap,1)];
    end
    
    % scale data so requested data lims fill colormap window
    % get appropriate 'clim' values 
    clim = newclim(cmapwin(1),cmapwin(2),...
                   dl(1),dl(2), size(cmap,1));
    
    % clip data to requested range (to avoid artifacts in neighboring
    % areas of the cmap)
    if ~a.noclip
      a.data(a.data > dl(2)) = dl(2);
      a.data(a.data < dl(1)) = dl(1);
    end
    
    if ~diff(clim)
      warning('clim values are non-increasing, setting second value to Inf');
      clim(2) = Inf;
    end
    
    set(a.ax,'clim',clim);
    
    if a.contour,
      % draw just like surf with 'interp', using centers
      [contours h] = contour(a.ax,...
                          a.xhistctrs,... 
                          a.yhistctrs,...
                          a.data,...
                          a.contourvals); %#ok dummy var     
      hs = [hs; h];
      
    else
      
      switch a.shading,
        
       case 'flat'
        
        if ~a.ylog && ~a.xlog && ~a.forcesurf,
          % we can use the much faster 'imagesc', if we don't need log y-axis, or interp
          % shading.

          % image likes to reset some properties, maybe a built-in? Undo these
          oldaxparams = {'visible', 'box', 'xdir', 'ydir'};
          oldaxvals = get(a.ax, oldaxparams);

          h = imagesc(a.xhistendctrs, a.yhistendctrs, a.data,...
                      'parent',a.ax,...
                      clim);

          set(a.ax, oldaxparams, oldaxvals);
          
          hs = [hs; h];
          
        else
          % to draw an image with a log y-axis scale, we have to use 'surf'

          % 'flat' shading draws the block using the color of the corner with the
          % smallest indices (i.e. the bottom-left corner; cf. 'help
          % shading'). We therefore want to shift the grid by 1/2 an image bin
          % down and to the left so that blocks are *centered* on the requested
          % data.
          %
          % We can do this by using the {x,y}histedges instead of 'ctrs',
          % but we need to add a 'dummy' row and column at the top and right
          % of the array so that the last data point will be drawn
          % correctly. (note this dummy data will not be used for plotting
          % with 'flat' shading);
          
          % add in dummy data (zeros)
          padded_data = zeros(size(a.data) + [1 1]); % faster than 2 appends (maybe?)
          padded_data(1:end-1,1:end-1) = a.data;

          % surf will gracefully handle negative data, with a 'Warning:
          % Negative data ignored', so we don't need to test for it  
          
% $$$           if (a.xlog && any(a.xhistedges <= 0)) || ...
% $$$                 (a.ylog && any(a.yhistedges <= 0)),
% $$$             warning('Can''t plot values at 0 or below on a log scale, clipping');
          
          h = surf(a.xhistedges,...
                   a.yhistedges,...
                   zeros(size(padded_data)),...
                   padded_data,...
                   'parent',a.ax,...
                   'facecolor', 'flat',...
                   'edgecolor', 'none');       
          
          hs = [hs; h];

        end
        
       case 'interp'
        
        % have to use 'surf' plot to get interp shading, whether or not we use the
        % log scale. Bin centers are appropriate in this case, since interp
        % plots the appropriate color *at* the requested data point.
        
        h = surf(a.xhistctrs,...
                 a.yhistctrs,...
                 a.data,...
                 'parent',a.ax,...
                 'facecolor', 'interp',...
                 'edgecolor', 'none');
        
        hs = [hs; h];
        
      end

    end % if a.contour
  end
  
  % required for surfs and for displaying of images        
  set(a.ax,'View',[0 90]);

  % imagesc defaults to 'ydir' = reverse
  set(a.ax, 'ydir', a.ydir);
  
  % set log scaling in y, if requested
  if a.ylog,
    set(a.ax,'yscale','log');
  else
    set(a.ax,'yscale','linear');
  end

  % set log scaling in y, if requested
  if a.xlog,
    set(a.ax,'xscale','log');
  else
    set(a.ax,'xscale','linear');
  end

  
  % set lims. xlim tight, ylim either tight or user-specified
  if ~isempty(a.ylim) && ~diff(a.ylim) == 0
    ylim(a.ax, a.ylim);
  end
  
  % restore grid state
  set(a.ax, {'XGrid' 'XMinorGrid' 'YGrid' 'YMinorGrid'}, oldgrid);
  
  
function [histedges histctrs histendedges histendctrs uniform] = ...
      subf_getallhist(histedges, histctrs, histendedges, histendctrs, ...
                      datasz)
  % given any of a number of ways of specifying the edge vectors of an
  % array, interconvert them all sensibly.
  
  %% redundant, but what the hell
  if sum([~isempty(histedges)...
          ~isempty(histctrs)...
          ~isempty(histendedges)...
          ~isempty(histendctrs)]) ~= 1,
    error('only one of edges/ctrs/endctrs may be provided');
  end
    
  % test for uniformity before mungeing inputs
  uniform = false;
  if ~isempty(histendctrs) || ~isempty(histendedges),...
    % when given first/last only, assume uniformity
    uniform = true;
  elseif size(histctrs,1) == 2 || size(histedges,1) == 2,
    uniform = true;
  elseif ~isempty(histctrs) &&...
        all(abs(diff(histctrs,2)) <= 4*eps(max(histctrs)));
    uniform = true;
  elseif ~isempty(histedges) &&...
        all(abs(diff(histedges,2)) <= 4*eps(max(histedges)));
    uniform = true;
  end
  
  if ~isempty(histendctrs)
    histctrs = (linspace(histendctrs(1), histendctrs(2), datasz))';
  end
  
  if ~isempty(histendedges)
    histedges = (linspace(histendedges(1), histendedges(2), datasz+1))';
  end
  
  if isempty(histedges)
    % if histedges is empty, then we def have histctrs by now

    % Assume that edges are at the mean of the neighboring bins. Also assume
    % that first bin is as wide as the second, and that the last is as wide
    % as the second-to-last. (works for uniform/nonuniform)
    ctrsoffset = diff(histctrs)./2;
    histedges = histctrs - reshape((ctrsoffset([1 1:end])),[],1); %column-ize      
    histedges = [histedges; histctrs(end)+ctrsoffset(end)];

  else
    % we have histedges, and therefore not histctrs yet
    
    histctrs = histedges(1:end-1) + (diff(histedges)./2);
  
  end
  
  % we have histctrs/histedges now; get histendctrs, histendedges
  histendctrs = histctrs([1 end]);
  histendedges = histedges([1 end]);
  