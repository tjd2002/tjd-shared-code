function [fh avih] = slp2dmovie(varargin)
% TODO
%  - base it on sla.plots
%  - rat's position as triangle?
%  - are we off by 0.5 bins? should we use start/middle of bins?
%  - display linearized position as grey dot on track? or as blue/red dot
%  depending on which direction he's going?
%  - framerate (for eeg, rat position updates) different from posest
%  binsize. I.e. don't update track every time, but keep movie smooth
%  - second/third (nmaps) tracks for reverse (in blue?)
%  - more generally: replot spikeslpos, but in time. Have spikeslpos just
%  spit out the pfmaps, slraster(no?), and contdata, and use as is? (sld
%  already has fpfmap, etc.), multiple contdatas. See sltodo
%  - use patches instead of lines, get rid of gaps
%  - resizfn that changes size of eeg trace, thickness of track line
%
% LATER/MAYBE
%  - DVD frame import/alignment?
%  - output to avi?/DVD (im2frame, mpgwrite,...)

a = struct(...
    'argstruct',[],...
    'e',{{}},...
    'movfigh',[],...
    'pfmaps',[],...
    'autoscale',[],...
    'scale',[],...
    'sharpen',[],...
    'distbinsize',[],...
    'lintrack',false,...
    'trackwidth',5,...
    'cmap',hot,...
    'contdata',[],...
    'contwindow',[-0.5 0.5],...
    'timebinsize',[],...
    'framerate',[],...
    'timewin',[],...
    'avifilename',[]);


a = parseArgs(varargin,a);


if isstruct(a.argstruct),
  
  matchi = strcmpi(varargin, 'argstruct');
  if sum(matchi) > 1, error([mfilename ':BadArgs'], 'Only one ''argstruct'' may be provided'); end
  if matchi(1) ~= 1, error([mfilename ':BadArgs'], '''argstruct'' must be first argument, if provided'); end

  % use argstruct as defaults, let additional params override.
  a = a.argstruct;
  a = parseArgs(varargin(3:end),a);

  % don't allow for recursive argstructs
  if ~isempty(a.argstruct);
    warning([mfilename ':BadArgs'],'Recursive argstructs not allowed--ignoring');
    a.argstruct = [];
  end

else
  if any(strcmpi(varargin, 'argstruct')),
    error([mfilename ':BadArgs'], '''argstruct'' does not specify a struct');
  end
end


if ~isempty(a.movfigh) && ishandle(a.movfigh) && strcmp(get(a.movfigh,'type'),'figure'),
  fh = a.movfigh;
else
  %bring to front
  fh = gcf;

  % don't accidentally overwrite a spikeslpos figure
  if ~isempty(getappdata(fh,'slargs')),
    fh = figure;
  end
end

set(0,'currentfigure',fh);

% if no framerate use realtime
if isempty(a.framerate),
  a.framerate = 1/a.timebinsize;
end


% are we making a movie file?
makeavi = false;
avih = [];
if ~isempty(a.avifilename),
  % apparently faster
  set(fh,'visible','off');
  makeavi = true;
  avih = avifile(a.avifilename,...
                 'FPS',a.framerate);
end




%%% figure stuff

%no menu bar
set (fh, 'menubar', 'none');

clf(fh);

%under opengl this is the trackmap bg color, too.
set (fh, 'color', [0.2 0.2 0.2]);
set (fh, 'color', 'k');


%set (fh, 'renderer', 'painters');
% seems to make no diff;
%set (fh, 'doublebuffer', 'on');
% hardware opengl:wicked fast, but glitchy
%set (fh, 'renderer', 'opengl');
set (fh, 'renderer', 'zbuffer');
%set (fh, 'renderer', 'painters');

% $$$ clf;
%%% setup subaxes

if ~isempty(a.contdata),
  if a.lintrack,
    trackax = subplot(3,1,1);
    contax = subplot(3,1,2:3);
  else
    trackax = subplot(4,1,1:3);
    contax = subplot(4,1,4);
  end
  axlist = [trackax contax];
else
  trackax = gca;
  axlist = trackax;
end

for ax = axlist;
%  set (ax, 'layer', 'top');
  set (ax, 'xcolor', 'w');
  set (ax, 'ycolor', 'w');
end

% setup track movie window

if ~a.lintrack,
  % equivalent to axis ij, use video co-ordinates
  set(trackax,'ydir','reverse');
  set(trackax,'dataaspectratio',[1 1 1]);
end

set(trackax,'dataaspectratio',[1 1 1]);
%set(trackax,'tickdir','out');
set(trackax, 'color', [0.2 0.2 0.2]);
set(trackax, 'box', 'on');



%%% select posdata from timewin
posflds = { 'time', 'centx', 'centy', 'headdir' 'x1' 'y1' 'x2' 'y2' ...
            'l' 'ldir' 'offtrackd'};
posdat = getepd(a.e.pos, posflds{:});

startpos = find(posdat(:,1) > a.timewin(1),1,'first');
endpos = find(posdat(:,1) < a.timewin(2),1,'last');

posdat = posdat(startpos:endpos,:);


%%% select contdata from timewin

if ~isempty(a.contdata),
  
  set(contax, 'color', 'k');
  
  set(contax,'units','pixels');
  axwidpix = get(contax, 'position');
  axwidpix = axwidpix(3);
  
  contstart = a.timewin(1) + a.contwindow(1); % this is -ve
  contend = a.timewin(2) + a.contwindow(2);
  contsamplerate = size(a.contdata,1) / (a.contdata(end,1) - ...
                                         a.contdata(1,1)); %#ok
  
  
  contwin(1) = max(1, floor((contstart - a.contdata(1,1)) * ...
                            contsamplerate));
  
  contwin(2) = min(size(a.contdata,1), ...
                   ceil((contend - a.contdata(1,1)) ...
                        * contsamplerate));

  % plot it all
  if all(contwin) && (contwin(1) < contwin(2)),

    
    % subsample to speed plotting
    sampsperpix = 5;
    
    nsamps = contsamplerate * diff(a.contwindow);

    if nsamps < sampsperpix*axwidpix
      % plot all, no subsample
      step = 1;
    else
      % step
      step = floor(nsamps/axwidpix/sampsperpix);
    end
    
    line(a.contdata(contwin(1):step:contwin(2),1), ...
         a.contdata(contwin(1):step:contwin(2),2), ...
         'parent',contax,...
         'color', [1 0 0]);

    % constant xticks scrolling by every 2 seconds
    contxticks = round(a.contdata(contwin(1),1)):...
        max([1 round(diff(a.contwindow)/2)]):...
        round(a.contdata(contwin(2),1));
    
  end

  % we'll scroll with xlim, want ylim to stay constant
  set(contax,'xlimmode','manual');
%  set(contax,'ylimmode','manual');
  ylim ([min(a.contdata(contwin(1):contwin(2),2)) ...
         max(a.contdata(contwin(1):contwin(2),2))]);
    
end

if a.lintrack, %track is going to be at x=l; y=0
  bins_lxy = [0:a.distbinsize:a.e.track.length a.e.track.length]';
  bins_lxy = [bins_lxy bins_lxy zeros(length(bins_lxy),1)];
  
   
else % 2d track:
  bins_lxy = [0:a.distbinsize:a.e.track.length a.e.track.length]';
  bins_lxy = [bins_lxy interp1q(a.e.track.l, a.e.track.splinevals, ...
                                bins_lxy)];
end

xlim(trackax, [min(bins_lxy(:,2))-20 max(bins_lxy(:,2))+20]);
ylim(trackax, [min(bins_lxy(:,3))-20 max(bins_lxy(:,3))+20]);



nsegs = size(bins_lxy, 1)-1;
nframes = size(a.pfmaps, 2);

lineh = zeros(1,nsegs);
rath = [];
eegnowlineh = [];


drawtime = zeros(1,nframes);
frametime = zeros(1,nframes);

if ~isempty(a.autoscale) && a.autoscale > 0,
  pfmax = prctile(a.pfmaps(:),a.autoscale*100);
else
%  pfmax = max(a.pfmaps(:));
  pfmax = 1;
end

pfmin = 0; % for now, could select a range if we wanted

%compute colormap index of each pt (then don't have to do it in the loop)

colscale = (length(a.cmap)-1)/(pfmax-pfmin);

%pfmapsi = zeros(size(pfmaps)); % 65K values, big enough!
pfmapsi = round((a.pfmaps-pfmin).*colscale) + 1;

% get rid of values out of cmap range (<1, >cmap length)
pfmapsi(pfmapsi<1) = 1;
pfmapsi(pfmapsi>length(a.cmap)) = length(a.cmap);

% bring fig to front (do it now to avoid flickering from clf, etc);
set(fh,'windowstyle','modal');

t_movstart = clock;
t_framestart = clock;
for j = 1:nframes,

  curtime = a.timewin(1) + ((j-0.5)*a.timebinsize);

  xlabel(trackax, num2str(curtime));
  
  for k = 1:nsegs

    % first run, init plot:
    if j == 1,
      
      % plot at z=-1 so that other stuff shows up on top
      lineh(k) = line(bins_lxy([k k+1],2), ...
                      bins_lxy([k k+1],3), ...
                      [-1 -1],...
                      'parent', trackax,...
                      'linewidth', dist2pts(a.trackwidth,trackax,'Y'), ...
                      'erasemode', 'normal');


    end
    % update color of line from pfmapsi (eventually use indexed color
    % with 'patch')
    set(lineh(k),'color',a.cmap(pfmapsi(k,j),:));
  end


  
  %%%  plottez le rat!

  posi = find((posdat(:,1) > a.timewin(1) + j*a.timebinsize) & ...
              (posdat(:,1) < a.timewin(1) + (j+1)*a.timebinsize));
  
  if ~isempty(posi),
    % for now just draw one position per frame
    posi = posi(1);
    if isempty(rath),
      rath(1) = line('EraseMode','normal',...
                     'Parent',trackax,...
                     'Color','w',...
                     'linewidth', 0.5, ...
                     'LineStyle','-',...
                     'MarkerEdgeColor','none',...
                     'MarkerFaceColor','g',...
                     'Marker','o',...
                     'MarkerSize',5);

      if ~a.lintrack,
        rath(2) = line('EraseMode','normal',...
                       'Parent',trackax,...
                       'LineStyle','none',...
                       'MarkerEdgeColor','none',...
                       'MarkerFaceColor','r',...
                       'Marker','o',...
                       'MarkerSize',5);
      end
    end
    if a.lintrack, % single point at l(i),offtrackd(i) for now
      set(rath(1),'xdata',posdat(posi,[9]),...
                  'ydata',posdat(posi,[11]),...
                  'zdata',1); %#ok 
    
    else % line between diode1/diode2
      set(rath(1),'xdata',posdat(posi,[5 7]),...
                  'ydata',posdat(posi,[6 8]),...
                  'zdata',[1 1]);
      set(rath(2),'xdata',posdat(posi,[5]),... 
                  'ydata',posdat(posi,[6]),...
                  'zdata',[2]); %#ok
    end
  end
  
  %%% scroll the eeg
  if ~isempty(a.contdata),
    if isempty(eegnowlineh)
      eegnowlineh = line([curtime curtime], ylim, 'color','w');
    else
      set(eegnowlineh,'xdata',[curtime curtime]);
    end
    xlim(contax,curtime + a.contwindow);
    set (contax,'xtick', sort([contxticks curtime]));
  end
  
  if makeavi,
    avih = addframe(avih,fh);
  end
  
  %%% GETFRAME IS DOG SLOW FOR NO GOOD REASON!!!
  % im2frame is lightning fast, but I'll have to draw into it
  % myself (as in slrasters).
  % -addframe now 10fps after getting hardware openGl working. Problem
  % with getframe may be that it is using 'noanimate', and 'addframe',
  % and switching the renderer to opengl/zbuffer and back every on every
  % frame.
  %
  % cf. matlab code: im2frame ( a builtin); hardcopy; noanimate...
  %
  % ***use screengrabber!***, make system call using:
  % !import -window ? filename.jpeg
  % !scrot filename.png
  % xwininfo/xwd/xwud/GIMP (fast, but only into xwd/ImageMagick format)
  % fbgrab (framebuffer only, faster??)
  % ... screenshot utilities
  % then drawnow() 
  % 
   drawnow;

  drawtime(j) = etime(clock,t_framestart);
  
  % calculate deviation from expected time, so errors don't accumulate
  % (means we try to catch up--is this good?)
  %wait_t = j/a.framerate - etime(clock,t_movstart);
  wait_t = 1/a.framerate - etime(clock,t_framestart);
  pause(wait_t);
  % time between frames, including drawing and waiting
  frametime(j) = etime(clock,t_framestart);
  t_framestart = clock;

  if makeavi
    disp(j)
  end
% $$$   if wait_t > 0,
% $$$     pause(wait_t);
% $$$   end
  
end
disp(prctile(drawtime,[0 10 50 90 100]));
disp(prctile(frametime,[0 10 50 90 100]));

if makeavi,
  avih = close(avih);
  set(fh,'visible','on');
end

setappdata(fh,'slmovargs',a);
setappdata(fh,'slmovdata',[]);
set(fh,'keypressfcn', 'slmovkeypress(gcbf)');
set(fh,'windowstyle','normal');
