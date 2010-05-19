function cache = mkcache(varargin)
% MKCACHE make/update a cache of 'obj's 
  
  a = struct(...
      'cache',[],...
      'add_obj',[],...
      'limit',[]);

  a = parseArgsLite(varargin, a);

  if ~isempty(a.cache),
    cache = a.cache; 
    cache.cache = [];
  else
    cache = a;
  end
  
  % default: keep 20 of each kind of object in cache
  if ~isempty(a.limit)
    cache.limit = a.limit;
  elseif isempty(cache.limit),
    cache.limit = 20;
  end
  
  if ~isempty(a.add_obj);
    
    if a.add_obj.cache_hit,
      return;
      %      warning('adding object with ''cache_hit'' == true');
    end
    
    % cache is only 1 level deep, please
    a.add_obj.cache = []; 
    
   otype = a.add_obj.type;

   if ~isfield(cache, otype),
     cache.(otype) = [];
   end
   
   if length(cache.(otype)) > cache.limit - 1;
     cache.(otype) = [a.add_obj cache.(otype)(1:cache.limit-1)];
   else
     cache.(otype) = [a.add_obj cache.(otype)];
   end
  end
  
  cache.add_obj = [];