function obj = obj_cleanup(obj)
% OBJ_CLEANUP
  
% just save e.desc (for comparison), not whole exp struct
  if all(isfieldmult(obj, {'e' 'edesc'})) && ...
        ~isempty(obj.e),
    
    obj.edesc = obj.e.desc;
    obj.e = [];
  end

  % just save contdata.name (for comparison), not whole exp struct
  if all(isfieldmult(obj, {'contdata' 'contvar' 'contname'})) && ...
        ~isempty(obj.contdata),
    obj.contname = obj.contdata.name;
    
    % don't remove contdata from 'cont' obj, but do remove it from,
    % e.g. specgramdata
    if ~strcmp(obj.type, 'cont'),
      obj.contdata = [];
    end
  end
  
  obj.cache = [];
  obj.template = [];
