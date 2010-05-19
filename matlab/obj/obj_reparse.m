function obj = obj_reparse(obj, arglist)
% OBJ_REPARSE parse 'template' args, do some other setup, for OBJs
  
  obj = parseArgsLite(arglist,obj);
  
  % parse 'template'
  if ~isempty(obj.template),
    
    % copy input args from obj.template to input arg list
    for field = obj.inputargs,
      field = field{:};
      obj.(field) = obj.template.(field);
    end
    
    % update with additional args 
    obj = parseArgs(arglist,obj);
    
  end
  
  % ensure we populate edesc
  if isfield(obj, 'e') && ~isempty(obj.e),
    obj.edesc = obj.e.desc;
  end

