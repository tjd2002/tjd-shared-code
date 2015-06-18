function equiv = structcmp(s1,s2,fields)
% STRUCTCMP compare structs (numeric, string, logical, cell arrays of
% strings, function handles and struct fields accepted)
%
% if 'fields' is given, compare only those fields of the two structures. If
% either argument is not a struct, or does not contain all the fields,
% return false.
    
    equiv = true;

    %%% both structs?

    if ~isstruct(s1) || ~isstruct(s2),
      equiv = false; return;
    end
    
    
    %%% structs have same fields?    

    % call it once here, rather than many times (esp in 'isfield')
    s1fnames = fieldnames(s1);
    s2fnames = fieldnames(s2);

    if nargin <=  2,

      if length(s1fnames) ~= length(s2fnames) || ...
            ~all(strcmp(s1fnames,s2fnames)),
        equiv = false; return;
      end
      
    else
      
      % only test selected fields
      
      % convert 'fields' arg to a cell, if user provided a string
      if ischar(fields),
        fields = {fields};
      end
      
      for f = fields(:)';
        f = f{:};
        
        if ~any(strcmp(s1fnames,f)) || ~any(strcmp(s2fnames,f)),
          equiv = false; return;
        end

        news1.(f) = s1.(f);
        news2.(f) = s2.(f);
      end
      s1 = news1;
      s2 = news2;
      s1fnames = fields';
%     s2fnames = fields';
      
    end

    %%% if array of structs, same size?
    if ~all(size(s1) == size(s2)),
      equiv = false; return
    end
    
    %%%  fields have same values?

    % iterate over arrays of structs (must have same order to be equiv)
    for k = 1:numel(s1);
      s1k = s1(k);
      s2k = s2(k);
      
      for f = s1fnames'
        f = f{:};
        if ~all(size(s1k.(f)) == size(s2k.(f))),
          equiv = false; return;
        end
        
        % compare numerics before testing class to allow for equivalency
        % across numeric classes (e.g. 2 == single(2))
        if isnumeric(s1k.(f)),
          if ~isnumeric(s2k.(f)) ||...
                ~(all(s1k.(f)(:) == s2k.(f)(:))),
            equiv = false; return;
          end
        else
          
          if ~strcmp(class(s1k.(f)), class(s2k.(f))),
            equiv = false; return;
          end
          
          switch class(s1k.(f)),

           case 'struct',
            % recursive call on struct args
            if ~structcmp(s1k.(f), s2k.(f)),
              equiv = false; return;
            end
            
           case 'char',
            if ~strcmp(s1k.(f), s2k.(f))
              equiv = false; return;
            end
            
           case 'logical',
            if ~all(s1k.(f) == s2k.(f))
              equiv = false; return;
            end
            
           case 'cell',
            if iscellstr(s1k.(f)),
              if ~all(strcmp(s1k.(f), s2k.(f))),
                equiv = false; return;
              end
            else
              for l = 1:numel(s1k.(f))
                s1kfl = s1k.(f){l};
                s2kfl = s2k.(f){l};
                if ~(isnumeric(s1kfl) || islogical(s1kfl)) ||...
                      ~(isnumeric(s2kfl) || islogical(s2kfl)),
                  error(['cell arrays that are not logical, numeric or ' ...
                         'strings are not supported']);
                end
                if ~all(size(s1kfl) == size(s2kfl)),
                  equiv = false; return;
                end
                if ~all(s1kfl == s2kfl),
                  equiv = false; return;
                end
              end
              
% $$$               % cell array of non-strings: ugly hack to convert cell array to struct,
% $$$               % which we can compare recursively using structcmp
% $$$               tmp_fnames = char(ceil(rand(numel(s1k.(f)),12).*26) + 64);
% $$$               if ~structcmp(cell2struct(s1k.(f)(:),tmp_fnames,1),...
% $$$                            cell2struct(s2k.(f)(:),tmp_fnames,1))
% $$$                 equiv = false; return;
% $$$               end
            end

           case 'function_handle',
            % only returns true if it's a pointer to the same fnhandle,
            % usu because they were set to be '='
            if ~isequal(s1k.(f), s2k.(f)),
              equiv = false; return;
            end
            
           otherwise
            error ('unsupported field data type');
          end
        end
      end
    end