function tf = isfieldmult(s,fields)
%ISFIELDMULT array of logicals if field is in structure array.
%   F = ISFIELDMULT(S,'field'/{fields}) returns true if 'field' is the name of a field
%   in the structure array S. 

if ~iscell(fields),
  fields = {fields};
end

%preallocate tf
tf = false(size(fields));

if isa(s,'struct')
  fnames = fieldnames(s);
  for k = 1:length(fields),
    f = fields{k};
    tf(k) = any(strcmp(fnames,f));
  end
end

