function varargout = inseg( A,B,varargin )
%INSEG test whether events/segments are contained in a set of reference segments
%
%  b=INSEG(s,t) where t can be a vector of events or a list of
%  segments. This function returns true for each event/segment t that lies
%  within at least one of the segments in s and false otherwise.
%
%  [b,n]=INSEG(s,t) also returns for each segment s the number of
%  segments/events in t that it contains.
%
%  [b,n,rnk]=INSEG(s,t) also returns for each event/segment in t the
%  relative order in the containing segment. A value of 0 indicates that
%  event/segment t is not contained in any reference segment s. If
%  segments in s overlap, it is possible that events/segments in t are
%  contained in multiple segments s. In that case the values in rnk are
%  invalid and a warning will be given.
%
%  [b,n,rnk,si]=INSEG(s,t) also returns for each event/segment in t the
%  index of the first segment s that contains that event/segment (0 if no
%  segment contains the event/segment).
%
%  [u,i,rnk]=INSEG(...,'expand') returns cell arrays with for
%  each segment in s the events/segments in t that it contains (u),
%  their index (i) and their relative order (rnk).
%
%  [...]=INSEG(...,'partial') only when t is a segment list, this also
%  returns true for segment that only partially overlap with any of the
%  reference segments in s.
%

%  Copyright 2005-2008 Fabian Kloosterman

%two input arguments required
if nargin<2
  help(mfilename)
  return
end

%check whether B is a list of events or segments
if (isvector(B) || isempty(B)) && size(B,2)~=2
  Bseg = false;
  B=B(:);   %make events a column vector
else
  Bseg = true;
end

%check for valid positive length (incl. 0 ) segments
if ( ~isnumeric(A) || ndims(A)~=2 || size(A,2)~=2 || any(diff(A,1,2)<0) ) || ...
      (Bseg && ( ~isnumeric(B) || ndims(B)~=2 || size(B,2)~=2 || any(diff(B,1,2)<0))),    
  error('seginseg:invalidArgument', 'Invalid segments')
end

%initialization
nA = size(A,1);    %number of segments
nB = size(B,1);    %number of events/segments
C = cell(nA,1);    %for storing events/segments
Cidx = cell(nA,1); %for storing indices
n = zeros(nA,1);   %for storing number of events/segments contained in
                   %reference segment
nthB = cell(nA,1); %for storing relative order statistic

option = false;
partial = false;

if ismember( 'expand', varargin )
  option = true;
end
if ismember( 'partial', varargin )
  partial = true;
end

containedBy = zeros( nB,1 );

%loop through all segments in A
for k = nA:-1:1
  
  %find all events/segments contained in reference segment
  if partial
    %partial overlap
    idx = find( B(:,1)<=A(k,2) & B(:,end)>=A(k,1) );    
  else
    %complete overlap
    idx = find( B(:,1)>=A(k,1) & B(:,end)<=A(k,2) );
  end

  C{k} = B(idx,:); %save events/segments
  Cidx{k} = idx;   %save indices

  n(k) = length(idx); %save number of events/segments found

  %compute relative order statistic
  [dummy, nthB{k}] = sort( B(idx,1) ); %#ok

  if ~option
    %get index of first A containing each B
    containedBy(idx) = k;
  end  
end

%construct output
if option
  varargout = {C, Cidx, nthB};
else
  
  [Cidx, tmp] = unique( cat( 1, Cidx{:} ) );
  varargout{1} = false(nB,1);
  varargout{1}(Cidx)=true;
  
  varargout{2} = n;
  
  varargout{3} = zeros(nB,1);
  nthB=cat(1,nthB{:});
  varargout{3}(Cidx)=nthB(tmp);
   
  varargout{4} = containedBy;
end
