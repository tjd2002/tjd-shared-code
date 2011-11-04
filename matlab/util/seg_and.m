function seg = seg_and(varargin)
%SEG_AND logical AND on segment lists
%
%  seg=SEG_AND(seg1,seg2,...) performs AND operation on segments,
%  returning the overlap between segments.
%

%  Copyright 2005-2008 Fabian Kloosterman

%check input arguments
narg = length(varargin);
if (narg<2)
    help(mfilename)
    return
end

s = {};

for i = 1:narg
    if ~isempty(varargin{i}) && (size(varargin{i}, 2)~=2) % | (size(varargin{i}, 1)<1)
        error('seg_and:invalidArgument', ['Expecting a nx2 matrix of segment start and end times (n>0). Error in argument ' num2str(i)])
    else
        s = [s ; varargin(i)];
    end
end


%take the first segment list
seg_stack = s{1};

%and find overlap with the other segments lists
for i = 2:length(s)

    overlap = zeros(0, 2);
    
    for j = 1:size(seg_stack, 1)
        %find overlapping segments
        ind = find( (seg_stack(j,1)<=s{i}(:,1) & seg_stack(j,2)>s{i}(:,1)) | (s{i}(:,1)<=seg_stack(j,1) & s{i}(:,2)>seg_stack(j,1)) );
        if ~isempty(ind)
            overlap = [overlap; [max(seg_stack(j,1), s{i}(ind,1)) min(seg_stack(j,2), s{i}(ind,2))]];
        end
    end
    
    if isempty(overlap)
        seg = [];
        return
    end
    
    seg_stack = overlap;
end
       
seg = overlap;
