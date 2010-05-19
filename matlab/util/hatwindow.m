function h1 = hatwindow(s,W)
% HATWINDOW create a 1-D 'mexican-hat' window
%
% s = sigma, very roughly the width of the central peak at 1/2 height
% W = width,in units of sigma at which to truncate
%
% brazenly stolen without understanding from:
% http://www.cs.ubc.ca/~woodham/cpsc505/examples/log.html

  if ~exist('W','var') || isempty(W),
    W = 4;
  end
  
x = -s*W : 1 : s*W;
n = length(x);

ss = s*s;

h1 = zeros(n,1);
h2 = h1;

for k = 1:n
 t = x(k);
 a = (t*t)/ss;
 b = exp (- a / 2);
 h1(k) = (1 - a) * b;
 h2(k) = b;
end

% Force h1 coefficients to sum to zero
h1 = h1 - sum(h1(:))/prod(size(h1));

