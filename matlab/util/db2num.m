function num = db2num(db)
% DB2NUM convert value in decibels to fraction
%
%  see also: db  (signal toolbox)
  
  num = 10.^(db./20);
  