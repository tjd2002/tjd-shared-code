function dateStringLocal = utc_to_local_time( dateVectorUtc, format )
%
% utc_to_local_time -- convert UTC to local time/date.
%
% dateStringLocal = utc_to_local_time( dateVectorUtc ) converts a vector of UTC
%    date/time information to local dates and times.  The format of dateVectorLocal must
%    be a vector of Matlab date vectors (ie, N x 6, for example [2008 10 22 13 41 00 ;
%    2008 10 22 13 41 12]).  The time/date values will be converted to local time and
%    returned in 'dd-mmm-yyyy HH:MM:SS' format.
%
% dateStringlocal = utc_to_local_time( dateVectorUtc, format ) returns the local time / 
%    date values in a format specified by numeric argument format.  The format values
%    match those used by the Matlab datestr function.  Supported formats are as follows:
%
%    Number           String                   Example
%    ===========================================================================
%       0             'dd-mmm-yyyy HH:MM:SS'   01-Mar-2000 15:45:17 
%      13             'HH:MM:SS'               15:45:17     
%      15             'HH:MM'                  15:45        
%      21             'mmm.dd,yyyy HH:MM:SS'   Mar.01,2000 15:45:17 
%      30 (ISO 8601)  'yyyymmddTHHMMSS'        20000301T154517 
%      31             'yyyy-mm-dd HH:MM:SS'    2000-03-01 15:45:17 
%
% Version date:  2008-November-12.
%

% Modification History:
%
%    2008-November-12, PT:
%        switch to use of Java methods and eliminate conditional on unix (since Java runs
%        correctly anywhere).
%
%=========================================================================================

% check that the format number is supported

  if ( nargin == 1)
      format = 0 ;
  end
  supportedFormats = [0 13 15 21 30 31] ;
  if ~ismember(format,supportedFormats)
      error('common:utcToLocalTime:unsupportedFormat' , ...
          'utc_to_local_time:  unsupported format requested' ) ;
  end
  
% call the converter vector

  dateVectorLocal = utc_date_vector_to_local( dateVectorUtc ) ;
  
% convert the date vectors to the desired format

  dateStringLocal = datestr(dateVectorLocal,format) ;
  stringLength = size(dateStringLocal,2) ;
  nStrings = size(dateStringLocal,1) ;
    
return  
  
% and that's it!

%
%
%

%=========================================================================================

% subfunction which converts a Matlab date vector in UTC to one in local time

function dateVectorLocal = utc_date_vector_to_local( dateVectorUtc )
      
  nDates = size(dateVectorUtc,1) ;
  dateVectorLocal = zeros(nDates,6) ;
  
% import the Java classes needed for the conversion and construct an appropriate
% SimpleDateFormat object in local time zone

  import java.text.SimpleDateFormat ;
  import java.util.Date ;
  import java.util.Calendar ;
  import java.util.TimeZone ;
  
  localFormatObject = SimpleDateFormat('yyyy-MM-dd HH:mm:ss') ;
  
% construct an appropriate calendar instance in UTC

  utcCalendarObject = Calendar.getInstance(TimeZone.getTimeZone('UTC')) ;
      
% loop over date strings 

  for iDate = 1:nDates

      dateVec = dateVectorUtc(iDate,:) ;
      
%     construct a Java Date object in UTC from the UTC dateVec object.  Note that Java's
%     month is zero-based.  Note further that, to define a date in UTC, we actually need
%     to use a Calendar object instance defined above

      utcCalendarObject.set( dateVec(1), dateVec(2)-1, dateVec(3), ...
          dateVec(4), dateVec(5), dateVec(6) ) ;
      utcDateObject = utcCalendarObject.getTime() ;
      
%     construct a local time string from the UTC date object

      dateStringNew = char(localFormatObject.format(utcDateObject)) ;
          
%     pick through the resulting string and extract the data we want, converting to
%     numbers as we go

      dateVectorLocal(iDate,1) = str2num(dateStringNew(1:4)) ;
      dateVectorLocal(iDate,2) = str2num(dateStringNew(6:7)) ;
      dateVectorLocal(iDate,3) = str2num(dateStringNew(9:10)) ;
      dateVectorLocal(iDate,4) = str2num(dateStringNew(12:13)) ;
      dateVectorLocal(iDate,5) = str2num(dateStringNew(15:16)) ;
      dateVectorLocal(iDate,6) = str2num(dateStringNew(18:19)) ;

  end % loop over dates
            
return