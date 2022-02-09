;+
; NAME:
;    convert_time_format
;
; PURPOSE:
;    This function converts a time vector from one format to another.
;
; CATEGORY:
;    Calendar
;
; CALLING SEQUENCE:
;    time_out = convert_time_format( time_in, format_in, format_out )
;
; INPUTS:
;    TIME_IN:  The required time vector, in units of "format_in", to be 
;        converted.  Of length N_TIME.  This can be a vector of strings, 
;        integers, or (single or double) floating point numbers, depending on 
;        the format.  For efficiency, values are assumed to be in increasing 
;        order.
;    FORMAT_IN:  The required description of the format of time values in 
;        TIME_IN.  Of type string.  Possible values are:
;        * '<units> since <origin>':  Where "<units>" is 'days', 
;          and "<origin>" is of the format "yyyy-mm-dd", "yyyy-mm-dd-hh", 
;          "yyyy-mm-dd-hh-mm", or "yyyy-mm-dd-hh-mm-ss".  The "-" the sub-daily 
;          components can be a ":" instead, and the "-" preceding the "hh" can 
;          be a space.  This is commonly used in NetCDF files.  'months' and 
;          'years' units are not fully supported yet.
;        * 'decimal year':  For example, [2000.04167,2000.1250].
;        * 'yyyy', 'yyyymm', 'yyyymmdd', 'yyyymmddhh', 'yyyymmddhhmm', or 
;          'yyyymmddhhmmss'.  This can be input either as a string or an 
;          integer (which will be converted to a string below).
;    FORMAT_OUT:  The required description of the format of time values to be 
;        output in TIME_OUT.  Of type string.  Supported values are:
;        * 'days since <origin>':  Where "<origin>" is of the format 
;          "yyyy-mm-dd", "yyyy-mm-dd-hh", "yyyy-mm-dd-hh-mm", or 
;          "yyyy-mm-dd-hh-mm-ss".
;        * 'decimal year':  For example, [2000.04167,2000.1250].
;        * 'yyyy', 'yyyymm', 'yyyymmdd', 'yyyymmddhh', 'yyyymmddhhmm', or 
;          'yyyymmddhhmmss'.  This can be input either as a string or an 
;          integer (which will be converted to a string below).
;
; KEYWORD PARAMETERS:
;    CALENDAR:  The optional name of the calendar to use.  Supported values are 
;        '360_day', '365_day' (or 'noleap' ), 'gregorian'.  The default is 
;        'gregorian'.
;
; USES:
;    correct_date.pro
;    str.pro
;    var_type.pro
;
; PROCEDURE:
;    This function first converts the inputed time values into an explicit date 
;    format ([yyyy,mm,dd,hh,mm,ss]), and then converts from that intermediate 
;    format to the requested output format.
;
; EXAMPLE:
;    ; Define the middle of all days between 1 January 2016 and 31 December 
;    ; 2018, using the "days since" format referenced to 1 November 2015
;    ; (using the Gregorian calendar):
;      time_in = 31 + 30 + 0.5 + findgen( 366 + 365 + 365 )
;      format_in = 'days since 2015-11-01T00:00:00'
;    ; Convert to "decimal year" format:
;      format_out = 'decimal year'
;      result = convert_time_format( time_in, format_in, format_out, calendar='gregorian' )
;
; MODIFICATION HISTORY:
;    Written by:  Daithi A. Stone (dstone@lbl.gov), 2017-04-28
;    Modified:  DAS, 2017-06-12 (completed code for when input units are 
;        months;  fixed bug in final adjustment of input hours)
;    Modified:  DAS, 2017-06-16 (fixed bug in incrementing months when there 
;        are surplus days with a Gregorian calendar and the "days since" input 
;        format;  added implementation of 'decimal' input type;  completed 
;        'yyyy...' output format implementation)
;    Modified:  DAS, 2017-10-10 (added 'yyyy...' to 'days since...' conversion)
;    Modified:  DAS, 2018-02-22 (relaxed requirements on format of 
;        'days since' date;  added capability for negative input time values 
;        under 'days since' format;  corrected bug in accounting for leap years 
;        when the first time value under the 'days since' format is multiple 
;        years after the time origin)
;    Modified:  DAS, 2018-03-01 (added some robustness to reading of time 
;        origin;  fixed bug when first time element is after a leap day and the 
;        origin is before the leap day;  added capability for negative output 
;        values)
;    Modified:  DAS, 2018-09-03 (Removed "-1" term when adding in hours and 
;        minutes for the 'days since' output format, I cannot remember why they 
;        were there but they were messing up the output)
;    Modified:  DAS, 2018-11-12 (Corrected bug in reading non-"-" spacers in 
;        the output origin specification)
;-

;***********************************************************************

FUNCTION convert_time_format, $
    TIME_IN, FORMAT_IN, FORMAT_OUT, $
    CALENDAR=calendar

;***********************************************************************
; Constants and check on input

; Confirm that required values are input
n_time = n_elements( time_in )
if n_time eq 0 then stop
if not( keyword_set( format_in ) ) then stop
if n_elements( format_in ) ne 1 then stop
if not( keyword_set( format_out ) ) then stop
if n_elements( format_out ) ne 1 then stop

; The default calendar
if not( keyword_set( calendar ) ) then calendar = 'gregorian'

; Confirm legible format for inputs
format_in = strlowcase( format_in )
format_out = strlowcase( format_out )
calendar = strlowcase( calendar )

; Determine the numbers of days in each calendar month (0=January)
temp = [ '365_day', 'noleap', 'gregorian', 'proleptic_gregorian' ]
if max( calendar eq temp ) eq 1 then begin
  dinm = [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]
endif else if calendar eq '360_day' then begin
  dinm = 30 + intarr( 12 )
endif else begin
  stop
endelse

;***********************************************************************
; Convert from input format to explicit format

; Initialise vector for explicit time format
; ([ 'yyyy', 'mm', 'dd', 'hh', 'mm', 'ss' ])
time_explicit = intarr( 6, n_time )
; Initialise check to confirm that requested input format was found below
check_in = 0

; Convert from "yyyy..." format
if strpos( format_in, 'yyyy' ) eq 0 then begin
  ; Ensure time_in is of type string, for easy reading
  if var_type( time_in ) ne 7 then time_in = strtrim( string( time_in ), 2 )
  ; Extract years
  time_explicit[0,*] = strmid( time_in, 0, 4 )
  ; Extract months and smaller units
  if strpos( format_in, 'yyyymm' ) eq 0 then begin
    time_explicit[1,*] = strmid( time_in, 4, 2 )
    ; Extract days and smaller units
    if strpos( format_in, 'yyyymmdd' ) eq 0 then begin
      time_explicit[2,*] = strmid( time_in, 6, 2 )
      ; Extract hours and smaller units
      if strpos( format_in, 'yyyymmddhh' ) eq 0 then begin
        time_explicit[3,*] = strmid( time_in, 8, 2 )
        ; Extract minutes and smaller units
        if strpos( format_in, 'yyyymmddhhmm' ) eq 0 then begin
          time_explicit[4,*] = strmid( time_in, 10, 2 )
          ; Extract seconds
          if strpos( format_in, 'yyyymmddhhmmss' ) eq 0 then begin
            time_explicit[5,*] = strmid( time_in, 12, 2 )
          endif
        endif
      endif
    endif
  endif
  ; Note that we have converted the values
  check_in = 1
endif

; Convert from "decimal" format
if format_in eq "decimal" then begin
  ; Iterate through time
  for i_time = 0, n_time - 1 do begin
    ; Intialise temporary date vector for current time value
    time_current = intarr( 6 )
    ; Copy current time value and adjust to assist numerical precision
    temp_time_in = time_in[i_time] + 1. / 365. / 24. / 60. / 60. / 4.
    ; Determine the year
    time_current[0] = floor( temp_time_in )
    ; Determine the month and days
    temp_dinm = dinm
    if strpos( calendar, 'gregorian' ) ge 0 then begin
      if ( ( time_current[0] mod 4 ) eq 0 ) $
           and ( ( ( time_current[0] mod 100 ) ne 0 ) $
           or ( ( time_current[0] mod 400 ) eq 0 ) ) then begin
        temp_dinm[1] = 29
      endif
    endif
    temp_resid = ( temp_time_in - time_current[0] ) * total( temp_dinm )
    temp_dinm_cumulative = total( [ 0, temp_dinm ], cumulative=1 )
    id = where( temp_resid gt temp_dinm_cumulative, n_id )
    time_current[1] = id[n_id-1] + 1
    temp_resid = temp_resid - temp_dinm_cumulative[id[n_id-1]]
    time_current[2] = floor( temp_resid ) + 1
    ; Determine hours, minutes, seconds
    temp_resid = ( temp_resid - time_current[2] + 1 ) * 24
    time_current[3] = floor( temp_resid )
    temp_resid = ( temp_resid - time_current[3] ) * 60
    time_current[4] = floor( temp_resid )
    temp_resid = ( temp_resid - time_current[4] ) * 60
    time_current[5] = floor( temp_resid )
    ; Record values
    time_explicit[*,i_time] = time_current
  endfor
  ; Note that we have converted the values
  check_in = 1
endif

; Convert from "<units> since <origin>" format
if strpos( format_in, ' since ' ) gt 0 then begin
  ; Determine the units and origin
  temp = strsplit( format_in, ' ', extract=1, count=n_temp )
  if ( n_temp lt 3 ) or ( n_temp gt 4 )  then stop
  units = temp[0]
  origin = temp[2]
  if n_temp eq 4 then origin = origin + ' ' + temp[3]
  ; Confirm that units are supported
  temp = [ 'days' ]
  if max( units eq temp ) ne 1 then stop
  ; Parse origin into components and confirm that it is legal
  origin_parsed = lonarr( 6 )
  temp = strsplit( origin, ' -:tz', extract=1, count=n_temp )
  if ( n_temp lt 3 ) or ( n_temp gt 6 ) then stop
  origin_parsed[0] = fix( temp[0] )
  origin_parsed[1] = fix( temp[1] )
  if ( origin_parsed[1] lt 1 ) or ( origin_parsed[1] gt 12 ) then stop
  if n_temp ge 3 then begin
    origin_parsed[2] = fix( temp[2] )
    if origin_parsed[2] lt 1 then stop
    if origin_parsed[2] gt dinm[origin_parsed[1]-1] then begin
      if strpos( calendar, 'gregorian' ) ge 0 then begin
        if origin_parsed[1] ne 2 then begin
          stop
        endif else if origin_parsed[2] ne 29 then begin
          stop
        endif else if origin_parsed[0] mod 4 ne 0 then begin
          stop
        endif else if ( origin_parsed[0] mod 100 eq 0 ) $
            and ( origin_parsed[0] mod 400 ne 0 ) then begin
          stop
        endif
      endif else begin
        stop
      endelse
    endif
  endif
  if n_temp ge 4 then begin
    origin_parsed[3] = fix( temp[3] )
    if ( origin_parsed[3] lt 0 ) or ( origin_parsed[3] gt 23 ) then stop
  endif
  if n_temp ge 5 then begin
    origin_parsed[4] = fix( temp[4] )
    if ( origin_parsed[4] lt 0 ) or ( origin_parsed[4] gt 59 ) then stop
  endif
  if n_temp ge 6 then begin
    origin_parsed[5] = fix( temp[5] )
    if ( origin_parsed[5] lt 0 ) or ( origin_parsed[5] gt 59 ) then stop
  endif
  ; Ensure that all values are positive (this code cannot handle negative time 
  ; values)
  if min( time_in ) lt 0 then begin
    ;; Parse output origin if it will be needed
    ;if strpos( format_out, 'days since ' ) ge 0 then begin
    ;  ; Determine the origin
    ;  temp = strsplit( format_out, ' ', extract=1, count=n_temp )
    ;  if n_temp ne 3 then stop
    ;  origin_out = temp[2]
    ;  ; Parse origin into components and confirm that it is legal
    ;  origin_parsed_out = intarr( 6 )
    ;  temp = strsplit( origin_out, '-:tz', extract=1, count=n_temp )
    ;  if ( n_temp lt 3 ) or ( n_temp gt 6 ) then stop
    ;  if strlen( temp[0] ) ne 4 then stop
    ;  origin_parsed_out[0] = fix( temp[0] )
    ;  if ( strlen( temp[1] ) lt 1 ) or ( strlen( temp[1] ) gt 2 ) then stop
    ;  origin_parsed_out[1] = fix( temp[1] )
    ;  if ( origin_parsed_out[1] lt 1 ) or ( origin_parsed_out[1] gt 12 ) then $
    ;      stop
    ;  if n_temp ge 3 then begin
    ;    if ( strlen( temp[2] ) lt 1 ) or ( strlen( temp[2] ) gt 2 ) then stop
    ;    origin_parsed_out[2] = fix( temp[2] )
    ;    if origin_parsed_out[2] lt 1 then stop
    ;    if origin_parsed_out[2] gt dinm[origin_parsed_out[1]-1] then begin
    ;      if strpos( calendar, 'gregorian' ) ge 0 then begin
    ;        if origin_parsed_out[1] ne 2 then begin
    ;          stop
    ;        endif else if ( origin_parsed_out[0] mod 4 ) ne 0 then begin
    ;          stop
    ;        endif else if ( origin_parsed_out[0] mod 100 eq 0 ) $
    ;            and ( origin_parsed_out[0] mod 400 ne 0 ) then begin
    ;          stop
    ;        endif
    ;      endif else begin
    ;        stop
    ;      endelse
    ;    endif
    ;  endif
    ;  if n_temp ge 4 then begin
    ;    if strlen( temp[3] ) ne 2 then stop
    ;    origin_parsed_out[3] = fix( temp[3] )
    ;    if ( origin_parsed_out[3] lt 0 ) or ( origin_parsed_out[3] gt 23 ) $
    ;        then stop
    ;  endif
    ;  if n_temp ge 5 then begin
    ;    if strlen( temp[4] ) ne 2 then stop
    ;    origin_parsed_out[4] = fix( temp[4] )
    ;    if ( origin_parsed_out[4] lt 0 ) or ( origin_parsed_out[4] gt 59 ) $
    ;        then stop
    ;  endif
    ;  if n_temp ge 6 then begin
    ;    if strlen( temp[5] ) ne 2 then stop
    ;    origin_parsed_out[5] = fix( temp[5] )
    ;    if ( origin_parsed_out[5] lt 0 ) or ( origin_parsed_out[5] gt 59 ) $
    ;        then stop
    ;  endif
    ;endif
    ; Shift to earlier origin to ensure positive values
    if max( units eq 'days' ) eq 1 then begin
      temp_years = ceil( abs( min( time_in ) ) / 366. )
      temp_days = 0
      ;if keyword_set( origin_parsed_out ) then temp_days_out = 0
      for i_years = 0, temp_years - 1 do begin
        temp_dinm = dinm
        ;if keyword_set( origin_parsed_out ) then temp_dinm_out = temp_dinm
        if strpos( calendar, 'gregorian' ) ge 0 then begin
          if origin_parsed[1] le 2 then begin
            temp_origin_year = origin_parsed[0] - i_years - 1
          endif else begin
            temp_origin_year = origin_parsed[0] - i_years
          endelse
          if ( ( temp_origin_year mod 4 ) eq 0 ) $
              and ( ( ( temp_origin_year mod 100 ) ne 0 ) $
              or ( ( temp_origin_year mod 400 ) eq 0 ) ) then begin
            temp_dinm[1] = 29
          endif
          ;if keyword_set( origin_parsed_out ) then begin
          ;  if origin_parsed_out[1] le 2 then begin
          ;    temp_origin_year_out = origin_parsed_out[0] - i_years - 1
          ;  endif else begin
          ;    temp_origin_year_out = origin_parsed_out[0] - i_years
          ;  endelse
          ;  if ( ( temp_origin_year_out mod 4 ) eq 0 ) $
          ;      and ( ( ( temp_origin_year_out mod 100 ) ne 0 ) $
          ;      or ( ( temp_origin_year_out mod 400 ) eq 0 ) ) then begin
          ;    temp_dinm_out[1] = 29
          ;  endif
          ;endif
        endif
        temp_days = temp_days + total( temp_dinm )
        ;if keyword_set( origin_parsed_out ) then begin
        ;  temp_days_out = temp_days_out + total( temp_dinm_out )
        ;endif
      endfor
      ;if keyword_set( origin_parsed_out ) then begin
      ;  if temp_days_out ne temp_days then stop
      ;  origin_parsed_out[0] = origin_parsed_out[0] - temp_years
      ;endif
      origin_parsed[0] = origin_parsed[0] - temp_years
      time_in = time_in + temp_days
    endif else begin
      stop
    endelse
  endif
  ; Iterate through values in the time vector
  for i_time = 0, n_time - 1 do begin
    ; Initialise the current explicit value with the previous value (automatic) 
    ; or the origin
    if i_time eq 0 then time_current = origin_parsed
    ; Determine the increment value
    if i_time eq 0 then begin
      delta_time = time_in[i_time]
    endif else begin
      delta_time = time_in[i_time] - time_in[i_time-1]
    endelse
    ; Increment amount specified by units
    temp_units = units
    if temp_units eq 'years' then begin
      time_current[0] = time_current[0] + floor( delta_time )
      delta_time_old = delta_time
      delta_time = ( delta_time - floor( delta_time ) ); * ???
      if delta_time_old eq delta_time - 1 then begin
        temp_factor = total( dinm[time_current[1]-1:11] )
        if time_current[1] eq 2 then begin
          if strpos( calendar, 'gregorian' ) ge 0 then begin
            if ( ( ( time_current[0] - 1 ) mod 4 ) eq 0 ) $
                and ( ( ( ( time_current[0] - 1 ) mod 100 ) ne 0 ) $
                or ( ( ( time_current[0] - 1 ) mod 400 ) eq 0 ) ) then begin
              temp_factor = temp_factor + 1
            endif
          endif
        endif
      endif else if delta_time_old eq delta_time then begin
        temp_factor = total( dinm[0:time_current[1]-1:11] )
        if time_current[1] eq 2 then begin
          if strpos( calendar, 'gregorian' ) ge 0 then begin
            if ( ( ( time_current[0] - 1 ) mod 4 ) eq 0 ) $
                and ( ( ( ( time_current[0] - 1 ) mod 100 ) ne 0 ) $
                or ( ( ( time_current[0] - 1 ) mod 400 ) eq 0 ) ) then begin
              temp_factor = temp_factor + 1
            endif
          endif
        endif
      endif else begin
        stop
      endelse
      if delta_time ne 0 then stop ; Need to determine conversion factor
      if delta_time ne 0 then temp_units = 'years'
    endif
    if temp_units eq 'months' then begin
      time_current[1] = time_current[1] + floor( delta_time )
      delta_time = ( delta_time - floor( delta_time ) )
      if delta_time ne 0 then begin
        dinm_use = dinm[time_current[1]]
        if ( time_current[1] eq 2 ) $
            and ( strpos( strlowcase( calendar ), 'gregorian' ) ge 0  ) $
            then begin
          if ( ( ( time_current[0] - 1 ) mod 4 ) eq 0 ) $
              and ( ( ( ( time_current[0] - 1 ) mod 100 ) ne 0 ) $
              or ( ( ( time_current[0] - 1 ) mod 400 ) eq 0 ) ) then begin
            dinm_use = 29
          endif
        endif
        delta_time = delta_time * dinm_use
        temp_units = 'days'
      endif
    endif
    if temp_units eq 'days' then begin
      time_current[2] = time_current[2] + floor( delta_time )
      delta_time = ( delta_time - floor( delta_time ) ) * 24
      if delta_time ne 0 then temp_units = 'hours'
    endif
    if temp_units eq 'hours' then begin
      time_current[3] = time_current[3] + floor( delta_time )
      delta_time = ( delta_time - floor( delta_time ) ) * 60
      if delta_time ne 0 then temp_units = 'minutes'
    endif
    if temp_units eq 'minutes' then begin
      time_current[4] = time_current[4] + floor( delta_time )
      delta_time = ( delta_time - floor( delta_time ) ) * 60
      if delta_time ne 0 then temp_units = 'seconds'
    endif
    if temp_units eq 'seconds' then begin
      time_current[5] = time_current[5] + round( delta_time )
    endif
    ; Ensure legal seconds entry
    while( time_current[5] lt 0 ) do begin
      time_current[4] = time_current[4] - 1
      time_current[5] = time_current[5] + 60
    endwhile
    while( time_current[5] gt 59 ) do begin
      time_current[4] = time_current[4] + 1
      time_current[5] = time_current[5] - 60
    endwhile
    ; Ensure legal minutes entry
    while( time_current[4] lt 0 ) do begin
      time_current[3] = time_current[3] - 1
      time_current[4] = time_current[4] + 60
    endwhile
    while( time_current[4] gt 59 ) do begin
      time_current[3] = time_current[3] + 1
      time_current[4] = time_current[4] - 60
    endwhile
    ; Ensure legal hours entry
    while( time_current[3] lt 0 ) do begin
      time_current[2] = time_current[2] - 1
      time_current[3] = time_current[3] + 24
    endwhile
    while( time_current[3] gt 23 ) do begin
      time_current[2] = time_current[2] + 1
      time_current[3] = time_current[3] - 24
    endwhile
    ; Ensure legal months entry
    while( time_current[1] lt 1 ) do begin
      time_current[0] = time_current[0] - 1
      time_current[1] = time_current[1] + 12
    endwhile
    while( time_current[1] gt 12 ) do begin
      time_current[0] = time_current[0] + 1
      time_current[1] = time_current[1] - 12
    endwhile
    ; Ensure legal days entry
    while( time_current[2] lt 1 ) do begin
      time_current[1] = time_current[1] - 1
      if time_current[1] ge 0 then begin
        temp_month = time_current[1] mod 12
      endif else begin
        temp_month = 12 - ( abs( time_current[1] ) mod 12 )
      endelse
      temp_dinm = dinm[temp_month-1]
      if ( strpos( calendar, 'gregorian' ) ge 0 ) $
          and ( temp_month eq 2 ) then begin
        if ( ( time_current[0] mod 4 ) eq 0 ) $
            and ( ( ( time_current[0] mod 100 ) ne 0 ) $
            or ( ( time_current[0] mod 400 ) eq 0 ) ) then begin
          temp_dinm = 29
        endif
      endif
      time_current[2] = time_current[2] + temp_dinm
    endwhile
    flag_done = 0
    while flag_done eq 0 do begin
      if time_current[1] gt 0 then begin
        temp_month = ( time_current[1] - 1 ) mod 12
      endif else begin
        temp_month = 11 - ( abs( time_current[1] ) mod 12 )
      endelse
      temp_dinm = dinm[temp_month]
      if ( strpos( calendar, 'gregorian' ) ge 0 ) $
          and ( temp_month eq 1 ) then begin
        if ( ( time_current[0] mod 4 ) eq 0 ) $
            and ( ( ( time_current[0] mod 100 ) ne 0 ) $
            or ( ( time_current[0] mod 400 ) eq 0 ) ) then begin
          temp_dinm = 29
        endif
      endif
      if time_current[2] gt temp_dinm then begin
        time_current[1] = time_current[1] + 1
        time_current[2] = time_current[2] - temp_dinm    
        ; Correct months entry if necessary
        if time_current[1] eq 13 then begin
          time_current[0] = time_current[0] + 1
          time_current[1] = time_current[1] - 12
        endif
      endif else begin
        flag_done = 1
      endelse
    endwhile
    ; Record the entries
    time_explicit[*,i_time] = time_current
  endfor
  ; Note that we have converted the values
  check_in = 1
endif

; Determine if we did indeed convert the values
if check_in eq 0 then stop

;***********************************************************************
; Convert from explicit format to requested format

; Initialise check to confirm that requested output format was found below
check_out = 0

; Convert to "yyyy..." format
if strpos( format_out, 'yyyy' ) eq 0 then begin
  ; Initialise output (string array)
  time_out = strarr( n_time )
  ; Add years
  time_out = time_out $
      + str( reform( time_explicit[0,*] ), 0, length=4, filler='0' )
  ; Add months
  if strpos( format_out, 'yyyymm' ) eq 0 then begin
  time_out = time_out $
      + str( reform( time_explicit[1,*] ), 0, length=2, filler='0' )
    ; Add days
    if strpos( format_out, 'yyyymmdd' ) eq 0 then begin
      time_out = time_out $
          + str( reform( time_explicit[2,*] ), 0, length=2, filler='0' )
      ; Add hours
      if strpos( format_out, 'yyyymmddhh' ) eq 0 then begin
        time_out = time_out $
            + str( reform( time_explicit[3,*] ), 0, length=2, filler='0' )
        ; Add minutes
        if strpos( format_out, 'yyyymmddhhmm' ) eq 0 then begin
          time_out = time_out $
              + str( reform( time_explicit[4,*] ), 0, length=2, filler='0' )
          ; Add seconds
          if strpos( format_out, 'yyyymmddhhmmss' ) eq 0 then begin
            time_out = time_out $
                + str( reform( time_explicit[5,*] ), 0, length=2, filler='0' )
          endif
        endif
      endif
    endif
  endif
  ; Note that we have converted the values
  check_out = 1
endif

; Convert to "days since <origin>" format
if strpos( format_out, 'days since ' ) ge 0 then begin
  ; Confirm that time values are in ascending order
  if n_time gt 1 then begin
    temp = double( time_explicit )
    temp = temp[0,*] + ( temp[1,*] + ( temp[2,*] + ( temp[3,*] + ( temp[4,*] $
        + temp[5,*] / 60.d ) / 60.d ) / 24.d ) / 31.d ) / 12.d
    temp = reform( temp )
    temp = temp[1:n_time-1] - temp[0:n_time-2]
    if min( temp ) le 0 then stop
    temp = 0
  endif
  ; Determine the origin
  temp = strsplit( format_out, ' ', extract=1, count=n_temp )
  if n_temp eq 4 then begin
    ; This bit ensures that a space between the yyyy-mm-dd and hh:mm:ss parts 
    ; of the origin will not confuse the code
    if ( n_elements( strsplit( temp[2], '-' ) ) eq 3 ) $
        and ( n_elements( strsplit( temp[3], ':' ) ) ge 2 ) then begin
      temp[2] = temp[2] + '-' + temp[3]
      n_temp = 3
      temp = temp[0:2]
    endif
  endif
  if n_temp ne 3 then stop
  origin = temp[2]
  ; Parse origin into components and confirm that it is legal.
  ;; This may have been done already in a time vector adjustment exercise
  ;if keyword_set( origin_parsed_out ) then begin
  ;  origin_parsed = origin_parsed_out
  ; Otherwise parse the origin
  ;endif else begin
    origin_parsed = intarr( 6 )
    temp = strsplit( origin, '-:tz', extract=1, count=n_temp )
    if ( n_temp lt 3 ) or ( n_temp gt 6 ) then stop
    if strlen( temp[0] ) ne 4 then stop
    origin_parsed[0] = fix( temp[0] )
    if ( strlen( temp[1] ) lt 1 ) or ( strlen( temp[1] ) gt 2 ) then stop
    origin_parsed[1] = fix( temp[1] )
    if ( origin_parsed[1] lt 1 ) or ( origin_parsed[1] gt 12 ) then stop
    if n_temp ge 3 then begin
      if ( strlen( temp[2] ) lt 1 ) or ( strlen( temp[2] ) gt 2 ) then stop
      origin_parsed[2] = fix( temp[2] )
      if origin_parsed[2] lt 1 then stop
      if origin_parsed[2] gt dinm[origin_parsed[1]-1] then begin
        if strpos( calendar, 'gregorian' ) ge 0 then begin
          if origin_parsed[1] ne 2 then begin
            stop
          endif else if ( origin_parsed[0] mod 4 ) ne 0 then begin
            stop
          endif else if ( origin_parsed[0] mod 100 eq 0 ) $
              and ( origin_parsed[0] mod 400 ne 0 ) then begin
            stop
          endif
        endif else begin
          stop
        endelse
      endif
    endif
    if n_temp ge 4 then begin
      if strlen( temp[3] ) ne 2 then stop
      origin_parsed[3] = fix( temp[3] )
      if ( origin_parsed[3] lt 0 ) or ( origin_parsed[3] gt 23 ) then stop
    endif
    if n_temp ge 5 then begin
      if strlen( temp[4] ) ne 2 then stop
      origin_parsed[4] = fix( temp[4] )
      if ( origin_parsed[4] lt 0 ) or ( origin_parsed[4] gt 59 ) then stop
    endif
    if n_temp ge 6 then begin
      if strlen( temp[5] ) ne 2 then stop
      origin_parsed[5] = fix( temp[5] )
      if ( origin_parsed[5] lt 0 ) or ( origin_parsed[5] gt 59 ) then stop
    endif
  ;endelse
  ;; Confirm that all time values are positive (this code cannot yet handle 
  ;; negative values)
  ;if min( time_explicit[0,*] ) lt origin_parsed[0] then stop
  ;id = where( time_explicit[0,*] eq origin_parsed[0], n_id )
  ;if n_id gt 0 then begin
  ;  temp_time = time_explicit[*,id]
  ;  if min( temp_time[1,*] ) lt origin_parsed[1] then stop
  ;  id = where( temp_time[1,*] eq origin_parsed[1], n_id )
  ;  if n_id gt 0 then begin
  ;    temp_time = temp_time[*,id]
  ;    if min( temp_time[2,*] ) lt origin_parsed[2] then stop
  ;    id = where( temp_time[2,*] eq origin_parsed[2], n_id )
  ;    if n_id gt 0 then begin
  ;      temp_time = temp_time[*,id]
  ;      if min( temp_time[3,*] ) lt origin_parsed[3] then stop
  ;      id = where( temp_time[3,*] eq origin_parsed[3], n_id )
  ;      if n_id gt 0 then begin
  ;        temp_time = temp_time[*,id]
  ;        if min( temp_time[4,*] ) lt origin_parsed[4] then stop
  ;        id = where( temp_time[4,*] eq origin_parsed[4], n_id )
  ;        if n_id gt 0 then begin
  ;          temp_time = temp_time[*,id]
  ;          if min( temp_time[5,*] ) lt origin_parsed[5] then stop
  ;        endif
  ;      endif
  ;    endif
  ;  endif
  ;endif
  ; If we are working with negative values
  temp_time_explicit = str( time_explicit[0,*], length=4, filler='0' ) $
      + str( time_explicit[1,*], length=2, filler='0' ) $
      + str( time_explicit[2,*], length=2, filler='0' ) $
      + str( time_explicit[3,*], length=2, filler='0' ) $
      + str( time_explicit[4,*], length=2, filler='0' ) $
      + str( time_explicit[5,*], length=2, filler='0' )
  temp_time_explicit = min( temp_time_explicit )
  temp_origin_parsed = str( origin_parsed[0,*], length=4, filler='0' ) $
      + str( origin_parsed[1,*], length=2, filler='0' ) $
      + str( origin_parsed[2,*], length=2, filler='0' ) $
      + str( origin_parsed[3,*], length=2, filler='0' ) $
      + str( origin_parsed[4,*], length=2, filler='0' ) $
      + str( origin_parsed[5,*], length=2, filler='0' )
  if temp_time_explicit lt temp_origin_parsed then begin
    ; Substitute the lowest value as the origin, recording a delta to add back 
    ; later
    origin_parsed = strmid( temp_time_explicit, 0, 4 ) + '-' $
        + strmid( temp_time_explicit, 4, 2 ) + '-' $
        + strmid( temp_time_explicit, 6, 2 ) + '-' $
        + strmid( temp_time_explicit, 8, 2 ) + ':' $
        + strmid( temp_time_explicit, 10, 2 ) + ':' $
        + strmid( temp_time_explicit, 12, 2 )
    time_out_shift = convert_time_format( temp_origin_parsed, $
        'yyyymmddhhmmss', 'days since ' + origin_parsed, calendar=calendar )
    time_out_shift = time_out_shift[0]
    origin_parsed = fix( strsplit( origin_parsed, '-:', extract=1 ) )
  endif
  ; Initialise output vector
  time_out = fltarr( n_time )
  ; Iterate through time vectors
  for i_time = 0, n_time - 1 do begin
    ; Determine number of units between previous value (or origin) and this 
    ; value
    if i_time eq 0 then begin
      time_previous = origin_parsed
    endif else begin
      time_previous = time_explicit[*,i_time-1]
    endelse
    time_current = time_explicit[*,i_time]
    ; If the two steps are separated by seconds within the same minute
    if ( time_current[0] eq time_previous[0] ) $
        and ( time_current[1] eq time_previous[1] ) $
        and ( time_current[2] eq time_previous[2] ) $
        and ( time_current[3] eq time_previous[3] ) $
        and ( time_current[4] eq time_previous[4] ) then begin
      ; Tally any seconds from the same minute
      time_out[i_time] = time_out[i_time] $
          + ( time_current[5] - time_previous[5] ) / 60. / 60. / 24.
   ; Otherwise
    endif else begin
      ; Tally seconds from the end of the previous minute and the start of the 
      ; current minute
      if time_previous[5] ne 0 then begin
        time_out[i_time] = time_out[i_time] $
            + ( 60. - time_previous[5] ) / 60. / 60. / 24.
        time_previous[4] = time_previous[4] + 1
        if time_previous[4] eq 60 then begin
          time_previous[4] = 0
          time_previous[3] = time_previous[3] + 1
          if time_previous[3] eq 24 then begin
            time_previous[3] = 0
            time_previous[2] = time_previous[2] + 1
            if time_previous[2] gt 28 then begin
              temp = correct_date( time_previous[0], time_previous[1], $
                  time_previous[2], calendar=calendar )
              time_previous[0] = strmid( strtrim( string( temp ), 2 ), 0, 4 )
              time_previous[1] = strmid( strtrim( string( temp ), 2 ), 4, 2 )
              time_previous[2] = strmid( strtrim( string( temp ), 2 ), 6, 2 )
            endif
          endif
        endif
      endif
      time_out[i_time] = time_out[i_time] + time_current[5] / 60. / 60. / 24.
      ; If the two steps are separated by minutes within the same hour
      if ( time_current[0] eq time_previous[0] ) $
          and ( time_current[1] eq time_previous[1] ) $
          and ( time_current[2] eq time_previous[2] ) $
          and ( time_current[3] eq time_previous[3] ) then begin
        ; Tally any minutes from the same hour
        ;time_out[i_time] = time_out[i_time] $
        ;    + ( time_current[4] - time_previous[4] - 1 ) / 60. / 24.
        time_out[i_time] = time_out[i_time] $
            + ( time_current[4] - time_previous[4] ) / 60. / 24.
      ; Otherwise
      endif else begin
        ; Tally minutes from the end of the previous hour and the start of the 
        ; current hour
        if time_previous[4] ne 0 then begin
          time_out[i_time] = time_out[i_time] $
            + ( 60. - time_previous[4] ) / 60. / 24.
          time_previous[3] = time_previous[3] + 1
          if time_previous[3] eq 24 then begin
            time_previous[3] = 0
            time_previous[2] = time_previous[2] + 1
            if time_previous[2] gt 28 then begin
              temp = correct_date( time_previous[0], time_previous[1], $
                  time_previous[2], calendar=calendar )
              time_previous[0] = strmid( strtrim( string( temp ), 2 ), 0, 4 )
              time_previous[1] = strmid( strtrim( string( temp ), 2 ), 4, 2 )
              time_previous[2] = strmid( strtrim( string( temp ), 2 ), 6, 2 )
            endif
          endif
        endif
        time_out[i_time] = time_out[i_time] + time_current[4] / 60. / 24.
        ; If the two steps are separated by hours within the same day
        if ( time_current[0] eq time_previous[0] ) $
            and ( time_current[1] eq time_previous[1] ) $
            and ( time_current[2] eq time_previous[2] ) then begin
          ; Tally any hours from the same day
          ;time_out[i_time] = time_out[i_time] $
          ;    + ( time_current[3] - time_previous[3] - 1 ) / 24.
          time_out[i_time] = time_out[i_time] $
              + ( time_current[3] - time_previous[3] ) / 24.
        ; Otherwise
        endif else begin
          ; Tally hours from the end of the previous day and the start of the 
          ; current day
          if time_previous[3] ne 0 then begin
            time_out[i_time] = time_out[i_time] $
                + ( 24. - time_previous[3] ) / 24.
            time_previous[2] = time_previous[2] + 1
            if time_previous[2] gt 28 then begin
              temp = correct_date( time_previous[0], time_previous[1], $
                  time_previous[2], calendar=calendar )
              time_previous[0] = strmid( strtrim( string( temp ), 2 ), 0, 4 )
              time_previous[1] = strmid( strtrim( string( temp ), 2 ), 4, 2 )
              time_previous[2] = strmid( strtrim( string( temp ), 2 ), 6, 2 )
            endif
          endif
          time_out[i_time] = time_out[i_time] + time_current[3] / 24.
          ; If the two steps are separated by days within the same month
          if ( time_current[0] eq time_previous[0] ) $
              and ( time_current[1] eq time_previous[1] ) then begin
            ; Tally any days from the same month
            time_out[i_time] = time_out[i_time] $
                + time_current[2] - time_previous[2]
          ; Otherwise
          endif else begin
            ; Tally days from the end of the previous month
            time_out[i_time] = time_out[i_time] + dinm[time_previous[1]-1] $
                - time_previous[2] + 1
            if time_previous[1] eq 2 then begin
              if ( strpos( calendar, 'gregorian' ) ge 0 ) $
                  and ( time_previous[2] ne 29 ) then begin
                if ( ( time_previous[0] mod 4 ) eq 0 ) $
                    and ( ( ( time_previous[0] mod 100 ) ne 0 ) $
                    or ( ( time_previous[0] mod 400 ) eq 0 ) ) then begin
                  time_out[i_time] = time_out[i_time] + 1
                endif
              endif
            endif
            ; Tally any days from the start of the current month
            time_out[i_time] = time_out[i_time] + time_current[2] - 1
            ; If the two steps are separated by months within the same year
            if ( time_current[0] eq time_previous[0] ) $
                and ( time_current[1] ge time_previous[1] + 2 ) then begin
              ; Tally any days from interim months of the same year
              time_out[i_time] = time_out[i_time] $
                  + total( dinm[time_previous[1]:time_current[1]-2] )
              if strpos( calendar, 'gregorian' ) ge 0 then begin
                if ( time_previous[1] lt 2 ) and ( time_current[1] gt 2 ) $
                    then begin
                  if ( ( time_previous[0] mod 4 ) eq 0 ) $
                      and ( ( ( time_previous[0] mod 100 ) ne 0 ) $
                      or ( ( time_previous[0] mod 400 ) eq 0 ) ) then begin
                    time_out[i_time] = time_out[i_time] + 1
                  endif
                endif
              endif
            ; Otherwise if the years are different
            endif else if time_current[0] gt time_previous[0] then begin
              ; Tally any months from the end of the previous year
              if time_previous[1] ne 12 then begin
                time_out[i_time] = time_out[i_time] $
                    + total( dinm[time_previous[1]:11] )
                if time_previous[1] eq 1 then begin
                  if strpos( calendar, 'gregorian' ) ge 0 then begin
                    if ( ( time_previous[0] mod 4 ) eq 0 ) $
                        and ( ( ( time_previous[0] mod 100 ) ne 0 ) $
                        or ( ( time_previous[0] mod 400 ) eq 0 ) ) then begin
                      time_out[i_time] = time_out[i_time] + 1
                    endif
                  endif
                endif
              endif
              ; Tally any months from the start of the current year
              if time_current[1] ne 1 then begin
                time_out[i_time] = time_out[i_time] $
                    + total( dinm[0:time_current[1]-2] )
                if strpos( calendar, 'gregorian' ) ge 0 then begin
                  if time_current[1] gt 2 then begin
                    if ( ( time_current[0] mod 4 ) eq 0 ) $
                        and ( ( ( time_current[0] mod 100 ) ne 0 ) $
                        or ( ( time_current[0] mod 400 ) eq 0 ) ) then begin
                      time_out[i_time] = time_out[i_time] + 1
                    endif
                  endif
                endif
              endif
              ; Tally any intervening years
              if time_current[0] ge time_previous[0] + 2 then begin
                time_out[i_time] = time_out[i_time] $
                    + total( dinm, integer=1 ) $
                    * ( time_current[0] - time_previous[0] - 1 )
                if strpos( calendar, 'gregorian' ) ge 0 then begin
                  temp = indgen( time_current[0] - time_previous[0] - 1 ) $
                      + time_previous[0] + 1
                  id = where( ( ( temp mod 4 ) eq 0 ) $
                      and ( ( ( temp mod 100 ) ne 0 ) $
                      or ( ( temp mod 400 ) eq 0 ) ), n_id )
                  time_out[i_time] = time_out[i_time] + n_id
                endif
              endif
            endif
          endelse
        endelse
      endelse
    endelse
    ; Add time-since-previous time step to previous time steps time-since-origin
    if i_time gt 0 then time_out[i_time] = time_out[i_time] + time_out[i_time-1]
  endfor
  ; Add any shift needed to account for using a placeholder origin
  if keyword_set( time_out_shift ) then time_out = time_out - time_out_shift
  ; Note that we have converted the values
  check_out = 1
endif

; Convert to decimal format
if format_out eq 'decimal year' then begin
  ; Initialise output vector
  time_out = dblarr( n_time )
  ; Iterate through time
  for i_time = 0, n_time - 1 do begin
    ; Copy the year
    time_out[i_time] = reform( time_explicit[0,i_time] )
    ; Determine if this is a leap year
    leap_value = 0
    if strpos( calendar, 'gregorian' ) ge 0 then begin
      temp_year = reform( time_explicit[0,i_time] )
      if ( ( temp_year mod 4 ) eq 0 ) $
          and ( ( ( temp_year mod 100 ) ne 0 ) $
          or ( ( temp_year mod 400 ) eq 0 ) ) then begin
        leap_value = 1
      endif
    endif
    ; Determine days in preceding months of this year
    temp_month = reform( time_explicit[1,i_time] )
    if temp_month gt 1 then begin
      temp_dinm_before = total( dinm[0:temp_month-2] )      
      if temp_month gt 2 then temp_dinm_before = temp_dinm_before + leap_value
    endif else begin
      temp_dinm_before = 0.
    endelse
    ; Determine days in this year
    temp_dina = total( dinm ) + leap_value
    ; Add days from preceding months of this year
    time_out[i_time] = time_out[i_time] + temp_dinm_before / temp_dina
    ; Add days, etc. from this month
    time_out[i_time] = time_out[i_time] $
        + ( time_explicit[2,i_time] - 1. $
        + ( time_explicit[3,i_time] $
        + ( time_explicit[4,i_time] + time_explicit[5,i_time] / 60. ) / 60. ) $
        / 24. ) / temp_dina
  endfor
  ; Note that we have converted the values
  check_out = 1
endif

; Determine if we did indeed convert the values
if check_out eq 0 then stop

;***********************************************************************
; The end

return, time_out
END
