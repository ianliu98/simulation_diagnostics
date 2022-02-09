;+
; NAME:
;    shengzwiers.pro
;
; PURPOSE:
;    This function applies the Sheng and Zwiers method for adjusting averaged  
;    data such that the original averaged values are re-obtained upon averaging 
;    of higher frequency data themselves calculated through linear 
;    interpolation from the adjusted data.  See Sheng and Zwiers (1998) and 
;    http://www-pcmdi.llnl.gov/projects/amip/AMIP2EXPDSN/BCS/amip2bcs.php for 
;    more details.
;
; CATEGORY:
;    Time Series Analysis
;
; CALLING SEQUENCE:
;    data_adj = shengzwiers( data_in, time )
;
; INPUTS:
;    DATA_IN:  A required floating point array of size N_SPACE*N_TIME*N_VAR 
;        containing the data to be adjusted.  The adjustmust is performed on 
;        the second (TIME) dimension, for each element in the first (SPACE) and 
;        third (VAR) dimensions.
;    TIME:  A required floating point vector of size N_TIME containing the 
;        times of the mid-points of the time steps, in units of time since a 
;        reference (e.g. "days since 2017-01-01T00:00:00".  Note that accurate 
;        calculations require that the time units are absolute (e.g. "days") 
;        and not variable (e.g. "months", which vary from 28 to 31 days in 
;        length).
;
; KEYWORD PARAMETERS:
;    DOUBLE:  If set then the inversion calculation is performed in double 
;        precision.  The default is the format of the data.
;    EXTEND_CYCLE:  If this is set to a positive integer then the time-ends of 
;        the data are padded with the average cycle of length EXTEND_CYCLE 
;        estimated along the time dimension in DATA.  For instance if 
;        EXTEND_CYCLE=12 and the data contain monthly averages, then the 
;        average annual cycle is calculated and used to pad the ends of the 
;        time series.  The default is 0, i.e. no padding.  Note that this needs 
;        to be 12 if padding monthly data, or 365 if padding 365-day-calendar 
;        daily data.
;
; OUTPUTS:
;    DATA_ADJ:  Returns a floating point array of identical size to DATA_IN 
;        containing the adjusted data.
;
; USES:
;    shengzwiers_weights.pro
;    dimension.pro
;
; PROCEDURE:
;    This function applies the method described in:
;      Sheng, J., and F. Zwiers.  1998.  An improved scheme for time-dependent 
;      boundary conditions in atmospheric general circulation models.  Climate 
;      Dynamics, 14, 609-613.
;
; EXAMPLE:
;    data_in = reform( [ 0., 1., 0., -1., 0. ], 1, 5, 1 )
;    data_out = shengzwiers( data_in, findgen( 5 ) )
;    ; data_out should be [ -0.195122, 1.36585, 0.00000, -1.36585, 0.195122 ].
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone (dastone@runbox.com), 2016-10-25, as 
;        sheng_zwiers.pro.
;    Modified:  DAS, 2017-10-10 (Branched to shengzwiers.pro;  standardised 
;        documentation;  removed requirement for monthly data)
;    Modified:  DAS, 2017-10-28 (Fixed bug in implementation of EXTEND_CYCLE)
;    Modified:  DAS, 2018-09-13 (Fixed additional bug in implementation of 
;        EXTEND_CYCLE and added more documentation about it)
;-

;***********************************************************************

FUNCTION SHENGZWIERS, $
    DATA_IN, TIME, $
    EXTEND_CYCLE=extend_cycle, $
    DOUBLE=double_opt

;***********************************************************************
; Constants

; The missing data flag
nan = !values.f_nan

; Confirm that we have a number of dimensions we can work with
temp = dimension( data_in )
if ( temp ne 2 ) and ( temp ne 3 ) then stop
; Determine the dimensions of the data set
n_space = n_elements( data_in[*,0,0] )
n_time = n_elements( data_in[0,*,0] )
n_var = n_elements( data_in[0,0,*] )

; Default no padding
if not( keyword_set( extend_cycle ) ) then extend_cycle = 0

; Option for double precision arithmetic
double_opt = keyword_set( double_opt )

;***********************************************************************
; Pad ends with average cycle

; If padding is requested
if extend_cycle gt 0 then begin
  ; Calculate average cycle values
  data_mean = fltarr( n_space, extend_cycle, n_var )
  for i_cycle = 0, extend_cycle - 1 do begin
    index = i_cycle + indgen( $
        + extend_cycle * ceil( ( 1. * n_time - i_cycle ) / extend_cycle ) )
    data_mean[*,i_cycle,*] = reform( $
        mean( data_in[*,index,*], dimension=2, nan=1 ), $
        n_space, 1, n_var )
  endfor
  ; Add this average cycle to data_in
  index = indgen( extend_cycle )
  temp = n_time mod extend_cycle
  if temp ne 0 then index = [ index[temp:extend_cycle-1], index[0:temp-1] ]
  data_use = [ [ data_mean ], [ data_in ], [ data_mean[*,index,*] ] ]
  ; (Note these times assume some regularity in the time stepping)
  time_use = [ time[0:extend_cycle-1] - ( time[extend_cycle] - time[0] ), $
      time, $
      time[n_time-extend_cycle:n_time-1] $
      + ( time[n_time-1] - time[n_time-1-extend_cycle] ) ]
  ; Clear memory
  data_mean = 0
; If padding is not requested
endif else begin
  ; Copy the input data
  data_use = data_in
  time_use = time
endelse
; The length of the (possibly revised) time vector
n_time_use = n_elements( time_use )

; Determine the lengths of each time step
d_time = time_use[1:n_time_use-1] - time_use[0:n_time_use-2]
; (Assume symmetry of time step lengths surrounding the first and last time 
; step)
time_len_use = [ d_time[0] * 2, $
    d_time[0:n_time_use-3] + d_time[1:n_time_use-2], d_time[n_time_use-2] * 2 ]
if double_opt eq 1 then time_len_use = double( time_len_use )
time_len_use = time_len_use / 2.

;***********************************************************************
; Do adjustments to produce mean-preserving interpolation

; Initialise output data array
if double_opt eq 1 then begin
  data_adj = double( data_use )
endif else begin
  data_adj = data_use
endelse

; Set up a linear system matrix for the full time series
if max( time_len_use, min=temp ) eq temp then begin
  linear_sys_full = shengzwiers_weights( n_month=n_time_use, double=double_opt )
endif else begin
  linear_sys_full = shengzwiers_weights( time_len_use, double=double_opt )
endelse
; Invert the matrix
linear_sys_full_inv = invert( linear_sys_full, double=double_opt )
linear_sys_full = 0

; Iterate through the var and space dimensions
for i_var = 0, n_var - 1 do begin
  for i_space = 0, n_space - 1 do begin
    ; Extract the time series
    temp_data = reform( data_use[i_space,*,i_var] )
    if double_opt eq 1 then temp_data = double( temp_data )
    ; Copy time series for determining segments
    temp_data_seg = temp_data
    ; Only take elements with valid data
    id = where( finite( temp_data_seg ) eq 1, n_id )
    temp_data_seg = 0
    ; If there is no data then mark as having no segments
    if n_id eq 0 then begin
      n_seg = 0
    ; If there is only one value then count it as a lone segment
    endif else if n_id eq 1 then begin
      n_seg = 1
      id_seg_start = id
      id_seg_end = id
    ; Otherwise identify the multiple segments
    endif else begin
      id_seg_start = id[0]
      id_seg_end = id[n_id-1]
      temp = id[1:n_id-1] - id[0:n_id-2]
      id_seg = where( temp gt 1, n_seg )
      if n_seg gt 0 then begin
        id_seg_start = [ id_seg_start, id[1+id_seg] ]
        id_seg_end = [ id[id_seg], id_seg_end ]
      endif
      n_seg = n_seg + 1
    endelse
    ; Iterate through the segments
    for i_seg = 0, n_seg - 1 do begin
      ; Check if these segments are identical to the previous round
      check = 0
      if n_elements( id_seg_start_old ) ne 0 then begin
        if ( id_seg_start[i_seg] eq id_seg_start_old ) $
            and ( id_seg_end[i_seg] eq id_seg_end_old ) then check = 1
      endif
      ; Otherwise check if we can adopt the prepared matrix for the full time
      ; series
      if check eq 0 then begin
        if ( id_seg_start[i_seg] eq 0 ) $
            and ( id_seg_end[i_seg] eq n_time_use - 1 ) then begin
          linear_sys_inv = linear_sys_full_inv
          check = 1
        endif
      endif
      ; If we need to set up a new system of equations
      if check eq 0 then begin
        ; Set up linear system matrix
        if id_seg_end[i_seg] eq id_seg_start[i_seg] then begin
          linear_sys = [ 1. ]
          if double_opt eq 1 then linear_sys = double( linear_sys )
        endif else if $
            max( time_len_use[id_seg_start[i_seg]:id_seg_end[i_seg]], $
            min=temp ) $
            eq temp then begin
          len_seg = id_seg_end[i_seg] - id_seg_start[i_seg] + 1
          linear_sys = shengzwiers_weights( n_month=len_seg, double=double_opt )
        endif else begin
          linear_sys = shengzwiers_weights( $
              time_len_use[id_seg_start[i_seg]:id_seg_end[i_seg]], $
              double=double_opt )
        endelse
        ; Invert the matrix
        linear_sys_inv = invert( linear_sys, double=double_opt )
      endif
      ; Calculate adjusted values
      temp_data[id_seg_start[i_seg]:id_seg_end[i_seg]] = linear_sys_inv $
          # temp_data[id_seg_start[i_seg]:id_seg_end[i_seg]]
      ; Record segment positions
      id_seg_start_old = id_seg_start[i_seg]
      id_seg_end_old = id_seg_end[i_seg]
    endfor
    ; Record adjusted time series
    data_adj[i_space,*,i_var] = temp_data
  endfor
endfor

;***********************************************************************
; Trim any padding

; If we did pad the ends
if keyword_set( extend_cycle ) then begin
  ; Remove the paddings
  data_adj = data_adj[*,extend_cycle:n_time_use-extend_cycle-1,*]
endif

;***********************************************************************
; The end

return, data_adj
END
