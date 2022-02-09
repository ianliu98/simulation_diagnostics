;+
; NAME:
;    c20c_dtos_v2_adjust_sic
;
; PURPOSE:
;    This procedure adjusts the sea ice coverage from an input data set by 
;    interpreting the effect of a delta-sea-surface-temperature field on input 
;    coverage and sea surface temperature.
;
; CATEGORY:
;    C20C dtos v2
;
; CALLING SEQUENCE:
;    c20c_dtos_v2_adjust_sic, fit_method=fit_method, in_tos_data=in_tos_data, in_sic_data=in_sic_data, out_sic_data=out_sic_data
;
; INPUTS:
;    FIT_METHOD, FIT_PERIOD, FIT_SIC_FILE, FIT_TOS_FILE, IN_TOS_DATA, 
;      IN_SIC_DATA, IN_LON, IN_LAT, IN_TIME, DELTA_TOS_DATA, OCEAN_MASK_DATA, 
;      TOS_FREEZE
;
; KEYWORD PARAMETERS:
;    DELTA_TOS_DATA:  An optional floating-point array of size 
;        N_IN_LON*N_IN_LAT*N_IN_TIME containing the delta field which should be 
;        added on to the IN_TOS_DATA input data.  This needs to be the same 
;        size as IN_SIC_DATA and IN_TOS_DATA.  If it is input, then code 
;        adjusts the IN_SIC_DATA sea ice concentration values in order to be 
;        consistent with the new sea surface temperatures.  If it is not input, 
;        then the code ignores IN_SIC_DATA and instead directly estimates sea 
;        ice concentration values that match the sea surface temperatures.
;    FIT_METHOD:  An required string vector providing details about the fitting 
;        method.  Possible values are:
;        * 'pall2007':  The linear fit used in Pall (2007) and Pall et alii 
;          (2011), calculated as a line through the freezing/full-ice point and 
;          the mean of the part-ice cells.  This is calculated separately for 
;          the Northern and Southern Hemispheres.  See 
;          c20c_dtos_v2_adjust_sic_pall2007.pro for more details.
;        * 'stonepall2018':  The linear fit used in Stone and Pall (2018), 
;          based on a linear fit to the freezing/full-ice point and the mean of
;          the median TOS values for all part-ice percentile bins.  Further 
;          options must be specified, as described in 
;          c20c_dtos_v2_adjust_sic_stonepall2018.pro.
;    FIT_PERIOD:  An optional two-element string vector defining the beginning 
;        and end dates, respectively, of the data to use for calculating the 
;        functional relation.  The format is ['yyyymmdd','yyyymmdd'].  The 
;        default is ['20010101','20101231'].  This only applies when then 
;        'stonepall2018' method is requested.
;    FIT_SIC_FILE:  In the 'stonepall2018' method is requested, then this is a 
;        required string scalar or array listing the NetCDF files from which to 
;        read the sea ice concentration data which will be used to calculate 
;        the functional relation with sea surface temperature.  This can 
;        include partial or full directory paths, and also can include 
;        wildcards.  The dimensions of the data must be the same as for the 
;        data specified in FIT_TOS_FILE.
;    FIT_TOS_FILE:  In the 'stonepall2018' method is requested, then this is a 
;        required string scalar or array listing the NetCDF files from which to 
;        read the sea surface temperature data which will be used to calculate 
;        the functional relation with sea ice concentration.  This can 
;        include partial or full directory paths, and also can include 
;        wildcards.  The dimensions of the data must be the same as for the 
;        data specified in FIT_SIC_FILE.
;    IN_LAT:  An optional floating-point array of length N_IN_LAT containing 
;        the latitude dimension for IN_TOS_DATA, IN_SIC_DATA, and 
;        DELTA_TOS_DATA.  This is required if the fitting method produces 
;        different fits as a function of latitude.
;    IN_LON:  An optional floating-point array of length N_IN_LON containing 
;        the longitude dimension for IN_TOS_DATA, IN_SIC_DATA, and 
;        DELTA_TOS_DATA.  This is required if the fitting method produces 
;        different fits as a function of longitude.
;    IN_SIC_DATA:  A required float-point array of size 
;        N_IN_LON*N_IN_LAT*N_IN_TIME containing the baseline sea ice 
;        concentration data to be adjusted.  This is ignored if DELTA_TOS_DATA 
;        is not provided.  In units of '%' (i.e. between 0 and 100).
;    IN_TIME:  An optional string array of length N_IN_TIME containing the time 
;        dimension for IN_TOS_DATA, IN_SIC_DATA, and DELTA_TOS_DATA.  The 
;        format should be 'yyyymmdd'.  This is required if the fitting method 
;        produces different fits as a function of time-of-year (season).
;    IN_TOS_DATA:  A required float-point array of size 
;        N_IN_LON*N_IN_LAT*N_IN_TIME containing the baseline sea surface 
;        temperature data.  The DELTA_TOS_DATA array will be added to this data.
;    OCEAN_MASK_DATA:  An optional floating-point array defining the fraction 
;        the area of each grid cell that is ocean.  Of size N_IN_LON*N_IN_LAT, 
;        where these are the same sizes of the spatial dimensions of 
;        IN_TOS_DATA and IN_SIC_DATA.  Values should be 1 over ocean, 0 over 
;        land, and between 0 and 1 for partially-oceanic cells.
;    OUT_SIC_DATA:  Returns a floating-point array containing the adjusted 
;        sea ice concentration data.  Of same size as IN_SIC_DATA.  In units of 
;        '%' (i.e. between 0 and 100).
;    TOS_FREEZE:  An optional floating-point scalar defining the freezing point 
;        in Kelvin.  The default is 271.35K.
;    V1: If set then the routine reproduces the output of the v1-0 code.  The 
;        difference is that the v1-0 set ice coverage to full as soon as the 
;        freezing point was passed, whereas the v2-0 code ignores the freezing 
;        point while shifting along the algorithm.  There is also a small 
;        difference in the selection of bins used to calculate the SST-SIC 
;        function, and v1-0 uses a time varying land-sea mask to for the 
;        calculation of the function whereas v2 uses a time-invariant mask.
;
; OUTPUTS:
;    OUT_SIC_DATA
;
; USES:
;    add_dim.pro
;    c20c_dtos_v2_adjust_sic_diag.pro
;    c20c_dtos_v2_adjust_sic_pall2007.pro
;    c20c_dtos_v2_adjust_sic_stonepall2018.pro
;    str.pro
;
; PROCEDURE:
;    If DELTA_TOS_DATA is input, then this procedure adjusts the sea ice 
;    concentration values in IN_SIC_DATA according to the change in 
;    temperatures dictated by adding DELTA_TOS_DATA onto IN_TOS_DATA.  If 
;    DELTA_TOS_DATA in not input, then instead this procedure directly 
;    estimates sea ice concentration values that match the IN_TOS_DATA sea 
;    surface temperature values.  Both options use a functional relationship 
;    between sea surface temperature and sea ice concentration.
;
; EXAMPLE:
;    See c20c_dtos_v2_make_tossic.pro.
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone (dastone@runbox.com), 2017-06-18, as 
;        c20c_adjust_sic_pall.pro.
;    Modified:  DAS, 2017-10-18 (Branched from c20c_adjust_sic.pro;  $
;        standardised documentation and code;  added to IDL routine library)
;    Modified:  DAS, 2018-05-08 (Added V1 keyword in call to 
;        c20c_dtos_v2_adjust_sic_stonepall2018.pro)
;    Modified:  DAS, 2018-11-26 (Added capability to deal with IN_TOS_DATA in 
;        units of degrees Celsius and IN_SIC_DATA in units of fraction)
;-

;***********************************************************************

PRO C20C_DTOS_V2_ADJUST_SIC, $
    FIT_SIC_FILE=fit_sic_file, FIT_TOS_FILE=fit_tos_file, $
    FIT_METHOD=fit_method, $
    FIT_PERIOD=fit_period, $
    IN_TOS_DATA=in_tos_data, IN_SIC_DATA=in_sic_data, IN_LON=in_lon_data, $
      IN_LAT=in_lat_data, IN_TIME=in_time_data, $
    DELTA_TOS_DATA=delta_tos_data, $
    OCEAN_MASK_DATA=ocean_mask_data, $
    TOS_FREEZE=tos_freeze_0, $
    OUT_SIC_DATA=out_sic_data, $
    V1=v1_opt

;***********************************************************************
; Constants

; Missing value
nan = !values.f_nan
; The freezing temperature
if n_elements( tos_freeze_0 ) eq 1 then begin
  tos_freeze = tos_freeze_0
  if tos_freeze lt 100 then tos_freeze = tos_freeze + ctok
endif else begin
  tos_freeze = 273.15 - 1.8
endelse
; No and full ice coverage
sic_none = 0.
sic_full = 100.

; Default instruction set
if not( keyword_set( instruction ) ) then instruction = ''

; Determine whether a delta field is being used
delta_opt = keyword_set( delta_tos_data )
if delta_opt eq 0 then begin
  print, 'NOTE:  No delta field being used (c20c_dtos_v2_adjust_sic.pro).'
endif

; Default periods
if not( keyword_set( fit_period ) ) then fit_period = [ '20010101', '20101231' ]

; The default fitting method
if not( keyword_set( fit_method ) ) then begin
  ; This setting was used for the C20C+ D&A project's Nat-Hist/CMIP5-est1 
  ; scenario estimate
  fit_method = [ 'stonepall2018', 'lon global', 'lat hemispheres', $
      'time annual', 'linear anchored' ]
endif

; The option to reproduce the v1 code output
v1_opt = keyword_set( v1_opt )

;***********************************************************************
; Check input data

; The size of the input data
n_in_lon = n_elements( in_lon_data )
n_in_lat = n_elements( in_lat_data )
n_in_time = n_elements( in_time_data )

; Check dimensions
if n_in_lon eq 0 then begin
  n_in_lon = n_elements( in_tos_data[*,0,0] )
  if n_in_lon gt 1 then stop
endif
if n_in_lat eq 0 then begin
  n_in_lat = n_elements( in_tos_data[0,*,0] )
  if n_in_lat gt 1 then stop
endif
if n_in_time eq 0 then begin
  if max( fit_method eq 'time seasonal' ) eq 1 then stop
  n_in_time = n_elements( in_tos_data[0,0,*] )
endif

; Check in_tos_data satisfies dimensions
if n_elements( in_tos_data[*,0,0] ) ne n_in_lon then stop
if n_elements( in_tos_data[0,*,0] ) ne n_in_lat then stop
if n_elements( in_tos_data[0,0,*] ) ne n_in_time then stop
; Check in_sic_data satisfies dimensions
if n_elements( in_sic_data[*,0,0] ) ne n_in_lon then stop
if n_elements( in_sic_data[0,*,0] ) ne n_in_lat then stop
if n_elements( in_sic_data[0,0,*] ) ne n_in_time then stop
; Check delta_tos_data satisfies dimensions
if delta_opt eq 1 then begin
  if n_elements( delta_tos_data[*,0,0] ) ne n_in_lon then stop
  if n_elements( delta_tos_data[0,*,0] ) ne n_in_lat then stop
  if n_elements( delta_tos_data[0,0,*] ) ne n_in_time then stop
endif
; Check if ocean mask satisfies dimensions
if keyword_set( ocean_mask_data ) then begin
  if n_elements( ocean_mask_data[*,0] ) ne n_in_lon then stop
  if n_elements( ocean_mask_data[0,*] ) ne n_in_lat then stop
endif

; Confirm legal sea surface temperature values (in K)
if max( finite( in_tos_data ) ) eq 0 then stop
if min( in_tos_data, nan=1, max=temp ) lt 271.15 then begin
  ; Convert from Celsius if appropriate
  temp = min( in_tos_data, nan=1 )
  if ( temp gt -1.81 ) and ( temp lt -1.79 ) then begin
    in_tos_data = in_tos_data + 273.15
    ; Flag conversion back to Celsius at end of procedure
    convert_to_celsius_opt = 1
  endif else begin
    stop
  endelse
endif
if temp gt 323.15 then stop
; Confirm legal sea ice concentration values (in %)
if max( finite( in_sic_data ) ) eq 0 then stop
if min( in_sic_data, nan=1, max=temp ) lt sic_none - 0.01 then stop
if ( temp lt sic_full * 0.999 ) or ( temp gt sic_full * 1.001) then begin
  ; Convert from fraction if appropriate
  if ( sic_full eq 100 ) and ( temp gt 0.999 ) and ( temp lt 1.001 ) then begin
    in_sic_data = in_sic_data * sic_full
    ; Flag conversion back to fraction at end of procedure
    convert_to_fraction_opt = 1
  endif else begin
    stop
  endelse
endif

; Confirm that data_in_tos and data_in_sic are consistent for missing values
id_tos = where( finite( in_tos_data ) eq 1, n_id_tos )
id_sic = where( finite( in_sic_data ) eq 1, n_id_sic )
if n_id_sic ne n_id_tos then stop
if max( abs( id_tos - id_sic ) ) gt 0 then stop
if delta_opt eq 1 then begin
  if min( finite( delta_tos_data[id_tos] ) ) eq 0 then stop
endif
id_tos = -1
id_sic = -1

;***********************************************************************
; Determine the relation of sea surface temperature and sea ice concentration

; If the Pall (2007) method is requested
if fit_method[0] eq 'pall2007' then begin
  c20c_dtos_v2_adjust_sic_pall2007, fit_tos_data=fit_tos_data, $
      fit_sic_data=fit_sic_data, fit_lon=fit_lon_data, fit_lat=fit_lat_data, $
      fit_time=fit_time_data, fit_type=fit_type
  ; Ensure that we have the freezing point defined
  if not( keyword_set( tos_freeze ) ) then tos_freeze = fit_tos_data[0,0,0,0]
; If the Stone and Pall (2018) method is requested
endif else if fit_method[0] eq 'stonepall2018' then begin
  n_fit_method = n_elements( fit_method )
  temp_fit_method = fit_method[1:n_fit_method-1]
  c20c_dtos_v2_adjust_sic_stonepall2018, sic_file=fit_sic_file, $
      tos_file=fit_tos_file, fit_method=temp_fit_method, $
      ocean_mask_data=ocean_mask_data, period=fit_period, $
      tos_freeze=tos_freeze_0, fit_tos_data=fit_tos_data, $
      fit_sic_data=fit_sic_data, fit_type=fit_type, fit_lon=fit_lon_data, $
      fit_lat=fit_lat_data, fit_time=fit_time_data, v1=v1_opt
; Otherwise the method has not yet been implemented
endif else begin
  stop
endelse

; Check that data was output
if not( keyword_set( fit_type ) ) then stop
if not( keyword_set( fit_sic_data ) ) then stop
if not( keyword_set( fit_tos_data ) ) then stop
;; Ensure sea ice concentration data is in %
;temp = max( fit_sic_data, nan=1 )
;if ( temp gt 0.99 ) and ( temp lt 1.01 ) then begin
;  fit_sic_data = fit_sic_data * 100.
;endif

; Determine the dimensions of the relationship function
n_fit_point = n_elements( fit_tos_data[*,0,0,0] )
n_fit_lon = n_elements( fit_lon_data )
if n_fit_lon eq 0 then begin
  n_fit_lon = n_elements( fit_tos_data[0,*,0,0] )
  if n_fit_lon ne 1 then stop
endif else begin
  if n_elements( fit_tos_data[0,*,0,0] ) ne n_fit_lon then stop
endelse
n_fit_lat = n_elements( fit_lat_data )
if n_fit_lat eq 0 then begin
  n_fit_lat = n_elements( fit_tos_data[0,0,*,0] )
  if n_fit_lat ne 1 then stop
endif else begin
  if n_elements( fit_tos_data[0,0,*,0] ) ne n_fit_lat then stop
endelse
n_fit_time = n_elements( fit_time_data )
if n_fit_time eq 0 then begin
  n_fit_time = n_elements( fit_tos_data[0,0,0,*] )
  if n_fit_time ne 1 then stop
endif else begin
  if n_elements( fit_tos_data[0,0,0,*] ) ne n_fit_time then stop
endelse

;***********************************************************************
; Apply SIC-TOS relationship to adjusted sea surface temperatures

; Copy input dimension vectors
out_lon_data = in_lon_data
out_lat_data = in_lat_data
out_time_data = in_time_data
; Create output SIC data array
out_sic_data = in_sic_data

; Iterate through time
for i_time = 0, n_in_time - 1 do begin
  ; Identify the relevant time-of-year in the relational function array
  if n_fit_time eq 1 then begin
    id_fit_time = 0
  endif else begin
    id_fit_time = where( strmid( in_time_data, 4, 4 ) eq fit_time_data, $
        n_id_fit_time )
    if n_id_fit_time ne 1 then stop
    id_fit_time = id_fit_time[0]
  endelse
  ; Iterate through latitude
  for i_lat = 0, n_in_lat - 1 do begin
    ; Identify the relevant latitude in the relational function array
    if n_fit_lat eq 1 then begin
      id_fit_lat = 0
    endif else begin
      temp = abs( fit_lat_data - in_lat_data[i_lat] )
      id_fit_lat = where( temp eq min( temp ) )
      id_fit_lat = id_fit_lat[0]
    endelse
    ; Iterate through longitude
    for i_lon = 0, n_in_lon - 1 do begin
      ; Extract input values
      temp_in_sic = in_sic_data[i_lon,i_lat,i_time]
      temp_in_tos = in_tos_data[i_lon,i_lat,i_time]
      ; Proceed if this is not land
      if ( finite( temp_in_sic ) + finite( temp_in_tos ) eq 2 ) then begin
        ; Identify the relevant longitude in the relational function array
        if n_fit_lon eq 1 then begin
          id_fit_lon = 0
        endif else begin
          temp = abs( fit_lon_data - in_lon_data[i_lon] )
          id_fit_lon = where( temp eq min( temp ) )
          id_fit_lon = id_fit_lon[0]
        endelse
        ; Extract the output, delta, and fit data
        if delta_opt eq 1 then begin
          temp_delta_tos = delta_tos_data[i_lon,i_lat,i_time]
        endif else begin
          temp_delta_tos = 0.
        endelse
        temp_fit_tos = fit_tos_data[*,id_fit_lon,id_fit_lat,id_fit_time]
        temp_fit_sic = fit_sic_data[*,id_fit_lon,id_fit_lat,id_fit_time]
        ; If we are working with a linear function within the partial-ice regime
        if fit_type eq 'straight line' then begin
          ; Set completion flag
          flag_done = 0
          ; If we are ending below the freezing point
          ; (only used when reproducing v1 of the code, as in 
          ; Nat-Hist/CMIP5-est1/v1-0)
          if ( v1_opt eq 1 ) $
              and ( temp_in_tos + temp_delta_tos le tos_freeze + 0.001 ) $
               then begin
            ; Cover with full ice
            temp_out_sic = sic_full
            ; Flag completion
            flag_done = 1
          endif
          ; If we are starting below the no-ice point
          if ( flag_done eq 0 ) and ( temp_in_tos le max( temp_fit_tos ) ) $
              then begin
            ; Shift the concentration according to the function
            if delta_opt eq 1 then begin
              temp_out_sic = ( temp_fit_sic[1] - temp_fit_sic[0] ) $
                  / ( temp_fit_tos[1] - temp_fit_tos[0] ) * temp_delta_tos $
                  + temp_in_sic
            ; Or define the concentration according to the function
            endif else begin
              temp_out_sic = ( temp_fit_sic[1] - temp_fit_sic[0] ) $
                  / ( temp_fit_tos[1] - temp_fit_tos[0] ) $
                  * ( temp_in_tos - temp_fit_tos[1] ) + temp_fit_sic[1]
            endelse
            ; Flag completion
            flag_done = 1
          endif
          ; If we are above the no-ice point
          if ( flag_done eq 0 ) and ( temp_in_tos gt max( temp_fit_tos ) ) $
              then begin
            ; If we are not shifting
            if delta_opt eq 0 then begin
              temp_out_sic = sic_none
            ; Shift by the portion of the cooling in the partial-ice range of 
            ; the function or by the full warming
            endif else begin
              if temp_delta_tos lt 0 then begin
                temp_delta_tos = temp_delta_tos $
                    + ( temp_in_tos - temp_fit_tos[1] )
                if temp_delta_tos lt 0 then begin
                  temp_out_sic = ( temp_fit_sic[1] - temp_fit_sic[0] ) $
                      / ( temp_fit_tos[1] - temp_fit_tos[0] ) * temp_delta_tos $
                      + temp_in_sic
                endif else begin
                  temp_out_sic = temp_in_sic
                endelse
              endif else begin
                temp_out_sic = ( temp_fit_sic[1] - temp_fit_sic[0] ) $
                    / ( temp_fit_tos[1] - temp_fit_tos[0] ) * temp_delta_tos $
                    + temp_in_sic
              endelse
            endelse
            ; Flag completion
            flag_done = 1
          endif
          ; Keep sea ice concentration in the legal range
          if temp_out_sic gt sic_full then begin
            temp_out_sic = sic_full
          endif else if temp_out_sic lt sic_none then begin
            temp_out_sic = sic_none
          endif
          ; Record new value
          out_sic_data[i_lon,i_lat,i_time] = temp_out_sic
        ; If the function is a nonlinear one following the bin values
        endif else if fit_type eq 'follow bins' then begin
          ; This is not yet implemented
          stop
        endif else begin
          stop
        endelse
      endif
    endfor
  endfor
endfor

;***********************************************************************
; Plot diagnostics if requested

; If this is requested
diag_opt = 1
if keyword_set( diag_opt ) then begin
  ; Iterate through time fits
  for i_fit_time = 0, n_fit_time - 1 do begin
    if i_fit_time gt 0 then stop ; Not yet implemented
    n_id_time = n_in_time
    id_time = indgen( n_id_time )
    ; Iterate through latitude fits
    for i_fit_lat = 0, n_fit_lat - 1 do begin
      ; Identify latitudes using this fit
      id_lat = -1
      for i_lat = 0, n_in_lat - 1 do begin
        temp = abs( in_lat_data[i_lat] - fit_lat_data )
        if temp[i_fit_lat] eq min( temp ) then id_lat = [ id_lat, i_lat ]
      endfor
      n_id_lat = n_elements( id_lat ) - 1
      if n_id_lat eq 0 then stop
      id_lat = id_lat[1:n_id_lat]
      ; Iterate through longitude fits
      for i_fit_lon = 0, n_fit_lon - 1 do begin
        ; Identify longitudes using this fit
        id_lon = -1
        for i_lon = 0, n_in_lon - 1 do begin
          temp = abs( in_lon_data[i_lon] - fit_lon_data )
          if temp[i_fit_lon] eq min( temp ) then id_lon = [ id_lon, i_lon ]
        endfor
        n_id_lon = n_elements( id_lon ) - 1
        if n_id_lon eq 0 then stop
        id_lon = id_lon[1:n_id_lon]
        ; Extract input and output data
        in_sic_data_plot = in_sic_data[*,*,id_time]
        in_sic_data_plot = in_sic_data_plot[*,id_lat,*]
        in_sic_data_plot = in_sic_data_plot[id_lon,*,*]
        in_tos_data_plot = in_tos_data[*,*,id_time]
        in_tos_data_plot = in_tos_data_plot[*,id_lat,*]
        in_tos_data_plot = in_tos_data_plot[id_lon,*,*]
        out_sic_data_plot = out_sic_data[*,*,id_time]
        out_sic_data_plot = out_sic_data_plot[*,id_lat,*]
        out_sic_data_plot = out_sic_data_plot[id_lon,*,*]
        out_tos_data_plot = in_tos_data[*,*,id_time] $
            + delta_tos_data[*,*,id_time]
        out_tos_data_plot = out_tos_data_plot[*,id_lat,*]
        out_tos_data_plot = out_tos_data_plot[id_lon,*,*]
        id = where( out_tos_data_plot lt tos_freeze, n_id )
        if n_id gt 0 then out_tos_data_plot[id] = tos_freeze
        ocean_mask_data_plot = ocean_mask_data[*,id_lat]
        ocean_mask_data_plot = ocean_mask_data_plot[id_lon,*]
        id = where( ocean_mask_data_plot eq 0, n_id )
        if n_id gt 0 then ocean_mask_data_plot[id] = !values.f_nan
        id = where( finite( ocean_mask_data_plot ) eq 1, n_id )
        ocean_mask_data_plot[id] = 1.
        ocean_mask_data_plot = add_dim( ocean_mask_data_plot, 2, n_id_time )
        in_sic_data_plot = in_sic_data_plot * ocean_mask_data_plot
        in_tos_data_plot = in_tos_data_plot * ocean_mask_data_plot
        out_sic_data_plot = out_sic_data_plot * ocean_mask_data_plot
        out_tos_data_plot = out_tos_data_plot * ocean_mask_data_plot
        ; Extract input and output dimension
        out_lon_data_plot = out_lon_data[id_lon]
        out_lat_data_plot = out_lat_data[id_lat]
        out_time_data_plot = out_time_data[id_time]
        ; Extract fit data
        fit_sic_data_plot = fit_sic_data[*,i_fit_lon,i_fit_lat,i_fit_time]
        fit_tos_data_plot = fit_tos_data[*,i_fit_lon,i_fit_lat,i_fit_time]
        ; Define the file name
        if fit_time_data[i_fit_time] eq '' then begin
          temp_fit_time_data = '00000000'
        endif else begin
          temp_fit_time_data = fit_time_data[i_fit_time]
        endelse
        file_diag = 'c20c_dtos_v2_adjust_sic_diag_' $
            + str( fit_lon_data[i_fit_lon], 1, length=3, filler='0' ) + '_' $
            + str( fit_lat_data[i_fit_lat], 1, length=2, filler='0' ) + '_' $
            + temp_fit_time_data + '.ps'
        ; Generate plots
        c20c_dtos_v2_adjust_sic_diag, fit_sic_data=fit_sic_data_plot, $
            fit_tos_data=fit_tos_data_plot, in_sic_data=in_sic_data_plot, $
            in_tos_data=in_tos_data_plot, adj_sic_data=out_sic_data_plot, $
            adj_tos_data=out_tos_data_plot, lon_data=out_lon_data_plot, $
            lat_data=out_lat_data_plot, time_data=out_time_data_plot, $
            file_name=file_diag
      endfor
    endfor
  endfor
endif

; Convert back to Celcius if required
if keyword_set( convert_to_celsius_opt ) then in_tos_data = in_tos_data - 273.15
; Convert back to fraction if required
if keyword_set( convert_to_fraction_opt ) then begin
  in_sic_data = in_sic_data / sic_full
  out_sic_data = out_sic_data / sic_full
endif

;***********************************************************************
; The end

;stop
return
END
