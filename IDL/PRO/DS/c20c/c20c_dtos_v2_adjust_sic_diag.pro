;+
; NAME:
;    c20c_dtos_v2_adjust_sic_diag
;
; PURPOSE:
;    This procedure produces a postscript set of plots which report on the 
;    sea ice adjustments performed in c20c_dtos_adjust_sic.pro.
;
; CATEGORY:
;    C20C dtos v2
;
; CALLING SEQUENCE:
;    c20c_dtos_v2_adjust_sic_diag
;
; INPUTS:
;    ADJ_SIC_DATA, ADJ_TOS_DATA, FILE_NAME, FIT_SIC_DATA, FIT_TOS_DATA, 
;      IN_SIC_DATA, IN_TOS DATA, LON_DATA, LAT_DATA, TIME_DATA
;
; KEYWORD PARAMETERS:
;    ADJ_SIC_DATA:  An optional floating-point array containing the adjusted  
;        sea ice concentration data, i.e. after adjustment.  Of size 
;        N_LON*N_LAT*N_TIME.  Required with ADJ_TOS_DATA for some plots, and 
;        also with IN_SIC_DATA and IN_TOS_DATA for some of those plots.  In 
;        units of %.
;    ADJ_TOS_DATA:  An optional floating-point array containing the adjusted 
;        sea surface temperature data, i.e. after adjustment.  Of size 
;        N_LON*N_LAT*N_TIME.  Required with ADJ_SIC_DATA for some plots, and 
;        also with IN_SIC_DATA and IN_TOS_DATA for some of those plots.  In 
;        units of K.
;    FILE_NAME:  An optional string containing the file name for the output 
;        postscript file.  The default is defined in ps_open.pro.
;    FIT_SIC_DATA:  An optional floating-point vector containing the SIC data 
;        in the SST-SIC function used to calculate the sea ice adjustment.  The 
;        array is of size N_FIT_POINT, the number of sea surface 
;        temperature/sea ice concentration points described by this array and 
;        FIT_TOS_DATA.  In units of Kelvin.  This vector must be the same size 
;        as FIT_TOS_DATA.  It will be ignored if FIT_TOS_DATA is not input, and 
;        the relevant plots will be omitted.
;    FIT_TOS_DATA:  An optional floating-point vector containing the SST data 
;        in the SST-SIC function used to calculate the sea ice adjustment.  The 
;        array is of size N_FIT_POINT, the number of sea surface 
;        temperature/sea ice concentration points described by this array and 
;        FIT_SIC_DATA.  In units of % (0 to 100)..  This vector must be the 
;        same size as FIT_SIC_DATA.  It will be ignored if FIT_SIC_DATA is not 
;        input, and the relevant plots will be omitted.
;    IN_SIC_DATA:  An optional floating-point array containing the original 
;        sea ice concentration data, before adjustment.  Of size 
;        N_LON*N_LAT*N_TIME.  Required with IN_TOS_DATA for some plots, and 
;        also with ADJ_SIC_DATA and ADJ_TOS_DATA for some of those plots.  In 
;        units of %.
;    IN_TOS_DATA:  An optional floating-point array containing the original 
;        sea surface temperature data, before adjustment.  Of size 
;        N_LON*N_LAT*N_TIME.  Required with IN_SIC_DATA for some plots, and 
;        also with ADJ_SIC_DATA and ADJ_TOS_DATA for some of those plots.  In 
;        units of K.
;    LAT_DATA:  An optional floating-point vector containing the latitude 
;        dimension in ADJ_SIC_DATA, ADJ_TOS_DAT, IN_SIC_DATA, and IN_TOS_DATA.  
;        LON_DATA, LAT_DATA, and TIME_DATA are required for plotting maps.  Of 
;        length N_LAT.
;    LON_DATA:  An optional floating-point vector containing the longitude 
;        dimension in ADJ_SIC_DATA, ADJ_TOS_DAT, IN_SIC_DATA, and IN_TOS_DATA.  
;        LON_DATA, LAT_DATA, and TIME_DATA are required for plotting maps.  Of 
;        length N_LON.
;    TIME_DATA:  An optional floating-point vector containing the time 
;        dimension in ADJ_SIC_DATA, ADJ_TOS_DAT, IN_SIC_DATA, and IN_TOS_DATA.  
;        LON_DATA, LAT_DATA, and TIME_DATA are required for plotting maps.  Of 
;        length N_TIME.
;
; OUTPUTS:
;    <FILE_NAME>
;
; USES:
;    arrows.pro
;    colorbrewer.pro
;    convert_time_format.pro
;    month_name.pro
;    ps_close.pro
;    ps_open.pro
;
; PROCEDURE:
;    This procedure generates various plots depending on provided input data.
;
; EXAMPLES:
;    See c20c_dtos_v2_adjust_sic.pro.
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone (dastone@runbox.com), 2017-06-18, as 
;        part of c20c_adjust_sic_pall.pro.
;    Modified:  DAS, 2017-11-22 (Extracted from c20c_adjust_sic_pall.pro into 
;        c20c_dtos_v2_adjust_sic_diag.pro)
;    Modified:  DAS, 201808-15 (Added to IDL routine library)
;-

;***********************************************************************

PRO C20C_DTOS_V2_ADJUST_SIC_DIAG, $
    FIT_SIC_DATA=fit_sic_data, FIT_TOS_DATA=fit_tos_data, $
    IN_SIC_DATA=in_sic_data, IN_TOS_DATA=in_tos_data, $
    ADJ_SIC_DATA=adj_sic_data, ADJ_TOS_DATA=adj_tos_data, $
    LON_DATA=lon_data, LAT_DATA=lat_data, TIME_DATA=time_data, $
    FILE_NAME=file_name

;***********************************************************************
; Constants

; The random number seed
seed = 1

; The size of the subset of original data to plot
n_subset = 3000l

; The sic plotting dimension
sic_none = 0.
sic_full = 100.
range_sic = [ sic_none, sic_full ]
style_sic = 2
title_sic = 'Sea ice coverage (%)'
; The tos plotting dimension
range_tos = [ 270., 280. ]
style_tos = 1
title_tos = 'Sea surface temperature (K)'
; Other plotting properties
charsize = 1.
font = 0
thick = 3
color_fit = 2
psym = 3

; Earth's radius (Mm) and surface area (Mm2)
r_earth = 6.37095
a_earth = 4. * !pi * r_earth ^ 2.

; Determine the dimensions of the relationship function
if keyword_set( fit_sic_data ) and keyword_set( fit_tos_data ) then begin
  n_fit_point = n_elements( fit_sic_data[*,0,0,0] )
  if n_elements( fit_sic_data ) ne n_fit_point then stop
  if n_elements( fit_tos_data ) ne n_fit_point then stop
endif else begin
  n_fit_point = 0
endelse

; Determine the existence of and dimensions of the input and data
if keyword_set( in_sic_data ) or keyword_set( adj_sic_data ) then begin
  n_lon = n_elements( lon_data )
  if n_lon eq 0 then begin
    if keyword_set( in_sic_data ) then begin
      n_lon = n_elements( in_sic_data[*,0,0] )
    endif else if keyword_set( adj_sic_data ) then begin
      n_lon = n_elements( adj_sic_data[*,0,0] )
    endif
  endif
  if n_lon gt 0 then begin
    if keyword_set( in_sic_data ) then begin
      if n_elements( in_sic_data[*,0,0] ) ne n_lon then stop
      if n_elements( in_tos_data[*,0,0] ) ne n_lon then stop
    endif
    if keyword_set( adj_sic_data ) then begin
      if n_elements( adj_sic_data[*,0,0] ) ne n_lon then stop
      if n_elements( adj_tos_data[*,0,0] ) ne n_lon then stop
    endif
  endif
  n_lat = n_elements( lat_data )
  if n_lat eq 0 then begin
    if keyword_set( in_sic_data ) then begin
      n_lat = n_elements( in_sic_data[0,*,0] )
    endif else if keyword_set( adj_sic_data ) then begin
      n_lat = n_elements( adj_sic_data[0,*,0] )
    endif
  endif
  if n_lat gt 0 then begin
    if keyword_set( in_sic_data ) then begin
      if n_elements( in_sic_data[0,*,0] ) ne n_lat then stop
      if n_elements( in_tos_data[0,*,0] ) ne n_lat then stop
    endif
    if keyword_set( adj_sic_data ) then begin
      if n_elements( adj_sic_data[0,*,0] ) ne n_lat then stop
      if n_elements( adj_tos_data[0,*,0] ) ne n_lat then stop
    endif
  endif
  n_time = n_elements( time_data )
  if n_time eq 0 then begin
    if keyword_set( in_sic_data ) then begin
      n_time = n_elements( in_sic_data[0,0,*] )
    endif else if keyword_set( adj_sic_data ) then begin
      n_time = n_elements( adj_sic_data[0,0,*] )
    endif
  endif
  if n_time gt 0 then begin
    if keyword_set( in_sic_data ) then begin
      if n_elements( in_sic_data[0,0,*] ) ne n_time then stop
      if n_elements( in_tos_data[0,0,*] ) ne n_time then stop
    endif
    if keyword_set( adj_sic_data ) then begin
      if n_elements( adj_sic_data[0,0,*] ) ne n_time then stop
      if n_elements( adj_tos_data[0,0,*] ) ne n_time then stop
    endif
  endif
endif else begin
  n_lon = 0
  n_lat = 0
  n_time = 0
endelse

;***********************************************************************
; Plot fitting process

; Open postscript output
ps_open, color=1, filename=file_name
device, helvetica=1, bold=1
tek_color

; Select a random subset of points/months for the plotting in the first three 
; panels
if ( keyword_set( in_sic_data ) and keyword_set( in_tos_data ) ) $
    or ( keyword_set( adj_sic_data ) and keyword_set( adj_tos_data ) ) $
    then begin
  ; Select the random subset
  if keyword_set( in_sic_data ) and keyword_set( adj_sic_data ) then begin
    id_plot = where( finite( in_sic_data ) + finite( in_tos_data ) $
        + finite( adj_sic_data ) + finite( adj_tos_data ) eq 4, n_id_plot )
  endif else if keyword_set( in_sic_data ) then begin
    id_plot = where( finite( in_sic_data ) + finite( in_tos_data ) eq 2, $
        n_id_plot )
  endif else if keyword_set( adj_sic_data ) then begin
    id_plot = where( finite( adj_sic_data ) + finite( adj_tos_data ) eq 2, $
        n_id_plot )
  endif else begin
    stop
  endelse
  if n_id_plot eq 0 then stop
  if keyword_set( in_sic_data ) and keyword_set( adj_sic_data ) then begin
    id = where( ( ( in_tos_data[id_plot] ge range_tos[0] ) $
        and ( in_tos_data[id_plot] le range_tos[1] ) ) $
        or ( ( adj_tos_data[id_plot] ge range_tos[0] ) $
        and ( adj_tos_data[id_plot] le range_tos[1] ) ), n_id_plot )
  endif else if keyword_set( in_sic_data ) then begin
    id = where( ( in_tos_data[id_plot] ge range_tos[0] ) $
        and ( in_tos_data[id_plot] le range_tos[1] ), n_id_plot )
  endif else if keyword_set( adj_sic_data ) then begin
    id = where( ( adj_tos_data[id_plot] ge range_tos[0] ) $
        and ( adj_tos_data[id_plot] le range_tos[1] ), n_id_plot )
  endif else begin
    stop
  endelse
  if n_id_plot eq 0 then stop
  id_plot = id_plot[id]
  if n_subset lt n_id_plot then begin
    id = floor( randomu( seed, n_subset ) * n_id_plot )
    id_plot = id_plot[id]
    n_id_plot = n_subset
  endif
endif

; Set the window for the first plot (original data)
!p.multi = [ 0, 3, 3 ]
; Proceed if there is data for this plot
if keyword_set( in_sic_data ) and keyword_set( in_tos_data ) then begin
  ; Prepare to plot original data for this domain
  plot, range_tos, range_sic, xthick=thick, ythick=thick, yrange=range_sic, $
      ystyle=style_sic, xrange=range_tos, xstyle=style_tos, charsize=charsize, $
      font=font, xtitle=title_tos, ytitle=title_sic, xticklen=0.04, $
      title='Before adjustment', nodata=1
  ; Plot a random subset of original data
  oplot, in_tos_data[id_plot], in_sic_data[id_plot], psym=3
  ; Plot fit line
  if n_fit_point gt 0 then begin
    oplot, fit_tos_data, fit_sic_data, color=color_fit, thick=thick
  endif
endif

; Set the window for the second plot (adjustment)
!p.multi[0] = 3 * 3 - 1
; Proceed if there is data for this plot
if ( n_fit_point gt 0 ) $
    or ( keyword_set( in_sic_data ) and keyword_set( adj_sic_data ) $
    and keyword_set( in_tos_data ) and keyword_set( adj_tos_data ) ) then begin
  ; Prepare to plot adjustment data for this domain
  plot, range_tos, range_sic, xthick=thick, ythick=thick, yrange=range_sic, $
      ystyle=style_sic, xrange=range_tos, xstyle=style_tos, charsize=charsize, $
      font=font, xtitle=title_tos, ytitle=title_sic, xticklen=0.04, $
      title='Adjustment', nodata=1
  ; Select a subset of the subset
  id = id_plot[0:n_subset/10-1]
  ; Plot a random subset of original data
  if keyword_set( in_sic_data ) and keyword_set( in_tos_data ) $
      and keyword_set( adj_sic_data ) and keyword_set( adj_tos_data ) then begin
    arrows, in_tos_data[id], in_sic_data[id], headx=adj_tos_data[id],$
        heady=adj_sic_data[id], absheadlength=0.25, solid=0, noclip=0
  endif
  ; Plot fit line
  if n_fit_point gt 0 then begin
    oplot, fit_tos_data, fit_sic_data, color=color_fit, thick=thick
  endif
endif

; Set the window for the third plot (adjusted)
!p.multi[0] = 3 * 3 - 2
; Proceed if there is data for this plot
if keyword_set( adj_sic_data ) and keyword_set( adj_tos_data ) then begin
  ; Prepare to plot adjusted data for this domain
  plot, range_tos, range_sic, xthick=thick, ythick=thick, yrange=range_sic, $
      ystyle=style_sic, xrange=range_tos, xstyle=style_tos, charsize=charsize, $
      font=font, xtitle=title_tos, ytitle=title_sic, xticklen=0.04, $
      title='After adjustment', nodata=1
  ; Plot a random subset of adjusted data
  oplot, adj_tos_data[id_plot], adj_sic_data[id_plot], psym=3
  ; Plot fit line
  if n_fit_point gt 0 then begin
    oplot, fit_tos_data, fit_sic_data, color=color_fit, thick=thick
  endif
endif

; Plot time series of monthly mean ice cover
if keyword_set( time_data ) and keyword_set( lat_data ) $
    and ( keyword_set( in_sic_data ) or keyword_set( adj_sic_data ) ) then begin
  ; Determine the latitude spacing, assuming a regular polar grid
  temp = lat_data[1:n_lat-1] - lat_data[0:n_lat-2]
  d_lat = mean( temp )
  if ( max( temp / d_lat ) gt 1.05 ) or ( min( temp / d_lat ) lt 0.95 ) $
      then stop
  ; Calculate surface area
  temp_weight = 0. * fltarr( n_lon, n_lat, n_time )
  for i_lat = 0, n_lat - 1 do begin
    temp_weight[*,i_lat,*] = cos( !pi / 180. * lat_data[i_lat] )
  endfor
  temp_weight = a_earth * temp_weight / total( temp_weight[*,*,0] ) $
      * n_lat * d_lat / 180.
  temp_weight = reform( temp_weight, n_lon * n_lat, n_time )
  ; Extract input and adjusted SICs
  if keyword_set( in_sic_data ) then begin
    in_sic_data_mean = reform( in_sic_data, n_lon * n_lat, n_time ) / sic_full
    in_sic_data_mean = total( in_sic_data_mean * temp_weight, 1, nan=1 )
  endif
  if keyword_set( adj_sic_data ) then begin
    adj_sic_data_mean = reform( adj_sic_data, n_lon * n_lat, n_time ) / sic_full
    adj_sic_data_mean = total( adj_sic_data_mean * temp_weight, 1, nan=1 )
  endif
  ; Define plotting time vector
  time_plot = convert_time_format( time_data, 'yyyymmdd', 'decimal year' )
  ; Determine plotting range
  yrange_sic_mean = [ 0, 0 ]
  if keyword_set( in_sic_data_mean ) then begin
    yrange_sic_mean[1] = max( [ yrange_sic_mean[1], max( in_sic_data_mean ) ] )
  endif
  if keyword_set( adj_sic_data_mean ) then begin
    yrange_sic_mean[1] = max( [ yrange_sic_mean[1], max( adj_sic_data_mean ) ] )
  endif
  ; Set up plotting window
  !p.multi[0] = 3 * 2
  temp_xtitle = 'Input=dotted, adjusted=solid'
  plot, time_plot, 0*time_plot, xthick=thick, ythick=thick, $
      yrange=yrange_sic_mean, ystyle=0, xstyle=1, charsize=charsize, $
      font=font, xtitle=temp_xtitle, ytitle='Mm!E2!N', xticklen=0.04, $
      title='Monthly total sea ice', nodata=1
  ; Plot time series
  if keyword_set( in_sic_data_mean ) then begin
    oplot, time_plot, in_sic_data_mean, thick=thick, color=4, linestyle=1
  endif
  if keyword_set( adj_sic_data_mean ) then begin
    oplot, time_plot, adj_sic_data_mean, thick=thick, color=4, linestyle=0
  endif
endif

; Prepare to plot maps of SIC for one month
if keyword_set( time_data ) and keyword_set( lat_data ) $
    and keyword_set( lon_data ) $
    and ( keyword_set( in_sic_data ) or keyword_set( adj_sic_data ) ) then begin
  ; Select date
  id_select = floor( randomu( seed, 1 ) * n_time )
  temp = fix( strmid( time_data[id_select], 4, 2 ) ) - 1
  title_date = month_name( temp ) + ' ' + strmid( time_data[id_select], 0, 4 )
  ; Define colours and levels
  color = 0
  n_level = 11
  levels = [ findgen( n_level ) / ( n_level - 1. ) ] * sic_full
  c_colors = [ 2 + n_level - 1, 2 + indgen( n_level - 1 ) ]
  colorbrewer, 'YlGnBu', n_level - 2, start_index=2, $
      add_rgb=[[100,100,100],[255,255,255]]
  ; Determine the centre of the map
  if ( min( lat_data, max=temp ) gt -1. ) and ( temp gt 70. ) then begin
    north_opt = 1
    centre_map = [ 0., 90. ]
  endif else if ( max( lat_data, min=temp ) lt 1. ) and ( temp lt -70. ) $
      then begin
    south_opt = 1
    centre_map = [ 0., -90. ]
  endif else if ( max( lon ) - min( lon ) gt 350. ) $
      and ( max( lat ) - min( lat ) gt 170. ) then begin
    centre_map = [ 10., 10. ]
    hammer_opt = 1
  endif else begin
    temp = fltarr( n_lon, n_lat, 2 )
    for i_lon = 0, n_lon - 1 do temp[i_lon,*,0] = lon_data[i_lon]
    for i_lat = 0, n_lat - 1 do temp[*,i_lat,1] = lat_data[i_lat]
    temp = transpose( reform( temp, n_lon * n_lat, 2 ) )
    centre_map = geo_mean( temp, minimise_maxdist=1 )
    temp = 0
    hammer_opt = 1
  endelse
  ; Plot legend
  !p.multi[0] = 3 * 2 - 1
  contour_legend, [1.03,1.07], [0.2,0.8], vertical=1, levels=levels, $
      c_colors=c_colors, color=color, thick=thick, charsize=charsize/2., $
      font=0, normal=1, subtitle='%'
endif

; Plot instantaneous map of input SIC
if keyword_set( time_data ) and keyword_set( lat_data ) $
    and keyword_set( lon_data ) and keyword_set( in_sic_data ) then begin
  ; Set up plotting window
  !p.multi[0] = 3 * 2 - 1
  title = 'Input SIC for ' + title_date
  ; Plot map
  temp_data = in_sic_data[*,*,id_select]
  id = where( finite( temp_data ) eq 0, n_id )
  if n_id gt 1 then temp_data[id] = 2.5 * sic_full
  contour_world, temp_data, lon_data, lat_data, levels=levels, $
      c_colors=c_colors, color=color, charcolor=color, thick=thick, $
      charsize=charsize, font=0, blockcell=1, noerase=1, title=title, $
      north=north_opt, south=south_opt, centre=centre_map, horizon=1, $
      nolines=1, noborder=1
endif

; Plot instantaneous map of adjusted SIC
if keyword_set( time_data ) and keyword_set( lat_data ) $
    and keyword_set( lon_data ) and keyword_set( adj_sic_data ) then begin
  ; Set up plotting window
  !p.multi[0] = 3 * 2 - 2
  title = 'Adjusted SIC for ' + title_date
  ; Plot map
  temp_data = adj_sic_data[*,*,id_select]
  id = where( finite( temp_data ) eq 0, n_id )
  if n_id gt 1 then temp_data[id] = 2.5 * sic_full
  contour_world, temp_data, lon_data, lat_data, levels=levels, $
      c_colors=c_colors, color=color, charcolor=color, thick=thick, $
      charsize=charsize, font=0, blockcell=1, noerase=1, title=title, $
      north=north_opt, south=south_opt, centre=centre_map, horizon=1, $
      nolines=1, noborder=1, hammer=hammer_opt
endif

; Close postscript output
ps_close

;***********************************************************************
; The end

return
END
