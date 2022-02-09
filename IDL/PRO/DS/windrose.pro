;+
; NAME:
;    windrose
;
; PURPOSE:
;    This function plots a wind rose given two-dimensional data.
;
; CATEGORY:
;    Graphics
;
; CALLING SEQUENCE:
;    result = windrose( DATA )
;
; INPUTS:
;    DATA:  A required floating point array providing the two-dimensional data 
;        to be analysed.  Of size N_DATA*2, where N_DATA is the total number of 
;        data points, with 2 values per data point.  If the XY option is set 
;        then these two values are the Cartesian x- and y-dimensions, 
;        otherwise they are assumed to be angle and magnitude.
;    C_COLORS, CENTRE, FILTER_WIDTH, FILTER_WINDOW, HIST_LP, LEVELS, N_PETAL, 
;      PETAL, SCALE_FACTOR, SCALE_RADIUS, WIDTH_PETAL
;
; KEYWORD PARAMETERS:
;    AREA_NORMALISE:  If set, then the plotted density is normalised such that 
;        the area covered by a specific level is proportional only to the 
;        actual data density, and not to the radius.  The default is to plot 
;        such that the length from the centre of the wind rose is proportional 
;        to the data density.
;    C_COLORS:  An optional integer vector of length N_LEVELS defining the 
;        colour indices corresponding to the contour levels listed in LEVELS.  
;        The default is 2+indgen(N_LEVELS).
;    CENTRE:  An optional two-element floating point vector defining the x- and 
;        y-coordinate of the centre of the wind rose.  The default is [0,0].
;    FILTER_WIDTH:  An optional scalar integer specifying the width of the 
;        smoothing filter to apply to the wind rose in the angular direction.  
;        The default is to have no smoothing.  See filter.pro for more details.
;    FILTER_WINDOW:  An optional scalar string specifying the smoothing 
;        function to use when applying a smoothing filter in the angular 
;        direction.  This is ignored unless FILTER_WIDTH is set.  See 
;        filter.pro for more details, including the default.
;    HIST_LP:  An optional floating point array of size N_LEVELS*N_PETAL 
;        directly providing the wind rose histogram to be plotted.  This may 
;        have been the output RESULT of a previous call to windrose.pro, for 
;        instance.
;    LEVELS:  An optional floating point vector defining the contour levels 
;        along the petals.  The default is generated using choose_levels.pro.
;    N_PETAL:  An optional integer describing the number of petals to draw.  
;        The default is the length of PETAL, if input, or otherwise 8 
;        (i,e. N, NE, E, SE, S, SW, W, NW).
;    NO_PLOT:  If set then no plot is produced.  The default is to plot the 
;        output to the plotting device.
;    NOCLIP:  The IDL noclip plotting keyword option.  Set to 0 to ensure 
;        clipping.  The default is no clipping (i.e. NOCLIP=1).
;    NORMAL:  If set, then CENTRE and SCALE_RADIUS are interpreted as being in 
;        normal coordinates.  The default is to assume data coordinates.
;    OVERPLOT:  If set then the wind rose is plotted on top of the current 
;        plot, using device coordinates.  The default is to set up a new 
;        plotting window.
;    PETAL:  An optional floating point vector defining the angular edges of 
;        the petals.  Of length N_PETAL.
;    SCALE_FACTOR:  An optional floating point scalar setting the scaling 
;        factor for converting the wind rose histogram density into petal 
;        lengths.  All petal lengths are directly multiplied by this factor.  A 
;        default can be estimated by the code.
;    SCALE_RADIUS:  An optional floating point scalar defining the reference 
;        target radius for scaling of the petal lengths.  This would be the 
;        radius of a circular rose.  The default is half of 
;        !x.crange[1]-!x.crange[0], or if NORMAL is set then 0.5.
;    UNIFY_FIRST_LEVEL:  If set, then the first level is represented by a round 
;        disk in the centre of the wind rose with no directional indication.  
;        The default is to include the first level as part of the petals.
;    WIDTH_PETAL:  An optional floating point scalar between 0 and 1 
;        specifying the angular width of the plotted petal in units of the full 
;        angular range coverd by the petal.  The default is 1 (full range).
;    XY:  If set, then the two values for each data point in DATA are 
;        interpreted as the Cartesian x- and y-dimensions.
;
; OUTPUTS
;    RESULT:  Returns a two-dimensional floating point array containing the 
;        counts for DATA for each of the N_LEVELS levels and N_PETAL petals of 
;        the wind rose.  Of size N_LEVELS*N_PETAL.  If HIST_LP is input, then 
;        RESULT simply returns HIST_LP.
;
; USES:
;    choose_levels.pro
;    filter.pro
;
; PROCEDURE:
;    This function calculates two-dimensional angle-radius histograms and then 
;    plots those histograms using colour density levels.
;
; EXAMPLE:
;    n_point = 10000l
;    data_ang = randomu( 1, n_point ) * 2. * !pi
;    data_mag = abs( randomn( 2, n_point ) + 5. )
;    wr_levels = findgen( 11 ) * 2.
;    n_wr_petal = 160
;    wr_petal = ( findgen( n_wr_petal ) + 0.5 ) / float( n_wr_petal ) * 2. * !pi
;    colorbrewer, 'PiYG', 10, start_index=2
;    temp = windrose( [[data_ang],[data_mag]], levels=wr_levels, petal=wr_petal, centre=[0.5,0.5] )
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone (dastone@runbox.com), 2018-09-28
;    Modified:  DAS, 2018-10-08 (Added SCALE_FACTOR and WIDTH_PETAL keyword 
;        inputs and UNIFY_FIRST_LEVEL option;  Completed documentation)
;-

;***********************************************************************

FUNCTION WINDROSE, $
    DATA, $
    CENTRE=centre, $
    FILTER_WIDTH=filter_width, FILTER_WINDOW=filter_window, $
    LEVELS=levels, C_COLORS=c_colors, $
    PETAL=petal, N_PETAL=n_petal, WIDTH_PETAL=petal_width, $
    SCALE_FACTOR=scale_factor, SCALE_RADIUS=scale_radius, $
    HIST_LP=hist_lp, $
    NO_PLOT=no_plot_opt, OVERPLOT=overplot_opt, $
    NOCLIP=noclip_opt, $
    NORMAL=normal_opt, $
    UNIFY_FIRST_LEVEL=unify_first_level_opt, $
    AREA_NORMALISE=area_normalise_opt, $
    XY=xy_opt

;***********************************************************************
; Constants and Options

; The default centre
if not( keyword_set( centre ) ) then begin
  if keyword_set( normal_opt ) then begin
    centre = [ 0.5, 0.5 ]
  endif else begin
    centre = [ 0., 0. ]
  endelse
endif
if n_elements( centre ) ne 2 then stop

; The number of levels to use
n_levels = n_elements( levels )

; The default number of petals
n_petal = n_elements( petal )
if n_petal eq 0 then n_petal = 8
; The default petal boundary angles
if not( keyword_set( petal ) ) then begin
  petal_use = ( findgen( n_petal ) + 0.5 ) / n_petal * 2. * !pi
endif else begin
  petal_use = petal
endelse
; Add an extra wrapping value for convenience
petal_use = [ petal_use[n_petal-1], petal_use ]

; The default fractional petal width
if not( keyword_set( petal_width ) ) then petal_width = 1.

; Ensure the data array is input
if n_elements( data ) eq 0 then begin
  ; If it is not then check whether the pre-calculated histogram is input
  if n_elements( hist_lp ) eq 0 then begin
    stop
  endif else begin
    ; Note that there are no data points
    n_data = 0
    ; Ensure that the dimensions are input
    if n_levels eq 0 then stop
    if n_levels ne n_elements( hist_lp[*,0] ) then stop
    if n_petal eq 0 then stop
    if n_petal ne n_elements( hist_lp[0,*] ) then stop
  endelse
endif else begin
  if n_elements( data[0,*] ) ne 2 then stop
  ; Number of data points
  n_data = n_elements( data[*,0] )
endelse

; Option to normalise radii by area
area_normalise_opt = keyword_set( area_normalise_opt )
; Option to unify the first level from all petals into a disk
unify_first_level_opt = keyword_set( unify_first_level_opt )
; Note that the unifying option is not currently implemented with 
; area_normalise_opt
if area_normalise_opt + unify_first_level_opt eq 2 then stop

;***********************************************************************
; Calculate 2-dimensional histogram

; If this is required
if n_data ne 0 then begin

  ; If we first need to convert to angle-magnitude format
  if keyword_set( xy_opt ) then begin
    data_angle = atan( data[*,1], data[*,0] )
    data_mag = sqrt( ( data[*,0] ^ 2. ) + ( data[*,1] ^ 2. ) )
  ; Otherwise copy the input data
  endif else begin
    data_angle = data[*,0]
    data_mag = data[*,1]
  endelse
  ; Ensure angles in the 0 to 2pi range
  id = where( data_angle lt 0, n_id )
  if n_id gt 0 then data_angle[id] = data_angle[id] + 2. * !pi
  while max( data_angle ) gt 2. * !pi do begin
    id = where( data_angle gt 2. * !pi, n_id )
    if n_id gt 0 then data_angle[id] = data_angle[id] - 2. * !pi
  endwhile

  ; Determine default levels
  if not( keyword_set( levels ) ) then begin
    levels_use = choose_levels( [0., data_mag ] )
    if ( min( levels_use ) lt 0. ) and ( min( data_mag ) ge 0. ) then begin
      id = where( levels_use gt 0. )
      levels_use = [ 0., levels_use[id] ]
    endif
  endif else begin
    levels_use = levels
  endelse
  n_levels = n_elements( levels_use )
  ; Add a maximum value for convenience
  if max( data_mag ) ge levels_use[n_levels-1] then begin
    levels_use = [ levels_use, max( data_mag ) + 1. ]
  endif else begin
    levels_use = [ levels_use, levels_use[n_levels-1] + 1. ]
  endelse

  ; Initialise histogram array
  hist = fltarr( n_levels, n_petal )
  ; Iterate through petals
  for i_petal = 0, n_petal - 1 do begin
    ; Find all values within the angle covered by this petal
    if petal_use[i_petal] lt petal_use[i_petal+1] then begin
      id_petal = where( ( data_angle ge petal_use[i_petal] ) $
          and ( data_angle lt petal_use[i_petal+1] ), n_id_petal )
    endif else begin
      id_petal = where( ( data_angle ge petal_use[i_petal] ) $
          or ( data_angle lt petal_use[i_petal+1] ), n_id_petal )
    endelse
    if n_id_petal gt 0 then begin
      ; Iterate through levels
      for i_levels = 0, n_levels - 1 do begin
        ; Determine how many values are in this angle and magnitude bin
        id_level = where( ( data_mag[id_petal] ge levels_use[i_levels] ) $
            and ( data_mag[id_petal] lt levels_use[i_levels+1] ), n_id_level )
        hist[i_levels,i_petal] = n_id_level
      endfor
    endif
  endfor

; Otherwise copy input histogram
endif else begin
  hist = hist_lp
  levels_use = levels
endelse

; Apply smoothing filter if requested
if keyword_set( filter_width ) then begin
  ; The default filter window
  if not( keyword_set( filter_window ) ) then filter_window = 'boxcar'
  ; Iterate through levels
  for i_levels = 0, n_levels - 1 do begin
    hist[i_levels,*] = reform( filter( reform( hist[i_levels,*] ), $
        filter_width, filter_window, wrap_edges=1 ), 1, n_petal )
  endfor
endif

;***********************************************************************
; Plot wind rose

; Do not do this if requested not to
if not( keyword_set( no_plot_opt ) ) then begin

  ; The default colour indices for the levels
  if not( keyword_set( c_colors ) ) then c_colors = 2 + indgen( n_levels )

  ; If we need to set up the plotting window
  if not( keyword_set( overplot_opt ) ) then begin
    plot, [-1,1], [-1,1], xstyle=5, ystyle=5, nodata=1, isotropic=1
  endif

  ; Determine the reference radius for scaling the rose
  if not( keyword_set( scale_factor ) ) then begin
    if area_normalise_opt eq 1 then begin
      density_max = sqrt( total( hist ) / !pi ) * 2.
    endif else begin
      density_max = total( hist ) / n_petal * 2.
    endelse
    if not( keyword_set( scale_radius ) ) then begin
      if keyword_set( normal_opt ) then begin
        scale_radius = 0.5
      endif else begin
        scale_radius = ( !x.crange[1] - !x.crange[0] ) / 2.
      endelse
    endif
    scale_factor = scale_radius / density_max
  endif

  ; Iterate through petals
  for i_petal = 0, n_petal - 1 do begin
    ; Calculate the angle covered by this petal
    petal_angle = petal_use[i_petal+1] - petal_use[i_petal]
    if petal_angle lt 0 then petal_angle = petal_angle + 2. * !pi
    ; Iterate through levels
    for i_levels = 0, n_levels - 2 do begin
      ; Determine the area inside of this section
      if i_levels eq 0 then begin
        area_inside = 0
      endif else if unify_first_level_opt eq 0 then begin
        area_inside = total( hist[0:i_levels-1,i_petal] )
      endif else if i_levels eq 1 then begin
        area_inside = total( hist[0,0:n_petal-1] ) / n_petal
      endif else begin
        area_inside = total( hist[1:i_levels-1,i_petal] ) $
            + total( hist[0,0:n_petal-1] ) / n_petal
      endelse
      ; Determine the area inside of and including this section
      if ( i_levels eq 0 ) and ( unify_first_level_opt eq 1 ) then begin
        area_include = total( hist[0,0:n_petal-1] ) / n_petal
      endif else if unify_first_level_opt eq 0 then begin
        area_include = total( hist[0:i_levels,i_petal] )
      endif else begin
        area_include = total( hist[1:i_levels,i_petal] ) $
            + total( hist[0,0:n_petal-1] ) / n_petal
      endelse
      ; Determine the radii corresponding to these areas
      if area_normalise_opt eq 1 then begin
        radius_inside = sqrt( 2. * area_inside / petal_angle )
        radius_outside = sqrt( 2. * area_include / petal_angle )
      endif else begin
        radius_inside = area_inside
        radius_outside = area_include
      endelse
      ; Normalise the radii
      radius_inside = radius_inside * scale_factor
      radius_outside = radius_outside * scale_factor
      ; Determine the petal width
      temp_petal = petal_use[i_petal:i_petal+1]
      if ( petal_width ne 1 ) $
          and ( ( unify_first_level_opt ne 1 ) or ( i_levels ne 0 ) ) $
          then begin
        if temp_petal[0] gt temp_petal[1] then begin
          temp_petal[0] = temp_petal[0] - 2. * !pi
        endif
        temp_petal = ( temp_petal[0] + temp_petal[1] ) / 2. $
            + [ -0.5, 0.5 ] * petal_width $
            * ( temp_petal[1] - temp_petal[0] )
      endif
      ; Plot the section
      temp_x = [ radius_inside * cos( temp_petal[0] ), $
          radius_inside * cos( temp_petal[1] ), $
          radius_outside * cos( temp_petal[1] ), $
          radius_outside * cos( temp_petal[0] ), $
          radius_inside * cos( temp_petal[0] ) ] $
          + centre[0]
      temp_y = [ radius_inside * sin( temp_petal[0] ), $
          radius_inside * sin( temp_petal[1] ), $
          radius_outside * sin( temp_petal[1] ), $
          radius_outside * sin( temp_petal[0] ), $
          radius_inside * sin( temp_petal[0] ) ] $
          + centre[1]
      polyfill, temp_x, temp_y, color=c_colors[i_levels], noclip=noclip_opt, $
          normal=normal_opt
    endfor
  endfor

endif

;***********************************************************************
; The end

;stop
return, hist
END
