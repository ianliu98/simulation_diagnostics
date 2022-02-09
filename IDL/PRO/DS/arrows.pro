;+
; NAME:
;    ARROWS
;
; PURPOSE:
;    This function plots arrows onto the current window.
;
; CATEGORY:
;    Graphics
;
; CALLING SEQUENCE:
;    arrows, TAIL_X, TAIL_Y, BODY_LENGTH, POINT_ANGLE
;
; INPUTS:
;    TAIL_X:  The x-coordinate of the tail of the arrow, of type floating 
;        point.  This can be a vector.
;    TAIL_Y:  The y-coordinate of the tail of the arrow, of type floating 
;        point.  This can be a vector.
;    BODY_LENGTH:  The length of the arrow, of type floating point.  This can 
;        be a vector.  This input is ignored if the HEADX and HEADY keyword 
;        parameters are inputted.
;    POINT_ANGLE:  The angle of the arrow direction in counterclockwise degrees 
;        from the positive x-axis, of type floating point.  This input is 
;        ignored if the HEADX and HEADY keyword parameters are inputted.  This 
;        can be a vector.
;
; KEYWORD PARAMETERS:
;    ABSHEADLENGTH:  The length of the head fins in the units of the x-axis.  
;        If given this overrides HEADLENGTH.  The default is not given.
;    COLOR:  An index value for the arrow color.  This can be a vector.
;    HEADANGLE:  The angle of the head fins in degrees, of type floating 
;        point.  The default is 45 degrees.
;    HEADLENGTH:  The length of the head fins as a fraction of arrow length, of 
;        type floating point.  The default is 0.25.  Can be overriden by 
;        ABSHEADLENGTH.
;    HEADX:  The x-coordinate of the tail of the arrow, of type floating 
;        point.  Use of HEADX and HEADY overrides the BODY_LENGTH and 
;        POINT_ANGLE inputs.  This can be a vector.
;    HEADY:  The y-coordinate of the tail of the arrow, of type floating
;        point.  Use of HEADX and HEADY overrides the BODY_LENGTH and 
;        POINT_ANGLE inputs.  This can be a vector.
;    LINESTYLE:  The line style index for drawing the arrow.
;    NOCLIP:  The IDL graphics NOCLIP option.
;    NORMAL:  If set, values are taken in normal coordinates.  The default is 
;        data coordinates.
;    SOLID:  If set, the head fin is filled.  The default is no fill.
;    TAIL:  If set, tail fins are also plotted, according to the parameters set 
;        for the head fins.
;    THICK:  The thickness of the lines forming the arrow, of type integer.
;        This can be a vector.
;
; USES:
;    constants.pro
;
; PROCEDURE:
;    This function overplots an arrow over the current plot.
;
; EXAMPLE:
;    Make a plot.
;      plot, [0,1], /nodata
;    Plot an arrow from the middle of the plot to halfway to the upper
;    right corner.
;      arrows, 0.5, 0.5, 0.25, 45.
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone (stoned@atm.ox.ac.uk), 2003-06-06.
;    Modified:  DAS, 2012-07-18 (removed bug which demanded LENGTH when HEADX 
;        and HEADY are input;  switched to default of clipping;  altered 
;        formating style;  re-named X to TAIL_X and Y to TAIL_Y, LENGTH to 
;        BODY_LENGTH, ANGLE to POINT_ANGLE)
;    Modified:  DAS, 2018-09-06 (Fixed error in setting default head fin length)
;-

;***********************************************************************

PRO ARROWS, $
    TAIL_X, TAIL_Y, BODY_LENGTH, POINT_ANGLE, $
    ABSHEADLENGTH=head_length_abs, HEADANGLE=head_angle, $
      HEADLENGTH=head_length_rel, $
    COLOR=color, $
    HEADX=head_x, HEADY=head_y, $
    LINESTYLE=linestyle, $
    NOCLIP=no_clip_opt, $
    NORMAL=normal_opt, $
    SOLID=solid_opt, $
    TAIL=tail_opt, $
    THICK=thick

;***********************************************************************
; Constants and Options

; Absolute constants
constants, degrad=degrad

; Option for solid filling head (and tail) fin
solid_opt = keyword_set( solid_opt )

; Option for plotting tail fins
tail_opt = keyword_set( tailopt )

; Number of arrows to plot
n_arrow = n_elements( tail_x )
if n_arrow ne n_elements( tail_y ) then stop

; Set head fin length
if not( keyword_set( head_length_abs ) ) $
    and not( keyword_set( head_length_rel ) ) then begin
  head_length_rel = 0.25
endif
; Set head fin angle
if n_elements( head_angle ) eq 0 then head_angle = 45.

; Determine the ratio of the x-y axis ranges in device units
xy_ratio = ( !y.crange[1] - !y.crange[0] ) / ( !p.clip[3] - !p.clip[1] ) $
    / ( ( !x.crange[1] - !x.crange[0] ) / ( !p.clip[2] - !p.clip[0] ) )

;***********************************************************************
; Initialise Variables

; Reform BODY_LENGTH into vector (vectors are denoted with the "vec" prefix)
if n_elements( body_length ) ne 0 then begin
  if n_arrow eq n_elements( body_length ) then begin
    vec_length = body_length
  endif else begin
    vec_length = body_length[0] + fltarr( n_arrow )
  endelse
endif

; Reform POINT_ANGLE into vector
if n_elements( point_angle ) ne 0 then begin
  if n_arrow eq n_elements( point_angle ) then begin
    vec_angle = point_angle
  endif else begin
    vec_angle = point_angle[0] + fltarr( n_arrow )
  endelse
endif

; Calculate head positions
if ( n_elements( head_x ) eq 0 ) or ( n_elements( head_y ) eq 0 ) then begin
  vec_head_x = tail_x + vec_length * cos( degrad * vec_angle )
  vec_head_y = tail_y + vec_length * sin( degrad * vec_angle ) * xy_ratio
endif else begin
  ; Reform HEAD_X into vector
  if n_arrow eq n_elements( head_x ) then begin
    vec_head_x = head_x
  endif else begin
    vec_head_x = head_x[0] + fltarr( n_arrow )
  endelse
  ; Reform HEAD_Y into vector
  if n_arrow eq n_elements( head_y ) then begin
    vec_head_y = head_y
  endif else begin
    vec_head_y = head_y[0] + fltarr( n_arrow )
  endelse
  ; Calculate vector length
  vec_length = sqrt( ( vec_head_x - tail_x ) ^ 2 $
      + ( ( vec_head_y - tail_y ) / xy_ratio ) ^ 2 )
  ; Calculate vector angle
  vec_angle = atan( ( vec_head_y - tail_y ) / xy_ratio, vec_head_x - tail_x ) $
      / degrad
endelse

; Reform COLOR into vector
if n_elements( color ) ne 0 then begin
  if n_arrow eq n_elements( color ) then begin
    vec_color = color
  endif else begin
    vec_color = color[0] + intarr( n_arrow )
  endelse
endif else begin
  vec_color = !p.color + intarr( n_arrow )
endelse

; Reform LINESTYLE into vector
if n_elements( linestyle ) ne 0 then begin
  if n_arrow eq n_elements( linestyle ) then begin
    vec_linestyle = linestyle
  endif else begin
    vec_linestyle = linestyle[0] + intarr( n_arrow )
  endelse
endif else begin
  vec_linestyle = !p.linestyle + intarr( n_arrow )
endelse

; Reform THICK into vector (vectors are denoted with the "v" prefix)
if n_elements( thick ) ne 0 then begin
  if n_arrow eq n_elements( thick ) then begin
    vec_thick = thick
  endif else begin
    vec_thick = thick[0] + fltarr( n_arrow )
  endelse
endif else begin
  vec_thick = !p.thick + intarr( n_arrow )
endelse

; Calculate head fin end points
if keyword_set( head_length_abs ) then begin
  head_length = head_length_abs
endif else begin
  head_length = head_length_rel * vec_length
endelse
vec_head_x_1 = vec_head_x + head_length $
           * cos( degrad * ( vec_angle - 180 - head_angle ) )
vec_head_y_1 = vec_head_y + head_length $
           * sin( degrad * ( vec_angle - 180 - head_angle ) ) * xy_ratio
vec_head_x_2 = vec_head_x + head_length $
           * cos( degrad * ( vec_angle - 180 + head_angle ) )
vec_head_y_2 = vec_head_y + head_length $
           * sin( degrad * ( vec_angle - 180 + head_angle ) ) * xy_ratio

; Calculate tail fin end points
if tail_opt then begin
  vec_tail_x_1 = tail_x + head_length $
             * cos( degrad * ( vec_angle - 180 - head_angle ) )
  vec_tail_y_1 = tail_y + head_length $
             * sin( degrad * ( vec_angle - 180 - head_angle ) ) * xy_ratio
  vec_tail_x_2 = tail_x + head_length $
             * cos( degrad * ( vec_angle - 180 + head_angle ) )
  vec_tail_y_2 = tail_y + head_length $
             * sin( degrad * ( vec_angle - 180 + head_angle ) ) * xy_ratio
endif

;***********************************************************************
; Plot 

; Iterate through arrows
for i_arrow = 0, n_arrow - 1 do begin
  ; Plot main line
  plots, [tail_x[i_arrow],vec_head_x[i_arrow]], $
      [tail_y[i_arrow],vec_head_y[i_arrow]], color=vec_color[i_arrow], $
      thick=vec_thick[i_arrow], normal=normal_opt, $
      linestyle=vec_linestyle[i_arrow], noclip=no_clip_opt
  ; Plot head fins
  if solid_opt then begin
    temp_x = [ vec_head_x[i_arrow], vec_head_x_1[i_arrow], $
        vec_head_x_2[i_arrow], vec_head_x[i_arrow] ]
    temp_y = [ vec_head_y[i_arrow], vec_head_y_1[i_arrow], $
        vec_head_y_2[i_arrow], vec_head_y[i_arrow] ]
    polyfill, temp_x, temp_y, color=vec_color[i_arrow], normal=normal_opt, $
        noclip=no_clip_opt
  endif else begin
    plots, [vec_head_x[i_arrow],vec_head_x_1[i_arrow]], $
        [vec_head_y[i_arrow],vec_head_y_1[i_arrow]], color=vec_color[i_arrow], $
        thick=vec_thick[i_arrow], normal=normal_opt, noclip=no_clip_opt
    plots, [vec_head_x[i_arrow],vec_head_x_2[i_arrow]], $
        [vec_head_y[i_arrow],vec_head_y_2[i_arrow]], color=vec_color[i_arrow], $
         thick=vec_thick[i_arrow], normal=normal_opt, noclip=no_clip_opt
  endelse
  ; Plot tail fins
  if tail_opt then begin
    if solid_opt then begin
      temp_x = [ tail_x[i_arrow], vec_tail_x_1[i_arrow], $
          vec_tail_x_2[i_arrow], tail_x[i_arrow] ]
      temp_y = [ tail_y[i_arrow], vec_tail_y_1[i_arrow], $
          vec_tail_y_2[i_arrow], tail_y[i_arrow] ]
      polyfill, temp_x, temp_y, color=vec_color[i_arrow], normal=normal_opt, $
          noclip=no_clip_opt
    endif else begin
      plots, [tail_x[i_arrow],vec_tail_x_1[i_arrow]], $
          [tail_y[i_arrow],vec_tail_y_1[i_arrow]], color=vec_color[i_arrow], $
          thick=vec_thick[i_arrow], normal=normal_opt, noclip=no_clip_opt
      plots, [tail_x[i_arrow],vec_tail_x_2[i_arrow]], $
          [tail_y[i_arrow],vec_tail_y_2[i_arrow]], color=vec_color[i_arrow], $
          thick=vec_thick[i_arrow], normal=normal_opt, noclip=no_clip_opt
    endelse
  endif
endfor

;***********************************************************************
; The end

return
END
