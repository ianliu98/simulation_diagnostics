;+
; NAME:
;    GEO_MEAN
;
; PURPOSE:
;    This function returns the coordinates of the central geographical point
;    of the input points.
;
; CATEGORY:
;    Geographical
;
; CALLING SEQUENCE:
;    Result = geo_mean( POINTS )
;
; INPUTS:
;    POINTS:  A floating point array of size 2*N_POINTS containing the 
;        longitude and latitude coordinates of the N_POINTS points.  Longitude 
;        values are in elements [0,*] and latitude values are in elements [1,*].
;
; KEYWORD PARAMETERS:
;    BUFFER_LIMIT:  A floating-point scalar containing size, in degrees, of a 
;        buffer area to add to the LIMIT rectangle.
;    LIMIT:  If set, then this returns the LIMIT keyword in the 8-element 
;        format used by map_set.pro.  The values are defined as the smallest 
;        rectangle that contains all of the points in POINTS with the 
;        restriction that the North Pole is directly above the central point.
;    MAXDIST:  If MINIMISE_MAXDIST is set, or if 1 is input to MAXDIST, then 
;        this returns the maximum distance between the central point and the 
;        points input in POINTS.  Of type floating point.  Units are km.
;    MEDIAN:  If set then the function calculates the median point.
;    MINIMISE_MAXDIST:  If set then the function estimates the central point 
;        that minimises the distance to the most distant point.
;    WEIGHT:  A floating point vector of length N_POINTS containing weights 
;        for the inputed POINTS.  This can include area-weighting if the points 
;        represent cells on a regular polar grid.
;
; OUTPUTS:
;    Result:  Returns a vector of [ longitude, latitude ] coordinates of the 
;        central point.  The default is to return the mean point.  See the 
;        MEDIAN and MINIMISE_MAXDIST keywords to use different definitions of 
;        the central point.
;    MAXDIST
;
; USES:
;    constants.pro
;    rotate_pole.pro
;
; PROCEDURE:
;    The function uses double precision floating point arithmetic to
;    calculate values.
;
; EXAMPLE:
;    Calculate the midpoint between Ottawa, ON (45.25N, 75.43W) and Victoria, 
;    BC (48.25N, 123.22W).
;      result = geo_mean( [ [ -75.43, 45.25 ], [ -123.22, 48.25 ] ] )
;    Compare to result from geo_mid.pro.
;      result = geo_mid( -75.43, 45.25, -123.22, 48.25 )
;
; MODIFICATION HISTORY:
;     Written by:  Daithi A. Stone (dstone@lbl.gov), 2015-01-15.
;-

;***********************************************************************

FUNCTION GEO_MEAN, $
    POINTS, $
    WEIGHT=weight, $
    LIMIT=limit, BUFFER_LIMIT=buffer_limit, $
    MEDIAN=median_opt, MINIMISE_MAXDIST=minimise_maxdist_opt

;***********************************************************************
; Define Constants

; Earth radius, degrees in a circle, degrees to radians conversion coefficient
constants, ndegree=n_degree, degrad=degrad, r_earth=r_earth
r_earth = r_earth / 1000.

; The number of inputed points
if n_elements( points[*,0] ) ne 2 then stop
n_points = n_elements( points[0,*] )

; Extract longitude and latitude values
lon = reform( points[0,*] )
lat = reform( points[1,*] )

; Determine the method
median_opt = keyword_set( median_opt )
minimise_maxdist_opt = keyword_set( minimise_maxdist_opt )
mean_opt = 1 - minimise_maxdist_opt

; Weighting option
weight_opt = keyword_set( weight )
if weight_opt eq 1 then weight_use = weight / total( weight )
; Weighting cannot currently be used with the median method
if ( weight_opt eq 1 ) and ( median_opt eq 1 ) then stop
; Weighting is not possible with the minimise_maxdist method
if ( weight_opt eq 1 ) and ( minimise_maxdist_opt eq 1 ) then stop

; The initial distance for the searching algorithm used for the $
; minimise_maxdist method (in degrees)
delta_try = 20.
; The threshold distance at which we say the search algorithm has done its 
; work (in degrees)
delta_try_thresh = 0.01

;***********************************************************************
; Calculations

; If we only have one point, the return that point
if n_points eq 1 then return, points[*,0]

; For the mean or median method
if ( mean_opt eq 1 ) or ( median_opt eq 1 ) then begin
  ; Convert to Cartesian coordinates
  x_val = cos( lon * degrad ) * cos( lat * degrad )
  y_val = sin( lon * degrad ) * cos( lat * degrad )
  z_val = sin( lat * degrad )
  ; Calculate Cartesian midpoint using the median method
  if median_opt eq 1 then begin
    ; If no weighting has been defined
    if weight_opt eq 0 then begin
      x_mid = median( x_val )
      y_mid = median( y_val )
      z_mid = median( z_val )
    ; If weighting has been defined
    endif else begin
      ; Not yet implemented
      stop
    endelse
  ; Calculate Cartesian midpoint using the mean method
  endif else if mean_opt eq 1 then begin
    ; If no weighting has been defined
    if weight_opt eq 0 then begin
      x_mid = mean( x_val )
      y_mid = mean( y_val )
      z_mid = mean( z_val )
    ; If weighting has been defined
    endif else begin
      x_mid = total( weight_use * x_val )
      y_mid = total( weight_use * y_val )
      z_mid = total( weight_use * z_val )
    endelse
  ; In case we missed a method
  endif else begin
    stop
  endelse
  ; Determine distance for Earth's centre
  dist = sqrt( x_mid ^ 2. + y_mid ^ 2. + z_mid ^ 2. )
  ; If we happen to be centred on the core of the planet
  if dist lt 0.01 then begin
    ; Use the lon=0, lat=0 point
    lon_mid = 0.
    lat_mid = 0.
  ; Otherwise
  endif else begin
    ; Convert to spherical midpoint
    lon_mid = atan( y_mid / x_mid ) / degrad
    if x_mid lt 0 then lon_mid = lon_mid - n_degree / 2.
    ;if lon_mid lt 0 then lon_mid = lon_mid + n_degree
    lat_mid = asin( z_mid / dist ) / degrad
  endelse
  ; Set to output
  point_mid = [ lon_mid, lat_mid ]
  ; Find the maximum distance to the central point to the input points
  if keyword_set( maxdist ) then begin
    maxdist = dblarr( n_points )
    for i_points = 0, n_points - 1 do begin
      maxdist[i_points] = geo_dist( point_mod[0], point_mid[1], $
          points[0,i_points], points[1,i_points] )
    endfor
    maxdist = max( maxdist )
  endif
endif

; Calculate Cartesian midpoint using the minimise_maxdist method
if minimise_maxdist_opt eq 1 then begin
points = double( points )
  ; Find the mean midpoint, to be used as a first guess
  point_guess = geo_mean( points, weight=weight_use )
  ; Keep doing this while we are still substantially revising our estimated 
  ; midpoint
  while delta_try gt delta_try_thresh do begin
    ; Define points to the north, east, south, and west of the this guess
    point_try = [ [ point_guess + [ 0, delta_try ] ], $
        [ point_guess + [ delta_try, 0 ] ], $
        [ point_guess + [ 0, -delta_try ] ], $
        [ point_guess + [ -delta_try, 0 ] ], $
        [ point_guess + [ delta_try, delta_try ] ], $
        [ point_guess + [ -delta_try, delta_try ] ], $
        [ point_guess + [ -delta_try, -delta_try ] ], $
        [ point_guess + [ delta_try, -delta_try ] ] ]
    n_try = n_elements( point_try[0,*] )
    id = where( point_try[1,*] gt n_degree / 4., n_id )
    for i_id = 0, n_id - 1 do begin
      point_try[0,id[i_id]] = -point_try[0,id[i_id]]
      point_try[1,id[i_id]] = 2. * n_degree / 4. - delta_try - point_guess[1]
    endfor
    id = where( point_try[1,*] lt -n_degree / 4., n_id )
    for i_id = 0, n_id - 1 do begin
      point_try[0,id[i_id]] = -point_try[0,id[i_id]]
      point_try[1,id[i_id]] = 2. * n_degree / 4. + delta_try - point_guess[1]
    endfor
    ; Calculate the maximum distance between the first guess and all of the 
    ; input points
    point_guess_maxdist = dblarr( n_points )
    for i_points = 0, n_points - 1 do begin
      point_guess_maxdist[i_points] = geo_dist( point_guess[0], $
          point_guess[1], points[0,i_points], points[1,i_points] )
    endfor
    point_guess_maxdist = max( point_guess_maxdist )
    ; Calculate the maximum distance between the search points and all of the 
    ; input points
    point_try_maxdist = dblarr( n_points, n_try )
    for i_try = 0, n_try - 1 do begin
      for i_points = 0, n_points - 1 do begin
        point_try_maxdist[i_points,i_try] = geo_dist( point_try[0,i_try], $
            point_try[1,i_try], points[0,i_points], points[1,i_points] )
       endfor
      point_try_maxdist[0,i_try] = max( point_try_maxdist[*,i_try] )
    endfor
    point_try_maxdist = reform( point_try_maxdist[0,*] )
    ; Make a new guess in the most promising direction
    id = where( point_try_maxdist eq min( point_try_maxdist ) )
    temp_weight = [ point_try_maxdist[id[0]], point_guess_maxdist ]
    point_guess = geo_mean( [ [ point_guess ], [ point_try[*,id[0]] ] ], $
        weight=temp_weight )
    ; Refine the search distance
    delta_try = abs( point_guess_maxdist - point_try_maxdist[id[0]] ) / degrad $
        / r_earth * 2.
  endwhile
  ; Output the estimate
  point_mid = point_guess
  ; Take the maximum distance from the central point
  maxdist = min( point_try_maxdist )
endif

; Ensure the longitude of the central point fits within standards
if point_mid[0] lt -n_degree / 2. then point_mid[0] = point_mid[0] + n_degree

; Estimate the smallest rectangle (with the North Pole north of the central 
; point) that encompasses all of the points
if keyword_set( limit ) then begin
  ; Determine the longitude of the new pole
  if point_mid[1] le 0 then begin
    lon_pole = point_mid[0]
  endif else begin
    lon_pole = point_mid[0] - n_degree / 2.
    if lon_pole lt -n_degree / 2. then lon_pole = lon_pole + n_degree
  endelse
  ; Determine the latitude of the new pole
  if point_mid[1] lt 0 then begin
    lat_pole = point_mid[1] + n_degree / 4.
  endif else begin
    lat_pole = n_degree / 4. - point_mid[1]
  endelse
  ; Calculate coordinates on rotated system
  pole_coord = [ lon_pole, lat_pole ]
  points_rotated = rotate_pole( points, pole_coord, double=1 )
  if point_mid[1] gt 0 then begin
    points_rotated[0,*] = points_rotated[0,*] + n_degree / 2.
  endif
  id = where( points_rotated[0,*] gt n_degree / 2., n_id )
  if n_id gt 0 then points_rotated[0,id] = points_rotated[0,id] - n_degree
  ; Find the longitude and latitude range
  range_lon_rotated = max( abs( points_rotated[0,*] ) ) * [ -1., 1. ]
  range_lat_rotated = max( abs( points_rotated[1,*] ) ) * [ -1., 1. ]
  if keyword_set( buffer_limit ) then begin
    range_lon_rotated = range_lon_rotated + [ -1., 1. ] * buffer_limit
    if range_lon_rotated[1] gt n_degree / 2. then begin
      range_lon_rotated = [ -1., 1. ] * n_degree / 2.
    endif
    range_lat_rotated = range_lat_rotated + [ -1., 1. ] * buffer_limit
    if range_lat_rotated[1] gt n_degree / 4. then begin
      range_lat_rotated = [ -1., 1. ] * n_degree / 4.
    endif
  endif
  ; If we are almost covering the full planet
  temp = round( range_lon_rotated[1] - range_lon_rotated[0] ) $
      + round( range_lat_rotated[1] - range_lat_rotated[0] )
  if temp ge 1.125 * n_degree then begin
    ; Then just make it a standard global map with the central point back on
    ; the Equator
    limit = [ 0., point_mid[0] - n_degree / 2., n_degree / 4., point_mid[0], $
        0, point_mid[0] + n_degree / 2., -n_degree / 4., point_mid[0] ]
    point_mid = [ point_mid[0], 0. ]
  ; Otherwise determine the best rectangle to cover
  endif else begin
    ; Convert these edges back to standard geographical coordinates
    limit_rotate = [ [ range_lon_rotated[0], 0. ], $
        [ 0., range_lat_rotated[1] ], [ range_lon_rotated[1], 0. ], $
        [ 0., range_lat_rotated[0] ] ]
    ;old_pole_rotated = rotate_pole( [ 0.d, n_degree / 4.d ], $
    ;    pole_coord, double=1 )
    old_pole_rotated = [ 0, n_degree / 4. - point_mid[1] ]
    limit = rotate_pole( limit_rotate, old_pole_rotated, double=1 )
    ;limit[0,*] = limit[0,*] - n_degree / 2. + point_mid[0]
    limit[0,*] = limit[0,*] + point_mid[0]
    id = where( limit[0,*] lt -n_degree / 2., n_id )
    if n_id gt 0 then limit[0,id] = limit[0,id] + n_degree
    temp_limit = limit
    limit[0,*] = temp_limit[1,*]
    limit[1,*] = temp_limit[0,*]
    limit = reform( limit, 8 )
  endelse
endif

;***********************************************************************
; The End

return, point_mid
END
