;+
; NAME:
;    ROTATE_POLE
;
; PURPOSE:
;    Convert geographical coordinates to a new polar coordinate system with a 
;    different North Pole.
;
; CATEGORY:
;    Geographical
;
; CALLING SEQUENCE:
;    Result = rotate_pole( POINTS_COORD, POLE_COORD )
;
; INPUT:
;    POINTS_COORD:  A floating-point array of size 2*N_POINTS containing the 
;        geographical coordinates in the standard longitude-latitude format.  
;        Elements [0,*] contain the longitude values and elements [1,*] contain 
;        the latitude values.  Units are in degrees.
;    POLE_COORD:  A 2-element floating-point array containing the longitude 
;        and latitude coordinates of the new North Pole, respectively.
;
; KEYWORD INPUTS:
;    DOUBLE:  If set, then calculations are performed in double precision.
;
; OUTPUT:
;    Result:  A floating-point array of size 2*N_POINTS containing the 
;        spherical coordinates under the new coordinate system.  In the same 
;        format as POINTS_COORD.
;
; EXAMPLES:
;    -
;
; MODIFICATION HISTORY:
;    Based on geo2mag.pro by:  Pascal Saint-Hilaire (Saint-Hilaire@astro.phys.ethz.ch), 2002-05
;    Written by:  Daithi A. Stone (dstone@lbl.gov), 2015-01-14
;-

;***********************************************************************

FUNCTION rotate_pole, $
    POINTS_COORD, $
    POLE_COORD, $
    DOUBLE=double_opt

;***********************************************************************
; Constanst and options

; Ensure proper format in the inputs
if n_elements( pole_coord ) ne 2 then stop
if n_elements( points_coord[*,0] ) ne 2 then stop
n_points = n_elements( points_coord[0,*] )

; Option for double precision
double_opt = keyword_set( double_opt )
if double_opt eq 1 then begin
  points_coord_use = double ( points_coord )
  pole_coord_use = double( pole_coord )
endif else begin
  points_coord_use = points_coord
  pole_coord_use = pole_coord
endelse

; Constants
if double_opt eq 1 then begin
  pi = !dpi
endif else begin
  pi = !pi
endelse
n_degree = 360.
degrad = 2. * pi / n_degree

;***********************************************************************
; Calculate coordinates on rotated grid

; Convert from degrees to radians
points_coord_use = points_coord_use * degrad
pole_coord_use = pole_coord_use * degrad

; Convert to Cartesian coordinates
points_cart = fltarr( 3, n_points )
if double_opt eq 1 then points_cart = double( points_cart )
points_cart[0,*] = cos( points_coord_use[0,*] ) * cos( points_coord_use[1,*] )
points_cart[1,*] = sin( points_coord_use[0,*] ) * cos( points_coord_use[1,*] )
points_cart[2,*] = sin( points_coord_use[1,*] )

; Construct first rotation matrix
; (Rotation around plane of the equator, from the Greenwich meridian to the 
; meridian containing the North Pole in the new coordinates)
rotate_matrix = fltarr( 3, 3 )
if double_opt eq 1 then rotate_matrix = double( rotate_matrix )
;temp = pi - pole_coord_use[0]
temp = pole_coord_use[0]
rotate_matrix[0,0] = cos( temp )
rotate_matrix[0,1] = sin( temp )
rotate_matrix[1,0] = -sin( temp )
rotate_matrix[1,1] = cos( temp )
rotate_matrix[2,2] = 1.
; Perform this rotation around the original pole
points_cart = rotate_matrix # points_cart

; Construct second rotation matrix
; (Rotation in the plane of the current meridian from geographic pole to the 
; North Pole in the new coordinates)
rotate_matrix[*,*] = 0.
temp = pole_coord_use[1] - pi / 2.
rotate_matrix[0,0] = cos( temp )
rotate_matrix[0,2] = sin( temp )
rotate_matrix[2,0] = -sin( temp )
rotate_matrix[2,2] = cos( temp )
rotate_matrix[1,1] = 1.
; Perform this rotation along the rotated meridian
points_cart = rotate_matrix # points_cart

; Convert back to longitude-latitude format (under the new coordinate system)
new_coord = 0. * points_coord_use
new_coord[0,*] = atan( points_cart[1,*], points_cart[0,*] ) / degrad
new_coord[1,*] = atan( points_cart[2,*], $
    sqrt( points_cart[0,*] ^ 2. + points_cart[1,*] ^ 2. ) ) / degrad

;***********************************************************************

return, new_coord
END
