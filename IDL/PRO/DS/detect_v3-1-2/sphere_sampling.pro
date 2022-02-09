;+
; NAME:
;    SPHERE_SAMPLING
;
; COPYRIGHT:
;    Copyright (2011) Daithi Stone under contract to the U.S. Department of 
;    Energy's Office of Science, Office of Biological and Environmental 
;    Research and the U.S. National Oceanic and Atmospheric Administration's 
;    Climate Program Office, via the International Detection and Attribution 
;    Group.
;
; PURPOSE:
;    This function determines the number of points to sample along the surface 
;    of a multi-dimensional sphere based on a suggestion.
;
; CATEGORY:
;    Optimal Detection Package v3.1.2
;
; CALLING SEQUENCE:
;    Result = sphere_sampling( N_DIM )
;
; INPUTS:
;    N_DIM:  A required scalar integer defining the dimension of the sphere.
;
; KEYWORD PARAMETERS:
;    COORD:  If SAMPLE is set, this returns a N_DIM*N_POINT floating point 
;        array containing the sampling of polar coordinates for the sphere.
;    DOUBLE:  If set then calculations are done in double precision arithmetic. 
;        The default is single precision.
;    FRAC_AREA:  If set then returns a vector containing the fraction of the 
;        total surface taken by each of the N_POINT points.  Of length 
;        N_POINT.  The total area is total(FRAC_AREA)=1.
;    N_PERDIM:  An optional input suggesting the number of points to sample per 
;        dimension along the surface of a multi-dimensional sphere.  Of type 
;        integer or long integer.  If N_POINT is set then it overrides N_PERDIM.
;        Returns a revised value based on the revised value of N_POINT given in 
;        Result.
;    N_POINT:  An optional input suggesting the number of points to sample 
;        along the surface of a multi-dimensional sphere.  Of type integer or 
;        long integer.
;    SAMPLE:  If set then sampling of the polar coordinates on the sphere is 
;        performed and returned in COORD.
;
; OUTPUTS:
;    Result:  The default or revised number of sampling points on the sphere's 
;        surface.
;    COORD, FRAC_AREA, N_PERDIM
;
; USES:
;    -
;
; PROCEDURE:
;    The function revises the value of N_POINT or N_PER_DIM^N_DIM down to 2 if 
;    N_DIM=1, or otherwise to N_PERDIM^(N_DIM-1) where N_POINT_PER_DIM is the 
;    largest integer such that the revised N_POINT is less than or equal to the 
;    original value.
;    When sampling, it first samples along longitude and then along the 
;    latitudes of additional dimensions.
;
; EXAMPLE:
;    See cpar_ols.pro and cpar_tls.pro.
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone (stoned@csag.uct.ac.za), 2011-10-14 (Based on 
;        cpar_ols.pro and cpar_tls.pro by Myles Allen;  fixed bug sampling 
;        south pole singularity;  v3.1.0)
;-

;***********************************************************************

FUNCTION SPHERE_SAMPLING, $
    N_DIM, $
    N_POINT=n_point, N_PERDIM=n_perdim, $
    SAMPLE=sample_opt, $
    DOUBLE=double_opt, $
    FRAC_AREA=frac_area, COORD=X_COORD

;***********************************************************************
; Interpret inputs

; Convert N_PERDIM suggestion to N_POINT suggestion
if keyword_set( n_perdim ) and not( keyword_set( n_point ) ) then begin
  n_point = n_perdim ^ ( n_dim - 1 )
endif

; Option for double precision
double_opt = keyword_set( double_opt )
one = 1.
if double_opt eq 1 then begin
  one = double( one )
  pi = !dpi
endif else begin
  pi = !pi
endelse

;***********************************************************************
; Set or revise the sampling density

; If we have a one-dimensional sphere
if n_dim eq 1 then begin
  if n_point ne 2 then begin
    print, 'sphere_sampling.pro:  Only 2 points required for sampling a line.'
    result = 2
  endif else begin
    result = n_point
  endelse
  n_perdim = n_point
; If we have a two-dimensional sphere
endif else if n_dim eq 2 then begin
  result = n_point
  n_perdim = result
; If we have a multi-dimensional sphere
endif else if n_dim gt 2 then begin
  n_perdim = fix( n_point ^ (1. / ( n_dim - 1 ) ) )
  result = ( one * float( n_perdim ) ) ^ ( n_dim - 1 )
  if result ne n_point then begin
    print, 'sphere_sampling.pro:  Number of points revised to ', result
  endif
endif

;***********************************************************************
; Generate array of sampling coordinates on the sphere

; If this is requested
if keyword_set( sample_opt ) then begin

  ; Initialise array of sampling points
  x_coord = one * fltarr( n_dim, result )
  ; Initialise fractional area array
  if keyword_set( frac_area ) then frac_area = replicate( one, result )

  ; If we only have one dimension
  if n_dim eq 1 then begin

    ; The we only need two points
    x_coord[0,*] = [ -one, one ]

  ; If we have multiple dimensions
  endif else begin

    ; Define vectors for generating latitudes
    z_lat_real = sin( pi * ( findgen( n_perdim ) + 0.5 ) / n_perdim ) 
    z_lat_imag = -cos( pi * ( findgen( n_perdim ) + 0.5 ) / n_perdim )

    ; Sample longitude in two dimensions (i.e. along an circle)
    ; (Note longitude is currently sampled at half the equatorial resolution of 
    ; the latitudes.)
    x_coord[0,0:n_perdim-1] = cos( 2. * pi * findgen( n_perdim ) / n_perdim )
    x_coord[1,0:n_perdim-1] = sin( 2. * pi * findgen( n_perdim ) / n_perdim )

    ; Iterate through additional dimensions
    for i_dim = 2, n_dim - 1 do begin
      ; Determine how many sampling points have been determined so far
      n_sofar = one * float( n_perdim ) ^ ( i_dim - one )
      id_sofar = indgen( n_sofar )
      ; Copy these previously sampled points
      x_coord_sofar = x_coord[0:i_dim-1,id_sofar]
      ; Iterate through latitudes in this dimension
      for i_perdim = 0, n_perdim - 1 do begin
        ; Determine the points on this latitude
        id_arc = i_perdim * n_sofar + id_sofar
        ; Set the coordinate of this latitude
        x_coord[i_dim,id_arc] = z_lat_imag[i_perdim]
        ; Weight other dimensions to maintain unit sphere
        x_coord[0:i_dim-1,id_arc] = z_lat_real[i_perdim] $
            * x_coord_sofar
        ; Calculate area corresponding to points on this latitude
        if keyword_set( frac_area ) then begin
          frac_area[id_arc] = abs( z_lat_real[i_perdim] ) * frac_area[id_sofar]
        endif
      endfor
    endfor

  endelse

  ; Make the elements of frac_area sum to unity
  if keyword_set( frac_area ) then frac_area = frac_area / total( frac_area )

endif

;***********************************************************************
; The end

return, result
END
