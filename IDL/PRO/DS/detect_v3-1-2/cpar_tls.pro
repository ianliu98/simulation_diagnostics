;+
; NAME:
;    CPAR_TLS
;
; COPYRIGHT:
;    Copyright (1999) Myles Allen, Space Science Department, Rutherford 
;    Appleton Laboratory.
;    Prepared under contract to the Hadley Centre for Climate Prediction and 
;    Research.
;    Revised (2011) Daithi Stone under contract to the U.S. Department of 
;    Energy's Office of Science, Office of Biological and Environmental 
;    Research and the U.S. National Oceanic and Atmospheric Administration's 
;    Climate Program Office, via the International Detection and Attribution 
;    Group.
;
; PURPOSE:
;    This function computes 1-D and m-D parametric confidence intervals from 
;    total least squares regression.
;
; CATEGORY:
;    Optimal Detection Package v3.1.2
;
; NOTES:
;    v3.1.0 is no longer compatible with PV-WAVE.  v2.1 is the most recent 
;    version tested on PV-WAVE.
;
; CALLING SEQUENCE:
;    Result = cpar_tls( BETA_TL, EIGLOAD, EIGVAL, T_CRIT, N_POINT )
;
; INPUTS:
;    BETA_TL:  A vector of best fit regression coefficients.  Of length N_COEF.
;    EIGLOAD:  A floating point matrix of size (N_COEF+1)*(N_COEF+1) containing 
;        the N_COEF+1 loadings corresponding to EIGVAL.
;    EIGVAL:  A floating point vector of size (N_COEF+1) containing the 
;        N_COEF+1 eigenvalues of 
;        transpose([[DATA_INDEP],[DATA_DEP]])#[[DATA_INDEP],[DATA_DEP]].  Here 
;        DATA_INDEP is a N_DATA*N_COEF matrix of N_COEF independent variables 
;        and DATA_DEP is a vector containing the dependent variable of length 
;        N_DATA.
;    N_POINT:  The number of points to sample along an isopleth of probability 
;        when estimating the multivariate probability density function of the 
;        scaling factors.  The function revises the value down to 2 if 
;        N_COEF=1, where N_COEF is the number of regression coefficients, or 
;        otherwise to N_PERDIM^(N_COEF-1) where N_PERDIM is the largest integer 
;        such that the revised N_POINT is less than or equal to the original 
;        value.  An isopleth of probability is an N_COEF dimensional ellipsoid 
;        and is sampled by the N_POINT points in a polar coordinate system.
;    T_CRIT:  A vector of critical points on one dimensional to N_THRESH 
;        dimensional t-distributions at which to estimate the confidence 
;        intervals on the estimates of the regression coefficients (betas) and 
;        attributable components (ONEDIM_PROJECT_CONF).  Of length N_THRESH.
;
; KEYWORD PARAMETERS:
;    DIST_AREA:  Returns a vector containing the fraction of the total surface 
;        of an isopleth taken by each of the N_POINT points on the isopleth.  
;        Of length N_POINT.  The total area of the isopleth is 
;        total(DIST_AREA)=1.
;    DOUBLE:  If set then calculations are done in double precision arithmetic. 
;        The default is single precision.
;    ONEDIM_BETA_DIST:  Returns the locations of the quantiles of the 
;        one-dimensional probability distributions of the regression 
;        coefficients (betas).  The t-values of the quantiles to sample are 
;        defined in ONEDIM_BETA_T.  Returns an array of size 
;        [N_ONEDIM_BETA,N_COEF] where N_COEF is the number of regression 
;        coefficients and N_ONEDIM_BETA is the size of ONEDIM_BETA_T.  Input 
;        from ONEDIM_BETA_T and ONEDIM_INDEP is required.
;    ONEDIM_BETA_T:  A vector of size N_ONEDIM_BETA containing the 
;        t-distribution values of the quantiles at which to sample the 
;        one-dimensional probability distributions of the regression 
;        coefficients (betas).
;    ONEDIM_INDEP:  An array of parameter values for which to return (in 
;        ONEDIM_PROJECT_CONF) the one dimensional confidence intervals of the 
;        attributable component.  Of size N_COEF*N_ONEDIM where N_COEF is the 
;        number of parameters and N_ONEDIM is the number of parameter 
;        combinations to consider.
;    ONEDIM_PROJECT_CONF:  An array of size 3*N_ONEDIM containing the best 
;        ([0,*]), lower confidence range ([1,*]), and upper confidence range 
;        ([2,*]) estimates of the one dimensional attributable components given 
;        the hypothetical parameter combinations in ONEDIM_INDEP.  N_ONEDIM is 
;        the number of parameter combinations to examine.
;    ONEDIM_PROJECT_DIST:  If set, then this returns an array describing the 
;        likelihood functions of the one dimensional attributable component 
;        given the parameter combinations in ONEDIM_INDEP.  Of size 
;        N_ONEDIM_BETA*ONEDIM_INDEP.  The attributable values of t-distribution 
;        value ONEDIM_BETA_T[i] are given in ONEDIM_PROJECT_DIST[i,*].  If 
;        ONEDIM_BETA_T equals T_crit, then ONEDIM_PROJECT_DIST equals 
;        ONEDIM_PROJECT_CONF.
;    Z_DATA:  An array containing the values of the independent and dependent 
;        variables.  Of size N_DATA*(N_COEF+1) where N_DATA is the length of 
;        the dependent variable and N_COEF is the number of independent 
;        variables (and so regression coefficients).  Elements [*,0:N_COEF-1] 
;        contain the values of the independent variables, while elements 
;        [*,N_COEF] contain the values of the dependent variable.
;    Z_POS:  If Z_DATA is input, this returns an array containing the 
;        estimated locations (values) of the probability density estimates of 
;        the values of the estimate noise free independent and dependent 
;        variables.  Of size [N_DATA,N_COEF+1,N_POINT,N_THRESH], where N_POINT 
;        is the number of points to use to sample along an isopleth of 
;        probability and N_THRESH is the number of isopleths as definied in 
;        T_CRIT.  Elements [*,0:N_COEF-1,J,K] correspond to the estimated 
;        values of the independent variables X_DATA corresponding to 
;        Result[*,J,K], that is the Jth sampled point on the (K+1) dimensional 
;        confidence interval surface.  Elements [*,N_COEF,J,K] correspond to 
;        the estimated values of the independent variables DATA_INDEP 
;        corresponding to Result[*,J,K].  Also see Result, N_Point, and T_CRIT.
;
; OUTPUTS:
;    Result:  Returns an N_COEF*N_POINT*N_THRESH array containing the 
;        coordinates of the isopleth of probability surfaces of the estimates 
;        of the regression coefficients.  N_COEF is the number of regression 
;        coefficients, N_POINT is the sample size on a surface, and N_THRESH is 
;        the number of isopleth surfaces.  Value [I,J,K] gives the coordinate 
;        along the direction of regression coefficient I of the sampling point 
;        J on the (K+1) dimensional confidence interval surface.  The 
;        confidence intervals to use are defined in T_crit.  If Result[*,J,K]=0 
;        then the confidence interval is unbounded.
;    DIST_AREA, ONEDIM_BETA_DIST, ONEDIM_PROJECT_DIST, Z_POS
;
; USES:
;    sphere_sampling.pro
;
; PROCEDURE:
;
; REFERENCES:
;    Allen, M. R., and P. A. Stott.  2003.  Estimating signal amplitudes in 
;        optimal fingerprinting. Part I: theory.  Climate Dynamics, 21, 477-491.
;    Stott, P. A., M. R. Allen, and G. S. Jones.  2003.  Estimating signal 
;        amplitudes in optimal fingerprinting. Part II: application to general 
;        circulation models.  Climate Dynamics, 21, 493-500.
;
; EXAMPLE:
;    See demo_gendetec.pro for a demonstration of the use of the optimal 
;    detection routines.
;
; MODIFICATION HISTORY:
;    Written by:  Myles Allen (m.r.allen@rl.ac.uk), 1999-05-15 (v1.0)
;    Modified:  MRA, 1999-08-02 (Revise to return n_point as no. of points on 
;        C-intvls;  v1.1)
;    Modified:  MRA, 1999-08-02 (Simplify computation of C-intvls;  v1.2)
;    Modified:  MRA, 1999-08-09 (Include Z_poss keyword;  v1.3)
;    Modified:  Daithi Stone (stoned@atm.ox.ac.uk), 2004-06-28 (Documentation 
;        for inclusion in routine library)
;    Modified:  DAS, 2005-03-13 (Added PV_WAVE keyword;  updated compliance 
;        with svdpvw.pro;  updated documentation)
;    Modified:  DAS, 2005-09-01 (Added B1DIM, BT1DIM, DOUBLE keywords;  updated 
;        documentation;  v3.0.0)
;    Modified:  DAS, 2008-04-08 (Fixed bug in rotation of B1DIM)
;    Modified:  DAS, 2009-07-29 (Fixed treatment of open-endedness in 
;        calculation of B1DIM)
;    Modified:  DAS, 2010-02-11 (Added H1DIM keyword;  edited documentation 
;        formatting;  added PV_WAVE note;  unset test options)
;    Modified:  DAS, 2011-10-09 (Discontinued PV-WAVE compatibility;  removed 
;        PV_WAVE keyword;  outsourced to sphere_sampling.pro;  removed pre-v3.0 
;        header documentation;  changed use of total_1d.pro to IDL's total;  
;        altered code variable names;  introduced Z_DATA to take over former 
;        input role of Z_POS;  removed T1DIM keyword;  changed input and output 
;        variable names (see translation notes below);  v3.1.0)
;    Modified:  DAS, 2012-01-24 (Edited for compliance with GDL;  v3.1.1)
;    Modified:  DAS, 2012-09-26 (Set to output ONEDIM_PROJECT_CONF only if 
;        N_ONEDIM > N_COEF;  v3.1.2)
;-

;***********************************************************************

FUNCTION CPAR_TLS, $
    BETA_TILDE, EIGLOAD, EIGVAL, T_CRIT, N_POINT, $
    DIST_AREA=dist_area, $
    ONEDIM_BETA_DIST=onedim_beta_dist, ONEDIM_BETA_T=onedim_beta_t, $
      ONEDIM_INDEP=onedim_indep, ONEDIM_PROJECT_CONF=onedim_project_conf, $
      ONEDIM_PROJECT_DIST=onedim_project_dist, $
    Z_DATA=z_data, Z_POS=z_pos, $
    DOUBLE=double_opt

; Changes in input/output variable names from v3.0.0 to v3.1.0
;   changed b1dim to onedim_beta_dist
;   changed betatl to beta_tilde
;   changed bt1dim to onedim_beta_t
;   changed c1dim to onedim_project_conf
;   changed cintvl to beta_confsurf
;   changed d1dim to onedim_indep
;   changed evalus to eigval
;   changed evects to eigload
;   changed h1dim to onedim_project_dist
;   changed npoint to n_point
;   changed p_area to dist_area
;   changed z_poss to z_pos

;***********************************************************************
; Options and checks

; Option for double precision
double_opt = keyword_set( double_opt )
one = 1.
if double_opt eq 1 then one = double( one )

; Warning that input role of Z_POS is no done by Z_DATA
if not( keyword_set( z_data ) ) and keyword_set( z_pos ) then begin
  stop, 'cpar_tls.pro:  Z_DATA now performs input role formely done by Z_POS.'
endif

;***********************************************************************
; Constants

; The number of regression coefficients inputted
n_coef = n_elements( beta_tilde )

; Check that the number of points requested for sampling is consistent with 
; dimensionality of surface
n_point = sphere_sampling( n_coef, n_point=n_point, frac_area=dist_area, $
    sample=1, coord=x_coord_unit )

; Find the number of thresholds for t-tests
n_thresh = n_elements( t_crit )

; The number of data points in the dependent variable
if keyword_set( z_data ) then n_data = n_elements( z_data ) / ( n_coef + 1 )

; The number of parameter values for which to return one-dimensional confidence 
; intervals
n_onedim = n_elements( onedim_indep ) / n_coef

;***********************************************************************

; Initialise working and output arrays
weight = one * fltarr( n_coef+1, n_point )
beta_confsurf = one * fltarr( n_coef, n_point, n_thresh )
if n_onedim gt 0 then onedim_project_conf = one * fltarr( 3, n_onedim )

; Initialise array z_pos if z_data is input
if keyword_set( z_data ) then begin
  ; Initialise z_pos (bug fix in last dimension by DAS)
  z_pos = one * fltarr( n_data, n_coef+1, n_point, n_thresh )
endif

; Compute n_thresh n_coef-dimensional confidence intervals
; (simplified to work in eigenvector coordinates).
; Iterate through dimensions (number of dimensions = 1+i_thresh)
for i_thresh = 0, n_thresh - 1 do begin
  ; Re-scale coordinates on the unit sphere to radius t_crit
  x_coord_use = x_coord_unit * t_crit[i_thresh]
  ; Compute weights on first n_coef eigenvectors by dividing x_coord_use by 
  ; square root of delta-eigenvalue
  for i_coef = 0, n_coef - 1 do begin
    weight[i_coef,*] = x_coord_use[i_coef,*] $
        / sqrt( eigval[i_coef] - eigval[n_coef] )
  endfor
  x_coord_use = 0
  ; Compute sum squared weights
  weight_sum_squared = total( weight[0:n_coef-1,*] ^ 2, 1 )
  ; Check if the sum is greater than unity (indicates open-ended interval)
  id_imag = where( weight_sum_squared gt 1., n_id_imag )
  if n_id_imag gt 0 then begin
    weight_sum_squared[id_imag] = 1.
    print, 'cpar_tls.pro:  Open-ended confidenc interval:  ', 1 + i_thresh
  endif
  id_real = where( weight_sum_squared lt 1., n_id_real )
  ; Do not compute anything if no bounds on confidence interval
  if n_id_real eq 0 then begin
    print, 'cpar_tls.pro:  No bounds on confidence interval:  ', 1 + i_thresh
  ; Otherwise compute stuff
  endif else begin
    ; Compute weight on (n_coef+1)th eigenvector from normalisation constraint.
    ; Use positive square root to give maximum projection onto last eigenvector.
    weight[n_coef,*] = sqrt( 1. - weight_sum_squared )
    weight_sum_squared = 0
    ; Transform to normal coordinates
    eigload_normal = eigload # weight
    ; Convert to conventional regression-coefficient-like scaling parameters
    for i_coef = 0, n_coef - 1 do begin
      beta_confsurf[i_coef,id_real,i_thresh] = -eigload_normal[i_coef,id_real] $
          / eigload_normal[n_coef,id_real]
    endfor
    ; Compute one-dimensional confidence limits
    if ( i_thresh eq 0 ) and keyword_set( onedim_indep ) then begin
      for i_onedim = 0, n_onedim - 1 do begin
        ; Project one-dimensional directions onto ellipsoids
        onedim_project_conf[0,i_onedim] $
            = transpose( [ onedim_indep[*,i_onedim] ] ) # beta_tilde
        temp = transpose( [ onedim_indep[*,i_onedim] ] ) $
            # reform( beta_confsurf[*,id_real,i_thresh] )
        onedim_project_conf[1,i_onedim] = min( temp )
        onedim_project_conf[2,i_onedim] = max( temp )
      endfor
      ; Clear memory
      temp_orthog = 0
      temp_orthog_weight = 0
    endif
    ; If we want noise-free estimates of the distributions of the variables
    if keyword_set( z_data ) then begin
      ; Iterate over real values (and allow long integer counter)
      if n_id_real le 2^15-1 then n_id_real = fix( n_id_real )
      for i_real = 0 * n_id_real, n_id_real - 1 do begin
        z_pos[*,*,id_real[i_real],i_thresh] = z_data $
            - ( z_data # eigload_normal[*,id_real[i_real]] ) $
            # transpose( eigload_normal[*,id_real[i_real]] )
      endfor
    endif
  endelse
endfor

; Compute one-dimensional likelihood distributions on the regressed 
; components.  The procedure is identical to above.
if keyword_set( onedim_project_dist ) then begin
  ; Initialise distribution arrays
  n_onedim_beta = n_elements( onedim_beta_t )
  beta_confsurf_project = one * fltarr( n_coef, n_point )
  if n_onedim gt n_coef then begin
    onedim_project_dist = one * fltarr( n_onedim_beta, n_onedim - n_coef )
  endif
  ; Iterate through quantiles
  for i_onedim_beta = 0, n_onedim_beta - 1 do begin
    ; Re-scale coordinates on the unit sphere to radius 
    ; onedim_beta_t[i_onedim_beta]
    x_coord_use = x_coord_unit * onedim_beta_t[i_onedim_beta]
    ; Compute weights on first n_coef eigenvectors by dividing x_coord_use by 
    ; square root of delta-eigenvalue
    for i_coef = 0, n_coef - 1 do begin
      weight[i_coef,*] = x_coord_use[i_coef,*] / sqrt( eigval[i_coef] $
          - eigval[n_coef] )
    endfor
    x_coord_use = 0
    ; Compute sum squared weights
    weight_sum_squared = total( weight[0:n_coef-1,*] ^ 2, 1 )
    ; Check if the sum is greater than unity (indicates open-ended interval)
    id_imag = where( weight_sum_squared gt 1., n_id_imag)
    if n_id_imag gt 0 then begin
      weight_sum_squared[id_imag] = 1.
      print, 'cpar_tls.pro:  ' $
          + 'Open-ended confidence interval for ONEDIM_PROJECT_DIST.'
    endif
    id_real = where( weight_sum_squared lt 1., n_id_real )
    ; Do not compute anything if no bounds on C-intvl
    if n_id_real eq 0 then begin
      print, 'cpar_tls.pro:  ' $
         + 'No bounds on confidence interval for ONEDIM_PROJECT_DIST.'
    endif else begin
      ; Compute weight on (n_coef+1)th eigenvector from normalisation 
      ; constraint.
      ; Use positive square root to give maximum projection onto last 
      ; eigenvector.
      weight[n_coef,*] = sqrt( 1. - weight_sum_squared )
      weight_sum_squared = 0
      ; Transform to normal coordinates
      eigload_normal = eigload # weight
      ; Convert to conventional regression-coefficient-like scaling parameters
      for i_coef = 0, n_coef - 1 do begin
        beta_confsurf_project[i_coef,id_real] $
            = -eigload_normal[i_coef,id_real] / eigload_normal[n_coef,id_real]
      endfor
      eigload_normal = 0
      ; Iterate through scenario states
      for i_onedim = 0, n_onedim - 1 - n_coef do begin
        ; Project one-dimensional directions onto ellipsoids
        temp = transpose( [ onedim_indep[*,i_onedim] ] ) $
            # reform( beta_confsurf_project[*,id_real] )
        if onedim_beta_t[i_onedim_beta] le 0. then begin
          onedim_project_dist[i_onedim,i_onedim] = max( temp )
        endif else begin
          onedim_project_dist[i_onedim,i_onedim] = min( temp )
        endelse
      endfor
    endelse
  endfor
  ; Clear memory
  beta_confsurf_project = 0
  d1d = 0
endif

; Determine one-dimensional likelihood functions for the regression coefficients
if keyword_set( onedim_beta_t ) then begin
  ; Initialise array of the locations of quantiles in the likelihood functions
  onedim_beta_dist $
      = reform( one * fltarr( n_elements( onedim_beta_t ), n_coef ), $
      n_elements( onedim_beta_t ), n_coef )
  ; Initialise a temporary work variable
  temp = one * fltarr( n_coef, n_point )
  ; Iterate through quantiles of the distribution
  for i_onedim_beta = 0, n_elements( onedim_beta_t ) - 1 do begin
    ; Re-scale coordinates on the unit sphere to radius 
    ; onedim_beta_t[i_onedim_beta]
    x_coord_use = x_coord_unit * onedim_beta_t[i_onedim_beta]
    ; Compute weights on first n_coef eigenvectors by dividing x_coord_use by 
    ; square root of delta-eigenvalue
    for i_coef = 0, n_coef - 1 do begin
      weight[i_coef,*] = x_coord_use[i_coef,*] $
          / sqrt( eigval[i_coef] - eigval[n_coef] )
    endfor
    x_coord_use = 0
    ; Compute sum squared weights
    weight_sum_squared = total( weight[0:n_coef-1,*] ^ 2, 1 )
    ; Check if the sum is greater than unity (indicates open-ended interval)
    id_imag = where( weight_sum_squared gt 1., n_id_imag )
    if n_id_imag gt 0  then begin
      weight_sum_squared[id_imag] = 1.
      print, 'cpar_tls.pro:  Open-ended confidence interval on ' $
          + 'one-dimensional likelihood functions.'
    endif
    id_real = where( weight_sum_squared lt 1., n_id_real )
    ; Do not compute anything if no bounds on beta_confsurf
    if n_id_real gt 0 then begin
      ; Compute weight on (n_coef+1)th eigenvector from normalisation 
      ; constraint.
      ; Use positive square root to give maximum projection onto last 
      ; eigenvector.
      weight[n_coef,*] = sqrt( 1. - weight_sum_squared )
      weight_sum_squared = 0
      ; Transform to normal coordinates
      eigload_normal = eigload # weight
      ; Convert to conventional regression-coefficient-like scaling parameters
      for i_coef = 0, n_coef - 1 do begin
        temp[i_coef,id_real] = -eigload_normal[i_coef,id_real] $
            / eigload_normal[n_coef,id_real]
      endfor
      eigload_normal = 0
      temp = temp ## transpose( [ onedim_indep[*,n_onedim-n_coef:n_onedim-1] ] )
      ; Compute location of current quantile
      if onedim_beta_t[i_onedim_beta] lt 0. then begin
        for i_coef = 0, n_coef - 1 do begin
          onedim_beta_dist[i_onedim_beta,i_coef] = max( temp[i_coef,id_real] )
        endfor
      endif else begin
        for i_coef = 0, n_coef - 1 do begin
          onedim_beta_dist[i_onedim_beta,i_coef] = min( temp[i_coef,id_real] )
        endfor
      endelse
    endif
  endfor
endif
; Clear memory
temp = 0
weight = 0

;***********************************************************************
; The end

return, beta_confsurf
END
