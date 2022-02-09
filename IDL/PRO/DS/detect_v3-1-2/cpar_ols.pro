;+
; NAME:
;    CPAR_OLS
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
;    ordinary least squares regression.
;
; CATEGORY:
;    Optimal Detection Package v3.1.2
;
; NOTES:
;    v3.1.0 is no longer compatible with PV-WAVE.  v2.1 is the most recent 
;    version tested on PV-WAVE.
;
; CALLING SEQUENCE:
;    Result = cpar_ols( BETA_TILDE, COV_NOISE, T_CRIT, N_POINT )
;
; INPUTS:
;    BETA_TILDE:  A vector of best fit regression coefficients.  Of length 
;        N_COEF.
;    COV_NOISE:  An array containing the estimated covariance.  Of size 
;        N_COEF*N_COEF.
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
;    ONEDIM_BETA_DIST:  Returns the locations of the quantiles of the one 
;        dimensional probability distributions of the regression coefficients 
;        (betas).  The t-values of the quantiles to sample are defined in 
;        ONEDIM_BETA_T.  Returns an array of size [N_ONEDIM_BETA,N_COEF] where 
;        N_COEF is the number of regression coefficients and N_ONEDIM_BETA is 
;        the size of ONEDIM_BETA_T.  Input from ONEDIM_BETA_T and is required.  
;        This uses the last N_COEF of the N_ONEDIM parameter combinations in 
;        ONEDIM_INDEP to define the directions of the betas to examine, 
;        otherwise it assumes unrotated directions.
;    ONEDIM_BETA_T:  A vector of size N_ONEDIM_BETA containing the 
;        t-distribution values of the quantiles at which to sample the 
;        one-dimensional probability distributions of the regression 
;        coefficients (betas).
;    ONEDIM_INDEP:  An array of parameter values for which to return (in 
;        ONEDIM_PROJECT_CONF and optionally ONEDIM_PROJECT_DIST) the one 
;        dimensional confidence intervals of the attributable component.  Of 
;        size N_COEF*N_ONEDIM where N_COEF is the number of parameters and 
;        N_ONEDIM is the number of parameter combinations to consider.
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
;        ONEDIM_BETA_T equals T_CRIT, then ONEDIM_PROJECT_DIST equals 
;        ONEDIM_PROJECT_CONF.
;    Z_DATA:  An array containing the values of the independent and dependent 
;        variables.  Of size N_DATA*(N_COEF+1) where N_DATA is the length of 
;        the dependent variable and N_COEF is the number of independent 
;        variables (and so regression coefficients).  Elements [*,0:N_COEF-1] 
;        contain the values of the independent variables, while elements 
;        [*,N_COEF] contain the values of the dependent variable.
;    Z_POS:  If Z_DATA is input, this returns an array containing the 
;        estimated locations (values) of the probability density estimates of 
;        the values of the estimated noise free independent and dependent 
;        variables.  Of size [N_DATA,N_COEF+1,N_POINT,N_THRESH], where N_POINT 
;        is the number of points to use to sample along an isopleth of 
;        probability and N_THRESH is the number of isopleths as definied in 
;        T_CRIT.  Elements [*,0:N_COEF-1,J,K] correspond to the estimated 
;        values of the independent variables X_DATA corresponding to 
;        RESULT[*,J,K], that is the Jth sampled point on the (K+1) dimensional 
;        confidence interval surface.  Elements [*,N_COEF,J,K] correspond to 
;        the estimated values of the independent variable Y_DATA corresponding 
;        to RESULT[*,J,K].  Also see Result, N_POINT, and T_CRIT.
;
; OUTPUTS:
;    Result:  Returns an N_COEF*N_POINT*N_THRESH array containing the 
;        coordinates of the isopleth of probability surfaces of the estimates 
;        of the regression coefficients.  N_COEF is the number of regression 
;        coefficients, N_POINT is the sample size on a surface, and N_THRESH is 
;        the number of isopleth surfaces.  Value [I,J,K] gives the coordinate 
;        along the direction of regression coefficient I of the sampling point 
;        J on the (K+1) dimensional confidence interval surface.  The 
;        confidence intervals to use are defined in T_CRIT.
;    DIST_AREA, ONEDIM_BETA_DIST, ONEDIM_PROJECT_DIST, Z_POS
;
; USES:
;    sphere_sampling.pro
;
; PROCEDURE:
;    This function calculates the singular value decomposition (SVD) of the 
;    estimated covariance matrix COV_NOISE.
;
; REFERENCES:
;    Allen, M. R., and S. F. B. Tett.  1999.  Checking internal consistency in 
;        optimal fingerprinting.  Climate Dynamics, 15, 419-434.
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
;    Modified:  MRA, 1999-08-02 (Revise to return npoint as no. of points on 
;        C-intvls; v1.1)
;    Modified:  MRA, 1999-08-09 (Include Z_POSS keyword; v1.3)
;    Modified:  MRA, 2000-08-03 (Allow threshold-type prior constraints; v2.0)
;    Modified:  Daithi Stone (stoned@atm.ox.ac.uk), 2004-06-28 (Documentation 
;        for inclusion in routine library)
;    Modified:  DAS, 2005-03-13 (Added PV_WAVE keyword;  updated compliance 
;        with svdpvw.pro;  updated documentation)
;    Modified:  DAS, 2005-04-06 (Allowed Npoint to be a long integer)
;    Modified:  DAS, 2005-09-01 (Added B1DIM, BT1DIM, DOUBLE keywords;  fixed 
;        bug in Z_POSS initialisation;  updated documentation;  v3.0.0)
;    Modified:  DAS, 2008-04-08 (Fixed bug in rotation of B1DIM)
;    Modified:  DAS, 2010-02-11 (Added H1DIM keyword;  edited documentation 
;        formatting;  added PV_WAVE note)
;    Modified:  DAS, 2011-11-15 (Discontinued PV-WAVE compatibility;  removed 
;        pre-v3.0 header documentation;  outsourced to sphere_sampling.pro;  
;        switched used of svdpvw.pro to pca.pro;  changed BETATL, Covb2s, and 
;        NPOINT inputs to BETA_TILDE, COV_NOISE, and N_POINT;  altered code 
;        variable names (see translation notes below);  introduced Z_DATA to 
;        take over former input role of Z_POS;  v3.1.0)
;    Modified:  DAS, 2012-01-24 (Edited for compliance with GDL;  v3.1.1)
;    Modified:  DAS, 2012-02-14 (Edited to allow case where N_ONEDIM equals 1;  
;        v3.1.2)
;-

;***********************************************************************

FUNCTION CPAR_OLS, $
    BETA_TILDE, COV_NOISE, T_CRIT, N_POINT,$
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
;   changed covb2s to cov_noise
;   changed d1dim to onedim_indep
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
  stop, 'cpar_ols.pro:  Z_DATA now performs input role formerly done by Z_POS.'
endif

;***********************************************************************
; Constants

; The number of regression coefficients inputted
n_coef = n_elements( beta_tilde )

; Check that the number of points requested for sampling is consistent with 
; dimensionality of surface
n_point_old = n_point
dist_area = 1.
n_point = sphere_sampling( n_coef, n_point=n_point_old, frac_area=dist_area, $
    sample=1, coord=x_coord )

; Find the number of thresholds for t-tests
n_thresh = n_elements( t_crit )

; The number of data points in the dependent variable
if keyword_set( z_data ) then n_data = n_elements( z_data ) / ( n_coef + 1 )

; The number of parameter values for which to return one-dimensional confidence 
; intervals
n_onedim = n_elements( onedim_indep ) / n_coef

;***********************************************************************

; Initialise output array
beta_confsurf = one * fltarr( n_coef, n_point, n_thresh )

; Initialise array Z_POS if Z_DATA is input
if keyword_set( z_data ) then begin
  ; Initialise Z_POS (bug fix in last dimension by DAS)
  z_pos = one * fltarr( n_data, n_coef+1, n_point, n_thresh )
endif

; compute 1-dimensional confidence intervals
if keyword_set( onedim_indep ) then begin
  ; Initialise array of one dimensional confidence intervals on output variable
  onedim_project_conf = reform( one * fltarr( 3, n_onedim ), 3, n_onedim )
  ; Iterate through input parameter sets
  for i_onedim = 0, n_onedim - 1 do begin
    ; Compute best-estimate value
    onedim_project_conf[0,i_onedim] $
        = transpose( [ onedim_indep[*,i_onedim] ] ) # beta_tilde
    ; Compute uncertainty range
    half_range = t_crit[0] * sqrt( transpose( [ onedim_indep[*,i_onedim] ] ) $
        # cov_noise # onedim_indep[*,i_onedim] )
    ; Compute min and max range values
    onedim_project_conf[1:2,i_onedim] = onedim_project_conf[0,i_onedim] $
        + [ -1, 1 ] * half_range[0]
  endfor
  ; Compute one-dimensional likelihood distributions on the regressed components
  if keyword_set( onedim_project_dist ) then begin
    ; Initialise output array of the distributions
    n_onedim_beta = n_elements( onedim_beta_t )
    onedim_project_dist = one * fltarr( n_onedim_beta, n_onedim - n_coef )
    ; Iterate through variable combinations
    for i = 0, n_onedim - 1 - n_coef do begin
      ; Compute the best-estimate value
      temp = ( transpose( [ onedim_indep[*,i] ] ) # beta_tilde )[0]
      ; Estimate the distribution for this variable combination
      temp1 = onedim_beta_t * ( sqrt( transpose( [ onedim_indep[*,i] ] ) $
          # cov_noise # onedim_indep[*,i] ) )[0]
      onedim_project_dist[*,i] = temp - temp1
    endfor
  endif
endif

; Compute one-dimensional distributions
if keyword_set( onedim_beta_t ) then begin
  ; Initialise output of distribution values
  onedim_beta_dist = one * fltarr( n_elements( onedim_beta_t ), n_coef )
  ; Iterate through independent variables
  for i_coef = 0, n_coef - 1 do begin
    ; Compute locations of the quantiles of the distribution
    temp_1 = transpose( [ onedim_indep[*,n_onedim-n_coef+i_coef] ] ) $
        # beta_tilde
    temp_2 = sqrt( transpose( [ onedim_indep[*,n_onedim-n_coef+i_coef] ] ) $
        # cov_noise # onedim_indep[*,n_onedim-n_coef+i_coef] )
    onedim_beta_dist[*,i_coef] = temp_1[0] - onedim_beta_t * temp_2[0]
  endfor
endif

; Compute n_thresh n_indep-dimensional confidence intervals.
; Perform principal component analysis on the noise covariance matrix
svdc, reform( cov_noise, n_coef, n_coef ), temp_eigval, temp_eigvec, temp, $
    double=double_opt
; Transfer variance to eigenvectors
for i_coef = 0, n_coef - 1 do begin
  temp_eigvec[*,i_coef] = temp_eigvec[*,i_coef] * sqrt( temp_eigval[i_coef] )
endfor
temp_eigval = 0
; Iterate through probability surfaces
for i_thresh = 1, n_thresh do begin
  ; Calculate the coordinates of the this probability surface
  beta_confsurf[*,*,i_thresh-1] = t_crit[i_thresh-1] $
      * ( temp_eigvec # x_coord ) + rebin( beta_tilde, n_coef, n_point )
  ; If we want noise-free estimates of the distributions of the variables
  if keyword_set( z_pos ) then begin
    ; Iterate through sampling points of the distribution
    for i_point = 0l, n_point - 1l do begin
      ; Calculate noise-free estimates of the independent variables
      z_pos[*,0:n_coef-1,i_point,i_thresh-1] = z_data[*,0:n_coef-1]
      ; Calculate noise-free estimates of the dependent variable
      z_pos[*,n_coef,i_point,i_thresh-1] = z_data[*,0:n_coef-1] $
          # beta_confsurf[*,i_point,i_thresh-1]
    endfor
  endif
endfor
; Clear memory
temp_eigvec = 0

;***********************************************************************
; The end

return, beta_confsurf
END
