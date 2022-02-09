;+
; NAME:
;    LINMOD
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
;    This function performs multiple linear regression given an optional 
;    prewhitening operator, externally specified noise realisation, and 
;    optional noise on the independent variables.
;
; CATEGORY:
;    Optimal Detection Package v3.1.2
;
; NOTES:
;    v3.1.0 and later are no longer compatible with PV-WAVE.  v2.1 is the most 
;    recent version tested on PV-WAVE.  v3.1.1 and later are compatible with 
;    GDL.
;
; CALLING SEQUENCE:
;    Result = linmod( DATA_INDEP, DATA_DEP )
;
; INPUTS:
;    DATA_INDEP:  A floating point matrix of size N_DATA*N_INDEP array 
;        containing N_DATA values of N_INDEP independent variables.
;    DATA_DEP:  A floating point vector of size N_DATA containing N_DATA values 
;        of the dependent variable.
;
; KEYWORD PARAMETERS:
;    ADDVAR_DEP:  An optional floating point input array of size N_DATA*N_DATA 
;        containing the covariance matrix of the additional variance in the 
;        dependent variable beyond that which is expected due to limited 
;        sampling.  Thus element [i,j] contains the covariance between data 
;        points i and j.  For example, if the dependent variable comes from 
;        observations, then element [i,j] will represent the covariance between 
;        space-time points i and j coming from errors in observational 
;        measurements.  Note that the covariance from limited sampling should 
;        not be included in ADDVAR_DEP.  If neither ADDVAR_INDEP nor ADDVAR_DEP 
;        are input then EIV is redundant with TLS.
;    ADDVAR_INDEP:  An optional floating point input array of size 
;        N_DATA*N_DATA*N_INDEP containing the inter-estimate covariance 
;        matrices of the additional variance in the independent variable 
;        patterns beyond that which is expected due to limited sampling.  Thus 
;        element [i,j,k] contains the covariance between data points i and j 
;        for independent variable k.  For example, if we are comparing model 
;        simulations against observations, then element [i,j,k] will represent 
;        the covariance across models between space-time points i and j in the 
;        estimation of pattern k.  Note that the covariance from limited 
;        sampling should not be included in ADDVAR_INDEP.  If neither 
;        ADDVAR_INDEP nor ADDVAR_DEP are input then EIV is redundant with TLS.
;    BETA_DIST:  If set then this returns an array containing the location 
;        (beta values) of the probability density estimates of the beta scaling 
;        parameters.  Of size [N_INDEP,N_POINT,N_DIST] where N_INDEP is the 
;        number of scenarios, N_POINT is the number of points sampled along an 
;        isopleth of probability, and N_DIST is the number of isopleths sampled 
;        as definied in T_DIST.  So elements [*,j,k] are the coordinates of the 
;        jth point on the kth isopleth of probability.  Initialise to 1 to 
;        ensure output.  See N_POINT, DIST_AREA, and T_DIST for more 
;        information.
;    BETA_NOISE:  An N_INDEP*N_NOISE matrix containing estimates of the N_INDEP 
;        coefficients estimated from each of the N_NOISE different noise 
;        realisations in DATA_NOISE added to the dependent variable DATA_DEP.
;    BETA_CONFSURF:  Returns an N_INDEP*N_POINT*N_INDEP array containing the 
;        coordinates of the isopleth of probability surfaces of the estimates 
;        of the regression coefficients.  N_INDEP is the number of regression 
;        coefficients and N_POINT is the sample size on a surface.  Value 
;        [I,J,K] gives the coordinate along the direction of regression 
;        coefficient I of the sampling point J on the (K+1) dimensional 
;        confidence interval surface.  The confidence intervals to use are 
;        defined in T_CRIT.  If Result[*,J,K]=0 then the confidence interval is 
;        unbounded (possible with TLS).
;    BLEACH:  A matrix containing a pre-whitening operator for the input 
;        invariables.  Of size N_DOF*N_DATA, where N_DOF is the number of 
;        degrees of freedom (truncation).
;    COV_NOISE:  Returns an N_INDEP*N_INDEP matrix containing the estimated 
;        covariance of the noise realisations in DATA_NOISE in the directions 
;        of the N_INDEP scenario patterns found in the N_INDEP independent 
;        variables of DATA_INDEP.
;    DATA_NOISE:  A floating point matrix of size N_DATA*N_NOISE containing 
;        N_NOISE independent noise realisations of the dependent variable 
;        DATA_DEP for use in estimating confidence intervals.
;    DIST_AREA:  Returns a vector containing the fraction of the total surface 
;        of an isopleth taken by each of the N_POINT points on the isopleth.  
;        Of length N_POINT.  The total area of the isopleth is 
;        total(DIST_AREA)=1.  See BETA_DIST and N_POINT for more information.
;    DOUBLE:  If set then calculations are done in double precision arithmetic. 
;        The default is single precision.
;    FINGER_TRANS:  If set then this returns the transpose of the matrix of 
;        fingerprint patterns.  Note that Result=FINGER_TRANS#DATA_DEP, to 
;        within the EOF space covered to the truncation N_DOF.  Of size 
;        N_INDEP*N_DATA.
;    FRAC_NOISE_VAR:  A floating point array of length N_INDEP.  It contains 
;        the fraction of the variance of the noise in DATA_DEP contained in 
;        each of the N_INDEP variables in DATA_INDEP.  Required for the TLS and 
;        EIV methods.
;    INV_BLEACH:  Returns the pseudo-inverse of BLEACH.
;    N_DOF:  Returns the number of degrees of freedom of the pre-whitened 
;        data.  This corresponds to the truncation in the EOF space used.  If 
;        BLEACH is input then its first dimension is of size N_DOF.
;    N_POINT:  The number of points N_POINT to sample along an isopleth of 
;        probability when estimating the multivariate probability density 
;        function of the scaling factors.  See BETA_DIST, DIST_AREA, and T_DIST 
;        for more information. The default (forced) is 2 if N_INDEP=1 (the 
;        problem is one-dimensional), otherwise the default is 
;        N_PERDIM^(N_INDEP-1) when N_INDEP>1, where N_INDEP is the number of 
;        independent variables and N_PERDIM=fix(N_POINT^(1./N_INDEP-1)).  An 
;        isopleth of probability is an N_INDEP-dimensional ellipsoid and is 
;        sampled by the N_POINT points in a polar coordinate system.
;    OLS:  If set the function uses the ordinary least squares (OLS) algorithm. 
;        The default is to use the total least squares (TLS) algorithm.  This 
;        option overrides the use of FRAC_NOISE_VAR.  This keyword is 
;        superceded by the TYPE keyword and overruled by it.
;    ONEDIM_BETA_DIST:  Returns the locations of the quantiles of the one 
;        dimensional probability distributions of the betas (the regression 
;        coefficients) in the direction of the last N_INDEP parameter 
;        combinations defined in ONEDIM_INDEP.  The t-values of the quantiles 
;        to sample are defined in ONEDIM_BETA_T.  Returns an array of size 
;        [N_ONEDIM_BETA,N_INDEP] where N_INDEP is the number of scenarios and 
;        N_ONEDIM_BETA is the size of ONEDIM_BETA_T.  Input from ONEDIM_BETA_T 
;        and ONEDIM_INDEP is required.
;    ONEDIM_BETA_T:  A vector of size N_ONEDIM_BETA containing the 
;        t-distribution values of the quantiles at which to sample the 
;        one-dimensional probability distributions of the betas (the regression 
;        coefficients).
;    ONEDIM_INDEP:  An array of parameter values (as in DATA_INDEP) for which 
;        to return (in ONEDIM_PROJECT_CONF and optionally ONEDIM_PROJECT_DIST) 
;        the one dimensional confidence intervals of the attributable component 
;        of DATA_DEP scaled according to the beta scaling parameters.  Of size 
;        N_INDEP*N_ONEDIM where N_INDEP is the number of parameters and 
;        N_ONEDIM is the number of parameter combinations to consider.  Values 
;        can be within or outside of the range of values in DATA_INDEP.
;    ONEDIM_PROJECT_CONF:  Returns an array of size 3*N_ONEDIM containing the 
;        best ([0,*]), lower confidence range ([1,*]), and upper confidence 
;        range ([2,*]) estimates of the one dimensional attributable components 
;        of DATA_DEP given the hypothetical parameter combinations in 
;        ONEDIM_INDEP.  N_ONEDIM is the number of parameter combinations to 
;        examine.
;    RESID_SUMSQ:  Returns the sum of squares of the residual differences 
;        between the prewhitened dependent variable DATA_DEP and the best 
;        regression estimate.  If DATA_NOISE is given, then the normalised 
;        residuals are used, with the normalisation done for each element of 
;        DATA_DEP according to the standard deviation of the corresponding 
;        N_NOISE realisations of that element in DATA_NOISE.
;    SCALE_COV_NOISE:  If set, then the covariance estimates are scaled by 
;        dividing by the residual sum of squares (RESID_SUMSQ).
;    T_CRIT:  A vector of critical points on one dimensional to N_INDEP-
;        dimensional t-distributions at which to estimate the confidence 
;        intervals on the estimates of the betas (regression coefficients).  Of 
;        length N_INDEP, the number of independent scenarios.
;    T_DIST:  A vector of length N_DIST containing the T-statistic values on 
;        the isopleths of probability to sample when estimating the probability 
;        density of the regression coefficients.  See BETA_DIST, N_POINT, and 
;        DIST_AREA for more information.
;    TILDE_DEP:  Returns a vector of length N_DATA containing the best estimate 
;        values of noise free DATA_DEP.  While the estimations are done in 
;        truncated EOF space, the results are projected back into normal space 
;        in TILDE_DEP.  But because we only use the leading N_DOF EOFs, the 
;        estimates will look progressively worse the smaller N_DOF is compared 
;        to N_DATA.
;    TILDE_INDEP:  Returns an array containing the best estimate values of 
;        noise free DATA_INDEP.  The array is of size N_DATA*N_INDEP.  While 
;        the estimations are done in truncated EOF space, the results are 
;        projected back into normal space in TILDE_INDEP.  But because we only 
;        use the leading N_DOF EOFs, the estimates will look progressively 
;        worse the smaller N_DOF is compared to N_DATA.
;    TYPE:  An optional string describing the type of regression algorithm to 
;        use.  Possibilities are:
;          'OLS' for ordinary least squares (OLS),
;          'TLS' for total least squares (TLS),
;          'EIV' for error in variables (EIV).
;        If not input then the function makes a guess based on the various 
;        inputs provided.
;    WEIGHT:  An optional floating point vector of length N_DATA defining 
;        weights on the N_DATA elements in the input data arrays.
;    Z_DIST:  Returns an array containing the estimated values of the estimated 
;        noise-free DATA_INDEP and DATA_DEP on the isopleths of probability 
;        defined in T_DIST.  Of size N_DATA*(N_INDEP+1)*N_POINT*N_DIST.  
;        Elements [*,0:N_INDEP-1,j,k] correspond to the estimates of 
;        BETA_DIST[*,j,k] for the independent variables in DATA_INDEP, while 
;        elements [*,N_INDEP,j,k] are for the dependent variable in DATA_DEP.
;    Z_POS:  Returns an array containing the estimated locations (values) of 
;        the probability density estimates of the values of the estimated 
;        noise-free DATA_INDEP and DATA_DEP.  Of size 
;        N_DATA*(N_INDEP+1),N_POINT,N_INDEP] where N_POINT is the number of 
;        points to use to sample probability distribution.  Elements 
;        [*,0:N_INDEP-1,J,K] correspond to the estimated values of DATA_INDEP 
;        corresponding to BETA_CONFSURF[*,J,K], that is the Jth sampled point 
;        on the (K+1) dimensional confidence interval surface.  Elements 
;        [*,N_INDEP,J,K] correspond to the estimated values of DATA_DEP 
;        corresponding to BETA_CONFSURF[*,J,K].  Also see BETA_CONFSURF and 
;        N_POINT.
;
; OUTPUTS:
;    Result:  A vector of length N_INDEP containing the best estimate of the 
;        regression coefficients (betas), where N_INDEP is the number of 
;        independent variables.
;    BETA_CONFSURF, BETA_DIST, BETA_NOISE, COV_NOISE, DIST_AREA, FINGER_TRANS, 
;      INV_BLEACH, N_DOF, ONEDIM_BETA_DIST, ONEDIM_PROJECT_CONF, RESID_SUMSQ, 
;      TILDE_DEP, TILDE_INDEP, Z_POS
;
; USES:
;    cpar_ols.pro
;    cpar_tls.pro
;    pca.pro
;    regols.pro
;    regtls.pro
;
; PROCEDURE:
;    This function fits the model
;      y = (X+U_x)b - u_y,
;    where
;      \expect{u_y u_y^T} = C_N
;      C_N^{-1} = P^T P,
;    where P is a "pre-whitening" operator, so
;      \expect{P u_y u_y^T P^T} = I_n'
;    (i.e. P u_y = z is unit variance white noise) and
;      \expect{P U_x S^{-2} U_x^T P^T} = m * I_n',
;    where
;      S[k,k]^2 = frac_noise_var[k].
;
; REFERENCES:
;    Allen, M. R., and S. F. B. Tett.  1999.  Checking internal consistency in 
;        optimal fingerprinting.  Climate Dynamics, 15, 419-434.
;    Allen, M. R., and P. A. Stott.  2003.  Estimating signal amplitudes in 
;        optimal fingerprinting. Part I: theory.  Climate Dynamics, 21, 477-491.
;    Huntingford, C., P. A. Stott, M. R. Allen, and F. H. Lambert.  2006.  
;        Incorporating model uncertainty into attribution of observed 
;        temperature change.  Geophysical Research Letters, 33, L05710, 
;        10.1029/2005GL024831.
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
;    Modified:  MRA, 1999-08-02 (Bug-fix on rssq MRA, bug spotted by PAS; v1.1)
;    Modified:  MRA, 1999-08-09 (Z_poss keyword implemented; v1.3)
;    Modified:  MRA, 2000-08-03 (B_dist keyword implemented: explicit output 
;        of PDFs;  v2.0)
;    Modified:  Daithi Stone (stoned@atm.ox.ac.uk), 2004-06-28 
;        (Documentation for inclusion in routine library)
;    Modified:  DAS, 2005-03-13 (Updated documentation;  added PV_WAVE 
;        keyword;  updated compliance with svdpvw.pro;  cpar_ols.pro;  
;        cpar_tls.pro)
;    Modified:  DAS, 2005-05-13 (Updated documentation)
;    Modified:  DAS, 2005-09-01 (Added B1DIM, BT1DIM, DOUBLE keywords;  
;        changed the default NPOINT; updated documentation;  v3.0.0)
;    Modified;  DAS, 2010-02-12 (Added H1DIM, TYPE keywords;  edited 
;        documentation formatting;  added PV_WAVE note)
;    Modified:  DAS, 2010-03-23 (Implemented MODVAR and Z_DIST keywords)
;    Modified:  DAS, 2011-10-14 (Discontinued PV-WAVE compatibility;  removed 
;        pre-v3.0 header documentation;  switched used of svdpvw.pro to 
;        pca.pro;  removed obsolete NDOFN and NONPAR keywords, and STATUS 
;        keyword;  clarified in documentation than N_DOF is an output only;  
;        altered code variable names;  changed input and output variable names 
;        (see translation notes below);  v3.1.0)
;    Modified:  DAS, 2012-01-24 (Switched to SVD option with pca.pro;  blocked 
;        use of EIV approach;  v3.1.1)
;    Modified:  DAS, 2012-02-15 (Added WEIGHT keyword;  v3.1.2)
;-

;***********************************************************************

FUNCTION LINMOD, $
    DATA_INDEP, DATA_DEP, $
    DATA_NOISE=data_noise, $
    ADDVAR_DEP=addvar_dep, ADDVAR_INDEP=addvar_indep, $
    BETA_CONFSURF=beta_confsurf, BETA_DIST=beta_dist, $
    BETA_NOISE=beta_noise, $
    BLEACH=bleach, INV_BLEACH=bleach_inv, $
    COV_NOISE=cov_noise, $
    DIST_AREA=dist_area, $
    FRAC_NOISE_VAR=frac_noise_var, $
    T_CRIT=t_crit, $
    FINGER_TRANS=finger_trans, $
    N_DOF=n_dof, $
    N_POINT=n_point, $
    ONEDIM_BETA_DIST=onedim_beta_dist, ONEDIM_BETA_T=onedim_beta_t, $
      ONEDIM_INDEP=onedim_indep, ONEDIM_PROJECT_CONF=onedim_project_conf, $
      ONEDIM_PROJECT_DIST=onedim_project_dist, $
    RESID_SUMSQ=resid_sumsq, $
    T_DIST=t_dist, $
    TILDE_DEP=tilde_dep, TILDE_INDEP=tilde_indep, $
    TYPE=type, OLS=ols_opt, $
    WEIGHT=weight, $
    Z_DIST=z_dist, Z_POS=z_pos, $
    DOUBLE=double_opt, $
    SCALE_COV_NOISE=scale_cov_noise_opt

; Changes in input/output variable names from v3.0.0 to v3.1.0
;   changed b_dist to beta_dist
;   changed b1dim to onedim_beta_dist
;   changed betatl to beta_tilde
;   changed betaun to beta_noise
;   changed bt1dim to onedim_beta_t
;   changed c1dim to onedim_project_conf
;   changed cintvl to beta_confsurf
;   changed colour to inv_bleach
;   changed covb2s to cov_noise
;   changed d1dim to onedim_indep
;   changed estvar to scale_cov_noise
;   changed ftrans to finger_trans
;   changed h1dim to onedim_project_dist
;   changed modvar to addvar_indep
;   changed ndofd to n_dof
;   changed obsvar to addvar_dep
;   changed p_area to dist_area
;   changed rssq to resid_sumsq
;   changed unoise to data_noise
;   changed x to data_indep
;   changed xnoise to frac_noise_var
;   changed xtilde to tilde_indep
;   changed y to data_dep
;   changed ytilde to tilde_dep
;   changed z_poss to z_pos

;***********************************************************************
; Options

; Determine whether we are in GDL
gdl_opt = !prompt eq 'GDL> '

; Default regression type
if not( keyword_set( type ) ) then begin
  ; Default to old ols_opt request
  if keyword_set( ols_opt ) then begin
    type = 'OLS'
  ; Otherwise check if frac_noise_var is given
  endif else if keyword_set( frac_noise_var ) then begin
    ; Check if addvar_indep is given
    if keyword_set( addvar_indep ) then begin
      print, 'linmod.pro:  Assuming EIV method desired.'
      type = 'EIV'
    endif else begin
      print, 'linmod.pro:  Assuming TLS method desired.'
      type = 'TLS'
    endelse
  endif else begin
    print, 'linmod.pro:  Assuming OLS method desired.'
    type = 'OLS'
  endelse
endif
type = strupcase( type )
; Ensure supported regression type
if max( type eq [ 'OLS', 'TLS', 'EIV' ] ) ne 1 then begin
  stop, 'linmod.pro:  Unsupported regression type'
endif

; Block use of EIV approach
if type eq 'EIV' then begin
  stop, 'linmod.pro:  ' $
      + 'Use of the EIV approach has been disabled in this version.'
  ; Weighting has not been tested with EIV yet
  if keyword_set( weight ) then begin
    stop, 'linmod.pro:  Weighting has not been tested with EIV yet.'
  endif
endif

; Option for double precision
double_opt = keyword_set( double_opt )
one = 1.
if double_opt eq 1 then one = double( one )

; Option to scale the covariance estimates
scale_cov_noise_opt = keyword_set( scale_cov_noise_opt )

;***********************************************************************
; Set constants and variables

; Number of data points
n_data = n_elements( data_dep )
; Number of independent variables
n_indep = n_elements( data_indep[0,*] )
; Number of noise realisations
n_noise = n_elements( data_noise ) / n_data

; Initialise addvar_dep
if not( keyword_set( addvar_dep ) ) then addvar_dep = 0.

; Set default t_crit
if not( keyword_set( t_crit ) ) then t_crit = sqrt( findgen( n_indep ) + 1. )

; Set default number of points on beta_confsurf
if not( keyword_set( n_point ) ) then begin
  if n_indep eq 1 then begin
    n_point = 2
  endif else if n_indep eq 2 then begin
    n_point = 32767
  endif else begin
    n_point = fix( 32767. ^ ( 1. / ( n_indep - 1 ) ) )
  endelse
endif

; Check that all elements of frac_noise_var are non-zero (i.e. not infinite 
; sample size)
if keyword_set( frac_noise_var ) then begin
  id = where( frac_noise_var lt 1.e-5, n_id )
  if n_id gt 0 then begin
    print, 'linmod.pro:  frac_noise_var input with zeros.  Resetting to 1e-5.'
    frac_noise_var[id] = 1.e-5
  endif
endif
; Calculate standard deviation of sampling noise
if keyword_set( frac_noise_var ) then frac_noise_std = sqrt( frac_noise_var )

; Determine the number of degrees of freedom in input data
if keyword_set( bleach ) then begin
  ; Set to rank of bleach
  n_dof = n_elements( bleach[*,0] )
  if n_dof gt n_data then stop, 'linmod.pro:  Bleach operator too large'
endif else begin
  ; Set to size of input data
  n_dof = n_data
endelse

;***********************************************************************
; Initialise output arrays

; Flag optional outputs for output if requested
tilde_indep_bleach = keyword_set( tilde_indep )
tilde_dep_bleach = keyword_set( tilde_dep )
finger_trans_bleach = keyword_set( finger_trans )

; If pre-whitening is required
if keyword_set( bleach ) then begin
  ; Perform pre-whitening transformation on data arrays
  data_dep_bleach = bleach # data_dep
  data_indep_bleach = bleach # data_indep
  if keyword_set( data_noise ) then begin
    data_noise_bleach = bleach # data_noise
  endif else begin
    data_noise_bleach = 0.
  endelse
  ; Perform pre-whitening transformation on addvar_dep and addvar_indep
  if keyword_set( addvar_dep ) then begin
    addvar_dep_bleach = bleach # addvar_dep # transpose( bleach )
  endif
  if keyword_set( addvar_indep ) then begin
    addvar_indep_bleach = fltarr( n_dof, n_dof, n_indep )
    for i_indep = 0, n_indep - 1 do begin
      addvar_indep_bleach[*,*,i_indep] = bleach # addvar_indep[*,*,i_indep] $
          # transpose( bleach )
    endfor
  endif
  ; Calculate the pseudo-inverse of the pre-whitening transformation of needed
  if keyword_set( tilde_indep ) or keyword_set( tilde_dep ) $
      or keyword_set( z_pos ) or keyword_set( bleach_inv ) then begin
    pca, bleach, evalue=bleach_eigval, evector=bleach_eigvec, $
        pc=bleach_eigload, double=double_opt, normalise_pc=1, no_anomaly=1, $
        svd=gdl_opt
    bleach_eigval = sqrt( n_data * bleach_eigval )
    id = where( bleach_eigval le 0., n_id )
    if n_id gt 0 then stop, 'linmod.pro:  Rank-deficient pre-whitening operator'
    for i_dof = 0, n_dof - 1 do begin
      bleach_eigvec[*,i_dof] = bleach_eigvec[*,i_dof] / bleach_eigval[i_dof]
    endfor
    bleach_eigval = 0
    bleach_inv = bleach_eigload # transpose( bleach_eigvec )
    bleach_eigvec = 0
    bleach_eigload = 0
  endif
; Otherwise copy inputs to working arrays
endif else begin
  ; Copy data arrays
  data_dep_bleach = data_dep
  data_indep_bleach = data_indep
  if keyword_set( data_noise ) then begin
    data_noise_bleach = data_noise
  endif else begin
    data_noise_bleach = 0.
  endelse
  ; Copy variance adjustment arrays
  if keyword_set( addvar_dep ) then addvar_dep_bleach = addvar_dep
  if keyword_set( addvar_indep ) then addvar_indep_bleach = addvar_indep
endelse

; Weight input data arrays
if keyword_set( weight ) then begin
  ; Normalise weighting
  weight_use = weight / total( weight )
  ; Weight input data arrays
  data_dep_bleach = data_dep_bleach * weight_use
  for i_indep = 0, n_indep - 1 do begin
    data_indep_bleach[*,i_indep] = weight_use * data_indep_bleach[*,i_indep]
  endfor
  if keyword_set( data_noise_bleach ) then begin
    for i_noise = 0, n_noise - 1 do begin
      data_noise_bleach[*,i_noise] = weight_use * data_noise_bleach[*,i_noise]
    endfor
  endif
endif

;***********************************************************************
; Perform regression

; For TLS or EIV
if max( type eq [ 'TLS', 'EIV' ] ) eq 1 then begin

  ; Interate through independent variables
  for i_indep = 0, n_indep - 1 do begin
    ; Scale independent data by frac_noise_std
    data_indep_bleach[*,i_indep] = data_indep_bleach[*,i_indep] $
        / frac_noise_std[i_indep]
    ; Scale estimates of noise in patterns by frac_noise_std
    if keyword_set( addvar_indep ) then begin
      addvar_indep_bleach[*,*,i_indep] = addvar_indep_bleach[*,*,i_indep] $
          / frac_noise_var[i_indep]
    endif
    ; Scale directions for 1-D confidence intervals by frac_noise_std
    if keyword_set( onedim_indep ) then begin
      onedim_indep[i_indep,*] = onedim_indep[i_indep,*] $
          / frac_noise_std[i_indep]
    endif
  endfor
  ; Perform regression on prewhitened data and scaled i.v.s
  beta_tilde = regtls( data_indep_bleach, data_dep_bleach, $
      data_noise=data_noise_bleach, addvar_dep=addvar_dep_bleach, $
      addvar_indep=addvar_indep_bleach, scale_cov_noise=scale_cov_noise_opt, $
      tilde_indep=tilde_indep_bleach, tilde_dep=tilde_dep_bleach, $
      resid_sumsq=resid_sumsq, finger_trans=finger_trans_bleach, $
      beta_noise=beta_noise, cov_noise=cov_noise, eigload=eigload, $
      eigval=eigval, double=double_opt )
  ; If we need to compute confidence intervals
  if keyword_set( beta_confsurf ) or keyword_set( onedim_indep ) $
      or keyword_set( z_pos ) or keyword_set( onedim_beta_t ) then begin
    ; Create array containing both dependent and independent data if necessary
    if keyword_set( z_pos ) then begin
      z_data = transpose( [ transpose( data_indep_bleach ), $
          transpose( data_dep_bleach ) ] )
    endif
    ; Estimate the confidence intervals of the regression coefficients
    beta_confsurf = cpar_tls( beta_tilde, eigload, eigval, t_crit, n_point, $
        onedim_indep=onedim_indep, onedim_project_conf=onedim_project_conf, $
        z_data=z_data, z_pos=z_pos, onedim_beta_dist=onedim_beta_dist, $
        onedim_beta_t=onedim_beta_t, onedim_project_dist=onedim_project_dist, $
        double=double_opt )
    z_data = 0
  endif
  ; If likelihood distributions of the regression coefficients have been 
  ; requested
  if keyword_set( beta_dist ) then begin
    ; Create array containing both dependent and independent data if necessary
    if keyword_set( z_dist ) then begin
      z_data = transpose( [ transpose( data_indep_bleach ), $
          transpose( data_dep_bleach ) ] )
    endif
    ; Estimate the likelihood distributions of the regression coefficients
    beta_dist = cpar_tls( beta_tilde, eigload, eigval, t_dist, n_point, $
        dist_area=dist_area, z_data=z_data, z_pos=z_dist, double=double_opt )
    z_data = 0
  endif
  ; Rescale everything by frac_noise_std
  for i_indep = 0, n_indep - 1 do begin
    ; Multiply reconstruction by frac_noise_std since original data was divided 
    ; by it
    if keyword_set( tilde_indep ) then begin
      tilde_indep_bleach[*,i_indep] = tilde_indep_bleach[*,i_indep] $
          * frac_noise_std[i_indep]
    endif
    if keyword_set( z_pos ) then begin
      z_pos[*,i_indep,*,*] = z_pos[*,i_indep,*,*] * frac_noise_std[i_indep]
    endif
    if keyword_set( z_dist ) then begin
      z_dist[*,i_indep,*,*] = z_dist[*,i_indep,*,*] * frac_noise_std[i_indep]
    endif
    if keyword_set( onedim_indep ) then begin
      onedim_indep[i_indep,*] = onedim_indep[i_indep,*] $
          * frac_noise_std[i_indep]
    endif
    ; Divide regression coefficients by frac_noise_std, to correspond to $
    ; rescaled data
    beta_tilde[i_indep] = beta_tilde[i_indep] / frac_noise_std[i_indep]
    if keyword_set( beta_dist ) then begin
      beta_dist[i_indep,*,*] = beta_dist[i_indep,*,*] / frac_noise_std[i_indep]
    endif
    if keyword_set( beta_noise ) then begin
      beta_noise[i_indep,*] = beta_noise[i_indep,*] / frac_noise_std[i_indep]
    endif
    ; Divide rows and columns of covariance estimates by frac_noise_std, to 
    ; correspond to rescaled data
    if keyword_set( cov_noise ) then begin
      cov_noise[*,i_indep] = cov_noise[*,i_indep] / frac_noise_std[i_indep]
      cov_noise[i_indep,*] = cov_noise[i_indep,*] / frac_noise_std[i_indep]
    endif
    ; Divide rows of beta_confsurf matrix by frac_noise_std (scaled as 
    ; beta_tilde)
    if keyword_set( beta_confsurf ) then begin
      beta_confsurf[i_indep,*,*] = beta_confsurf[i_indep,*,*] $
          / frac_noise_std[i_indep]
    endif
  endfor

; For OLS regression
endif else if type eq 'OLS' then begin

  ; Perform standard OLS regression on prewhitened data
  beta_tilde = regols( data_indep_bleach, data_dep_bleach, $
      data_noise=data_noise_bleach, addvar_dep=addvar_dep_bleach, $
      scale_cov_noise=scale_cov_noise_opt, tilde_dep=tilde_dep_bleach, $
      resid_sumsq=resid_sumsq, finger_trans=finger_trans_bleach, $
      beta_noise=beta_noise, cov_noise=cov_noise, double=double_opt )
  tilde_indep_bleach = data_indep_bleach
  if keyword_set( frac_noise_var ) then begin
    ; Inflate rows and columns of cov_noise to account for noise in 
    ; data_indep_bleach (approximate)
    cov_noise = cov_noise * ( sqrt( 1. + frac_noise_var ) $
        # transpose( sqrt( 1. + frac_noise_var ) ) )
    ; Deflate the residual by the mean noise in data_indep_bleach (approximate)
    ; (Tests are applied as if this noise is not present)
    resid_sumsq = resid_sumsq / ( 1. + total( frac_noise_var ) / n_indep )
  endif
  ; If we need to compute confidence intervals
  if keyword_set( beta_confsurf ) or keyword_set( onedim_indep ) $
      or keyword_set( onedim_beta_t ) then begin
    ; Create array containing both dependent and independent data if necessary
    if keyword_set( z_pos ) then begin
      z_data = transpose( [ transpose( data_indep_bleach ), $
          transpose( data_dep_bleach ) ] )
    endif
    ; Estimate the confidence intervals of the regression coefficients
    beta_confsurf = cpar_ols( beta_tilde, cov_noise, t_crit, n_point, $
        onedim_indep=onedim_indep, onedim_project_conf=onedim_project_conf, $
        z_data=z_data, z_pos=z_pos, onedim_beta_dist=onedim_beta_dist, $
        onedim_beta_t=onedim_beta_t, onedim_project_dist=onedim_project_dist, $
        double=double_opt )
    z_data = 0
  endif
  ; If likelihood distributions of the regression coefficients have been 
  ; requested
  if keyword_set( beta_dist ) then begin
    ; Create array containing both dependent and independent data if necessary
    if keyword_set( z_dist ) then begin
      z_data = transpose( [ transpose( data_indep_bleach ), $
          transpose( data_dep_bleach ) ] )
    endif
    ; Estimate the likelihood distributions of the regression coefficients
    beta_dist = cpar_ols( beta_tilde, cov_noise, t_dist, n_point, $
        dist_area=dist_area, z_data=z_data, z_pos=z_dist )
    z_data = 0
  endif
endif

; Compute best-fit data_indep, data_dep and fingerprint matrix in original 
; (unbleached) coordinates if bleach was given
if keyword_set( bleach ) then begin
  if keyword_set( tilde_dep ) then tilde_dep = bleach_inv # tilde_dep_bleach
  if keyword_set( tilde_indep ) then begin
    tilde_indep = bleach_inv # tilde_indep_bleach
  endif
  if keyword_set( finger_trans ) then begin
    finger_trans = finger_trans_bleach # bleach
  endif
  if keyword_set( z_pos ) then begin
    z_pos = reform( temporary( z_pos ), n_dof, (n_indep+1)*n_point*n_indep )
    z_pos = temporary( bleach_inv # z_pos )
    z_pos = reform( temporary( z_pos ), n_data, n_indep+1, n_point, n_indep )
  endif
  if keyword_set( z_dist ) then begin
    z_dist = reform( z_dist, n_dof, (n_indep+1)*n_point*n_dist )
    z_dist = bleach_inv # z_dist
    z_dist = reform( z_dist, n_data, n_indep+1, n_point, n_dist )
  endif
; Or just copy the variables if bleach was not given
endif else begin
  if keyword_set( tilde_dep ) then tilde_dep = tilde_dep_bleach
  if keyword_set( tilde_indep ) then tilde_indep = tilde_indep_bleach
  if keyword_set( finger_trans ) then finger_trans = finger_trans_bleach
endelse

;***********************************************************************
; The end

return, beta_tilde
END
