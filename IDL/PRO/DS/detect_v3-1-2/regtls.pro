;+
; NAME:
;    REGTLS
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
;    This function performs multiple total least squares (TLS) regression, 
;    including the extension of the error-in-variables (EIV) approach.
;
; CATEGORY:
;    Optimal Detection Package v3.1.2
;
; NOTES:
;    v3.1.0 is no longer compatible with PV-WAVE.  v2.1 is the most recent 
;    version tested on PV-WAVE.
;
; CALLING SEQUENCE:
;    Result = regtls( DATA_INDEP, DATA_DEP )
;
; INPUTS:
;    DATA_INDEP:  A floating point matrix of size N_DATA*N_INDEP containing the 
;        N_DATA values of N_INDEP independent variables.
;    DATA_DEP:  A floating point vector of size N_DATA containing the N_DATA 
;        values of the dependent variable.
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
;    BETA_NOISE:  Returns an N_INDEP*N_NOISE matrix containing estimates of the 
;        N_INDEP coefficients estimated from each of the N_NOISE different 
;        noise realisations in DATA_NOISE added to the dependent variable 
;        DATA_DEP.
;    COV_NOISE:  If set then returns an N_INDEP*N_INDEP matrix containing the 
;        estimated approximate covariance of the noise realisations in 
;        DATA_NOISE in the directions of the N_INDEP scenario patterns found in 
;        the N_INDEP independent variables of DATA_INDEP.
;    DATA_NOISE:  A floating point matrix of size N_DATA*N_NOISE containing 
;        N_NOISE independent noise realisations of the dependent variable 
;        DATA_DEP for use in estimating confidence intervals.
;    DOUBLE:  If set then calculations are done in double precision arithmetic. 
;        The default is single precision.
;    EIGLOAD:  Returns a floating point matrix of size (N_INDEP+1)*(N_INDEP+1) 
;        containing the N_INDEP+1 loadings corresponding to the eigenvectors of 
;        transpose([[DATA_INDEP],[DATA_DEP]])#[[DATA_INDEP],[DATA_DEP]].
;    EIGVAL:  Returns a floating point vector of size N_INDEP+1 containing the 
;        N_INDEP+1 eigenvalues of 
;        transpose([[DATA_INDEP],[DATA_DEP]])#[[DATA_INDEP],[DATA_DEP]].
;    FINGER_TRANS:  If set then this returns the transpose of the matrix of 
;        fingerprint patterns.  Of size N_INDEP*N_DATA.
;    RESID_SUMSQ:  Returns the sum of squares of the residual differences 
;        between the prewhitened dependent variable Y and the best regression 
;        estimate.  If DATA_NOISE is given, then the normalised residuals are 
;        used, with the normalisation done for each element of DATA_DEP 
;        according to the standard deviation of the corresponding N_NOISE 
;        realisations of that element in DATA_NOISE.
;    SCALE_COV_NOISE:  If set, then the covariance estimates are scaled by 
;        dividing by the residual sum of squares (RESID_SUMSQ).
;    STATUS:  Returns 0 if the routine has completed without problems, 1 
;        otherwise.
;    TILDE_DEP:  Returns a floating point vector of size N_DATA containing the 
;        N_DATA predicted values of the DATA_DEP variable from the best fit 
;        model.
;    TILDE_INDEP:  Returns a floating point matrix of size N_DATA*N_INDEP 
;        containing the N_DATA predicted values of the N_INDEP variables in 
;        DATA_INDEP from the best fit model.
;
; OUTPUTS:
;    Result:  A floating point vector of length N_INDEP containing the 
;        estimated regression coefficients for the N_INDEP variables in 
;        DATA_INDEP.
;    BETA_NOISE, COV_NOISE, EIGLOAD, EIGVAL, FINGER_TRANS, RESID_SUMSQ, STATUS, 
;        TILDE_DEP, TILDE_INDEP
;
; USES:
;    pca.pro
;    shuffle.pro
;    sign.pro
;
; PROCEDURE:
;    This function fits to the model y = (X+U_x)b - u_y
;    where \expect{u_y u_y^T} = I
;    and \expect{U_x U_x^T} = m * I
;    by applying singular value decomposition on the data matrix Z=[X,y].
;
; REFERENCES:
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
;    Modified:    Daithi Stone (stoned@atm.ox.ac.uk), 2004-06-28 (Documentation 
;        for inclusion in routine library)
;    Modified:    DAS, 2005-09-01 (Added DOUBLE keyword;  fixed bug in 
;        calculation of ZN;  updated documentation;  v3.0)
;    Modified:  DAS, 2011-11-06 (Added EIV capability based on code by Chris 
;        Huntingford;  edited documentation formating;  discontinued PV-WAVE 
;        compatibility and added note;  changed use of total_1d.pro to IDL's 
;        total;  switched use of svdpvw.pro to pca.pro;  removed pre-v3.0 
;        header documentation;  added use of shuffle.pro and sign.pro;  imposed 
;        non-replacement in random selection in non-parametric estimate of 
;        parameter spread;  altered code variable names;  made return of 
;        EIGLOAD and EIGVAL automatic;  changed input and output variable names 
;        (see translation notes below);  corrected eigen terminology;  v3.1.0)
;    Modified:  DAS, 2012-01-24 (Switched to SVD option with pca.pro;  edited 
;        for compliance with GDL;  blocked use of EIV approach;  v3.1.1)
;    Modified:  DAS, 2012-09-26 (Fixed disabling of EIV approach;  v3.1.2)
;-

;***********************************************************************

FUNCTION REGTLS, $
    DATA_INDEP, DATA_DEP, $
    DATA_NOISE=data_noise, $
    ADDVAR_DEP=addvar_dep, ADDVAR_INDEP=addvar_indep, $
    BETA_NOISE=beta_noise, COV_NOISE=cov_noise, $
    EIGLOAD=eigload, EIGVAL=eigval, $
    FINGER_TRANS=finger_trans, $
    RESID_SUMSQ=resid_sumsq, $
    STATUS=status, $
    TILDE_INDEP=tilde_indep, TILDE_DEP=tilde_dep, $
    DOUBLE=double_opt, $
    SCALE_COV_NOISE=scale_cov_noise_opt

; Changes in input/output variable names from v3.0.0 to v3.1.0
;   changed betatl to beta_tilde
;   changed betaun to beta_noise
;   changed covb2s to cov_noise
;   changed estvar to scale_cov_noise
;   changed evalus to eigval
;   changed evects to eigload
;   changed ftrans to finger_trans
;   changed modvar to addvar_indep
;   changed obsvar to addvar_dep
;   changed rssq to resid_sumsq
;   changed unoise to data_noise
;   changed x to data_indep
;   changed xtilde to tilde_indep
;   changed y to data_dep
;   changed ytilde to tilde_dep

;***********************************************************************
; Options

; Determine whether we are in GDL
gdl_opt = !prompt eq 'GDL> '

; Option for EIV
if keyword_set( addvar_dep ) or keyword_set( addvar_indep ) then begin
  eiv_opt = 1
endif else begin
  eiv_opt = 0
endelse

; Block use of EIV approach
if eiv_opt eq 1 then begin
  stop, 'regtls.pro:  ' $
      + 'Use of the EIV approach has been disabled in this version.'
endif

; Option for double precision
double_opt = keyword_set( double_opt )
one = 1.
if double_opt eq 1 then one = double( one )

;***********************************************************************
; Set constants and variables

; Number of data points
n_data = n_elements( data_dep )
; Number of independent variables
n_indep = n_elements( data_indep[0,*] )
; Number of noise realisations
n_noise = n_elements( data_noise ) / n_data

; Create data=[data_indep,data_dep], remembering IDL's indexing convention
data_z = transpose( [ transpose( data_indep ), transpose( data_dep ) ] )

; Default addvar_dep and addvar_indep
if not( keyword_set( addvar_dep ) ) then addvar_dep = 0
if not( keyword_set( addvar_indep ) ) then addvar_indep = 0

;***********************************************************************
; Estimate regression coefficients

; Diagonalise data_z
pca, data_z, evalue=eigval, evector=eigvec, pc=eigload, double=double_opt, $
    normalise_pc=1, no_anomaly=1, svd=gdl_opt
eigval = ( n_indep + 1 ) * eigval
; Flip the signs to make eigload[n_indep,*] all negative
; (signs are arbitrary and it helps book-keeping later)
id = where( sign( eigload[n_indep,*] ) eq -1, n_id )
if n_id gt 0 then begin
  eigvec[*,id] = -eigvec[*,id]
  eigload[*,id] = -eigload[*,id]
endif

; Compute best-fit regression coefficients from the mininum eigenvector
y_coef = eigload[n_indep,n_indep]
if y_coef eq 0. then begin
  stop, 'regtls.pro:  Zero weight on observations in regtls:  ', $
      eigload[*,n_indep]
endif
beta_tilde = -reform( eigload[0:n_indep-1,n_indep] ) / y_coef

; Compute best-fit model from all singular values orthogonal to 
; eigload[*,n_indep]
if keyword_set( tilde_indep ) or keyword_set( tilde_dep ) then begin
  tilde_z = data_z $
      - data_z # eigload[*,n_indep] # transpose( eigload[*,n_indep] )
  tilde_indep = tilde_z[*,0:n_indep-1]
  tilde_dep = reform( tilde_z[*,n_indep] )
endif

; Set finger_trans to the operator which extracts beta_tilde from data_z
if keyword_set( finger_trans ) then begin
  finger_trans = transpose( eigvec[*,n_indep] )
endif

; Copy eigenvalues for use in residual calculations
eigval_tls = eigval

; If noise samples are given
if keyword_set( data_noise ) then begin
  ; Revise variance in each var_noise[i] to reflect independent noise estimate
  var_noise = total( data_noise ^ 2, 2 ) / n_noise
  ; Normalise eigenvalues for residual calculations with respect to projection 
  ; onto noise
  for i_indep = 0, n_indep do begin
    eigval_tls[i_indep] = eigval_tls[i_indep] $
        / ( transpose( eigvec[*,i_indep] ^ 2 ) # var_noise )
  endfor
  ; If this is EIV
  if eiv_opt eq 1 then begin
    ;; Revise eigen-decomposition given additional variance estimates
    ;; NOTE THE EIV OPTION IS CURRENTLY DISABLED.  IF IT IS REINSTATED THEN 
    ;; eivtls.pro NEEDS TO BE INCLUDED IN THE "USES" LIST ABOVE.
    ;temp_eigload = eivtls( data_z=data_z, tilde_z=tilde_z, eigvec=eigvec, $
    ;    eigload=eigload, eigval=eigval, addvar_dep=addvar_dep, $
    ;    addvar_indep=addvar_indep, double=double_opt )
    eigload = temporary( temp_eigload )
  endif
  ; Normalise eigenvalues with respect to projection onto noise
  for i_indep = 0, n_indep do begin
    eigval[i_indep] = eigval[i_indep] $
        / ( transpose( eigvec[*,i_indep] ^ 2 ) # var_noise )
  endfor
endif

; Compute sum-squared residual
resid_sumsq = eigval_tls[n_indep]

; Revise variance estimates to residuals if scale_cov_noise is set
if keyword_set( scale_cov_noise_opt ) then begin
  eigval = eigval / ( resid_sumsq / ( n_data - n_indep ) )
endif

; Check for degenerate eigenvalues for the independent variables
eigval_diff = eigval[0:n_indep-1] - eigval[n_indep]
id_indep = where( eigval_diff lt 0., n_id_indep )
if n_id_indep gt 0 then begin
  print, 'regtls.pro:  degenerate eigenvalues.  Resetting to 0.'
  status = 1.
  eigval[id_indep] = eigval[n_indep]
endif
; Identify the non-degenerate eigenvalues of the independent variables
id_indep = where( eigval_diff ge 0., n_id_indep )

; Estimate of covariance
if keyword_set( cov_noise ) then begin
  ; Initialise new eigenvector matrix
  temp_eigload = one * fltarr( n_indep+1, n_indep+1 )
  ; Normalise eigenvectors by sqrt(eigenvalue)
  for i_indep = 0, n_id_indep - 1 do begin
    temp_eigload[*,id_indep[i_indep]] = eigload[*,id_indep[i_indep]] $
        / sqrt( eigval_diff[id_indep[i_indep]] )
  endfor
  ; Compute covariance of PCs
  cov_eigload = temp_eigload # transpose( temp_eigload )
  ; Compute second-order approximation to covariance of beta_tilde
  eigload_min = reform( eigload[*,n_indep] )
  eigload_dep = eigload_min[n_indep]
  cov_noise = cov_eigload[0:n_indep-1,0:n_indep-1] / eigload_dep ^ 2 $
      - cov_eigload[0:n_indep-1,n_indep] $
      # transpose( [ eigload[0:n_indep-1] ] ) / eigload_dep ^ 3 $
      - eigload[0:n_indep-1] # cov_eigload[n_indep,0:n_indep-1] $
      / eigload_dep ^ 3 $
      + eigload[0:n_indep-1] # transpose( [ eigload[0:n_indep-1] ] ) $
      * cov_eigload[n_indep,n_indep] / eigload_dep ^ 4
  ; Clear memory
  temp_eigload = 0
  cov_eigload = 0
  eigload_min = 0
  eigload_dep = 0
endif

; Non-parametric estimate of parameter spread
if keyword_set( beta_noise ) then begin
  ; Create new eigenvector matrix for noise
  eigload_noise = one * fltarr( n_indep+1, n_indep+1 )
  ; Weight values according to sqrt(eigenvalues)
  for i_indep = 0, n_id_indep - 1 do begin
    eigload_noise[*,id_indep[i_indep]] = eigload_noise[*,id_indep[i_indep]] $
        * sqrt( eigval_diff[id_indep[i_indep]] )
  endfor
  ; Initialise estimate of parameter spread
  beta_noise = one * fltarr( n_indep, n_noise )
  ; Reconstruct noise-reduced data
  tilde_noise = eigvec # transpose( eigload_noise )
  ; Iterate through noise realisations
  for i_noise = 0, n_noise-1 do begin
    ; Copy data
    temp_tilde_noise = tilde_noise
    ; Contaminate observations, temp_tilde_noise_noise[*,n_indep], with noise 
    ; segment
    ; (Bug noted by Gareth Jones, fixed by DAS)
    temp_tilde_noise[*,n_indep] = temp_tilde_noise[*,n_indep] $
        + data_noise[*,i_noise]
    ; Extract independent noise segments
    ; (non-replacement imposed by DAS)
    index = indgen( n_noise )
    id = where( index ne i_noise )
    index = index[id]
    index = shuffle( index, seed=seed )
    ; Contaminate independent variables, temp_tilde_noise[0:n_indep-1,*], with 
    ; independent random noise segments
    temp_tilde_noise[0:n_indep-1,*] = temp_tilde_noise[0:n_indep-1,*] $
        + data_noise[*,index[0:n_indep-1]]
    ; Compute estimates of regression coefficients (beta_tilde) for noise 
    ; realisations
    pca, temp_tilde_noise, pc=temp_eigload, double=double_opt, normalise_pc=1, $
        no_anomaly=1, svd=gdl_opt
    beta_noise[*,i_noise] = temp_eigload[0:n_indep-1,n_indep] $
        / temp_eigload[n_indep,n_indep]
  endfor
  ; Clear memory
  temp_tilde_noise = 0
  eigload_noise = 0
  tilde_noise = 0
  temp_eigload = 0
endif

;***********************************************************************
; The end

return, beta_tilde
END
