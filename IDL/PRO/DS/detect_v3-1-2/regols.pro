;+
; NAME:
;    REGOLS
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
;    This function performs multiple ordinary least squares regression.
;
; CATEGORY:
;    Optimal Detection Package v3.1.2
;
; NOTES:
;    v3.1.0 is no longer compatible with PV-WAVE.  v2.1 is the most recent 
;    version tested on PV-WAVE.
;
; CALLING SEQUENCE:
;    Result = regols( DATA_INDEP, DATA_DEP )
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
;        measurements.
;    BETA_NOISE:  Returns an N_INDEP*N_NOISE matrix containing estimates of the 
;        N_INDEP coefficients estimated from each of the N_NOISE different 
;        noise realisations in DATA_NOISE added to the dependent variable 
;        DATA_DEP.
;    COV_NOISE:  If set then returns an N_INDEP*N_INDEP matrix containing the 
;        estimated covariance of the noise realisations in DATA_NOISE in the 
;        directions of the N_INDEP scenario patterns found in the N_INDEP 
;        independent variables of DATA_INDEP.  If SCALE_COV_NOISE is set then 
;        the values are scaled by dividing by the residual sum of squares 
;        (RESID_SUMSQ).
;    DATA_NOISE:  A floating point matrix of size N_DATA*N_NOISE containing 
;        N_NOISE independent noise realisations of the dependent variable 
;        DATA_DEP for use in estimating confidence intervals.
;    DOUBLE:  If set then calculations are done in double precision arithmetic. 
;        The default is single precision.
;    FINGER_TRANS:  Returns the transpose of the matrix of fingerprint 
;        patterns.  Note that Result=FINGER_TRANS#DATA_DEP.  Of size 
;        N_INDEP*N_DATA.
;    RESID_SUMSQ:  Returns the sum of squares of the residual differences 
;        between the prewhitened dependent variable DATA_DEP and the best 
;        regression estimate.  If DATA_NOISE is given, then the normalised 
;        residuals are used, with the normalisation done for each element of 
;        DATA_DEP according to the standard deviation of the corresponding 
;        N_NOISE realisations of that element in DATA_NOISE.
;    SCALE_COV_NOISE:  If set, then the covariance estimates in COV_NOISE are 
;        scaled by the residual sum of squares (RESID_SUMSQ).
;    TILDE_DEP:  Returns a floating point vector of size N_DATA containing the 
;        N_DATA predicted values of the DATA_DEP variable from the best fit 
;        model.
;
; OUTPUTS:
;    Result:  A floating point vector of length N_INDEP containing the 
;        estimated regression coefficients for the N_INDEP variables in 
;        DATA_INDEP.
;    BETA_NOISE, COV_NOISE, FINGER_TRANS, RESID_SUMSQ, TILDE_DEP
;
; USES:
;    -
;
; PROCEDURE:
;    This function computes the Hessian explicitly and then inverts it.
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
;    Modified:    MRA, 2000-01-11 (Modified name of invert1 to invert1k;  v1.1)
;    Modified:    Daithi Stone (stoned@atm.ox.ac.uk), 2004-06-28 
;        (Documentation for inclusion in routine library)
;    Modified:    DAS, 2005-09-01 (Added DOUBLE keyword;  updated 
;        documentation;  v3.0.0)
;    Modified:    DAS, 2010-02-12 (Edited documentation formatting;  added 
;        PV_WAVE note)
;    Modified:  DAS, 2011-11-06 (Discontinued PV-WAVE compatibility;  changed 
;        use of total_1d.pro and invert1k.pro to IDL's total and invert;  
;        removed pre-v3.0 header documentation;  removed STATUS keyword;  
;        altered code variable names;  changed input and output variable names 
;        (see translation notes below);  v3.1.0)
;-

;***********************************************************************

FUNCTION REGOLS, $
    DATA_INDEP, DATA_DEP, $
    DATA_NOISE=data_noise, $
    ADDVAR_DEP=addvar_dep, $
    BETA_NOISE=beta_noise, COV_NOISE=cov_noise, $
    FINGER_TRANS=finger_trans, $
    RESID_SUMSQ=resid_sumsq, $
    SCALE_COV_NOISE=scale_cov_noise_opt, $
    TILDE_dep=tilde_dep, $
    DOUBLE=double_opt

; Changes in input/output variable names from v3.0.0 to v3.1.0
;   changed betatl to beta_tilde
;   changed betaun to beta_noise
;   changed covb2s to cov_noise
;   changed estvar to scale_cov_noise
;   changed ftrans to finger_trans
;   changed obsvar to addvar_dep
;   changed rssq to resid_sumsq
;   changed unoise to data_noise
;   changed x to data_indep
;   changed y to data_dep
;   changed ytilde to tilde_dep

;***********************************************************************
; Options

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
; n_noise=number of noise realisations
n_noise = n_elements( data_noise ) / n_data

; Perform standard OLS regression on input data
; (Forcing double precision substantially increases accuracy.)
cov_noise = invert( transpose( data_indep ) # data_indep, status, double=1 )
if double_opt eq 0 then cov_noise = float( cov_noise )
; Check that the inversion has gone okay
if status ne 0 then begin
  stop, 'regols.pro:  ' $
      + 'Failure in matrix inversion for estimating regression coefficients.'
endif

; Compute F^T, the matrix of fingerprints
finger_trans = cov_noise # transpose( data_indep )

; Compute best-fit regression coefficients
beta_tilde = finger_trans # data_dep

; Compute predicted noise-free observations
tilde_dep = data_indep # beta_tilde

; Compute residual
z_tilde = data_dep - tilde_dep
if n_noise ne 0 then begin
  ; Compute estimates of beta from independent noise realisations
  beta_noise = finger_trans # data_noise
  ; Re-estimate cov_noise from covariance of independent estimates
  cov_noise = ( beta_noise # transpose( beta_noise ) ) / n_noise
  ; Re-estimate noise variance in each element of z_tilde
  var_noise = total( data_noise ^ 2, 2 ) / n_noise
  ; Add estimated observation error
  if keyword_set( addvar_dep ) then begin
    for i_data = 0, n_data - 1 do begin
      var_noise[i_data] = var_noise[i_data] + addvar_dep[i_data,i_data]
    endfor
  endif
  ; Re-normalise residuals
  z_tilde = z_tilde / sqrt( var_noise )
  ; Add beta_tilde to beta_noise, to centre noise on estimate
  beta_noise = beta_noise + rebin( beta_tilde, n_indep, n_noise )
endif

; Compute sum squared residual
resid_sumsq = total( z_tilde ^ 2 )

; Add estimated observation error to cov_noise if supplied
if keyword_set( addvar_dep ) then begin
  cov_noise = cov_noise + finger_trans # addvar_dep # transpose( finger_trans )
endif

; Re-estimate covariances using resid_sumsq if required
if keyword_set( scale_cov_noise_opt ) then begin
  cov_noise = cov_noise * resid_sumsq / ( n_data - n_indep )
endif

;***********************************************************************
; The end

return, beta_tilde
END
