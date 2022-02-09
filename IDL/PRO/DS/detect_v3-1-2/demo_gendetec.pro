;+
; NAME:
;    DEMO_GENDETEC
;
; COPYRIGHT:
;    Copyright (2004) Daithi Stone, University of Oxford.
;    Revised (2011) Daithi Stone under contract to the U.S. Department of 
;    Energy's Office of Science, Office of Biological and Environmental 
;    Research and the U.S. National Oceanic and Atmospheric Administration's 
;    Climate Program Office, via the International Detection and Attribution 
;    Group.
;
; PURPOSE:
;    This function is a simple example demonstration of how to use the Optimal 
;    Detection Package.  All inputs are optional.
;
; CATEGORY:
;    Optimal Detection Package v3.1.2
;
; CALLING SEQUENCE:
;    demo_gendetec, [N_SCEN]
;
; INPUTS:
;    N_SCEN:  An optional integer defining the number of scenarios (independent 
;        variables in the regression) to use.  The default is 3.
;
; KEYWORD PARAMETERS:
;    AMP_EPSILON:  An optional floating point scalar containing the amplitude 
;        of the Gaussian sampling noise.  See code for the random default value.
;    AMP_PERT:  An optional floating point scalar containing the amplitude 
;        of the uncertainty in the model parameter across models.  This is only 
;        used by EIV.  Ignored if FRAC_PERT_MODEL is given.  See code for the 
;        random default value.
;    AMP_SIGNAL:  An optional floating point vector of length N_SCEN 
;        containing scalings by which to modify the N_SCEN signals.  See 
;        code for the random default values.
;    BETA_REAL:  An optional floating point vector of length N_SCEN containing 
;        the actual regression coefficient values.  See code for the random 
;        default values.
;    DIST_PROB:  An optional integer scalar giving the number of isopleths 
;        of probability to sample for the estimation of likelihood functions 
;        for the regression coefficients and other output of gendetec.pro.  See 
;        code for the default value.
;    MAX_N_SCEN_SAMP:  An optional integer scalar containing the maximum number 
;        of time series samples to use for estimating the scenario functions 
;        from the models for the random default setting of N_SCEN_SAMP.  See 
;        code for the default value.  Ignored if N_SCEN_SAMP is set.
;    MIN_N_SCEN_SAMP:  An optional integer scalar containing the minimum number 
;        of time series samples to use for estimating the scenario functions 
;        from the models for the random default setting of N_SCEN_SAMP.  See 
;        code for the default value.  Ignored if N_SCEN_SAMP is set.
;    N_POINT:  An optional integer scalar giving the number of points to sample 
;        along an isopleth of probability when estimating the multivariate 
;        probability density function of the scaling factors.  See gendetec.pro 
;        for more information.  See code for the default value.
;    N_HIST:  An optional integer defining the number of bins to use in 
;        histogram plots.
;    N_MODEL:  An optional integer giving the number of different models of the 
;        time series to sample when estimating the scenarios when using EIV.  
;        See code for the default value.
;    FRAC_PERT_MODEL:  An optional floating point vector of length N_MODEL 
;        containing the perturbations to the model parameter for use in 
;        sampling across models.  Zero implies no perturbation.  For use with 
;        EIV only.  See code for the default random value.  If set then this 
;        overrides AMP_PERT.
;    N_TEST:  An optional integer scalar containing the number of random times 
;        to run the gendetec code, with the resulting ensemble scanning 
;        uncertainty in the inputs.  If non-zero, additional plots are 
;        returned based on this output.  The default is zero.
;    N_NOISE_1:  An optional integer scalar giving the number of random time 
;        series to use for estimating the structure of the noise.  These 
;        samples are used to estimate the optimal fingerprints.  See code for 
;        the default value.
;    N_NOISE_2:  An optional integer scalar giving the number of random time 
;        series to use for estimating the structure of the noise.  These 
;        samples are used for significance testing.  A value of zero means that 
;        no random time series will be generated.  See code for the default 
;        value.
;    N_SCEN_SAMP:  An optional integer vector of length N_SCEN (for OLS and 
;        TLS) or array of size N_SCEN*N_MODEL (for EIV) giving the number of 
;        time series samples to use for estimating the scenario functions from 
;        the models.  For instance, for EIV if value [i,j]=4 gives then 4 time 
;        series are averaged to estimate the underlying function of scenario i 
;        for model j.  See code for the default random values.  If set then 
;        MIN_N_SCEN_SAMP and MAX_N_SCEN_SAMP are ignored. 
;    N_TIME:  An optional integer defining the length of the regressed vectors.
;    NO_OPTIMISE:  If set then no optimisation of the fingerprints is used.  
;        The default is to use optimisation.
;    P_LIMIT:  An optional floating point number containing the
;        two-sided p-value for confidence interval estimation.  The default is 
;        0.1.
;    SEED:  An optional integer defining the seed for the random number 
;        generator.
;    TRUNC:  An optional integer scalar defining the number of leading EOFs of 
;        the first noise sample to retain for estimating the fingerprints.  See 
;        gendetec.pro for more information and the default.
;    TYPE:  An optional string describing the type of regression algorithm to 
;        use.  Possibilities are:
;          'OLS' for ordinary least squares (OLS),
;          'TLS' for total least squares (TLS),
;          'EIV' for error in variables (EIV).
;        The default is TLS.
;    FILENAME_HTML:  An optional filename for an HTML file in which to output 
;        the results of the procedure.  The default is to output the results 
;        to the screen.
;
; OUTPUTS:
;    Some keyword parameters return their default value if none was input.
;
; USES:
;    betas_hist.pro
;    gendetec.pro
;    line_legend.pro  (see http://web.csag.uct.ac.za/~daithi/idl_lib)
;    pdf.pro  (see http://web.csag.uct.ac.za/~daithi/idl_lib)
;    ps_close.pro  (see http://web.csag.uct.ac.za/~daithi/idl_lib)
;    ps_open.pro  (see http://web.csag.uct.ac.za/~daithi/idl_lib)
;    str.pro  (see http://web.csag.uct.ac.za/~daithi/idl_lib)
;    string_from_vector.pro  (see http://web.csag.uct.ac.za/~daithi/idl_lib)
;
; PROCEDURE:
;    This function creates a multiple linear regression problem and then has 
;    gendetec.pro solve it.  Let us say we have an observed quantity DATA_OBS.  
;    We think that DATA_OBS is the weighted linear sum of N_SCEN scenarios plus 
;    some Gaussian "noise".  So 
;      DATA_OBS = DATA_SCEN # BETA_REAL + AMP_EPSILON * EPSILON.  
;    DATA_OBS is our observations of length N_TIME.  BETA_REAL is a vector 
;    containing the N_SCEN regression coefficients.  DATA_SCEN is a matrix of 
;    size N_TIME*N_SCEN containing the N_TIME values for each scenario.  
;    AMP_EPSILON * EPSILON is the residual from this regression model, which we 
;    suppose to be distributed as a Gaussian EPSILON of standard deviation 
;    AMP_EPSILON.  Now if our scenarios DATA_SCEN are estimated from some other 
;    model, then there is noise there too which we can reduce by sampling 
;    N_SCEN_SAMP times and averaging.
;
; REFERENCES:
;    Allen, M. R., and S. F. B. Tett.  1999.  Checking for model consistency in 
;        optimal fingerprinting. Climate Dynamics, 15, 419-434. 
;    Allen, M. R., and P. A. Stott.  2003.  Estimating signal amplitudes in 
;        optimal fingerprinting. Part I: theory.  Climate Dynamics, 21, 477-491.
;    Stott, P. A., M. R. Allen, and G. S. Jones.  2003.  Estimating signal 
;        amplitudes in optimal fingerprinting. Part II: application to general 
;        general circulation models.  Climate Dynamics, 21, 493-500.
;
; EXAMPLE:
;    OLS regression against one independent variable with randomly defined 
;    regression properties:
;      demo_gendetec, 1, type='OLS'
;    EIV regression against two independent variables with regression 
;    coefficients 0.8 and 1.3, unmodifed signal amplitude, noise standard 
;    deviation of 0.1, no optimisation, standard deviation in random 
;    perturbation of model parameter of 0.01, 100 random test samples, and 
;    output sent to file demo_gendetec.html:
;      demo_gendetec, 2, amp_epsilon=0.1, amp_signal=[1.,1.], $
;          beta_real=[0.8,1.3], n_test=100, no_optimise=1, type='EIV', $
;          filename_html='demo_gendetec.html', amp_pert=0.01
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone (stoned@atm.ox.ac.uk), 2004-06-28
;    Modified:  DAS, 2005-03-15 (added PDF estimation)
;    Modified:  DAS, 2005-09-01 (majorly expanded PDF plotting;  added N_Scen 
;        input;  v3.0.0 of Optimal Detection Package)
;    Modified:  DAS, 2010-02-11 (changed reporting of P_RESID to be a 
;        two-sided test;  edited documentation formatting)
;    Modified:  DAS, 2011-11-06 (standardised plotting window;  altered code 
;        variable names;  added AMP_EPSILON, AMP_PERT, AMP_SIGNAL, BETA_REAL, 
;        DIST_PROB, N_POINT, N_HIST, N_MODEL, FRAC_PERT_MODEL, N_NOISE_1, 
;        N_NOISE_2, N_SCEN_SAMP, MIN_N_SCEN_SAMP, MAX_N_SCEN_SAMP, N_TIME, 
;        NO_OPTIMISE, P_LIMIT, SEED, TRUNC, TYPE keywords;  added random 
;        sampling test via N_TEST keyword;  added HTML output via FILENAME_HTML 
;        keyword;  v3.1.0)
;    Modifed:  DAS, 2012-01-24 (Blocked use of EIV approach;  v3.1.1)
;-

;***********************************************************************

PRO DEMO_GENDETEC, $
    N_SCEN, $
    AMP_EPSILON=amp_epsilon, AMP_PERT=amp_pert, AMP_SIGNAL=amp_signal, $
    BETA_REAL=beta_real, $
    DIST_PROB=dist_prob, N_POINT=n_point, $
    N_HIST=n_hist, $
    N_MODEL=n_model, FRAC_PERT_MODEL=frac_pert_model, $
    N_TEST=n_test, $
    N_NOISE_1=n_noise_1, N_NOISE_2=n_noise_2, $
    N_SCEN_SAMP=n_scen_samp, MIN_N_SCEN_SAMP=n_scen_samp_min, $
        MAX_N_SCEN_SAMP=n_scen_samp_max, $
    N_TIME=n_time, $
    NO_OPTIMISE=no_optimise_opt, $
    P_LIMIT=p_limit, $
    SEED=seed, $
    TRUNC=trunc, $
    TYPE=type, $
    FILENAME_HTML=filename_html

;***********************************************************************
; Constants

; Assume TLS if no type specified
if not( keyword_set( type ) ) then type = 'TLS'

; Block use of EIV approach
if type eq 'EIV' then begin
  stop, 'demo_gendetec.pro:  ' $
      + 'Use of the EIV approach has been disabled in this version.'
endif

; Initialise the random seed (this default allows reproducibility)
if n_elements( seed ) eq 0 then seed = 2

; Length of space-time series
if not( keyword_set( n_time ) ) then n_time = 100
; Number of scenarios (polynomial or sinusoidal functions)
if not( keyword_set( n_scen ) ) then n_scen = 3
; The number of models used to estimate the scenarios
if type eq 'EIV' then begin
  if not( keyword_set( n_model ) ) then n_model = 5
endif else begin
  ; Set a default n_model of 1 (simplifies coding later)
  n_model = 1
endelse
; Ensure legal input
if keyword_set( frac_pert_model ) then begin
  if n_elements( frac_pert_model ) ne n_model then begin
    stop, 'demo_gendetec.pro:  frac_pert_model needs to have n_model values.'
  endif
endif
; Number of samples of each scenario
if not( keyword_set( n_scen_samp ) ) then begin
  if not( keyword_set( n_scen_samp_min ) ) then n_scen_samp_min = 3
  if not( keyword_set( n_scen_samp_max ) ) then begin
    n_scen_samp_max = n_scen_samp_min + 4
  endif
  n_scen_samp = n_scen_samp_min $
      + round( randomu( seed, n_scen, n_model ) $
      * ( n_scen_samp_max - n_scen_samp_min ) )
endif else begin
  ; Ensure no zero-sized samples
  if min( n_scen_samp ) lt 1 then begin
    stop, 'demo_gendetec.pro:  n_scen_samp sets zero-sized samples.'
  endif
  ; Ensure proper size
  if n_elements( n_scen_samp ) ne n_scen * n_model then begin
    stop, 'demo_gendetec.pro:  n_scen_samp must have n_scn*n_model values'
  endif
endelse
n_scen_samp = reform( n_scen_samp, n_scen, n_model )
; Number of control samples for estimating noise
if not( keyword_set( n_noise_1 ) ) then n_noise_1 = 90
; Number of independent control samples for use in significance testing
if n_elements( n_noise_2  ) eq 0 then n_noise_2 = n_noise_1

; Actual regression coefficient values
if keyword_set( beta_real ) then begin
  if n_elements( beta_real ) ne n_scen then begin
    stop, 'demo_gendetec.pro:  beta_real must have n_scen values.'
  endif
endif else begin
  beta_real = 3. * randomu( seed, n_scen )
endelse
; Amplitude of noise in the variable
if n_elements( amp_epsilon ) eq 0 then amp_epsilon = 0.1 + 0.2 * randomu( seed )
; Amplitude of the parameter uncertainty in the models
if ( type eq 'EIV' ) and ( n_elements( amp_pert ) eq 0 ) then begin
  amp_pert = 0.1 + 0.2 * randomu( seed )
endif
; Amplitude of the signals in the variable
if n_elements( amp_signal ) eq 0 then begin
  amp_signal = 1. + randomu( seed, n_scen )
endif else begin
  if n_elements( amp_signal ) ne n_scen then begin
    stop, 'demo_gendetec.pro:  amp_signal must have n_scen values.'
  endif
endelse

; The p-value for significance tests
if not( keyword_set( p_limit ) ) then p_limit = 0.1
; The number of points (~ 1/resolution) to plot in the 1-D PDFs
; (must be odd)
if not( keyword_set( n_hist ) ) then n_hist = 101
; Initialise probability density function output
beta_dist = 1.
dist_weight = 1
; The number of isopleths of probability to sample
if not( keyword_set( dist_prob ) ) then begin
  if n_scen eq 1 then begin
    dist_prob = 1000
  endif else begin
    dist_prob = 100
  endelse
endif
; The number of points to sample along each isopleth of probability.
if not( keyword_set( n_point ) ) then n_point = 2500l
; The number of quantiles to evenly sample along the one dimensional 
; distributions of the regression coefficients
onedim_beta_axis = n_hist

; The partitioning of the noise between the scenarios.
if type ne 'OLS' then scen_noise = 1. / n_scen_samp

; The default number of test samples for comparing distributions
if not( keyword_set( n_test ) ) then n_test = 0

; Option to output to html
html_opt = keyword_set( filename_html )
; Directory of html output
if html_opt eq 1 then begin
  pos = strpos( filename_html, '/', reverse_search=1 )
  if pos eq -1 then begin
    dirname_html = filename_html
  endif else begin
    dirname_html = strmid( filename_html, 0, pos + 1 )
  endelse
endif
; Plot settings
tek_color
!p.multi = 0
charsize = 1.
if html_opt eq 1 then font = 0
if html_opt eq 1 then begin
  thick = 5
endif else begin
  thick = 1
endelse
color_scen = 2 + indgen( n_scen )
color_obs = 1 - html_opt
color_fit = 14
linestyle_anal = 0
linestyle_test = 1
linestyle_real = 2
pos_legend = [ 0.14, 0.6 ]

;***********************************************************************
; Initialise Variables

; Initialise scenarios (data_scen_0 is noise-free, data_scen has noise)
data_scen_0 = fltarr( n_time, n_scen )
; Iterate through scenarios
for i_scen = 0, n_scen - 1 do begin
  ; Create sinusoidal signal of wavenumber (i_scen+1)
  data_scen_0[*,i_scen] = amp_signal[i_scen] $
      * sin( findgen( n_time ) / n_time * 2. * !pi * ( i_scen + 1. ) )
endfor

; Create non-free observations
data_obs_0 = data_scen_0 # beta_real

; Create a noise free version of data_obs which we will use as an 
; "another" scenario
onedim_indep = transpose( data_scen_0 )

; Initialise result arrays for test samples
if n_test gt 0 then begin
  beta_est_test = fltarr( n_scen, 3, n_test )
  resid_sumsq_test = fltarr( n_test )
  p_resid_test = fltarr( n_test )
endif

;***********************************************************************
; Create sample data and estimate regression coefficients

; Iterate through the base sample and further test samples
for i_test = -1, n_test - 1 do begin

  ; If we are using EIV then create uncertain scenarios from model limitation
  if type eq 'EIV' then begin
    ; Set the standard deviation of the variation in the frequency of the 
    ; scenarios
    if n_elements( frac_pert_model ) eq 0 then begin
      frac_pert_model_use = 1. * amp_pert * randomn( seed, n_model )
    endif else begin
      frac_pert_model_use = frac_pert_model
    endelse

    ; Initialise the array containing pure model scenarios
    data_scen_model = fltarr( n_time, n_scen, n_model )
    for i_model = 0, n_model - 1 do begin
      for i_scen = 0, n_scen - 1 do begin
        ;; These differ in wavelength
        ;data_scen_model[*,i_scen,i_model] = amp_signal[i_scen] $
        ;    * sin( findgen( n_time ) / n_time * 2. * !pi $
        ;    * ( ( i_scen + 1. ) * ( 1. + frac_pert_model_use[i_model] ) ) )
        ; These differ in phase
        data_scen_model[*,i_scen,i_model] = amp_signal[i_scen] $
            * sin( findgen( n_time ) / n_time * 2. * !pi $
            * ( i_scen + 1. ) + frac_pert_model_use[i_model] )
      endfor
    endfor
  ; If no model limitation considered then just copy actual scenario data
  endif else begin
    data_scen_model = data_scen_0
  endelse

  ; Add scenario's noise as sum of noise from n_scen_samp samples
  data_scen_use = fltarr( n_time, n_scen, n_model )
  for i_model = 0, n_model - 1 do begin
    for i_scen = 0, n_scen - 1 do begin
      data_scen_use[*,i_scen,i_model] = data_scen_model[*,i_scen,i_model] $
          + amp_epsilon * randomn( seed, n_time ) $
          / sqrt( n_scen_samp[i_scen,i_model] )
      data_scen_use[*,i_scen,i_model] = data_scen_use[*,i_scen,i_model] $
          - mean( data_scen_use[*,i_scen,i_model] )
    endfor
  endfor

  ; Add sampling noise to observations
  data_obs_use = data_obs_0 + amp_epsilon * randomn( seed, n_time )
  data_obs_use = data_obs_use - mean( data_obs_use )

  ; Create control time series for estimating noise
  data_noise_1 = amp_epsilon * randomn( seed, n_time, n_noise_1 )
  for i_noise = 0, n_noise_1 - 1 do begin
    data_noise_1[*,i_noise] = data_noise_1[*,i_noise] $
        - mean( data_noise_1[*,i_noise] )
  endfor
  ; Create more control samples for use in significance testing
  if n_noise_2 gt 0 then begin
    data_noise_2 = amp_epsilon * randomn( seed, n_time, n_noise_2 )
    for i_noise = 0, n_noise_1 - 1 do begin
      data_noise_2[*,i_noise] = data_noise_2[*,i_noise] $
          - mean( data_noise_2[*,i_noise] )
    endfor
  endif else begin
    data_noise_2 = 0
  endelse

  ; If this is the base sample then request distribution outputs
  if i_test eq -1 then begin
    temp_onedim_project_dist = 1
    temp_beta_dist = beta_dist
    temp_dist_weight = dist_weight
    temp_dist_prob = dist_prob
    temp_onedim_beta_axis = onedim_beta_axis
    temp_z_pos = 1
  endif

  ; Run gendetec.pro to estimate the regression coefficients from the data we
  ; have produced.
  ; The output is delivered as beta_est.  beta_est[*,0] are the best-estimate 
  ; values.  The p_limit/2 and (1-p_limit/2) confidence limits are given by 
  ; beta_est[*,1] and beta_est[*,2] respectively. 
  ; Because we expect no covariance between noise estimates, let's not include 
  ; optimisation (no_optimise=1).
  gendetec, data_obs_use, data_scen_use, data_noise_1, temp_beta_est, $
      p_limit=p_limit, frac_noise_var=scen_noise, data_noise_2=data_noise_2, $
      trunc=trunc, type=type, n_point=n_point, beta_dist=temp_beta_dist, $
      dist_prob=temp_dist_prob, dist_weight=temp_dist_weight, $
      onedim_beta_axis=temp_onedim_beta_axis, $
      onedim_beta_dist=temp_onedim_beta_dist, p_resid=temp_p_resid, $
      no_optimise=no_optimise_opt, z_best=temp_z_best, z_pos=temp_z_pos, $
      double=1, onedim_indep=onedim_indep, $
      onedim_project_conf=temp_onedim_project_conf, $
      onedim_project_dist=temp_onedim_project_dist
  ; Estimate root of the sum of the squares of the residuals between data_obs 
  ; and the reconstruction
  temp_resid_sumsq = sqrt( $
      total( ( reform( temp_onedim_project_conf[0,*] ) - data_obs_0 ) ^ 2. ) )

  ; If this is the base sample
  if i_test eq -1 then begin
    ; Record basic results
    beta_est = temp_beta_est
    resid_sumsq = temp_resid_sumsq
    p_resid = temp_p_resid[0]
    ; Record inputs
    data_obs = data_obs_use
    data_scen = data_scen_use
    ; Record distribution results and clear temporary variables
    onedim_beta_dist = temporary( temp_onedim_beta_dist )
    onedim_beta_axis = temporary( temp_onedim_beta_axis )
    beta_dist = temporary( temp_beta_dist )
    dist_prob = temporary( temp_dist_prob )
    dist_weight = temporary( temp_dist_weight )
    onedim_project_dist = temporary( temp_onedim_project_dist )
    onedim_project_conf = temp_onedim_project_conf
    z_best = temporary( temp_z_best )
    z_pos = temporary( temp_z_pos )
  ; If this is a test sample
  endif else begin
    ; Record basic results
    beta_est_test[*,*,i_test] = temp_beta_est
    resid_sumsq_test[i_test] = temp_resid_sumsq
    p_resid_test[i_test] = temp_p_resid[0]
  endelse

endfor

;***********************************************************************
; Compare Results to Truth

; Print output
print, 'The actual regression coefficients are:'
print, beta_real
print, 'The estimated coefficients are:'
print, beta_est[*,0]
print, 'The estimated ' + str( p_limit ) + '-level confidence intervals are:'
print, beta_est[*,1]
print, beta_est[*,2]
if n_test gt 0 then begin
  print, 'The sampled ' + str( p_limit ) + '-level confidence intervals are:'
  id_conf = floor( [ p_limit / 2., 1. - p_limit / 2. ] * n_test )
  temp_conf = fltarr( n_scen, 2 )
  for i_scen = 0, n_scen - 1 do begin
    id_sort = sort( beta_est_test[i_scen,0,*] )
    temp_conf[i_scen,*] = beta_est_test[i_scen,0,id_sort[id_conf]]
  endfor
  print, temp_conf
endif
print, 'The number of samples used to estimate each scenario were:'
print, n_scen_samp
print, 'The p-value the residuals in the regression is:  ' + str( p_resid )

;***********************************************************************
; Plot the PDFs of the First Two Scaling Parameters

; Plot the 2-D multivariate PDF of the first two regression coefficients.  
; Note this may be messy if N_SCEN is large because the N_POINT samples we are 
; taking of the N_SCEN dimensional probability volume may not be a good enough 
; sampling.
if n_scen gt 1 then begin
  ; Open postscript output
  if html_opt eq 1 then ps_open, color=1, filename='temp.ps'
  ; Use a single plotting window
  !p.multi = 0
  ; Estimate the two-dimensional likelihood function
  beta_2dpdf = betas_hist( beta_dist[0:1,*,*], dist_weight, $
      hist_axis=beta_2dpdf_axis, n_hist=n_hist )
  ; Determine a set of plotting levels
  ; (because of numerical issues when estimating quantities based on 
  ; beta_dist, automatic levels can overestimate the useful range)
  levels = max( beta_2dpdf, dimension=1 )
  temp = levels * findgen( n_hist ) / total( levels )
  temp_mean = total( temp )
  temp_std = sqrt( total( ( temp - temp_mean ) ^ 2 ) ) / n_hist
  temp = 0
  id = round( temp_mean + [ -1, 1 ] * temp_std )
  levels = 1.4 * mean( levels[id] )
  n_levels = 20
  levels = ( findgen( n_levels ) + 1. ) / n_levels * levels
  ; Contour plot the result
  contour, beta_2dpdf, beta_2dpdf_axis[*,0], beta_2dpdf_axis[*,1], $
      levels=levels, xtitle='beta1', ytitle='beta2', isotropic=1, $
      title='2-D PDF of first two betas', charsize=charsize, font=font, $
      thick=thick, xthick=thick, ythick=thick, color=color_obs, $
      xrange=beta_2dpdf_axis[[0,n_hist-1],0], $
      yrange=beta_2dpdf_axis[[0,n_hist-1],1]
  ; Overplot test sample results
  if n_test gt 0 then begin
    ; Estimate the smoothed 2-D PDF of samples
    pdf, beta_est_test[0,0,*], beta_est_test[1,0,*], noplot=1, $
        pdf=beta_2dpdf_test, xid=x_beta_2dpdf_test, yid=y_beta_2dpdf_test, $
        xrange=beta_2dpdf_axis[[0,n_hist-1],0], $
        yrange=beta_2dpdf_axis[[0,n_hist-1],1]
    ; Plot the result
    contour, beta_2dpdf_test, x_beta_2dpdf_test, y_beta_2dpdf_test, $
        levels=levels, overplot=1, color=2, thick=thick, $
        c_linestyle=linestyle_test, xrange=beta_2dpdf_axis[[0,n_hist-1],0], $
        yrange=beta_2dpdf_axis[[0,n_hist-1],1]
  endif
  ; Close postscript output and incorporate to html
  if html_opt eq 1 then begin
    ps_close
    spawn, 'ps2epsi temp.ps temp.epsi'
    spawn, 'convert temp.ps -background white -flatten ' + dirname_html $
        + 'beta_pdf_2d.gif'
    spawn, 'rm temp.ps temp.epsi'
  endif
  ; Wait for prompt to continue
  if html_opt eq 0 then begin
    temp_str = ''
    read, 'Press enter to continue...', temp_str
  endif
endif

;***********************************************************************
; Plot one-dimensional likelihood functions of the scaling parameters 
; (regression coefficients)
; Two different ways of doing this are shown.

; Calculate this from the higher level output.
; Initialise PDF array information
beta_1dpdf_high = fltarr( n_hist, n_scen )
beta_1dpdf_high_axis = fltarr( n_hist, n_scen )
; Iterate through scenarios
for i = 0, n_scen - 1 do begin
  ; Estimate the one-dimensional PDF
  temp_1 = 0
  temp = betas_hist( beta_dist[i,*,*], dist_weight, hist_axis=temp_1, $
      n_hist=n_hist )
  beta_1dpdf_high[*,i] = temp
  beta_1dpdf_high_axis[*,i] = temp_1
endfor

; Unfortunately, as N_SCEN gets larger then N_POINT starts not becoming enough 
; to fully sample the N_SCEN dimensional space and we end up with some oddly 
; shaped estimated PDFs.  An alternative is to get the optimal fingerprinting 
; package to do all of these calculations at a lower level.
beta_1dpdf_low_axis = ( onedim_beta_dist[0:n_hist-2,*] $
    + onedim_beta_dist[1:n_hist-1,*] ) / 2.
beta_1dpdf_low = onedim_beta_dist[0:n_hist-2,*]
for i = 0, n_scen - 1 do begin
  beta_1dpdf_low[*,i] = ( onedim_beta_axis[1:n_hist-1] $
      - onedim_beta_axis[0:n_hist-2] ) $
      / ( onedim_beta_dist[1:n_hist-1,i] - onedim_beta_dist[0:n_hist-2,i] )
endfor

; Estimate test sample PDFs for reference
if n_test gt 0 then begin
  beta_1dhist_test = fltarr( n_hist, n_scen )
  beta_1dhist_test_axis = fltarr( n_hist, n_scen )
  for i_scen = 0, n_scen - 1 do begin
    pdf, beta_est_test[i_scen,0,*], noplot=1, pdf=temp_pdf, xid=temp_axis, $
        npdf=n_hist
    beta_1dhist_test[*,i_scen] = temp_pdf $
        / total( temp_pdf * ( temp_axis[n_hist-1] - temp_axis[0] ) / n_hist )
    beta_1dhist_test_axis[*,i_scen] = temporary( temp_axis )
    temp_pdf = 0
    temp_axis = 0
  endfor
endif

; Open postscript output
if html_opt eq 1 then ps_open, color=1, filename='temp.ps'
; Use two plotting windows
!p.multi = [ 0, 1, 2 ]

; Determine the plotting range and axis properties
xrange = [ min( [ beta_1dpdf_high_axis, beta_1dpdf_low_axis ], max=temp ), $
    temp ]
yrange = [ 0, max( [ beta_1dpdf_high, beta_1dpdf_low ] ) ]
xtitle = 'Coefficient value'
ytitle = 'Probability'
; Plot the distributions from the high-level method
title = 'Distributions of coefficients estimated from beta_dist'
plot, beta_1dpdf_high_axis[*,0], beta_1dpdf_high[*,0], nodata=1, $
    yrange=yrange, xrange=xrange, title=title, ytitle=ytitle, xtitle=xtitle, $
    charsize=charsize, font=font, xthick=thick, ythick=thick, color=color_obs
for i_scen = 0, n_scen - 1 do begin
  ; Plot the distribution for this coefficient
  oplot, beta_1dpdf_high_axis[*,i_scen], beta_1dpdf_high[*,i_scen], $
      color=color_scen[i_scen], thick=thick, linestyle=linestyle_anal
  ; Mark the true value
  oplot, beta_real[[i_scen,i_scen]], !y.crange, color=color_scen[i_scen], $
      linestyle=linestyle_real, thick=thick
  ; Plot the test sample distribution
  if n_test gt 0 then begin
    oplot, beta_1dhist_test_axis[*,i_scen], beta_1dhist_test[*,i_scen], $
        color=color_scen[i_scen], linestyle=linestyle_test, thick=thick
  endif
endfor
; Plot the distributions from the low-level method
title = 'Distributions of coefficients estimated from onedim_beta_dist'
plot, beta_1dpdf_low_axis[*,0], beta_1dpdf_low[*,0], nodata=1, yrange=yrange, $
    xrange=xrange, title=title, ytitle=ytitle, xtitle=xtitle, $
    charsize=charsize, font=font, xthick=thick, ythick=thick, color=color_obs
for i_scen = 0, n_scen - 1 do begin
  ; Plot the distribution for this coefficient
  oplot, beta_1dpdf_low_axis[*,i_scen], beta_1dpdf_low[*,i_scen], $
      color=color_scen[i_scen], thick=thick, linestyle=linestyle_anal
  ; Mark the true value
  oplot, beta_real[[i_scen,i_scen]], !y.crange, color=color_scen[i_scen], $
      thick=thick, linestyle=linestyle_real
  ; Plot the test sample distribution
  if n_test gt 0 then begin
    oplot, beta_1dhist_test_axis[*,i_scen], beta_1dhist_test[*,i_scen], $
        color=color_scen[i_scen], linestyle=linestyle_test, thick=thick
  endif
endfor
; Plot line legends
label = 'Estimated distribution'
if n_test gt 0 then label = [ label, 'Distribution from test samples' ]
label = [ label, 'Actual value' ] 
if n_test gt 0 then begin
  linestyle = [ linestyle_anal, linestyle_test, linestyle_real ]
endif else begin
  linestyle = [ linestyle_anal, linestyle_real ]
endelse
line_legend, pos_legend, label, linestyle=linestyle, charsize=2.*charsize, $
    font=font, thick=thick, noborder=1, length=0.04
label = 'Scenario #' + str( 1 + indgen( n_scen ) )
line_legend, pos_legend-[0,1], label, charsize=2.*charsize, $
    color=color_scen, font=font, thick=thick, noborder=1

; Close postscript output and incorporate to html
if html_opt eq 1 then begin
  ps_close
  spawn, 'ps2epsi temp.ps temp.epsi'
  spawn, 'convert temp.ps -background white -flatten ' + dirname_html $
      + 'beta_pdf_1d.gif'
  spawn, 'rm temp.ps temp.epsi'
endif
; Wait for prompt to continue
if html_opt eq 0 then begin
  temp_str = ''
  read, 'Press enter to continue...', temp_str
endif

;***********************************************************************
; Plot residuals

; This can only be done with the test sample
if n_test gt 0 then begin

  ; Open postscript output
  if html_opt eq 1 then ps_open, color=1, filename='temp.ps'
  ; Use two plotting windows
  !p.multi = [ 0, 1, 2 ]

  ; Estimate distribution of residuals from test samples
  temp_axis = 0
  xrange = [ min( resid_sumsq_test, max=temp ), temp ]
  pdf, resid_sumsq_test, noplot=1, pdf=temp_pdf, xid=temp_axis, npdf=n_hist, $
      xrange=xrange
  ; Plot this distribution
  !p.multi = [ 0, 1, 2 ]
  title = 'Distribution of residuals'
  plot, temp_axis, temp_pdf, xrange=xrange, xstyle=1, title=title, $
      linestyle=linestyle_test, thick=thick, charsize=charsize, font=font, $
      xthick=thick, ythick=thick, color=color_obs, yticks=1, ytickname=[' ',' ']
  temp_pdf = 0
  temp_axis = 0
  ; Mark residual from base sample
  oplot, [1,1]*resid_sumsq, !y.crange, linestyle=linestyle_anal, $
      color=color_obs, thick=thick
  ; Mark expected residual
  oplot, [1,1]*amp_epsilon, !y.crange, linestyle=linestyle_real, $
      color=color_obs, thick=thick
  ; Plot legend
  linestyle = [ linestyle_test, linestyle_anal, linestyle_real ]
  label = [ 'Distribution from test samples', 'Value from original sample', $
      'Expected value' ]
  line_legend, pos_legend+[0,0.65], label, linestyle=linestyle, $
      charsize=2.*charsize, font=font, thick=thick, noborder=1, length=0.04

  ; Estimate distribution of p-values of residuals from test samples
  n_p_resid_hist = n_hist
  if n_p_resid_hist gt n_test / 5 then begin
    if n_test mod 5 eq 0 then n_p_resid_hist = n_test / 5
  endif
  p_resid_test_hist_axis = ( findgen( n_p_resid_hist ) + 0.5 ) / n_p_resid_hist
  temp = floor( p_resid_test * n_p_resid_hist )
  p_resid_test_hist = intarr( n_p_resid_hist )
  for i_hist = 0, n_p_resid_hist - 1 do begin
    id = where( temp eq i_hist, n_id )
    p_resid_test_hist[i_hist] = n_id
  endfor
  ; Plot this histogram
  id = [ 0, indgen( n_p_resid_hist ), n_p_resid_hist-1 ]
  title = 'Histogram of p-values of residuals'
  plot, p_resid_test_hist_axis[id], p_resid_test_hist[id], psym=10, $
      xrange=[0,1], xstyle=1, title=title, linestyle=linestyle_test, $
      thick=thick, charsize=charsize, font=font, xthick=thick, ythick=thick, $
      color=color_obs
  ; Mark the p-value of the residual from the base sample
  oplot, [1,1]*p_resid, !y.crange, linestyle=linestyle_anal, color=color_obs, $
      thick=thick
  ; Plot legend
  linestyle = [ linestyle_test, linestyle_anal ]
  label = [ 'Distribution from test samples', 'Value from original sample' ]
  line_legend, pos_legend-[0,1.], label, linestyle=linestyle, $
      charsize=2.*charsize, font=font, thick=thick, noborder=1, length=0.04

  ; Close postscript output and incorporate to html
  if html_opt eq 1 then begin
    ps_close
    spawn, 'ps2epsi temp.ps temp.epsi'
    spawn, 'convert temp.ps -background white -flatten ' + dirname_html $
        + 'residual_pdf.gif'
    spawn, 'rm temp.ps temp.epsi'
  endif
  ; Wait for prompt to continue
  if html_opt eq 0 then begin
    temp_str = ''
    read, 'Press enter to continue...', temp_str
  endif

endif

;***********************************************************************
; Plot the Time Series

; Open postscript output
if html_opt eq 1 then ps_open, color=1, filename='temp.ps'
; Use two plotting windows
!p.multi = [ 0, 1, 2 ]

; Initialise the plot
yrange = [ min( [ [ data_obs ], $
    [ reform( data_scen, n_time, n_scen * n_model ) ], [ z_best ] ], $
    max=temp ), temp ]
plot, data_obs, nodata=1, yrange=yrange, xstyle=1, title='The time series', $
    thick=thick, charsize=charsize, font=font, xthick=thick, ythick=thick
; The data_scen series are in colour
for i_scen = 0, n_scen - 1 do begin
  oplot, data_scen[*,i_scen], color=color_scen[i_scen], thick=thick
endfor
; Plot observed series
oplot, data_obs, thick=thick, color=color_obs
; Plot best estimate of data_obs
oplot, z_best[*,n_scen], color=color_fit, thick=thick
; Plot line legend
label = [ 'Dependent sample', $
    'Scenario #' + str( 1 + indgen( n_scen ) ) + ' sample', $
    'Best regression fit' ]
color = [ color_obs, color_scen, color_fit ]
line_legend, pos_legend+[0,0.65], label, color=color, charsize=2.*charsize, $
    font=font, thick=thick, noborder=1

; Plot the confidence interval on the estimate of data_obs.
; Plot the observed data
plot, data_obs, yrange=yrange, xstyle=1, title='Confidence intervals on fit', $
    thick=thick, charsize=charsize, font=font, xthick=thick, ythick=thick
; Plot the confidence interval
temp = fltarr( n_time, 2 )
for i_time = 0, n_time - 1 do begin
  temp[i_time,0] = min( z_pos[i_time,n_scen,*,0] )
  temp[i_time,1] = max( z_pos[i_time,n_scen,*,0] )
endfor
oplot, temp[*,0], color=color_fit, thick=thick
oplot, temp[*,1], color=color_fit, thick=thick
; Plot the confidence interval from a noise free version of data_obs
oplot, onedim_project_conf[1,*], color=color_fit, linestyle=1, thick=thick
oplot, onedim_project_conf[2,*], color=color_fit, linestyle=1, thick=thick
; Plot line legend
label = [ 'Dependent sample', 'Confidence interval', $
    'Noise-free confidence interval' ]
color = [ color_obs, color_fit, color_fit ]
linestyle = [ 0, 0, 1 ]
line_legend, pos_legend-[0,1.35], label, color=color, charsize=2.*charsize, $
    font=font, thick=thick, noborder=1, linestyle=linestyle

; Close postscript output and incorporate to html
if html_opt eq 1 then begin
  ps_close
  spawn, 'ps2epsi temp.ps temp.epsi'
  spawn, 'convert temp.ps -background white -flatten ' + dirname_html $
      + 'series.gif'
  spawn, 'rm temp.ps temp.epsi'
endif

;***********************************************************************
; Create HTML output page

; Proceed only if requested
if html_opt eq 1 then begin
  ; Open html output file
  openw, 1, filename_html
  ; Write header
  printf, 1, '<HTML>'
  printf, 1, '<HEAD>'
  printf, 1, '<STYLE TYPE="text/css">A { text-decoration: none }</STYLE>'
  printf, 1, '</HEAD>'
  printf, 1, '<BODY BGCOLOR="#bbbb99" LINK="#0033aa" VLINK="#0000aa" ' $
      + 'ALINK="blue" TEXT="#222222" FACE="helvetica" OLOR="#aa0000" SIZE="1">'
  ; Write metadata
  temp_date = strsplit( systime(), ' ', extract=1 )
  temp_date = temp_date[3] + ' on ' + temp_date[0] + ', ' + temp_date[2] + ' ' $
      + temp_date[1] + ' ' + temp_date[4]
  printf, 1, '<B>This output file was produced at ' + temp_date $
      + ' by demo_gendetec.pro.</B><BR><BR>'
  ; Write table of sample sizes
  printf, 1, 'The number of samples used to estimate each scenario<BR>'
  printf, 1, '<TABLE BORDER="1">'
  printf, 1, '  <TR>'
  printf, 1, '    <TD><B>Scenario</B></TD>' $
      + string_from_vector( '<TD><B>' + str( 1 + indgen( n_scen ) ) $
      + '</B></TD>', nospace=1, spacer='' )
  printf, 1, '  </TR>'
  for i_model = 0, n_model - 1 do begin
    printf, 1, '  <TR>'
    printf, 1, '    <TD><B>Model ' + str( 1 + i_model ) + '</B></TD>' $
        + string_from_vector( $
        '<TD>' + str( n_scen_samp[*,i_model] ) + '</TD>', nospace=1, spacer='' )
    printf, 1, '  </TR>'
  endfor
  printf, 1, '</TABLE>'
  printf, 1, '<BR>'
  ; Write p-level of residual
  printf, 1, 'The p-value of the residuals in the regression is:  ' $
       + str( p_resid ) + '<BR><BR>'
  ; Write table of regression coefficients
  printf, 1, 'Actual and estimated regression coefficients<BR>'
  printf, 1, '<TABLE BORDER="1">'
  printf, 1, '  <TR>'
  printf, 1, '    <TD><B>Scenario</B></TD>' $
      + string_from_vector( '<TD><B>' + str( 1 + indgen( n_scen ) ) $
      + '</B></TD>', nospace=1, spacer='' )
  printf, 1, '  </TR>'
  printf, 1, '  <TR>'
  printf, 1, '    <TD><B>Actual</B></TD>' $
      + string_from_vector( '<TD>' + str( beta_real, 2 ) + '</TD>', nospace=1, $
      spacer='' )
  printf, 1, '  </TR>'
  printf, 1, '  <TR>'
  printf, 1, '    <TD><B>Estimated</B></TD>' $
      + string_from_vector( '<TD>' + str( beta_est[*,0], 2 ) + '</TD>', $
      nospace=1, spacer='' )
  printf, 1, '  </TR>'
  printf, 1, '  <TR>'
  printf, 1, '    <TD><B>Estimated confidence interval</B></TD>' $
      + string_from_vector( '<TD>' + str( beta_est[*,1], 2 ) + ', ' $
      + str( beta_est[*,2], 2 ) + '</TD>', nospace=1, spacer='' )
  printf, 1, '  </TR>'
  if n_test gt 0 then begin
    printf, 1, '  <TR>'
    id_conf = floor( [ p_limit / 2., 1. - p_limit / 2. ] * n_test )
    temp_conf = fltarr( n_scen, 2 )
    for i_scen = 0, n_scen - 1 do begin
      id_sort = sort( beta_est_test[i_scen,0,*] )
      temp_conf[i_scen,*] = beta_est_test[i_scen,0,id_sort[id_conf]]
    endfor
    printf, 1, '    <TD><B>Sampled confidence interval</B></TD>' $
      + string_from_vector( '<TD>' + str( temp_conf[*,0], 2 ) + ', ' $
      + str( temp_conf[*,1], 2 ) + '</TD>', nospace=1, spacer='' )
    printf, 1, '  </TR>'
  endif
  printf, 1, '</TABLE>'
  printf, 1, '<BR>'
  ; Plot images
  if n_scen gt 1 then begin
    printf, 1, 'Estimated two-dimensional multivariate PDFs of the first two ' $
        + 'regression coefficients<BR>'
    printf, 1, '<IMG SRC="beta_pdf_2d.gif"><BR><BR>'
  endif
  printf, 1, 'Estimated one-dimensional PDFs of the regression coefficients<BR>'
  printf, 1, '(The legends apply to both plots.)<BR>'
  printf, 1, '<IMG SRC="beta_pdf_1d.gif"><BR><BR>'
  if n_test gt 1 then begin
    printf, 1, 'Residuals from the regressed fit<BR>'
    printf, 1, '<IMG SRC="residual_pdf.gif"><BR><BR>'
  endif
  printf, 1, 'Time series of original data and fits<BR>'
  printf, 1, '<IMG SRC="series.gif"><BR><BR>'
  ; Write end of file
  printf, 1, '</BODY>'
  printf, 1, '</HTML>'
  ; Close html output file
  close, 1
endif

;***********************************************************************
; The End

!p.multi = 0

stop, 'demo_gendetec.pro is finished but has stopped before returning so as ' $
    + 'to allow manual exploration of the results.'

return
END
