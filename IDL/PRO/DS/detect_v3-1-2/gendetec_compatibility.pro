;+
; NAME:
;    GENDETEC_COMPATIBILITY
;
; COPYRIGHT:
;    Copyright (2011) Daithi Stone under contract to the U.S. Department of 
;    Energy's Office of Science, Office of Biological and Environmental 
;    Research and the U.S. National Oceanic and Atmospheric Administration's 
;    Climate Program Office, via the International Detection and Attribution 
;    Group.
;
; PURPOSE:
;    This procedure is a wrapper for gendetec.pro v3.1.0 which allows 
;    simultaneous use of the v3.0.0 and v3.1.0 input/output keyword parameter 
;    names.  
;
; CATEGORY:
;    Optimal Detection Package v3.1.2
;
; CALLING SEQUENCE:
;    gendetec, OBS, SCN, CTL [, BOBS] [, RSSQ]
;
; INPUTS:
;    See gendetec.pro v3.0.0 and v3.1.0.
;
; KEYWORD PARAMETERS:
;    See gendetec.pro v3.0.0 and v3.1.0.
;
; OUTPUTS:
;    See gendetec.pro v3.0.0 and v3.1.0.
;
; USES:
;    gendetec.pro
;
; PROCEDURE:
;    This routine translates only input/output names and then runs gendetec.pro.
;
; EXAMPLE:
;    See demo_gendetec.pro for a demonstration of the use of the optimal 
;    detection routines.
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone (stoned@csag.uct.ac.za), 2011-11-06
;    Modified:  DAS, 2012-03-15 (Fixed bug in call to gendetec.pro;  removed 
;        warning about WEIGHT;  v.3.1.2)
;-

;***********************************************************************

PRO GENDETEC_COMPATIBILITY, $
    OBS, SCN, CTL, BOBS, RSSQ, $
    PV_WAVE=pv_wave, $
    CTLIND=ctlind, WEIGHT=weight, POBS=pobs, COVB=Covb, $
    DOFCTR=dofctr, PRESID=presid, PREWHI=prewhi, $
    PLIMIT=plimit, CINTVL=cintvl, NPOINT=npoint, ICOM=icom, $
    Z_POSS=z_poss, $
    MODVAR=modvar, OBSVAR=obsvar, $
    P_DIST=p_dist, B_DIST=b_dist, W_DIST=w_dist, $
    CTLSVL=ctlsvl, XNOISE=xnoise, D1DIM=d1dim, $
    C1dim=C1dim,ry=ry,rX=rX, $
    B1DIM=b1dim, BAXIS1DIM=baxis1dim, H1DIM=h1dim, $
    ESTVAR=estvar, $
    DATA_NOISE_2=data_noise_2, $
    ADDVAR_OBS=addvar_obs, ADDVAR_SCEN=addvar_scen, $
    BETA_CONFSURF=beta_confsurf, $
    BETA_DIST=beta_dist, DIST_PROB=dist_prob, DIST_WEIGHT=dist_weight, $
    CORR_OBS=corr_obs, CORR_SCEN=corr_scen, $
    COV_NOISE_2=cov_noise_2, $
    DOF_NOISE=dof_noise, $
    FINGER_TRANS=finger_trans, $
    FRAC_NOISE_VAR=frac_noise_var, $
    N_POINT=n_point, $
    NOISE_SINGVAL=noise_singval, $
    ONEDIM_INDEP=onedim_indep, ONEDIM_BETA_AXIS=onedim_beta_axis, $
      ONEDIM_BETA_DIST=onedim_beta_dist, $
      ONEDIM_PROJECT_DIST=onedim_project_dist, $
      ONEDIM_PROJECT_CONF=onedim_project_conf, $
    P_LIMIT=p_limit, $
    P_RESID=p_resid, $
    TRANSFORM=transform, $
    TRUNC=trunc, $
    TYPE=type, OLS=ols_opt, $
    Z_BEST=z_best, $
    DOUBLE=double_opt, $
    NO_OPTIMISE=no_optimise_opt, $
    SCALE_COV_NOISE=scale_cov_noise_opt

;***********************************************************************
; Check for obsolete input

; Stop if PREWHI input
if keyword_set( prewhi ) then begin
  stop, 'gendetec_compatibility.pro:  keyword PREWHI not supported in v3.1.0.'
endif
; Stop if PV_WAVE option set
if keyword_set( pv_wave ) then begin
  stop, 'gendetec_compatibility.pro:  PV_WAVE option not supported in v3.1.0.'
endif

;***********************************************************************
; Translate old names to new, checking for duplication

; Copy BAXIS1DIM input
if keyword_set( baxis1dim ) then begin
  if keyword_set( onedim_beta_axis ) then begin
    stop, $
        'gendetec_compatibility.pro:  input only BAXIS1DIM or ONEDIM_BETA_AXIS.'
  endif else begin
    onedim_beta_axis_use = baxis1dim
  endelse
endif else if keyword_set( onedim_beta_axis ) then begin
  onedim_beta_axis_use = onedim_beta_axis
endif
; Copy CTLIND input
if keyword_set( ctlind ) then begin
  if keyword_set( data_noise_2 ) then begin
    stop, 'gendetec_compatibility.pro:  input only CTLIND or DATA_NOISE_2.'
  endif else begin
    data_noise_2_use = ctlind
  endelse
endif else if keyword_set( data_noise_2 ) then begin
  data_noise_2_use = data_noise_2
endif
; Copy D1DIM input
if keyword_set( d1dim ) then begin
  if keyword_set( onedim_indep ) then begin
    stop, 'gendetec_compatibility.pro:  input only D1DIM or ONEDIM_INDEP.'
  endif else begin
    onedim_indep_use = d1dim
  endelse
endif else if keyword_set( onedim_indep ) then begin
  onedim_indep_use = onedim_indep
endif
; Copy DOFCTR input
if keyword_set( dofctr ) then begin
  if keyword_set( dof_noise ) then begin
    stop, 'gendetec_compatibility.pro:  input only DOFCTR or DOF_NOISE.'
  endif else begin
    dof_noise_use = dofctr
  endelse
endif else if keyword_set( dof_noise ) then begin
  dof_noise_use = dof_noise
endif
; Copy ESTVAR input
if keyword_set( estvar ) then begin
  if keyword_set( scale_cov_noise_opt ) then begin
    stop, 'gendetec_compatibility.pro:  input only ESTVAR or SCALE_COV_NOISE.'
  endif else begin
    scale_cov_noise_opt_use = estvar
  endelse
endif else if keyword_set( scale_cov_noise_opt ) then begin
  scale_cov_noise_opt_use = scale_cov_noise_opt
endif
; Copy ICOM input
if keyword_set( icom ) then begin
  if keyword_set( transform ) then begin
    stop, 'gendetec_compatibility.pro:  input only ICOM or TRANSFORM.'
  endif else begin
    transform_use = icom
  endelse
endif else if keyword_set( transform ) then begin
  transform_use = transform
endif
; Copy MODVAR input
if keyword_set( modvar ) then begin
  if keyword_set( addvar_scen ) then begin
    stop, 'gendetec_compatibility.pro:  input only MODVAR or ADDVAR_SCEN.'
  endif else begin
    addvar_scen_use = modvar
  endelse
endif else if keyword_set( addvar_scen ) then begin
  addvar_scen_use = addvar_scen
endif
; Copy OBSVAR input
if keyword_set( obsvar ) then begin
  if keyword_set( addvar_obs ) then begin
    stop, 'gendetec_compatibility.pro:  input only OBSVAR or ADDVAR_OBS.'
  endif else begin
    addvar_obs_use = obsvar
  endelse
endif else if keyword_set( addvar_obs ) then begin
  addvar_obs_use = addvar_obs
endif
; Copy NPOINT input
if keyword_set( npoint ) then begin
  if keyword_set( n_point ) then begin
    stop, 'gendetec_compatibility.pro:  input only NPOINT or N_POINT.'
  endif else begin
    n_point_use = npoint
  endelse
endif else if keyword_set( n_point ) then begin
  n_point_use = n_point
endif
; Copy P_DIST input
if keyword_set( p_dist ) then begin
  if keyword_set( dist_prob ) then begin
    stop, 'gendetec_compatibility.pro:  input only P_DIST or DIST_PROB.'
  endif else begin
    dist_prob_use = p_dist
  endelse
endif else if keyword_set( dist_prob ) then begin
  dist_prob_use = dist_prob
endif
; Copy PLIMIT input
if keyword_set( plimit ) then begin
  if keyword_set( p_limit ) then begin
    stop, 'gendetec_compatibility.pro:  input only PLIMIT or P_LIMIT.'
  endif else begin
    p_limit_use = plimit
  endelse
endif else if keyword_set( p_limit ) then begin
  p_limit_use = p_limit
endif
; Copy XNOISE input
if keyword_set( xnoise ) then begin
  if keyword_set( frac_noise_var ) then begin
    stop, 'gendetec_compatibility.pro:  input only XNOISE or FRAC_NOISE_VAR.'
  endif else begin
    frac_noise_var_use = xnoise
  endelse
endif else if keyword_set( frac_noise_var ) then begin
  frac_noise_var_use = frac_noise_var
endif

; Copy H1DIM request
if keyword_set( h1dim ) then onedim_project_dist = 1

;***********************************************************************
; Run gendetec.pro

; Call gendetec with hybrid variable names
gendetec, obs, scn, ctl, bobs, rssq, $
    data_noise_2=data_noise_2_use, $
    addvar_obs=addvar_obs_use, addvar_scen=addvar_scen_use, $
    beta_confsurf=beta_confsurf, $
    beta_dist=beta_dist, dist_prob=dist_prob_use, dist_weight=dist_weight, $
    corr_obs=corr_obs, corr_scen=corr_scen, $
    cov_noise_2=cov_noise_2, $
    dof_noise=dof_noise_use, $
    finger_trans=finger_trans, $
    frac_noise_var=frac_noise_var_use, $
    n_point=n_point_use, $
    noise_singval=noise_singval, $
    onedim_indep=onedim_indep_use, onedim_beta_axis=onedim_beta_axis_use, $
      onedim_beta_dist=onedim_beta_dist, $
      onedim_project_dist=onedim_project_dist, $
      onedim_project_conf=onedim_project_conf, $
    p_limit=p_limit_use, $
    p_resid=p_resid, $
    transform=transform_use, $
    trunc=trunc, $
    type=type, ols=ols_opt, $
    weight=weight, $
    z_best=z_best, z_pos=z_poss, $
    double=double_opt, $
    no_optimise=no_optimise_opt, $
    scale_cov_noise=scale_cov_noise_opt_use

;***********************************************************************
; Copy outputs to old parameter names

; Copy BETA_DIST to B_DIST
if keyword_set( beta_dist ) then b_dist = beta_dist
; Copy ONEDIM_BETA_DIST to B1DIM
if keyword_set( onedim_beta_dist ) then b1dim = onedim_beta_dist
; Copy ONEDIM_PROJECT_CONF to C1DIM
if keyword_set( onedim_project_conf ) then c1dim = onedim_project_conf
; Copy BETA_CONFSURF to CINTVL
if keyword_set( beta_confsurf ) then cintvl = beta_confsurf
; Copy COV_NOISE_2 to COVB
if keyword_set( cov_noise_2 ) then covb = cov_noise_2
; Copy NOISE_SINGVAL to CTLSVL
if keyword_set( noise_singval ) then ctlsvl = noise_singval
; Copy ONEDIM_PROJECT_DIST to H1DIM
if keyword_set( onedim_project_dist ) then h1dim = onedim_project_dist
; Copy FINGER_TRANS to POBS
if keyword_set( finger_trans ) then pobs = finger_trans
; Copy P_RESID to PRESID
if keyword_set( p_resid ) then presid = p_resid
; Copy CORR_SCEN to RX
if keyword_set( corr_scen ) then rx = corr_scen
; Copy CORR_OBS to RY
if keyword_set( corr_obs ) then ry = corr_obs
; Copy DIST_WEIGHT to W_DIST
if keyword_set( dist_weight ) then w_dist = dist_weight

;***********************************************************************
; The end

return
END
