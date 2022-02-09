;+
; NAME:
;    GENDETEC
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
;    This procedure computes confidence intervals on undetermined model 
;    parameters by multiple regression using observational constraints.  
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
;    gendetec, DATA_OBS, DATA_SCEN, DATA_NOISE_1 [, BETA_EST] [, RESID_SUMSQ]
;
; INPUTS:
;    DATA_NOISE_1:  An array of dimensions N_DATA by N_NOISE_1 containing the 
;        N_DATA spatio-temporal model data for each of the N_NOISE_1 control 
;        simulations.  This is used for estimating a pre-whitening 
;        transformation.  The spatio-temporal format is identical to DATA_OBS 
;        and DATA_SCEN.
;    DATA_OBS:  A vector of length N_DATA containing the spatio-temporal 
;        observations, with no missing data allowed.
;    DATA_SCEN:  An array of dimensions N_DATA by N_SCEN containing the N_DATA 
;        spatio-temporal model data for each of the N_SCEN scenario 
;        simulations.  The spatio-temporal format is identical to DATA_OBS and 
;        DATA_NOISE_1.
;
; KEYWORD PARAMETERS:
;    ADDVAR_OBS:  An optional floating point input array of size N_DATA*N_DATA 
;        containing the covariance matrix of the additional variance in the 
;        dependent variable beyond that which is expected due to limited 
;        sampling.  Thus element [i,j] contains the covariance between data 
;        points i and j.  For example, if the dependent variable comes from 
;        observations, then element [i,j] will represent the covariance between 
;        space-time points i and j coming from errors in observational 
;        measurements.  Note that the covariance from limited sampling should 
;        not be included in ADDVAR_OBS.
;    ADDVAR_SCEN:  An optional floating point input array of size 
;        N_DATA*N_DATA*N_SCEN containing the inter-estimate covariance 
;        matrices of the additional variance in the scenario patterns beyond 
;        that which is expected due to limited sampling.  Thus element [i,j,k] 
;        contains the covariance between data points i and j for scenario k.  
;        Element [i,j,k] represents the covariance across models between 
;        space-time points i and j in the estimation of pattern k.  Note that 
;        the covariance from limited sampling should not be included in 
;        ADDVAR_SCEN.  If ADDVAR_SCEN is not input but EIV is requested then 
;        ADDVAR_SCEN is calculated from DATA_SCEN.
;    BETA_CONFSURF:  Returns an N_TRANSFORM*N_POINT*N_SCEN array containing the 
;        coordinates of the isopleth of probability surface of the estimates of 
;        the regression coefficients in various dimensions.  N_SCEN is the 
;        number of scenarios and N_POINT is the sampling size on the isopleth 
;        surface defined in P_LIMIT (two sided).  Element [I,J,K] gives the 
;        coordinate along the direction of regression coefficient I of the 
;        sampling point J on the (K+1) dimensional confidence interval 
;        surface.  If BETA_CONFSURF[*,J,K]=0 then the confidence interval is 
;        unbounded (possible with TLS/EIV).
;    BETA_DIST:  If set then this returns an array containing the location 
;        (beta values) of the probability density estimates of the beta scaling 
;        parameters.  Of size [N_TRANSFORM,N_POINT,N_DIST], where N_TRANSFORM 
;        is the number of linear combinations of scenarios requested in 
;        TRANSFORM if given or otherwise N_TRANSFORM=N_SCEN where N_SCEN is the 
;        number of input scenarios, N_POINT is the number of points sampled 
;        along an isopleth of probability, and N_DIST is the number of 
;        isopleths sampled.  So elements [*,j,k] are the coordinates of the jth 
;        point on the kth isopleth of probability.  Initialise to 1 to ensure 
;        output, also initialise DIST_WEIGHT and input DIST_PROB.  See N_POINT, 
;        DIST_PROB, and DIST_WEIGHT for more information.
;    CORR_OBS:  Returns the correlation coefficient between DATA_OBS and the 
;        projection of DATA_OBS on the control simulation EOFs.
;    CORR_SCEN:  Returns the correlation coefficients between the signal 
;        patterns and the projection of the signal patterns onto the control 
;        simulation EOFs.  A vector of length N_SCEN (the number of scenarios).
;    COV_NOISE_2:  If set then returns an N_TRANSFORM*N_TRANSFORM matrix 
;        containing the estimated covariance of the noise realisations in 
;        DATA_NOISE_2 in the directions of the N_TRANSFORM scenario patterns in 
;        from the linearly transformed DATA_SCEN (or else the N_SCEN patterns 
;        from DATA_SCEN).  Approximate if TLS is used.  Currently returned 
;        always if OLS is used.
;    DATA_NOISE_2:  An optional array of dimensions N_DATA by N_NOISE_2 
;        containing N_DATA spatio-temporal model data for each of N_NOISE_2 
;        unforced control simulations.  This is for hypothesis testing, 
;        independent of the original noise estimate.  The spatio-temporal 
;        format is identical to DATA_OBS, DATA_NOISE_1, and DATA_SCEN.
;    DIST_PROB:  A vector of length N_DIST containing the isopleths of 
;        probability to sample when estimating the probability density of the 
;        scaling parameters.  If a scalar value is entered then that number of 
;        equally spaced isopleths will be sampled over the (0,1) interval.  See 
;        BETA_DIST, N_POINT, and DIST_WEIGHT for more information.
;    DIST_WEIGHT:  Returns an array containing the weighting estimates of the 
;        beta scaling parameters at the locations in BETA_DIST.  Note these are 
;        not actual density estimates because BETA_DIST contains an irregular 
;        sampling (polar).  Of size [N_POINT,N_DIST], so elements [j,k] are the 
;        fraction of the PDF at the jth point on the kth isopleth of 
;        probability.  Initialise to 1 to ensure output, also initialise 
;        BETA_DIST and input DIST_PROB.  See BETA_DIST, N_POINT, and DIST_PROB 
;        for more information.
;    DOF_NOISE:  The number of degrees of freedom in the estimate of the noise 
;        covariance to be used when estimating probability distributions.  The 
;        default is N_NOISE_2, the number of independent noise realisations in 
;        DATA_NOISE_2.
;    DOUBLE:  If set then calculations are done in double precision arithmetic. 
;        The default is single precision.
;    FINGER_TRANS:  If set, then this returns the transpose of the matrix 
;        of fingerprint patterns.  Note that 
;        BETA_EST[*,0]=FINGER_TRANS#DATA_OBS, to within the EOF space covered 
;        to the truncation TRUNC.  Of size N_SCEN*N_DATA where N_SCEN is the 
;        number of scenarios in DATA_SCEN and N_DATA is the length of 
;        DATA_OBS.
;    FRAC_NOISE_VAR:  A floating point array of length N_SCEN.  It contains the 
;        fraction of the variance of the noise in DATA_OBS contained in each of 
;        the N_SCEN variables in DATA_SCEN.  Required for the TLS method, so if 
;        not defined then the procedure defaults to the OLS method.
;    N_POINT:  The number of points N_POINT to sample along an isopleth of 
;        probability when estimating the multivariate probability density 
;        function of the scaling factors.  An isopleth of probability is an 
;        N_SCEN dimensional ellipsoid and is sampled by the N_POINT points in a 
;        polar coordinate system.  See BETA_DIST, DIST_PROB, and DIST_WEIGHT 
;        for more information.
;    NO_OPTIMISE:  If set, the procedure does not optimise the fingerprints, 
;        but it still projects onto the noise EOFs procedurally.  The default 
;        is to optimise the fingerprints.
;    NOISE_SINGVAL:  Returns a vector of length N_NOISE_1 containing the 
;        singular values of the input matrix DATA_NOISE_1 (equivalently the 
;        square roots of the eigenvalues of the covariance matrix of 
;        DATA_NOISE_1).
;    OLS:  If set the procedure uses the ordinary least squares (OLS) 
;        estimators.  The default is to use the total least squares (TLS) 
;        algorithm.  This keyword is superceded by the TYPE keyword and 
;        overruled by it.
;    ONEDIM_BETA_AXIS:  A vector of size N_ONEDIM_BETA containing the quantiles 
;        at which to sample the one dimensional probability distributions of 
;        the betas (the regression coefficients).  Values must be in the (0,1) 
;        range.
;    ONEDIM_BETA_DIST:  Returns the locations of the quantiles of the one 
;        dimensional probability distributions of the betas (the regression 
;        coefficients).  The quantiles to sample are defined in 
;        ONEDIM_BETA_AXIS.  Returns an array of size [N_ONEDIM_BETA,N_SCEN] 
;        where N_SCEN is the number of scenarios and N_ONEDIM_BETA is the size 
;        of ONEDIM_BETA_AXIS.  Input from ONEDIM_BETA_AXIS is required.
;    ONEDIM_INDEP:  An array of scenarios (as in DATA_SCEN) for which to return 
;        (in ONEDIM_PROJECT_CONF and optionally ONEDIM_PROJECT_DIST) the one 
;        dimensional confidence intervals of the attributable component of 
;        DATA_OBS scaled by the scalings in BETA_EST.  Of size N_SCEN*N_ONEDIM 
;        where N_SCEN is the number of scenarios and N_ONEDIM is the number of 
;        conditions of scenario combinations to consider.  Values can be within 
;        or outside of the range of values in DATA_SCEN.
;    ONEDIM_PROJECT_CONF:  Returns an array of size 3*N_ONEDIM containing the 
;        best ([0,*]), lower confidence range ([1,*]), and upper confidence 
;        range ([2,*]) estimates of the one dimensional attributable components 
;        given the hypothetical scenarios in ONEDIM_INDEP.  N_ONEDIM is the 
;        number of scenario combinations to examine.
;    ONEDIM_PROJECT_DIST:  If set, then this returns an array describing the 
;        likelihood functions of the one dimensional attributable component 
;        given the parameter combinations in ONEDIM_INDEP.  Of size 
;        N_ONEDIM_BETA*ONEDIM_INDEP.  The attributable values of quantile 
;        ONEDIM_BETA_AXIS[i] are given in ONEDIM_PROJECT_DIST[i,*].  If 
;        ONEDIM_BETA_AXIS equals [0.5,P_LIMIT/2,1-P_LIMIT/2] then 
;        ONEDIM_PROJECT_DIST equals ONEDIM_PROJECT_CONF.
;    P_LIMIT:  An optional floating point number containing the
;        two-sided p-value for confidence interval estimation.  The default is 
;        0.1.
;    P_RESID:  Returns the P-value of the F or chi-squared test on the 
;        residuals between DATA_OBS and the best construction of DATA_OBS 
;        obtained from the regression.  If TRUNC is set to a positive integer 
;        then P_RESID is a vector of length 1 containing the value 
;        corresponding to the truncation TRUNC;  otherwise P_RESID is a vector 
;        of length NTRUNC+1 where P_RESID[0] is the value corresponding to the 
;        maximum allowed truncation, while P_RESID[1:NTRUNC] are the values 
;        corresponding to the N_TRUNC truncations defined in TRUNC.
;    SCALE_COV_NOISE:  If set, then the covariance estimates are scaled by 
;        dividing by the residual sum of squares (RESID_SUMSQ).
;    SCAN_BETA_EST:  If a range of TRUNC values are scanned before a final 
;        value is selected and used, then this returns an N_TRANSFORM*3*N_TRUNC 
;        array containing the BETA_EST values estimated for all of the scanned 
;        truncations.
;    SCAN_P_RESID:  If a range of TRUNC values are scanned before a final value 
;        is selected and used, then this returns an N_TRUNC vector containing 
;        the P_RESID values estimated for all of the scanned truncations.
;    SCAN_RESID_SUMSQ:  If a range of TRUNC values are scanned before a final 
;        value is selected and used, then this returns an N_TRUNC vector 
;        containing the RESID_SUMSQ values estimated for all of the scanned 
;        truncations.
;    TRANSFORM:  An optional string vector of length N_TRANSFORM containing 
;        linear transformation instructions for the coefficients on output.  
;        The format is of the form '{sign}{index}' or 
;        '{scaling factor}*{index}'.  {sign} can be '+', '-', or '', 
;        {scaling factor} can be any integer or floating point number, and 
;        {index} is the index number of the requested independent variable.  
;        Sequences of these command formats are allowed too.  So 
;        ['0.5*0','-1+0.5*0'] would return coefficients for half of the first 
;        variable and the difference of half the first variable and all the 
;        second variable.  Output from other variables is similarly rotated, 
;        except ONEDIM_PROJECT_CONF and ONEDIM_BETA_DIST.  
;        The old integer input format is still supported.
;        If TRANSFORM is not given, then for output purposes, N_TRANSFORM can 
;        be considered equal to N_SCEN and TRANSFORM can be considered to have 
;        been an identity matrix.
;    TRUNC:  If set to a positive integer, then this is the number of EOFs of 
;        the control simulation to retain in the noise model, and N_TRUNC=1.  
;        If TRUNC is not set then N_TRUNC=N_BLEACH-N_SCEN, where N_BLEACH is 
;        the degrees of freedom of DATA_NOISE_1, and the routine will scan 
;        through all truncations from N_SCEN+1 to N_TRUNC to find the maximum 
;        truncation for which P_RESID>P_LIMIT, and then use that truncation in 
;        the calculations.  If TRUNC is set to a negative value, then 
;        N_TRUNC=abs(TRUNC) and the maximum truncation in the range 
;        [N_SCEN+1,N_TRUNC] for which P_RESID>P_LIMIT is calculated and used.  
;        If TRUNC set to an integer vector of length N_TRUNC, then only 
;        suggested values are scanned and the maximum for which 
;        P_RESID>P_LIMIT is used.
;        TRUNC returns an integer containing the truncation value that was used.
;    TYPE:  An optional string describing the type of regression algorithm to 
;        use.  Possibilities are:
;          'OLS' for ordinary least squares (OLS),
;          'TLS' for total least squares (TLS),
;          'EIV' for error in variables (EIV).
;    WEIGHT:  An optional floating point vector of length N_DATA defining 
;        weights on the N_DATA elements in the input data arrays.
;    Z_BEST:  Returns an array containing the best estimate values of noise 
;        free (linearly transformed) scenarios and observations.  The array is 
;        of size [N_DATA,N_TRANSFORM+1].  Elements [*,0:N_TRANSFORM-1] contain 
;        the estimates for the (transformed) scenarios, while elements 
;        [*,N_TRANSFORM] contains the estimates for the observations.  While 
;        the estimations are done in truncated EOF space, the results are 
;        projected back into normal space in Z_BEST.  But because we only use 
;        the leading TRUNC EOFs, the estimates will look progressively worse 
;        the smaller TRUNC is compared to N_DATA.
;    Z_DIST:  Returns an array containing the estimated values of the estimated 
;        noise-free DATA_SCEN after linear transformation according to 
;        TRANSFORM (if given) and DATA_OBS on the isopleths of probability 
;        defined in DIST_PROB.  Of size N_DATA*(N_TRANSFORM+1)*N_POINT*N_DIST.  
;        Elements [*,0:N_TRANSFORM-1,j,k] correspond to the estimates of 
;        BETA_DIST[*,j,k] for the (transformed) independent variables from 
;        DATA_SCEN, while elements [*,N_TRANSFORM,j,k] are for the dependent 
;        variable in DATA_OBS.
;    Z_POS:  If set then returns an array containing the estimated locations 
;        (values) of the probability density estimates of the values of the 
;        estimated noise-free (linearly transformed) DATA_SCEN and DATA_OBS.  
;        Of size [N_DATA,N_TRANSFORM+1,N_POINT,N_SCEN] where N_DATA is the 
;        number of values in DATA_OBS, N_TRANSFORM is the number of linear 
;        transformed scenarios from TRANSFORM and DATA_SCEN (or is just 
;        N_SCEN), and N_POINT is the number of points to use to sample the 
;        probability distribution.  Elements [*,0:N_TRANSFORM-1,i,k] correspond 
;        to the estimated values of (the linearly transformed) DATA_SCEN 
;        corresponding to BETA_CONFSURF[*,j,k], that is the jth sampled point 
;        on the (k+1)-dimensional confidence interval surface.  Elements 
;        [*,N_TRANSFORM,j,k] correspond to the estimated values of DATA_OBS 
;        corresponding to BETA_CONFSURF[*,j,k].  Also see BETA_CONFSURF and 
;        N_POINT.
;
; OUTPUTS:
;    BETA_EST:  An array of size N_SCEN*3 containing the best ([*,0]), lower 
;        confidence range ([*,1]), and upper confidence range ([*,2]) estimates 
;        of the projection amplitudes of the observations DATA_OBS onto the 
;        various scenarios DATA_SCEN.
;    RESID_SUMSQ:  Returns the sum of squares of the normalised residual 
;        differences between the prewhitened DATA_OBS and the best regression 
;        estimate.  The normalisation is done for each element of DATA_OBS 
;        according to the standard deviation of the corresponding N_NOISE_2 
;        realisations of that element in DATA_NOISE_2.  An adjustment for the 
;        bias in RESID_SUMSQ is made if OLS is used in order to account for 
;        bias.  If scanning over multiple truncations (see TRUNC) then this 
;        returns a vector with RESID_SUMSQ[0] containing the value when using 
;        the maximum allowed truncation and RESID_SUMSQ[N_SCEN+1:abs(TRUNC)] 
;        containing the values for truncations of N_SCEN+1 to TRUNC.
;    BETA_CONFSURF, BETA_DIST, CORR_OBS, CORR_SCEN, COV_NOISE_2, DIST_WEIGHT, 
;      FINGER_TRANS, NOISE_SINGVAL, ONEDIM_BETA_DIST, ONEDIM_PROJECT_CONF, 
;      ONEDIM_PROJECT_DIST, P_RESID, SCAN_BETA_EST, SCAN_P_RESID, 
;      SCAN_RESID_SUMSQ, TRUNC, Z_BEST, Z_POS
;
; USES:
;    algebra_to_matrix.pro
;    linmod.pro
;    pca.pro
;    regtls.pro
;
; PROCEDURE:
;    This routine is primarily a driver for linmod.pro, but computes 
;    prewhitening transformations and transformations on input/output variables 
;    if uncertainties on parameter-combinations are required, and scans for 
;    maximum allowable truncation if not specified.
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
;    Modified:  MRA, 1999-08-09 (Z_POSS keyword implemented; v1.3)
;    Modified:  MRA, 1999-08-13 (Name revised for DOS file transfers;  v1.4)
;    Modified:  MRA, 1999-01-11 (Bug in use of weights fixed; v1.5)
;    Modified:  MRA, 2000-08-03 (Explicit PDFs computed; v2.0)
;    Modified:  Daithi Stone (stoned@atm.ox.ac.uk), 2004-06-28 
;        (Documentation for inclusion in routine library)
;    Modified:  DAS, 2005-03-13 (updated documentation, compliance with 
;        svdpvw.pro and linmod.pro)
;    Modified:  DAS, 2005-05-13 (Updated documentation)
;    Modified:  DAS, 2005-09-01 (Added B1DIM, BAXIS1DIM, DOUBLE keywords and 
;        related calculations;  updated documentation;  v3.0.0)
;    Modified:  DAS, 2005-11-01 (Fixed bug which insisted on B_DIST input)
;    Modified:  DAS, 2005-12-16 (Fixed bug which did not allow D1DIM1 to be 
;        scalar;  added H1DIM keyword)
;    Modified:  DAS, 2008-04-01 (Added string input capability to ICOM;  
;        added capability for ICOM to be of different size than N_SCN)
;    Modified:  DAS, 2008-07-09 (Noted apparent bug in presid calculation 
;        when CTLIND is not given)
;    Modified:  DAS, 2010-02-11 (Added H1DIM keyword;  edited documentation 
;        formatting;  added PV-WAVE notes)
;    Modified:  DAS, 2011-01-03 (Fixed bug in returning N_ICOM values of BOBS)
;    Modified:  DAS, 2011-11-08 (Discontinued PV-WAVE compatibility;  
;        outsourced TRANSFORM parsing to algebra_to_matrix.pro;  removed 
;        pre-v3.0 header documentation;  changed default p_limit to 0.1;  
;        removed use of fcdf.pro;  switched use of svdpvw.pro to pca.pro;  
;        altered code variable names;  changed input and output variable names 
;        (see translation notes below);  removed PREWHI, STATUS, and WEIGHT 
;        keywords;  v3.1.0)
;    Modified:  DAS, 2012-01-23 (Switched to SVD option with pca.pro;  blocked 
;        use of EIV approach;  v3.1.1)
;    Modified:  DAS, 2012-03-15 (Re-introduced WEIGHT keyword; introduced 
;        SCAN_BETA_EST, SCAN_P_RESID, and SCAN_RESID_SUMSQ keywords;  
;        introduced option of vector input for TRUNC;  added further checks on 
;        TRUNC input;  corrected when N_TRANSFORM should be used in the 
;        documentation instead of N_SCEN;  v3.1.2)
;-

;***********************************************************************

PRO GENDETEC, $
    DATA_OBS, DATA_SCEN, DATA_NOISE_1, BETA_EST, RESID_SUMSQ, $
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
    TRUNC=trunc, SCAN_BETA_EST=scan_beta_est, SCAN_P_RESID=scan_p_resid, $
      SCAN_RESID_SUMSQ=scan_resid_sumsq, $
    TYPE=type, OLS=ols_opt, $
    WEIGHT=weight, $
    Z_BEST=z_best, Z_POS=z_pos, $
    DOUBLE=double_opt, $
    NO_OPTIMISE=no_optimise_opt, $
    SCALE_COV_NOISE=scale_cov_noise_opt

; Changes in input/output variable names from v3.0.0 to v3.1.0
;   changed b_dist to beta_dist
;   changed b1dim to onedim_beta_dist
;   changed baxis1dim to onedim_beta_axis
;   changed bobs to beta_est
;   changed c1dim to onedim_project_conf
;   changed cintvl to beta_confsurf
;   changed covb to cov_noise_2
;   changed ctl to data_noise_1
;   changed ctlind to data_noise_2
;   changed ctlsvl to noise_singval
;   changed d1dim to onedim_indep
;   changed dofctr to dof_noise
;   changed estvar to scale_cov_noise
;   changed h1dim to onedim_project_dist
;   changed icom to transform
;   changed modvar to addvar_scen
;   changed no_opt to no_optimise
;   changed npoint to n_point
;   changed obs to data_obs
;   changed obsvar to addvar_obs
;   changed p_dist to dist_prob
;   changed pobs to finger_trans
;   changed plimit to p_limit
;   changed presid to p_resid
;   changed rssq to resid_sumsq
;   changed rx to corr_scen
;   changed ry to corr_obs
;   changed scn to data_scn
;   changed w_dist to dist_weight
;   changed xnoise to frac_noise_var
;   changed z_poss to z_pos

;***********************************************************************
; Options

; Determine whether we are in GDL
gdl_opt = !prompt eq 'GDL> '

; Option for no optimisation
no_optimise_opt = keyword_set( no_optimise_opt )

; Option for double precision
double_opt = keyword_set( double_opt )
one = 1.
if double_opt eq 1 then one = double( one )

; Default regression type
if not( keyword_set( type ) ) then begin
  ; Default to old ols_opt request
  if keyword_set( ols_opt ) then begin
    type = 'OLS'
  ; Otherwise check if frac_noise_var is given
  endif else if keyword_set( frac_noise_var ) then begin
    ; Check if addvar_scen is given
    if keyword_set( addvar_scen ) then begin
      print, 'Assuming EIV method desired.'
      type = 'EIV'
    endif else begin
      print, 'Assuming TLS method desired.'
      type = 'TLS'
    endelse
  endif else begin
    print, 'Assuming OLS method desired.'
    type = 'OLS'
  endelse
endif
type = strupcase( type )

; Block use of EIV approach
if type eq 'EIV' then begin
  stop, 'gendetec.pro:  ' $
      + 'Use of the EIV approach has been disabled in this version.'
endif

;***********************************************************************
; Set constants

; The default two-sided p-level for significance testing
if not( keyword_set( p_limit ) ) then p_limit = 0.1

; Determine sizes of inputs.
; Number of points in observed sample
n_data = n_elements( data_obs[*,0] )
; Number of signals/scenarios
n_scen = n_elements( data_scen[0,*,0] )
; Number of models
n_model = n_elements( data_scen[0,0,*] )
; Number of samples of noise for pre-whitening estimation
n_noise_1 = n_elements( data_noise_1[0,*] )
; Number of samples of noise for hypothesis-testing
n_noise_2 = n_elements( data_noise_2 ) / n_data

; Count number of linear transformation instructions
n_transform = n_elements( transform )
; Create linear transformation matrix according to instructions
if n_transform ne 0 then begin
  transform_matrix = one * algebra_to_matrix( transform, n_var=n_scen )
  transform_matrix = reform( transform_matrix, n_scen, n_transform )
; Otherwise instruct for direct output
endif else begin
  n_transform = n_scen
  transform_matrix = identity( n_transform )
endelse

; Check number of scenario values for processing
n_onedim = n_elements( onedim_indep )
if n_onedim ne 0 then begin
  n_onedim = n_elements( onedim_indep[0,*] )
  if n_elements( onedim_indep[*,0] ) ne n_scen then begin
    stop, 'gendetec.pro:   ' $
        + 'The first dimension of ONEDIM_INDEP must have N_SCEN elements.'
  endif
endif
; Set last n_scen values equal to transform_matrix so that it all gets 
; processed together
if n_onedim ne 0 then begin
  onedim_indep_out = transpose( [ transpose( [onedim_indep] ), $
      transpose( transform_matrix ) ] )
endif else begin
  onedim_indep_out = transform_matrix
endelse

; The default degrees of freedom for significance tests
if not( keyword_set( dof_noise ) ) then begin
  ; The rank of the second input noise variable if provided
  if n_noise_2 gt 0 then begin
    dof_noise = n_noise_2
  ; Otherwise we will substitute the first input noise variable later
  endif else begin
    dof_noise = n_noise_1
  endelse
endif

; The default isopleths of probability to sample (defined by density)
if keyword_set( beta_dist ) then begin
  ; Count the requested number of input isopleth sampling requests
  n_dist = n_elements( dist_prob )
  ; Default dist_prob as density of 1./100.
  if n_dist eq 0 then dist_prob = 100
  ; If dist_prob is supplied as a density
  if n_dist eq 1 then begin
    n_dist = dist_prob[0]
    dist_prob = ( 0.5 + one * findgen( n_dist ) ) / n_dist
  ; If specific dist_prob values are supplied
  endif else begin
    ; dist_prob must be in ascending order
    id = sort( dist_prob )
    dist_prob = one * dist_prob[id]
    ; Check for acceptable values
    if ( dist_prob[0] le 0. ) or ( dist_prob[n_dist-1] ge 1. ) then begin
      stop, 'gendetec.pro:  Impossible p-values:', dist_prob[0], $
          dist_prob[n_dist-1]
    endif
  endelse
  ; Compute critical values of t-statistic corresponding to dist_prob
  t_dist = one * fltarr( n_dist )
  for i_dist = 0, n_dist - 1 do begin
    t_dist[i_dist] $
        = sqrt( n_scen * f_cvf( one * dist_prob[i_dist], n_scen, dof_noise ) )
  endfor
endif

; Compute critical values of the t- or T-statistic corresponding to the 
; values in onedim_beta_axis
if keyword_set( onedim_beta_axis ) then begin
  ; If only the size of the axis vector is given then create the axis vector
  if n_elements( onedim_beta_axis ) eq 1 then begin
    onedim_beta_axis = findgen( onedim_beta_axis ) / onedim_beta_axis $
        + 1. / onedim_beta_axis / 2.
  endif
  ; Initialise the vector of critical values
  onedim_beta_t = one * fltarr( n_elements( onedim_beta_axis ) )
  ; Iterate through values in onedim_beta_axis
  for i = 0, n_elements( onedim_beta_axis ) - 1 do begin
    if onedim_beta_axis[i] gt 0.5 then begin
      onedim_beta_t[i] = -sqrt( f_cvf( 2.*one*(1.-onedim_beta_axis[i]), 1, $
          dof_noise ) )
    endif else begin
      onedim_beta_t[i] = sqrt( f_cvf( 2.*one*onedim_beta_axis[i], 1, $
          dof_noise ) )
    endelse
  endfor
endif

; Compute critical values of the t- or T-statistic corresponding to p_limit
t_crit = one * fltarr( n_scen )
for i_scen = 1, n_scen do begin
  t_crit[i_scen-1] = sqrt( i_scen * f_cvf( one * p_limit, i_scen, dof_noise ) )
endfor

;***********************************************************************
; Evaluate inputs

; Ensure supported regression type
if max( type eq [ 'OLS', 'TLS', 'EIV' ] ) ne 1 then begin
  stop, 'gendetec.pro:  Unsupported regression type.'
endif

; Ensure data_obs is 1-dimensional
if ( size( data_obs ) )[0] ne 1 then begin
  stop, 'gendetec.pro:  data_obs must be one-dimensional..'
endif

; Ensure data_scen is 2- or 3-dimensional and of correct size
dim_scen = ( size( data_scen ) )[0]
if ( max( type eq ['OLS','TLS'] ) eq 1 ) and ( dim_scen ne 2 ) then begin
  ; Resize to two-dimensions if one-dimensional
  if dim_scen lt 2 then begin
    data_scen = reform( data_scen, n_elements( data_scen[*,0] ), n_scen )
  ; Otherwise report error
  endif else begin
    stop, 'gendetec.pro:  data_scen must be two-dimensional for OLS/TLS.'
  endelse
endif else if ( type eq 'EIV' ) and ( dim_scen ne 3 ) then begin
  ; Resize to three-dimensions if one- or two-dimensional
  if dim_scen lt 3 then begin
    data_scen = reform( data_scen, n_elements( data_scen[*,0,0] ), n_scen, $
        n_model )
  ; Otherwise report error
  endif else begin
    stop, 'gendetec.pro:  data_scen must be three-dimensional for EIV.'
  endelse
endif
if n_elements( data_scen[*,0,0] ) ne n_data then begin
  stop, 'gendetec.pro:  data_obs and data_scen have incompatible sizes.'
endif

; Ensure data_noise_1 is 2-dimensional and of correct size
if ( size( data_noise_1 ) )[0] ne 2 then begin
  stop, 'gendetec.pro:  data_noise_1 must be two-dimensional.'
endif
if n_elements( data_noise_1[*,0] ) ne n_data then begin
  stop, 'gendetec.pro:  data_obs and data_noise_1 have incompatible sizes.'
endif
; Ensure data_noise_2 is 2-dimensional and of correct size
if n_noise_2 gt 0 then begin
  if ( size( data_noise_2 ) )[0] ne 2 then begin
    stop, 'gendetec.pro:  data_noise_2 must be two-dimensional.'
  endif
  if n_elements( data_noise_2[*,0] ) ne n_data then begin
    stop, 'gendetec.pro:  data_obs and data_noise_2 have incompatible sizes.'
  endif
endif

; Ensure weight is 1-dimensional and of correct size
if keyword_set( weight ) then begin
  if ( size( weight ) )[0] ne 1 then begin
    stop, 'gendetec.pro:  weight must be one-dimensional.'
  endif
  if n_elements( weight[*,0] ) ne n_data then begin
    stop, 'gendetec.pro:  data_obs and weight have incompatible sizes.'
  endif
  ; Also ensure that no_optimise is set
  if no_optimise_opt eq 0 then begin
    stop, 'gendetec.pro:  ' $
        + 'weight input makes no difference in no_optimise is not set.'
  endif
endif

;***********************************************************************
; Set up variables

; Copy inputs to working arrays
data_obs_use = one * data_obs
data_scen_use = one * data_scen
data_noise_1_use = one * data_noise_1
; Use independent estimate of noise for significance estimation if input
if n_noise_2 gt 0 then begin
  data_noise_2_use = one * data_noise_2
; Otherwise use the original estimate
endif else begin
  data_noise_2_use = data_noise_1_use
  n_noise_2 = n_noise_1
endelse
; Copy vector of noise contribution to signals
if keyword_set( frac_noise_var ) then frac_noise_var_use = frac_noise_var

; Copy weight and ensure it sums to unity
if keyword_set( weight ) then begin
  weight_use = one * weight / total( weight )
; Otherwise take even weighting
endif else begin
  weight_use = one / n_data + 0. * data_obs
endelse

; Take anomalies from input data arrays
data_obs_use = data_obs_use - mean( data_obs_use )
for i_model = 0, n_model - 1 do begin
  for i_scen = 0, n_scen - 1 do begin
    data_scen_use[*,i_scen,i_model] = data_scen_use[*,i_scen,i_model] $
        - total( weight_use * data_scen_use[*,i_scen,i_model] )
  endfor
endfor
for i_noise = 0, n_noise_1 - 1 do begin
  data_noise_1_use[*,i_noise] = data_noise_1_use[*,i_noise] $
      - total( weight_use * data_noise_1_use[*,i_noise] )
endfor
for i_noise = 0, n_noise_2 - 1 do begin
  data_noise_2_use[*,i_noise] = data_noise_2_use[*,i_noise] $
      - total( weight_use * data_noise_2_use[*,i_noise] )
endfor

; Estimate inter-estimate covariance matrix for EIV.
; Do this only for EIV
if type eq 'EIV' then begin
  ; Weighting not yet supported for EIV
  if keyword_set( weight ) then begin
    stop, 'gendetec.pro:  weight not yet supported with EIV.'
  endif
  ; Infer sample sizes from frac_noise_var
  n_sample_scen = 1. / frac_noise_var
  ; Estimate model-mean scenario
  data_scen_use_mean = one * fltarr( n_data, n_scen )
  for i_scen = 0, n_scen - 1 do begin
    data_scen_use_mean[*,i_scen] = data_scen_use_mean[*,i_scen] $
        + ( reform( data_scen_use[*,i_scen,*] ) $
        # ( reform( n_sample_scen[i_scen,*] ) $
        / total( n_sample_scen[i_scen,*] ) ) )
  endfor
  ; If addvar_scen is not already set
  if not( keyword_set( addvar_scen ) ) then begin
    ; Estimate TLS fits of model scenario responses to model-mean responses
    fit_est_model = one * fltarr( n_data, n_scen, n_model )
    for i_model = 0, n_model - 1 do begin
      for i_scen = 0, n_scen - 1 do begin
        temp_beta = regtls( data_scen_use[*,i_scen,i_model], $
            data_scen_use_mean[*,i_scen] )
        fit_est_model[*,i_scen,i_model] = data_scen_use[*,i_scen,i_model] $
            * temp_beta[0]
      endfor
    endfor
    ; Initialise variance array
    addvar_scen = one * fltarr( n_data, n_data, n_scen )
    ; Iterate through scenarios
    for i_scen = 0, n_scen - 1 do begin
      ; Iterate through data points
      for i_data_2 = 0, n_data - 1 do begin
        temp_2 = reform( fit_est_model[i_data_2,i_scen,*] ) $
            - data_scen_use_mean[i_data_2,i_scen]
        ; Iterate through data points again
        for i_data_1 = 0, n_data - 1 do begin
          temp_1 = reform( fit_est_model[i_data_1,i_scen,*] )
          addvar_scen[i_data_1,i_data_2,i_scen] = total( $
              ( temp_1 - data_scen_use_mean[i_data_1,i_scen] ) * ( temp_2 ) ) $
              / ( n_model - 1. )
        endfor
      endfor
    endfor
  endif
  ; Substitute model-mean response
  data_scen_use = temporary( data_scen_use_mean )
  ; Set effective model-mean frac_noise_var
  frac_noise_var_use = total( frac_noise_var, 2 ) / n_model
endif

; Estimate the pre-whitening operator from data_noise_1
; (The singular values^2 are n_noise_1*eigenvalues of 
; data_noise_1#transpose(data_noise_1).)
if keyword_set( weight ) then temp_weight = weight_use
pca, data_noise_1_use, evalue=noise_singval, evector=bleach, $
    double=double_opt, normalise_pc=1, no_anomaly=1, svd=gdl_opt, $
    weight_point=temp_weight, unweight=1
temp_weight = 0
noise_singval = sqrt( noise_singval )
n_bleach = n_elements( noise_singval )
; Remove eigenvalues that are effectively zero
id = where( noise_singval[1:n_bleach-1] / noise_singval[0:n_bleach-2] gt 0.01, $
    n_bleach )
noise_singval = noise_singval[0:n_bleach-1]
bleach = bleach[*,0:n_bleach-1]
; Apply inverse noise weighting if no_optimise option is not set
if no_optimise_opt eq 0 then begin
  for i_bleach = 0, n_bleach - 1 do begin
    bleach[*,i_bleach] = bleach[*,i_bleach] / noise_singval[i_bleach]
  endfor
; Otherwise apply weighting if provided
endif else if keyword_set( weight ) then begin
  for i_bleach = 0, n_bleach - 1 do begin
    bleach[*,i_bleach] = bleach[*,i_bleach] * weight_use / mean( weight_use )
  endfor
endif

; Count the number of requested truncation values
n_trunc = n_elements( trunc )
; If no truncation requests have been made
if n_trunc eq 0 then begin
  ; Set to scan truncations to find the maximum for which p_resid is greater 
  ; than p_limit
  n_trunc = n_bleach - n_scen
  trunc = n_scen + indgen( n_trunc ) + 1
; If a list of truncations is given
endif else if n_trunc gt 1 then begin
  ; Ensure truncations within permitted range
  id = where( ( trunc gt n_scen + 1 ) and ( trunc le n_bleach ), n_id )
  if n_id eq 0 then stop, 'gendetec.pro:  Impossible truncations specified.'
  if n_id lt n_trunc then begin
    print, 'gendetec.pro:  Some impossible truncations specified.  ' $
        + 'Restricting to possible truncation values only.'
    trunc = trunc[id]
    n_trunc = n_id
  endif
; If a negative truncation value is given
endif else if trunc[0] lt 0 then begin
  ; Set to scan truncations to find the maximum for which p_resid is greater 
  ; than p_limit
  n_trunc = abs( trunc[0] ) - n_scen
  if n_trunc le 0 then stop, 'gendetec.pro:  Impossible truncation specified.'
  if n_trunc gt n_bleach - n_scen then begin
    print, 'gendetec.pro:  Maximum truncation is too large.  ' $
        + 'Resetting to maximum possible value.'
    n_trunc = n_bleach - n_scen
  endif
  trunc = n_scen + indgen( n_trunc ) + 1
; If a single truncation value is given
endif else begin
  ; Ensure truncation within permitted range
  if trunc[0] le n_scen then begin
    stop, 'gendetec.pro:  Impossible truncation specified.'
  endif
  if trunc[0] gt n_bleach then begin
    print, 'gendetec.pro:  Truncation is too large.  ' $
        + 'Resetting to maximum possible value.'
    trunc = n_bleach
  endif
endelse

; Always request the residual sum of squares
temp_resid_sumsq = one

; If we need to determine the optimal truncation then set for a scanning pass 
; of linmod.pro
if n_trunc gt 1 then begin
  i_pass_0 = 0
; Otherwise skip the scanning pass
endif else begin
  i_pass_0 = 1
endelse
; Iterate through passes of linmod.pro
; (i_pass=0 is an optional scanning pass, i_pass=1 is the calculation pass)
for i_pass = i_pass_0, 1 do begin
  ; Initialise vectors for recording resid_sumsq, p_resid, and best_est values
  resid_sumsq = one * fltarr( n_trunc )
  p_resid = fltarr( n_trunc ) + one
  beta_est = one * fltarr( n_transform, 3, n_trunc )
  ; Request reconstructions and confidence intervals on second pass
  if i_pass eq 1 then begin
    tilde_obs = 1
    tilde_scen = 1
    temp_finger_trans = keyword_set( finger_trans )
    temp_onedim_indep = onedim_indep_out
    temp_onedim_project_dist = keyword_set( onedim_project_dist )
    beta_confsurf = 1
    cov_noise_2 = 1
    temp_z_pos = keyword_set( z_pos )
    temp_beta_dist = keyword_set( beta_dist )
    if keyword_set( beta_dist ) then dist_area = 1.
    if keyword_set( onedim_beta_t ) then temp_onedim_beta_t = onedim_beta_t
  ; Otherwise request the minimum necessary for getting the scanning output
  endif else begin
    temp_onedim_indep = transform_matrix
  endelse
  ; Iterate through requested truncations
  for i_trunc = 0, n_trunc - 1 do begin
    ; Determine the pre-whitening transformation
    dof_trunc = trunc[i_trunc]
    temp_bleach = transpose( bleach[*,0:dof_trunc-1] )
    ; Perform the regression
    temp = linmod( data_scen_use, data_obs_use, bleach=temp_bleach, $
        data_noise=data_noise_2_use, addvar_dep=addvar_obs, $
        scale_cov_noise=scale_cov_noise_opt, $
        frac_noise_var=frac_noise_var_use, tilde_dep=tilde_obs, $
        tilde_indep=tilde_scen, resid_sumsq=temp_resid_sumsq, $
        finger_trans=temp_finger_trans, t_crit=t_crit, $
        beta_confsurf=beta_confsurf, n_point=n_point, z_pos=temp_z_pos, $
        beta_dist=temp_beta_dist, t_dist=t_dist, dist_area=dist_area, $
        cov_noise=cov_noise_2, ols=ols_opt, type=type, $
        onedim_indep=temp_onedim_indep, $
        onedim_project_conf=onedim_project_conf, inv_bleach=bleach_inv, $
        onedim_beta_dist=temp_onedim_beta_dist, $
        onedim_beta_t=temp_onedim_beta_t, $
        onedim_project_dist=onedim_project_dist, double=double_opt, $
        addvar_indep=addvar_scen, z_dist=z_dist )
    ; Record resid_sumsq value
    resid_sumsq[i_trunc] = temp_resid_sumsq
    ; Evaluate the probability of obtaining an resid_sumsq this large if the 
    ; noise model is adequate
    if keyword_set( data_noise_2 ) then begin
      p_resid[i_trunc] = 1. - f_pdf( temp_resid_sumsq $
          / ( dof_trunc - n_scen ), $
          dof_trunc - n_scen, dof_noise )
    endif else begin
      p_resid[i_trunc] = 1. - chisqr_pdf( temp_resid_sumsq, dof_trunc - n_scen )
    endelse
    ; Record best_est value
    beta_est[*,*,i_trunc] = transpose( onedim_project_conf[*,0:n_transform-1] )
  endfor
  ; Clear memory
  if i_pass eq 1 then begin
    temp_onedim_indep = 0
  endif
  ; Record results of scan on scanning pass
  if i_pass eq 0 then begin
    scan_beta_est = beta_est
    scan_p_resid = p_resid
    scan_resid_sumsq = resid_sumsq
  endif
  ; On scanning pass, find maximum allowable truncation
  if i_pass eq 0 then begin
    id = where( p_resid gt p_limit, n_id )
    if n_id eq 0 then begin
      trunc = trunc[n_trunc-1]
    endif else begin
      trunc = trunc[id[n_id-1]]
    endelse
    n_trunc = 1
  endif
  ; Record results on the final pass
  if i_pass eq 1 then begin
    finger_trans = temporary( temp_finger_trans )
    if keyword_set( onedim_beta_axis ) then begin
      onedim_beta_dist= temporary( temp_onedim_beta_dist )
    endif
    z_pos = temporary( temp_z_pos )
    ; Retain this pre-whitening transformation
    bleach = transpose( bleach[*,0:dof_trunc-1] )
    ; Record regression coefficient distribution array
    if keyword_set( temp_beta_dist ) then begin
      beta_dist = temporary( temp_beta_dist )
    endif
  endif
endfor

; Output correlation coefficient between original and projected data
corr_obs = transpose( data_obs_use ) # bleach_inv # bleach # data_obs_use $
    / sqrt( transpose( data_obs_use ) # bleach_inv # bleach $
    # bleach_inv # bleach # data_obs_use $
    * ( transpose( data_obs_use ) # data_obs_use ) )
corr_scen = transpose( data_scen_use ) # bleach_inv # bleach # data_scen_use $
    / sqrt( transpose( data_scen_use ) # bleach_inv # bleach $
    # bleach_inv # bleach # data_scen_use $
    * ( transpose( data_scen_use ) # data_scen_use ) )
; Take diagonal elements
id = indgen( n_scen )
corr_scen = corr_scen[id,id]

; Combine the best estimates of the noise-free observations and scenarios for 
; output
z_best = transpose( [ transpose( tilde_scen ), transpose( tilde_obs ) ] )

; Extract parameter estimates and confidence limits from onedim_project_conf
if keyword_set( onedim_indep ) then begin
  onedim_project_conf = onedim_project_conf[*,0:n_onedim-1]
endif

; Convert dist_area to dist_weight
if keyword_set( beta_dist ) then begin
  p_weight = one * fltarr( n_dist )
  p_weight[0] = ( dist_prob[1] + dist_prob[0] ) / 2.
  p_weight[1:n_dist-2] = ( dist_prob[2:n_dist-1] - dist_prob[0:n_dist-3] ) / 2.
  p_weight[n_dist-1] = 1. - ( dist_prob[n_dist-2] + dist_prob[n_dist-1] ) / 2.
  dist_weight = dist_area # transpose( p_weight )
  p_weight = 0
endif

; Linearly transform output if requested
if keyword_set( transform ) then begin
  ; Calculate the inverse of the transformation matrix
  ;transform_matrix_inv $
  ;    = transpose( invert( transform_matrix, status, double=double_opt ) $
  ;    ## transpose ( transform_matrix ) )
  ;temp = transpose( transform_matrix ) ## transform_matrix
  ;transform_matrix_inv = transpose( invert( temp, status, double=double_opt ) $
  ;    ## transpose ( transform_matrix ) )
  ;transform_matrix_inv = reform( transform_matrix_inv, n_scen, n_transform )
  ;temp = 0
  ;if status ne 0 then begin
  ;  stop, 'gendetec.pro:  Failed inversion of transform_matrix.'
  ;endif
  ; Transform output
  cov_noise_2 = transpose( transform_matrix ) # cov_noise_2 # transform_matrix
  beta_confsurf = reform( beta_confsurf, n_scen, n_point*n_scen )
  beta_confsurf = transpose( transform_matrix ) # beta_confsurf
  beta_confsurf = reform( beta_confsurf, n_transform, n_point, n_scen )
  if keyword_set( beta_dist ) then begin
    beta_dist = reform( beta_dist, n_transform, n_point*n_dist )
    beta_dist = transpose( transform_matrix ) # beta_dist
    beta_dist = reform( beta_dist, n_transform, n_point, n_dist )
  endif
  ;z_best[*,0:n_transform-1] = z_best[*,0:n_scen-1] # transform_matrix_inv
  ;if keyword_set( z_pos ) then begin
  ;  for i_scen = 0, n_scen - 1 do begin
  ;    for i_point = 0, n_point - 1 do begin
  ;      z_pos[*,0:n_transform-1,i_point,i_scen] $
  ;          = reform( z_pos[*,0:n_scen-1,i_point,i_scen] ) $
  ;          # transform_matrix_inv
  ;    endfor
  ;  endfor
  ;endif
  ;if keyword_set( z_dist ) then begin
  ;  for i_dist = 0, n_dist - 1 do begin
  ;    for i_point = 0, n_point - 1 do begin
  ;      z_dist[*,0:n_transform-1,i_point,i_dist] $
  ;          = reform( z_dist[*,0:n_scen-1,i_point,i_dist] ) $
  ;          # transform_matrix_inv
  ;    endfor
  ;  endfor
  ;endif
endif

;***********************************************************************
; The end

return
END
