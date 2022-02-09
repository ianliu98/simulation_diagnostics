;+
; NAME:
;    BETAS_HIST
;
; COPYRIGHT:
;    Copyright (2005) Daithi Stone, University of Oxford.
;    Revised (2011) Daithi Stone under contract to the U.S. Department of 
;    Energy's Office of Science, Office of Biological and Environmental 
;    Research and the U.S. National Oceanic and Atmospheric Administration's 
;    Climate Program Office, via the International Detection and Attribution 
;    Group.
;
; PURPOSE:
;    This function produces n-dimensional histograms of the beta scaling 
;    factors from the gendetec Optimal Detection Package.
;
; CATEGORY:
;    Optimal Detection Package v3.1.2
;
; CALLING SEQUENCE:
;    Result = betas_hist( BETA_DIST, DIST_WEIGHT )
;
; INPUTS:
;    BETA_DIST:  An array containing the location (beta values) of the 
;        probability density estimates of the beta scaling parameters.  Of size 
;        [N_SCEN,N_POINT,N_PROB] where N_SCEN is the number of scenarios, 
;        N_POINT is the number of points sampled along an isopleth of 
;        probability, and N_PROB is the number of isopleths sampled.  So 
;        elements [*,j,k] are the coordinates of the jth point on the kth 
;        isopleth of probability.
;    DIST_WEIGHT:  An array containing the weighting estimates of the beta 
;        scaling parameters at the locations in BETA_DIST.  Note these are not 
;        actual density estimates because BETA_DIST contains an irregular 
;        sampling (polar).  Of size [N_POINT,N_PROB], so elements [j,k] are the 
;        fraction of the PDF at the jth point on the kth isopleth of 
;        probability.
;
; KEYWORD PARAMETERS:
;    D_HIST:  A vector containing the resolution of the histogram along each 
;        scenario.  Of length [N_SCEN] (the number of scenarios).  A default 
;        can be used and returned by the function.
;    HIST_AXIS:  An array containing the values to sample for each scenario.  
;        Of size [N_HIST,N_SCEN].  A default can be used and returned by the 
;        function.
;    HIST_RANGE:  An array containing the range of values to sample for each 
;        scenario.  Of size [2,N_SCEN], where the first element is the minimum 
;        and the second element is the maximum.  A default can be used and 
;        returned by the function.
;    N_HIST:  The number of values to sample along each scenario.  The default 
;        is 100.
;
; OUTPUTS:
;    Result:  An array containing the histogram estimator of the probability 
;        density of the scaling factors.  Of N_SCEN dimensions, each with 
;        N_HIST elements.  The corresponding axes are in HIST_AXIS.
;    D_HIST, HIST_AXIS, HIST_RANGE, N_HIST
;
; USES:
;    ---
;
; PROCEDURE:
;    This function converts the information in BETA_DIST into coordinates in 
;        the RESULT array, and copies the DIST_WEIGHT values into these 
;        coordinates in the RESULT array.
;
; EXAMPLE:
;    See demo_gendetec.pro.
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone (stoned@atm.ox.ac.uk), 2005-03-13
;    Modified:  DAS, 2005-03-16 (implemented faster 1-D algorithm, fixed 
;        normalisation bug)
;    Modified:  DAS, 2005-03-18 (implemented a much, much faster algorithm;  
;        v3.0 of Optimal Detection Package)
;    Modified:  DAS, 2011-10-14 (Edited documentation formatting;  altered code 
;        variable names;  changed input and output variable names (see 
;        translation notes below);  v3.1.0)
;    Modified:  DAS, 2012-01-24 (Edited for compliance with GDL;  v3.1.1)
;-

;***********************************************************************

FUNCTION BETAS_HIST, $
    BETA_DIST, DIST_WEIGHT, $
    D_HIST=d_hist, HIST_AXIS=hist_axis, HIST_RANGE=hist_range, N_HIST=n_hist

; Changes in input/output variable names from v3.0.0 to v3.1.0
;   changed dhist to d_hist
;   changed histaxis to hist_axis
;   changed histrange to hist_range
;   changed nhist to n_hist
;   changed w_dist to dist_weight

;***********************************************************************
; Constants

; Number of scenarios (dimensions) included in the input
n_scen = n_elements( beta_dist[*,0,0] )
; Number of sampling points
n_point = n_elements( beta_dist[0,*,0] )
; The number of isopleths of probability to sample
n_prob = n_elements( beta_dist[0,0,*] )
; The total number of points on all isopleths of probability
n_dist = n_point * n_prob

; Range to sample for each scenario.
; Check if ranges were inputted
if keyword_set( hist_range ) then begin
  ; Check that the right number of ranges were given
  if n_elements( hist_range ) ne n_scen * 2 then stop
  ; Reform to range-scenario format
  hist_range = reform( hist_range, 2, n_scen )
  ; A reminder to not update this range later
  check = 0
; If ranges were not inputted but hist_axis was input
endif else if keyword_set( hist_axis ) then begin
  ; Initialise range array
  hist_range = fltarr( 2, n_scen )
  ; Iterate through scenarios
  for i_scen = 0, n_scen - 1 do begin
    ; Determine a range for this scenario
    hist_range[*,i_scen] = [ min( hist_axis, max=temp ), temp ]
  endfor
  ; A reminder to not update this range later
  check = 0
endif else begin
  ; Initialise range array
  hist_range = fltarr( 2, n_scen )
  ; Iterate through scenarios
  for i_scen = 0, n_scen - 1 do begin
    ; Determine a range for this scenario
    hist_range[*,i_scen] = [ min( beta_dist[i_scen,*,*], max=temp ), temp ]
  endfor
  ; A reminder to update this range later
  check = 1
endelse
hist_range = reform( hist_range, 2, n_scen )
; The default number of points in each dimension of our histogram of the betas
if not( keyword_set( n_hist ) ) then begin
  if keyword_set( hist_axis ) then begin
    n_hist = n_elements( hist_axis[*,0] )
  endif else begin
    n_hist = 100
  endelse
endif
; The beta sampling step size for each scenario
if not( keyword_set( d_hist ) ) then begin
  d_hist = reform( hist_range[1,*] - hist_range[0,*] ) / 1. / n_hist
endif
; An update on the histogram axis ranges if automatic
if check eq 1 then begin
  ; Update to ensure extreme points are included (numerically)
  hist_range[0,*] = hist_range[0,*] - d_hist / 2.
  hist_range[1,*] = hist_range[1,*] + d_hist / 2.
endif
; Create an array of beta values to sample in our histogram
if not( keyword_set( hist_axis ) ) then begin
  hist_axis = reform( fltarr( n_hist, n_scen ), n_hist, n_scen )
  for i_scen = 0, n_scen - 1 do begin
    hist_axis[*,i_scen] = findgen( n_hist ) / n_hist $
        * ( hist_range[1,i_scen] - hist_range[0,i_scen] ) $
        + d_hist[i_scen] / 2. + hist_range[0,i_scen]
  endfor
endif

;***********************************************************************
; Estimate the Histogram

; Reform the n_scen dimensional distributional information so that all 
; distributional information is in a single dimension for each scenario
beta_dist_use = transpose( reform( beta_dist, n_scen, n_dist ) )
dist_weight = reform( dist_weight, n_dist )

; Initialise our histogram
hist = fltarr( n_hist ^ n_scen )

; Iterate through scenarios
for i_scen = 0, n_scen - 1 do begin
  ; Convert beta_dist_use values to indices
  beta_dist_use[*,i_scen] = ( beta_dist_use[*,i_scen] - hist_axis[0,i_scen] ) $
      / ( hist_axis[n_hist-1,i_scen] - hist_axis[0,i_scen] ) * ( n_hist - 1 )
endfor
; Retain only points inside of the axis range
for i_scen = 0, n_scen - 1 do begin
  id = where( ( beta_dist_use[*,i_scen] gt -0.5 ) $
      and ( beta_dist_use[*,i_scen] lt n_hist - 0.5 ) )
  beta_dist_use = beta_dist_use[id,*]
endfor
n_dist = n_elements( beta_dist_use[*,0] )
; Convert beta_dist_use values to indices of hist
beta_dist_use = round( beta_dist_use )
beta_dist_use = ( n_hist ^ indgen( n_scen ) ) ## beta_dist_use
; Iterate through distributional values
for i_dist = 0l, n_dist - 1l do begin
  ; Add this weighting to our histogram
  id = beta_dist_use[i_dist]
  hist[id] = hist[id] + dist_weight[i_dist]
endfor

; Reform histogram to n_scen-dimensional array
hist = reform( hist, intarr( n_scen ) + n_hist )
dist_weight = reform( dist_weight, n_point, n_prob )
; Normalise the histogram
hist = hist / product( d_hist )

;***********************************************************************
; The End

return, hist
END
