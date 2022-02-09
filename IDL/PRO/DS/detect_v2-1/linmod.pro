;+
; NAME:
;	LINMOD
;
; COPYRIGHT:
;	Copyright (1999) Myles Allen, Space Science Department, 
;	Rutherford Appleton Laboratory.
; 	Prepared under contract to the Hadley Centre for Climate Prediction 
;	and Research.
;
; PURPOSE:
;	This function performs multiple linear regression given an optional 
;	prewhitening operator, externally specified noise realisation, and 
;	optional noise on the independent variables.
;
; CATEGORY:
;	Optimal Detection Package v2.1
;
; CALLING SEQUENCE:
;	Result = linmod( X, Y )
;
; MODIFICATION HISTORY:
;	Written by:	Myles R. Allen (m.r.allen@rl.ac.uk), 1999-05-15 (v1.0)
;	Modified:	MRA, 1999-08-02 (Bug-fix on rssq MRA, bug spotted by 
;			PAS; v1.1)
;	Modified:	MRA, 1999-08-09 (Z_poss keyword implemented; v1.3)
;	Modified:	MRA, 2000-08-03 (B_dist keyword implemented: explicit 
;			output of PDFs; v2.0)
;	Modified:	Daithi A. Stone (stoned@atm.ox.ac.uk), 2004-06-28 
;			(Documentation for inclusion in routine library)
;-

function linmod,X,y,Bleach=Bleach,Unoise=Unoise,Obsvar=Obsvar,$
 estvar=estvar,xnoise=xnoise,nonpar=nonpar,ndofn=ndofn,ndofd=ndofd,$
 Xtilde=Xtilde,ytilde=ytilde,rssq=rssq,Ftrans=Ftrans,$
 T_crit=T_crit,Cintvl=Cintvl,npoint=npoint,Z_poss=Z_poss,$
 B_dist=B_dist,T_dist=T_dist,P_area=P_area,$
 BetaUn=BetaUn,covb2s=covb2s,ols=ols,d1dim=d1dim,C1dim=C1dim,$
 status=status,colour=colour

; Copyright (1999) Myles Allen, Space Science Department, Rutherford Appleton Laboratory
; Prepared under contract to the Hadley Centre for Climate Prediction and Research
;+
; Name: function linmod
;
; Description:
; Performs multiple linear regression given an optional
; prewhitening operator, externally-specified
; noise realisation and optional noise on indep. variables.
; Uniform weighting with no missing values
; (absorb weighting in prewhitening operator if desired)
;
; Method: fits the model
; y = (X+U_x)b - u_y where \expect{u_y u_y^T} = C_N
; C_N^{-1} = B^T B, so \expect{B u_y u_y^T B^T} = I_n'
; i.e. B u_y = z is unit variance white noise
; and \expect{B U_x S^{-2} U_x^T B^T} = m * I_n' where S(k,k)^2=xnoise(k)
; using either regols.pro if xnoise(*)=0. or regtls.pro otherwise
; No regression constant: if required, input 1s in jth column of X
; and set xnoise(j)=0. or (numerically better) subtract all means
; before calling and work it out manually.
; Indices follow standard maths notation: i.e. (row,column)
;
; History:
; Vers.	Date		Comment			Author
; ----  -------- 	--------		--------
; 1.0   15/05/99 	Original code 	Myles Allen m.r.allen@rl.ac.uk
; 1.1   02/08/99	Bug-fix on rssq MRA, bug spotted by PAS
; 1.3   09/08/99	Z_poss keyword implemented
; 2.0   03/08/00	B_dist keyword implemented: explicit output of PDFs
;
; Code Description: IDL / PV-WAVE
;
; Category: 		Function
;
; Classification keywords: detection, pre-whitening
;
; Calling sequence: betatl=linmod(X,y) performs standard ols regression
;
; Example call: 	bobs1=linmod(X,y,Bleach=Bleach,Unoise=Unoise,obsvar=obsvar,$
; 					 estvar=estvar,xnoise=xnoise,nonpar=nonpar,ndofn=ndofn,ndofd=ndofd,$
; 					 Xtilde=Xtilde,ytilde=ytilde,rssq=rssq,Ftrans=Ftrans,$
; 					 T_crit=T_crit,Cintvl=Cintvl,npoint=npoint,$
; 					 BetaUn=BetaUn,Covb2s=Covb2s,ols=ols,d1dim=d1dim1,C1dim=C1dim1,status=status)
; Inputs:
; 		arg1:		X = l*m array of independent variable values
; 		arg2:		y = l-rank vector of dependent variable values
;
; Optional Inputs:	None
;
; Keywords:
; Bleach = prewhitening operator (e.g. W_N^{-1} V_N^T where W and V
;          are noise EOFs and singular values respectively
; Unoise = l*n array of p independent noise realisations for conf. ints.
; estvar = if set, use residual ssq to scale all variance estimates
; xnoise = if set, xnoise(k) = expected variance(BX_k)/variance(By)
; nonpar = if set, use non-parametric confidence intervals based on Unoise
; ndofn  = no. degrees of freedom of noise estimate (-> inf. if not set)
; ndofd  = no. degrees of freedom of prewhitened data (rank of Bleach)
; Xtilde = l*m array, predicted values of X in best-fit model
; ytilde = l-rank vector, predicted values of y in best-fit model
; rssq   = prewhitened residual sum of squares
; Ftrans = Transposed fingerprint matrix, or operator used to obtain betatl
; T_crit = m-rank vector of critical points on 1-D to m-D confidence intervals
; Cintvl = m*npoint*m array of confidence intervals
; Cintvl(*,j,k) = coordinates of jth point on (k+1)-D interval
; npoint = number of points on ellipsoids (defaults to 2 if m=1, np_fix^(m-1) otherwise)
; Z_poss = returned as possible values of noise-free Z
;          ldp*(m+1)*(npoint^(m-1))*m array
;          Z_poss(*,0:m-1,j,k)=best-fit values of predictors, X,
;          corresponding to Cintvl(*,j,k)
;          Z_poss(*,m,j,k)=best-fit values of observations, y,
;          corresponding to Cintvl(*,j,k)
; B_dist=fltarr(nscn,npoint,ndist) = PDF of estimators
; B_dist(*,j,k)=coordinates of jth point on kth isopleth of probability
; P_area=fltarr(npoint)
; P_area(j)=fraction of area on isopleths corresponding to jth point
; T_dist=fltarr(ndist)=T-statistic on requested isopleths of probability
; BetaUn = coefficient values estimated from Unoise
; covb2s = if xnoise not set, m*m array of covariances of coefficients
;        = if xnoise set, (m+1)*(m+1) array of covs of normalised coeffs
;          (for advanced users only -- see regrma documentation)
; ols    = use ordinary least squares algorithm even if xnoise is set
; d1dim  = m*n1dim array = parameter combins. on which uncertainties are reqd.
; C1dim  = 3*n1dim array = best, lower and upper estimates of composite parameters
; status = 0 if routine completed OK, 1 otherwise
; colour = pseudo-inverse of Bleach.
;
; Outputs:			Updated keywords & function value
;
; Optional Outputs: None
; Return Value:    	betatl = m-rank vector of regression coefficients
; Common Blocks: 	None
; Side Effects: 	None known
; Restrictions:
;-
;;; Just ignore the next four lines
author_name = '$Author: m.r.allen@rl.ac.uk $'
date_name = '$Date: 1999/05/15 10:30 $'
version_name = '$Revision: 1.0 $'

; ldp=number of datapoints
ldp=n_elements(X(*,0))
; miv=number of independent variables
miv=n_elements(X(0,*))
; nnr=number of noise realisations
if (keyword_set(Unoise)) then nnr=n_elements(Unoise)/ldp else nnr=0
; set default d.o.f. of Unoise to nnr
if (not keyword_set(ndofn)) then ndofn=nnr
; initialise estvar
if (not keyword_set(estvar)) then estvar=0.
; initialise Obsvar
if (not keyword_set(Obsvar)) then Obsvar=0.
; set default T_crit to sqrt(m)
if (not keyword_set(T_crit)) then T_crit=sqrt(findgen(miv)+1.)
; set default null-hypothesis to zero
if (not keyword_set(b_null)) then b_null=fltarr(miv)
; set default no. of points on Cintvl
if (not keyword_set(npoint)) then $
 if (miv eq 1) then npoint=2. else $
 if (miv eq 2) then npoint=100. else $
  npoint=fix(25000.^(1./(miv-1)))
; initialise BXtilde and bytilde if necessary
if (keyword_set(Xtilde)) then BXtilde=1.
if (keyword_set(ytilde)) then Bytilde=1.
if (keyword_set(Ftrans)) then BFtrans=1.
if (keyword_set(d1dim)) then n1dim=n_elements(d1dim)/miv else n1dim=0

; Prewhitening step:
if (keyword_set(Bleach)) then begin
; set number of degrees of freedom of data to rank of Bleach
 ndofd=n_elements(Bleach)/ldp
 if (ndofd gt ldp) then stop,'Bleach operator cannot be rank',ndofd
 By=Bleach#y
 BX=Bleach#X
 if (keyword_set(Unoise)) then BUnoise=Bleach#Unoise else BUnoise=0.
 if (keyword_set(Obsvar)) then BObsvar=Bleach#Obsvar#transpose(Bleach) $
                          else BObsvar=0.
 if (keyword_set(Xtilde) or keyword_set(ytilde) or keyword_set(Z_poss) or keyword_set(colour)) then begin
; Colour operator = pseudo-inverse of Bleach
  svdpvw,Bleach,Q,P,R,/unsort
  for k=0, ndofd-1 do begin
   if (Q(k) le 0.) then stop,'Rank-deficient prewhitening operator'
   P(*,k)=P(*,k)/Q(k)
  endfor
  Colour=R#transpose(P)
 endif
endif else begin
 ndofd=ldp
 By=y
 BX=X
 if (keyword_set(Unoise)) then BUnoise=Unoise else BUnoise=0.
endelse

if (keyword_set(xnoise) and not keyword_set(ols)) then begin

; Check all elements of xnoise are non-zero
 for m=0, miv-1 do begin
  if (xnoise(m) lt 1.e-5) then begin
   print,'Warning: xnoise input with zeros - resetting to 1e-5',xnoise
   xnoise(m)=1.e-5
  endif
  sx=sqrt(xnoise(m))
; scale i.v.s by xnoise
  BX(*,m)=BX(*,m)/sx
; scale directions for 1-D confidence intervals by xnoise
  if (keyword_set(d1dim)) then d1dim(m,*)=d1dim(m,*)/sx
 endfor
 if (keyword_set(Cintvl)) then begin
; request eigenvector and eigenvalue information from regrma
  Evects=1.
  Evalus=1.
 endif
; Perform RMA regression on prewhitened data and scaled i.v.s
 betatl=regtls(BX,By,Unoise=BUnoise,Obsvar=BObsvar,$
               estvar=estvar,Xtilde=BXtilde,ytilde=Bytilde,rssq=rssq,$
               Ftrans=BFtrans,BetaUn=BetaUn,covb2s=covb2s,$
               Evects=Evects,Evalus=Evalus,status=status)
; Compute confidence intervals using Evects and Evalus information
; Store original values in Z_poss if required
 if (keyword_set(Z_poss)) then Z_poss=transpose([transpose(BX),transpose(By)]) $
  else Z_poss=0.
 t1dim=0.
 if (keyword_set(Cintvl) or keyword_set(d1dim) or keyword_set(Z_poss)) then $
  Cintvl=cpar_tls(betatl,Evects,Evalus,T_crit,npoint,$
                  d1dim=d1dim,C1dim=C1dim,t1dim=t1dim,Z_poss=Z_poss)
 if (keyword_set(B_dist)) then $
  B_dist=cpar_tls(betatl,Evects,Evalus,T_dist,npoint,P_area=P_area)
; rescale everything by xnoise
 for m=0, miv-1 do begin
; recompute scaling factor
  sx=sqrt(xnoise(m))
; reconstruction multiplied by sx since original data divided by sx
  if (keyword_set(Xtilde)) then BXtilde(*,m)=BXtilde(*,m)*sx
  if (keyword_set(Z_poss)) then Z_poss(*,m,*,*)=Z_poss(*,m,*,*)*sx
; estimates divided by sx, to correspond to rescaled data
  betatl(m)=betatl(m)/sx
  if (keyword_set(B_dist)) then B_dist(m,*,*)=B_dist(m,*,*)/sx
  if (keyword_set(BetaUn)) then BetaUn(m,*)=BetaUn(m,*)/sx
; rows and columns of covariance estimates divided by sx
  if (keyword_set(covb2s)) then begin
   covb2s(*,m)=covb2s(*,m)/sx
   covb2s(m,*)=covb2s(m,*)/sx
  endif
; rows of Cintvl matrix divided by sx (scaled as betatl)
  if (keyword_set(Cintvl)) then Cintvl(m,*,*)=Cintvl(m,*,*)/sx
  if (keyword_set(d1dim)) then d1dim(m,*)=d1dim(m,*)*sx
 endfor

endif else begin

; Perform standard OLS regression on prewhitened data
 betatl=regols(BX,By,Unoise=BUnoise,Obsvar=BObsvar,$
               estvar=estvar,ytilde=Bytilde,rssq=rssq,$
               Ftrans=BFtrans,BetaUn=BetaUn,covb2s=covb2s,status=status)
 BXtilde=BX
 if (keyword_set(xnoise)) then begin
; inflate rows and columns of covb2s to account for noise in BX (fudge)
  covb2s=covb2s*(sqrt(1.+xnoise)#transpose(sqrt(1.+xnoise)))
; deflate residual by mean noise in BX (bigger fudge)
; tests are applied as if this noise is not present
  rssq=rssq/(1.+total(xnoise)/miv)
 endif
; stop
; compute confidence intervals
; Store original values in Z_poss if required
 if (keyword_set(Z_poss)) then Z_poss=transpose([transpose(BX),transpose(By)]) $
  else Z_poss=0.
 if (keyword_set(Cintvl) or keyword_set(d1dim)) then $
  Cintvl=cpar_ols(betatl,covb2s,T_crit,npoint,$
                  d1dim=d1dim,C1dim=C1dim,Z_poss=Z_poss)
 if (keyword_set(B_dist)) then $
  B_dist=cpar_ols(betatl,covb2s,T_dist,npoint,P_area=P_area)
endelse

; Compute best-fit X, y and Fingerprint matrix in original coordinates
if (keyword_set(Bleach)) then begin
 if (keyword_set(ytilde)) then ytilde=Colour#Bytilde
 if (keyword_set(Xtilde)) then Xtilde=Colour#BXtilde
 if (keyword_set(Ftrans)) then Ftrans=BFtrans#Bleach
 if (keyword_set(Z_poss)) then begin
  Z_poss=reform(Z_poss,ndofd,(miv+1)*npoint*miv)
  Z_poss=Colour#Z_poss
  Z_poss=reform(Z_poss,ldp,miv+1,npoint,miv)
 endif
endif else begin
 if (keyword_set(ytilde)) then ytilde=Bytilde
 if (keyword_set(Xtilde)) then Xtilde=BXtilde
 if (keyword_set(Ftrans)) then Ftrans=BFtrans
endelse

return,betatl
end
