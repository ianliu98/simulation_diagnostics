;+
; NAME:
;	GENDETEC
;
; COPYRIGHT:
;	Copyright (1999) Myles Allen, Space Science Department, 
;	Rutherford Appleton Laboratory.
; 	Prepared under contract to the Hadley Centre for Climate Prediction 
;	and Research.
;
; PURPOSE:
;	This procedure computes confidence intervals on undetermined model
;	parameters by multiple regression using observational constraints.  
;
; CATEGORY:
;	Optimal Detection Package v2.1
;
; CALLING SEQUENCE:
;	gendetec, Obs, Scn, Ctl [, Bobs] [, Rssq]
;
; MODIFICATION HISTORY:
;	Written by:	Myles R. Allen (m.r.allen@rl.ac.uk), 1999-05-15 (v1.0)
;	Modified:	MRA, 1999-08-09 (Z_POSS keyword implemented; v1.3)
;	Modified:	MRA, 1999-08-13 (Name revised for DOS file transfers; 
;			v1.4)
;	Modified:	MRA, 1999-01-11 (Bug in use of weights fixed; v1.5)
;	Modified:	MRA, 2000-08-03 (Explicit PDFs computed; v2.0)
;	Modified:	Daithi A. Stone (stoned@atm.ox.ac.uk), 2004-06-28 
;			(Documentation for inclusion in routine library)
;-

pro gendetec,obs,scn,ctl,bobs,rssq,pv_wave=pv_wave,$
 ctlind=ctlind,weight=weight,trunc=trunc,Pobs=Pobs,Covb=Covb,$
 dofctr=dofctr,Presid=Presid,prewhi=prewhi,ols=ols,$
 Plimit=Plimit,Cintvl=Cintvl,npoint=npoint,icom=icom,$
 Z_best=Z_best,Z_poss=Z_poss,$
 P_dist=P_dist,B_dist=B_dist,W_dist=W_dist,$
 no_opt=no_opt,ctlsvl=ctlsvl,xnoise=xnoise,d1dim=d1dim,$
 C1dim=C1dim,ry=ry,rX=rX

; Copyright (1999) Myles Allen, Space Science Department, Rutherford Appleton Laboratory
; Prepared under contract to the Hadley Centre for Climate Prediction and Research
;+
; Name: Program gendetec (was gendetect)
;
; Description:
; Given a vector, obs, of observations containing no missing data,
; a set of similar vectors, scn, from model simulations with changing forcing and a
; control set, ctl, of pure noise, computes confidence intervals on undetermined model
; parameters by multiple regression using either the ordinary (ols) or total (tls)
; least squares algorithm.
; Method: This routine is primarily a driver for linmod.pro, but computes prewhitening
; transformations and transformations on input/output variables if uncertainties
; on parameter-combinations are required, and scans for maximum allowable
; truncation if not specified
;
; History:
; Vers.	Date		Comment			Author
; ----  -------- 	--------		--------
; 1.0   15/05/99 	Original code 	Myles Allen m.r.allen@rl.ac.uk
; 1.3	   09/08/99	Z_poss keyword implemented
; 1.4	   13/08/99	Name revised for DOS file transfers
; 1.5   11/01/99	Bug in use of weights fixed
; 2.0   03/08/00	Explicit PDFs computed
;
; Code Description: IDL / PV-WAVE
;
; Category: 		Program
;
; Classification keywords: detection, pre-whitening
;
; Calling sequence: 5 arguments + keywords
;
; Example call: 	gendetect,obs,scn,ct1,bobs,rssq,$
; 					ctlind=ct2,weight=wgt,trunc=trunc,Pobs=Pobs,Covb=Covb,$
; 					dofctr=dofctr,Presid=Presid,prewhi=prewhi,ols=ols,$
; 					Plimit=Plimit,Cintvl=Cintvl,npoint=npoint,Z_poss=Z_poss,icom=icom,$
; 					no_opt=no_opt,ctlsvl=ctlsvl,xnoise=xnoise,d1dim=d1dim,C1dim=C1dim
;
; Inputs:
; 		arg1:		obs(imx)=observations
; 		arg2:		scn(imx,nscn)=runs with changing forcing: nscn experiments (scenarios)
; 		arg3:		ctl(imx,nct1)=control run for estimating pre-whitening transformation
;
; Optional Inputs:	None
;
; Keywords:
; pv_wave=use PV-wave stats routine names
; ctlind=fltarr(imx,nct2), if set, independent control run for hypothesis-testing
; weight=fltarr(imx) = vector of weights
; trunc=number of EOFs of the control run to retain in the noise model
;  if not set or negative, then max truncation for which Presid > Plimit used
;  if negative, then max truncation considered = abs(trunc), returned as truncation used
; Pobs=fltarr(imx,nscn)=Projection operators used to extract pattern amps from obs in ols
; Covb=fltarr(nscn,nscn)=estimated covariance of bobs (approximate unless ols is used)
; dofctr=degrees of freedom of noise covariance used for uncertainties (input only)
; Presid=P-value of F or chi-squared test on residuals
; Presid(0)=value for max allowed truncation (or only truncation if keyword_set(trunc)
; Presid(nscn+1:trmax)=values for truncations nscn+1:trmax if scanning over truncations
; prewhi=externally specified prewhitening transformation
; ols=use ordinary least squares estimators (still prewhitened but not tls)
; Plimit=P-value for confidence intervals
; Cintvl=fltarr(nscn,npoint,nscn):
;  Cintvl(*,j,k)=coords of jth point in (k+1)-D confidence interval
; npoint=no. of points on nscn-D confidence intervals
; Z_best = best-guess values of noise-free Z, ldp*(m+1) array
; Z_poss = returned as possible values of noise-free Z
;          ldp*(m+1)*npoint*m array
;          Z_poss(*,0:m-1,j,k)=best-fit values of predictors, X,
;          corresponding to Cintvl(*,j,k)
;          Z_poss(*,0:m-1,j,k)=best-fit values of observations, y,
;          corresponding to Cintvl(*,j,k)
; P_dist=fltarr(ndist)=requested isopleths of probability
;      can be entered as scalar = number of equal-spaced isopleths in [0,1] interval
; B_dist=fltarr(nscn,npoint,ndist) = PDF of estimators
; B_dist(*,j,k)=coordinates of jth point on kth isopleth of probability
; W_dist=fltarr(npoint,ndist) = array of weights for PDF
; W_dist(j,k)=fraction of total PDF corresponding to jth point on kth isopleth
; icom=combinations of coefficients to output
; e.g. icom=[12,2] -> output [bobs(0)+bobs(1),bobs(1)]
; no_opt=do not optimise fingerprints (but still project onto noise EOFs)
; ctlsvl=singular values of matrix of control segments
; xnoise=fltarr(nscn)=fraction variance of obs in nscn response-patterns
; d1dim(nscn,n1dim)=directions for 1-D confidence intervals (before combination in icom)
; C1dim(3,n1dim)=mid,min&max points on 1-D intervals
; ry=correlation coefficient between obs before and after projection on control EOFs.
; rX=vector of correlation coefficients of signal patterns before and after projection
; on control EOFs.
;
; Outputs:			Updated keywords &
; 		arg4:		bobs(nscn,3)=vector of best, lower and upper estimates of pattern-amplitudes
; 		arg5:		rssq=weighted residual sum of squares (trunc-nscn d.o.f.)
;       rssq(0)=value for max allowed truncation (or only truncation if keyword_set(trunc)
;       rssq(nscn+1:trmax)=values for truncations nscn+1:trmax if scanning over truncations
; Optional Outputs: None
; Return Value:    	N/A
; Common Blocks: 	None
; Side Effects: 	None known
; Restrictions:
;-
;;; Just ignore the next four lines
author_name = '$Author: m.r.allen@rl.ac.uk $'
date_name = '$Date: 1999/05/15 10:30 $'
version_name = '$Revision: 1.0 $'

if (keyword_set(pv_wave)) then begin
 @stat_startup
endif

if (not keyword_set(no_opt)) then no_opt=0

; use independent estimate of control if available
if (keyword_set(ctlind)) then ct2=ctlind else ct2=ctl
if (not keyword_set(weight)) then wgt=replicate(1.,n_elements(obs)) else wgt=weight

; Extract parameters of observational dataset, control and model simulations
s1=size(obs) & s2=size(scn) & s3=size(ctl) & s4=size(ct2) & s5=size(wgt)
if (s1(0) ne 1) then begin print,'obs must be 1-d in gendetect' & stop & end
if (s2(0) ne 2) then begin print,'scn must be 2-d in gendetect' & stop & end
if (s3(0) ne 2) then begin print,'ctl must be 2-d in gendetect' & stop & end
if (s4(0) ne 2) then begin print,'ct2 must be 2-d in gendetect' & stop & end
if (s5(0) ne 1) then begin print,'wgt must be 1-d in gendetect' & stop & end
; imx=number of datapoints in total in each field
imx=s1(1)
if (s2(1) ne imx) then begin print,'obs and scn incompatible' & stop & end
if (s3(1) ne imx) then begin print,'obs and ctl incompatible' & stop & end
if (s4(1) ne imx) then begin print,'obs and ct2 incompatible' & stop & end
if (s5(1) ne imx) then begin print,'obs and wgt incompatible' & stop & end
; nscn=number of experiments (ensembles) with changing forcing
nscn=s2(2)
; nct1=number of fields in the control model run for prewhitening
nct1=s3(2)
; nct2=number of fields in the control model run for hypothesis-testing
nct2=s4(2)

; multiply everything by wgt
obswgt=obs
scnwgt=scn
ct1wgt=ctl
ct2wgt=ct2
for i=0, imx-1 do begin
 obswgt(i)=obs(i)*wgt(i)
 scnwgt(i,*)=scn(i,*)*wgt(i)
 ct1wgt(i,*)=ctl(i,*)*wgt(i)
 ct2wgt(i,*)=ct2(i,*)*wgt(i)
endfor

; Find the EOFs of the control to estimate Bleach, the EOFs of ct1 weighted
; by inverse(sqrt(eigenvalue))
; Use IDL/PV-wave svdpvw program, which is the same as svdc with /column
; except that it transposes rather than zero-padding
; Remember the singular values^2 are nct1*eigenvalues of ct1#transpose(ct1)

if (keyword_set(prewhi)) then begin
; set U to be the left singular vectors of prewhi
 svdpvw,prewhi,W,U,V
 for k=0, n_elements(W)-1 do begin
  if (W(k) le 0.) then stop,'Singular prewhitening transformation supplied'
  W(k)=1./W(k)
 endfor
endif else begin
 svdpvw,ct1wgt,W,U,V
 W=W/sqrt(nct1)
endelse

; extract control singular values if requested
if (keyword_set(ctlsvl)) then ctlsvl=W

; Apply inverse noise weighting if (no_opt) not set
if (not keyword_set(no_opt)) then for k=0, n_elements(W)-1 do U(*,k)=U(*,k)/W(k)

; compute linear transformation on output
tr_out=fltarr(nscn,nscn)
icom_flag=0
if (keyword_set(icom)) then begin
 if (n_elements(icom) ne nscn) then stop,'Inconsistent icom keyword in gendetect:',icom
 for k=0, nscn-1 do begin
; unpack the integers in icom(k)
  for j=0, fix(alog10(icom(k))) do begin
   i=fix(icom(k)/10^j)-10*fix(icom(k)/10^(j+1))
   if (i ne k+1) then icom_flag=1
   if (i eq 0 or i gt nscn) then stop,'Inconsistent icom keyword in gendetect:',icom
   tr_out(i-1,k)=1.
  endfor
 endfor
endif else for k=0, nscn-1 do tr_out(k,k)=1.

; compute 1-D directions in which confidence intervals are required
; set last nscn dirns equal to tr_out
if (keyword_set(d1dim)) then begin
 n1dim=n_elements(d1dim(0,*))
 ndirn=n_elements(d1dim(*,0))
 if (ndirn ne nscn) then begin
  d1dim1=fltarr(nscn,n1dim)
  d1dim1(0:ndirn-1,*)=d1dim
 endif else d1dim1=d1dim
 d1dim1=transpose([transpose(d1dim1),transpose(tr_out)])
endif else begin
 n1dim=0
 d1dim1=tr_out
endelse

; Perform the regression
if (not keyword_set(nonpar)) then nonpar=0.

; Fix bug in v1_4 which didn't pick up the weighted version
X=scnwgt
y=obswgt
if (keyword_set(ctlind)) then Unoise=ct2wgt else Unoise=0.
if (keyword_set(obsvar)) then obsvar=obsvar else obsvar=0.
if (keyword_set(estvar)) then estvar=estvar else estvar=0.
if (keyword_set(xnoise)) then xnoise=xnoise else xnoise=0.
if (keyword_set(dofctr)) then ndofn=dofctr else ndofn=nct2

if (keyword_set(Plimit)) then Plimit=Plimit else Plimit=0.05
; compute critical values of the t- or T-statistic corresponding to Plimit
T_crit=fltarr(nscn)
for m=1, nscn do begin
 if (keyword_set(pv_wave)) then $
  T_crit(m-1)=sqrt(m*Fcdf((1.-Plimit),m,ndofn,/inverse)) $
 else $
  T_crit(m-1)=sqrt(m*f_cvf(Plimit,m,ndofn))
endfor

if (keyword_set(B_dist)) then begin
; compute critical values of T-statistic corresponding to P_dist
 if (not keyword_set(P_dist)) then P_dist=100
 if (n_elements(P_dist) eq 1) then begin
  ndist=total(P_dist)
  P_dist=(0.5+findgen(ndist))/ndist
 endif else begin
; P_dist must be in ascending order
  P_dist=P_dist(sort(P_dist))
  ndist=n_elements(P_dist)
  if (P_dist(0) le 0. or P_dist(ndist-1) ge 1.) then stop,'Impossible P-values:',P_dist(0),P_dist(ndist-1)
 endelse
 T_dist=fltarr(ndist)
 for k=0, ndist-1 do begin
  if (keyword_set(pv_wave)) then $
  T_dist(k)=sqrt(nscn*Fcdf((1.-P_dist(k)),nscn,ndofn,/inverse)) $
 else $
  T_dist(k)=sqrt(nscn*f_cvf(P_dist(k),nscn,ndofn))
 endfor
endif

; if trunc not set or negative, scan to determine maximum truncation for which
; Presid gt Plimit
; set maximum truncation to rank of pre-whitening operator unless specified otherwise
if (not keyword_set(trunc)) then trmax=n_elements(W) else $
 if (trunc lt 0) then trmax=abs(trunc)<n_elements(W) else trmax=0
; stop if max truncation less than nscn is requested
if (trmax ne 0 and trmax le nscn) then stop,'Impossible truncation specified:',trmax

; initialise arrays to store rssq and Presid
rssq=fltarr(trmax+1)
Presid=fltarr(trmax+1)+1.

; comment this line in to use default for number of points on ellipses
; npoint=0
; always request the residual sum of squares
rssq_k=1.

; make two passes, the first to determine the truncation, only if trmax ne 0
if (trmax gt 0) then np0=0 else np0=1
for np=np0, 1 do begin
 if (np eq 0) then begin
; do not request reconstructions and confidence intervals
  Xtilde=0.
  ytilde=0.
  Ftrans=0.
  d1dimk=0.
  C1dim1=0.
  Cintvl=0.
  BetaUn=0.
  Covb2s=0.
  Z_poss=0.
  trunc_k=nscn+indgen(trmax-nscn)+1
  P_area=0.
 endif else begin
; request reconstructions and confidence intervals
  Xtilde=1.
  ytilde=1.
  Ftrans=1.
  d1dimk=d1dim1
  C1dim1=1.
  Cintvl=1.
  BetaUn=1.
  Covb2s=1.
  Z_poss=1.
  trunc_k=[trunc]
  if (keyword_set(B_dist)) then P_area=1.
 endelse

 for k=0, n_elements(trunc_k)-1 do begin
; Compute the pre-whitening transformation
  Bleach=transpose(U(*,0:trunc_k(k)-1))
  ndofd=trunc_k(k)
  bobs1=linmod(X,y,Bleach=Bleach,Unoise=Unoise,obsvar=obsvar,$
   estvar=estvar,xnoise=xnoise,nonpar=nonpar,ndofn=ndofn,ndofd=ndofd,$
   Xtilde=Xtilde,ytilde=ytilde,rssq=rssq_k,Ftrans=Ftrans,$
   T_crit=T_crit,Cintvl=Cintvl,npoint=npoint,Z_poss=Z_poss,$
   B_dist=B_dist,T_dist=T_dist,P_area=P_area,$
   BetaUn=BetaUn,Covb2s=Covb2s,ols=ols,d1dim=d1dimk,C1dim=C1dim1,$
   status=status,colour=colour)
; on first pass, store results in rssq(trunk_k(k))
; on second pass, store in rssq(0) and Presid(0)
  if (np eq 0) then k_t=trunc_k(k) else k_t=k
  rssq(k_t)=rssq_k
; evaluate probability of obtaining an rssq this large if noise model adequate
  if (keyword_set(ctlind)) then begin
   if (keyword_set(pv_wave)) then $
    Presid(k_t)=1.-fcdf(rssq_k/(ndofd-nscn),ndofd-nscn,ndofn) else $
    Presid(k_t)=1.-f_pdf(rssq_k/(ndofd-nscn),ndofd-nscn,ndofn)
  endif else begin
   if (keyword_set(pv_wave)) then $
    Presid(k_t)=1.-chisqcdf(rssq_k,ndofd-nscn) else $
    Presid(k_t)=1.-chisqr_pdf(rssq_k,ndofd-nscn)
  endelse
 endfor
; on first pass, find maximum allowable truncation
 if (np eq 0) then begin
  trunc=max(where(Presid gt Plimit,count))
  if (count eq 0) then trunc=trmax
 endif
endfor
;
;NPG - output correlation coefficient between original and projected data.
y=reform(y,imx,1)
ry=transpose(y)#Colour#Bleach#y/$
  sqrt(transpose(y)#Colour#Bleach#Colour#Bleach#y*(transpose(y)#y))
rX=transpose(X)#Colour#Bleach#X/$
  sqrt(transpose(X)#Colour#Bleach#Colour#Bleach#X*(transpose(X)#X))
;Take diagonal elements
rX=rX(indgen(nscn),indgen(nscn))


Pobs=Ftrans
Covb=Covb2s
Bctl=BetaUn
Z_best=transpose([transpose(Xtilde),transpose(ytilde)])

; extract parameter estimates and limits from C1dim1
if (keyword_set(d1dim)) then C1dim=C1dim1(*,0:n1dim-1)
; transformation already included in d1dim
bobs=transpose(C1dim1(*,n1dim:n1dim+nscn-1))

; convert P_area to W_dist, weights on distribution
if (keyword_set(B_dist)) then begin
 P_wgt=fltarr(ndist)
 P_wgt(0)=(P_dist(1)+P_dist(0))/2.
 for n=1, ndist-2 do P_wgt(n)=(P_dist(n+1)-P_dist(n-1))/2.
 P_wgt(ndist-1)=1.-(P_dist(ndist-2)+P_dist(ndist-1))/2.
 W_dist=P_area#transpose(P_wgt)
endif

; transform output if necessary
if (icom_flag) then begin
 Covb=transpose(tr_out)#Covb#tr_out
 sz=size(Cintvl)
 Cintvl=reform(Cintvl,sz(1),sz(2)*sz(3))
 Cintvl=transpose(tr_out)#Cintvl
 Cintvl=reform(Cintvl,sz(1),sz(2),sz(3))
 if (keyword_set(B_dist)) then begin
  sz=size(B_dist)
  B_dist=reform(B_dist,sz(1),sz(2)*sz(3))
  B_dist=transpose(tr_out)#B_dist
  B_dist=reform(B_dist,sz(1),sz(2),sz(3))
 endif
 tr_inv=transpose(invert1k(tr_out,status=status))
 Z_best(*,0:nscn-1)=Z_best(*,0:nscn-1)#tr_inv
 if (keyword_set(Z_poss)) then begin
  for m=0, nscn-1 do begin
   for j=0, npoint-1 do begin
    Z_poss(*,0:nscn-1,j,m)=reform(Z_poss(*,0:nscn-1,j,m))#tr_inv
   endfor
  endfor
 endif
endif

return
end
