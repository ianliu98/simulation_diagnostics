;+
; NAME:
;	CPAR_OLS
;
; COPYRIGHT:
;	Copyright (1999) Myles Allen, Space Science Department, 
;	Rutherford Appleton Laboratory.
; 	Prepared under contract to the Hadley Centre for Climate Prediction 
;	and Research.
;
; PURPOSE:
;	This function computes 1-D and m-D parametric confidence intervals 
;	from ordinary least squares regression.
;
; CATEGORY:
;	Optimal Detection Package v2.1
;
; CALLING SEQUENCE:
;	Result = cpar_ols( Betatl, Covb2s, T_crit, Npoint )
;
; MODIFICATION HISTORY:
;	Written by:	Myles R. Allen (m.r.allen@rl.ac.uk), 1999-05-15 (v1.0)
;	Modified:	MRA, 1999-08-02 (Revise to return npoint as no. of 
;			points on C-intvls; v1.1)
;	Modified:	MRA, 1999-08-09 (Include Z_poss keyword; v1.3)
;	Modified:	MRA, 2000-08-03 (Allow threshold-type prior 
;			constraints; v2.0)
;	Modified:	Daithi A. Stone (stoned@atm.ox.ac.uk), 2004-06-28 
;			(Documentation for inclusion in routine library)
;-

function cpar_ols,betatl,covb2s,T_crit,npoint,$
 ndofn=ndofn,d1dim=d1dim,C1dim=C1dim,Z_poss=Z_poss,P_area=P_area

; Copyright (1999) Myles Allen, Space Science Department, Rutherford Appleton Laboratory
; Prepared under contract to the Hadley Centre for Climate Prediction and Research
;+
; Name: function cpar_ols
;
; Description:
; computes 1 and m-dimensional parametric conf. intervals from ols regression
;
; Method:
; svd of estimated covariance matrix, covb2s
; Indices follow standard maths notation: i.e. (row,column)
;
; History:
; Vers.	Date		Comment			Author
; ----  -------- 	--------		--------
; 1.0   15/05/99 	Original code 	Myles Allen m.r.allen@rl.ac.uk
; 1.1	02/08/99	Revise to return npoint as no. of points on C-intvls
; 1.3	09/08/99	Include Z_poss keyword
; 2.0	03/08/00	Allow threshold-type prior constraints
;
; Code Description: IDL / PV-WAVE
;
; Category: 		Function
;
; Classification keywords: regression
;
; Calling sequence: Cintvl=cpar_ols(betatl,covb2s,T_crit,npoint)
;
; Example call: 	Cintvl=cpar_ols(betatl,covb2s,T_crit,npoint,$
; 					 ndofn=ndofn,d1dim=d1dim,C1dim=C1dim)
;
; Inputs:
; 		arg1:		betatl = m-rank vector of best-fit reg. coeffs
; 		arg2:		covb2s = mxm array estimated covariance
;		arg3:		T_crit = m-rank vector of critical points on 1-D to
;                            m-D intervals
; 		arg4:		npoint = no. of points required on ellipsoids
;							 returned as np_fix^(miv-1), np_fix integer
;
; Optional Inputs:	None
;
; Keywords:
; ndofn  = d.o.f. of noise used to estimate covb2s (0 for infinite)
; npoint = no. of points to plot on confidence ellipses
; d1dim  = m*n1dim array: directions in which 1-D confidence intervals are required
; C1dim  = 3*n1dim array: best-guess, low and high scaling factors in dirn. d1dim
; Z_poss = entered as values of predictors and observations, ldp*(m+1) array
;          Z_poss(*,0:m-1) = predictors (independent variables)
; 		   Z_poss(*,m) = observations (dependent variable)
;		   returned as possible values of noise-free Z
;          ldp*(m+1)*(npoint^(m-1))*m array
;          Z_poss(*,0:m-1,j,k)=independent variables, X
;          Z_poss(*,m,j,k)=reconstructed obs, y_poss=X#Cintvl(*,j,k)
; P_area = fltarr(npoint)
;          P_area(j) = fraction of hypersurface area corresponding to jth point
;
; Outputs:			Updated keywords & function value
;
; Optional Outputs: None
; Return Value:    	Cintvl = m*(npoint^(m-1))*m array of points on confidence intervals
;                   Cintvl(*,j,k) = coordinates of jth point on (k+1)-D confidence interval
; Common Blocks: 	None
; Side Effects: 	None known
; Restrictions:
;-
;;; Just ignore the next four lines
author_name = '$Author: m.r.allen@rl.ac.uk $'
date_name = '$Date: 1999/05/15 10:30 $'
version_name = '$Revision: 1.0 $'

miv=n_elements(betatl)

; generate points on the surface of a unit miv-sphere
; first check no. of points requested is consistent with dimensionality of surface
if (miv eq 1) then begin
 if (npoint ne 2) then print,'Only 2 points required for 1-D intervals - resetting'
 npoint=2
 np_fix=npoint
endif else if (miv eq 2) then begin
 np_fix=npoint
endif else if (miv gt 2) then begin
 np_fix=fix(npoint^(1./(miv-1)))
 nprevd=float(np_fix)^(miv-1)
 if (nprevd ne npoint) then print,'No. of points revised to:',nprevd
 npoint=nprevd
endif

x=fltarr(miv,npoint)
if (keyword_set(P_area)) then P_area=replicate(1.,npoint)

; the next few lines generates a set of lat-long points on a miv-sphere
if (miv eq 1) then begin
 x(0,*)=[-1.,1.]
endif else begin
 np=fltarr(miv) & for m=1, miv do np(m-1)=float(np_fix)^(m-1.)
 z1=exp(complex(0,1)*1.*!pi*findgen(np_fix)/np_fix)
 z1(np_fix/2:np_fix-1)=-z1(np_fix/2:np_fix-1)
 z2=exp(complex(0,1)*2.*!pi*findgen(np_fix)/np_fix)
 x0=transpose([[float(z2)],[imaginary(z2)]])
 x(0:1,0:np_fix-1)=x0
 for m=2, miv-1 do begin
  indx0=indgen(np(m-1))
  for n=1, np_fix-1 do begin
   indxn=n*np(m-1)+indx0
   x(m,indxn)=imaginary(z1(n))
   x(0:m-1,indxn)=total(float(z1(n)))*x(0:m-1,indx0)
   if (keyword_set(P_area)) then $
    P_area(indxn)=abs(total(float(z1(n))))*P_area(indx0)
  endfor
 endfor
endelse

; make the elements of P_area sum to unity
if (keyword_set(P_area)) then P_area=P_area/total(P_area)

; find the number of thresholds
nthold=n_elements(T_crit)

; initialise output array
Cintvl=fltarr(miv,npoint,nthold)

; initialise array Z_poss if keyword set
if (keyword_set(Z_poss)) then begin
 ldp=n_elements(Z_poss)/(miv+1)
; store input Z_poss (original values) in Z_orig
 Z_orig=reform(Z_poss,ldp,miv+1)
 Z_poss=fltarr(ldp,miv+1,npoint,miv)
endif

; compute 1-dimensional confidence intervals
if (keyword_set(d1dim)) then begin
 n1dim=n_elements(d1dim)/miv
 C1dim=fltarr(3,n1dim)
 for n=0, n1dim-1 do begin
; compute best-guess value
  C1dim(0,n)=transpose(d1dim(*,n))#betatl
; compute uncertainty range
  d1=T_crit(0)*sqrt(transpose(d1dim(*,n))#covb2s#d1dim(*,n))
; compute min and max
  C1dim(1,n)=C1dim(0,n)-d1
  C1dim(2,n)=C1dim(0,n)+d1
 endfor
endif

; compute nthold m-dimensional confidence intervals
svdpvw,covb2s,w,u,v,/unsort
; Set u to sqrt(covb2s)
for k=0, miv-1 do u(*,k)=u(*,k)*sqrt(w(k))
for m=1, nthold do begin
 d2=T_crit(m-1)
 Cintvl(*,*,m-1)=d2*(u#x)+rebin(betatl,miv,npoint)
 if (keyword_set(Z_poss)) then begin
  for j=0, npoint-1 do begin
   Z_poss(*,0:miv-1,j,m-1)=Z_orig(*,0:miv-1)
   Z_poss(*,miv,j,m-1)=Z_orig(*,0:miv-1)#Cintvl(*,j,m-1)
  endfor
 endif
endfor

return,Cintvl
end
