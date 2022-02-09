;+
; NAME:
;	CPAR_TLS
;
; COPYRIGHT:
;	Copyright (1999) Myles Allen, Space Science Department, 
;	Rutherford Appleton Laboratory.
; 	Prepared under contract to the Hadley Centre for Climate Prediction 
;	and Research.
;
; PURPOSE:
;	This function computes 1-D and m-D parametric confidence intervals 
;	from total least squares regression.
;
; CATEGORY:
;	Optimal Detection Package v2.1
;
; CALLING SEQUENCE:
;	Result = cpar_tls( Betatl, Evects, Evalus, T_crit, Npoint )
;
; MODIFICATION HISTORY:
;	Written by:	Myles R. Allen (m.r.allen@rl.ac.uk), 1999-05-15 (v1.0)
;	Modified:	MRA, 1999-08-02 (Revise to return npoint as no. of 
;			points on C-intvls; v1.1)
;	Modified:	MRA, 1999-08-02 (Simplify computation of C-intvls; 
;			v1.2)
;	Modified:	MRA, 1999-08-09 (Include Z_poss keyword; v1.3)
;	Modified:	Daithi A. Stone (stoned@atm.ox.ac.uk), 2004-06-28 
;			(Documentation for inclusion in routine library)
;-

function cpar_tls,betatl,Evects,Evalus,T_crit,npoint,$
 d1dim=d1dim,C1dim=C1dim,t1dim=t1dim,Z_poss=Z_poss,P_area=P_area

; Copyright (1999) Myles Allen, Space Science Department, Rutherford Appleton Laboratory
; Prepared under contract to the Hadley Centre for Climate Prediction and Research
;+
; Name: function cpar_tls
;
; Description:
; computes 1-D and m-D parametric confidence intervals from TLS regression
;
; Method:
; see techreport
; Indices follow standard maths notation: i.e. (row,column)
;
; History:
; Vers.	Date		Comment			Author
; ----  -------- 	--------		--------
; 1.0   15/05/99 	Original code 	Myles Allen m.r.allen@rl.ac.uk
; 1.1   02/08/99	Revise to return npoint as no. of points on C-intvls
; 1.2   02/08/99	Simplify computation of C-intvls
; 1.3	09/08/99	Include Z_poss keyword
; 2.0	future		Allow threshold-type prior constraints
;
; Code Description: IDL / PV-WAVE
;
; Category: 		Function
;
; Classification keywords: regression
;
; Calling sequence: Cintvl=cpar_tls(betatl,Evects,Evalus,T_crit,npoint)
;
; Example call: 	Cintvl=cpar_ols(betatl,covb2s,T_crit,npoint,$
; 					 d1dim=d1dim,C1dim=C1dim,t1dim=t1dim,Z_poss=Z_poss)
;
; Inputs:
; 		arg1:		betatl = m-rank vector of best-fit reg. coeffs
; 		arg2:		Evects = eigenvectors of transpose([X,y])#[X,y]
; 		arg3:		Evalus = eigenvalues  of ditto
;		arg3:		T_crit = m-rank vector of critical points on 1-D to
;                            m-D intervals
; 		arg4:		npoint = no. of points required on ellipsoids
;							 returned as np_fix^(miv-1), np_fix integer
;
; Optional Inputs:	None
;
; Keywords:
; d1dim  = m*n1dim array: directions in which 1-D confidence intervals are required
; C1dim  = 3*n1dim array: best-guess, low and high scaling factors in dirn. d1dim
; t1dim  = n1dim vector: threshold values of t-statistic, returned as sqrt(rmin)
; Z_poss = entered as values of predictors and observations, ldp*(m+1) array
;          Z_poss(*,0:m-1) = predictors (independent variables)
; 		   Z_poss(*,m) = observations (dependent variable)
;		   returned as possible values of noise-free Z
;          ldp*(m+1)*(npoint^(m-1))*m array
;          Z_poss(*,*,j,k)=best-fit values of predictors and observations
;          corresponding to Cintvl(*,j,k)
; P_area = fltarr(npoint)
;          P_area(j) = fraction of hypersurface area corresponding to jth point
;
; Outputs:			Updated keywords & function value
;
; Optional Outputs: None
; Return Value:    	Cintvl = m*(npoint^(m-1))*m array of points on confidence intervals
; 					Cintvl(*,j,k) = coordinates of jth point on (k+1)-D confidence interval
; 					if (Cintvl(*,j,k) eq 0.) then confidence interval unbounded
; Common Blocks: 	None
; Side Effects: 	None known
; Restrictions:
;-
;;; Just ignore the next four lines
author_name = '$Author: m.r.allen@rl.ac.uk $'
date_name = '$Date: 1999/08/09 10:30 $'
version_name = '$Revision: 1.3 $'

test=1

miv=n_elements(betatl)

; generate points on the surface of a sphere of radius 1.
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

; initialise array of points
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

; initialise working and output arrays
b=fltarr(miv+1,npoint)
Cintvl=fltarr(miv,npoint,nthold)
if (keyword_set(d1dim)) then begin
 n1dim=n_elements(d1dim)/miv
 C1dim=fltarr(3,n1dim)
endif else n1dim=0

if (keyword_set(t1dim)) then begin
 if (n1dim eq 0) then stop,'d1dim required in cpar_tls if t1dim set'
 if (n_elements(t1dim) ne n1dim) then t1dim=replicate(0.,n1dim)
endif

; initialise array Z_poss if keyword set
if (keyword_set(Z_poss)) then begin
 ldp=n_elements(Z_poss)/(miv+1)
; store input Z_poss (original values) in Z_orig
 Z_orig=reform(Z_poss,ldp,miv+1)
 Z_poss=fltarr(ldp,miv+1,npoint,nthold)
endif

; compute nthold miv-D confidence intervals -- simplified to work in Evects coordinates
; m=1 corresponds to 1-D confidence interval etc.
for m=1, nthold do begin
; rescale the unit sphere to radius T_crit
 a=x*T_crit(m-1)
; compute weights on 1st miv e-vectors by dividing a by sqrt(delta-evalus)
 for i=0, miv-1 do b(i,*)=a(i,*)/sqrt(Evalus(i)-Evalus(miv))
; compute sum squared weights
 f=total_1d(b(0:miv-1,*)^2,1)
; check if the sum if greater than unity (indicates open-ended interval)
 iimag=where(f gt 1.,nimag)
 if (nimag gt 0) then begin
  f(iimag)=1.
  print,'Warning: open-ended C-intvl:',m
 endif
 ireal=where(f lt 1.,nreal)
; Don't compute anything if no bounds on C-intvl
 if (nreal eq 0) then print,'Warning: no bounds on C-intvl:',m else begin
; compute weight on (miv+1)-th e-vector from normalisation constraint
; use positive sqrt to give max projection onto last e-vector
  b(miv,*)=sqrt(1.-f)
; transform to normal coordinates
  v=Evects#b
; convert to conventional regression-coefficient-like scaling parameters
  for i=0, miv-1 do Cintvl(i,ireal,m-1)=-v(i,ireal)/v(miv,ireal)
; compute 1-D confidence limits
  if (m eq 1 and (keyword_set(d1dim) or keyword_set(t1dim))) then begin
   for n=0, n1dim-1 do begin
; project 1-D directions onto ellipsoids
    C1dim(0,n)=transpose(d1dim(*,n))#betatl
    d1d=transpose(d1dim(*,n))#reform(Cintvl(*,ireal,m-1))
    C1dim(1,n)=min(d1d)
    C1dim(2,n)=max(d1d)
    test=1
    if (keyword_set(t1dim)) then begin
     if (test) then t1dim(n)=C1dim(1,n)
; define the vector c (see notes)
     c=[d1dim(*,n),t1dim(n)]
; define a complete set of vectors orthogonal to c
     D=-c#transpose(c)/total(c^2)
     for k=0, miv do D(k,k)=D(k,k)+1.
     svdpvw,D,Q,P,R
     D2=transpose(Evects)#P(*,0:miv-1)
; weight the vectors by the Evalu difference
     for i=0, miv-1 do D2(i,*)=D2(i,*)*sqrt(Evalus(i)-Evalus(miv))
     v2=D2#x
     t1dim(n)=sqrt(min(total_1d(v2^2,1)))
     if (test) then print,t1dim(n),T_crit(m-1)
    endif
   endfor
  endif
; compute Z_poss if requested
  if (keyword_set(Z_poss)) then begin
   for j=0, nreal-1 do begin
    Z_poss(*,*,ireal(j),m-1)=Z_orig - $
     (Z_orig#v(*,ireal(j)))#transpose(v(*,ireal(j)))
   endfor
  endif
 endelse
endfor

return,Cintvl
end

; code to check [a,b2] is a solution
;  for n=0, np_fix-1 do begin
;   if (f(n) lt 1.) then begin
;    vsoln=a(n)*Evects(*,miv)
;    vsoln(i)=vsoln(i)+b2(0,n)
;    vsoln(j)=vsoln(j)+b2(1,n)
;    print,'transpose(vsoln)#vsoln=',transpose(vsoln)#vsoln
;    print,'transpose(vsoln)#E_miv=',transpose(vsoln)#Evects(*,miv)
;    vswgt=transpose(Evwgts)#vsoln
;    print,'transpose(vswgt)#vswgt=',transpose(vswgt)#vswgt
;   endif
;  endfor

; code to extract 2-D from miv-D confidence intervals
; else begin
;  for i=0, miv-1 do begin
;   for j=i+1, miv-1 do begin
;    plane2d=cellip([i,j],*)
;    plane2d=plane2d(*,i_def)
; just extract 2-D curve if miv eq 2
;    if (miv eq 2) then cell2d=plane2d else begin
; use x-coord/y-coord to avoid x=0
;     pratio=plane2d(0,*)/plane2d(1,*)
;     xratio=x0(0,*)/x0(1,*)
;     rxsign=pratio*xratio
;     cell2d=fltarr(2,np_fix)
;     for n=0, np_fix-1 do begin
;      idirn=where(pratio gt xratio(n) and pratio lt xratio(n+1) and rxsign gt 0.,ndirn)
;      if (ndirn gt 0) then begin
;       sub_plane=plane2d(*,idirn)
;       dist1d=total_1d(sub_plane^2,1)
;       cell2d(*,n)=plane1d(*,where(dist1d eq max(dist1d)))
;       defd2d(n)=1.
;      endif
;     endfor
;    endelse
;    Cintvl(i,j,i_def)=cell2d(0,*)
;    Cintvl(j,i,i_def)=cell2d(1,*)
;   endfor
;  endfor
; endelse
