;+
; NAME:
;	REGTLS
;
; COPYRIGHT:
;	Copyright (1999) Myles Allen, Space Science Department, 
;	Rutherford Appleton Laboratory.
; 	Prepared under contract to the Hadley Centre for Climate Prediction 
;	and Research.
;
; PURPOSE:
;	This function performs multiple total least squares regression.
;
; CATEGORY:
;	Optimal Detection Package v2.1
;
; CALLING SEQUENCE:
;	Result = regtls( X, Y )
;
; MODIFICATION HISTORY:
;	Written by:	Myles R. Allen (m.r.allen@rl.ac.uk), 1999-05-15 (v1.0)
;	Modified:	Daithi A. Stone (stoned@atm.ox.ac.uk), 2004-06-28 
;			(Documentation for inclusion in routine library)
;-

function regtls,X,y,Unoise=Unoise,Obsvar=Obsvar,$
 estvar=estvar,Xtilde=Xtilde,ytilde=ytilde,rssq=rssq,$
 Ftrans=Ftrans,BetaUn=BetaUn,covb2s=covb2s,$
 Evects=Evects,Evalus=Evalus,status=status

; Copyright (1999) Myles Allen, Space Science Department, Rutherford Appleton Laboratory
; Prepared under contract to the Hadley Centre for Climate Prediction and Research
;+
; Name: function regtls
;
; Description:
; Performs multiple total-least-squares regression given an optional
; externally-specified noise realisation and equal noise on indep. variables.
; Updated 17-Nov-1998, uniform weighting
; Model: y = (X+U_x)b - u_y where \expect{u_y u_y^T} = I
; and \expect{U_x U_x^T} = m * I
;
; Method:
; SVD on data matrix, Z=[X;y]
; No regression constant: if required, input 1s in jth column of X
; Indices follow standard maths notation: i.e. (row,column)
;
; History:
; Vers.	Date		Comment			Author
; ----  -------- 	--------		--------
; 1.0   15/05/99 	Original code 	Myles Allen m.r.allen@rl.ac.uk
;
; Code Description: IDL / PV-WAVE
;
; Category: 		Function
;
; Classification keywords: regression
;
; Calling sequence: b=regols(X,y) performs standard ols regression
;
; Example call: 	betatl=regtls(BX,By,Unoise=BUnoise,Obsvar=BObsvar,$
; 					 estvar=estvar,Xtilde=BXtilde,ytilde=Bytilde,rssq=rssq,$
; 					 Ftrans=BFtrans,BetaUn=BetaUn,covb2s=covb2s,$
; 					 Evects=Evects,Evalus=Evalus,status=status)
;
; Inputs:
; 		arg1:		X = l*m array of independent variable values
; 		arg2:		y = l-rank vector of dependent variable values
;
; Optional Inputs:	None
;
; Keywords:
; Unoise = l*n array of p independent noise realisations for conf. ints.
; Obsvar = observation error covariance (not supported)
; estvar = if set, use residual ssq to scale all variance estimates
; Xtilde = l*m array, predicted values of X in best-fit model
; ytilde = l-rank vector, predicted values of y in best-fit model
; rssq   = prewhitened residual sum of squares
; Ftrans = Transposed fingerprint matrix, or LSV with min e-value
; BetaUn = coefficient values estimated from Unoise
; covb2s   = approximate m*m array of covariances of coefficients
; Evects = eigenvectors of transpose([X,y])#[X,y]
; Evalus = eigenvalues  of ditto
; status = 0 if routine completed OK, 1 otherwise
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
if (keyword_set(Obsvar)) then stop,'Obsvar keyword not yet supported in regtls'

; Create Z=[X,y], remembering IDL's indexing convention
Z=transpose([transpose(X),transpose(y)])
; diagonalise
svdpvw,Z,Q,P,R
; compute best-fit regression coefficients from min RSV
ycoeff=R(miv,miv)
if (ycoeff eq 0.) then stop,'Zero weight on observations in regtls',R(*,miv)
betatl=-reform(R(0:miv-1,miv))/ycoeff
if (keyword_set(Xtilde) or keyword_set(ytilde)) then begin
; compute best-fit model from all SVs orthogonal to R(*,miv)
 Ztilde=Z-Z#R(*,miv)#transpose(R(*,miv))
 Xtilde=Ztilde(*,0:miv-1)
 ytilde=reform(Ztilde(*,miv))
endif
; set Ftrans to the operator which extracts betatl from Z
if (keyword_set(Ftrans)) then Ftrans=transpose(P(*,miv))
if (keyword_set(Unoise)) then begin
; Revise variance in each W(k) to reflect independent noise estimate
 varUnoise=total_1d(Unoise^2,2)/nnr
 for k=0, miv do Q(k)=Q(k)/sqrt(transpose(P(*,k)^2)#varUnoise)
endif
; compute sum-squared residual
rssq=total(Q(miv)^2)
; Revise variance estimates to residuals if estvar is set
if (keyword_set(estvar)) then Q=Q/sqrt(rssq/(ldp-miv))
; Store eigenvectors and eigenvalues if requested
if (keyword_set(Evects)) then Evects=R
if (keyword_set(Evalus)) then Evalus=Q^2
; initialise new eigenvectors matrices if required
if (keyword_set(covb2s)) then R1=fltarr(miv+1,miv+1)
if (keyword_set(BetaUn)) then R2=fltarr(miv+1,miv+1)
for k=0, miv-1 do begin
 q2dif=total(Q(k)^2-Q(miv)^2)
 if (q2dif lt 0.) then begin
  print,'Warning: degenerate evalues in regtls -- resetting to 0'
  Q(k)=Q(miv)
  status=1.
 endif else begin
  if (keyword_set(covb2s)) then R1(*,k)=R(*,k)/sqrt(q2dif)
  if (keyword_set(BetaUn)) then R2(*,k)=R(*,k)*sqrt(q2dif)
 endelse
endfor

; Estimate of covariance
if (keyword_set(covb2s)) then begin
; Compute covariance of Rmin
 covRmin=R1#transpose(R1)
; Compute second-order approximation to covariance of betatl
 rmin=reform(R(*,miv)) & rmin_miv=total(rmin(miv))
 covb2s=covRmin(0:miv-1,0:miv-1)/rmin_miv^2 - $
      covRmin(0:miv-1,miv)#transpose(rmin(0:miv-1))/rmin_miv^3 - $
      rmin(0:miv-1)#covRmin(miv,0:miv-1)/rmin_miv^3 + $
      rmin(0:miv-1)#transpose(rmin(0:miv-1))*covRmin(miv,miv)/rmin_miv^4
endif

; Non-parametric estimate of parameter spread
if (keyword_set(BetaUn)) then begin
 nnr=n_elements(Unoise)/ldp
 BetaUn=fltarr(miv,nnr)
; reconstruct noise-reduced data
 Z2=P#transpose(R2)
; loop over noise realisations
 for n=0, nnr-1 do begin
  ZN=Z2
; contaminate observations, ZN(*,miv), with nth noise segment
  ZN(*,miv)=ZN(*,miv)+Unoise(*,k)
; generate an nnr-length randomized sequence
  index=sort(randomu(seed,nnr))
; eliminate element n
  index=index(where(index ne n))
; contaminate indep. vars., ZN(0:miv-1,*), with random noise segments
  ZN(0:miv-1,*)=ZN(0:miv-1,*)+Unoise(*,index(0:miv-1))
; compute estimates of betatl for noise realisations
  svdpvw,ZN,Q,P,R
  BetaUn(*,n)=R(0:miv-1,miv)/total(R(miv,miv))
 endfor
endif

return,betatl
end
