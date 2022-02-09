;+
; NAME:
;	REGOLS
;
; COPYRIGHT:
;	Copyright (1999) Myles Allen, Space Science Department, 
;	Rutherford Appleton Laboratory.
; 	Prepared under contract to the Hadley Centre for Climate Prediction 
;	and Research.
;
; PURPOSE:
;	This function performs multiple ordinary least squares regression.
;
; CATEGORY:
;	Optimal Detection Package v2.1
;
; CALLING SEQUENCE:
;	Result = regols( X, Y )
;
; MODIFICATION HISTORY:
;	Written by:	Myles R. Allen (m.r.allen@rl.ac.uk), 1999-05-15 (v1.0)
;	Modified:	MRA, 2000-01-11 (Modified name of invert1 to invert1k; 
;			v1.1)
;	Modified:	Daithi A. Stone (stoned@atm.ox.ac.uk), 2004-06-28 
;			(Documentation for inclusion in routine library)
;-

function regols,X,y,Unoise=Unoise,Obsvar=Obsvar,$
                estvar=estvar,ytilde=ytilde,rssq=rssq,$
                Ftrans=Ftrans,BetaUn=BetaUn,covb2s=covb2s,status=status

; Copyright (1999) Myles Allen, Space Science Department, Rutherford Appleton Laboratory
; Prepared under contract to the Hadley Centre for Climate Prediction and Research
;+
; Name: function regols
;
; Description:
; Performs ordinary least squares regression given an optional
; externally-specified noise realisation
; Indices follow standard maths notation: i.e. (row,column)
;
; Method:
; Computes Hessian explicitly and inverts
; No regression constant: if required, input 1s in jth column of X
; Indices follow standard maths notation: i.e. (row,column)
;
; History:
; Vers.	Date		Comment			Author
; ----  -------- 	--------		--------
; 1.0   15/05/99 	Original code 	Myles Allen m.r.allen@rl.ac.uk
; 1.1   11/01/00	Modified name of invert1 to invert1k
;
; Code Description: IDL / PV-WAVE
;
; Category: 		Function
;
; Classification keywords: regression
;
; Calling sequence: b=regols(X,y) performs standard ols regression
;
; Example call: 	betatl=regols(BX,By,Unoise=BUnoise,Obsvar=BObsvar,$
; 					 estvar=estvar,ytilde=Bytilde,rssq=rssq,$
; 					 Ftrans=BFtrans,BetaUn=BetaUn,covb2s=covb2s,status=status)
; Inputs:
; 		arg1:		X = l*m array of independent variable values
; 		arg2:		y = l-rank vector of dependent variable values
;
; Optional Inputs:	None
;
; Keywords:
; Unoise = l*n array of p independent noise realisations for conf. ints.
; estvar = if set, use residual ssq to scale all variance estimates
; ytilde = l-rank vector, predicted values of y in best-fit model
; rssq   = prewhitened residual sum of squares
; BetaUn = coefficient values estimated from Unoise
; covb2s = m*m array of 2nd order estimate of covariances of coefficients
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

if (keyword_set(Obsvar)) then print,'Warning: Obsvar not yet supported'

; ldp=number of datapoints
ldp=n_elements(X(*,0))
; miv=number of independent variables
miv=n_elements(X(0,*))
; nnr=number of noise realisations
if (keyword_set(Unoise)) then nnr=n_elements(Unoise)/ldp else nnr=0

; Perform standard OLS regression on input data
covb2s=invert1k(transpose(X)#X,status=status)
; check inversion has gone OK
if (status ne 0) then return,fltarr(miv)
; Compute F^T, the matrix of fingerprints
Ftrans=covb2s#transpose(X)
; Compute best-fit regression coefficients
betatl=Ftrans#y
; Compute predicted noise-free observations
ytilde=X#betatl
; Compute residual
ztilde=y-ytilde
if (keyword_set(Unoise)) then begin
; compute estimates of betatl from independent noise realisations
 BetaUn=Ftrans#Unoise
; re-estimate covb2s from covariance of independent estimates
 covb2s=(BetaUn#transpose(BetaUn))/n_elements(BetaUn(0,*))
; re-estimate noise variance in each element of ztilde
 varUnoise=total_1d(Unoise^2,2)/nnr
 if (keyword_set(Obsvar)) then $
  for k=0, ldp-1 do varUnoise(k)=varUnoise(k)+Obsvar(k,k)
; re-normalise residuals
 ztilde=ztilde/sqrt(varUnoise)
; add betatl to BetaUn, to centre noise on estimate
 BetaUn=BetaUn+rebin(betatl,miv,nnr)
endif
; compute sum squared residual
rssq=total(ztilde^2)
; add estimated observation error to covb2s if supplied
if keyword_set(Obsvar) then covb2s=covb2s+Ftrans#Obsvar#transpose(Ftrans)
; re-estimate covariances using rssq if required
if (keyword_set(estvar)) then covb2s=covb2s*rssq/(ldp-miv)

return,betatl
end




