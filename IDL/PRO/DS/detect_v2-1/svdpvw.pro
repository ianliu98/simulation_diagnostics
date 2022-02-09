;+
; NAME:
;	SVDPVW
;
; COPYRIGHT:
;	Copyright (1999) Myles Allen, Space Science Department, 
;	Rutherford Appleton Laboratory.
; 	Prepared under contract to the Hadley Centre for Climate Prediction 
;	and Research.
;
; PURPOSE:
;	This procedure is a singular value decomposition (SVD) routine for 
;	PV-Wave.
;
; CATEGORY:
;	Optimal Detection Package v2.1
;
; CALLING SEQUENCE:
;	SVDPVW, A, W, U, V
;
; MODIFICATION HISTORY:
;	Written by:	Myles R. Allen (m.r.allen@rl.ac.uk), 1999-05-15 (v1.0)
;	Modified:	Daithi A. Stone (stoned@atm.ox.ac.uk), 2004-06-28 
;			(Documentation for inclusion in routine library)
;-

pro svdpvw,a,w,u,v,unsort=unsort

; Copyright (1999) Myles Allen, Space Science Department, Rutherford Appleton Laboratory
; Prepared under contract to the Hadley Centre for Climate Prediction and Research
;+
; Name: function svdpvw
;
; Description:
; svd routine for PV-wave which transposes as necessary
; output is in standard math notation
;
; Method:
; uses the svd,a,w,u,v command after reforming and transposing
; as necessary
; Indices follow standard maths notation: i.e. (row,column)
;
; History:
; Vers.	Date		Comment			Author
; ----  -------- 	--------		--------
; 1.0   15/05/99 	Original code 	Myles Allen m.r.allen@rl.ac.uk
;
; Code Description: IDL / PV-WAVE
;
; Category: 		program
;
; Classification keywords: general linear algebra
;
; Calling sequence: svdpvw,a,w,u,v
;
; Example call: 	svdpvw,a,w,u,v
;
; Inputs:
; 		arg1:		a = real mxn array
;
; Optional Inputs:	None
;
; Keywords:
; unsort = if set and positive, do not sort singular values (slightly faster)
;          if set and negative, sort singular values in ascending order
;
; Outputs:			k = smaller of m and n
; 		arg2:		w = k-rank vector of singular values
;		arg3:		u = m x k matrix of LSVs, stored columnwise
; 		arg4:		v = n x k matrix of RSVs, stored columnwise
;
; Optional Outputs: None
; Return Value:    	N/A
; Common Blocks: 	None
; Side Effects: 	None known
; Restrictions: 	Not tested for complex matrices
;-
;;; Just ignore the next four lines
author_name = '$Author: m.r.allen@rl.ac.uk $'
date_name = '$Date: 1999/05/15 10:30 $'
version_name = '$Revision: 1.0 $'

if (not keyword_set(unsort)) then unsort=0
sa=size(a)
if (sa(0) eq 0) then stop,'svdpvw not set up for scalar input'
if (sa(0) eq 1) then sa(2)=1
if (sa(0) gt 2) then stop,'svdpvw not set up for >2-D arrays'
if (sa(1) gt sa(2)) then $
 svd,reform(a,sa(1),sa(2)),w,u,v else $
 svd,reform(transpose(a),sa(2),sa(1)),w,v,u
if (unsort le 0) then begin
; find ordering of singular values
 index=sort(w)
; default is to sort in descending order
 if (unsort eq 0) then index=reverse(index)
 w=w(index)
 su=size(u)
 sv=size(v)
 for k=0, su(1)-1 do u(k,*)=u(k,index)
 for k=0, sv(1)-1 do v(k,*)=v(k,index)
endif
return
end
