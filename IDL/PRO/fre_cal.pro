; new method for instantaneous frequency calculation
; 

Function  Fre_cal, by, bz, dt, AMPLITUDE=amplitude

  tlen = n_elements( by )
  period = Fltarr(tlen, /nozero)
  freq = Fltarr(tlen, /nozero)
  
  rad = Atan(bz[*], by[*])
  unrad = unwrap(rad) / (2.d * !DPI)
  
  idx = 0
  idx_prec = 0.0
  
  for it=1,tlen-1 do begin
    tmp_cycle = unrad[it]
    target_cycle = tmp_cycle + 1
    
    if ((unrad[it] - unrad[it-1]) ge 0) then begin
      for idx=idx,tlen-1 do begin
        if (unrad[idx] gt target_cycle) then begin
		  idx_prec = (target_cycle - unrad[idx-1]) / (unrad[idx] - unrad[idx-1]) + idx - 1.0
          break
        endif
      endfor
    endif else begin
      for idx=idx,0,-1 do begin
        if (unrad[idx] le target_cycle) then begin
		  idx_prec = (target_cycle - unrad[idx]) / (unrad[idx+1] - unrad[idx]) + idx
          break
        endif
      endfor
    endelse
    
    period[it-1] = idx_prec - it
    
    if (idx ge tlen-1) then begin
      break
    endif
  endfor
  
  period = period * dt
  freq = 2.d * !DPI / period
 
  if KEYWORD_SET(amplitude) then begin
	amplitude = sqrt( by[*]^2 + bz[*]^2 )
	amplitude[0] = !values.f_nan
	return, [ [freq[*]], [amplitude[*]]]
  endif else begin
  	return, freq[*]
  endelse
END
