Function MOVMEAN, a, k

n = n_elements(a)
m = fltarr(n)
for i=0,n-1 do begin
	il = max([ceil(i - k/2.d), 0])
	ih = min([ceil(i + k/2.d - 1), n-1])
	m[i] = mean(a[il:ih])
endfor

return, m

END
