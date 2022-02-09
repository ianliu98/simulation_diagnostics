; routine to cgplot
; loadct, num  in advance

Function myplot, data, minv, maxv, interp, axes, bart

cgwindow, wxsize=960, wysize=960
cgplot, data, pos=[.1,.1,.9,.9], interpolate=interp, $
	minvalue=minv, maxvalue=maxv, /axis, axkeywords=axes, /window

cgcolorbar, /vertical, /right, pos=[.91,.1,.92,.9], $
	range=[minv,maxv], title=bartitle, /addcmd

End
