;+
; NAME:
;    LINE_LEGEND
;
; PURPOSE:
;    This procedure plots a legend for line plots.
;
; CATEGORY:
;    Graphics
;
; CALLING SEQUENCE:
;    LINE_LEGEND, Position, Label
;
; INPUTS:
;    Position:  A 2-element vector containing the position ([x,y]) in normal 
;        coordinates of the lower left corner of the box containing the legend.
;    Label:  A vector, of type string, containing the labels for each line.
;
; OUTPUTS:
;    LIMIT
;
; KEYWORD PARAMETERS:
;    BACKGROUND:  The colour index of a background for the legend box.  If not 
;        set, no background is drawn.
;    BORDER_COLOR:  An input integer giving the colour index of the border 
;        line.  The default is !p.color.
;    CHARSIZE:  The size of the label characters, of type floating point.  The 
;        default is the IDL default (!p.charsize).
;    COLOR:  A vector, of type integer, containing the colour index values of 
;        the lines.  The default is the IDL default (!p.color).
;    FONT:  An integer specifying the graphics font to use.  The default is the 
;        IDL default (!p.font).
;    LENGTH:  The length of the lines in normal coordinates.
;    LIMIT:  An output floating point vector containing the coordinates of the 
;        edges of the legend box.  The format is [ left, bottom, right, top ].
;    LINESTYLE:  A vector, of type integer, containing the linestyle index 
;        value for each line.  The default is the IDL default (!p.linestyle).
;    NOBORDER:  If set then no border line is drawn around the legend.  The 
;        default is for a border box of colour BORDER_COLOR.
;    PSYM:  A vector, of type integer, containing the symbol codes.  The 
;        default is the IDL default (!p.psym).
;    THICK:  A vector, of type integer, containing the line thickness value for 
;        each line.  The default is the IDL default (!p.thick).
;    TITLE:  A string containing the title of the legend.
;
; USES:
;    -
;
; PROCEDURE:
;    This procedure uses the input values to construct an appropriate box 
;    containing the legend for a line plot.
;
; EXAMPLE:
;    Create a legend for a two-line plot (colours red and green).
;      tek_color
;      line_legend, [0.2,0.2], ['Red','Green'], color=[2,3], title='Legend'
;
; MODIFICATION HISTORY:
;    Written by:  Daithi A. Stone (stoned@atm.ox.ac.uk), 2000-09-18.
;    Modified:  DAS, 2003-02-05 (added FONT and PSYM keywords, converted LINE 
;        keyword to LINESTYLE).
;    Modified:  DAS, 2003-06-25 (added BACKGROUND keyword).
;    Modified:  DAS, 2005-04-12 (fixed output bugs in CHARSIZE, COLOR, 
;        LINESTYLE, PSYM, THICK keywords)
;    Modified:  DAS, 2009-10-01 (fixed so only one symbol plotted per label 
;        instead of two when PSYM given)
;    Modified:  DAS, 2010-02-19 (fixed bug which always set PSYM)
;    Modified:  DAS, 2010-09-16 (added NOBORDER keyword option, LIMIT keyword 
;        output;  removed use of var_type.pro;  modified documentation format)
;    Modified:  DAS, 2010-09-26 (fixed bug in PSYM keyword implementation;  
;        added BORDER_COLOR keyword)
;    Modified:  DAS, 2013-02-01 (set text to print in BORDER_COLOR colour)
;    Modified:  DAS, 2015-01-08 (fixed bug in drawing border when BACKGROUND
;        is input)
;-

;***********************************************************************

PRO LINE_LEGEND, $
    Position, $
    Label, $
    BACKGROUND=background, BORDER_COLOR=border_color, NOBORDER=noborder_opt, $
    CHARSIZE=charsize, $
    COLOR=color, $
    FONT=font, $
    LENGTH=length, $
    LINESTYLE=linestyle, $
    PSYM=psym, $
    THICK=thick, $
    TITLE=title, $
    LIMIT=limit

;***********************************************************************
; Variables and Options

; Number of lines
nlabels = n_elements( label )

; Line colours (COLOR keyword)
if not( keyword_set( color ) ) then begin
  color0 = !p.color + intarr( nlabels )
endif else begin
  color0 = color
  ; Assure sufficient number defined
  if n_elements( color0 ) lt nlabels then begin
    color0 = color0[0] + intarr( nlabels )
  endif
endelse

; Graphics font (FONT keyword)
if n_elements( font ) eq 0 then font = !p.font

; Line styles (LINESTYLE keyword)
if not( keyword_set( linestyle ) ) then begin
  linestyle0 = !p.linestyle + intarr( nlabels )
endif else begin
  linestyle0 = linestyle
  ; Assure sufficient number defined
  if n_elements( linestyle0 ) lt nlabels then begin
    linestyle0 = linestyle0[0] + intarr( nlabels )
  endif
endelse

; Line thicknesses (THICK keyword)
if not( keyword_set( thick ) ) then begin
  thick0 = !p.thick + intarr( nlabels )
endif else begin
  thick0 = thick
  ; Assure sufficient number defined
  if n_elements( thick0 ) lt nlabels then begin
    thick0 = thick0[0] + intarr( nlabels )
  endif
endelse

; Symbols (PSYM keyword)
if not( keyword_set( psym ) ) then begin
  psym0 = !p.psym + intarr( nlabels )
endif else begin
  psym0 = psym
  ; Assure sufficient number defined
  if n_elements( psym0 ) lt nlabels then psym0 = psym0[0] + intarr( nlabels )
endelse

; Line lengths (LENGTH keyword)
if not( keyword_set( length ) ) then length = 0.02

; Determine the legend box corner positions
xpos = position[[0,0]]
ypos = position[[1,1]]
if !p.multi[1] ne 0 then begin
  if !p.multi[0] eq 0 then pmulti0 = !p.multi[1] * !p.multi[2] $
                      else pmulti0 = !p.multi[0]
  xpos = ( pmulti0 + !p.multi[1] - 1 ) / !p.multi[1] $
         - pmulti0 / 1. / !p.multi[1] + xpos / !p.multi[1]
endif
if !p.multi[2] ne 0 then begin
  if !p.multi[0] eq 0 then pmulti0 = !p.multi[1] * !p.multi[2] $
                      else pmulti0 = !p.multi[0]
  ypos = ( ypos + ( pmulti0 - 1 ) / !p.multi[1] ) / 1. / !p.multi[2]
endif

; Character scale
if not( keyword_set( charsize ) ) then begin
  charsize0 = !p.charsize
endif else begin
  charsize0 = charsize
endelse
if charsize0 eq 0. then charsize0 = 1.
if !p.multi[2] gt 1 then charsize0 = charsize0 / !p.multi[2]
spacing = 1. * charsize0 * !d.y_ch_size / !d.y_size

;***********************************************************************
; Plot Legend

; Re-determine the legend box corner positions
ypos[1] = ypos[0] + ( 1.5 * nlabels + 0.5 ) * spacing
if n_elements( background ) ne 0 then begin
  id = where( strlen( label ) eq max( strlen( label ) ) )
  xyouts, xpos[0]+length/2., ypos[0]+spacing/2., label[id[0]], $
      /normal, charsize=charsize0, width=strwidth, font=font
  xpos[1] = xpos[0] + 2.5 * length + strwidth
  ; Clear area for legend box
  polyfill, xpos[[0,1,1,0,0]], ypos[[0,0,1,1,0]]+[0,0,1,1,0]*2.*spacing, $
      color=background, /normal
  ; Draw the legend border
  if not( keyword_set( noborder_opt ) ) then begin
    ;plots, xpos[[0,1,1,0,0]], ypos[[0,0,1,1,0]], color=!p.color, /normal
    plots, xpos[[0,1,1,0,0]], ypos[[0,0,1,1,0]], color=border_color, normal=1
  endif
endif

; Draw the bullets
for i = 0, nlabels - 1 do begin
  j = nlabels - i - 1
  psym_opt = keyword_set( psym0[j] )
  plots, xpos[0]+length+[-1,1]*length/2.*(1-psym_opt), $
      ypos[0]+(1.5*i+1.)*spacing+[0,0], /normal, color=color0[j], $
      linestyle=linestyle0[j], thick=thick0[j], psym=psym0[j]
endfor

; Label the bullets and determine the longest label
strwidth = 0.
for i = 0, nlabels - 1 do begin
  xyouts, xpos[0]+2*length, ypos[0]+(3.*i+1)*spacing/2., (reverse(label))[i], $
      width=strwidth1, /normal, charsize=charsize0, alignment=0, font=font, $
      color=border_color
  if strwidth1 gt strwidth then strwidth = strwidth1
endfor

; Legend box
if n_elements( background ) eq 0 then begin
  ; Re-determine the legend box corner positions
  xpos[1] = xpos[0] + 2.5 * length + strwidth
  ; Draw the legend border
  if not( keyword_set( noborder_opt ) ) then begin
    ;plots, xpos[[0,1,1,0,0]], ypos[[0,0,1,1,0]], color=!p.color, /normal
    plots, xpos[[0,1,1,0,0]], ypos[[0,0,1,1,0]], color=border_color, normal=1
  endif
endif
; Pass edges to output
limit = [ xpos[0], ypos[0], xpos[1], ypos[1] ]

; Title
if keyword_set( title ) then begin
  xyouts, (xpos[1]+xpos[0])/2., ypos[1]+spacing/2., title, /normal, $
          charsize=charsize0, alignment=0.5, font=font, color=border_color
endif

;***********************************************************************
; The End

return
END
