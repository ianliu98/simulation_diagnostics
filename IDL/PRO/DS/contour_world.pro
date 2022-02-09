;+
; NAME:
;    CONTOUR_WORLD
;
; PURPOSE:
;    This procedure draws a colour contour of a data field over a world map.
;
; CATEGORY:
;    Graphics
;
; CALLING SEQUENCE:
;    contour_world, Data, Lon, Lat
;
; INPUTS:
;    DATA:  The data field array in [longitude,latitude] format, or a vector of 
;        [space] format.  Of type integer or floating point.
;    LON:  A vector of the longitude coordinates of the data array, of type 
;        integer or floating point.  Of format [longitude] or [space], 
;        depending on the format of Data.
;    LAT:  A vector of the latitude coordinates of the data array, of type 
;        integer or floating point.  Of format [longitude] or [space], 
;        depending on the format of Data.
;
; KEYWORD PARAMETERS:
;    ARROWEND:  If set, then pointed end levels are plotted, indicating that 
;        the end levels extend beyond the range of the plot.  The default is 
;        rectangular end levels, like all the other levels.  A 2 element vector 
;        can also be input, where the first element gives the setting (0 or 1) 
;        for the left/bottom end and the second for the right/top end.
;    AXIS:  If set axis labels are printed.  This only works with the 
;        cylindrical projection at the moment.
;    BOLD:  If set then text is written in bold font, if set to zero then 
;        regular font is used.  The default is to be set.
;    BLOCKCELL:  If set then cells in a regular grid are filled with blocks.  
;        The default is to use IDL's interpolating contour procedure.  This is 
;        automatically set to 0 if irregularly gridded data has been input, 
;        i.e. if Data is a vector.
;    BUFFER_LIMIT:  An optional scalar number defining the buffer size to add 
;        to the region defined by LIMIT when LIMIT is being calculated by the 
;        procedure (i.e. when LIMIT=1 is input).
;    C_COLORS:  A vector of colour indices for the contoured levels.  The 
;        default may not look great with certain colour tables.
;    CENTRE:  The coordinates of the central point of the plot, of format 
;        [longitude,latitude].  If set to 1 then the coordinates are determined 
;        automatically.
;    CHARCOLOR:  The colour index of the map title and border, and the legend 
;        title and border.  The default is COLOR.
;    CHARSIZE:  The size of the text characters.  The default is 1.
;    CLIP:  If set to 0 then no clipping is performed.  See map_set for more 
;        information.
;    COASTS:  If set, coastlines, islands, and lakes are drawn.
;    COLOR:  The colour index of the continent outlines.  The default is set in 
;        !P.COLOR.
;    COUNTRIES:  If set, national borders are drawn.
;    FILL_CONTINENTS_COLOR:  Set to the colour index for filling the continents 
;        according to the FILL_CONTINENTS keyword of map_continents.
;    FONT:  An integer specifying the graphics font to use.  The default is 
;        Helvetica bold.
;    HIRES:  If set, a high resolution map of boundaries is used.
;    HORIZON:  If set then a horizon is drawn, if the projection permits.  See 
;        map_set for more information.
;    JUSTIFY_TITLE:  An optional scalar defining the maximum line length in the 
;        title string input in TITLE.  Hard returns will be inserted to ensure 
;        this maximum line length is satisfied (subject to single words not 
;        exceeding it).
;    LEGEND:  If set, the procedure plots a colour legend at the bottom of the 
;        window.  The default is no legend.  This can also be a 2-element 
;        vector containing the position ([x,y]) in normal coordinates of the 
;        lower left corner of the legend box.
;    LEVELS:  A vector of values for the contour levels.  This can be 
;        determined automatically from the data.
;    LIMIT:  A four-element vector specifying the limiting coordinates of the 
;        mapping area, in the format of [minlat,minlon,maxlat,maxlon].  If set 
;        to 1 then the limit is determined automatically.
;    NOERASE:  If set, the procedure does not erase the screen before 
;        plotting.  The default is to erase before plotting.
;    NOBORDER:  If set, then no window border line is drawn.  The default is to 
;        draw a border.
;    NOLINES:  If set, then no outlines are drawn.  The default is not set.
;    NTICKS:  The number of ticks and labels to be used in the legend.  This 
;        can be determined automatically.
;    OUTLINE_COLOR:  The colour index of the OUTLINE_DATA contour lines.  The 
;        default is COLOR.
;    OUTLINE_DATA:  An optional data array of the same size as Data.  Contour 
;        lines for OUTLINE_DATA are plotted on top of the colour-contours of 
;        Data.
;    OUTLINE_LEVELS:  Like LEVELS but corresponding to OUTLINE_DATA.  The 
;        default is LEVELS.
;    OUTLINE_POINT:  If contour lines are being drawn, then setting this 
;        keyword causes perpendicular ticks to be drawn on the contours lines.  
;        If set to -1 then the ticks point downhill on OUTLINE_DATA, if set to 
;        1 then they point uphill.  The 1 setting does not work of BLOCKCELL=0.
;    [X,Y]RANGE:  A 2-element vector containing the minimum and maximum 
;        [x,y]-coordinates to be plotted.
;    REGION:  Obsolete keyword retained for continuity.
;    RIVERS:  If set rivers are drawn.
;    THICK:  The line thickness used for drawing.
;    [X,Y]THICK:  The line thickness used for drawing the axes.  This only 
;        works if the AXIS keyword is set.
;    TICKNAME:  A vector of strings to place on the tick marks in the legend.  
;        This can be determined automatically.
;    TITLE:  A string containing the title of the plot.
;    [X,Y]TITLE:  A string containing the title of the [X,Y] axis.  This only 
;        works if the AXIS keyword is set.
;    UNITS:  A string containing the title or units subtitle for the legend.
;    WRAP:  If set, the plot connects the westernmost and eastermost points to 
;        produce a continuous plot.  The default is to do this only if the data 
;        appears global in longitude.
;    CYLINDRICAL:  If set a cylindrical projection is used.  This is the 
;        default.
;    HAMMER:  If set, the hammer projection is used.
;    MERCATOR:  If set a Mercator projection is used.
;    NORTH:  If set an orthographic projection is used, centred on the North 
;        Pole and covering the Northern Hemisphere.
;    ORTHOGRAPHIC:  If set an orthographic projection is used.
;    SATELLITE:  If set a satellite projection is used.  See map_set for more 
;        information.
;    SAT_P:  Input vector required if SATELLITE is set.  See map_set for more 
;        information.
;    SOUTH:  If set an orthographic projection is used, centred on the South 
;        Pole and covering the Southern Hemisphere.
;
; USES:
;    block_cell.pro
;    choose_levels.pro
;    contour_legend.pro
;    dimension.pro
;    geo_mean.pro
;    odd.pro
;    sign.pro
;
; PROCEDURE:
;    This procedure plots a world map using map_set.pro and map_continents.pro 
;    and then plots a colour contour on top using contour.pro.
;
; EXAMPLE:
;    Contour an increasing-value 20*10 array over a map of Earth.
;      arr = findgen(20,10)
;      loadct, 5
;      contour_world, arr
;
; MODIFICATION HISTORY:
;    Written by:  Daithi A. Stone (stoned@atm.ox.ac.uk), 2000-09-21.
;    Modified:  DAS, 2000-09-28 (level-picking stuff).
;    Modified:  DAS, 2000-11-16 (fixed polar projection bug).
;    Modified:  DAS, 2000-11-29 (added NTICKS and TICKNAME keywords).
;    Modified:  DAS, 2001-04-11 (added ability to handle irregular data).
;    Modified:  DAS, 2001-05-16 (added REGION, WRAP, LIMIT keywords).
;    Modified:  DAS, 2002-04-10 (switched STD.pro to STDDEV.pro)
;    Modified:  DAS, 2004-10-26 (improved documentation;  added AXIS, COASTS, 
;        COUNTRIES, HIRES, RIVERS, THICK, [X,Y]THICK keywords)
;    Modified:  DAS, 2006-02-14 (added BLOCKCELL keyword;  removed some WRAP 
;        value selection stuff;  made REGION option automatic and therefore 
;        obsolete)
;    Modified:  DAS, 2008-03-14 (fixed bug with LIMIT keyword wrapping 
;        longitude when BLOCKCELL set)
;    Modified:  DAS, 2009-02-12 (added NOLINES keyword)
;    Modified:  DAS, 2009-09-02 (fixed bug using coordinate ordering with 
;        irregular data;  added use of BLOCKCELL with irregular data;  ensured 
;        LEGEND option does not overwrite map plotting coordinates)
;    Modified:  DAS, 2009-11-09 (added BOLD keyword)
;    Modified:  DAS, 2009-12-08 (added ARROWEND keyword)
;    Modified:  DAS, 2009-12-09 (added fix to ensure longitudes do not span 
;        more than 360 degrees)
;    Modified:  DAS, 2009-12-15 (fixed bug in dealing with values outside the 
;        range of LEVELS)
;    Modified:  DAS, 2009-12-23 (added OUTLINE contour feature)
;    Modified:  DAS, 2010-01-06 (added OUTLINE_POINT keyword;  edited format;  
;        fixed OUTLINE_COLOR bug)
;    Modified:  DAS, 2010-05-02 (fixed bug in OUTLINE with non-monotonic 
;        coordinates)
;    Modified:  DAS, 2010-05-21 (fixed bug with shifting !p.multi values after 
;        call to contour_legend.pro;  added modification to CHARSIZE in 
;        contour_legend.pro if !p.multi is used)
;    Modified:  DAS, 2010-09-16 (added FONT keyword)
;    Modified:  DAS, 2010-11-19 (added CLIP, HORIZON, NOBORDER, SATELLITE, 
;        SAT_P keywords;  added OUTLINE_DATA capability when BLOCKCELL=1 and 
;        DATA is irregularly gridded)
;    Modified:  DAS, 2011-03-16 (fixed bug with legend tick labels when 
;        !p.multi=0)
;    Modified:  DAS, 2012-12-07 (added the FILL_CONTINENTS_COLOR keyword)
;    Modified:  DAS, 2015-01-15 (added the HAMMER projection keyword;  modified
;        clipping behaviour when BLOCKCELL=1 and LIMIT is input;  improved
;        performance when LIMIT is of size 8)
;    Modified:  DAS, 2015-05-15 (added options for automatic CENTRE and LIMIT 
;        selection;  added BUFFER_LIMIT keyword;  added JUSTIFY_TITLE keyword;  
;        set CHARCOLOR default to COLOR;  modified formating)
;    Modified:  DAS, 2016-10-10 (Outsourced regular BLOCKCELL plotting to 
;        block_cell.pro;  fixed bug in ORTHOSCOPIC, NORTH, and SOUTH options;  
;        made automatic estimation of LIMIT more robust to numerical 
;        imprecision)
;    Modified:  DAS, 2018-06-07 (Modified the default LON vector)
;    Modified:  DAS, 2018-08-24 (Added capability to run BLOCK_CELL option for 
;        a cylindrical projection under GDL)
;    Modified:  DAS, 2018-10-08 (Added plotting fix when latitude limits are on 
;        the poles)
;-

;***********************************************************************

PRO CONTOUR_WORLD, $
    DATA_0, LON_0, LAT_0, $
    ARROWEND=arrowend_opt, $
    AXIS=axis_opt, $
    BOLD=bold_opt, FONT=font, $
    BLOCKCELL=blockcell_opt, $
    C_COLORS=c_colors, $
    CENTRE=centre, $
    CHARCOLOR=charcolor, CHARSIZE=charsize, $
    CLIP=clip_opt, $
    COASTS=coasts_opt, COUNTRIES=countries_opt, HIRES=hires_opt, $
      RIVERS=rivers_opt, FILL_CONTINENTS_COLOR=fill_continents_color, $
    COLOR=color, $
    CYLINDRICAL=cylindrical_opt, HAMMER=hammer_opt, MERCATOR=mercator_opt, $
      ORTHOGRAPHIC=orthographic_opt, NORTH=north_opt, SOUTH=south_opt, $
      SATELLITE=satellite_opt, SAT_P=sat_p, $
    HORIZON=horizon_opt, $
    LEGEND=legend, $
    LEVELS=levels, NLEVELS=nlevels, $
    LIMIT=limit, BUFFER_LIMIT=buffer_limit, $
    XMARGIN=xmargin, YMARGIN=ymargin, $
    NOBORDER=noborder_opt, $
    NOERASE=noerase, $
    NOLINES=nolines_opt, $
    OUTLINE_COLOR=outline_color, OUTLINE_DATA=outline_data_0, $
      OUTLINE_LEVELS=outline_levels, OUTLINE_POINT=outline_point, $
    NTICKS=nticks, TICKNAME=tickname, $
    THICK=thick, xthick=xthick, ythick=ythick, $
    TITLE=title, XTITLE=xtitle, YTITLE=ytitle, JUSTIFY_TITLE=justify_title, $
    UNITS=units, $
    WRAP=wrap_opt, $
    REGION=region_opt

;***********************************************************************
; Constants and Options

; Physical and mathematical constants
n_degree = 360
degrad = !pi / n_degree * 2.

;; Equivalent to zero threshold
;zthresh = 0.00001

; Legend option
legend_opt = keyword_set( legend )

; Wrapping option
wrap_opt = keyword_set( wrap_opt )

; Option to print axes
axis_opt = keyword_set( axis_opt )

; Map projections
; Cylindrical projection
cylindrical_opt = keyword_set( cylindrical_opt )
; Hammer projection
hammer_opt = keyword_set( hammer_opt )
; Mercator projection
mercator_opt = keyword_set( mercator_opt )
; Orthoscopic projection
orthographic_opt = keyword_set( orthographic_opt )
north_opt = keyword_set( north_opt )
south_opt = keyword_set( south_opt )
if north_opt or south_opt then begin
  orthographic_opt = 1
endif else begin
  orthographic_opt = 0
endelse
; The satellite projection
satellite_opt = keyword_set( satellite_opt )
if satellite_opt eq 0 then satellite_opt = keyword_set( sat_p )

; Irregular gridding
if dimension( data_0 ) eq 1 then begin
  irregular_opt = 1
endif else begin
  irregular_opt = 0
endelse

; Option to plot block cells rather than smooth contours
if keyword_set( blockcell_opt ) then begin
  blockcell_opt = 1
endif else begin
  blockcell_opt = 0
endelse
;; This cannot work with an irregular grid
;if irregular_opt eq 1 then blockcell_opt = 0

; Arrowed end option
if keyword_set( arrowend_opt ) then begin
  if n_elements( arrowend_opt ) eq 1 then begin
    arrowend_opt = [ 0, 0 ] + arrowend_opt
  endif
endif else begin
  arrowend_opt = [ 0, 0 ]
endelse

;***********************************************************************
; Default Settings

; Printing settings
if !d.name eq 'PS' then begin
  if not( keyword_set( thick ) ) then thick = 3
  if not( keyword_set( xthick ) ) then xthick = thick
  if not( keyword_set( ythick ) ) then ythick = xthick
  if n_elements( bold_opt ) eq 0 then bold_opt = 1
  if n_elements( font ) eq 0 then begin
    device, helvetica=1, bold=bold_opt
    font = 0
  endif
  !p.font = font
  unset_font_opt = 1
endif

; Legend settings
if legend_opt then begin
  if not( keyword_set( ymargin ) ) then ymargin = [6,3]
  if n_elements( legend ) ne 2 then legend = [0.3,0.09]
endif

; Default cylindrical map projection
if ( mercator_opt eq 0 ) and ( orthographic_opt eq 0 ) $
    and ( satellite_opt eq 0 ) and ( hammer_opt eq 0 ) then begin
  cylindrical_opt = 1
endif

; Copy input to temporary plotting arrays
data = data_0
if keyword_set( lon_0 ) then lon = lon_0
if keyword_set( lat_0 ) then lat = lat_0
n_limit = n_elements( limit )
if max( n_limit eq [ 0, 1, 4, 8 ] ) eq 0 then stop
if keyword_set( outline_data_0 ) then outline_data = outline_data_0
; Ensure positive and increasing longitude
if ( irregular_opt eq 0 ) and keyword_set( lon ) then begin
  ; Note this may be reversed later if limit is provided
  id = where( lon lt 0, nid )
  if nid ne 0 then lon[id] = lon[id] + n_degree
  ; We redo the ordering of LON later in this section
  id = sort( lon )
  lon = lon[id]
  data = data[id,*]
  if keyword_set( outline_data ) then outline_data = outline_data[id,*]
endif
; Ensure increasing latitude
if ( irregular_opt eq 0 ) and keyword_set( lat ) then begin
  id = sort( lat )
  lat = lat[id]
  data = data[*,id]
  if keyword_set( outline_data ) then outline_data = outline_data[*,id]
endif

; Longitude and latitude dimensions
if keyword_set( lon ) then begin
  n_lon = n_elements( lon )
endif else begin
  n_lon = n_elements( data[*,0] )
  lon = ( findgen( n_lon ) + 0.5 ) / n_lon * n_degree
endelse
if keyword_set( lat ) then begin
  n_lat = n_elements( lat )
endif else begin
  n_lat = n_elements( data[0,*] )
  lat = ( findgen( n_lat ) - ( n_lat - 1. ) / 2. ) / n_lat * 2. * n_degree / 4.
endelse

; Find convenient longitude order
if irregular_opt eq 0 then begin
  ; Find a gap in longitude
  fd = abs( lon[1:n_lon-1] - lon[0:n_lon-2] )
  fd = round( fd / min( fd ) )
  id = ( where( fd gt 1, nid ) )[0]
  ; Re-organise longitude if we have found a gap
  if nid ne 0 then begin
    id1 = indgen( id + 1 )
    id2 = indgen( n_lon - 1 - id ) + id + 2
    lon[id2] = lon[id2] - n_degree
    lon = lon[[id2,id1]]
    data = data[[id2,id1],*]
    if keyword_set( outline_data ) then outline_data = outline_data[[id2,id1],*]
  endif
endif
; Set the default mapping limits
if n_limit eq 0 then begin
  if irregular_opt eq 1 then begin
    limit = [ min( lat ), min( lon ), max( lat ), max( lon ) ]
  endif else begin
    limit = [ lat[0] - ( lat[1] - lat[0] ) / 2., $
        lon[0] - ( lon[1] - lon[0] ) / 2., $
        lat[n_lat-1] + ( lat[n_lat-1] - lat[n_lat-2] ) / 2., $
        lon[n_lon-1] + ( lon[n_lon-1] - lon[n_lon-2] ) / 2. ]
    if limit[3] - limit[1] gt 360. then begin
      limit[3] = limit[3] - 0.001 * ( lon[1] - lon[0] )
      if limit[3] - limit[1] gt 360. then begin
        limit[1] = limit[1] + 0.001 * ( lon[1] - lon[0] )
      endif
      if limit[3] - limit[1] gt 360. then stop
    endif
  endelse
  n_limit = 4
; Or if we are to determine the limits
endif else if n_limit eq 1 then begin
  ; The buffer to retain around the region, in degrees
  if n_elements( buffer_limit ) eq 0 then buffer_limit = 10
  ; Determine the coordinates of the defined points within the map
  id_points = where( finite( data ) eq 1, n_id_points )
  if n_id_points eq 0 then stop
  temp_points = fltarr( 2, n_id_points )
  if irregular_opt eq 1 then begin
    temp_points[0,*] = lon[id_points]
    temp_points[1,*] = lat[id_points]
  endif else begin
    id = id_points mod n_lon
    temp_points[0,*] = lon[id]
    id = id_points / n_lon
    temp_points[1,*] = lat[id]
  endelse
  ; Determine the limit and centre
  if n_id_points eq 1 then begin
    temp_points = [ [ temp_points ], [ temp_points ] ]
    temp_centre = geo_mean( temp_points, median=1, limit=limit, $
        buffer_limit=buffer_limit )
  endif else begin
    temp_centre = geo_mean( temp_points, minimise_maxdist=1, limit=limit, $
        buffer_limit=buffer_limit )
  endelse
  n_limit = n_elements( limit )
  ; Take calculated centre point if requested
  if n_elements( centre ) eq 1 then centre = temporary( temp_centre )
  ; Clear memory
  id_points = 0
  temp_points = 0
; Otherwise ensure increasing longitude limits
endif else if n_limit eq 4 then begin
  if limit[3] lt limit[1] then limit[1] = limit[1] - n_degree
endif
; Ensure legal ranges
if ( cylindrical_opt + mercator_opt eq 1 ) and ( n_limit eq 4 ) then begin
  temp = min( lat[1:n_lat-1] - lat[0:n_lat-2] )
  if temp eq 0 then begin
    temp = 0.001
  endif else begin
    temp = min( [ 0.001, temp ] )
  endelse
  if limit[0] lt -n_degree / 4. then limit[0] = -n_degree / 4. + temp
  if limit[2] gt n_degree / 4. then limit[2] = n_degree / 4. - temp
  if limit[1] gt limit[3] then limit[1] = limit[3]
  if limit[3] - n_degree - limit[1] gt 0 then limit[3] = limit[1] + n_degree
endif

; Ensure that the longitude values are within the mapping limits
if ( dimension( data ) eq 2 ) and ( n_limit eq 4 ) then begin
  id_limit = [ 1, 3 ]
  ;if n_limit eq 8 then id_limit = [ id_limit, 5, 7 ]
  id = where( lon - n_degree ge min( limit[id_limit] ), n_id_1 )
  if n_id_1 gt 0 then lon[id] = lon[id] - n_degree
  id = where( lon + n_degree le max( limit[id_limit] ), n_id_2 )
  if n_id_2 gt 0 then lon[id] = lon[id] + n_degree
  id = sort( lon )
  if max( abs( id[1:n_lon-1] - id[0:n_lon-2] ) ) gt 1 then begin
    lon = lon[id]
    data = data[id,*]
    if keyword_set( outline_data ) then outline_data = outline_data[id,*]
  endif
endif

;; Choice to wrap the data around the map sides
;if not( wrap_opt ) then begin
;  fd = lon[1:n_lon-1] - lon[0:n_lon-2]
;  if stddev( fd[0:n_lon-2], nan=1 ) lt zthresh then begin
;    if round( fd[0] - fd[n_lon-1] ) eq n_degree then wrap_opt = 1
;  endif
;endif

; Central point of the plot
if keyword_set( centre ) then begin
  centlon = centre[0]
  centlat = centre[1]
endif else if keyword_set( limit ) then begin
  if n_limit eq 4 then begin
    centlon = ( limit[1] + limit[3] ) / 2.
    if south_opt then begin
      centlat = -90.
    endif else if north_opt then begin
      centlat = 90.
    endif else begin
      centlat = ( limit[0] + limit[2] ) / 2.
    endelse
  endif else begin
    centlon = ( limit[1] + limit[5] ) / 2.
    if south_opt then begin
      centlat = -90.
    endif else if north_opt then begin
      centlat = 90.
    endif else begin
      centlat = ( limit[2] + limit[6] ) / 2.
    endelse
  endelse
  ;centlat = 0.
endif else begin
  stop
endelse

; Contour levels
if not( keyword_set( levels ) ) then begin
  if not( keyword_set( nlevels ) ) then nlevels = 29
  levels = choose_levels( data )
  levels = levels[0] + findgen( nlevels+1 ) / nlevels $
      * ( levels[n_elements(levels)-1] - levels[0] )
endif
nlevels = n_elements( levels )
; Ensure proper behaviour at end levels
levels_contour = levels
if keyword_set( arrowend_opt ) then begin
  if arrowend_opt[0] eq 1 then begin
    temp = min( data, nan=1 )
    if temp lt levels_contour[0] then begin
      if temp lt 0 then begin
        temp = temp * 1.1
      endif else begin
        temp = temp / 1.1
      endelse
      levels_contour[0] = temp
    endif
  endif
  if arrowend_opt[1] eq 1 then begin
    temp = max( data, nan=1 )
    if temp gt levels_contour[nlevels-1] then begin
      if temp lt 0 then begin
        temp = temp / 1.1
      endif else begin
        temp = temp * 1.1
      endelse
      levels_contour[nlevels-1] = temp
    endif
  endif
endif

; Colour scale
if not( keyword_set( c_colors ) ) then c_colors = indgen( nlevels + 1 ) + 2

; The defaults for the option of plotting contour lines
if keyword_set( outline_data ) then begin
  if not( keyword_set( outline_color ) ) then begin
    if keyword_set( color ) then begin
      outline_color = color
    endif else begin
      outline_color = !p.color
    endelse
  endif
  if not( keyword_set( outline_levels ) ) then outline_levels = levels
  if keyword_set( outline_point ) then downhill_opt = -1 * outline_point
endif

; The default character color
if ( n_elements( charcolor ) eq 0 ) and ( n_elements( color ) eq 1 ) then begin
  charcolor = color
endif

; Option to write title in justified text format
if keyword_set( justify_title ) and keyword_set( title ) then begin
  ; Split the title into its constituent words
  temp = strsplit( title, ' ', extract=1, count=n_temp )
  temp_len = strlen( temp ) + 1
  temp_len[0] = temp_len[0] - 1
  temp_len = total( temp_len, cumulative=1, integer=1 )
  ; Iterate through lines
  check = 0
  id_line = 0
  while check eq 0 do begin
    ; Identify words to retain in this line
    id = max( where( temp_len[id_line:n_temp-1] lt justify_title ) )
    ; If the first word is longer than justify_title then take it anyway
    if id eq -1 then id = 0
    ; Add a hard return if this is not the last line
    if id_line + id ne n_temp - 1 then begin
      temp[id_line+id] = temp[id_line+id] + '!C'
    endif
    ; Move on to the next line
    id_line = id_line + id + 1
    temp_len = temp_len - temp_len[id_line-1]
    ; Exit loop if this is the last line
    if id_line gt n_temp - 1 then check = 1
  endwhile
  ; Re-build the title
  title_use = string_from_vector( temp, spacer=' ', nospace=1 )
  temp = 0
; Otherwise just copy the titlee
endif else if keyword_set( title ) then begin
  title_use = title
endif

;***********************************************************************
; Plot Map

; Record !p.multi settings (MAP_SET messes them up)
pmulti = !p.multi

; Set up map
gdl_opt = !prompt eq 'GDL> '
if ( gdl_opt eq 1 ) and ( cylindrical_opt eq 1 ) and ( n_limit eq 4 ) then begin
  xstyle = 1
  ystyle = 1
  if keyword_set( noborder_opt ) then begin
    xstyle = 5
    ystyle = 5
  endif
  if n_elements( clip_opt ) eq 0 then clip_opt = 1
  plot, limit[[1,3]], limit[[0,2]], nodata=1, isotropic=1, noerase=noerase, $
      title=title_use, charsize=charsize, xmargin=xmargin, ymargin=ymargin, $
      xstyle=xstyle, ystyle=ystyle, clip=clip_opt, color=charcolor
endif else begin
  map_set, centlat, centlon, 0, cylindrical=cylindrical_opt, $
      mercator=mercator_opt, orthographic=orthographic_opt, $
      satellite=satellite_opt, sat_p=sat_p, color=charcolor, $
      isotropic=1, noerase=noerase, title=title_use, charsize=charsize, $
      xmargin=xmargin, ymargin=ymargin, limit=limit, horizon=horizon_opt, $
      noborder=noborder_opt, clip=clip_opt, hammer=hammer_opt
endelse

; Draw land coverage
if n_elements( fill_continents_color ) eq 1 then begin
  map_continents, hires=hires_opt, color=fill_continents_color, $
      fill_continents=1
endif

; Determine the longitude and latitude extent
check = 1
if clip_opt eq 1 then begin
  if ( blockcell_opt eq 1 ) and ( n_limit ge 4 ) then check = 0
endif else begin
  check = 0
endelse
if check eq 0 then begin
  limit_lon_min = centlon - n_degree / 2.
  limit_lon_max = centlon + n_degree / 2.
  limit_lat_min = -n_degree / 4.
  limit_lat_max = n_degree / 4.
endif else begin
  if n_limit eq 4 then begin
    limit_lon_min = limit[1]
    limit_lon_max = limit[3]
    limit_lat_min = limit[0]
    limit_lat_max = limit[2]
  endif else begin
    limit_lon_min = min( limit[[1,3,5,7]] )
    limit_lon_max = max( limit[[1,3,5,7]] )
    limit_lat_min = min( limit[[0,2,4,6]] )
    limit_lat_max = max( limit[[0,2,4,6]] )
  endelse
endelse

; Contour over map.
; If we have irregularly gridded data for a smooth contour style of plot
if ( irregular_opt eq 1 ) and ( blockcell_opt eq 0 ) then begin
  contour, data, lon, lat, levels=levels_contour, c_colors=c_colors, $
      cell_fill=1, overplot=1, irregular=1
; If we have irregularly gridded data for a block cell style of plot
endif else if ( irregular_opt eq 1 ) and ( blockcell_opt eq 1 ) then begin
  ; Set the polygon size
  n_vertex = 4
    ; Take one pass for the data and one for the outline contours
    for i_pass = 0, keyword_set( outline_data ) do begin
    ; Iterate through elements in the data
    for i = 0, n_lon - 1 do begin
      ; Find n_vertex neighbouring points to this one
      temp = ( lon - lon[i] ) ^ 2 + ( lat - lat[i] ) ^ 2
      id_data = sort( temp )
      id_data = id_data[1:n_vertex]
      lon_near = lon[id_data]
      lat_near = lat[id_data]
      if i_pass eq 1 then outline_data_near = outline_data[id_data]
      ; Sort vertices in counterclockwise order from the smallest angle
      ang = atan( lat_near - lat[i], lon_near - lon[i] )
      id = sort( ang )
      lon_near = lon_near[id]
      lat_near = lat_near[id]
      ang = ang[id]
      if i_pass eq 1 then outline_data_near = outline_data_near[id]
      ; If the largest angle is greater than pi
      temp = [ ang[1:n_vertex-1], ang[0] + 2. * !pi ] - ang
      if max( temp ) gt 0.9 * !pi then begin
        ; Retain the closest points
        temp_1 = ( lon_near - lon[i] ) ^ 2 + ( lat_near - lat[i] ) ^ 2
        id = sort( temp_1 )
        ; If more than 1.5*pi is missing then retain only half the points
        if max( temp ) gt 1.4 * !pi then begin
          id = id[0:(n_vertex+1)/2-1]
        ; Otherwise retain 3/4 of the points
        endif else begin
          id = id[0:3*n_vertex/4-1]
        endelse
        lon_near = lon_near[id]
        lat_near = lat_near[id]
        if i_pass eq 1 then outline_data_near = outline_data_near[id]
        ; Create fake points reflected across our point
        lon_near = [ lon_near, lon[i] + ( lon[i] - lon_near ) ]
        lat_near = [ lat_near, lat[i] + ( lat[i] - lat_near ) ]
        if i_pass eq 1 then begin
          outline_data_near = [ outline_data_near, outline_data_near ]
        endif
        ; Restrict to n_vertex points
        n_temp = n_elements( lon_near )
        if n_temp gt n_vertex then begin
          ang = atan( lat_near - lat[i], lon_near - lon[i] )
          id = sort( ang )
          lon_near = lon_near[id]
          lat_near = lat_near[id]
          ang = ang[id]
          if i_pass eq 1 then outline_data_near = outline_data_near[id]
          temp = [ ang[1:n_temp-1], ang[0] + 2. * !pi ] - ang
          ; Remove redundant points
          id = where( temp gt !pi / 10., n_id )
          lon_near = lon_near[id]
          lat_near = lat_near[id]
          ang = ang[id]
          if i_pass eq 1 then outline_data_near = outline_data_near[id]
          temp = temp[id]
          if n_id gt n_vertex then begin
            id = sort( temp )
            id = id[0:n_vertex-1]
            lon_near = lon_near[id]
            lat_near = lat_near[id]
            if i_pass eq 1 then outline_data_near = outline_data_near[id]
          endif
        endif
      endif
      ; Determine vertices of polygon
      lon_vertex = fltarr( n_vertex )
      lat_vertex = fltarr( n_vertex )
      lon_vertex[0] = ( lon_near[n_vertex-1] + lon_near[0] ) / 2.
      lat_vertex[0] = ( lat_near[n_vertex-1] + lat_near[0] ) / 2.
      for j = 1, n_vertex - 1 do begin
        lon_vertex[j] = ( lon_near[j-1] + lon_near[j] ) / 2.
        lat_vertex[j] = ( lat_near[j-1] + lat_near[j] ) / 2.
      endfor
      ; If we are drawing filled contours
      if i_pass eq 0 then begin
        ; Determine colour
        id = max( where( levels_contour - data[i] lt 0 ) )
        if id eq -1 then id = 0
        ; Plot cell
        polyfill, lon_vertex, lat_vertex, color=c_colors[id]
      ; If we are drawing outline contours
      endif else begin
        ; Determine the contour level of this cell
        id_level = max( where( outline_levels - outline_data[i] lt 0 ) )
        ; Iterate through neighbouring cells
        for i_near = 0, n_vertex - 1 do begin
          ; Determine the contour level of the cell to the east
          id_level_1 = max( where( $
              outline_levels - outline_data_near[i_near] lt 0 ) )
          ; If the contour levels are different then plot a contour line
          if id_level_1 ne id_level then begin
            if i_near eq n_vertex - 1 then begin
              temp_lon = lon_vertex[[n_vertex-1,0]]
              temp_lat = lat_vertex[[n_vertex-1,0]]
            endif else begin
              temp_lon = lon_vertex[[i_near,i_near+1]]
              temp_lat = lat_vertex[[i_near,i_near+1]]
            endelse
            plots, temp_lon, temp_lat, color=outline_color, thick=thick
            ; Plot a down-/up-hill tick if requested
            if keyword_set( outline_point ) then begin
              d_lon_tick = [ 0, ( lon_near[i_near] - mean( lon_vertex ) ) / 4. ]
              d_lat_tick = [ 0, ( lat_near[i_near] - mean( lat_vertex ) ) / 4. ]
              lon_tick = [ 0, 0 ] + mean( temp_lon )
              lat_tick = [ 0, 0 ] + mean( temp_lat )
              temp = outline_point * sign( id_level_1 - id_level )
              oplot, lon_tick+temp*d_lon_tick, lat_tick+temp*d_lat_tick, $
                  color=outline_color, thick=thick
            endif
          endif
        endfor
      endelse
    endfor
  endfor
; If we have regularly gridded data for a block cell style of plot
endif else if blockcell_opt eq 1 then begin
  ; Plot filled contours
  if n_limit eq 4 then begin
    temp_limit = limit[[1,0,3,2]]
  endif else if n_limit eq 8 then begin
    temp_limit = limit[[1,0,3,2,5,4,7,6]]
  endif
  block_cell, data, lon, lat, levels=levels_contour, c_colors=c_colors, $
      color=color, limit=temp_limit
  ; Plot outline contours
  if keyword_set( outline_data ) then begin
    block_cell, outline_data, lon, lat, outline=1, levels=outline_levels, $
        color=outline_color, thick=thick, limit=temp_limit, $
        point_ticks=outline_point
  endif
; If we have regularly gridded data for a smooth contour style of plot
endif else begin
  if wrap_opt then begin
    ; Wrapped global data
    contour, [data,data[0,*]], [lon,lon[0]], lat, levels=levels_contour, $
        c_colors=c_colors, cell_fill=1, overplot=1
    ; Draw outline contours if requested
    if keyword_set( outline_data ) then begin
      contour, [outline_data,outline_data[0,*]], [lon,lon[0]], lat, $
          levels=outline_levels, overplot=1, color=outline_color, thick=thick, $
          downhill=downhill_opt
    endif
  endif else begin
oplot,[0,1],[0,1]
    ; Regional data
    contour, data, lon, lat, levels=levels_contour, c_colors=c_colors, $
        cell_fill=1, overplot=1
    ; Draw outline contours if requested
    if keyword_set( outline_data ) then begin
      contour, outline_data, lon, lat, levels=outline_levels, overplot=1, $
          color=outline_color, thick=thick, downhill=downhill_opt
    endif
  endelse
endelse

; Draw coastlines
if not( keyword_set( nolines_opt ) ) then begin
  map_continents, hires=hires_opt, coasts=coasts_opt, rivers=rivers_opt, $
      countries=countries_opt, color=color, thick=thick
endif

; Draw axes.
; This is simple only if a cylindrical projection is used
if cylindrical_opt and axis_opt then begin
  xrange = [ min( lon ), max( lon ) ] $
      + 0.01 * [ -1., 1. ] * ( max( lon ) - min( lon ) )
  yrange = [ min( lat ), max( lat ) ] $
      + 0.01 * [ -1., 1. ] * ( max( lat ) - min( lat ) )
  axis, xrange[0], 0, yaxis=0, xmargin=xmargin, ymargin=ymargin, $
      ytitle=ytitle, ythick=ythick, charsize=charsize, $
      yrange=yrange, ystyle=1, noerase=1, font=font
  axis, xrange[1], 0, yaxis=1, xmargin=xmargin, ymargin=ymargin, $
      ythick=ythick, charsize=charsize, yrange=yrange, ystyle=1, noerase=1, $
      ytickname=strarr(30)+' ', font=font
  axis, 0, yrange[0], xaxis=0, xmargin=xmargin, ymargin=ymargin, $
      xtitle=xtitle, xthick=xthick, charsize=charsize, $
      xrange=xrange, xstyle=1, noerase=1, font=font
  axis, 0, yrange[1], xaxis=1, xmargin=xmargin, ymargin=ymargin, $
      xthick=xthick, charsize=charsize, xrange=xrange, xstyle=1, noerase=1, $
      xtickname=strarr(30)+' ', font=font
endif

; Restore !p.multi settings
!p.multi = pmulti

;***********************************************************************
; Plot Legend

if legend_opt then begin
  ; Set legend charsize
  ;if n_elements( pmulti ) ge 3 then begin
  if max( pmulti ) ge 2 then begin
    temp_charsize = 1. * charsize / ( min( pmulti[1:2] ) )
  endif else begin
    temp_charsize = charsize
  endelse
  ; Plot legend
  contour_legend, legend[0]+[0.,0.4], legend[1]+[0.,0.03], $
      levels=levels, c_colors=c_colors, color=charcolor, $
      charsize=temp_charsize, normal=1, horizontal=1, subtitle=units, $
      nticks=nticks, tickname=tickname, arrowend=arrowend_opt, font=font
  ; Ensure that plotting the legend has not forced us to a new sub-window
  ;if pmulti[1] + pmulti[2] gt 2 then begin
  ;  if !p.multi[0] eq 0 then !p.multi[0] = !p.multi[1] * !p.multi[2]
  ;  !p.multi[0] = !p.multi[0] - 1
  ;endif
  !p.multi = pmulti
  ; Ensure that mapping coordinates are retained in the plotting window
  map_set, centlat, centlon, 0, cylindrical=cylindrical_opt, $
      mercator=mercator_opt, orthographic=orthographic_opt, $
      satellite=satellite_opt, sat_p=sat_p, color=charcolor, $
      isotropic=1, noerase=1, title=title_use, charsize=charsize, $
      xmargin=xmargin, ymargin=ymargin, limit=limit, horizon=horizon_opt, $
      noborder=noborder_opt, clip=clip_opt, hammer=hammer_opt
endif

;***********************************************************************
; Setdown

; Undo printing settings
if keyword_set( unset_font_opt ) then !p.font = -1

;***********************************************************************
; The End

return
END
