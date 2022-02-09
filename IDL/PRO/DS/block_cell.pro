;+
; NAME:
;    BLOCK_CELL
;
; PURPOSE:
;    This procedure draws a boxed contour plot, where boxes are assigned a 
;    level instead of calculating smooth boundaries for the levels.
;
; CATEGORY:
;    Graphics
;
; CALLING SEQUENCE:
;    block_cell, data, xvec, yvec
;
; INPUTS:
;    DATA:  The data field array of size [N_XVEC,N_YVEC], or a vector of length 
;        N_SPACE.  Of type integer or floating point.
;    XVEC:  A vector or array of the X-axis coordinates of the data array, of 
;        type integer or floating point.  Of length N_XVEC or N_SPACE, 
;        depending on the format of Data.
;    YVEC:  A vector or array of the Y-axis coordinates of the data array, of 
;        type integer or floating point.  Of length N_YVEC or N_SPACE, 
;        depending on the format of Data.
;
; KEYWORD PARAMETERS:
;    BORDER:  If RETURN_BORDER is set, then this returns a 2-dimensional 
;        floating point number array containing the coordinates of the borders 
;        of the shapes.  The array is of the form [3,N_POINT], where N_POINT is 
;        the total number of points in all of the shapes.  Values in [0,*] and 
;        [1,*] are the coordinates in the X-axis and Y-axis respectively.  
;        Values in [2,*] are the indices of the shapes, running from 0 to 
;        N_SHAPE-1 where N_SHAPE is the total number of shapes.  If 
;        RETURN_BORDER is not set then nothing is returned.
;    C_COLORS:  A vector of colour indices for the contoured levels.  The 
;        default may not look great with certain colour tables.
;    CLIP:  If set to 0 then no clipping is performed.
;    COLOR:  The colour index of the various lines.  The default is set in 
;        !P.COLOR.
;    LEVELS:  A vector of values for the contour levels.  This can be 
;        determined automatically from the data.
;    LIMIT:  A four-element or eight-element vector specifying the limiting 
;        values of the [x,y] dimensions.  For a four-element vector, the format 
;        is [xmin,ymin,xmax,ymax].  For an eight-element vector the format is 
;        [x-left,y-left,x-top,y-top,x-right,y-right,x-bottom,y-bottom], for 
;        points on the left, top, right, and bottom edges of the plot, 
;        respectively.  If not set then the values of [X,Y]RANGE are used, or 
;        are taken as the minimum/maximum value minus/plus half the distance to 
;        the second-smallest/largest value.
;    NO_PLOT:  If set then no plotting is performed.  If RETURN_BORDER is set 
;        then the procedure can serve as a shape calculation tool.  The default 
;        is for plotting.
;    OUTLINE:  If set then contour lines for DATA are plotted instead of the 
;        filled colour-contours.
;    POINT_TICKS:  If contour lines are being drawn, then setting this 
;        keyword causes perpendicular ticks to be drawn on the contours lines.  
;        If set to -1 then the ticks point downhill on DATA, if set to 1 then 
;        they point uphill.
;    [X,Y]RANGE:  A 2-element vector containing the minimum and maximum 
;        [x,y]-coordinates to be plotted.
;    RETURN_BORDER:  If set then output is returned in BORDER.  The default 
;        is for no output to be returned in BORDER.
;    THICK:  The line thickness used for drawing.
;
; OUTPUTS:
;    BORDER
;
; USES:
;    choose_levels.pro
;    dimension.pro
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
;    Written by:  Daithi Stone (dstone@lbl.gov), 2012-08-19 (code extracted 
;        from contour_world.pro).
;    Modified:  DAS, 2016-10-20 (Extracted from contour_world.pro)
;    Modified:  DAS, 2017-01-24 (Added BORDER output.  Added NO_PLOT, 
;        RETURN_BORDER options)
;    Modified:  DAS, 2017-06-08 (Corrected error in calculation of minimum 
;        longitude deltas;  sorted plotting from larger shapes to smaller 
;        shapes in order to avoid overplotting)
;    Modified:  DAS, 2018-01-06 (Corrected bug in which no amalgamation of 
;        shapes was done in the x-dimension... not sure when that one appeared)
;    Modified:  DAS, 2018-05-02 (Fixed issue in which NaN-defined shapes were 
;        being lumped together and thus only the first shape was being plotted)
;    Modified:  DAS, 2018-09-25 (Fixed automatic selection of contour levels;  
;        Fixed errors in calculating node distances and limiting numbers of 
;        nearby nodes with irregular gridss)
;
; TASKS:
;    Enact CLIP option
;    Enact POINT_TICKS option
;    Implement efficient method for irregular grids
;-

;***********************************************************************

PRO BLOCK_CELL, $
    DATA, XVEC, YVEC, $
    C_COLORS=c_colors, $
    CLIP=clip_opt, $
    COLOR=color, $
    LEVELS=levels, $
    LIMIT=limit, $
    OUTLINE=outline_opt, $
    POINT_TICKS=point_ticks, $
    XRANGE=xrange, YRANGE=yrange, $
    THICK=thick, $
    NO_PLOT=no_plot_opt, $
    RETURN_BORDER=return_border_opt, $
    BORDER=border_out

;***********************************************************************
; Constants and Options

; Option of no plot output
no_plot_opt = keyword_set( no_plot_opt )

; Option for line contours
outline_opt = keyword_set( outline_opt )

; Option to return border coorinates for the shapes
return_border_opt = keyword_set( return_border_opt )

; The default number of vertices for boxes when given irregularly gridded data
if outline_opt eq 1 then begin
  n_vertex_0 = 4
endif else begin
  n_vertex_0 = 4
endelse

; Determine if the data is gridded or irregular
if dimension( data ) eq 1 then begin
  irregular_opt = 1
endif else begin
  irregular_opt = 0
endelse

; Determine the input array dimensions
n_xvec = n_elements( xvec )
n_yvec = n_elements( yvec )
n_data = n_elements( data )
if irregular_opt eq 1 then begin
  if n_xvec ne n_data then stop
  if n_yvec ne n_data then stop
endif else begin
  if n_data ne n_xvec * n_yvec then stop
endelse

; Default contour levels
if not( keyword_set( levels ) ) then begin
  ;if not( keyword_set( n_levels ) ) then n_levels = 29
  levels = choose_levels( data )
  n_levels = n_elements( levels )
  ;levels = levels[0] + findgen( n_levels + 1 ) / n_levels $
  ;    * ( levels[n_levels-1] - levels[0] )
endif else begin
  n_levels = n_elements( levels )
endelse

; Colour scale
if not( keyword_set( c_colors ) ) and ( outline_opt eq 0 ) then begin
  c_colors = indgen( n_levels + 1 ) + 2
endif

;***********************************************************************
; Prepare data for analysis

; Copy input arrays
data_use = data
xvec_use = xvec
yvec_use = yvec

; Ensure monotonically increasing dimensions for regularly gridded data
if irregular_opt eq 0 then begin
  id = sort( xvec_use )
  if max( id[1:n_xvec-1] - id[0:n_xvec-2] ) gt 1 then begin
    xvec_use = xvec_use[id]
    data_use = data_use[id,*]
  endif
  id = sort( yvec_use )
  if max( id[1:n_yvec-1] - id[0:n_yvec-2] ) gt 1 then begin
    yvec_use = yvec_use[id]
    data_use = data_use[*,id]
  endif
  ; Get the minimum deltas
  d_xvec_use_min = min( xvec_use[1:n_xvec-1] - xvec_use[0:n_xvec-2] )
  d_yvec_use_min = min( yvec_use[1:n_yvec-1] - yvec_use[0:n_yvec-2] )
endif

; Set the default plotting limits
if not( keyword_set( limit ) ) then begin
  limit = fltarr( 4 )
  if keyword_set( xrange ) then begin
    limit[[0,2]] = xrange
  endif else begin
    id = sort( xvec )
    limit[[0,2]] = [ xvec[id[0]] - ( xvec[id[1]] - xvec[id[0]] ) / 2., $
        xvec[id[n_xvec-1]] + ( xvec[id[n_xvec-1]] - xvec[id[n_xvec-2]] ) / 2. ]
    id = -1
  endelse
  if keyword_set( yrange ) then begin
    limit[[1,3]] = yrange
  endif else begin
    id = sort( yvec )
    limit[[1,3]] = [ yvec[id[0]] - ( yvec[id[1]] - yvec[id[0]] ) / 2., $
        yvec[id[n_yvec-1]] + ( yvec[id[n_yvec-1]] - yvec[id[n_yvec-2]] ) / 2. ]
    id = -1
  endelse
endif
n_limit = n_elements( limit )

;***********************************************************************
; Plot Map

; Determine the longitude and latitude extent
if n_limit eq 4 then begin
  limit_xvec_min = limit[0]
  limit_xvec_max = limit[2]
  limit_yvec_min = limit[1]
  limit_yvec_max = limit[3]
endif else begin
  limit_xvec_min = min( limit[[0,2,4,6]] )
  limit_xvec_max = max( limit[[0,2,4,6]] )
  limit_yvec_min = min( limit[[1,3,5,7]] )
  limit_yvec_max = max( limit[[1,3,5,7]] )
endelse

; Contour over map.
; If we have irregularly gridded data
if irregular_opt eq 1 then begin
  ; Iterate through elements in the data
  for i_data = 0l, n_data - 1l do begin
    ; Calculate dimensional distances from this point
    near_xvec_delta = xvec_use - xvec_use[i_data]
    near_yvec_delta = yvec_use - yvec_use[i_data]
    ; Find n_vertex neighbouring points to this one
    near_delta = near_xvec_delta ^ 2 + near_yvec_delta ^ 2
    id_data = sort( near_delta )
    n_vertex = n_vertex_0
    id_data = id_data[1l:n_vertex]
    near_xvec_delta = near_xvec_delta[id_data]
    near_yvec_delta = near_yvec_delta[id_data]
    near_delta = near_delta[id_data]
    if outline_opt eq 1 then near_data = data_use[id_data]
    ; Sort vertices in counterclockwise order from the smallest angle
    near_ang = atan( near_yvec_delta, near_xvec_delta )
    id_near = sort( near_ang )
    near_xvec_delta = near_xvec_delta[id_near]
    near_yvec_delta = near_yvec_delta[id_near]
    near_ang = near_ang[id_near]
    near_delta = near_delta[id_near]
    if outline_opt eq 1 then near_data = near_data[id_near]
    ; If the largest angle is near or greater than pi
    temp_near_ang_diff = [ near_ang[1:n_vertex-1], near_ang[0] + 2. * !pi ] $
        - near_ang
    if max( temp_near_ang_diff ) gt 0.9 * !pi then begin
      ; Retain the closest points
      id_near = sort( near_delta )
      ; If more than 1.5*pi is missing then retain only half the points
      if max( temp_near_ang_diff ) gt 1.4 * !pi then begin
        n_vertex = n_vertex / 2
      ; Otherwise retain 3/4 of the points
      endif else begin
        n_vertex = 3 * n_vertex / 4
      endelse
      id_near = id_near[0:n_vertex-1]
      near_xvec_delta = near_xvec_delta[id_near]
      near_yvec_delta = near_yvec_delta[id_near]
      near_ang = near_ang[id_near]
      near_delta = near_delta[id_near]
      if outline_opt eq 1 then near_data = near_data[id_near]
      ; Create fake points reflected across our point
      near_xvec_delta = [ near_xvec_delta, -near_xvec_delta ]
      near_yvec_delta = [ near_yvec_delta, -near_yvec_delta ]
      near_delta = [ near_delta, near_delta ]
      near_ang = atan( near_yvec_delta, near_xvec_delta )
      if outline_opt eq 1 then begin
        ; Set data values according to nearest available point, rather than the 
        ; reflected point
        near_data = [ near_data, near_data ]
        for i_near = 0, n_vertex - 1 do begin
          temp = ( xvec_use - near_xvec_delta[n_vertex+i_near] $
              - xvec_use[i_data] ) ^ 2 $
              + ( yvec_use - near_yvec_delta[n_vertex+i_near] $
              - yvec_use[i_data] ) ^ 2
          id_sort = sort( temp )
          if id_sort[0] eq i_data then begin
            near_data[n_vertex+i_near] = data_use[id_sort[1]]
          endif else begin
            near_data[n_vertex+i_near] = data_use[id_sort[0]]
          endelse
        endfor
      endif
      n_vertex = n_vertex * 2
      ; Restrict to n_vertex_0 points
      if n_vertex gt n_vertex_0 then begin
        id_near = sort( near_ang )
        near_xvec_delta = near_xvec_delta[id_near]
        near_yvec_delta = near_yvec_delta[id_near]
        near_ang = near_ang[id_near]
        near_delta = near_delta[id_near]
        if outline_opt eq 1 then near_data = near_data[id_near]
        temp_near_ang_diff $
            = [ near_ang[1:n_vertex-1], near_ang[0] + 2. * !pi ] - near_ang
        id_near = sort( temp_near_ang_diff )
        id_near = id_near[n_vertex-n_vertex_0:n_vertex-1]
        near_xvec_delta = near_xvec_delta[id_near]
        near_yvec_delta = near_yvec_delta[id_near]
        near_ang = 0
        near_delta  = 0
        if outline_opt eq 1 then near_data = near_data[id_near]
        n_vertex = n_vertex_0
      endif
    endif
    ; Determine vertices of polygon
    vertex_xvec = fltarr( n_vertex )
    vertex_yvec = fltarr( n_vertex )
    vertex_xvec[0] = ( near_xvec_delta[n_vertex-1] + near_xvec_delta[0] ) / 2.
    vertex_yvec[0] = ( near_yvec_delta[n_vertex-1] + near_yvec_delta[0] ) / 2.
    for i_vertex = 1, n_vertex - 1 do begin
      vertex_xvec[i_vertex] = ( near_xvec_delta[i_vertex-1] $
          + near_xvec_delta[i_vertex] ) / 2.
      vertex_yvec[i_vertex] = ( near_yvec_delta[i_vertex-1] $
          + near_yvec_delta[i_vertex] ) / 2.
    endfor
    vertex_xvec = vertex_xvec + xvec_use[i_data]
    vertex_yvec = vertex_yvec + yvec_use[i_data]
    ; If we are drawing filled contours
    if ( outline_opt eq 0 ) and ( no_plot_opt eq 0 ) then begin
      ; Determine colour
      id = max( where( levels - data_use[i_data] lt 0 ) )
      if id eq -1 then id = 0
      ; Plot cell
      polyfill, vertex_xvec, vertex_yvec, color=c_colors[id]
    ; If we are drawing outline contours
    endif else begin
      ; Determine the contour level of this cell
      id_level = max( where( levels - data_use[i_data] lt 0 ) )
      ; Iterate through neighbouring cells
      for i_near = 0, n_vertex - 1 do begin
        ; Determine the contour level of the cell to the east
        id_level_1 = max( where( levels - near_data[i_near] lt 0 ) )
        ; If the contour levels are different then plot a contour line
        if id_level_1 ne id_level then begin
          if i_near eq n_vertex - 1 then begin
            temp_xvec = vertex_xvec[[n_vertex-1,0]]
            temp_yvec = vertex_yvec[[n_vertex-1,0]]
          endif else begin
            temp_xvec = vertex_xvec[[i_near,i_near+1]]
            temp_yvec = vertex_yvec[[i_near,i_near+1]]
          endelse
          if no_plot_opt eq 0 then begin
            plots, temp_xvec, temp_yvec, color=color, thick=thick
          endif
          ; Plot a down-/up-hill tick if requested
          if keyword_set( point_ticks ) and ( no_plot_opt eq 0 ) then begin
            d_tick_xvec = [ 0, $
                ( near_xvec_delta[i_near] $
                - ( mean( vertex_xvec ) - xvec_use[i_data] ) ) / 4. ]
            d_tick_yvec = [ 0, $
                ( near_xvec_delta[i_near] $
                - ( mean( vertex_yvec ) - yvec_use[i_data] ) ) / 4. ]
            tick_xvec = [ 0, 0 ] + mean( temp_xvec )
            tick_yvec = [ 0, 0 ] + mean( temp_yvec )
            temp = point_ticks * sign( id_level_1 - id_level )
            oplot, tick_xvec+temp*d_tick_xvec, tick_yvec+temp*d_tick_yvec, $
                color=color, thick=thick
          endif
        endif
      endfor
    endelse
  endfor
; If we have regularly gridded data
endif else begin
  ; Initialise array identifying shape membership
  index_shape = lonarr( n_xvec, n_yvec )
  index_shape_max = 1
  ; Initialise the vector remembering removed (unused) shape counter
  id_shape_unused = 0l
  n_id_shape_unused = n_elements( id_shape_unused ) - 1l
  ; Iterate through rows
  for i_yvec = 0l, n_yvec - 1l do begin
    ; If this row is within the mapping range
    if ( yvec_use[i_yvec] ge limit_yvec_min ) $
        and ( yvec_use[i_yvec] le limit_yvec_max ) then begin
      ; Copy data for this row and the one below
      temp_data_use_row = data_use[*,i_yvec]
      temp_index_shape_row = index_shape[*,i_yvec]
      if i_yvec gt 0l then begin
        if yvec_use[i_yvec-1l] ge limit_yvec_min then begin
          temp_data_use_row_below = data_use[*,i_yvec-1l]
          temp_index_shape_row_below = index_shape[*,i_yvec-1l]
        endif
      endif
      ; Iterate through columns
      for i_xvec = 0l, n_xvec - 1l do begin
        ; If this column is within the mapping range
        if ( xvec_use[i_xvec] ge limit_xvec_min ) $
            and ( xvec_use[i_xvec] le limit_xvec_max ) then begin
          ; Copy data for this cell
          temp_data_use = temp_data_use_row[i_xvec]
          ; Determine if this cell has data
          cell_finite_opt = finite( temp_data_use )
          ; Note the level of this cell
          if cell_finite_opt eq 1 then begin
            id_level_here = max( where( levels - temp_data_use lt 0 ) )
            if id_level_here eq -1 then id_level_here = 0
          endif else begin
            id_level_here = -2
          endelse
          ; Initialise the value for this cell (which may already have been 
          ; calculated to something other than 0)
          temp_index = 0l
          id_level_left = -1
          id_level_below = -1
          ; Check if this matches the level of the cell to the left
          if i_xvec gt 0 then begin
            if xvec_use[i_xvec-1l] ge limit_xvec_min then begin
              temp_data_use_left = temp_data_use_row[i_xvec-1l]
              if cell_finite_opt eq 1 then begin
                if finite( temp_data_use_left ) eq 1 then begin
                  id_level_left = max( $
                      where( levels - temp_data_use_left lt 0 ) )
                  if id_level_left eq -1 then id_level_left = 0
                endif
              endif else begin
                if finite( temp_data_use_left ) eq 0 then id_level_left = -2
              endelse
              if id_level_left eq id_level_here then begin
                temp_index = temp_index_shape_row[i_xvec-1l]
              endif
            endif
          endif
          ; Check if this matches the level of the cell below
          ; (Check even if we already assigned to a shape above, in case we 
          ; find that we have now arrived at where the same shape was 
          ; assigned a different identifier in the previous row.)
          if i_yvec gt 0 then begin
            if yvec_use[i_yvec-1l] ge limit_yvec_min then begin
              temp_data_use_below = temp_data_use_row_below[i_xvec]
              if cell_finite_opt eq 1 then begin
                if finite( temp_data_use_below ) eq 1 then begin
                  id_level_below = max( $ 
                      where( levels - temp_data_use_below lt 0 ) )
                  if id_level_below eq -1 then id_level_below = 0
                endif
              endif else begin
                if finite( temp_data_use_below ) eq 0 then id_level_below = -2
              endelse
              if id_level_below eq id_level_here then begin
                ; Copy the index value for the cell below
                temp_index_below = temp_index_shape_row_below[i_xvec]
                ; If we should adopt the value for the cell below
                if temp_index eq 0 then begin
                  temp_index = temp_index_below
                ; If we need to revise the shape identifier on this cell
                endif else if temp_index ne temp_index_below then begin
                  ; If we need to revise the shape identifier on cells to the 
                  ; left as well
                  if id_level_left eq id_level_below then begin
                    if abs( temp_index_below ) gt abs( temp_index ) then begin
                      temp_index_use = temp_index
                      temp_index_discard = temp_index_below
                    endif else begin
                      temp_index_use = temp_index_below
                      temp_index_discard = temp_index
                    endelse
                    id = where( index_shape eq temp_index_discard, n_id )
                    index_shape[id] = temp_index_use
                    id = where( temp_index_shape_row eq temp_index_discard, $
                        n_id )
                    if n_id gt 0 then begin
                      temp_index_shape_row[id] = temp_index_use
                    endif
                    id = where( temp_index_shape_row_below $
                        eq temp_index_discard, n_id )
                    if n_id gt 0 then begin
                      temp_index_shape_row_below[id] = temp_index_use
                    endif
                    ; Change the shape identifier
                    temp_index = temp_index_use
                    ; Record the identifier that was going to be used
                    id_shape_unused = [ id_shape_unused, $
                        abs( temp_index_discard ) ]
                    if abs( temp_index ) eq index_shape_max then begin
                      index_shape_max = max( abs( index_shape ) )
                    endif
                  endif else begin
                    ; Change the shape identifier
                    temp_index = temp_index_below
                  endelse
                endif
              endif
            endif
          endif
          ; If we have not matched to a neighbouring cell then start a new
          ; shape
          if temp_index eq 0 then begin
            n_id_shape_unused = n_elements( id_shape_unused ) - 1l
            if n_id_shape_unused gt 0 then begin
              temp_index = id_shape_unused[n_id_shape_unused]
              id_shape_unused = id_shape_unused[0l:n_id_shape_unused-1l]
            endif else begin
              temp_index = max( [ 2, index_shape_max + 1 ] )
            endelse
            if temp_index gt index_shape_max then index_shape_max = temp_index
            ; Flag NaN shapes
            if cell_finite_opt eq 0 then temp_index = -abs( temp_index )
          endif
          ; Record the value
          index_shape[i_xvec,i_yvec] = temp_index
          temp_index_shape_row[i_xvec] = temp_index
        endif
      endfor
    endif
  endfor
  ; Count the number of shapes
  n_shape = max( abs( index_shape ) ) - 1
  ; Sort shapes by size (largest to smallest).  If a shape is nested in a 
  ; larger shape then it will be overwritten by the larger shape if the 
  ; larger is plotted afterward.
  area = fltarr( n_shape )
  for i_shape = 0l, n_shape - 1l do begin
    id = where( abs( index_shape ) eq i_shape + 2l, n_id )
    area[i_shape] = n_id
  endfor
  id_sort_shape = reverse( sort( area ) )
  ; Iterate through shapes
  for i_shape = 0l, n_shape - 1l do begin
    ; Find all cells in this shape
    temp_id_sort_shape = id_sort_shape[i_shape] + 2l
    id_shape = where( index_shape eq temp_id_sort_shape, n_id_shape )
    if n_id_shape eq 0 then begin
      temp_id_sort_shape = -1 * temp_id_sort_shape
      id_shape = where( index_shape eq temp_id_sort_shape, n_id_shape )
    endif
    if n_id_shape gt 0 then begin
      ; Determine the x and y coordinates of the first cell
      ; (Note we are assuming it is on the border, which it has to be.)
      id_shape_x = id_shape[0] mod n_xvec
      id_shape_y = id_shape[0] / n_xvec
      ; Initialise counter for border points
      i_point = 0l
      ; Initialise a vector for instructing the order of sides to examine in 
      ; the next cell.  Start with left, top, right, bottom ([0,1,2,3]).
      index_side = [ 0, 1, 2, 3 ]
      ; Iterate through cells until we return to the original spot
      flag_done_shape = 0
      while flag_done_shape eq 0 do begin
        ; Determine the shape values of this cell and the one to the left, top, 
        ; right, and bottom, as well as the distance to the edges with those 
        ; neighbouring cells from the centre of this cell
        temp_shape_here = index_shape[id_shape_x,id_shape_y]
        temp_shape_near = -1 + intarr( 4 )
        temp_shape_dist = fltarr( 4 )
        if id_shape_x ne 0 then begin
          temp_shape_near[0] = index_shape[id_shape_x-1,id_shape_y]
          ;id = id_shape_x + [ -1, 0 ]
          temp_shape_dist[0] = xvec_use[id_shape_x] $
              - ( xvec_use[id_shape_x] - xvec_use[id_shape_x-1] ) / 2.
        endif else begin
          temp_shape_dist[0] = xvec_use[0] - ( xvec_use[1] - xvec_use[0] ) / 2.
        endelse
        if id_shape_y ne n_yvec - 1 then begin
          temp_shape_near[1] = index_shape[id_shape_x,id_shape_y+1]
          temp_shape_dist[1] = yvec_use[id_shape_y] $
              + ( yvec_use[id_shape_y+1] - yvec_use[id_shape_y] ) / 2.
        endif else begin
          temp_shape_dist[1] = yvec_use[n_yvec-1] $
              + ( yvec_use[n_yvec-1] - yvec_use[n_yvec-2] ) / 2.
        endelse
        if id_shape_x ne n_xvec - 1 then begin
          temp_shape_near[2] = index_shape[id_shape_x+1,id_shape_y]
          temp_shape_dist[2] = xvec_use[id_shape_x] $
              + ( xvec_use[id_shape_x+1] - xvec_use[id_shape_x] ) / 2.
        endif else begin
          temp_shape_dist[2] = xvec_use[n_xvec-1] $
              + ( xvec_use[n_xvec-1] - xvec_use[n_xvec-2] ) / 2.
        endelse
        if id_shape_y ne 0 then begin
          temp_shape_near[3] = index_shape[id_shape_x,id_shape_y-1]
          temp_shape_dist[3] = yvec_use[id_shape_y] $
              - ( yvec_use[id_shape_y] - yvec_use[id_shape_y-1] ) / 2.
        endif else begin
          temp_shape_dist[3] = yvec_use[0] - ( yvec_use[1] - yvec_use[0] ) / 2.
        endelse
        ; Iterate through sides until we come to a non-border
        flag_done_sides = 0
        for i_side = 0, 3 do begin
          if flag_done_sides eq 0 then begin
            ; Determine the identity of the side to examine
            ; (0=left, 1=top, 2=right, 3=bottom)
            id_side = index_side[i_side]
            ; If this is on the border
            if temp_shape_here ne temp_shape_near[id_side] then begin
              ; If this is the first edge on the border then record the first 
              ; (counterclockwise-most) point on the edge
              if i_point eq 0l then begin
                id_before = ( id_side - 1 + 4 ) mod 4
                if ( id_side eq 0 ) or ( id_side eq 2 ) then begin
                  border = [ temp_shape_dist[id_side], $
                      temp_shape_dist[id_before] ]
                endif else begin
                  border = [ temp_shape_dist[id_before], $
                      temp_shape_dist[id_side] ]
                endelse
                i_point = 1l
              endif
              ; Record the second (clockwise-most) point on the edge
              border = [ [ border ], [ 0, 0 ] ]
              id_after = ( id_side + 1 ) mod 4
              if ( id_side eq 0 ) or ( id_side eq 2 ) then begin
                border[*,i_point] = [ temp_shape_dist[id_side], $
                   temp_shape_dist[id_after] ]
              endif else begin
                border[*,i_point] = [ temp_shape_dist[id_after], $
                   temp_shape_dist[id_side] ]
              endelse
              i_point = i_point + 1l
              ; Check if we have returned to the start point
              temp = abs( border[0,i_point-1l] - border[0,0] )
              if temp lt d_xvec_use_min / 2. then begin
                temp = abs( border[1,i_point-1l] - border[1,0] )
                if temp lt d_yvec_use_min / 2. then begin
                  ; If we have returned, then exit the loop for this shape
                  flag_done_sides = 1
                  flag_done_shape = 1
                endif
              endif
            ; If this is not a border then proceed to that cell
            endif else begin
              ; If this is the left border, then move to the left
              if id_side eq 0 then begin
                id_shape_x = id_shape_x - 1
                index_side = [ 3, 0, 1, 2 ]
              ; If this is the top border, then move up
              endif else if id_side eq 1 then begin
                id_shape_y = id_shape_y + 1
                index_side = [ 0, 1, 2, 3 ]
              ; If this is the right border, then move to the right
              endif else if id_side eq 2 then begin
                id_shape_x = id_shape_x + 1
                index_side = [ 1, 2, 3, 0 ]
              ; If this is the bottom border, then move down
              endif else if id_side eq 3 then begin
                id_shape_y = id_shape_y - 1
                index_side = [ 2, 3, 0, 1 ]
              endif
              ; Flag the exit from this cell
              flag_done_sides = 1
            endelse
          endif
        endfor
      endwhile
      ; If plotting is requested
      if no_plot_opt eq 0 then begin
        ; If we are only outlining the shape
        if ( outline_opt eq 1 ) and ( temp_id_sort_shape ge 0 ) then begin
          ; Plot the outline of the polygon
          plots, reform( border[0,*] ), reform( border[1,*] ), color=color, $
              thick=thick
        ; If we are filling
        endif else begin
          ; Determine the colour
          if temp_id_sort_shape ge 0 then begin
            id_shape_x = id_shape[0] mod n_xvec
            id_shape_y = id_shape[0] / n_xvec
            id_level = max( $
                where( levels - data_use[id_shape_x,id_shape_y] lt 0 ) )
            if id_level eq -1 then id_level = 0
            temp_color = c_colors[id_level]
          endif else begin
            temp_color = !p.background
          endelse
          ; Plot the filled polygon
          polyfill, reform( border[0,*] ), reform( border[1,*] ), $
              color=temp_color
        endelse
      endif
      ;print,i_shape,n_elements( border[0,*] )
      ; Record border coordinates for output, if requested
      if return_border_opt eq 1 then begin
        ; Determine the feature this belongs to (each feature potentially 
        ; containing multiple shapes)
        id_shape_x = id_shape[0] mod n_xvec
        id_shape_y = id_shape[0] / n_xvec
        id_level = max( $
            where( levels - data_use[id_shape_x,id_shape_y] lt 0 ) )
        if id_level eq -1 then id_level = 0
        ; Determine the number of vertices to add
        n_temp = n_elements( border[0,*] )
        ; Initialise the border information array for this shape, including 
        ; coordinates
        temp_border = [ border, fltarr( 2, n_temp ) ]
        ; Add feature identifier
        temp_border[2,*] = id_level
        ; Determine the shape identifier within the feature
        if not( keyword_set( border_out ) ) then begin
          id_subshape = 0
        endif else begin
          id = where( border_out[2,*] eq id_level, n_id )
          if n_id eq 0 then begin
            id_subshape = 0
          endif else begin
            id_subshape = max( border_out[3,id] ) + 1
          endelse
        endelse
        temp_border[3,*] = id_subshape
        ; Copy this shape to the output array
        if not( keyword_set( border_out ) ) then begin
          border_out = temp_border
        endif else begin
          border_out = [ [ border_out ], [ temp_border ] ]
        endelse
      endif
    endif
  endfor
endelse

;***********************************************************************
; The End

;stop
return
END
