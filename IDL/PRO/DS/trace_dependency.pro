;+
; NAME:
;    trace_dependency
;
; PURPOSE:
;    This procedure plots a diagram showing the dependency tree of the
;    requested IDL programs and returns the list of dependencies.
;
; CATEGORY:
;    Miscellaneous
;
; CALLING SEQUENCE:
;    trace_dependency, program
;
; INPUTS:
;    PROGRAM:  A required string vector listing the IDL routines to examine.
;
; KEYWORD PARAMETERS:
;    CODE_NAME:  A string output vector of length N_CODE listing the routines 
;        used by PROGRAM.  PROGRAM itself is included in the list.
;    CODE_DEPENDENCY:  A string output vector of length N_CODE listing the 
;        direct dependencies of each routine specified in N_CODE.  Each list is 
;        in comma-delimited format.  The routine listed in CODE_NAME is 
;        included in the dependency list.
;    CODE_COLOR:  An integer output of length N_CODE listing colour indices 
;        associated with the routines listed in CODE_NAME.  If a specified 
;        routine was found then the colour index is that specified in the 
;        COLOR_GOOD constant, otherwise the COLOR_MISSING constant value is 
;        used (see code for values).
;    IDL_PATH:  An optional input string vector listing the directories in $
;        which to search for the routines listed in PROGRAM.
;    NO_PLOT:  If set then no plot is produced, with the keyword output still 
;        output.  The default is to produce a plot.
;
; OUTPUTS:
;    CODE_COLOR, CODE_DEPENDENCY, CODE_NAME
;
; USES:
;    string_from_vector.pro
;
; PROCEDURE:
;    This code find the list of dependencies in the "USES:" header field of 
;    each file listed in PROGRAM, and then iteratively performs the same search 
;    on the dependencies until the end of the dependency tree is reached.
;
; EXAMPLE:
;    ; Run this is the directory containing the published IDL library code:
;      trace_dependency, 'trace_dependency.pro', no_plot=1, code_name=code_name
;    ; CODE NAME should return as [ 'trace_dependency.pro', 'string_from_vector.pro', 'str.pro', 'dimension.pro', 'var_type.pro' ]
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone (dastone@runbox.com), 2012-05-25
;    Modified:  DAS, 2018-08-10 (Added standard documentation;  Added 
;        IDL_PATH keyword input)
;-

;***********************************************************************

PRO TRACE_DEPENDENCY, $
    PROGRAM, $
    IDL_PATH=idl_path, $
    NO_PLOT=no_plot_opt, $
    CODE_NAME=code_name, CODE_DEPENDENCY=code_dependency, CODE_COLOR=code_color

;***********************************************************************
; Constants

; Confirm program is defined
if not( keyword_set( program ) ) then stop
if min( strlen( program ) ) eq 0 then stop
;print,program
;print

; The IDL path
if n_elements( idl_path ) eq 0 then idl_path = ''

; The number of programmes to examine
n_program = n_elements( program )

; The default colours
color_good = 4
color_missing = 2
color_link = 15
color_highlight = 3

; Option for no plotting
no_plot_opt = keyword_set( no_plot_opt )

; Default plot settings
charsize_0 = 1.
thick = 2

; Confirm output vectors are clear
code_name = ''
code_dependency = ''
code_color = 0

;***********************************************************************
; Iteratively find all files in these dependency trees

; Iterate through programs
for i_program = 0, n_program - 1 do begin
  ; Find the program
  temp_file = file_search( idl_path, program[i_program], count=temp_file_count )
  ;; If multiple copies exist
  ;if temp_file_count gt 1 then begin
  ;  ; We are lost
  ;  stop
  ; If the file could not be found
  ;endif else if temp_file_count eq 0 then begin
  if temp_file_count eq 0 then begin
    ; Set the colour to missing
    temp_color = color_missing
    ; Set blank dependency list
    temp_code = program[i_program]
  ; If the file exists
  endif else begin
    ; Set the colour to good
    temp_color = color_good
    ; Set blank dependency list
    temp_code = program[i_program]
    ; If this is an IDL routine
    if strpos( program[i_program], '.pro' ) gt 0 then begin
      ; Prepare to find the dependencies
      flag_uses = 0
      temp = ''
      ; Open the file for reading
      openr, 1, temp_file[0]
      ; Iterate through code lines
      while flag_uses ne 2 do begin
        ; Read the line
        readf, 1, temp
        ; If this is the beginning of the USES list
        if ( strpos( temp, '; USES:' ) eq 0 ) $
            or ( strpos( temp, ';USES:' ) eq 0 ) then begin
          ; Flag it
          flag_uses = 1
          ; If we are within the USES list
        endif else if flag_uses eq 1 then begin
          ; If this is the end of the USES list
          if ( strpos( temp, '; EXAMPLE' ) eq 0 ) or ( temp eq '' ) $
              or ( strpos( temp, '; PROCEDURE' ) eq 0 ) $
              or ( strpos( temp, '---' ) gt 0 ) $
              or ( strpos( temp, ';    -' ) eq 0 ) or ( temp eq ';' ) $
              or ( strpos( temp, ';     -' ) eq 0 ) $
              or ( strpos( temp, ';   -' ) eq 0 ) $
              or ( strpos( temp, ';  -' ) eq 0 ) $
              or ( strpos( temp, ';-' ) eq 0 ) then begin
            ; Flag it
            flag_uses = 2
          ; Otherwise this is a dependency
          endif else begin
            ; Remove any bracketed comments
            temp = ( strsplit( temp, '(', extract=1 ) )[0]
            ; Record the dependency
            temp_len = strlen( temp )
            temp = strtrim( strmid( strlowcase( temp ), 1, temp_len - 1 ), 2 )
            temp_code = [ temp_code, temp ]
          endelse
        ; Otherwise check if we are at the end of a documentation header 
        ; without a USES list
        endif else begin
          if strpos( temp, ';-' ) eq 0 then flag_uses = 2
        endelse
      endwhile
      ; Close the file
      close, 1
    endif
  endelse
  ; Combine dependencies into a list
  n_temp_code = n_elements( temp_code )
  if n_temp_code eq 1 then begin
    temp_code = [ temp_code, '' ]
  endif else begin
    temp_code = [ temp_code[0], $
        string_from_vector( temp_code[1:n_temp_code-1], nospace=1 ) ]
  endelse
  ; Record the details of this file
  if not( keyword_set( code_name ) ) then begin
    code_name = temp_code[0]
    code_dependency = temp_code[1]
    code_color = temp_color
  endif else begin
    code_name = [ code_name, temp_code[0] ]
    code_dependency = [ code_dependency, temp_code[1] ]
    code_color = temp_color
  endelse
  ; Find dependencies of the dependencies
  if n_temp_code gt 1 then begin
    ; Call trace_dependency.pro again but with these dependencies
    temp_program = strsplit( temp_code[1], ',', extract=1 )
    trace_dependency, temp_program, no_plot=1, code_name=temp_code_name, $
        code_dependency=temp_code_dependency, code_color=temp_code_color, $
        idl_path=idl_path
    ; Add these to the file list
    code_name = [ code_name, temp_code_name ]
    code_dependency = [ code_dependency, temp_code_dependency ]
    code_color = [ code_color, temp_code_color ] 
  endif
endfor

; Remove redundancies
n_code = n_elements( code_name )
for i_code = 0, n_code - 2 do begin
  if code_name[i_code] ne '' then begin
    id = where( code_name[i_code+1:n_code-1] eq code_name[i_code], n_id )
    if n_id gt 0 then code_name[i_code+1+id] = ''
  endif
endfor
id = where( code_name ne '', n_code )
code_name = code_name[id]
code_dependency = code_dependency[id]
code_color = code_color[id]

;***********************************************************************
; Create dependency trees

; Only do this if we are plotting
if no_plot_opt eq 0 then begin

  ; Determine rank of each program
  code_rank = n_code + 1 + intarr( n_code )
  ; Find lowest ranked programs
  id_low = where( code_dependency eq '', n_id_low )
  ; Rank these programs
  code_rank[id_low] = n_code
  ; Iterate through ranks from lowest upward
  for i_rank = 0, n_code - 1 do begin
    ; Iterate through programs of this rank
    id_here = where( code_rank eq n_code - i_rank, n_id_here )
    for i_here = 0, n_id_here - 1 do begin
      ; Determine if program is a dependency of something else
      temp = strpos( code_dependency, code_name[id_here[i_here]] )
      ;id = where( ( temp ge 0 ) and ( code_rank eq n_code + 1 ), n_id )
      id = where( temp ge 0, n_id )
      ; Assign rank to this higher level
      if n_id gt 0 then code_rank[id] = n_code - i_rank - 1
    endfor
  endfor
  ; Re-write rank starting from 0 downward
  rank_new = -1 + intarr( n_code )
  for i_rank = 0, n_code - 1 do begin
    id = where( code_rank eq i_rank + 1, n_id )
    if n_id gt 0 then rank_new[id] = max( rank_new ) + 1
  endfor
  code_rank = temporary( rank_new )
  ; Count the ranks
  n_rank = max( code_rank ) + 1

endif

;***********************************************************************
; Plot tree

; Only do this if we are plotting
if no_plot_opt eq 0 then begin

  ; Initialise plotting window
  plot, [0,1], [0,1], nodata=1, xstyle=4, ystyle=4, xmargin=[0,0], ymargin=[0,0]
  tek_color

  ; Determine the character size
  charsize = charsize_0 / n_rank * 5

  ; Determine position of program names.
  ; Initialise x-y position array
  code_pos = -1 + fltarr( n_code, 2 )
  ; Iterate through ranks
  for i_rank = 0, n_rank - 1 do begin
    ; Identify programs in this rank
    id_rank = where( code_rank eq i_rank, n_id_rank )
    ; If this is the top rank
    if i_rank eq 0 then begin
      ; Then randomly sort
      id_sort = indgen( n_id_rank )
    ; Otherwise select order best matched to higher rank
    endif else begin
      ; Initialise pseudo x-position vector
      temp_x = fltarr( n_id_rank )
      ; Find members of the higher ranks
      id_rank_high = where( code_rank le i_rank - 1, n_id_rank_high )
      ; Iterate through members of the lower rank
      for i_id_rank = 0, n_id_rank - 1 do begin
        ; Find the middle x-position for linking to higher ranks
        temp = strpos( code_dependency[id_rank_high], $
            code_name[id_rank[i_id_rank]] )
        id = where( temp ge 0 )
        temp_x[i_id_rank] = mean( code_pos[id_rank_high[id],0] )
      endfor
      ; Sort x-positions
      id_sort = -1 + intarr( n_id_rank )
      for i_sort = 0, n_id_rank - 1 do begin
        id = where( ( temp_x eq min( temp_x ) ) and ( temp_x lt 1.9 ), n_id )
        if n_id gt 0 then begin
          id_sort[id] = max( id_sort ) + 1 + indgen( n_id )
          temp_x[id] = 2
        endif
      endfor
    endelse
    ; Assign horizontal positions
    code_pos[id_rank,0] = float( id_sort ) / n_id_rank + 0.5 / n_id_rank $
        + randomn( seed, n_id_rank ) / n_id_rank / 10.
    ; Assign vertical positions
    code_pos[id_rank,1] = 1. - float( i_rank ) / n_rank - 0.5 / n_rank
  endfor

  ; Highlight top ranker programs.
  ; Iterate through the programs
  for i_code = 0, n_code - 1 do begin
    ; Check if it is a top ranker
    temp = max( strpos( code_dependency, code_name[i_code] ) )
    if temp lt 0 then begin
      ; Highlight its spot
      temp_x = code_pos[i_code,0] + [ -0.15, 0.15 ] * charsize
      temp_y = code_pos[i_code,1] + [ -0.02, 0.04 ] * charsize
      polyfill, temp_x[[0,1,1,0,0]], temp_y[[0,0,1,1,0]], color=color_highlight
    endif
  endfor

  ; Plot program links.
  ; Iterate through programs
  for i_code = 0, n_code - 1 do begin
    ; Iterate through dependencies
    dependency = strsplit( code_dependency[i_code], ',', extract=1, $
        count=n_dependency )
    for i_dependency = 0, n_dependency - 1 do begin
      ; Find the dependency
      id_dependency = where( code_name eq dependency[i_dependency] )
      ; Plot the link
      temp_x = [ code_pos[i_code,0], code_pos[id_dependency[0],0] ]
      temp_y = [ code_pos[i_code,1] - 0.01 * charsize, $
          code_pos[id_dependency[0],1] + 0.03 * charsize ]
      plots, temp_x, temp_y, thick=thick, color=color_link
    endfor
  endfor

  ; Print program names.
  ; Iterate through the programs
  for i_code = 0, n_code - 1 do begin
    ; Print the program name
    xyouts, code_pos[i_code,0], code_pos[i_code,1], code_name[i_code], $
        color=code_color[i_code], alignment=0.5, charsize=charsize, font=0
  endfor

endif

;***********************************************************************
; The end

return
END
