;+
; NAME:
;    SREAD
;
; PURPOSE:
;    This function extracts column delimited data (whitespace between columns) 
;    of data from a file and loads it into an string array.  
;
; CATEGORY:
;    Input/Output
;
; CALLING SEQUENCE:
;    Result = sread( FILENAME )
;
; INPUTS:
;    FILENAME:  The name of the file to open and read.
;    COMMENT_CHAR, NOCOMPRES, SEPARATOR, SKIP_LINES, TRUNCATE, VERBOSE
;
; KEYWORD PARAMETERS:
;    NO_COLUMNS:  If set then no column separation is performed.  The default 
;        is to separate the input into columns.
;    COMMENT_CHAR:  A string labeling the character used to denote comment 
;        within the file being read, and thus to be skipped.  If not input then 
;        ';' is assumed.  If set to '' then no comment lines are assumed.
;    NOCOMPRESS:  If set then whitespaces are not compressed to length 1.  The 
;        default is compression.
;    PRESERVE_NULL:  The PRESERVE_NULL keyword for the IDL strsplit function.
;    PROTECTOR:  A string which, when present in pairs, protects splitting 
;        within the pair.  For instance, if SEPARATOR=',' and the line read 
;        from the file is '"a,b","c,d"', then if PROTECTOR='"' the resulting 
;        split will be '"a,b"' and '"c,d"'.  The default is no protector.
;    SEPARATOR:  The string that separates columns of data.
;    SKIP_LINES: An integer number of lines to skip from the beginning of the 
;        file.
;    TRUNCATE:  If this keyword is set any (non comment) lines longer than the 
;        first line will be truncated.
;    VERBOSE:  Write some information about the reading process to the terminal.
;
; OUTPUTS:
;    Result:  A string array matching the content of the file is returned.  A 
;        null string ('') is returned if there is an error.
;
; USES:
;    odd.pro
;    str.pro
;
; RESTRICTIONS:
;    I am sure there are some but they are unknown.
;
; EXAMPLE:
;    text = sread( file, SKIP_LINES=5 )
;
; MODIFICATION HISTORY:
;    Written by:  Edward C. Wiebe, 2000-06-09.
;    Modified:  Edward C. Wiebe, 2002-02-07 (changed code to always return '' 
;        if there is an error)
;    Modified:  Edward C. Wiebe, 2002-03-05 (fixed a bug where the routine 
;        stopped with an error if the input filename was '' -- it now returns 
;        '') 
;    Modified:  Daithi A. Stone (stoned@atm.ox.ac.uk), 2005-03-10 (switched use 
;        of obsolete findfile to file_search)
;    Modified:  DAS, 2007-01-17 (added NOCOMPRESS keyword)
;    Modified:  DAS, 2009-09-29 (changed use of stringc to str.pro)
;    Modified:  DAS, 2014-03-31 (Added COMMENT_CHAR keyword)
;    Modified:  DAS, 2014-08-15 (Fixed bug in treatment of blank lines)
;    Modified:  DAS, 2015-02-27 (Added the PRESERVE_NULL keyword)
;    Modified:  DAS, 2017-10-09 (Added the NO_COLUMNS keyword;  fixed bug where 
;        blank lines could interfere with counting of columns first lines 
;        include COMMENT_CHAR)
;    Modified:  DAS, 2018-05-28 (Added PROTECTOR keyword input)
;-

;***********************************************************************

FUNCTION sread, $
    FILENAME, $
    COMMENT_CHAR=comment_char, $
    SKIP_LINES=skip_lines, $  
    VERBOSE=verbose, $
    PROTECTOR=protector, $
    SEPARATOR=separator, $
    TRUNCATE=truncate, $
    DEBUG=debug, $
    PRESERVE_NULL=preserve_null_opt, $
    NO_COLUMNS=no_columns_opt, $
    NOCOMPRESS=nocompressopt

;***********************************************************************

  if (filename eq '') then begin
    Message,/INFO,'Filename is a null string.'
    Return,''
  endif

; evaluate the keywords 
  if (not Keyword_Set(separator)) then sep = ' ' else sep = separator
  all = 0
  if (not Keyword_Set(verbose)) then verbose = 0 else verbose = 1
  if (not Keyword_Set(skip_lines)) then skip_lines = 0
  no_columns_opt = keyword_set( no_columns_opt )

; Default comment character flag
if n_elements( comment_char ) eq 0 then comment_char = ';'

; does the file exist?
  a = file_search( filename, COUNT = count )
  if (count eq 0) then begin
    Print,' SREAD ERROR: at read.pro, FILE NOT FOUND. '+filename
    Return, ''
  endif else begin
;   Open the file and read the first (valid) line
    OpenR,lun,filename,/GET_LUN
    a = FStat(lun)
    if (a.size gt 0) then begin
      s = ''
      if (Keyword_Set(skip_lines)) then begin
        skip_lines = Round(skip_lines)        
        for i=1,skip_lines do ReadF,lun,s
      endif
  
;     Get the first line (not a comment) and parse it to figure out how many 
;     columns there are
      stop = 0
      while ((not EOF(lun)) and (not stop)) do begin
        ReadF,lun,s
        if not( keyword_set( nocompressopt ) ) then begin
          s = StrCompress(StrTrim(s,2),REMOVE_ALL=all)
        endif
        if ( StrPos( s, comment_char ) lt 0 ) and ( s ne '' ) then begin
          stop = 1
        endif
      endwhile
;     if the first valid line has zero length then return a null string
      if (StrLen(s) eq 0) then begin
        if (verbose) then Print,'SREAD: no valid data in '+filename+'.'
        Return, ''
      endif

;     Begin the process of reading the data.

      if no_columns_opt eq 0 then begin
        s = strsplit( s, sep, extract=1, preserve_null=preserve_null_opt )  
  ;     get the number of elements in the first valid line  
        n         = (Size(s))[1]
      endif else begin
        n = 1
      endelse

      line      = ''
      start_len = Long(1000)
      inc_len   = Long(1000)
      tot_len   = start_len
      tmp       = StrArr(n,start_len)
      cnt       = Long(0)
    
      if (verbose) then Print,'Opening '+filename
   
      Point_Lun,lun,0  
      if (Keyword_Set(skip_lines)) then begin
;       Check that skip_lines is an integer 
        s = ''
        if (Var_Type(skip_lines) eq 2) then begin
          for i=1,skip_lines do ReadF,lun,s
        endif
      endif
   
      inc_cnt = 0
      wrote_message = 0
;     read to the end of the file
      while (not EOF(lun)) do begin
        s = ''
        ReadF,lun,line          
        if not( keyword_set( nocompressopt ) ) then begin
          line = StrCompress(StrTrim(line,2),REMOVE_ALL=all) 
        endif
        if ( StrPos( line, comment_char ) lt 0 ) then begin
          if no_columns_opt eq 0 then begin
            ; If there is a protector character defined
            if n_elements( protector ) ne 0 then begin
              ; Determine the locations of the protector character
              pos_protect = strsplit( line, protector, count=n_pos_protect, $
                  preserve_null=1 )
              ; If there are no occurences
              if n_pos_protect le 1 then begin
                ; The proceed to split normally into columns
                sepstring = strsplit( line, sep, extract=1, $
                    preserve_null=preserve_null_opt )
              ; If there are occurrences of the protector character
              endif else begin
                ; Remove the position markers for the start and end of the 
                ; string if the protector character does not actually occur 
                ; there
                if strmid( line, 0, 1 ) ne protector then begin
                  id = where( pos_protect ne 0, n_pos_protect )
                  if n_pos_protect eq 0 then stop
                  pos_protect = pos_protect[id]
                endif
                temp_len = strlen( line )
                if strmid( line, temp_len - 1, 1 ) ne protector then begin
                  id = where( pos_protect ne temp_len, n_pos_protect )
                  if n_pos_protect eq 0 then stop
                  pos_protect = pos_protect[id]
                endif
                ; Find the locations of the separator string
                pos_sep = strsplit( line, sep, $
                    preserve_null=preserve_null_opt, count=n_pos_sep )
                ; Iterate through columns
                for i_pos_sep = 0, n_pos_sep - 1 do begin
                  ; Determine if this is between two instances of the protector 
                  ; character
                  id_below = where( pos_protect lt pos_sep[i_pos_sep], $
                      n_id_below )
                  id_above = where( pos_protect gt pos_sep[i_pos_sep], $
                      n_id_above )
                  if ( n_id_below gt 0 ) and ( n_id_above gt 0 ) then begin
                    id_below = id_below[n_id_below-1]
                    id_above = id_above[0]
                    ; If the instances immediately below and above form a pair
                    if odd( id_below ) eq 0 then begin
                      if id_above eq id_below + 1 then begin
                        ; Then flag this separator location for removal
                        pos_sep[i_pos_sep] = -1
                      endif
                    endif
                  endif
                endfor
                ; Remove cases of protected separator strings
                id = where( pos_sep ne -1, n_pos_sep )
                pos_sep = [ pos_sep[id], strlen( line ) + 1 ]
                ; Extract columns
                sepstring = strmid( line, pos_sep[0:n_pos_sep-1], $
                    pos_sep[1:n_pos_sep] - pos_sep[0:n_pos_sep] - 1 )
              endelse
            ; Otherwise extract columns normally
            endif else begin
              sepstring = strsplit( line, sep, extract=1, $
                  preserve_null=preserve_null_opt )
            endelse
          endif else begin
            sepstring = line
          endelse
          if (keyword_Set(debug)) then Print,n,size(sepstring),cnt
          if (Keyword_Set(truncate)) then begin
            tmp[*,cnt] = sepstring[0:n-1]
          endif else begin
            if strlen( line ) eq 0 then begin
              nn = 1
            endif else begin
              nn = n_elements( sepstring )
            endelse
            if (nn gt n) then begin
;             Since we encountered a longer line (a line with more
;             elements than we expected) increase that dimension of
;             temp and don't forget n.      
              if (verbose) then Print,'SREAD:  Found a longer line;' $
                                     +' adjusting array.'
              tmp2 = StrArr(nn,tot_len)
              tmp2[0:n-1,0:cnt-1] = tmp[*,0:cnt-1]
              n    = nn
              tmp  = tmp2
            endif
            tmp[0:nn-1,cnt] = sepstring
          endelse
          cnt = cnt + 1
        endif
;       increase the size of the arrays to account for more lines 
        if (cnt gt tot_len-1) then begin 
          tot_len = tot_len + inc_len
          tmp2 = StrArr(n,tot_len)
          tmp2[0:n-1,0:cnt-1] = tmp
          tmp = tmp2
          inc_cnt = inc_cnt +1
          if (inc_cnt ge 10) then begin
            inc_len = inc_len*2
            inc_cnt = 0
            if (not wrote_message) then begin
              if (verbose) then Print,"This file appears to be quite long."
              wrote_message = 1
            endif
            ;if (verbose) then Print,"I'm increasing the read increment to " $
            ;                       +StringC(inc_len)+" lines."
            if verbose then begin
              Print, "I'm increasing the read increment to " + str( inc_len ) $
                  + " lines."
            endif
          endif
        endif
      endwhile
 
      Free_lun,lun

      if (verbose) then begin 
        ;Print,'Read '+StringC(cnt)+' lines from '+filename+'.' 
        Print, 'Read ' + str( cnt ) + ' lines from ' + filename + '.' 
      endif
      
;     Now reduce the array if we made it too big
      if (n eq 1) then begin
        tmp2 = StrArr(cnt)
        tmp2[*] = tmp[0,0:cnt-1]
      endif else begin
        tmp2 = StrArr(n,cnt)
        tmp2 = tmp[0:n-1,0:cnt-1]   
      endelse

      Return,tmp2
 
    endif else begin
      Print,'File '+filename+' has zero length.'
      Return, ''
    endelse

  endelse
end
