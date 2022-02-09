;+
; NAME:
;    IDL_LIB_HTML
;
; PURPOSE:
;    This procedure creates the web page listing of a library's IDL utilities.
;
; CATEGORY:
;    Miscellaneous
;
; CALLING SEQUENCE:
;    IDL_LIB_HTML
;
; KEYWORD PARAMETERS:
;    ADMIN_WEB:  The administrator's web page.  The default is defined in the 
;        Constants section.
;    ADMIN_NAME:  The administrator's name.  The default is defined in the 
;        Constants section.
;    ADMIN_EMAIL:  The administrator's e-mail address.  The default is defined 
;        in the Constants section.
;    CATEG:  An optionsl string vector listing the function and procedure 
;        categories to which to restrict.
;    COLOR_TABLE:  An optional string describing the colour table to use.  
;        Possible values are "black on beige" and "black on cyan and blue".
;    FIRST_UPDATE:  A vector of strings containing information on updates to be 
;        added before the list of routines.
;    HEAD_META:  An optional string array of size 2*N_HEAD_META containing the 
;        entries for the "<META" fields within the "<HEAD>" field.  The values 
;        in [0,*] denote the entries for the "NAME=" variable, while the [1,*] 
;        values denote the entries for the "CONTENT=" variable. 
;    INDIR:  A scalar or vector string containing the directories in which to 
;        search for files.  The default is the current directory.
;    INFILE:  The flag for filenames to include in the library.  The default is 
;        '*.pro'.
;    MODIFIED_AUTHOR:  If set then authors under the "Modified" heading are 
;        included as Contributors in the output webpage.
;    OUTFILE:  The filename for the output webpage.  The default is 
;        'idl_lib.html'.
;    PACKAGE_CATEG:  An optional string array of length N_PACKAGE providing a 
;        list of all categories to include the package collections listed in 
;        PACKAGE_NAME.  For each package, multiple categories can be included 
;        in a semi-colon-delimited format.  The code traces dependences in the 
;        routines included in these categories, so the final list of files 
;        included may include routines outside of the specified categories.
;    PACKAGE_INFO:  An optional string array of length N_PACKAGE providing 
;        descriptions of the packages listed in PACKAGE_NAME, to be included 
;        in the package listings on the web page.
;    PACKAGE_NAME:  An optional string array of length N_PACKAGE listing names 
;        of package collections.  These collections are listed at the top of 
;        the page with a description (see PACKAGE_INFO) and tar and zip 
;        archives of all necessary code for usage of the package.
;    REVERSE_CATEGORY:  Order the categories in the output webpage in reverse 
;        alphabetical order.  The default is alphabetical order.
;    SUPERTITLE_PAGE:  The optional super-title of the web page.
;    TABS_HREF:  An optional string vector of length N_TABS containing the URLs 
;        for the tabs appearing at the top of the page.
;    TABS_NAME:  An optional string vector of length N_TABS containing the 
;        names of the tabs appearing at the top of the page.
;    TABS_WIDTH:  An optional integer vector of length N_TABS containing the 
;        width of the tabs appearing at the top of the page.
;    TAR:  If set, then a .tar archive file is created containing all of the 
;        routines in the library.  This archive will be included on the web 
;        page.
;    TITLE_FIRST_UPDATE:  A string contain the title of the FIRST_UPDATE list 
;        of information.  The default is 'UPDATES'.
;    TITLE_PAGE:  The title of the web page.
;    TITLE_UPDATE:  A string contain the title of the UPDATE list of 
;        information.  The default is 'UPDATES'.
;    UPDATE:  A vector of strings containing information on updates to be added 
;        after the list of routines.
;    ZIP:  If set, then a .zip archive file is created containing all of the 
;        routines in the library.  This archive will be included on the web 
;        page.
;
; USES:
;    isin.pro
;    month_name.pro
;    str.pro
;    string_from_vector.pro
;    string_substitute.pro
;    trace_dependency.pro
;
; PROCEDURE:
;    This procedure reads the documentation of the IDL files in the group 
;    utility directory, and creates a web page that lists them.
;
; MODIFICATION HISTORY:
;    Written by:  Daithi A. Stone (stoned@atm.ox.ac.uk), 2002-02-08.
;    Modified:  DAS, 2002-04-10 (STD.pro removal update)
;    Modified:  DAS, 2002-08-09 (CONSTANTS.pro changes update)
;    Modified:  DAS, 2002-11-29 (altered directory structure)
;    Modified:  DAS, 2002-12-06 (removed dependency on me)
;    Modified:  DAS, 2003-05-07 (renamed, and added ADMIN* keywords)
;    Modified:  DAS, 2004-06-24 (added old modification date printing)
;    Modified:  DAS, 2005-09-01 (added FIRST_UPDATE, FUTITLE, INDIR, INFILE, 
;        MODIFIED_AUTHOR, OUTFILE, REVERSE_CATEGORY, TITLE, UPDATE, UTITLE 
;        keywords;  permitted multiple input directories;  removed use of 
;        constants.pro)
;    Modified:  DAS, 2006-03-23 (corrected self-reference;  increased 
;        modification highlighting to 90 days)
;    Modified:  DAS, 2007-05-24 (added TAR and ZIP options)
;    Modified:  DAS, 2008-12-01 (updated administrator details)
;    Modified:  DAS, 2011-11-06 (modified documentation;  altered keyword 
;        parameter names;  modified colours;  added CATEG keyword;  introduced 
;        removal of middle initials in names;  allowed multiple categories per 
;        file)
;    Modified:  DAS, 2013-07-19 (switched from use of findfile to file_search)
;    Modified:  DAS, 2018-08-17 (added COLOR_TABLE keyword input, switched to 
;        "black on cyan and blue" default colour table;  add HEAD_META, 
;        TABS_HREF, TABS_NAME, TABS_WIDTH, PACKAGE_* keywords)
;    Modified:  DAS, 2018-09-18 (Made modification date detection more robust;  
;        Completed documentation for PACKAGE_* keywords)
;-

;***********************************************************************

PRO IDL_LIB_HTML, $
    ADMIN_WEB=admin_web, ADMIN_NAME=admin_name, ADMIN_EMAIL=admin_email, $
    CATEG=categ_0, $
    COLOR_TABLE=COLOR_TABLE, $
    FIRST_UPDATE=first_update, TITLE_FIRST_UPDATE=title_first_update, $
    HEAD_META=head_meta, $
    INDIR=indir, INFILE=infile, $
    MODIFIED_AUTHOR=modified_author_opt, $
    OUTFILE=outfile, $
    PACKAGE_NAME=package_name, PACKAGE_CATEG=package_categ, $
      PACKAGE_INFO=package_info, $
    REVERSE_CATEGORY=reverse_category_opt, $
    TABS_HREF=tabs_href, TABS_NAME=tabs_name, TABS_WIDTH=tabs_width, $
    TAR=tar_opt, ZIP=zip_opt, $
    TITLE_PAGE=title_page, SUPERTITLE_PAGE=SUPERTITLE_PAGE, $
    UPDATE=update, TITLE_UPDATE=title_update

;***********************************************************************
; Constants and variables

; Input directory
if not( keyword_set( indir ) ) then indir = ''
if not( keyword_set( infile ) ) then infile = '*.pro'
; Output file
if not( keyword_set( outfile ) ) then outfile = 'idl_lib.html'

; The time for mentioning modifications (in days)
hist_delay = 90

; Maximum number of categories per routine
n_categ_max = 5

; The number of requested packages
n_package = n_elements( package_name )
if n_elements( package_categ ) ne n_package then stop
; The default empty package description
if n_package gt 0 then begin
  if n_elements( package_info ) eq 0 then package_info = strarr( n_package )
endif

; The default title for the update sections
if not( keyword_set( title_first_update ) ) then title_first_update = 'UPDATES'
if not( keyword_set( title_update ) ) then title_update = 'UPDATES'

; Default colours and text style
if not( keyword_set( color_table ) ) then color_table = 'black on cyan and blue'
; The "black on beige" colour table
if color_table eq 'black on beige' then begin
  ; Off-page background colour
  color_offpage = '#222222'
  ; Page background and border colours
  color_page_background = '#bbbb99'
  color_page_border = '#bbbb99'
  ; Text colour
  color_text = '#222222'
  ; Header colour
  color_header = '#aa0000'
  ; Link colour
  color_link = '#0033aa'
  ; Viewed link colour
  color_viewedlink = '#0000aa'
  ; Text size
  size_text = '"4"'
endif else if color_table eq 'black on cyan and blue' then begin
  ; Off-page background colour
  color_offpage = '#222222'
  ; Page background and border colours
  color_page_background = '#36649c'
  color_page_border = '#298f8f'
  ; Window background and border colours
  color_box_background = '#b7d7df'
  color_box_border = '#005959'
  ; Text colours
  color_text = '#000000'
  color_text_shadow = '#97b7bf'
  ; Header colours
  color_header = '#1f278a'
  color_header_shadow = '#678fbf'
  ; Title colours
  color_supertitle = ' #ffcd7f'
  color_title = '#eeac44'
  color_title_shadow = '#021832'
  ; Link colour
  color_link = '#4d2f00'
  ; Viewed link colour
  color_viewedlink = '#4d2f00'
  ; Tab colours
  color_tab_background = '#5ab5b5'
  color_tab_text = '#083061'
  ; Footer box colors
  color_footer_background = '#0c7878'
  color_footer_text = '#ffffff'
  color_footer_shadow = '#002e2e'
  ; Text font
  font_text = 'Arial,Helvetica'
  ; Header font
  font_header = font_text
  ; Text size
  size_text = '"4"'
endif else begin
  stop
endelse

; The default width of the page and boxes
width_page = '1000'
width_margin = '10'
width_box_head = '800'
width_head = str( fix( width_box_head ) - 4 * fix( width_margin ) )
width_box_text = str( fix( width_page ) - 2 * fix( width_margin ) )
width_text = str( fix( width_box_text ) - 4 * fix( width_margin ) )


; Months in a year
mina=12
; Date and time
temp_time = strsplit( systime(), ' ', extract=1 )
temp_month_abbrev = month_name( indgen( mina ), abbreviate=1 )
id_month = where( temp_month_abbrev eq temp_time[1], n_id )
if n_id ne 1 then stop
date_current_numeric = temp_time[4] $
    + str( id_month[0] + 1, length=2, filler='0' ) + temp_time[2]
date_current_abbrev = temp_time[4] + '-' $
    + str( id_month[0] + 1, length=2, filler='0' ) + '-' + temp_time[2]
date_current_written = str( fix( temp_time[2] ) ) + ' ' $
    + month_name( id_month ) + ' ' + temp_time[4]

;***********************************************************************
; Set-up for reading IDL files

; Find files
list = file_search( indir[0]+infile, count=n_pro )
n_indir = n_elements( indir )
if n_indir ne 1 then begin
  for i_indir = 1, n_indir - 1 do begin
    list = [ list, file_search( indir[i_indir]+infile, count=temp ) ]
    n_pro = n_pro + temp
  endfor
endif

; Program information vectors
pro_name = strarr( n_pro )
pro_purp = strarr( n_pro )
pro_categ = strarr( n_pro, n_categ_max )
pro_hist = strarr( n_pro )
pro_date = lonarr( n_pro )
pro_auth = strarr( n_pro )
;Number of variables to read
n_read = 4

;***********************************************************************
; Read IDL files

; Set temporary variable, counters
temp_str = ''
check_name = 0
check_purp = 0
check_categ = 0
check_hist = 0

; Iterate through routines
for i_pro = 0, n_pro-1 do begin
  ; Open file
  openr, 1, list[i_pro]
  ; Set counters
  check = 0
  ; Run loop until all information retrieved or end of documentation
  while check le n_read + 1 do begin

    ; Read line from file
    readf, 1, temp_str
    temp_str = strtrim( strcompress( temp_str ), 2 )
    if check eq 0 then begin
      if temp_str eq ';+' then check = 1
    endif

    ; If in documentation section
    if check ge 1 then begin
      ; Check for end of documentation
      if temp_str eq ';-' then check = n_read + 2

      ; Read name
      if check_name ne 0 then begin
        ;pro_name[i_pro] = strupcase( $
        ;     strtrim( strmid( temp_str, 1, strlen( temp_str ) - 1 ), 2 ) )
        pro_name[i_pro] $
             = strtrim( strmid( temp_str, 1, strlen( temp_str ) - 1 ), 2 )
        check_name = 0
        check = check + 1
        ; Remove ".pro"
        pos = strpos( pro_name[i_pro], '.PRO' )
        if pos ne -1 then pro_name[i_pro] = strmid( pro_name[i_pro], 0, pos )
      endif

      ; Read purpose
      if check_purp ne 0 then begin
        if strlen( temp_str ) gt 3 then begin
          pro_purp[i_pro] = pro_purp[i_pro] $
              + strmid( temp_str, 1, strlen( temp_str ) - 1 )
        endif else begin
          pro_purp[i_pro] = strtrim( strcompress( pro_purp[i_pro] ), 2 )
          check_purp = 0
          check = check + 1
        endelse
      endif

      ; Read category
      if check_categ ne 0 then begin
        ; Extract category
        temp_str = strupcase( $
            strtrim( strmid( temp_str, 1, strlen( temp_str ) - 1 ), 2 ) )
        ; Check for multiple categories
        temp = strtrim( strsplit( temp_str, ',', extract=1, count=n_temp ), 2 )
        pro_categ[i_pro,0:n_temp-1] = temp
        check_categ = 0
        check = check + 1
      endif

      ; Read history
      if check_hist ne 0 then begin
        if strlen( temp_str ) gt 3 then begin
          ; Get author
          pos = strpos( temp_str, 'Written by:' )
          if pos ne -1 then begin
            pos_1 = strpos( temp_str, ', ' )
            pos_2 = strpos( temp_str, ' (' )
            if pos_2 ne -1 then pos_1 = min( [ pos_1, pos_2 ] )
            pro_auth[i_pro] = strmid( temp_str, pos+12, pos_1-pos-12 )
          endif
          ; Get modifying authors if requested
          if keyword_set( modified_author_opt ) then begin
            pos = strpos( temp_str, 'Modified:' )
            if pos ne -1 then begin
              pos_1 = strpos( temp_str, ', ' )
              pos_2 = strpos( temp_str, ' (' )
              if pos_2 ne -1 then pos_1 = min( [ pos_1, pos_2 ] )
              ; Only save this author name if it is not initials because 
              ; initials are often used as a shorthand for the original author 
              ; who will be in PRO_AUTH anyway
              temp_1 = strmid( temp_str, pos+10, pos_1-pos-10 )
              if strupcase( temp_1 ) ne temp_1 then begin
                if keyword_set( mod_auth ) then begin
                  mod_auth = [ mod_auth, temp_1 ]
                endif else begin
                  mod_auth = temp_1
                endelse
              endif
            endif
          endif
          ; Get more recent date
          pos = strpos( temp_str, '-' )
          if pos ne -1 then begin
            pos_1 = strpos( strmid( temp_str, pos+1, strlen( temp_str )-pos ), $
                '-' )
            if pos_1 eq 2 then begin
              temp = long( strmid( temp_str, pos-4, 4 ) $
                  + strmid( temp_str, pos+1, 2 ) $
                  + strmid( temp_str, pos+4, 2 ) )
              if temp gt pro_date[i_pro] then begin
                pro_date[i_pro] = temp
                check_hist = check_hist + 1
              endif
            endif
          endif
        endif else begin
          ; Create modification alert
          if check_hist eq 2 then begin
            pro_hist[i_pro] = 'Written '
          endif else begin
            pro_hist[i_pro] = 'Modified '
          endelse
          pro_hist[i_pro] = pro_hist[i_pro] $
              + strmid( str( pro_date[i_pro] ), 0, 4 ) + '-' $
              + strmid( str( pro_date[i_pro] ), 4, 2 ) + '-' $
              + strmid( str( pro_date[i_pro] ), 6, 2 )
          check_hist = 0
          check = check + 1
        endelse
      endif

      ; Search for start of documentation fields
      if temp_str eq '; NAME:' then check_name = 1 
      if temp_str eq '; PURPOSE:' then check_purp = 1 
      if temp_str eq '; CATEGORY:' then check_categ = 1 
      if temp_str eq '; MODIFICATION HISTORY:' then check_hist = 1 

    endif
  endwhile

  ; Close file
  close, 1
endfor

;***********************************************************************
; Sort Programs

; Category list
categ_name = ['']
categ_name_use = ['']
n_categ = 1
; Category classification of programs
pro_categ_index = intarr( n_pro, n_categ_max )

; Build categories
; Go through programs
for i_pro = 0, n_pro-1 do begin
  ; Iterage through possible multiple categories
  for i_categ_mult = 0, n_categ_max - 1 do begin
    if pro_categ[i_pro,i_categ_mult] ne '' then begin
      ; Set counter
      check = 0
      ; Compare program category to category list
      for i_categ = 0, n_categ-1 do begin
        if pro_categ[i_pro,i_categ_mult] eq categ_name[i_categ] then begin
          ; Link to existing category
          pro_categ_index[i_pro,i_categ_mult] = i_categ
          check = 1
        endif
      endfor
      ; Create new category and link
      if check eq 0 then begin
        ; But only if it is consistent with restricted list
        if keyword_set( categ_0 ) then begin
          id = where( strupcase( categ_0 ) eq pro_categ[i_pro,i_categ_mult], $
              n_id )
          if n_id gt 0 then begin
            categ_name = [ categ_name, pro_categ[i_pro,i_categ_mult] ]
            categ_name_use = [ categ_name_use, categ_0[id[0]] ]
            n_categ = n_categ + 1
            pro_categ_index[i_pro,i_categ_mult] = n_categ - 1
          endif
        endif else begin
          categ_name = [ categ_name, pro_categ[i_pro,i_categ_mult] ]
          n_categ = n_categ + 1
          pro_categ_index[i_pro,i_categ_mult] = n_categ - 1
        endelse
      endif
    endif
  endfor
endfor

; Adopt favoured format for category name
if keyword_set( categ_name_use ) then categ_name = temporary( categ_name_use )
; Remove initialisation value
categ_name = categ_name[1:n_categ-1]
n_categ = n_categ - 1
pro_categ_index = pro_categ_index - 1
; Sort categories
id_categ_sort = sort( categ_name )
; Reverse order if requested
if keyword_set( reverse_category_opt ) then begin
  id_categ_sort = reverse( id_categ_sort )
endif

; Build author list
auth_name = ['']
n_auth = 1
; Go through programs
for i_pro = 0, n_pro-1 do begin
  ; Remove middle initial
  temp = strsplit( pro_auth[i_pro], ' ', extract=1, count=n_temp )
  if n_temp gt 2 then begin
    for i_temp = 1, n_temp - 2 do begin
      if strlen( temp[i_temp] ) eq 1 then temp[i_temp] = ''
      if ( strlen( temp[i_temp] ) eq 2 ) $
          and ( strmid( temp[i_temp], 1, 1 ) eq '.' ) then temp[i_temp] = ''
    endfor
    id = where( temp ne '' )
    temp = temp[id]
  endif
  temp = string_from_vector( temp, spacer=' ', nospace=1 )
  ; Set counter
  check = 0
  ; Compare program author to author list
  for i_auth = 0, n_auth-1 do begin
    if temp eq auth_name[i_auth] then check = 1
  endfor
  ; Add new author
  if check eq 0 then begin
    auth_name = [ auth_name, temp ]
    n_auth = n_auth + 1
  endif
endfor
; Go through modifying authors if requested
if keyword_set( modified_author_opt ) then begin
  for i_mod = 0, n_elements( mod_auth ) - 1 do begin
    ; Remove middle initial
    temp = strsplit( mod_auth[i_mod], ' ', extract=1, count=n_temp )
    if n_temp gt 2 then begin
      for i_temp = 1, n_temp - 2 do begin
        if strlen( temp[i_temp] ) eq 1 then temp[i_temp] = ''
        if ( strlen( temp[i_temp] ) eq 2 ) $
            and ( strmid( temp[i_temp], 1, 1 ) eq '.' ) then temp[i_temp] = ''
      endfor
      id = where( temp ne '' )
      temp = temp[id]
    endif
    temp = string_from_vector( temp, spacer=' ', nospace=1 )
    ; Set counter
    check = 0
    ; Compare program author to author list
    for i_auth = 0, n_auth-1 do begin
      if temp eq auth_name[i_auth] then check = 1
    endfor
    ; Add new author
    if check eq 0 then begin
      auth_name = [ auth_name, temp ]
      n_auth = n_auth + 1
    endif
  endfor
endif
; Remove the initial entry (used for initialisation only)
auth_name = auth_name[1:n_auth-1]
n_auth = n_auth - 1
; Sort alphabetically by surname
temp_surname = strarr( n_auth )
for i_auth = 0, n_auth - 1 do begin
  temp = strsplit( auth_name[i_auth], ' ', extract=1 )
  temp_surname[i_auth] = temp[1]
endfor
id = sort( temp_surname )
auth_name = auth_name[id]

;***********************************************************************
; Creat Archive Files

;; Create .tar archive if requested
;if keyword_set( tar_opt ) then begin
;  ; Create archive file
;  spawn, 'tar -cf ' + indir[0] + 'idl_lib.tar ' $
;      + string_from_vector( list, spacer=' ' )
;  ; Allow access to everyone
;  spawn, 'chmod a+r ' + indir[0] + 'idl_lib.tar'
;endif
;
;; Create .zip archives if requested
;if keyword_set( zip_opt ) then begin
;  ; Create archive file
;  spawn, 'zip ' + indir[0] + 'idl_lib.zip ' $
;      + string_from_vector( list, spacer=' ' )
;  ; Allow access to everyone
;  spawn, 'chmod a+r ' + indir[0] + 'idl_lib.zip'
;endif

; Sort out the packages
if n_package gt 0 then begin
  ; Initialise the list of archive files
  if keyword_set( tar_opt ) then package_file_tar = strarr( n_package )
  if keyword_set( zip_opt ) then package_file_zip = strarr( n_package )
  ; Initialise the lists of dependencies outside of the available library
  package_external = strarr( n_package )
  ; Iterate through packages
  for i_package = 0, n_package - 1 do begin
    ; Find files identified with the specified categories
    temp = isin( $
        strupcase( strsplit( package_categ[i_package], ';', extract=1 ) ), $
        strupcase( categ_name ) )
    id_categ = where( temp eq 1, n_id_categ )
    temp = isin( fix( id_categ ), fix( pro_categ_index ) )
    temp = reform( temp, n_pro, n_categ_max )
    if n_categ_max gt 1 then temp = max( temp, dimension=2 )
    id_pro = where( temp eq 1, n_id_pro )
    temp_file_list = list[id_pro]
    ; Remove redundancies
    id = uniq( temp_file_list )
    temp_file_list = temp_file_list[id]  
    n_temp_file_list = n_elements( temp_file_list )
    ; Determine where this collection should sit
    if n_elements( temp_file_list ) eq 1 then stop ; Not yet implemented
    ctr_same = -1
    check_flag = 0
    for i_char = 0, min( strlen( temp_file_list ) ) - 1 do begin
      if check_flag eq 0 then begin
        temp = strmid( temp_file_list, i_char, 1 )
        if n_elements( uniq( temp ) ) eq 1 then begin
          ctr_same = i_char
        endif else begin
          check_flag = 1
        endelse
      endif
    endfor
    if ctr_same eq 0 then begin
      dir_package = ''
    endif else begin
      dir_package = strmid( temp_file_list[0], 0, ctr_same + 1 )
      pos = strpos( dir_package, '/', reverse_search=1 )
      if pos eq -1 then begin
        dir_package = ''
      endif else begin
        dir_package = strmid( dir_package, 0, pos + 1 )
      endelse    
    endelse
    ; Split directory information from the file list
    temp_dir = temp_file_list
    temp_name = temp_file_list
    for i_file = 0, n_temp_file_list - 1 do begin
      temp = strsplit( temp_file_list[i_file], '/', extract=1, count=n_temp )
      if n_temp gt 1 then temp_dir[i_file] = strjoin( temp[0:n_temp-2], '/' )
      temp_name[i_file] = temp[n_temp-1]
    endfor
    ; Find dependencies to include
    for i_file = 0, n_temp_file_list - 1 do begin
      temp_code_name = ''
      trace_dependency, temp_name[i_file], no_plot=1, $
          code_name=temp_code_name, idl_path=temp_dir[i_file]
      id = where( isin( temp_name, temp_code_name ) eq 0, n_id )
      if n_id gt 0 then begin
        temp_name = [ temp_name, temp_code_name[id] ]
        temp_file_list = [ temp_file_list, strarr( n_id ) ]
      endif
    endfor
    n_temp_file_list = n_elements( temp_name )
    ; Locate the dependences within the library
    id_missing = where( temp_file_list eq '', n_id_missing )
    for i_missing = 0, n_id_missing - 1 do begin
      id_sub = where( $
          strpos( '/' + list, '/' + temp_name[id_missing[i_missing]] ) ge 0, $
          n_id_sub )
      ; If there are multiple routines with this name
      if n_id_sub gt 1 then begin
        ; Take the one in the category that appears first in the input category 
        ; list
        id_categ = where( isin( strupcase( pro_categ[id_sub] ), $
            strupcase( categ_name ) ) eq 1, n_id_categ )
        if n_id_categ ne n_id_sub then stop
        id = where( strupcase( pro_categ[id_sub] ) $
            eq strupcase( categ_name[id_categ[0]] ), n_id_sub )
        if n_id_sub ne 1 then stop
        id_sub = id_sub[id]
      endif
      if n_id_sub eq 1 then begin
        temp_file_list[id_missing[i_missing]] = list[id_sub]
      endif else if n_id_sub eq 0 then begin
        if package_external[i_package] eq '' then begin
          package_external[i_package] = temp_name[id_missing[i_missing]]
        endif else begin
          package_external[i_package] = package_external[i_package] + ', ' $
              + temp_name[id_missing[i_missing]]
        endelse
      endif else begin
        stop
      endelse
    endfor
    ; Remove any external dependences from the list
    id = where( temp_file_list ne '', n_id )
    if n_id eq 0 then stop
    temp_file_list = temp_file_list[id]
    ; Create the package file name
    temp_package_name = string_substitute( package_name[i_package], ' ', '_' )
    ; Create .tar archive if requested
    if keyword_set( tar_opt ) then begin
      ; Create archive file
      package_file_tar[i_package] = dir_package + temp_package_name + '.tar'
      spawn, 'tar -cf ' + package_file_tar[i_package] + ' ' $
          + string_from_vector( temp_file_list, spacer=' ' )
      ; Allow access to everyone
      spawn, 'chmod a+r ' + package_file_tar[i_package]
    endif
    ; Create .zip archives if requested
    if keyword_set( zip_opt ) then begin
      ; Create archive file
      package_file_zip[i_package] = dir_package + temp_package_name + '.zip'
      spawn, 'zip ' + package_file_zip[i_package] + ' ' $
          + string_from_vector( temp_file_list, spacer=' ' )
      ; Allow access to everyone
      spawn, 'chmod a+r ' + package_file_zip[i_package]
    endif
  endfor
endif

;***********************************************************************
; Create Web Page

; Open web page file
openw, 1, outfile

; Print comment header
printf, 1, '<!-- *********************************************** /-->'
printf, 1, '<!-- ********** ' + title_page + ' ********** /-->'
printf, 1, '<!-- * ' + admin_name + ' (' + admin_email + ') * /-->'
printf, 1, '<!-- * This file was produced by idl_lib_html.pro. * /-->'
printf, 1, '<!-- * ' + date_current_abbrev + ' * /-->'
printf, 1, '<!-- *********************************************** /-->'
printf, 1

; Web page header
printf, 1, '<HTML>'
printf, 1, '<HEAD>'
printf, 1, '  <TITLE>' + title_page + '</TITLE>'
if keyword_set( head_meta ) then begin
  for i_meta = 0, n_elements( head_meta[0,*] ) - 1 do begin
    printf, 1, '  <META NAME="' + head_meta[0,i_meta] + '" CONTENT="' $
        + head_meta[1,i_meta] + '">'
  endfor
endif
printf, 1, '  <STYLE TYPE="text/css">'
printf, 1, '    A { text-decoration: none }'
printf, 1, '  </STYLE>'
printf, 1, '</HEAD>'

; Colours and fonts
printf, 1, '<BODY BGCOLOR="' + color_offpage + '" LINK="' + color_link $
    + '" VLINK="' + color_viewedlink + '" ALINK="' + color_link + '" TEXT="' $
    + color_text + '">'

; Begin table containing web page
printf, 1, '<!-- Table containing web page /-->'
printf, 1, '<CENTER>'
printf, 1, '<TABLE WIDTH=' + width_page + ' STYLE="border:5px outset ' + color_page_border + ';" BGCOLOR="' + color_page_background + '">'

; Links tabs
if keyword_set( tabs_name ) then begin
  printf, 1, '  <!-- Links row /-->'
  printf, 1, '  <TR>'
  printf, 1, '    <TD VALIGN="top" WIDTH=' + width_page + ' ALIGN="center">'
  printf, 1, '      <TABLE>'
  printf, 1, '        <TR>'
  for i_tabs = 0, n_elements( tabs_name ) - 1 do begin
    temp = ''
    if keyword_set( tabs_width ) then temp = str( tabs_width[i_tabs] )
    printf, 1, '          <TD WIDTH=' + temp $
        + ' ALIGN="center" STYLE="border:5px outset ' + color_box_border + ';" BGCOLOR="' + color_tab_background + '">'
    printf, 1, '          <A HREF="' + tabs_href[i_tabs] $
        + '" style="font-family : ' + font_header + '; color : ' + color_tab_text + '; font-size : large; font-weight : bold; font-variant:small-caps; text-shadow: 1px 1px 2px ' + color_page_border + ';">'
    printf, 1, '          ' + tabs_name[i_tabs]
    printf, 1, '        </A>'
    printf, 1, '      </TD>'
  endfor
  printf, 1, '        </TR>'
  printf, 1, '      </TABLE>'
  printf, 1, '    </TD>'
  printf, 1, '  </TR>'
endif

; Print title row
printf, 1, '  <!-- Title row /-->'
printf, 1, '  <TR>'
printf, 1, '    <TD WIDTH=' + width_page + ' STYLE="font-size:normal;">'
printf, 1, '      <BR>'
printf, 1, '    </TD>'
printf, 1, '  </TR>'
if keyword_set( supertitle_page ) then begin
  printf, 1, '  <TR>'
  printf, 1, '    <TD WIDTH=' + width_page + ' ALIGN="right" style="font-family: ' + font_header + '; color : ' + color_supertitle + '; font-size : xx-large; font-weight : bold; text-shadow: 2px 2px 2px ' + color_title_shadow + ';">'
  printf, 1, '      ' + supertitle_page
  printf, 1, '    </TD>'
  printf, 1, '  </TR>'
endif
printf, 1, '  <TR>'
printf, 1, '    <TD WIDTH=' + width_page + ' ALIGN="right" style="font-family : ' + font_header + '; color : ' + color_title + '; font-size : xx-large; font-weight : bold; text-shadow: 2px 2px 2px ' + color_title_shadow + ';">'
printf, 1, '      ' + title_page + ' &nbsp &nbsp'
printf, 1, '    </TD>'
printf, 1, '  </TR>'
printf, 1, '  <TR>'
printf, 1, '    <TD WIDTH=' + width_page + ' STYLE="font-size:small;">'
printf, 1, '      <BR>'
printf, 1, '    </TD>'
printf, 1, '  </TR>'
printf, 1, ''

printf, 1, '  <!-- Table containing general information /-->'
printf, 1, '  <TR>'
printf, 1, '    <TD ALIGN="center">'
printf, 1, '      <TABLE WIDTH=' + width_box_head + ' ALIGN="center">'
printf, 1, '        <TR>'
printf, 1, '          <TD VALIGN="top" WIDTH=' + width_box_head + ' ALIGN="center" STYLE="border:5px outset ' + color_box_border + '" BGCOLOR="' + color_box_background + '">'
printf, 1, '            <TABLE WIDTH=' + width_head + ' ALIGN="center">'
printf, 1, '              <TR>'
printf, 1, '                <TD ALIGN="left" WIDTH=' + width_head + ' style="font-family:' + font_header + '; color:' + color_text + '; text-shadow: 1px 1px 1px ' + color_text_shadow + ';">'
printf, 1, '                  <BR>'
temp = '<FONT COLOR="' + color_header + '"><B>This page lists a compilation of IDL code made available for public use.</B></FONT><BR>'
printf, 1, '                  ' + temp
temp = 'Please report any bugs to <A HREF="mailto:' + admin_email + '"><B>' $
    + admin_email + '</B></A>.<BR>'
printf, 1, '                  ' + temp
printf, 1, 'We welcome any modified or new routines you would like to add.  '
temp = 'We would like to know if you have found any of these routines ' $
    + 'useful.<BR>'
printf, 1, '                  ' + temp + '<BR>'
temp = '<FONT COLOR="' + color_header + '"><B>Maintained by:</B></FONT>  ' $
    + admin_name + ' (' + '<A HREF="mailto:' + admin_email + '"><B>' $
    + admin_email + '</B></A>)<BR>'
printf, 1, '                  ' + temp
temp = '<FONT COLOR="' + color_header + '"><B>Last modified:</B></FONT>  ' $
    + date_current_written + ' by <A HREF="' $
    + admin_web + '/idl_lib/pro/idl_lib_html.pro"><B>idl_lib_html.pro</B></A>.<BR>'
printf, 1, '                  ' + temp
temp = '<FONT COLOR="' + color_header + '"><B>Number of routines:</B></FONT>  ' $
    + str( n_pro ) + '<BR>'
printf, 1, '                  ' + temp
temp = '<FONT COLOR="' + color_header + '"><B>Contributors:</B></FONT>  ' $
    + string_from_vector( auth_name ) + '<BR>'
printf, 1, '                  ' + temp
temp = '<FONT COLOR="' + color_header $
    + '"><B>Licence:</B></FONT> The IDL routines available from this page are free ' $
    + 'for non-commercial use under the terms of this ' $
    + '<A HREF="http://creativecommons.org/licenses/by-nc-sa/2.0/"><B>Creative ' $
    + 'Commons License</B></A> unless otherwise noted in the routine.<BR><BR>' 
printf, 1, '                  ' + temp
; Write information to page
;if keyword_set( tar_opt ) or keyword_set( zip_opt ) then begin
;  temp = 'The entire routine library can be downloaded as an archive file.<BR>'
;  printf, 1, '                  ' + temp
;  if keyword_set( tar_opt ) then begin
;    temp = '  Click <A HREF="' + indir[0] $
;        + 'idl_lib.tar">here</A> for a .tar archive.<BR>'
;  endif
;  if keyword_set( zip_opt ) then begin
;    temp = '  Click <A HREF="' + indir[0] $
;        + 'idl_lib.zip">here</A> for a .zip archive.<BR>'
;  endif
;  printf, 1, '                  ' + temp
;  printf, 1, '                  <BR>'
;endif
printf, 1, '                </TD>'
printf, 1, '              </TR>'
printf, 1, '            </TABLE>'
printf, 1, '          </TD>'
printf, 1, '        </TR>'
printf, 1, '      </TABLE>'
printf, 1, '    </TD>'
printf, 1, '  </TR>'

printf, 1, '  <!-- Table containing links to the category lists /-->'
printf, 1, '  <TR>'
printf, 1, '    <TD ALIGN="center">'
printf, 1, '      <TABLE WIDTH=' + width_box_head + ' ALIGN="center">'
printf, 1, '        <TR>'
printf, 1, '          <TD VALIGN="top" WIDTH=' + width_head + ' ALIGN="center" STYLE="border:5px outset ' + color_box_border + '" BGCOLOR="' + color_box_background + '">'
printf, 1, '            <TABLE WIDTH=' + width_head + ' ALIGN="center">'
printf, 1, '              <TR>'
printf, 1, '                <TD>'
printf, 1, '                  <P STYLE="font-size:small;">'
printf, 1, '                    <BR>'
printf, 1, '                  </P>'
printf, 1, '                </TD>'
printf, 1, '              </TR>'
printf, 1, '              <TR>'
printf, 1, '                <TD ALIGN="center" WIDTH=' + width_head + ' style="font-family:' + font_header + '; color:' + color_header + '; font-size:x-large; font-weight:bold; text-shadow: 1px 1px 2px ' + color_header_shadow + ';">'
printf, 1, '                  <B>Links to categories</B><BR>'
printf, 1, '                </TD>'
printf, 1, '              </TR>'
printf, 1, '              <TR>'
printf, 1, '                <TD ALIGN="center" WIDTH=' + width_head + ' style="font-family:' + font_text + '; color:' + color_text + '; text-shadow: 1px 1px 1px ' + color_text_shadow + ';">'
temp_done = ''
for i_categ = 0, n_categ-1 do begin
  temp = string_substitute( categ_name[id_categ_sort[i_categ]], ' ', '_' )
  temp_done = temp_done + temp
  if strlen( temp_done ) gt 42 then begin
    temp_done = ''
    printf, 1, '                  <BR><BR>'
  endif
  printf, 1, '                  <A HREF="#' + temp + '" style="font-family:' + font_text + '; color:' + color_link + '; text-shadow: 1px 1px 1px ' + color_text_shadow + '; border:3px outset ' + color_box_border + '; background-color:' + color_tab_background + ';"><B>&nbsp;' + categ_name[id_categ_sort[i_categ]] + '&nbsp;</B></A>'
endfor
printf, 1, '                  <BR><BR>'
printf, 1, '                </TD>'
printf, 1, '              </TR>'
printf, 1, '            </TABLE>'
printf, 1, '          </TD>'
printf, 1, '        </TR>'
printf, 1, '      </TABLE>'
printf, 1, '    </TD>'
printf, 1, '  </TR>'
printf, 1

if n_package gt 0 then begin
  printf, 1, '  <!-- Table containing package lists /-->'
  printf, 1, '  <TR>'
  printf, 1, '    <TD ALIGN="center">'
  printf, 1, '      <TABLE WIDTH=' + width_box_head + ' ALIGN="center">'
  printf, 1, '        <TR>'
  printf, 1, '          <TD VALIGN="top" WIDTH=' + width_head + ' ALIGN="center" STYLE="border:5px outset ' + color_box_border + '" BGCOLOR="' + color_box_background + '">'
  printf, 1, '            <TABLE WIDTH=' + width_head + ' ALIGN="center">'
  printf, 1, '              <TR>'
  printf, 1, '                <TD>'
  printf, 1, '                  <P STYLE="font-size:small;">'
  printf, 1, '                    <BR>'
  printf, 1, '                  </P>'
  printf, 1, '                </TD>'
  printf, 1, '              </TR>'
  printf, 1, '              <TR>'
  printf, 1, '                <TD ALIGN="center" WIDTH=' + width_head + ' style="font-family:' + font_header + '; color:' + color_header + '; font-size:x-large; font-weight:bold; text-shadow: 1px 1px 2px ' + color_header_shadow + ';">'
  printf, 1, '                  <B>Complete packages</B><BR>'
  printf, 1, '                </TD>'
  printf, 1, '              </TR>'
  printf, 1, '              <TR>'
  printf, 1, '                <TD ALIGN="left" WIDTH=' + width_head + ' style="font-family:' + font_text + '; color:' + color_text + '; text-shadow: 1px 1px 1px ' + color_text_shadow + ';">'
  printf, 1, '                  The following packages are available with all required files contained in a single archive file.<BR>'
  printf, 1, '                  <UL>'
  ; Iterate through packages
  for i_package = 0, n_package - 1 do begin
    printf, 1, '                    <LI>'
    printf, 1, '                      <A NAME="PACKAGE_' + string_substitute( package_name[i_package], ' ', '_' ) + '"></A><FONT style="font-family:' + font_header + '; color:' + color_header + '; text-shadow: 1px 1px 1px ' + color_text_shadow + '; font-size : large"><B>' + package_name[i_package] + '</B></FONT><BR>'
    if package_info[i_package] ne '' then begin
      printf, 1, '                        ' + package_info[i_package] + '<BR>'
    endif
    if keyword_set( package_file_tar ) then begin
      printf, 1, '                        Click <A HREF="' + package_file_tar[i_package] + '"><B>here</B></A> for a .tar archive.'
      if keyword_set( package_file_tar ) then begin
        printf, 1, '                        <BR>'
      endif
    endif
    if keyword_set( package_file_zip ) then begin
      printf, 1, '                        Click <A HREF="' + package_file_zip[i_package] + '"><B>here</B></A> for a .zip archive.'
    endif
    if package_external[i_package] ne '' then begin
      printf, 1, '                        <BR>'
      printf, 1, '                        These dependencies external to IDL are also required:  ' + package_external[i_package] + '.'
    endif
    printf, 1, '                    </LI>'
  endfor
  printf, 1, '                  <UL>'
  ;printf, 1, '                  <BR>'
  printf, 1, '                </TD>'
  printf, 1, '              </TR>'
  printf, 1, '            </TABLE>'
  printf, 1, '          </TD>'
  printf, 1, '        </TR>'
  printf, 1, '      </TABLE>'
  printf, 1, '    </TD>'
  printf, 1, '  </TR>'
  printf, 1
endif

printf, 1, '  <!-- Table containing new news and links to the category lists /-->'
printf, 1, '  <TR>'
printf, 1, '    <TD ALIGN="center">'
printf, 1, '      <TABLE WIDTH=' + width_box_head + ' ALIGN="center">'
printf, 1, '        <TR>'
printf, 1, '          <TD VALIGN="top" WIDTH=' + width_head + ' ALIGN="center" STYLE="border:5px outset ' + color_box_border + '" BGCOLOR="' + color_box_background + '">'
printf, 1, '            <TABLE WIDTH=' + width_head + ' ALIGN="center">'
printf, 1, '              <TR>'
printf, 1, '                <TD>'
printf, 1, '                  <P STYLE="font-size:small;">'
printf, 1, '                    <BR>'
printf, 1, '                  </P>'
printf, 1, '                </TD>'
printf, 1, '              </TR>'
if keyword_set( first_update ) or keyword_set( update ) then begin
  printf, 1, '              <TR>'
  printf, 1, '                <TD ALIGN="center" WIDTH=' + width_head + ' style="font-family:' + font_header + '; color:' + color_header + '; font-size:x-large; font-weight:bold; text-shadow: 1px 1px 2px ' + color_header_shadow + ';">'
  printf, 1, '                  <B>News</B><BR>'
  printf, 1, '                </TD>'
  printf, 1, '              </TR>'
  printf, 1, '              <TR>'
  printf, 1, '                <TD ALIGN="left" WIDTH=' + width_head + ' style="font-family:' + font_text + '; color:' + color_text + '; text-shadow: 1px 1px 1px ' + color_text_shadow + ';">'
  printf, 1, '                  <UL>'
  for i_update = 0, n_elements( first_update ) - 1 do begin
    printf, 1, '                    <LI>' + first_update[i_update] + '</LI>'
  endfor
  if keyword_set( update ) then begin
    printf, 1, '                    <LI>See <A HREF="#NEWS"><B>bottom</B></A> for more update information.</LI>'
  endif
  printf, 1, '                  </UL>'
  printf, 1, '                </TD>'
  printf, 1, '              </TR>'
endif
printf, 1, '              <TR>'
printf, 1, '                <TD ALIGN="center" WIDTH=' + width_head + ' style="font-family:' + font_header + '; color:' + color_header + '; font-size:x-large; font-weight:bold; text-shadow: 1px 1px 2px ' + color_header_shadow + ';">'
printf, 1, '                  <B>Recent updates</B><BR>'
printf, 1, '                </TD>'
printf, 1, '              </TR>'
printf, 1, '              <TR>'
printf, 1, '                <TD ALIGN="left" WIDTH=' + width_head + ' style="font-family:' + font_text + '; color:' + color_text + '; text-shadow: 1px 1px 1px ' + color_text_shadow + ';">'
pro_date_diff = lonarr( n_pro )
for i_pro = 0, n_pro - 1 do begin
  pro_date_diff[i_pro] = julday( $
      fix( strmid( str( date_current_numeric ), 4, 2 ) ), $
      fix( strmid( str( date_current_numeric ), 6, 2 ) ), $
      fix( strmid( str( date_current_numeric ), 0, 4 ) ) ) $
      - julday( fix( strmid( str( pro_date[i_pro] ), 4, 2 ) ), $
      fix( strmid( str( pro_date[i_pro] ), 6, 2 ) ), $
      fix( strmid( str( pro_date[i_pro] ), 0, 4 ) ) )
endfor
id_new = where( pro_date_diff lt hist_delay, n_id_new )
if n_id_new gt 0 then begin
  id_sort = sort( pro_date_diff[id_new] )
  printf, 1, '                <UL>'
  for i_new = 0, n_id_new - 1 do begin
    temp = str( pro_date[id_new[id_sort[i_new]]] )
    printf, 1, '                  <LI>'
    printf, 1, '                    <A HREF="' + strlowcase( list[id_new[id_sort[i_new]]] ) + '"><B>' + pro_name[id_new[id_sort[i_new]]] + '</B></A> (' $
        + strmid( temp, 0, 4 ) + '-' + strmid( temp, 4, 2 ) + '-' $
        + strmid( temp, 6, 2 ) + ')<BR>'
  printf, 1, '                  </LI>'
  endfor
  printf, 1, '                </UL>'
endif
;printf, 1, '                  <BR>'
printf, 1, '                </TD>'
printf, 1, '              </TR>'
printf, 1, '            </TABLE>'
printf, 1, '          </TD>'
printf, 1, '        </TR>'
printf, 1, '      </TABLE>'
printf, 1, '    </TD>'
printf, 1, '  </TR>'
printf, 1

; Write categories to page
for i_categ = 0, n_categ-1 do begin
  printf, 1, '  <!-- Table containing the ' + categ_name[id_categ_sort[i_categ]] + ' category /-->'
  printf, 1, '  <TR>'
  printf, 1, '    <TD ALIGN="center">'
  printf, 1, '      <TABLE WIDTH=' + width_box_text + ' ALIGN="center">'
  printf, 1, '        <TR>'
  printf, 1, '          <TD VALIGN="top" WIDTH=' + width_box_text + ' ALIGN="center" STYLE="border:5px outset ' + color_box_border + '" BGCOLOR="' + color_box_background + '">'
  printf, 1, '            <TABLE WIDTH=' + width_text + ' ALIGN="center">'
  printf, 1, '              <TR>'
  printf, 1, '                <TD>'
  printf, 1, '                  <P STYLE="font-size:small;">'
  printf, 1, '                    <BR>'
  printf, 1, '                  </P>'
  printf, 1, '                </TD>'
  printf, 1, '              </TR>'
  printf, 1, '              <TR>'
  ; Category title
  printf, 1, '                <TD ALIGN="center" WIDTH=' + width_text + ' style="font-family:' + font_header + '; color:' + color_header + '; font-size:x-large; font-weight:bold; text-shadow: 1px 1px 2px ' + color_header_shadow + ';">'
  temp = string_substitute( categ_name[id_categ_sort[i_categ]], ' ', '_' )
  printf, 1, '                  <A NAME="' + temp + '"></A>' + categ_name[id_categ_sort[i_categ]]
  printf, 1, '                </TD>'
  printf, 1, '              </TR>'
  printf, 1, '              <TR>'
  printf, 1, '                <TD ALIGN="left" WIDTH=' + width_text + ' style="font-family:' + font_text + '; color:' + color_text + '; text-shadow: 1px 1px 1px ' + color_text_shadow + ';">'
  ; Program list
  printf, 1, '                  <UL>'
  temp = max( pro_categ_index eq id_categ_sort[i_categ], dimension=2 )
  id = where( temp eq 1, n_id )
  id_sort = sort( strlowcase( pro_name[id] ) )
  id = id[id_sort]
  for i_id = 0, n_id-1 do begin
    printf, 1, '                    <LI>'
    temp = '<A HREF="' + strlowcase( list[id[i_id]] ) + '"><B>' $
        + pro_name[id[i_id]] + '</B></A>'
    printf, 1, '                      ' + temp + '<BR>'
    printf, 1, '                      ' + pro_purp[id[i_id]] + '<BR>'
    temp = julday( fix( strmid( str( date_current_numeric ), 4, 2 ) ), $
        fix( strmid( str( date_current_numeric ), 6, 2 ) ), $
        fix( strmid( str( date_current_numeric ), 0, 4 ) ) ) $
        - julday( fix( strmid( str( pro_date[id[i_id]] ), 4, 2 ) ), $
        fix( strmid( str( pro_date[id[i_id]] ), 6, 2 ) ), $
        fix( strmid( str( pro_date[id[i_id]] ), 0, 4 ) ) )
    if temp lt hist_delay then begin
      printf, 1, '<B>' + pro_hist[id[i_id]] + '</B>'
    endif else begin
      printf, 1, pro_hist[id[i_id]]
    endelse
    printf, 1, '                    </LI>'
  endfor
  printf, 1, '                  <UL>'
  printf, 1, '                </TD>'
  printf, 1, '              </TR>'
  printf, 1, '            </TABLE>'
  printf, 1, '          </TD>'
  printf, 1, '        </TR>'
  printf, 1, '      </TABLE>'
  printf, 1, '    </TD>'
  printf, 1, '  </TR>'
  printf, 1
endfor

; List update information unless already done earlier in the webpage
if keyword_set( update ) then begin
  printf, 1, '  <!-- Table containing old news /-->'
  printf, 1, '  <TR>'
  printf, 1, '    <TD ALIGN="center">'
  printf, 1, '      <TABLE WIDTH=' + width_box_head + ' ALIGN="center">'
  printf, 1, '        <TR>'
  printf, 1, '          <TD VALIGN="top" WIDTH=' + width_head + ' ALIGN="center" STYLE="border:5px outset ' + color_box_border + '" BGCOLOR="' + color_box_background + '">'
  printf, 1, '            <TABLE WIDTH=' + width_head + ' ALIGN="center">'
  printf, 1, '              <TR>'
  printf, 1, '                <TD ALIGN="center" WIDTH=' + width_head + ' style="font-family:' + font_header + '; color:' + color_header + '; font-size:x-large; font-weight:bold; text-shadow: 1px 1px 2px ' + color_header_shadow + ';">'
  printf, 1, '                  <A NAME="NEWS"></A><B>News</B><BR>'
  printf, 1, '                </TD>'
  printf, 1, '              </TR>'
  printf, 1, '              <TR>'
  printf, 1, '                <TD ALIGN="left" WIDTH=' + width_head + ' style="font-family:' + font_text + '; color:' + color_text + '; text-shadow: 1px 1px 1px ' + color_text_shadow + ';">'
  printf, 1, '                  <UL>'
  for i_update = 0, n_elements( update ) - 1 do begin
    printf, 1, '                    <LI>' + update[i_update]
  endfor
  printf, 1, '                  </UL>'
  printf, 1, '                </TD>'
  printf, 1, '              </TR>'
endif
printf, 1, '                  <BR>'
printf, 1, '                </TD>'
printf, 1, '              </TR>'
printf, 1, '            </TABLE>'
printf, 1, '          </TD>'
printf, 1, '        </TR>'
printf, 1, '      </TABLE>'
printf, 1, '    </TD>'
printf, 1, '  </TR>'
printf, 1

printf, 1, '<!-- Row and column containing contact details /-->'
printf, 1, '  <TR>'
printf, 1, '    <TD ALIGN="center">'
printf, 1, '      <TABLE WIDTH=' + width_box_head + '>'
printf, 1, '        <TR>'
printf, 1, '          <TD BGCOLOR="' + color_footer_background + '" STYLE="border:5px outset ' + color_box_border + ';">'
printf, 1, '            <TABLE WIDTH=' + width_head + '>'
printf, 1, '              <TR>'
printf, 1, '                <TD WIDTH=' + str( fix( width_head ) / 2 - 2 * fix( width_margin ) ) + ' ALIGN="center" VALIGN="top" STYLE="font-family : ' + font_header + '; color : ' + color_footer_text + '; font-size : small; font-weight : normal; text-shadow: 1px 1px 1px ' + color_footer_shadow + ';">'
printf, 1, '                  <BR>'
printf, 1, '                  Maintained by: ' + admin_email + '<BR>'
printf, 1, '                  Last updated:  ' + date_current_written + '<BR><BR>'
printf, 1, '                </TD>'
printf, 1, '                <TD WIDTH=' + str( fix( width_head ) / 2 - 2 * fix( width_margin ) ) + ' ALIGN="center" VALIGN="top" STYLE="font-family : ' + font_header + '; color : ' + color_footer_text + '; font-size : small; font-weight : normal; text-shadow: 1px 1px 1px ' + color_footer_shadow + ';">'
printf, 1, '                  <BR>'
printf, 1, '                    &copy; Copyright 2017-' + strmid( date_current_numeric, 0, 4 ) + ' ' + admin_name + '<BR>'
printf, 1, '                </TD>'
printf, 1, '              </TR>'
printf, 1, '            </TABLE>'
printf, 1, '          </TD>'
printf, 1, '        </TR>'
printf, 1, '      </TABLE>'
printf, 1, '    </TD>'
printf, 1, '  </TR>'
printf, 1, '</TABLE>'

; Web page ending
;printf, 1, '<BR>'
printf, 1
printf, 1, '</CENTER>'
printf, 1, '</BODY>'
printf, 1, '</HTML>'

; Close web page file
close, 1

;***********************************************************************
; The End

;stop
return
END
