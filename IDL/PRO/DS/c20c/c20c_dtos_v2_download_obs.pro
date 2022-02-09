;+
; NAME:
;    c20c_dtos_v2_download_obs
;
; PURPOSE:
;    This procedure downloads observational data from an ftp or html address 
;    and converts non-NetCDF data to NetCDF format.
;
; CATEGORY:
;    C20C dtos v2
;
; CALLING SEQUENCE:
;    c20c_dtos_v2_download_obs, url
;
; INPUTS:
;    URL:  A required vector string containing the list of URLs of files to 
;        download.  Each file may contain data for multiple data variables.
;    CDO_SHIFTTIME, DIR_DEST, VAR_LABEL
;
; KEYWORD PARAMETERS:
;    CDO_SHIFTTIME:  An optional scalar string providing the setting for the 
;        CDO shifttime option, which shifts the time vector by the specified 
;        setting.  For instance '+15days' will shift the time vector forward by 
;        15 days.
;    DIR_DEST:  An optional vector string listing the directories to which to 
;        send the downloaded data files.  If not input then then current 
;        directory is used.  If a single value is input then all files are 
;        files are put in that directory.  If VAR_LABEL is set and the number 
;        of DIR_DEST and VAR_LABEL entries is equal, then data for each 
;        VAR_LABEL[i] is put in the DIR_DEST[i] directory.
;    VAR_LABEL:  An optional vector string listing the labels of the data 
;        variables to retain/extract from the downloaded files.  For example 
;        ['tos','sic'] will retain/extract the sea surface temperature and sea 
;        ice concentration variables.  If not input then all variables are 
;        retained (and not extracted).
;
; OUTPUTS:
;    ---
;
; USES:
;    cdo
;    ncdump
;    ncks
;    wget
;
; PROCEDURE:
;    This procedure uses wget to fetch requested files from somewhere on the 
;    internet.
;
; EXAMPLES:
;    ; The following downloads the July 2018 sea surface temperature and sea 
;    ; ice concentration data from the NOAA OI.v2 observationally-based 
;    ; product, and puts each variable in a different directory.
;    c20c_dtos_v2_download_obs, 'ftp://ftp.emc.ncep.noaa.gov/cmb/sst/oimonth_v2/GRIB/oiv2mon.201807.grb', var_label=['tos','sic'], dir_dest=['tos','sic'], cdo_shifttime='+15days'
;
; MODIFICATION HISTORY:
;    Written by:  Daithi Stone (dastone@runbox.com), 2017-11-29.
;    Modified:  DAS, 2018-08-23 (Completed header documentation.  Fixed various 
;        bugs with DIR_DEST and CDO_SHIFTTIME implementation)
;    Modified:  DAS, 2018-11-12 (Removed deletion of GRIB parameter table file)
;-

PRO C20C_DTOS_V2_DOWNLOAD_OBS, $
    URL, $
    DIR_DEST=dir_dest, $
    VAR_LABEL=var_label, $
    CDO_SHIFTTIME=cdo_shifttime

;***********************************************************************
; Constants and options

; Confirm a URL or series of URLs are requested
n_url = n_elements( url )
if n_url eq 0 then stop

; Default directory ($PWD) to which to download
n_dir_dest = n_elements( dir_dest )
if n_dir_dest eq 0 then begin
  dir_dest = ''
; Checks on input directories
endif else begin
  ; Iterate through defined directories
  for i_dir = 0, n_dir_dest - 1 do begin
    if dir_dest[i_dir] ne '' then begin
      ; Check that the directory(ies) exists
      spawn, 'cd ' + dir_dest[i_dir] + ' 2> /dev/null', exit_status=temp_status
      if temp_status ne 0 then begin
        ; The make the missing directory
        spawn, 'mkdir -p ' + dir_dest[i_dir], exit_status=temp_status
        if temp_status ne 0 then begin
          print, 'ERROR c20c_dtos_v2_download_obs:  Directory "' $
              + dir_dest[i_dir] + '" does not exist and unable to create.'
          stop
        endif
      endif
      ; Ensure the directory label ends with a slash
      if strmid( dir_dest[i_dir], strlen( dir_dest[i_dir] ) - 1, 1 ) ne '/' $
          then begin
        dir_dest[i_dir] = dir_dest[i_dir] + '/'
      endif
    endif
  endfor
endelse

; The default (lack of) list of requested variables
n_var_label = n_elements( var_label )
if n_var_label eq 0 then var_label = ''
; Confirm that dir_dest is either scalar or, for moving extracted files, of a 
; consistent size
if ( n_dir_dest ne 1 ) and ( n_dir_dest ne n_var_label ) then begin
  print, 'ERROR c20c_dtos_v2_download_obs:  DIR_DEST must be absent, scalar, ' $
      + 'or the same size as VAR_LABEL.'
  stop
endif

;***********************************************************************
; Prepare for extraction of variables

; If we are downloading GRIB files
if min( strpos( url, '.grb' ) ) gt 0 then begin
  ; Create a parameter_table.txt file for use by CDO
  file_parameter_table = 'GRIB_parameter_table.txt'
  openw, 1, dir_dest[0] + file_parameter_table
  temp_ctr = 0
  if ( max( var_label eq 'tos' ) eq 1 ) or ( n_var_label eq 0 ) then begin
    printf, 1, '11' + string(9B) + 'tos' + string(9B) $
        + 'Sea surface temperature [K]'
    temp_ctr = temp_ctr + 1
  endif
  if ( max( var_label eq 'sic' ) eq 1 ) or ( n_var_label eq 0 ) then begin
    printf, 1, '91' + string(9B) + 'sic' + string(9B) $
        + 'Sea ice area fraction [fraction]'
    temp_ctr = temp_ctr + 1
  endif
  close, 1
endif

;***********************************************************************
; Download the files, convert formats

; Iterate through URLs (files)
for i_url = 0, n_url - 1 do begin
  ; Download the file
  if dir_dest[0] ne '' then begin
    temp_command = 'cd ' + dir_dest[0] + ' ; '
  endif else begin
    temp_command = ''
  endelse
  temp_command = temp_command + 'wget ' + url[i_url]
  spawn, temp_command, exit_status=temp_status
  ; If we have successfully downloaded the file
  if temp_status eq 0 then begin
    ; Determine the downloaded file name
    url_file_name = strsplit( url[i_url], '/', extract=1, count=n_temp )
    url_file_name = url_file_name[n_temp-1]
    ; Determine the file type and base name
    url_file_type = strsplit( url_file_name, '.', extract=1, count=n_temp )
    if n_temp eq 1 then stop
    url_file_type = url_file_type[n_temp-1]
    pos = strpos( url_file_name, '.' + url_file_type, reverse_search=1 )
    out_file_base = strmid( url_file_name, 0, pos )
    ; If this is a .gz compressed file
    if url_file_type eq 'gz' then begin
      ; Uncompress the file
      spawn, 'gunzip ' + dir_dest[0] + url_file_name, $
          exit_status=temp_status
      ; If we did not manage to uncompress the file it may be because it is 
      ; incorrectly labeled as compressed (this happens sometimes for instance 
      ; with the HadISST1 update files
      if temp_status ne 0 then begin
        ; If it may in fact be an uncompressed .nc file
        if strspos( out_file_base, '.nc' ) gt 0 then begin
          ; Confirm that it is an uncompressed .nc file
          spawn, 'ncdump -h ' + dir_dest[0] + url_file_name $
              + ' &> /dev/null', $
              exit_status=temp_status
          ; If so then relabel it
          if temp_status eq 0 then begin
            ; Re-determine the file type and base name
            url_file_type = strsplit( out_file_base, '.', extract=1, $
                count=n_temp )
            if n_temp eq 1 then stop
            url_file_type = url_file_type[n_temp-1]
            pos = strpos( url_file_name, '.' + url_file_type, reverse_search=1 )
            out_file_base = strmid( url_file_name, 0, pos )
            ; Re-name the file
            spawn, 'mv ' + dir_dest[0] + url_file_name + ' ' $
                + dir_dest[0] + out_file_base + '.nc', $
                exit_status=temp_status
            if temp_status ne 0 then stop
          endif else begin
            print, 'ERROR c20c_dtos_v2_download_obs:  File ' + url_file_name $
                + ' is neither a .gz or a .nc file, as it claims to be.  ' $
                + 'I do not know what to do.'
            stop
          endelse
        ; Otherwise the solution is not obvious
        endif else begin
          print, 'ERROR c20c_dtos_v2_download_obs:  File ' + url_file_name $
              + ' is neither a .gz or a .nc file.  I do not know what to do.'
          stop
        endelse
      endif else begin
        ; Document the file name and extension change
        url_file_type = strsplit( out_file_base, '.', extract=1, count=n_temp )
        if n_temp eq 1 then stop
        url_file_type = url_file_type[n_temp-1]
        pos = strpos( url_file_name, '.' + url_file_type, reverse_search=1 )
        out_file_base = strmid( url_file_name, 0, pos )
      endelse
    endif
    ; If this is a NetCDF file
    if max( url_file_type eq [ 'nc', 'nc4' ] ) eq 1 then begin
      ; If variables are requested for extraction
      if n_var_label gt 0 then begin
        ; Get the list of variables available in this file
        ncdf_fileinq, dir_dest[0] + url_file_name, temp
        temp_var_label = temp.vars[*].name
        temp = 0
        ; Determine which requested variables are provided by this file
        index = isin( var_label, temp_var_label )
        id_temp_var = where( index eq 1, n_id_temp_var )
        ; If there are requested variables in this file
        if n_id_temp_var gt 0 then begin
          ; Iterate through requested variables
          for i_var = 0, n_id_temp_var - 1 do begin
            ; Extract the variable
            temp_command = 'ncks -v ' + temp_var_label[id_temp_var[i_var]] $
                + ' ' + url_file_name + ' ' + out_file_base + '_' $
                + temp_var_label[id_temp_var[i_var]] + '.nc'
            if dir_dest[0] eq '' then begin
              spawn, temp_command, exit_status=temp_status
            endif else begin
              spawn, 'cd ' + dir_dest[0] + ' ; ' + temp_command, $
                  exit_status=temp_status
            endelse
            if temp_status ne 0 then begin
              print, 'ERROR c20c_dtos_v2_download_obs:  ' $
                  + 'Could not extract variable ' $
                  + temp_var_label[id_temp_var[i_var]] + ' from file ' $
                  + url_file_name + '.'
              stop
            endif
          endfor
        endif
        ; Now we can delete the original file
        spawn, 'rm ' + dir_dest[0] + url_file_name
      endif
    ; If this is a grib file
    endif else if url_file_type eq 'grb' then begin
      ; Confirm that we entered something in the table
      if temp_ctr eq 0 then begin
        print, 'ERROR c20c_dtos_v2_download_obs:  ' $
            + 'Unknown variables requested for extraction from GRIB file.'
        stop
      endif
      ; Convert GRIB to NetCDF
      if dir_dest[0] ne '' then begin
        temp_command = 'cd ' + dir_dest[0] + ' ; '
      endif else begin
        temp_command = ''
      endelse
      temp_command = temp_command + 'cdo -f nc splitname -setpartab,' $
          + file_parameter_table
      if keyword_set( cdo_shifttime ) then begin
        temp_command = temp_command + ' -shifttime,' + cdo_shifttime
      endif
      temp_command = temp_command + ' ' + url_file_name + ' ' $
          + out_file_base + '_'
      spawn, temp_command, exit_status=temp_status
      if temp_status ne 0 then begin
        print, 'ERROR c20c_dtos_v2_download_obs:  ' $
            + 'Unable to extract data from GRIB file ' + url_file_name + '.'
        stop
      endif
      ;; Remove the original file and parameter table
      ;spawn, 'rm ' + dir_dest[0] + url_file_name + ' ' $
      ;    + dir_dest[0] + file_parameter_table
      ; Ensure expected file and variable labels
      if dir_dest[0] ne '' then begin
        temp_command = 'cd ' + dir_dest[0] + ' ; '
      endif else begin
        temp_command = ''
      endelse
      if max( var_label eq 'sic' ) then begin
        temp = file_search( dir_dest[0] + out_file_base + '_var91.nc', $
            count=n_temp )
        if n_temp eq 1 then begin
          spawn, temp_command + 'mv ' + out_file_base + '_var91.nc ' $
              + out_file_base + '_sic.nc'
          spawn, temp_command + 'ncrename -O -v var91,sic ' + out_file_base $
              + '_sic.nc'
        endif
      endif
      if max( var_label eq 'tos' ) then begin
        temp = file_search( dir_dest[0] + out_file_base + '_var11.nc', $
            count=n_temp )
        if n_temp eq 1 then begin
          spawn, temp_command + 'mv ' + out_file_base + '_var11.nc ' $
              + out_file_base + '_tos.nc'
          spawn, temp_command + 'ncrename -O -v var11,tos ' + out_file_base $
              + '_tos.nc'
        endif
      endif
    ; Otherwise we do not know what to do with this file
    endif else begin
      stop
    endelse
    ; If specific variables were to be extracted to other directories
    if n_dir_dest gt 1 then begin
      ; Iterate through requested variables
      for i_var = 1, n_var_label - 1 do begin
        ; Determine if we have extracted this variable
        temp = dir_dest[0] + out_file_base + '_' + var_label[i_var] $
            + '.nc'
        out_file = file_search( temp, count=n_out_file )
        if n_out_file eq 1 then begin
          ; Move the file to the other directory
          temp = 'mv ' + dir_dest[0] + out_file_base + '_' $
              + var_label[i_var] + '.nc ' + dir_dest[i_var] $
              + out_file_base + '_' + var_label[i_var] + '.nc '
          spawn, temp, exit_status=temp_status
          if temp_status ne 0 then begin
            print, 'ERROR c20c_dtos_v2_download_obs:  Unable to move file ' $
                + out_file_base + '_' + var_label[i_var] + '.nc to directory ' $
                + dir_dest[i_var] + '.'
            stop
          endif
        endif
      endfor
    endif
  ; If we are unable to download the file
  endif else begin
    print, 'ERROR c20c_dtos_v2_download_obs:  Unable to download ' $
        + url[i_url] + '.'
    stop
  endelse
endfor

;***********************************************************************
; The end

return
END

