;
; PURPOSE:
;    return an array[column, row] of data that read from text file
;
; INPUT:
;    file: filename
;
; OUTPUT:
;    The double type array of data
;
; USE:
;    > a = getdata()
;
; Thanks:
;    coyote programming
;    http://www.idlcoyote.com/fileio_tips/read_columns.html
;
;                                      M.Hikishima  Mar 27, 2013@Sendai
;

FUNCTION  getdata
;PRO  getdata

   file = dialog_pickfile()

   ; Determine the number of rows in the file.
   rows = File_Lines(file)

   ; Determine the number of colums in the file by reading
   ; the first line and parsing it into column units.
   OpenR, lun, file, /Get_Lun
   line = ""
   ReadF, lun, line

   ; Find the number of columns in the line.
   cols = N_Elements(StrSplit(line, /RegEx, /Extract))

   ; Create a variable to hold the data.
;   data = FltArr(cols, rows)
   data = DblArr(cols, rows)

   ; if array dimension is [1, ??] -> [??]
   data = Reform(data)

   print, 'The array is : '  &  help, data

   ; Rewind the data file to its start.
   Point_Lun, lun, 0

   ; Read the data.
   ReadF, lun, data
   Free_Lun, lun

   Return, data

end
