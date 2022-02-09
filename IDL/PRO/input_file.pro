;+
; PURPOSE:
;    Give prefix names of input files
;
;                                                     M.Hikishima  Jun 1, 2012
;
; OUTPUT
;    jobname:  prefix of job name 
;    njob:     total numver of input jobs
;-


PRO input_file, jobname, filename_prefix, firnum, endnum, njob

print, ''
filename_prefix = ''     ; define file name as string
print, 'if test1.bx, input "test"'
read, 'prefix of file name (without number) ? : ', filename_prefix
read, 'first number of files ? : ', firnum
read, 'end number of files ? : ', endnum

njob = endnum - firnum + 1
fnum = String(Fix(firnum) + Indgen(njob))
jobname = Strcompress(filename_prefix + fnum, /remove_all)  ; delete some spaces
print, 'open files ... ', jobname


END
