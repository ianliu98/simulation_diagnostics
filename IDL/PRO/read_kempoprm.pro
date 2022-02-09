;+
; PURPOSE:
;    To obtain the following variables which are weitten in .input_idl as a structure
;
;                                                     M.Hikishima  May 31, 2012
; INPUT:
;    prmfile:  Name of .input_idl in dat directory
;
; OUTPUT
;    kmp:      Structure consisting of parameters 
;-



PRO read_kempoprm, prmfile, kmp


Get_lun, unit
Openr, unit, prmfile
Point_lun, unit, 0

; parameters to read
; see conv_input.f
prm_name = [  $
;
'ncpu',  $
;
'jobnum',  $
'nstep',  $
'ifdiag',  $
'ijdiag',  $
'ivdiag1',  $
'ivdiag2',  $
'ivdiag3',  $
'ipdig1',  $
'ipdig2',  $
'ips1',  $
'ips2',  $
'iwpidiag',  $
;
'wp1',  $
'wp2',  $
'wc',  $
'cv',  $
'dt',  $
'dr',  $
'nx',  $
'nspec',  $
'mltstp',  $
'juncan',  $
;
'np1',  $
'peth1',  $
'path1',  $
'np2',  $
'peth2',  $
'path2',  $
'amjz',  $
'omegjz',  $
'nxl',  $
'nxr',  $
'npin1',  $
'npin2',  $
'b02'  $
;
]

; all parameters are converted to double precision
nprm = N_elements(prm_name)
tmp = 0.0
for j=0, nprm-1 do begin
    Readf, unit, tmp
    if j eq 0 then  kmp = Create_struct(     prm_name[j], Double(tmp) )  $
              else  kmp = Create_struct(kmp, prm_name[j], Double(tmp) )
endfor

Free_lun, unit


END
