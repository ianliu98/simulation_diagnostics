; loscone

LOADCT, 39
GLBVAR, var

once = 0
tskip = 1
diag  = 256

INPUT_FILE, jobname, prefname, firnum, endnum, njob

for ijob = 0, njob-1 do begin    ; job loop
  datjob   = '../../dat/' + prefname + '/' + jobname[ijob]
  file     = datjob + '.los'
  prm_file = datjob + '.prm'
  READ_KEMPOPRM, prm_file, kmp
  nprev = (kmp.jobnum - 1) * kmp.nstep /kmp.jobnum
  ntstep = (kmp.nstep - nprev) /diag
  
  if (once eq 0) then begin
    ;+++++++++++++++++++++++++++++++++++++++
    ;++++++++++ preset for arrays ++++++++++
    ;+++++++++++++++++++++++++++++++++++++++
    total_tstep = ntstep * njob/ tskip
    time        = Dblarr(total_tstep)
    edis10k_r   = dblarr(90,total_tstep)
    edis10k_l   = edis10k_r
    edis100k_r  = edis10k_r
    edis100k_l  = edis10k_r
    edis1M_l    = edis10k_r
    edis1M_r    = edis10k_r
    
    time_tmp         = Dblarr(1)
    edis10k_r_tmp    = dblarr(90)
    edis10k_l_tmp    = edis10k_r_tmp
    edis100k_l_tmp   = edis10k_r_tmp
    edis100k_r_tmp   = edis10k_r_tmp
    edis1M_l_tmp     = edis10k_r_tmp
    edis1M_r_tmp     = edis10k_r_tmp
    jjt         = 0L

    once = 1
  endif
  
  print, 'Opening ', file   &   Openr, 1, file   &   Point_lun, 1, 0
  head  = var.head
  
  for jt=0L, total_tstep-1 do begin
    readU, 1, head  &  readU, 1, time_tmp, edis10k_r_tmp, edis10k_l_tmp, edis100k_r_tmp, edis100k_l_tmp, edis1M_r_tmp, edis1M_l_tmp   &  readU, 1, head
    time[jjt]          =  time_tmp
    edis10k_r[*, jjt]  = edis10k_r_tmp
    edis10k_l[*, jjt]  = edis10k_l_tmp
    edis100k_r[*, jjt] = edis100k_r_tmp
    edis100k_l[*, jjt] = edis100k_l_tmp
    edis1M_r[*, jjt]   = edis1M_r_tmp
    edis1M_l[*, jjt]   = edis1M_l_tmp
    jjt = jjt + 1
  endfor

  close,1
endfor


end
