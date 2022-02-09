
iediag = 256


time = Dblarr(1)
epe1 = time
epe2 = time
epa = time
bpe1 = time
bpe2 = time
bpa = time
ep1 = time
ep2 = time
ep3 = time
tskip = 1

; Get global variables
    GLBVAR, var

; Reading files
    INPUT_FILE, jobname, prefname, firnum, endnum, njob


    ps = 2
;    READ, 'Make PS file ?, y[1]/ n[2] : ', ps
    SET_DISP, ps, 'tx'

; Reading files
     jfrm = 0
     jjt = 0L
     jnt = 0L

     for ijob = 0, njob-1 do begin		; job loop

        file = '../../dat/' +prefname +'/' +jobName[ijob] +'.ene'

        ; Read input parameters in KEMPO from .input_idl file
        prm_file = '../../dat/' +prefname +'/'  + jobName(ijob) + '.prm'
        READ_KEMPOPRM, prm_file, kmp 

        nprev = (kmp.jobnum - 1) * kmp.nstep/ kmp.jobnum 
        ntime = (kmp.nstep - nprev)/ iediag 

        if (ijob eq 0) then begin
           total_tstep = ntime * njob / tskip
           time_ar = Dblarr(total_tstep)
           epe1_ar = time_ar
           epe2_ar = time_ar
           epa_ar = time_ar
           bpe1_ar = time_ar
           bpe2_ar = time_ar
           bpa_ar = time_ar
           ep1_ar = time_ar
           ep2_ar = time_ar
        endif 

        ;--- Open files ---
        print, 'Opening ', file   &   Openr, 1, file   &   Point_lun, 1, 0

        head = var.head

        for jt=0L, ntime-1 do begin

           readU, 1, head  &  readU, 1, time, epe1, epe2, epa, bpe1, bpe2, bpa, ep1, ep2, ep3  &  readU, 1, head

           time_ar[jjt] = time
           epe1_ar[jjt] = epe1
           epe2_ar[jjt] = epe2
           epa_ar[jjt] = epa
           bpe1_ar[jjt] = bpe1
           bpe2_ar[jjt] = bpe2
           bpa_ar[jjt] = bpa
           ep1_ar[jjt] = ep1
           ep2_ar[jjt] = ep2
           jjt = jjt + 1
        endfor
      close,1
    endfor


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Plot
    epe_ar = epe1_ar + epe2_ar + epa_ar
    bpe_ar = bpe1_ar + bpe2_ar ;+ bpa_ar
    eall_ar = epe_ar + bpe_ar + ep1_ar + ep2_ar

    !p.multi = [0,1,3]
    !P.charsize = 3     ; X

    color = ['red', 'blue', 'red', 'blue', 'purple']

    cgplot, time_ar, epe_ar, color=color[0], yrange=[1e-8,1e-1],/ylog
    cgplot, time_ar, bpe_ar, color=color[1], /overplot
    cglegend, alignment=2, colors=color[0], tcolors=color[0], location=[0.9, 0.8], titles='E field', thick=1
    cglegend, alignment=2, colors=color[1], tcolors=color[1], location=[0.9, 0.85], titles='B field', thick=1

    cgplot, time_ar, ep1_ar, yrange=[min(ep1_ar),max(ep1_ar)], color=color[2], title='P cold'
    cgplot, time_ar, ep2_ar, yrange=[min(ep2_ar),max(ep2_ar)], color=color[3], title='P hot', xtitle='time'

    Close, 1

    if (!d.name eq 'PS') then  device, /close

end
