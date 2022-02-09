GLBVAR, var
var.dir = '../../dat/'

INPUT_FILE, jobname, prefname, firnum, endnum, njob

ivzeta = 2 

FOR ijob = 0, njob-1 DO BEGIN	; job loop

	prm_file = var.dir + prefname + '/' + jobname(ijob) + '.prm'
	READ_KEMPOPRM, prm_file, kmp

	nprev = (kmp.jobnum - 1) * kmp.nstep / kmp.jobnum
	ntime = (kmp.nstep - nprev) / ivzeta	; time steps in one job

	ddt = kmp.dt * ivzeta

	print, 'ntime = ', ntime
	print, 'part 1 -> 500000 : 750000'
	print, 'part 2 -> 937500 : 1125000'

	ntime1 = 0L  &  ntime2 = 0L
	read, 'start step = ', ntime1
	read, 'end   step = ', ntime2

	center = 16384  &  width = 256  &  position1 = 0  &  position2 = 0
	print, 'position -> nx/2-128 : nx/2+128 -> 16256 : 16512'
	read, 'center(default:16384) = ', center
	read, 'width(default:256) = ', width
	position1 = center - width/2
	position2 = center + width/2

	ncpu = kmp.ncpu
	head = var.head	; 4 bytes

	pos  = lonarr(ncpu)  &  pos[*] = 0L    ; position of file pointer
	tim  = dblarr(ncpu)  &  tim[*] = 0.0
	data = dblarr(5)	; array to store data -> t, vx, vy, vz, x

	sav  = strcompress('./particles' + string(ijob) + '.sav', /REMOVE_ALL)
	openw, 1, sav  &  point_lun, 1, 0

	for jt = ntime1, ntime2-1 do begin	; time loop

		print, 'time loop ', jt, '...'

		for icpu = 0, ncpu-1 do begin	; cpu loop

			file = strcompress(var.dir + prefname + '/' + jobname(ijob) + '.wpia' + string(icpu), /REMOVE_ALL)
			tmp = 0

			openr, 2, file  &  point_lun, 2, pos[icpu]*48  &  pos[icpu] = pos[icpu] - 1

			print, 'icpu = ', icpu

			while (1) do begin	; read loop

				if (tmp eq 2) then break

				readu, 2, head  &  readu, 2, data  &  readu, 2, head
				pos[icpu] = pos[icpu] + 1

				;if (data[0] lt ntime1*ddt) then continue

				if (data[4] eq 0) then begin
					tmp = tmp + 1
					continue
				endif

				if ((data[4] ge position1)  &&  data[4] le position2) then begin
					writeu, 1, data
					tim[icpu] = data[0]
				endif

				if (EOF(2)) then break

			endwhile
	
			close, 2
		endfor

		;if (NOT (jt mod 4000)) then print, 'time loop:', jt, '...'

	endfor	; time loop

	close, 1

ENDFOR	; job loop


end
