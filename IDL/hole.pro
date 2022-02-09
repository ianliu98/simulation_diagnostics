glbvar,var

INPUT_FILE, jobname, prefname, firnum, endnum, njob

ijob = 0
datjob = '../../dat/' + prefname + '/' + jobname[ijob]
prm_file = datjob +'.prm'
READ_KEMPOPRM, prm_file, kmp

fname = datjob + '.wpia100'

openr, 1, fname  &  point_lun, 1, 0
head = var.head

num = 100000
data = dblarr(5,num)
tmp = dblarr(5)

for i=0,num-1 do begin
	readu, 1, head  &  readu, 1, tmp  &  readu, 1, head
	data[*,i] = tmp	
endfor


end
