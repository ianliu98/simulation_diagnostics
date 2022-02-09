f1 = './pti1_vx.csv'
aa = Fltarr(20480,6)
aa[*,0] = vxes[6,*]
aa[*,1] = vxes[10,*]
aa[*,2] = vxes[12,*]
aa[*,3] = vxes[13,*]
aa[*,4] = vxes[14,*]
aa[*,5] = vxes[15,*]
write_csv,f1,aa


f1 = './pti1_vy.csv'
aa = Fltarr(20480,6)
aa[*,0] = vyes[6,*]
aa[*,1] = vyes[10,*]
aa[*,2] = vyes[12,*]
aa[*,3] = vyes[13,*]
aa[*,4] = vyes[14,*]
aa[*,5] = vyes[15,*]
write_csv,f1,aa


f1 = './pti1_vz.csv'
aa = Fltarr(20480,6)
aa[*,0] = vzes[6,*]
aa[*,1] = vzes[10,*]
aa[*,2] = vzes[12,*]
aa[*,3] = vzes[13,*]
aa[*,4] = vzes[14,*]
aa[*,5] = vzes[15,*]
write_csv,f1,aa


f1 = './pti1_xe.csv'
aa = Fltarr(20480,6)
aa[*,0] = xxes[6,*]
aa[*,1] = xxes[10,*]
aa[*,2] = xxes[12,*]
aa[*,3] = xxes[13,*]
aa[*,4] = xxes[14,*]
aa[*,5] = xxes[15,*]
write_csv,f1,aa


end
