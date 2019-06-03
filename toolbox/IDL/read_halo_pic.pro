close,/all

openr,1,'/data/inspur_disk01/userdir/tczhang/work4/grid/halo_134/halo_pic_134.0'
ngrid=0L
readu,1,ngrid
print,"Ngrid= ", ngrid
lbox=0.0
readu,1,lbox
print,"Lbox= ", lbox
nhalo=0L
readu,1,nhalo
print,"Nhalo= ", nhalo

picture = {$
         sxy         :dblarr(ngrid,ngrid),$
         syz         :dblarr(ngrid,ngrid),$
         sxz         :dblarr(ngrid,ngrid)$
	       $
         }
halopic=replicate(picture,nhalo)

readu,1,halopic

close,1



set_plot, 'ps'
device, xsize=100, ysize=100,$
        filename='test.eps', font_size=12,/color
!P.BACKGROUND='FFFFFF'xl
!P.COLOR='000000’xl
!P.Font=-1
Mp=0.496298
loadct,1



ii=where(halopic[3].sxy eq 0,count)
halopic[3].sxy[ii]=0.0000001
;print,count/(ngrid^2)


halopic[3].sxy=alog10(halopic[3].sxy*Mp*10.^10./(lbox/ngrid)^3.0)
halopic[3].sxy=transpose(halopic[3].sxy)
cgimage,bytscl(halopic[3].sxy),pos=[0.1,0.1,0.9,0.9],/noerase,/KEEP_ASPECT_RATIO

device,/close
set_plot,'x'
end

