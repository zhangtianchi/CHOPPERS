close,/all

openr,1,'/data/inspur_disk01/userdir/tczhang/work4/grid/halo_051/halo_id_051.0'
Ngroups=0L
readu,1,Ngroups
print,'Nhalo=',Ngroups
Nid=0L
readu,1,Nid
print,'Nid=',Nid
Len=lonarr(Ngroups)
Off=lonarr(Ngroups)
ID=lon64arr(Nid)
readu,1, Len
readu,1, Off
print,'Len=',Len
print,'Offset=',Off
readu,1,ID
print,'PID=',ID(0:100)
close,1

end

