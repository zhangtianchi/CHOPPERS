close,/all




openr,1,'/data/halo_051/halo_051.0'
ng=0L
rhobin=0L
vcbin=0L
readu,1,ng
readu,1,rhobin
readu,1,vcbin



halo = {$
         pos          :fltarr(3), $
         vel          :fltarr(3), $
         sam          :fltarr(3), $    ; halo specific angular momentum
         cm           :fltarr(3), $    ; center of mass
         inertia      :fltarr(6), $    
         m200         :0.0, $
         n200         :0L, $
         r200         :0.0, $
         v200         :0.0, $
         vmax         :0.0, $
         rmax         :0.0, $
         rhalf        :0.0, $          ; half mass radius
         spinB        :0.0, $          ; Bullock 2001 spin parameter
         spinP        :0.0, $          ; Peebles 1969 spin Parameter
         rconv        :0.0, $          ; Power 2003 et al. convergence radius
         c            :0.0, $          ; fit NFW concentration
         cvmax        :0.0, $          ; vmax concentration
         vdisp        :0.0, $          ; velocity dispersion
         ea           :0.0, $          ; first axis of moment of inertia tensor
         eb           :0.0, $          ; second axis of moment of inertia tensor
         ec           :0.0, $          ; thrid axis of moment of inertia tensor
         pe           :0.0, $          ; specific potential energy
         ke           :0.0, $          ; specific kinetic energy
         fsub         :0.0, $          ; sub fraction
         soff         :0.0, $          ; center of mass displacement
         virial       :0.0, $          ;   2ke/|pe|
         denr         :fltarr(rhobin),$
         denrho       :fltarr(rhobin),$
         dennp        :lonarr(rhobin),$
         velr         :fltarr(vcbin),$
         velv         :fltarr(vcbin),$
         velnp        :lonarr(vcbin),$
         trackid      :0L, $
         birth        :0L $            ; birth snapshot
	 $
         }

tree=replicate(halo,ng)
readu,1,tree


print,'nhalo=',ng
print,'density_bin=',rhobin
print,'vc_bin=',vcbin

print, tree

;print,'pos=',tree.pos
;print,'trickid=',tree[0].trackid
;print,'m200=',tree[0].m200
;print,'lambda1=',tree.ea
;print,'lambda2=',tree.eb
;print,'lambda3=',tree.ec
;print,'vdisp=',tree.vdisp
;print,'rmax=',tree.rmax
;print,'Vmax=',tree.vmax
;print,'r200=',tree[0].r200
;print,'virial=',tree.virial
;print,'r/r200=',tree[0].denr
;print,'rho/rhoc=',tree[0].denrho
;print,tree.npvbin
;print,'c=',tree.c
;print,'cvmax=',tree.cvmax
close,1





end

