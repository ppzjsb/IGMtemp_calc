
;--------------------------------------------------------------------

;; Plots the temperature and ionisation fractions as a function of
;; redshift.

;--------------------------------------------------------------------


pro plot_tevol

;; Time in Gyr
timeGyr = (365.25*24.0*3600.0*1.0d9)
MPC     = 3.08568025e24
ELECTRONVOLT = 1.602176565d-12
PROTONMASS   = 1.672621777d-24
XH = 0.76d
yhelium = (1.0-XH)/(4.0*XH)


base = '../output/'

filename = ['puchwein2019.dat', $
            'puchwein2019_m8.0e-14_e5.0e-15.dat']
                  
   
print
for ii=0, n_elements(filename)-1 do begin
   
@read_files

   
   if(ii eq 0) then begin
      zval_all = dblarr(N,n_elements(filename))
      temp_all = dblarr(N,n_elements(filename))
      nH0_all  = dblarr(N,n_elements(filename))
      nHep_all = dblarr(N,n_elements(filename))
      
   endif
   
   zval_all[*, ii] = redshift
   temp_all[*, ii] = temp
   nH0_all[*, ii]  = nH0
   nHep_all[*, ii] = nHep
   
endfor
   
   
window,0,xsize=500,ysize=500
device,Retain=2,true_color=24,decomposed=0
!p.font=-1
plot,zval_all[*,0],alog10(temp_all[*,0]),xrange=[-0.1,15],yrange=[3.0,4.7],xstyle=1,ystyle=1,charsize=1.75,xtitle='redshift z',ytitle='log(T/K)'
loadct,5
oplot,zval_all[*,1],alog10(temp_all[*,1]),linestyle=2,color=100

;window,1,xsize=500,ysize=500
;device,Retain=2,true_color=24,decomposed=0
;!p.font=-1
;plot,zval_all[*,0],temp_all[*,1]/temp_all[*,0] ,xrange=[-0.1,15],yrange=[0.9,3.0],xstyle=1,ystyle=1,charsize=1.75,xtitle='redshift z',ytitle='T / Tdp'
;loadct,5

window,1,xsize=500,ysize=500
device,Retain=2,true_color=24,decomposed=0
!p.font=-1
plot,zval_all[*,0],alog10(nH0_all[*,0]),xrange=[0,15],yrange=[-7,0.5],xstyle=1,ystyle=1,charsize=1.75,xtitle='redshift z',ytitle='log(n/nH), [H0, He+]'
oplot,zval_all[*,1],alog10(nH0_all[*,1]),linestyle=2,color=100
oplot,zval_all[*,0],alog10(nHep_all[*,0]/yhelium),linestyle=0
oplot,zval_all[*,1],alog10(nHep_all[*,1]/yhelium),linestyle=2,color=100


end
