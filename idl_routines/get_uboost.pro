
;;----------------------------------------------------------------

;; Calculate the specific energy input associated with a given UV
;; background.  

;;----------------------------------------------------------------

pro get_uboost

base = '../output/'

ZMIN = 0.1d
ZMAX = 20.0d

NWEIGHTS = 10000


filename = ['puchwein2019.dat', $
            'puchwein2019_m8.0e-14_e5.0e-15.dat']
            


timeGyr      = 365.0 * 24.0 * 3600.0 * 1.0d9
MPC          = 3.08568025d24
ELECTRONVOLT = 1.602176565d-12
PROTONMASS   = 1.672621777d-24

u0_all = dblarr(n_elements(filename))


print
for ii=0, n_elements(filename)-1 do begin
   
@read_files
   
   if(ii eq 0) then begin
      zval_all = dblarr(N,n_elements(filename))
      temp_all = dblarr(N,n_elements(filename))
   endif
   
@gausslegendre
@compute_u0
   
   u0_all[ii]      = u0
   zval_all[*, ii] = redshift
   temp_all[*, ii] = temp
   
endfor


print
print,'Integration range',ZMAX,ZMIN
for ii=1, n_elements(filename)-1 do begin
print,filename[ii],': du [eV/mp]=',u0_all[ii]-u0_all[0]
endfor
print




end
