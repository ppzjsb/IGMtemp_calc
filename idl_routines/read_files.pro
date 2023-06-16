

f    = base+filename[ii]
openr,1,f

omegaM = 0.0d ;; matter (cdm+baryons)
omegaL = 0.0d ;; vacuum energy
omegaB = 0.0d ;; baryons
hubble = 0.0d ;; H0/100 km s^-1 Mpc^-1
Xh     = 0.0d ;; hydrogen fraction by mass
Delta  = 0.0d ;; normalised gas parcel density, rho/rho_crit
N      = 0L   ;; output rate

readu,1,omegaM,omegaL,omegaB,hubble,Xh,Delta,N

redshift = dblarr(N) ;; redshift
time     = dblarr(N) ;; time [s]
nH0      = dblarr(N) ;; H0 fraction, nH0/nH
nHep     = dblarr(N) ;; He+ fraction, nHe+/nH
nHepp    = dblarr(N) ;; He2+ fraction, nHe2+/nH
nelec    = dblarr(N) ;; electron fraction, ne/nH
temp     = dblarr(N) ;; temperature [K]
heat     = dblarr(N) ;; Heat/rho [erg s^-1 g^-1]
gH0      = dblarr(N) ;; HI ionisation rate  [s^-1]
gHe0     = dblarr(N) ;; HeI ionisation rate [s^-1]
gHep     = dblarr(N) ;; HeII ionisation rate [s^-1]
epsH0    = dblarr(N) ;; HI heating rate  [erg s^-1]
epsHe0   = dblarr(N) ;; HeI heating rate [erg s^-1]
epsHep   = dblarr(N) ;; HeII heating rate [erg s^-1]



readu,1,redshift,time,temp,nH0,nHep,nHepp,nelec,heat
readu,1,gH0,gHe0,gHep
readu,1,epsH0,epsHe0,epsHe
close,1

time /= timeGyr
heat *= (PROTONMASS/ELECTRONVOLT)
