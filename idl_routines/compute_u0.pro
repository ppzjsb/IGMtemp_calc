

u0  = dblarr(N) 
ind = where(heat gt 0.0 and redshift le ZMAX and redshift ge ZMIN) ;; ignore zero values

a = time[ind(0)]*timeGyr
b = time[max(ind)]*timeGyr
ab_plus  = 0.5d*(b+a)
ab_minus = 0.5d*(b-a)
eval   = ab_plus + ab_minus*absc[0: NWEIGHTS-1]
f_eval = interpol(heat,time*timeGyr,eval) 
u0     = total(ab_minus * weight[0 : NWEIGHTS-1] * f_eval) 
