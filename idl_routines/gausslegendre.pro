
;; ------------------------------------------------------------------------------

;; Gauss-Legendre quadrature.  Compute the weights and abscissa.

;; ------------------------------------------------------------------------------

absc   = dblarr(NWEIGHTS)
weight = dblarr(NWEIGHTS)

eps = 3.d-14
m   = (NWEIGHTS+1)/2
for i=1, m do begin
   z = cos(!dpi*(i - 0.25d)/(NWEIGHTS + 0.5d))
   
   noconv:
   
   p1 = Legendre(z, NWEIGHTS, /double)
   p2 = Legendre(z, NWEIGHTS-1, /double)
   pp = NWEIGHTS*(z*p1 - p2)/(z*z - 1.0d)
   z1 = z
   z  = z1-p1/pp
   
   if (abs(z-z1) gt eps) then goto, noconv
   
   absc(i-1)          = -1.0d * z
   absc(NWEIGHTS-i)   = z
   weight(i-1)        = 2.0d/((1.0d - z*z)*pp*pp)
   weight(NWEIGHTS-i) = weight(i-1)
endfor


;; ------------------------------------------------------------------------------
