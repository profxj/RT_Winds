pro calc_tau_igm, relvel, z=z, RVIR=rvir, TAU_IGM=tau_v

  c = x_constants()
  wrest = 1215.6701d

  ;;;;;;;;;;;;;;;;;
  ;; Cosmology
  if not keyword_set(z) then z = 3.
  h = 0.72
  Hubb = cosm_hubble(z, /w05map, /silen)
  nH_bar = 0.04 * (c.rhoc * h^2) * (1+z)^3 / c.mp ;; Should correct for He

  ;;
  nval = 1000L

  ;;;;;;;;;;;;;;;;;
  ;; Model the Halo with Barkana
  if not keyword_set(RFIL) then rfil = getenv('LYAP')+'/Analysis/IGM/Barkana/rho_12_z3.dat'
  if not keyword_set(VFIL) then vfil = getenv('LYAP')+'/Analysis/IGM/Barkana/vel_12_z3.dat'

  ;; From Michele
  if not keyword_set(RVIR) then rvir = 80.  ; pMpc [Should be redshift dependent!]
  if not keyword_set(VC) then vc = 230. ; km/s [This is virial velocity, not v_c]
  r = rvir + 100*rvir*findgen(nval)/(nval-1)

  ;;;;;;;;;;;;;;;;;
  ;; Model density
  readcol, RFIL, B_r_rv, B_rho_rhob, format='F,F', /silen
  rho_rhob = interpol(B_rho_rhob, B_r_rv, r/rvir)
;  x_splot, alog10(r/rvir), rho_rhob, /block

  nH_r = rho_rhob*nH_bar ; cm^-3

  ;;;;;;;;;;;;;;;;;
  ;; Velocity  (Following Fig 3 of Santos)
  ;; This includes Hubble flow and is *wrong* for z=3!!!
  readcol, VFIL, B_r_rv, B_velp_vc, format='F,F', /silen
  velp_vc = interpol(B_velp_vc, B_r_rv, r/rvir)
  velr = abs(velp_vc * vc) - Hubb*(r/1e3)  ; km/s
;  x_splot, r/rvir, velr, /block, xlog=1

  ;; dv/dr
  dr = r - shift(r,1)  
  dr[0] = dr[1]
  dv = velr - shift(velr,1)
  dv[0] = dv[1]
  dvdr = dv/dr  ; km/s/kpc

;  x_splot, r/rvir, dvdr, /block, xlog=1
;  stop

  ;;;;;;;;;;;;;;;;;
  ;; Ionization
  Gamma =  0.5d-12 ; s^-1  [Check!]
  T = 1e4 ; K -- Might be too cold
  beta = c.Ryd / (c.k * T)
  L = c.mp / (2*c.k*T)
  A = 2.1e-22 ; cm^-2
  alpha_A = (2*A)/sqrt(!pi*L) * beta * 0.5 * (1.735 + alog(beta) + 1./(6*beta))

  chi_r = Gamma / (nH_r * alpha_A)
  xHI_r = 1 + 0.5*chi_r - 0.5*sqrt(chi_r^2 + 4*chi_r)  ;; Santos 04

  nHI_r = nH_r * xHI_r ; cm^-3

  ;;;;;;;;;;;;;;;;;
  ;; Sobolev
  getfnam, wrest, fval, nam
  Kappa_l = !pi*c.e^2/c.me/c.c * fval

  tau_r = nHI_r*Kappa_l / (abs(dvdr)*1e5/c.kpc) * (wrest*1e-8) 
  ;x_splot, velr, tau_r, /block

  ;;;;;;;;;;;;
  ;; Interpolate
  npix = n_elements(relvel)
  tau_v = fltarr(npix)
  mn = min(velr, max=mx)
  gd = where(relvel GE MN and relvel LE mx, ngd)
  if ngd GT 0 then tau_v[gd] = interpol(tau_r, velr, relvel[gd])

  return

end
