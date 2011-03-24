pro calc_tau_igm, tau=tau, z=z, RVIR=rvir

  c = x_constants()
  wrest = 1215.6701d

  ;;;;;;;;;;;;;;;;;
  ;; Cosmology
  if not keyword_set(z) then z = 3.
  h = 0.72
  Hubb = cosm_hubble(z, /w05map)
  nH_bar = 0.04 * (c.rhoc * h^2) * (1+z)^3 / c.mp ;; Should correct for He

  ;;
  nval = 1000L

  ;;;;;;;;;;;;;;;;;
  ;; Model the Halo with Barkana
  if not keyword_set(RFIL) then rfil = getenv('LYAP')+'/Analysis/IGM/Barkana/rho_12_z3.dat'
  if not keyword_set(VFIL) then vfil = getenv('LYAP')+'/Analysis/IGM/Barkana/vel_12_z3.dat'

  if not keyword_set(RVIR) then rvir = 50.  ; pMpc [Should be redshift dependent!]
  if not keyword_set(VC) then vc = 200. ; km/s
  r = rvir + 100*rvir*findgen(nval)/(nval-1)

  ;;;;;;;;;;;;;;;;;
  ;; Model density
  readcol, RFIL, B_r_rv, B_rho_rhob, format='F,F'
  rho_rhob = interpol(B_rho_rhob, B_r_rv, r/rvir)
;  x_splot, alog10(r/rvir), rho_rhob, /block

  nH_r = rho_rhob*nH_bar ; cm^-3

  ;;;;;;;;;;;;;;;;;
  ;; Velocity  (Following Fig 3 of Santos)
  ;; This includes Hubble flow and is *wrong* for z=3!!!
  readcol, VFIL, B_r_rv, B_velp_vc, format='F,F'
  velp_vc = interpol(B_velp_vc, B_r_rv, r/rvir)
  velr = velp_vc * vc - Hubb*r*(1+z)

  x_splot, alog10(r/rvir), velr, /block
  stop

;  dvdr = abs(v_1rvir / r)  ; km/s/kpc

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
  dr = r

  getfnam, wrest, fval, nam
  Kappa_l = !pi*c.e^2/c.me/c.c * fval

  tau_r = nHI_r*Kappa_l / (dvdr*1e5/c.kpc) * (wrest*1e-8) 
  stop

  return
end
